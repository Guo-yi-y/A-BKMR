#' There are some internal functions
#'
#' @param y outcome
#' @param x covariates
#' @param Z exposures
#' @return objects for consquent computation
#' @keywords internal
# model fit ---------------------------------------------------------------


makeKpart <- function(r, Z1, Z2 = NULL) {
  Z1r <- sweep(Z1, 2, sqrt(r), "*")
  if (is.null(Z2)) {
    Z2r <- Z1r
  } else {
    Z2r <- sweep(Z2, 2, sqrt(r), "*")
  }
  Kpart <- fields::rdist(Z1r, Z2r)^2
  Kpart
}
a_makeVcomps <- function(X =X, r, lambda, Z, data.comps) {
  if (is.null(data.comps$knots)) {
    Kpart <- makeKpart(r, Z)
    V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
    if (data.comps$nlambda == 2) {
      V <- V + lambda[2]*data.comps$crossTT
    }
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)
    logdetVinv <- -2*sum(log(diag(cholV)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
  } else {## predictive process approach
    ## note: currently does not work with random intercept model
    nugget <- 0.001
    n0 <- nrow(Z)
    n1 <- nrow(data.comps$knots)
    nall <- n0 + n1

    K1 <- py$exp_neg_mat(makeKpart(r, data.comps$knots))
    K10 <- py$exp_neg_mat(makeKpart(r, data.comps$knots, Z))
    Q <- K1 + diag(nugget, n1, n1)
    R <- Q + lambda[1]*py$tcrossprod_py(K10)

    cholQ <- chol(Q)
    cholR <- chol(R)
    Qinv <- chol2inv(cholQ)
    Rinv <- chol2inv(cholR)

    XVinv = py$comp_XVinv(X, lambda[1], K10, Rinv)
    logdetVinv <- 2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR)))
    Vcomps <- list(lambda = lambda[1], XVinv = XVinv, logdetVinv = logdetVinv, cholR = cholR, Q = Q, K10 = K10, Qinv = Qinv, Rinv = Rinv)
  }
  Vcomps
}



makeVcomps <- function(r, lambda, Z, data.comps) {
  if (is.null(data.comps$knots)) {
    Kpart <- makeKpart(r, Z)
    V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
    if (data.comps$nlambda == 2) {
      V <- V + lambda[2]*data.comps$crossTT
    }
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)
    logdetVinv <- -2*sum(log(diag(cholV)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
  } else {## predictive process approach
    ## note: currently does not work with random intercept model
    nugget <- 0.001
    n0 <- nrow(Z)
    n1 <- nrow(data.comps$knots)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(Z, data.comps$knots))
    # Kall <- exp(-Kpartall)
    # K0 <- Kall[1:n0, 1:n0 ,drop=FALSE]
    # K1 <- Kall[(n0+1):nall, (n0+1):nall ,drop=FALSE]
    # K10 <- Kall[(n0+1):nall, 1:n0 ,drop=FALSE]
    K1 <- py$exp_neg_mat(makeKpart(r, data.comps$knots))
    K10 <- py$exp_neg_mat(makeKpart(r, data.comps$knots, Z))
    Q <- K1 + diag(nugget, n1, n1)
    R <- Q + lambda[1]*py$tcrossprod_py(K10)
    cholQ <- chol(Q)
    cholR <- chol(R)
    Qinv <- chol2inv(cholQ)
    Rinv <- chol2inv(cholR)
    Vinv <- diag(1, n0, n0) - py$compute_Vinv(lambda[1], K10, Rinv)
    logdetVinv <- 2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv, cholR = cholR, Q = Q, K10 = K10, Qinv = Qinv, Rinv = Rinv)
  }
  Vcomps
}


a_beta.update <- function(X, XVinv, y, sigsq.eps) {
  Vbeta <- chol2inv(chol(XVinv %*% X))
  cholVbeta <- chol(Vbeta)
  betahat <- Vbeta %*% XVinv %*% y
  #set.seed(s)
  n01 <- rnorm(ncol(X))
  betahat + crossprod(sqrt(sigsq.eps)*cholVbeta, n01)
}

a_sigsq.eps.update <- function(y, X, beta, Vcomps, a.eps=1e-3, b.eps=1e-3) {
  mu <- y - X%*%beta

  muVinv = py$compute_muVinv(mu, Vcomps$lambda, Vcomps$K10, Vcomps$Rinv)


  #set.seed(s)
  prec.y <- rgamma(1, shape=a.eps + nrow(X)/2, rate=b.eps + 1/2*muVinv)
  1/prec.y
}

#' @importFrom truncnorm rtruncnorm

ystar.update <- function(y, X, beta, h) {
  mu <-  drop(h + X %*% beta)
  lower <- ifelse(y == 1, 0, -Inf)
  upper <- ifelse(y == 0, 0,  Inf)
  #set.seed(s)
  samp <- truncnorm::rtruncnorm(1, a = lower, b = upper, mean = mu, sd = 1)
  drop(samp)
}
#' @importFrom tmvtnorm rtmvnorm
ystar.update.noh <- function(y, X, beta, Vinv, ystar) {
  mu <-  drop(X %*% beta)
  lower <- ifelse(y == 1, 0, -Inf)
  upper <- ifelse(y == 0, 0,  Inf)
  #set.seed(s)
  samp <- tmvtnorm::rtmvnorm(1, mean = mu, H = Vinv, lower = lower, upper = upper, algorithm = "gibbs", start.value = ystar)
  #samp <- truncnorm::rtruncnorm(1, a = lower, b = upper, mean = mu, sd = 1)
  drop(samp)
}

a_r.update <- function(r, whichcomp, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, data.comps, control.params, rprop.gen, rprop.logdens, rprior.logdens, ...) {
  # r.params <- set.r.params(r.prior = control.params$r.prior, comp = whichcomp, r.params = control.params$r.params)
  r.params <- bkmr:::make_r_params_comp(control.params$r.params, whichcomp)
  rcomp <- unique(r[whichcomp])
  if(length(rcomp) > 1) stop("rcomp should only be 1-dimensional")

  ## generate a proposal
  rcomp.star <- rprop.gen(current = rcomp, r.params = r.params)
  lambda.star <- lambda
  delta.star <- delta
  move.type <- NA

  ## part of M-H ratio that depends on the proposal distribution
  negdifflogproposal <- -rprop.logdens(rcomp.star, rcomp, r.params = r.params) + rprop.logdens(rcomp, rcomp.star, r.params = r.params)

  ## prior distribution
  diffpriors <- rprior.logdens(rcomp.star, r.params = r.params) - rprior.logdens(rcomp, r.params = r.params)

  r.star <- r
  r.star[whichcomp] <- rcomp.star

  ## M-H step
  return(a_MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps))
}

a_rdelta.comp.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, ztest, data.comps, control.params, rprop.gen2, rprop.logdens1, rprior.logdens, rprior.logdens2, rprop.logdens2, rprop.gen1,...) { ## individual variable selection
  r.params <- control.params$r.params
  a.p0 <- control.params$a.p0
  b.p0 <- control.params$b.p0
  delta.star <- delta
  r.star <- r

  move.type <- ifelse(all(delta[ztest] == 0), 1, sample(c(1,2),1))
  move.prob <- ifelse(all(delta[ztest] == 0), 1, 1/2)
  if(move.type == 1) {
    comp <- ifelse(length(ztest) == 1, ztest, sample(ztest, 1))
    r.params <- bkmr:::set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

    delta.star[comp] <- 1 - delta[comp]
    move.prob.star <- ifelse(all(delta.star[ztest] == 0), 1, 1/2)
    r.star[comp] <- ifelse(delta.star[comp] == 0, 0, rprop.gen1(r.params = r.params))

    diffpriors <- (lgamma(sum(delta.star[ztest]) + a.p0) + lgamma(length(ztest) - sum(delta.star[ztest]) + b.p0) - lgamma(sum(delta[ztest]) + a.p0) - lgamma(length(ztest) - sum(delta[ztest]) + b.p0)) + ifelse(delta[comp] == 1, -1, 1)*with(list(r.sel = ifelse(delta[comp] == 1, r[comp], r.star[comp])), rprior.logdens(x = r.sel, r.params = r.params))

    negdifflogproposal <- -log(move.prob.star) + log(move.prob) - ifelse(delta[comp] == 1, -1, 1)*with(list(r.sel = ifelse(delta[comp] == 1, r[comp], r.star[comp])), rprop.logdens1(x = r.sel, r.params = r.params))

  } else if(move.type == 2) {
    comp <- ifelse(length(which(delta == 1)) == 1, which(delta == 1), sample(which(delta == 1), 1))
    r.params <- bkmr:::set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

    r.star[comp] <- rprop.gen2(current = r[comp], r.params = r.params)

    diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[comp], r.params = r.params)

    negdifflogproposal <- -rprop.logdens2(r.star[comp], r[comp], r.params = r.params) + rprop.logdens2(r[comp], r.star[comp], r.params = r.params)
  }

  lambda.star <- lambda

  ## M-H step
  return(a_MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps))
}

a_rdelta.group.update <- function(r, delta, lambda, y, X, beta, sigsq.eps, Vcomps, Z, ztest, data.comps, control.params, rprop.gen1, rprior.logdens, rprop.logdens1, rprop.gen2, rprop.logdens2,...) { ## grouped variable selection
  r.params <- control.params$r.params
  a.p0 <- control.params$a.p0
  b.p0 <- control.params$b.p0
  groups <- control.params$group.params$groups
  sel.groups <- control.params$group.params$sel.groups
  neach.group <- control.params$group.params$neach.group
  delta.star <- delta
  r.star <- r

  # if(length(mu.r) == 1) mu.r <- rep(mu.r, nz)
  # if(length(sigma.r) == 1) sigma.r <- rep(sigma.r, nz)

  delta.source <- sapply(sel.groups, function(x) ifelse(any(delta[which(groups == groups[x])] == 1), 1, 0))
  delta.source.star <- delta.source

  ## randomly select move type
  if(all(delta.source == 0)) {
    move.type <- 1
    move.prob <- 1
  } else if(length(which(neach.group > 1 & delta.source == 1)) == 0) {
    move.type <- sample(c(1, 3), 1)
    move.prob <- 1/2
  } else {
    move.type <- sample(1:3, 1)
    move.prob <- 1/3
  }
  # move.type <- ifelse(all(delta.source == 0), 1, ifelse(length(which(neach.group > 1 & delta.source == 1)) == 0, sample(c(1, 3), 1), sample(1:3, 1)))

  # print(move.type)

  if(move.type == 1) { ## randomly select a source and change its state (e.g., from being in the model to not being in the model)

    source <- sample(seq_along(delta.source), 1)
    source.comps <- which(groups == source)

    # r.params <- set.r.params(r.prior = control.params$r.prior, comp = source.comps, r.params = r.params)

    delta.source.star[source] <- 1 - delta.source[source]
    delta.star[source.comps] <- rmultinom(1, delta.source.star[source], rep(1/length(source.comps), length(source.comps)))
    move.prob.star <- ifelse(all(delta.source.star == 0), 1, ifelse(length(which(neach.group > 1 & delta.source.star == 1)) == 0, 1/2, 1/3))

    ## which component got switched
    comp <- ifelse(delta.source[source] == 1, source.comps[which(delta[source.comps] == 1)], source.comps[which(delta.star[source.comps] == 1)])
    r.params <- bkmr:::set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

    r.star[comp] <- ifelse(delta.star[comp] == 0, 0, rprop.gen1(r.params = r.params))

    # diffpriors <- ifelse(delta.source[source] == 1, log(length(sel.groups) - sum(delta.source) + b.p0) - log(sum(delta.source.star) + a.p0), log(sum(delta.source) + a.p0) - log(length(sel.groups) - sum(delta.source.star) + b.p0)) + ifelse(delta.source[source] == 1, 1, -1)*log(length(source.comps)) + ifelse(delta.source[source] == 1, -1, 1)*with(list(r.sel = ifelse(delta.source[source] == 1, r[source.comps][which(delta[source.comps] == 1)], r.star[source.comps][which(delta.star[source.comps] == 1)])), rprior.logdens(x = r.sel, r.params = r.params))
    diffpriors <- ifelse(delta.source[source] == 1, log(length(sel.groups) - sum(delta.source) + b.p0) - log(sum(delta.source.star) + a.p0), log(sum(delta.source) + a.p0) - log(length(sel.groups) - sum(delta.source.star) + b.p0)) + ifelse(delta.source[source] == 1, 1, -1)*log(length(source.comps)) + ifelse(delta.source[source] == 1, -1, 1)*with(list(r.sel = ifelse(delta.source[source] == 1, r[comp], r.star[comp])), rprior.logdens(x = r.sel, r.params = r.params))

    # negdifflogproposal <- -log(move.prob.star) + log(move.prob) -ifelse(delta.source[source] == 1, 1, -1)*(log(length(source.comps)) - with(list(r.sel = ifelse(delta.source[source] == 1, r[source.comps][which(delta[source.comps] == 1)], r.star[source.comps][which(delta.star[source.comps] == 1)])), rprop.logdens1(x = r.sel, r.params = r.params)))
    negdifflogproposal <- -log(move.prob.star) + log(move.prob) -ifelse(delta.source[source] == 1, 1, -1)*(log(length(source.comps)) - with(list(r.sel = ifelse(delta.source[source] == 1, r[comp], r.star[comp])), rprop.logdens1(x = r.sel, r.params = r.params)))

  } else if(move.type == 2) { ## randomly select a multi-component source that is in the model and change which component is included

    tmp <- which(neach.group > 1 & delta.source == 1)
    source <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))
    source.comps <- which(groups == source)

    oldcomp <- source.comps[delta[source.comps] == 1]
    tmp <- source.comps[delta[source.comps] == 0]
    comp <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))

    r.params.oldcomp <- bkmr:::set.r.params(r.prior = control.params$r.prior, comp = oldcomp, r.params = r.params)
    r.params <- bkmr:::set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

    delta.star[oldcomp] <- 0
    delta.star[comp] <- 1

    r.star[oldcomp] <- 0
    r.star[comp] <- rprop.gen1(r.params = r.params)

    diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[oldcomp], r.params = r.params.oldcomp)

    negdifflogproposal <- -rprop.logdens1(r.star[comp], r.params = r.params) + rprop.logdens1(r[oldcomp], r.params = r.params.oldcomp)

  } else if(move.type == 3) { ## randomly select a component that is in the model and update it
    tmp <- which(delta == 1)
    comp <- ifelse(length(tmp) == 1, tmp, sample(tmp, 1))

    r.params <- bkmr:::set.r.params(r.prior = control.params$r.prior, comp = comp, r.params = r.params)

    r.star[comp] <- rprop.gen2(current = r[comp], r.params = r.params)

    diffpriors <- rprior.logdens(r.star[comp], r.params = r.params) - rprior.logdens(r[comp], r.params = r.params)

    negdifflogproposal <- -rprop.logdens2(r.star[comp], r[comp], r.params = r.params) + rprop.logdens2(r[comp], r.star[comp], r.params = r.params)
  }

  lambda.star <- lambda

  ## M-H step
  return(a_MHstep(r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps,))
}

a_lambda.update <- function(r, delta, lambda, whichcomp=1, y, X, Z = Z, beta, sigsq.eps, Vcomps, data.comps, control.params) {
  lambda.jump <- control.params$lambda.jump[whichcomp]
  mu.lambda <- control.params$mu.lambda[whichcomp]
  sigma.lambda <- control.params$sigma.lambda[whichcomp]
  lambdacomp <- lambda[whichcomp]

  ## generate a proposal
  #set.seed(s)
  lambdacomp.star <- rgamma(1, shape=lambdacomp^2/lambda.jump^2, rate=lambdacomp/lambda.jump^2)
  r.star <- r
  delta.star <- delta
  move.type <- NA

  ## part of M-H ratio that depends on the proposal distribution
  negdifflogproposal <- -dgamma(lambdacomp.star, shape=lambdacomp^2/lambda.jump^2, rate=lambdacomp/lambda.jump^2, log=TRUE) + dgamma(lambdacomp, shape=lambdacomp.star^2/lambda.jump^2, rate=lambdacomp.star/lambda.jump^2, log=TRUE)

  ## prior distribution
  diffpriors <- dgamma(lambdacomp.star, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE) - dgamma(lambdacomp, shape=mu.lambda^2/sigma.lambda^2, rate=mu.lambda/sigma.lambda^2, log=TRUE)

  lambda.star <- lambda
  lambda.star[whichcomp] <- lambdacomp.star

  ## M-H step
  return(a_MHstep( r=r, lambda=lambda, lambda.star=lambda.star, r.star=r.star, delta=delta, delta.star=delta.star, y=y, X=X, Z=Z, beta=beta, sigsq.eps=sigsq.eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move.type=move.type, data.comps=data.comps))
}

a_MHstep <- function( r, lambda, lambda.star, r.star, delta, delta.star, y, X, Z, beta, sigsq.eps, diffpriors, negdifflogproposal, Vcomps, move.type, data.comps) {
  ## compute log M-H ratio
  Vcomps.star <- a_makeVcomps(X, r.star, lambda.star, Z, data.comps)
  mu <- y - X%*%beta
  stainv1 = py$compute_stainv(mu, Vcomps.star$lambda, Vcomps.star$K10, Vcomps.star$Rinv)
  stainv2 = py$compute_stainv(mu, Vcomps$lambda, Vcomps$K10, Vcomps$Rinv)
  diffliks <- 1/2*Vcomps.star$logdetVinv - 1/2*Vcomps$logdetVinv - 1/2/sigsq.eps*(stainv1 - stainv2 )%*%mu
  logMHratio <- diffliks + diffpriors + negdifflogproposal
  logalpha <- min(0,logMHratio)

  ## return value
  acc <- FALSE
  if( log(runif(1)) <= logalpha ) {
    r <- r.star
    delta <- delta.star
    lambda <- lambda.star
    Vcomps <- Vcomps.star
    acc <- TRUE
  }
  return(list(r=r, lambda=lambda, delta=delta, acc=acc, Vcomps=Vcomps, move.type=move.type))
}

h.update <- function(lambda, Vcomps, sigsq.eps, y, X, beta, r, Z, data.comps) {
  if (is.null(Vcomps)) {
    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
  }
  if(is.null(Vcomps$Q)) {
    Kpart <- makeKpart(r, Z)
    K <- py$exp_neg_mat(Kpart)
    Vinv <- Vcomps$Vinv
    lambda <- lambda[1] ## in case with random intercept (randint==TRUE), where lambda is 2-dimensional
    lamKVinv <- lambda*K%*%Vinv
    h.postmean <- lamKVinv%*%(y-X%*%beta)
    ##h.postvar <- sigsq.eps*lamKVinv
    h.postvar <- sigsq.eps*lambda*(K - lamKVinv%*%K)
    h.postvar.sqrt <- try(chol(h.postvar), silent=TRUE)
    if(inherits(h.postvar.sqrt, "try-error")) {
      sigsvd <- svd(h.postvar)
      h.postvar.sqrt <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    hsamp <- h.postmean + crossprod(h.postvar.sqrt, rnorm(length(h.postmean)))
    hcomps <- list(hsamp = hsamp)
  } else {
    h.star.postvar.sqrt <- sqrt(sigsq.eps*lambda)*forwardsolve(t(Vcomps$cholR), Vcomps$Q)
    h.star.postmean <- lambda[1]*(Vcomps$Q %*% Vcomps$Rinv) %*% (Vcomps$K10 %*% (y - X %*% beta))
    hsamp.star <- h.star.postmean + crossprod(h.star.postvar.sqrt, rnorm(length(h.star.postmean)))
    hsamp <- t(Vcomps$K10) %*% (Vcomps$Qinv %*% hsamp.star)
    hcomps <- list(hsamp = hsamp, hsamp.star = hsamp.star)
  }
  hcomps
}

newh.update <- function(Z, Znew, Vcomps, lambda, sigsq.eps, r, y, X, beta, data.comps) {

  if(is.null(data.comps$knots)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(Z, Znew))
    # Kmat <- exp(-Kpartall)
    # Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
    # Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    # Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
    Kmat1 <- exp(-makeKpart(r, Znew))
    Kmat10 <- exp(-makeKpart(r, Znew, Z))

    if(is.null(Vcomps)) {
      Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
    }
    Vinv <- Vcomps$Vinv

    lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
    Sigma.hnew <- lambda[1]*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
    mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
    root.Sigma.hnew <- try(chol(Sigma.hnew), silent=TRUE)
    if(inherits(root.Sigma.hnew, "try-error")) {
      sigsvd <- svd(Sigma.hnew)
      root.Sigma.hnew <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    hsamp <- mu.hnew + crossprod(root.Sigma.hnew, rnorm(n1))
  } else {
    n0 <- nrow(data.comps$knots)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(data.comps$knots, Znew))
    # Kmat <- exp(-Kpartall)
    # Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
    # Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    # Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
    Kmat10 <- py$exp_neg_mat(makeKpart(r, Znew, data.comps$knots))

    if(is.null(Vcomps)) {
      Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
      h.star.postvar.sqrt <- sqrt(sigsq.eps*lambda[1])*forwardsolve(t(Vcomps$cholR), Vcomps$Q)
      h.star.postmean <- lambda[1]*Vcomps$Q %*% Vcomps$Rinv %*% Vcomps$K10 %*% (y - X %*% beta)
      Vcomps$hsamp.star <- h.star.postmean + crossprod(h.star.postvar.sqrt, rnorm(length(h.star.postmean)))
    }
    hsamp <- Kmat10 %*% Vcomps$Qinv %*% Vcomps$hsamp.star
  }

  hsamp
}

## function to obtain posterior samples of h(znew) from fit of Bayesian kernel machine regression
predz.samps <- function(fit, Znew, verbose = TRUE) {
  if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
  if(inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
  Z <- fit$Z
  if(ncol(Z) != ncol(Znew)) {
    stop("Znew must have the same number of columns as Z")
  }

  hnew.samps <- sapply(1:fit$nsamp, function(s) {
    if(s%%(fit$nsamp/10)==0 & verbose) print(s)
    newh.update(Z = Z, Znew = Znew, Vcomps = NULL, lambda = fit$lambda[s], sigsq.eps = fit$sigsq.eps[s], r = fit$r[s,], y = fit$y, X = fit$X, beta = fit$beta[s,], data.comps = fit$data.comps)
  })
  rownames(hnew.samps) <- rownames(Znew)
  t(hnew.samps)
}

## function to approximate the posterior mean and variance as a function of the estimated tau, lambda, beta, and sigsq.eps
newh.postmean <- function(fit, Znew, sel) {
  if(is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
  if(inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)

  Z <- fit$Z
  X <- fit$X
  y <- fit$y
  data.comps <- fit$data.comps
  lambda <- colMeans(fit$lambda[sel, ,drop = FALSE])
  sigsq.eps <- mean(fit$sigsq.eps[sel])
  r <- colMeans(fit$r[sel,])
  beta <- colMeans(fit$beta[sel, ,drop=FALSE])

  if(is.null(data.comps$knots)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    Kpartall <- makeKpart(r, rbind(Z, Znew))
    Kmat <- exp(-Kpartall)
    Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
    Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]

    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
    Vinv <- Vcomps$Vinv

    lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
    Sigma.hnew <- lambda[1]*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
    mu.hnew <- lamK10Vinv %*% (y - X%*%beta)
  } else {
    n0 <- nrow(data.comps$knots)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    Kpartall <- makeKpart(r, rbind(data.comps$knots, Znew))
    # Kmat <- exp(-Kpartall)
    # Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
    # Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
    # Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
    Kmat1 <- exp(-makeKpart(r, Znew))
    Kmat10 <- exp(-makeKpart(r, Znew, data.comps$knots))

    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)

    Sigma.hnew <- lambda[1]*sigsq.eps*Kmat10 %*% Vcomps$Rinv %*% t(Kmat10)
    mu.hnew <- lambda[1]*Kmat10 %*% Vcomps$Rinv %*% Vcomps$K10 %*% (y - X%*%beta)
  }

  ret <- list(postmean = drop(mu.hnew), postvar = Sigma.hnew)
  ret
}


set.r.MH.functions = function (r.prior) {
  if (r.prior == "gamma") {
    rprior.logdens <- function(x, r.params) {
      mu.r <- r.params$mu.r
      sigma.r <- r.params$sigma.r
      dgamma(x, shape = mu.r^2/sigma.r^2, rate = mu.r/sigma.r^2,
             log = TRUE)
    }
    rprop.gen1 <- function(r.params) {
      r.muprop <- r.params$r.muprop
      r.jump <- r.params$r.jump1
      #set.seed(s)
      rgamma(1, shape = r.muprop^2/r.jump^2, rate = r.muprop/r.jump^2)
    }
    rprop.logdens1 <- function(x, r.params) {
      r.muprop <- r.params$r.muprop
      r.jump <- r.params$r.jump1
      dgamma(x, shape = r.muprop^2/r.jump^2, rate = r.muprop/r.jump^2,
             log = TRUE)
    }
    rprop.gen2 <- function(current, r.params) {
      r.jump <- r.params$r.jump2
      #set.seed(s)
      rgamma(1, shape = current^2/r.jump^2, rate = current/r.jump^2)
    }
    rprop.logdens2 <- function(prop, current, r.params) {
      r.jump <- r.params$r.jump2
      dgamma(prop, shape = current^2/r.jump^2, rate = current/r.jump^2,
             log = TRUE)
    }
    rprop.gen <- function(current, r.params) {
      r.jump <- r.params$r.jump
      #set.seed(s)
      rgamma(1, shape = current^2/r.jump^2, rate = current/r.jump^2)
    }
    rprop.logdens <- function(prop, current, r.params) {
      r.jump <- r.params$r.jump
      dgamma(prop, shape = current^2/r.jump^2, rate = current/r.jump^2,
             log = TRUE)
    }
  }
  if (r.prior == "invunif") {
    rprior.logdens <- function(x, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      ifelse(1/r.b <= x & x <= 1/r.a, -2 * log(x) - log(r.b -
                                                          r.a), log(0))
    }
    rprop.gen1 <- function(r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      #set.seed(s)
      1/runif(1, r.a, r.b)
    }
    rprop.logdens1 <- function(x, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      ifelse(1/r.b <= x & x <= 1/r.a, -2 * log(x) - log(r.b -
                                                          r.a), log(0))
    }
    rprop.gen2 <- function(current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump2
      #set.seed(s)
      truncnorm::rtruncnorm(1, a = 1/r.b, b = 1/r.a, mean = current,
                            sd = r.jump)
    }
    rprop.logdens2 <- function(prop, current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump2
      #set.seed(s)
      log(truncnorm::dtruncnorm(prop, a = 1/r.b, b = 1/r.a,
                                mean = current, sd = r.jump))
    }
    rprop.gen <- function(current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump
      #set.seed(s)
      truncnorm::rtruncnorm(1, a = 1/r.b, b = 1/r.a, mean = current,
                            sd = r.jump)
    }
    rprop.logdens <- function(prop, current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump
      #set.seed(s)
      log(truncnorm::dtruncnorm(prop, a = 1/r.b, b = 1/r.a,
                                mean = current, sd = r.jump))
    }
  }
  if (r.prior == "unif") {
    rprior.logdens <- function(x, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      dunif(x, r.a, r.b, log = TRUE)
    }
    rprop.gen1 <- function(r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      #set.seed(s)
      runif(1, r.a, r.b)
    }
    rprop.logdens1 <- function(x, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      dunif(x, r.a, r.b, log = TRUE)
    }
    rprop.gen2 <- function(current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump2
      #set.seed(s)
      truncnorm::rtruncnorm(1, a = r.a, b = r.b, mean = current,
                            sd = r.jump)
    }
    rprop.logdens2 <- function(prop, current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump2
      #set.seed(s)
      log(truncnorm::dtruncnorm(prop, a = r.a, b = r.b,
                                mean = current, sd = r.jump))
    }
    rprop.gen <- function(current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump
      #set.seed(s)
      truncnorm::rtruncnorm(1, a = r.a, b = r.b, mean = current,
                            sd = r.jump)
    }
    rprop.logdens <- function(prop, current, r.params) {
      r.a <- r.params$r.a
      r.b <- r.params$r.b
      r.jump <- r.params$r.jump
      #set.seed(s)
      log(truncnorm::dtruncnorm(prop, a = r.a, b = r.b,
                                mean = current, sd = r.jump))
    }
  }
  list(rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1,
       rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2,
       rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen,
       rprop.logdens = rprop.logdens)
}




# output ------------------------------------------------------------------

a_ComputePostmeanHnew = function (fit, y = NULL, Z = NULL, X = NULL, Znew = NULL, sel = NULL,
                                   method = "approx", data.comps) {
  if (method == "approx") {
    res <- a_ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z,
                                         X = X, Znew = Znew, sel = sel, data.comps = data.comps)
  }
  else if (method == "exact") {
    res <- ComputePostmeanHnew.exact(fit = fit, y = y, Z = Z,
                                     X = X, Znew = Znew, sel = sel)
  }
  res
}


a_ComputePostmeanHnew.approx = function (fit, y = NULL, Z = NULL, X = NULL, Znew = NULL, sel = NULL, data.comps)
{

  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }
  if (!is.null(Znew)) {
    if (is.null(dim(Znew)))
      Znew <- matrix(Znew, nrow = 1)
    if (inherits(Znew, "data.frame"))
      Znew <- data.matrix(Znew)
  }
  if (is.null(dim(X)))
    X <- matrix(X, ncol = 1)
  ests <- ExtractEsts(fit, sel = sel)
  sigsq.eps <- ests$sigsq.eps[, "mean"]
  r <- ests$r[, "mean"]
  beta <- ests$beta[, "mean"]
  lambda <- ests$lambda[, "mean"]
  if (fit$family == "gaussian") {
    ycont <- y
  }else if (fit$family == "binomial") {
    ycont <- ests$ystar[, "mean"]
  }

  Vcomps_sum = a_makeVcomps(X = X, r = r, lambda = lambda, Z = Z, data.comps = data.comps)


  # Kpart <- makeKpart(r, Z)
  # K <- exp(-Kpart)
  # V <- diag(1, nrow(Z), nrow(Z)) + lambda[1] * K
  # cholV <- chol(V)
  # Vinv <- chol2inv(cholV)
  if (!is.null(Znew)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(Z, Znew))
    # Kmat <- exp(-Kpartall)
    # Kmat0 <- Kmat[1:n0, 1:n0, drop = FALSE]
    # Kmat1 <- Kmat[(n0 + 1):nall, (n0 + 1):nall, drop = FALSE]
    # Kmat10 <- Kmat[(n0 + 1):nall, 1:n0, drop = FALSE]
    #Kmat1 = exp_neg_mat(makeKpart(r, Znew))
    #Kmat10 = exp_neg_mat(makeKpart(r, Znew, Z))

    Kmat1 <- exp(-makeKpart(r, Znew))
    Kmat10 <- exp(-makeKpart(r, Znew, data.comps$knots))
    #lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
    #Vinv1 = diag(1, n0, n0) - lambda[1]*t(Vcomps_sum$K10) %*% Vcomps_sum$Rinv %*% Vcomps_sum$K10
    #lamK10Vinv <- lambda[1] * as.bigq(Kmat10) - lambda[1]*lambda[1] * as.bigq(Kmat10) %*% as.bigq(t(Vcomps_sum$K10)) %*% as.bigq(Vcomps_sum$Rinv) %*% as.bigq(Vcomps_sum$K10)
    #lamK10Vinv <- lambda[1] * Kmat10 - lambda[1]*lambda[1] * Kmat10 %*% t(Vcomps_sum$K10) %*% Vcomps_sum$Rinv %*% Vcomps_sum$K10

    #lamK10Vinv = comp_lamK10Vinv(exp_neg_mat(makeKpart(r, Znew, Z)), lambda[1], Vcomps_sum$K10, Vcomps_sum$Rinv)
    #postvar <- lambda[1] * sigsq.eps * (Kmat1 - lamK10Vinv %*% t(Kmat10))
    postvar <- lambda[1] * sigsq.eps * Kmat10 %*% Vcomps_sum$Rinv %*% t(Kmat10)
    #postmean <- lamK10Vinv %*% (ycont - X %*% beta)
    postmean = lambda[1]*Kmat10 %*% Vcomps_sum$Rinv %*% Vcomps_sum$K10 %*% (ycont - X%*%beta)
  }
  else {
    lamKVinv <- lambda[1] * K %*% Vinv
    postvar <- lambda[1] * sigsq.eps * (K - lamKVinv %*% K)

    postmean <- lamKVinv %*% (ycont - X %*% beta)
  }
  ret <- list(postmean = drop(postmean), postvar = postvar)
  ret
}



a_SingVarIntSummary = function (whichz = 1, fit, y = NULL, Z = NULL, X = NULL, data.comps, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), method = "approx", sel = NULL,
                                 ...) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }
  q.fixed <- qs.fixed[1]
  point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
  point2[whichz] <- quantile(Z[, whichz], qs.diff[2])
  point1[whichz] <- quantile(Z[, whichz], qs.diff[1])
  newz.q1 <- rbind(point1, point2)
  q.fixed <- qs.fixed[2]
  point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
  point2[whichz] <- quantile(Z[, whichz], qs.diff[2])
  point1[whichz] <- quantile(Z[, whichz], qs.diff[1])
  newz.q2 <- rbind(point1, point2)
  if (method %in% c("approx", "exact")) {
    preds.fun <- function(znew) a_ComputePostmeanHnew(fit = fit,
                                                       y = y, Z = Z, X = X, data.comps = data.comps, Znew = znew, sel = sel, method = method)
    interactionSummary <- bkmr:::interactionSummary.approx
  }
  else {
    stop("method must be one of c('approx', 'exact')")
  }
  interactionSummary(newz.q1, newz.q2, preds.fun, ...)
}






a_PredictorResponseUnivarVar = function (whichz = 1, fit, y, Z, X, data.comps = data.comps, method = "approx", ngrid = 50,
                                          q.fixed = 0.5, sel = NULL, min.plot.dist = Inf, center = TRUE,
                                          z.names = colnames(Z), ...) {
  if (ncol(Z) < 2)
    stop("requires there to be at least 2 predictor variables")
  if (is.null(z.names)) {
    colnames(Z) <- paste0("z", 1:ncol(Z))
  }
  else {
    colnames(Z) <- z.names
  }
  ord <- c(whichz, setdiff(1:ncol(Z), whichz))
  z1 <- seq(min(Z[, ord[1]]), max(Z[, ord[1]]), length = ngrid)
  z.others <- lapply(2:ncol(Z), function(x) quantile(Z[, ord[x]],
                                                     q.fixed))
  z.all <- c(list(z1), z.others)
  newz.grid <- expand.grid(z.all)
  colnames(newz.grid) <- colnames(Z)[ord]
  newz.grid <- newz.grid[, colnames(Z)]
  if (!is.null(min.plot.dist)) {
    mindists <- rep(NA, nrow(newz.grid))
    for (i in seq_along(mindists)) {
      pt <- as.numeric(newz.grid[i, colnames(Z)[ord[1]]])
      dists <- fields::rdist(matrix(pt, nrow = 1), Z[,
                                                     colnames(Z)[ord[1]]])
      mindists[i] <- min(dists)
    }
  }
  if (method %in% c("approx", "exact")) {
    preds <- a_ComputePostmeanHnew(fit = fit, y = y, Z = Z, data.comps = data.comps,
                                    X = X, Znew = newz.grid, sel = sel, method = method)
    preds.plot <- preds$postmean
    se.plot <- sqrt(diag(preds$postvar))
  }
  else {
    stop("method must be one of c('approx', 'exact')")
  }
  if (center)
    preds.plot <- preds.plot - mean(preds.plot)
  if (!is.null(min.plot.dist)) {
    preds.plot[mindists > min.plot.dist] <- NA
    se.plot[mindists > min.plot.dist] <- NA
  }
  res <- dplyr::tibble(z = z1, est = preds.plot, se = se.plot)
}

a_VarRiskSummary = function (whichz = 1, fit, y = NULL, Z = NULL, X = NULL, data.comps, qs.diff = c(0.25,
                                                                                                     0.75), q.fixed = 0.5, method = "approx", sel = NULL, ...) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }
  point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
  point2[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile,
                          qs.diff[2])
  point1[whichz] <- apply(Z[, whichz, drop = FALSE], 2, quantile,
                          qs.diff[1])

  cc <- c(-1, 1)
  newz <- rbind(point1, point2)
  preds <- a_ComputePostmeanHnew(fit = fit,
                                  y = y, Z = Z, X = X, data.comps = data.comps, Znew = newz, sel = sel, method = method)
  diff <- drop(cc %*% preds$postmean)
  diff.sd <- drop(sqrt(cc %*% preds$postvar %*% cc))
  return(c(est = diff, sd = diff.sd, est_e = preds$postmean[2], sd_e = preds$postvar[4]^0.5))



}



a_PredictorResponseBivarPair_gy = function (fit, y = NULL, Z = NULL, X = NULL, data.comps, whichz1 = 1, whichz2 = 2,
                                          whichz3 = NULL, method = "approx", prob = 0.5, q.fixed = 0.5,
                                          sel = NULL, ngrid = 50, min.plot.dist = 0.5, center = TRUE,
                                          ...) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }
  if (ncol(Z) < 3)
    stop("requires there to be at least 3 Z variables")
  if (is.null(colnames(Z)))
    colnames(Z) <- paste0("z", 1:ncol(Z))
  if (is.null(whichz3)) {
    ord <- c(whichz1, whichz2, setdiff(1:ncol(Z), c(whichz1,
                                                    whichz2)))
  }
  else {
    ord <- c(whichz1, whichz2, whichz3, setdiff(1:ncol(Z),
                                                c(whichz1, whichz2, whichz3)))
  }
  z1 <- seq(min(Z[, ord[1]]), max(Z[, ord[1]]), length = ngrid)
  z2 <- seq(min(Z[, ord[2]]), max(Z[, ord[2]]), length = ngrid)
  z3 <- quantile(Z[, ord[3]], probs = prob)
  z.all <- c(list(z1), list(z2), list(z3))
  if (ncol(Z) > 3) {
    z.others <- lapply(4:ncol(Z), function(x) quantile(Z[,
                                                         ord[x]], q.fixed))
    z.all <- c(z.all, z.others)
  }
  newz.grid <- expand.grid(z.all)
  z1save <- newz.grid[, 1]
  z2save <- newz.grid[, 2]
  colnames(newz.grid) <- colnames(Z)[ord]
  newz.grid <- newz.grid[, colnames(Z)]
  if (!is.null(min.plot.dist)) {
    mindists <- rep(NA, nrow(newz.grid))
    for (k in seq_along(mindists)) {
      pt <- as.numeric(newz.grid[k, c(colnames(Z)[ord[1]],
                                      colnames(Z)[ord[2]])])
      dists <- fields::rdist(matrix(pt, nrow = 1), Z[,
                                                     c(colnames(Z)[ord[1]], colnames(Z)[ord[2]])])
      mindists[k] <- min(dists)
    }
  }
  if (method %in% c("approx", "exact")) {
    preds <- a_ComputePostmeanHnew(fit = fit, y = y, Z = Z, data.comps = data.comps,
                                    X = X, Znew = newz.grid, sel = sel, method = method)
    preds.plot <- preds$postmean
    se.plot <- sqrt(diag(preds$postvar))
  }
  else {
    stop("method must be one of c('approx', 'exact')")
  }
  if (center)
    preds.plot <- preds.plot - mean(preds.plot)
  if (!is.null(min.plot.dist)) {
    preds.plot[mindists > min.plot.dist] <- NA
    se.plot[mindists > min.plot.dist] <- NA
  }
  res <- dplyr::tibble(z1 = z1save, z2 = z2save, est = preds.plot,
                       se = se.plot)
}

