#' advanced BKMR method
#'
#' This packages achieves acceleration of BKMR and provide quantitative estimates
#' It includes three parts: weighted sampling, approximate gaussian process, and quantization

#' @import dplyr
#' @import magrittr
#' @import fields
#' @import bkmr
#' @import statmod
#' @import glue
#' @import reticulate
#' @importFrom MASS mvrnorm



#' @title Simulate data

#' @description used for generating simulated data
#'
#' @param n number of observations
#' @param M number of exposures
#' @param sigsq.true Variance of normally distributed residual error
#' @param hfun An integer from 1 to 3 identifying which predictor-response function to generate
#' @param beta.true coefficient of covariates
#' @param Zgen Method for generating the matrix Z of exposure variables, taking one of the values c("unif", "norm", "corr", "realistic")
#' @param ind select which predictor(s) will be included in the h function; how many predictors that can be included will depend on which h function is being used
#' @param family 	a description of the error distribution and link function to be used in the model. Currently implemented for gaussian and binomial families
#' @export
sim_data = function (n = 100, M = 5, sigsq.true = 0.5, hfun = 3, beta.true = 2,
                     Zgen = "norm", ind = 1:2, family = "gaussian")
{
  stopifnot(n > 0, M > 0, sigsq.true >= 0, family %in% c("gaussian",
                                                         "binomial"))
  if (family == "binomial") {
    sigsq.true <- 1
  }
  if (hfun == 1) {
    HFun <- HFun1
  }
  else if (hfun == 2) {
    HFun <- HFun2
  }
  else if (hfun == 3) {
    HFun <- bkmr:::HFun3
  }
  else {
    stop("hfun must be an integer from 1 to 3")
  }
  if (Zgen == "unif") {
    Z <- matrix(runif(n * M, -2, 2), n, M)
  }
  else if (Zgen == "norm") {
    Z <- matrix(rnorm(n * M), n, M)
  }
  else if (Zgen == "corr") {
    if (M < 3) {
      stop("M must be an integer > 2 for Zgen = 'corr'")
    }
    Sigma <- diag(1, M, M)
    Sigma[1, 3] <- Sigma[3, 1] <- 0.95
    Sigma[2, 3] <- Sigma[3, 2] <- 0.3
    Sigma[1, 2] <- Sigma[2, 1] <- 0.1
    Z <- mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }
  else if (Zgen == "realistic") {
    VarRealistic <- structure(c(0.72, 0.65, 0.45, 0.48, 0.08,
                                0.14, 0.16, 0.42, 0.2, 0.11, 0.35, 0.1, 0.11, 0.65,
                                0.78, 0.48, 0.55, 0.06, 0.09, 0.17, 0.2, 0.16, 0.11,
                                0.32, 0.12, 0.12, 0.45, 0.48, 0.56, 0.43, 0.11, 0.15,
                                0.23, 0.25, 0.28, 0.16, 0.31, 0.15, 0.14, 0.48, 0.55,
                                0.43, 0.71, 0.2, 0.23, 0.32, 0.22, 0.29, 0.14, 0.3,
                                0.22, 0.18, 0.08, 0.06, 0.11, 0.2, 0.95, 0.7, 0.45,
                                0.22, 0.29, 0.16, 0.24, 0.2, 0.13, 0.14, 0.09, 0.15,
                                0.23, 0.7, 0.8, 0.36, 0.3, 0.35, 0.13, 0.23, 0.17,
                                0.1, 0.16, 0.17, 0.23, 0.32, 0.45, 0.36, 0.83, 0.24,
                                0.37, 0.2, 0.36, 0.34, 0.25, 0.42, 0.2, 0.25, 0.22,
                                0.22, 0.3, 0.24, 1.03, 0.41, 0.13, 0.39, 0.1, 0.1,
                                0.2, 0.16, 0.28, 0.29, 0.29, 0.35, 0.37, 0.41, 0.65,
                                0.18, 0.3, 0.18, 0.16, 0.11, 0.11, 0.16, 0.14, 0.16,
                                0.13, 0.2, 0.13, 0.18, 0.6, 0.18, 0.13, 0.08, 0.35,
                                0.32, 0.31, 0.3, 0.24, 0.23, 0.36, 0.39, 0.3, 0.18,
                                0.79, 0.42, 0.12, 0.1, 0.12, 0.15, 0.22, 0.2, 0.17,
                                0.34, 0.1, 0.18, 0.13, 0.42, 1.27, 0.1, 0.11, 0.12,
                                0.14, 0.18, 0.13, 0.1, 0.25, 0.1, 0.16, 0.08, 0.12,
                                0.1, 0.67), .Dim = c(13L, 13L))
    if (M > ncol(VarRealistic)) {
      stop("Currently can only generate exposure data based on a realistic correlation structure with M = 13 or fewer. Please set M = 13 or use Zgen = c('unif','norm'")
    }
    else if (M <= 13) {
      Sigma <- VarRealistic[1:M, 1:M]
    }
    Z <- mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }
  colnames(Z) <- paste0("z", 1:M)
  X1 <- cbind(3 * cos(Z[, 1]) + 2 * rnorm(n))
  X2 = rnorm(n)
  X3 = rbinom(n, 1, 0.3)
  eps <- rnorm(n, sd = sqrt(sigsq.true))
  h <- apply(Z, 1, HFun)
  mu <- X1 * beta.true +X2+1.2*X3+ h
  y <- drop(mu + eps)
  if (family == "binomial") {
    ystar <- y
    y <- ifelse(ystar > 0, 1, 0)
  }
  dat <- list(n = n, M = M, sigsq.true = sigsq.true, beta.true = beta.true,
              Z = Z, h = h, X = data.frame(X1 = X1, X2 = X2, X3 = X3), y = y, hfun = hfun, HFun = HFun,
              family = family)
  if (family == "binomial") {
    dat$ystar <- ystar
  }
  dat
}

# weighted sampling ------------------------------------------------------------

#' used for knot sampling
#' @param R exposures
#' @param nd number of the representative observations
#' @param num_nn number of nearest neighbor
#' @param P computation parameter
#' @param Q computation parameter
#' @param w whether the weight is used
#' @param max_loop max number of loop finding knots
#' @export


sam_py_r = function(R, nd, num_nn, P = -20, Q = 20, max_loop = 20, w = F){
  nd = as.integer(nd)
  num_nn = as.integer(num_nn)

  R_ic <- data.frame(R) %>%
    mutate(id = match(do.call(paste, as.list(.)), unique(do.call(paste, as.list(.)))))


  R_uni = data.frame(R_ic) %>% group_by(across(everything())) %>%
    summarise(count = n(), .groups = 'drop') %>% ungroup() %>% arrange(id)

  py <- get_py()
  if (w == F){
    id_ini = sample(unique(R_uni$id), nd)-1
    return(sam_py(R, nd, as.integer(id_ini), num_nn, P=-20, Q=20, max_loop=20))
  } else {
    id_ini = sample(unique(R_uni$id), nd, prob = R_uni$count/sum(R_uni$count))-1
    return(sam_py_w(R, nd, as.integer(id_ini), num_nn, P=-20, Q=20, max_loop=20))
  }


}

# model fit ---------------------------------------------------------------



#' advanced BKMR model fit
#' the parameters in function a_kmbayes is same as those in kmbayes in packages "bkmr"
#' @param y a vector of outcome data of length n
#' @param Z an n-by-M matrix of predictor variables to be included in the h function. Each row represents an observation and each column represents an predictor
#' @param X covariates
#' @param iter number of iterations to run the sampler
#' @param family a description of the error distribution and link function to be used in the model. Currently implemented for gaussian and binomial families
#' @param id optional vector (of length n) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting
#' @param Znew optional matrix of new predictor values at which to predict h, where each row represents a new observation. This will slow down the model fitting, and can be done as a post-processing step using SamplePred
#' @param starting.values list of starting values for each parameter. If not specified default values will be chosen.
#' @param control.params list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen
#' @param varsel TRUE or FALSE
#' @param groups optional vector (of length M) of group indicators for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed
#' @param knots optional matrix of knot locations for implementing the Gaussian predictive process of Banerjee et al. (2008). Currently only implemented for models without a random intercept
#' @param ztest optional vector indicating on which variables in Z to conduct variable selection (the remaining variables will be forced into the model)
#' @param rmethod for those predictors being forced into the h function, the method for sampling the r[m] values. Takes the value of 'varying' to allow separate r[m] for each predictor; 'equal' to force the same r[m] for each predictor; or 'fixed' to fix the r[m] to their starting values
#' @param est.h TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting
#' @export

a_kmbayes = function (y, Z, X = NULL, iter = 1000, family = "gaussian", id = NULL,
                       verbose = TRUE, Znew = NULL, starting.values = NULL, control.params = NULL,
                       varsel = FALSE, groups = NULL, knots = NULL, ztest = NULL,
                       rmethod = "varying", est.h = TRUE) {
  missingX <- is.null(X)
  if (missingX) X <- matrix(0, length(y), 1)
  hier_varsel <- !is.null(groups)

  ##Argument check 1, required arguments without defaults
  ##check vector/matrix sizes
  stopifnot (length(y) > 0, is.numeric(y), anyNA(y) == FALSE)
  if (!inherits(Z, "matrix"))  Z <- as.matrix(Z)
  stopifnot (is.numeric(Z), nrow(Z) == length(y), anyNA(Z) == FALSE)
  if (!inherits(X, "matrix"))  X <- as.matrix(X)
  stopifnot (is.numeric(X), nrow(X) == length(y), anyNA(X) == FALSE)

  ##Argument check 2: for those with defaults, write message and reset to default if invalid
  if (iter < 1) {
    message ("invalid input for iter, resetting to default value 1000")
    nsamp <- 1000
  } else {
    nsamp <- iter
  }
  if (!family %in% c("gaussian", "binomial")) {
    stop("family", family, "not yet implemented; must specify either 'gaussian' or 'binomial'")
  }
  if (family == "binomial") {
    message("Fitting probit regression model")
    if (!all(y %in% c(0, 1))) {
      stop("When family == 'binomial', y must be a vector containing only zeros and ones")
    }
  }
  if (rmethod != "varying" & rmethod != "equal" & rmethod != "fixed") {
    message ("invalid value for rmethod, resetting to default varying")
    rmethod <- "varying"
  }
  if (verbose != FALSE & verbose != TRUE) {
    message ("invalid value for verbose, resetting to default FALSE")
    verbose <- FALSE
  }
  if (varsel != FALSE & varsel != TRUE) {
    message ("invalid value for varsel, resetting to default FALSE")
    varsel <- FALSE
  }

  ##Argument check 3: the rest id (below) znew, knots, groups, ztest
  if (!is.null(id)) {
    stopifnot(length(id) == length(y), anyNA(id) == FALSE)
    if (!is.null(knots)) {
      message ("knots cannot be specified with id, resetting knots to null")
      knots<-NA
    }
  }
  if (!is.null(Znew)) {
    if (!inherits(Znew, "matrix"))  Znew <- as.matrix(Znew)
    stopifnot(is.numeric(Znew), ncol(Znew) == ncol(Z), anyNA(Znew) == FALSE)
  }
  if (!is.null(knots)) {
    if (!inherits(knots, "matrix"))  knots <- as.matrix(knots)
    stopifnot(is.numeric(knots), ncol(knots )== ncol(Z), anyNA(knots) == FALSE)
  }
  if (!is.null(groups)) {
    if (varsel == FALSE) {
      message ("groups should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(groups), length(groups) == ncol(Z), anyNA(groups) == FALSE)
    }
  }
  if (!is.null(ztest)) {
    if (varsel == FALSE) {
      message ("ztest should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(ztest), length(ztest) <= ncol(Z), anyNA(ztest) == FALSE, max(ztest) <= ncol(Z) )
    }
  }

  ## start JB code
  if (!is.null(id)) { ## for random intercept model
    randint <- TRUE
    id <- as.numeric(as.factor(id))
    nid <- length(unique(id))
    nlambda <- 2

    ## matrix that multiplies the random intercept vector
    TT <- matrix(0, length(id), nid)
    for (i in 1:nid) {
      TT[which(id == i), i] <- 1
    }
    crossTT <- tcrossprod(TT)
    rm(TT, nid)
  } else {
    randint <- FALSE
    nlambda <- 1
    crossTT <- 0
  }
  data.comps <- list(randint = randint, nlambda = nlambda, crossTT = crossTT)
  if (!is.null(knots)) data.comps$knots <- knots
  rm(randint, nlambda, crossTT)

  ## create empty matrices to store the posterior draws in
  chain <- list(h.hat = matrix(0, nsamp, nrow(Z)),
                beta = matrix(0, nsamp, ncol(X)),
                lambda = matrix(NA, nsamp, data.comps$nlambda),
                sigsq.eps = rep(NA, nsamp),
                r = matrix(NA, nsamp, ncol(Z)),
                acc.r = matrix(0, nsamp, ncol(Z)),
                acc.lambda = matrix(0, nsamp, data.comps$nlambda),
                delta = matrix(1, nsamp, ncol(Z))
  )
  if (varsel) {
    chain$acc.rdelta <- rep(0, nsamp)
    chain$move.type <- rep(0, nsamp)
  }
  if (family == "binomial") {
    chain$ystar <- matrix(0, nsamp, length(y))
  }

  ## components to predict h(Znew)
  if (!is.null(Znew)) {
    if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow=1)
    if (inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
    if (ncol(Z) != ncol(Znew)) {
      stop("Znew must have the same number of columns as Z")
    }
    ##Kpartall <- as.matrix(dist(rbind(Z,Znew)))^2
    chain$hnew <- matrix(0,nsamp,nrow(Znew))
    colnames(chain$hnew) <- rownames(Znew)
  }

  ## components if model selection is being done
  if (varsel) {
    if (is.null(ztest)) {
      ztest <- 1:ncol(Z)
    }
    rdelta.update <- a_rdelta.comp.update
  } else {
    ztest <- NULL
  }

  ## control parameters
  control.params.default <- list(lambda.jump = rep(10, data.comps$nlambda), mu.lambda = rep(10, data.comps$nlambda), sigma.lambda = rep(10, data.comps$nlambda), a.p0 = 1, b.p0 = 1, r.prior = "invunif", a.sigsq = 1e-3, b.sigsq = 1e-3, mu.r = 5, sigma.r = 5, r.muprop = 1, r.jump = 0.1, r.jump1 = 2, r.jump2 = 0.1, r.a = 0, r.b = 100)
  if (!is.null(control.params)){
    control.params <- modifyList(control.params.default, as.list(control.params))
    validateControlParams(varsel, family, id, control.params)
  } else {
    control.params <- control.params.default
  }

  control.params$r.params <- with(control.params, list(mu.r = mu.r, sigma.r = sigma.r, r.muprop = r.muprop, r.jump = r.jump, r.jump1 = r.jump1, r.jump2 = r.jump2, r.a = r.a, r.b = r.b))

  ## components if grouped model selection is being done
  if (!is.null(groups)) {
    if (!varsel) {
      stop("if doing grouped variable selection, must set varsel = TRUE")
    }
    rdelta.update <- a_rdelta.group.update
    control.params$group.params <- list(groups = groups, sel.groups = sapply(unique(groups), function(x) min(seq_along(groups)[groups == x])), neach.group = sapply(unique(groups), function(x) sum(groups %in% x)))
  }

  ## specify functions for doing the Metropolis-Hastings steps to update r
  e <- environment()
  rfn <- set.r.MH.functions(r.prior = control.params$r.prior)
  rprior.logdens <- rfn$rprior.logdens
  environment(rprior.logdens) <- e
  rprop.gen1 <- rfn$rprop.gen1
  environment(rprop.gen1) <- e
  rprop.logdens1 <- rfn$rprop.logdens1
  environment(rprop.logdens1) <- e
  rprop.gen2 <- rfn$rprop.gen2
  environment(rprop.gen2) <- e
  rprop.logdens2 <- rfn$rprop.logdens2
  environment(rprop.logdens2) <- e
  rprop.gen <- rfn$rprop.gen
  environment(rprop.gen) <- e
  rprop.logdens <- rfn$rprop.logdens
  environment(rprop.logdens) <- e
  rm(e, rfn)

  ## initial values
  starting.values0 <- list(h.hat = 1, beta = NULL, sigsq.eps = NULL, r = 1, lambda = 10, delta = 1)
  if (is.null(starting.values)) {
    starting.values <- starting.values0
  } else {
    starting.values <- modifyList(starting.values0, starting.values)
    validateStartingValues(varsel, y, X, Z, starting.values, rmethod)
  }
  if (family == "gaussian") {
    if (is.null(starting.values$beta) | is.null(starting.values$sigsq.eps)) {
      lmfit0 <- lm(y ~ Z + X)
      if (is.null(starting.values$beta)) {
        coefX <- coef(lmfit0)[grep("X", names(coef(lmfit0)))]
        starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
      }
      if (is.null(starting.values$sigsq.eps)) {
        starting.values$sigsq.eps <- summary(lmfit0)$sigma^2
      }
    }
  } else if (family == "binomial") {
    starting.values$sigsq.eps <- 1 ## always equal to 1
    if (is.null(starting.values$beta) | is.null(starting.values$ystar)) {
      probitfit0 <- try(glm(y ~ Z + X, family = binomial(link = "probit")))
      if (!inherits(probitfit0, "try-error")) {
        if (is.null(starting.values$beta)) {
          coefX <- coef(probitfit0)[grep("X", names(coef(probitfit0)))]
          starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
        }
        if (is.null(starting.values$ystar)) {
          #prd <- predict(probitfit0)
          #starting.values$ystar <- ifelse(y == 1, abs(prd), -abs(prd))
          starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
        }
      } else {
        starting.values$beta <- 0
        starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
      }
    }
  }

  ##print (starting.values)
  ##truncate vectors that are too long
  if (length(starting.values$h.hat) > length(y)) {
    starting.values$h.hat <- starting.values$h.hat[1:length(y)]
  }
  if (length(starting.values$beta) > ncol(X)) {
    starting.values$beta <- starting.values$beta[1:ncol(X)]
  }
  if (length(starting.values$delta) > ncol(Z)) {
    starting.values$delta <- starting.values$delta[1:ncol(Z)]
  }
  if (varsel==FALSE & rmethod == "equal" & length(starting.values$r) > 1) {
    starting.values$r <- starting.values$r[1] ## this should only happen if rmethod == "equal"
  } else if (length(starting.values$r) > ncol(Z)) {
    starting.values$r <- starting.values$r[1:ncol(Z)]
  }

  chain$h.hat[1, ] <- starting.values$h.hat
  chain$beta[1, ] <- starting.values$beta
  chain$lambda[1, ] <- starting.values$lambda
  chain$sigsq.eps[1] <- starting.values$sigsq.eps
  chain$r[1, ] <- starting.values$r
  if (varsel) {
    chain$delta[1,ztest] <- starting.values$delta
  }
  if (family == "binomial") {
    chain$ystar[1, ] <- starting.values$ystar
    chain$sigsq.eps[] <- starting.values$sigsq.eps ## does not get updated
  }
  if (!is.null(groups)) {
    ## make sure starting values are consistent with structure of model
    if (!all(sapply(unique(groups), function(x) sum(chain$delta[1, ztest][groups == x])) == 1)) {
      # warning("Specified starting values for delta not consistent with model; using default")
      starting.values$delta <- rep(0, length(groups))
      starting.values$delta[sapply(unique(groups), function(x) min(which(groups == x)))] <- 1
    }
    chain$delta[1,ztest] <- starting.values$delta
    chain$r[1,ztest] <- ifelse(chain$delta[1,ztest] == 1, chain$r[1,ztest], 0)
  }
  chain$est.h <- est.h

  ## components


  Vcomps <- a_makeVcomps(X = X, r = chain$r[1, ], lambda = chain$lambda[1, ], Z = Z, data.comps = data.comps)






  ## start sampling ####
  #chain$time1 <- Sys.time()
  s=2
  for (s in 2:nsamp) {

    ## continuous version of outcome (latent outcome under binomial probit model)
    if (family == "gaussian") {
      ycont <- y
    } else if (family == "binomial") {
      if (est.h) {
        chain$ystar[s,] <- ystar.update(y = y, X = X, beta = chain$beta[s - 1,], h = chain$h[s - 1, ])
      } else {
        chain$ystar[s,] <- ystar.update.noh(y = y, X = X, beta = chain$beta[s - 1,], Vinv = Vcomps$Vinv, ystar = chain$ystar[s - 1, ])
      }
      ycont <- chain$ystar[s, ]
    }

    ## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)

    ## beta
    if (!missingX) {
      chain$beta[s,] <- a_beta.update(X = X, XVinv = Vcomps$XVinv, y = ycont, sigsq.eps = chain$sigsq.eps[s - 1])
    }



    ## \sigma_\epsilon^2
    if (family == "gaussian") {
      chain$sigsq.eps[s] <- a_sigsq.eps.update(y = ycont, X = X, beta = chain$beta[s,], Vcomps = Vcomps, a.eps = control.params$a.sigsq, b.eps = control.params$b.sigsq)
    }

    ## lambda
    lambdaSim <- chain$lambda[s - 1,]
    for (comp in 1:data.comps$nlambda) {
      varcomps <- a_lambda.update(r = chain$r[s - 1,], delta = chain$delta[s - 1,], lambda = lambdaSim, whichcomp = comp, y = ycont, X = X, Z = Z, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, data.comps = data.comps, control.params = control.params)
      lambdaSim <- varcomps$lambda
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.lambda[s,comp] <- varcomps$acc
      }
    }
    chain$lambda[s,] <- lambdaSim

    ## r
    rSim <- chain$r[s - 1,]
    comp <- which(!1:ncol(Z) %in% ztest)
    if (length(comp) != 0) {
      if (rmethod == "equal") { ## common r for those variables not being selected
        varcomps <- a_r.update(r = rSim, whichcomp = comp, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
        rSim <- varcomps$r
        if (varcomps$acc) {
          Vcomps <- varcomps$Vcomps
          chain$acc.r[s, comp] <- varcomps$acc
        }
      } else if (rmethod == "varying") { ## allow a different r_m
        for (whichr in comp) {
          varcomps <- a_r.update(r = rSim, whichcomp = whichr, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
          rSim <- varcomps$r
          if (varcomps$acc) {
            Vcomps <- varcomps$Vcomps
            chain$acc.r[s, whichr] <- varcomps$acc
          }
        }
      }
    }
    ## for those variables being selected: joint posterior of (r,delta)
    if (varsel) {
      varcomps <- rdelta.update(r = rSim, delta = chain$delta[s - 1,], lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z, ztest = ztest, data.comps = data.comps, control.params = control.params, rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1, rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2, rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen, rprop.logdens = rprop.logdens)
      chain$delta[s,] <- varcomps$delta
      rSim <- varcomps$r
      chain$move.type[s] <- varcomps$move.type
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.rdelta[s] <- varcomps$acc
      }
    }
    chain$r[s,] <- rSim

    ###################################################
    ## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)

    if (est.h) {
      hcomps <- h.update(lambda = chain$lambda[s,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s], y = ycont, X = X, beta = chain$beta[s,], r = chain$r[s,], Z = Z, data.comps = data.comps)
      chain$h.hat[s,] <- hcomps$hsamp
      if (!is.null(hcomps$hsamp.star)) { ## GPP
        Vcomps$hsamp.star <- hcomps$hsamp.star
      }
      rm(hcomps)
    }

    ###################################################
    ## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)

    if (!is.null(Znew)) {
      chain$hnew[s,] <- newh.update(Z = Z, Znew = Znew, Vcomps = Vcomps, lambda = chain$lambda[s,], sigsq.eps = chain$sigsq.eps[s], r = chain$r[s,], y = ycont, X = X, beta = chain$beta[s,], data.comps = data.comps)
    }

    ###################################################
    ## print details of the model fit so far
    # opts <- bkmr:::set_verbose_opts(
    #   verbose_freq = control.params$verbose_freq,
    #   verbose_digits = control.params$verbose_digits,
    #   verbose_show_ests = control.params$verbose_show_ests
    # )
    # bkmr:::print_diagnostics(verbose = verbose, opts = opts, curr_iter = s, tot_iter = nsamp, chain = chain, varsel = varsel, hier_varsel = hier_varsel, ztest = ztest, Z = Z, groups = groups)
    #


  }
  control.params$r.params <- NULL
  chain$time2 <- Sys.time()
  chain$iter <- nsamp
  chain$family <- family
  chain$starting.values <- starting.values
  chain$control.params <- control.params
  chain$X <- X
  chain$Z <- Z
  chain$y <- y
  chain$ztest <- ztest
  chain$data.comps <- data.comps
  if (!is.null(Znew)) chain$Znew <- Znew
  if (!is.null(groups)) chain$groups <- groups
  chain$varsel <- varsel
  class(chain) <- c("bkmrfit", class(chain))
  chain
}


# output ------------------------------------------------------------------


#' @title Generate overall exposure-risk curves
#' @description Estimate joint exposure effects across a set of counterfactual scenarios (newz).
#'
#' @param fit      A `bkmrfit` object returned by `kmbayes()`
#' @param y        Outcome vector of length n
#' @param Z        n×M exposure matrix
#' @param X        n×K covariate matrix (no intercept column)
#' @param newz     Data frame of counterfactual exposure scenarios (rows are scenarios, columns are exposures)
#' @param point1   Reference scenario, defaults to column-wise medians of Z
#' @param method   One of "approx" (faster) or "exact" for posterior summaries
#' @param sel      Numeric vector indicating which MCMC iterations to use
#' @param data.comps A list of components for the approximate Gaussian process, generated by `a_kmbayes`
#'
#' @return A list containing:
#'   \itemize{
#'     \item preds_e: list with `postmean` and `postvar` of h(newz)
#'     \item preds_p: list with `postmean` and `postvar` of differences from reference
#'     \item risk_overall: data frame with newz, estimated effects and standard errors
#'   }
#' @export

a_OverallRiskSummaries_vary <- function(
    fit,
    y = NULL,
    Z = NULL,
    X = NULL,
    newz = NULL,
    method = "approx",
    sel = NULL,
    data.comps,
    point1 = NULL
) {
  # If fit is a bkmrfit object, extract y, Z, and X if not provided
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))   y <- fit$y
    if (is.null(Z))   Z <- fit$Z
    if (is.null(X))   X <- fit$X
  }

  # Generate default grid of quantiles from 25% to 75% in 1% increments
  if (is.null(newz)) {
    quants <- seq(0.25, 0.75, 0.01)
    newz <- lapply(seq_len(ncol(Z)), function(i) {
      df <- data.frame(quantile(Z[, i], quants, na.rm = TRUE))
      names(df) <- colnames(Z)[i]
      df
    }) %>% bind_cols()
  }

  # Default reference point is the median of each exposure
  if (is.null(point1)) {
    point1 <- data.frame(t(colMedians(Z, na.rm = TRUE)))
  }
  colnames(newz) <- colnames(Z)

  # Find if reference scenario is already in newz
  hits <- which(apply(newz, 1, function(row) all(row == point1)))

  # If not, prepend reference scenario to the grid
  if (length(hits) == 0) {
    newz <- rbind(point1, newz)
    which_ref <- 1L
  } else {
    which_ref <- hits[1]
  }

  # Compute posterior mean and variance for h(newz)
  preds <- aBKMR:::a_ComputePostmeanHnew(
    fit = fit, y = y, Z = Z, X = X,
    Znew = newz, sel = sel,
    method = method, data.comps = data.comps
  )

  n <- nrow(newz)
  # Create contrast matrix to compute differences from reference
  C <- diag(n)
  C[, which_ref] <- C[, which_ref] - 1

  # Compute differences and their variance
  diff_vals <- as.vector(C %*% preds$postmean)
  Sigma_diff <- C %*% preds$postvar %*% t(C)

  # Prepare output risk_overall without the reference row
  if (which_ref == 1) {
    risk <- newz[-1, , drop = FALSE]
    risk$est_p <- diff_vals[-1]
    risk$sd_p  <- sqrt(diag(Sigma_diff))[-1]
    risk$est_e <- preds$postmean[-1]
    risk$sd_e  <- sqrt(diag(preds$postvar))[-1]
  } else {
    risk <- newz
    risk$est_p <- diff_vals
    risk$sd_p  <- sqrt(diag(Sigma_diff))
    risk$est_e <- preds$postmean
    risk$sd_e  <- sqrt(diag(preds$postvar))
  }

  # Return list of posterior summaries and summary data frame
  list(
    preds_e = list(postmean = preds$postmean[-which_ref],
                   postvar  = preds$postvar[-which_ref, -which_ref]),
    preds_p = list(postmean = diff_vals[-which_ref],
                   postvar  = Sigma_diff[-which_ref, -which_ref]),
    risk_overall = risk
  )
}


#' used for multivariate interaction plot
#' @param fit An object containing the results returned by a the kmbayes function
#' @param y a vector of outcome data of length n
#' @param Z an n-by-M matrix of predictor variables to be included in the h function. Each row represents an observation and each column represents an predictor
#' @param X an n-by-K matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column
#' @param data.comps a list including objects used for approximate gaussian process, which is generated by function "a_kmbayes"
#' @param which.z vector indicating which variables (columns of Z) for which the summary should be computed
#' @param qs.diff vector indicating the two quantiles at which to compute the single-predictor risk summary
#' @param q.fixed vector indicating the two quantiles at which to fix all of the remaining exposures in Z
#' @param method method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets; see details
#' @param sel logical expression indicating samples to keep; defaults to keeping the second half of all samples
#' @param z.names optional vector of names for the columns of z
#' @param ... other arguments to pass on to the prediction function
#' @export


a_SingVarRiskSummaries = function (fit, y = NULL, Z = NULL, X = NULL, data.comps, which.z = 1:ncol(Z),
                                   qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.5, 0.75), method = "approx",
                                   sel = NULL, z.names = colnames(Z), ...) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }
  if (is.null(z.names))
    z.names <- paste0("z", 1:ncol(Z))

  para = expand.grid(qf = seq_along(q.fixed), wz = seq_along(which.z))

  df = bind_rows(lapply(1:nrow(para), function(i){
    whichz_temp = which.z[para$wz[i]]
    q_temp = q.fixed[para$qf[i]]

    point2 <- point1 <- apply(Z, 2, quantile, q_temp)
    point2[whichz_temp] <- apply(Z[, whichz_temp, drop = FALSE], 2, quantile,
                                 qs.diff[2])
    point1[whichz_temp] <- apply(Z[, whichz_temp, drop = FALSE], 2, quantile,
                                 qs.diff[1])

    cc <- c(-1, 1)
    newz <- rbind(point1, point2)
    preds <- aBKMR:::a_ComputePostmeanHnew(fit = fit,
                                   y = y, Z = Z, X = X, data.comps = data.comps, Znew = newz, sel = sel, method = method)
    diff <- drop(cc %*% preds$postmean)
    diff.sd <- drop(sqrt(cc %*% preds$postvar %*% cc))
    return( data.frame(q.fixed = q.fixed[para$qf[i]], variable = z.names[para$wz[i]], est_p = diff, sd_p = diff.sd, est_e = preds$postmean[2], sd_e = preds$postvar[4]^0.5))
  }))

  df <- dplyr::mutate_at(df, "variable", function(x) factor(x,
                                                            levels = z.names[which.z]))
  df <- dplyr::mutate_at(df, "q.fixed", function(x) as.factor(x))
  attr(df, "qs.diff") <- qs.diff
  df
}



#' used for univariate effect plot
#' @param fit An object containing the results returned by a the kmbayes function
#' @param y a vector of outcome data of length n
#' @param Z an n-by-M matrix of predictor variables to be included in the h function. Each row represents an observation and each column represents an predictor
#' @param X an n-by-K matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column
#' @param data.comps a list including objects used for approximate gaussian process, which is generated by function "a_kmbayes"
#' @param quants specify the quantiles of exposures
#' @param method method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets; see details
#' @param sel logical expression indicating samples to keep; defaults to keeping the second half of all samples
#' @export


a_PredictorResponseUnivar = function (fit, y = NULL, Z = NULL, X = NULL, quants = seq(0.25, 0.75, 0.01), method = "approx", sel = NULL,
                                   data.comps, point1 = NULL) {

  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }

  if(is.null(point1)){
    point1 = apply(Z, 2, median, na.rm = TRUE)
  }

  newz_temp = lapply(1:length(quants), function(i) point1) %>% bind_rows()


  s_res = lapply(1:length_z, function(i){
    newz = newz_temp
    newz[,i] = quantile(Z[, i], quants, na.rm = T)
    res = a_OverallRiskSummaries_vary(fit, y = NULL, Z = NULL, X = NULL, newz = newz , method = "approx", sel = NULL, data.comps, point1 = point1)
    res$est = com_je(postmean = res$preds_p$postmean, postvar = res$preds_p$postvar, X = quants*10, B = 1000)
    res$risk_overall$quants = quants

    res$risk_overall = res$risk_overall[, c("quants", glue("z{i}"), "est_p", "sd_p")]
    names(res$risk_overall) = c("quants", "z", "est", "se")
    res$risk_overall$vars = names(data.frame(Z))[i]

    return(res)
  })

}




#' used for bivariate interaction plot
#' @param fit An object containing the results returned by a the kmbayes function
#' @param y a vector of outcome data of length n
#' @param Z an n-by-M matrix of predictor variables to be included in the h function. Each row represents an observation and each column represents an predictor
#' @param X an n-by-K matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column
#' @param data.comps a list including objects used for approximate gaussian process, which is generated by function "a_kmbayes"
#' @param z.pairs grid dataset of z
#' @param z.names names of z
#' @param verbose whether to print details of running
#' @param method method for obtaining posterior summaries at a vector of new points. Options are "approx" and "exact"; defaults to "approx", which is faster particularly for large datasets; see details
#' @param q.fixed vector of quantiles at which to fix the remaining predictors in Z
#' @param sel logical expression indicating samples to keep; defaults to keeping the second half of all samples
#' @param ngrid number of grid points to cover the range of each predictor (column in Z)
#' @param center flag for whether to scale the exposure-response function to have mean zero
#' @param min.plot.dist specifies a minimum distance that a new grid point needs to be from an observed data point in order to compute the prediction; points further than this will not be computed
#' @param center flag for whether to scale the exposure-response function to have mean zero
#' @param ... other arguments to pass on to the prediction function
#' @export

a_PredictorResponseBivar = function (fit, y = NULL, Z = NULL, X = NULL, data.comps = data.comps, z.pairs = NULL,
                                      method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL,
                                      min.plot.dist = 0.5, center = TRUE, z.names = colnames(Z),
                                      verbose = TRUE, ...) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y))
      y <- fit$y
    if (is.null(Z))
      Z <- fit$Z
    if (is.null(X))
      X <- fit$X
  }
  if (is.null(z.names)) {
    z.names <- colnames(Z)
    if (is.null(z.names)) {
      z.names <- paste0("z", 1:ncol(Z))
    }
  }
  if (is.null(z.pairs)) {
    z.pairs <- expand.grid(z1 = 1:ncol(Z), z2 = 1:ncol(Z))
    z.pairs <- z.pairs[z.pairs$z1 < z.pairs$z2, ]
  }
  df <- dplyr::tibble()
  for (i in 1:nrow(z.pairs)) {
    compute <- TRUE
    whichz1 <- z.pairs[i, 1] %>% unlist %>% unname
    whichz2 <- z.pairs[i, 2] %>% unlist %>% unname
    if (whichz1 == whichz2)
      compute <- FALSE
    z.name1 <- z.names[whichz1]
    z.name2 <- z.names[whichz2]
    names.pair <- c(z.name1, z.name2)
    if (nrow(df) > 0) {
      completed.pairs <- df %>% dplyr::select_at(c("variable1",
                                                   "variable2")) %>% dplyr::distinct() %>% dplyr::transmute(z.pair = paste("variable1",
                                                                                                                           "variable2", sep = ":")) %>% unlist %>% unname
      if (paste(names.pair, collapse = ":") %in% completed.pairs |
          paste(rev(names.pair), collapse = ":") %in% completed.pairs)
        compute <- FALSE
    }
    if (compute) {
      if (verbose)
        message("Pair ", i, " out of ", nrow(z.pairs))
      res <- a_PredictorResponseBivarPair(fit = fit, y = y,
                                           Z = Z, X = X, data.comps = data.comps, whichz1 = whichz1, whichz2 = whichz2,
                                           method = method, ngrid = ngrid, q.fixed = q.fixed,
                                           sel = sel, min.plot.dist = min.plot.dist, center = center,
                                           z.names = z.names, ...)
      df0 <- res
      df0$variable1 <- z.name1
      df0$variable2 <- z.name2
      df0 %<>% dplyr::select_at(c("variable1", "variable2",
                                  "z1", "z2", "est", "se"))
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$variable1 <- factor(df$variable1, levels = z.names)
  df$variable2 <- factor(df$variable2, levels = z.names)
  df
}







# quantitative estimation -------------------------------------------------

#' used for joint effect estimation
#' @param postmean prrdicted values of y
#' @param postvar variance-covariance matrix of predicted y
#' @param B the number of bootstrap. By default, it is 1000
#' @export

com_je = function(postmean, postvar, X, B=1000){

  coef_df = lapply(1:B, function(i){

    y_star <- mvrnorm(n = 1, mu = postmean, Sigma = postvar)

    fit_b   <- lm(y_star ~ X)

    return(coef(fit_b))


  }) %>% bind_rows()

  return(data.frame(est = coef(lm(postmean~X))[2],
                    lower_quant = quantile(coef_df$X, 0.025),
                    upper_quant = quantile(coef_df$X, 0.975),
                    sd = sd(coef_df$X),
                    lower_norm = coef(lm(postmean~X))[2] - 1.96*sd(coef_df$X),
                    upper_norm = coef(lm(postmean~X))[2] + 1.96*sd(coef_df$X)

  ))

}



#' used for bivariate interaction estimation
#' @param var1 first variable
#' @param var2 second variable
#' @param bivar_obj the object obtained by a_PredictorResponseBivar
#' @export

com_2interaction = function(var1, var2, bivar_obj){

  z2_data = bivar_obj %>%
    filter(variable1 %in% c(var1, var2) & variable2 %in% c(var1, var2))



  z2_var1 = z2_data %>% filter(variable1 == var1)

  z2_var2 = z2_data %>% filter(variable1 == var2)

  z2_var1_dev = z2_var1 %>%
    group_by(quantile) %>%
    arrange(z1) %>%
    mutate(
      der_est = c(diff(est) / diff(z1), NA)
    )


  z2_var2_dev = z2_var2 %>%
    group_by(quantile) %>%
    arrange(z1) %>%
    mutate(
      der_est = c(diff(est) / diff(z1), NA)
    )

  var1_modified_by_var2 = compareGrowthCurves(as.character(z2_var1_dev$quantile), as.matrix(z2_var1_dev$der_est) , levels=NULL, nsim=100, fun=meanT, times=NULL,
                                              verbose=TRUE, adjust="holm", n0=0.5)


  var2_modified_by_var1 = compareGrowthCurves(as.character(z2_var2_dev$quantile), as.matrix(z2_var2_dev$der_est) , levels=NULL, nsim=100, fun=meanT, times=NULL,
                                              verbose=TRUE, adjust="holm", n0=0.5)

  return(data.frame(interaction_type = c("var1_modified_by_var2", "var2_modified_by_var1"),
                    stat = c(var1_modified_by_var2$Stat, var2_modified_by_var1$Stat),
                    P = c(var1_modified_by_var2$P.Value, var2_modified_by_var1$P.Value)    ))


}


#' used for multivariate interaction estimation
#' @param beta coefficient of a variable
#' @param se se of beta
#' @export

trend <- function(beta, se) {

  # beta must be a vector of estimates
  # se must be a vector of their standard errors


  z <- seq(1:length(beta))
  w <- 1 / (se^2)
  test <- (sum(w * beta *(z - (sum(z * w) / sum(w))))^2) / sum(w * ((z - sum(z * w) / sum(w))^2))
  p <- pchisq(test, df = 1, lower.tail = FALSE)
  return(list("test" = test, "p_trend" = p))
}
