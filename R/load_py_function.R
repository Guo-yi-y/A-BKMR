.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the 'reticulate' package first.")
  }

  # 只有在交互式会话且 Python 没有被初始化时，才去激活我们的 conda 环境
  if (interactive() && !reticulate::py_available(initialize = FALSE)) {
    env <- "aBKMR-env"
    # 如果用户自己创建了这个环境，就激活它
    if (env %in% reticulate::conda_list()$name) {
      reticulate::use_condaenv(env, required = TRUE)
      packageStartupMessage("Using conda env '", env, "'")
    } else if (nzchar(reticulate::miniconda_path())) {
      # 否则如果存在默认 miniconda，就激活它
      reticulate::use_miniconda(reticulate::miniconda_path(), required = TRUE)
      packageStartupMessage("Using Miniconda at ", reticulate::miniconda_path())
    }
  }
}
