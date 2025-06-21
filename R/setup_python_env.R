#' @title Initialize Python environment for aBKMR
#' @description Run this once in an interactive R session to create
#'   a conda env and install required Python packages.
#' @export
setup_python_env <- function(envname = "aBKMR-env") {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install reticulate first.")
  }
  if (!nzchar(reticulate::miniconda_path())) {
    reticulate::install_miniconda()
  }
  reticulate::conda_create(envname,
                           packages = c("numpy", "pandas", "scipy", "mpmath")
  )
  message("Conda environment '", envname,
          "' created with numpy,pandas,scipy,mpmath.")
}
