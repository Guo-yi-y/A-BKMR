# R/load_py_function.R

#' @import reticulate
.onLoad <- function(libname, pkgname) {
  # Ensure reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the 'reticulate' package first.")
  }

  # Only perform environment setup if Python hasn't been initialized yet
  if (!reticulate::py_available(initialize = FALSE)) {
    envname <- "aBKMR-env"

    # Install Miniconda if not present
    mc_path <- tryCatch(reticulate::miniconda_path(), error = function(e) "")
    if (!nzchar(mc_path) || !dir.exists(mc_path)) {
      packageStartupMessage("Miniconda not detected. Installing via reticulate::install_miniconda()...")
      tryCatch(
        reticulate::install_miniconda(),
        error = function(e) {
          packageStartupMessage("Miniconda installation failed. Please install it manually.")
        }
      )
      mc_path <- reticulate::miniconda_path()
    }

    # Create conda environment if missing
    envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))
    if (!(envname %in% envs)) {
      packageStartupMessage(sprintf("Creating conda environment '%s' with required Python packages...", envname))
      reticulate::conda_create(envname, packages = c("numpy", "pandas", "scipy", "mpmath"))
    }

    # Activate the persistent environment
    reticulate::use_condaenv(envname, required = TRUE)
  } else {
    packageStartupMessage("Python has already been initialized; skipping environment setup.")
  }

  # Source the Python functions file
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") {
    stop("Python file 'py_functions.py' not found in the installed package.")
  }
  reticulate::source_python(python_file)
}

# Avoid R CMD check notes for imported Python-side variables
utils::globalVariables(c(
  "ComputePostmeanHnew.exact", "K", "Vinv", "est",
  "variable", "variable1", "variable2", "z1", "."
))
