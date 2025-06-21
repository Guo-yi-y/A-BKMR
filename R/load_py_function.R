## R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 1. Ensure reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the 'reticulate' package first.")
  }

  # 2. Install Miniconda if missing
  if (!reticulate::miniconda_exists()) {
    packageStartupMessage("Miniconda not found. Installing Miniconda...")
    reticulate::install_miniconda()
  }

  # 3. Create (once) and activate a persistent conda environment
  envname <- "aBKMR-env"
  conda_info <- reticulate::conda_list()
  if (!(envname %in% conda_info$name)) {
    packageStartupMessage(sprintf("Creating conda environment '%s' and installing Python dependencies...", envname))
    reticulate::conda_create(envname)
    reticulate::conda_install(envname,
                              packages = c("numpy", "pandas", "scipy", "mpmath"),
                              channel = "conda-forge")
  }
  # Always use the persistent env
  reticulate::use_condaenv(envname, required = TRUE)
  packageStartupMessage(sprintf("Using conda environment '%s'", envname))

  # 4. Load the Python functions file
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") {
    stop("Python file 'py_functions.py' not found in inst/python/ folder.")
  }
  reticulate::source_python(python_file)
}

# Register global variables
utils::globalVariables(
  c("ComputePostmeanHnew.exact", "K", "Vinv", "est",
    "variable", "variable1", "variable2", "z1", ".")
)
