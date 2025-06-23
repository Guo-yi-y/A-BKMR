# R/load_py_function.R

# Automatically check Python environment, required packages, and source Python script on package load
.onLoad <- function(libname, pkgname) {
  # 1. Check for 'reticulate' package
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    packageStartupMessage(
      "⚠️ 'reticulate' package not found. Please install it via: install.packages('reticulate')"
    )
    return()
  }

  # 2. Verify Python availability
  py_ok <- reticulate::py_available()
  if (!py_ok) {
    packageStartupMessage(
      "⚠️ Python environment not detected. To install Miniconda, run in R console:\n",
      "    install.packages('reticulate'); reticulate::install_miniconda()"
    )
    return()
  }

  # 3. Print Python configuration and check required Python packages
  cfg <- reticulate::py_config()
  packageStartupMessage(
    sprintf("✅ Python detected at %s (version %s, %d-bit)",
            cfg$python, cfg$version, cfg$bitarchitecture)
  )
  required_pkgs <- c("numpy", "mpmath", "pandas", "scipy")
  for (pkg in required_pkgs) {
    if (reticulate::py_module_available(pkg)) {
      packageStartupMessage(sprintf("  ✅ Python package '%s' is installed.", pkg))
    } else {
      packageStartupMessage(
        sprintf(
          "  ❌ Python package '%s' not found; you can install it in R via: reticulate::py_install('%s')",
          pkg, pkg
        )
      )
    }
  }

  # 4. Source the Python script from inst/python
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") {
    stop("❌ Python script 'py_functions.py' not found. Please ensure the package is installed correctly.")
  }
  reticulate::source_python(python_file)
}




utils::globalVariables(c(
  "ComputePostmeanHnew.exact", "K", "Vinv", "est",
  "variable", "variable1", "variable2", "z1", "."
))
