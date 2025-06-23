# R/load_py_function.R

# Automatically check Python environment, required packages, and source Python script on package load
.onLoad <- function(libname, pkgname) {


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
