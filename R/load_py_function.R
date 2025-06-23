# R/load_py_function.R

.onLoad <- function(libname, pkgname) {


  # 4. 加载 Python 文件
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") {
    stop("Python file 'py_functions.py' not found.")
  }

  # 加载 Python 文件
  reticulate::source_python(python_file)
}

