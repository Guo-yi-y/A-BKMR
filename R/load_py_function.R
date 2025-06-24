# R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 1. 确保加载 reticulate 包
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the reticulate package first.")
  }
  library(reticulate)


  # 2. 声明 Python 依赖（只执行一次，后续会话重用同一环境）
  reticulate::py_require(
    packages = c("numpy", "pandas", "scipy")
  )



  # 3. 加载 Python 文件
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") {
    stop("Python file 'py_functions.py' not found.")
  }

  # 加载 Python 文件
  reticulate::source_python(python_file)
}
