# R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 0. 指定使用预先创建的 conda 环境
  reticulate::use_condaenv("r-reticulate", required = TRUE)

  # 1. 确保 reticulate 已加载
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Please install the 'reticulate' package first.")

  # 2. 检查 Python 可用性
  if (!inherits(try(reticulate::py_config(), silent = TRUE), "py_config"))
    packageStartupMessage("Warning: Could not detect Python. Continuing...")

  # 3. 确保依赖包已“require”
  for (pkg in c("numpy","mpmath","pandas","scipy")) {
    reticulate::py_require(pkg)
  }

  # 4. 加载 Python 代码
  python_file <- system.file("python/py_functions.py", package = pkgname)
  if (python_file == "") stop("Python file not found.")
  reticulate::source_python(python_file)
}


utils::globalVariables(c(
  "ComputePostmeanHnew.exact", "K", "Vinv", "est",
  "variable", "variable1", "variable2", "z1", "."
))
