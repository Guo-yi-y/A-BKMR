# R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 1. 确保加载 reticulate 包
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the reticulate package first.")
  }
  library(reticulate)
  # 2. 检查 Python 是否已安装
  python_installed <- tryCatch({
    reticulate::py_config()
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!python_installed) {
    message("Python is not installed or not detected.")

    # 尝试通过 Miniconda 安装 Python
    message("Attempting to install Python using reticulate::install_python()...")
    tryCatch({
      reticulate::install_python()
      message("Python has been successfully installed.")
    }, error = function(e) {
      message("Failed to install Python. Please install Python manually.")
      message("You can install Python via: https://www.python.org/downloads/")
    })
  } else {
    message("Python is installed. Proceeding with the package load...")
  }

  # 3. 声明 Python 依赖（只执行一次，后续会话重用同一环境）
  reticulate::py_require(c("numpy", "mpmath", "pandas", "scipy"))

  # 4. 延迟加载（delay_load）Python 模块，防止在 onLoad 时马上初始化 Python
  numpy  <<- reticulate::import("numpy",  delay_load = TRUE)
  mpmath <<- reticulate::import("mpmath", delay_load = TRUE)
  pandas <<- reticulate::import("pandas", delay_load = TRUE)
  scipy  <<- reticulate::import("scipy",  delay_load = TRUE)

  # 5. 加载 Python 文件
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") {
    stop("Python file 'py_functions.py' not found.")
  }

  # 加载 Python 文件
  reticulate::source_python(python_file)
}
