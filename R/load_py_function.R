## R/load_py_function.R
.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the reticulate package first.")
  }

  # 检查是否已安装 Python
  python_installed <- tryCatch({
    reticulate::py_config()
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!python_installed) {
    message("Python is not installed or not detected.")

    # 自动安装 Python（通过 Miniconda）
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

  # 检查并安装所需的 Python 包
  required_packages <- c("numpy", "mpmath", "pandas")

  for (pkg in required_packages) {
    # 检查 Python 包是否已安装
    installed <- tryCatch({
      reticulate::py_run_string(paste("import", pkg))
      TRUE
    }, error = function(e) {
      FALSE
    })

    if (!installed) {
      message(paste(pkg, "not found. Installing...", sep = " "))

      # 安装缺失的 Python 包
      tryCatch({
        reticulate::py_install(pkg)
        message(paste(pkg, "has been successfully installed.", sep = " "))
      }, error = function(e) {
        message(paste("Failed to install", pkg, ". Please install it manually.", sep = " "))
      })
    } else {
      message(paste(pkg, "is already installed.", sep = " "))
    }
  }
}
