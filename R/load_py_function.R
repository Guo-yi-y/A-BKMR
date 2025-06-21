# R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 1. 确保加载 reticulate 包
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the 'reticulate' package first.")
  }

  # 2. 检查 Python 是否已安装
  python_installed <- tryCatch({
    reticulate::py_config()
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!python_installed) {
    packageStartupMessage("Python is not installed or not detected.")
    packageStartupMessage("Attempting to install Python using reticulate::install_python()...")
    tryCatch({
      reticulate::install_python()
      packageStartupMessage("Python has been successfully installed.")
    }, error = function(e) {
      packageStartupMessage("Failed to install Python. Please install Python manually.")
      packageStartupMessage("You can install Python via: https://www.python.org/downloads/")
    })
  } else {
    packageStartupMessage("Python is installed. Proceeding with the package load...")
  }

  # 3. 检查并安装所需的 Python 包
  required_packages <- c("numpy", "mpmath", "pandas", "scipy")

  for (pkg in required_packages) {
    installed <- tryCatch({
      reticulate::py_run_string(paste("import", pkg))
      TRUE
    }, error = function(e) {
      FALSE
    })

    if (!installed) {
      packageStartupMessage(paste(pkg, "not found. Installing..."))
      tryCatch({
        reticulate::py_install(pkg)
        packageStartupMessage(paste(pkg, "has been successfully installed."))
      }, error = function(e) {
        packageStartupMessage(paste("Failed to install", pkg, ". Please install it manually."))
      })
    } else {
      packageStartupMessage(paste(pkg, "is already installed."))
    }
  }

  # 4. 加载 Python 文件
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") {
    stop("Python file 'py_functions.py' not found.")
  }

  reticulate::source_python(python_file)
}

utils::globalVariables(c(
  "ComputePostmeanHnew.exact", "K", "Vinv", "est",
  "variable", "variable1", "variable2", "z1", "."
))
