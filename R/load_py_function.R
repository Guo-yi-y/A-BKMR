.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the 'reticulate' package first.")
  }

  # 只在交互式会话里尝试安装 miniconda
  if (interactive()) {
    mc <- tryCatch(reticulate::miniconda_path(), error = function(e) "")
    if (!nzchar(mc) || !dir.exists(mc)) {
      packageStartupMessage("Miniconda not found; installing with reticulate::install_miniconda()...")
      tryCatch(reticulate::install_miniconda(),
               error = function(e) {
                 packageStartupMessage("Miniconda 安装失败，请手动安装：https://docs.conda.io/en/latest/miniconda.html")
               })
    }
    reticulate::use_miniconda(reticulate::miniconda_path(), required = TRUE)
  }

  # 然后只做基础的模块可用性检查，不再调用 py_install()
  missing <- Filter(function(pkg) !reticulate::py_module_available(pkg),
                    c("numpy","pandas","scipy","mpmath"))
  if (length(missing)) {
    packageStartupMessage(
      "缺少 Python 包：", paste(missing, collapse = ", "),
      "；请在 R 里运行 reticulate::py_install(missing) 或者自行在 Conda 环境中安装。"
    )
  }

  # 最后再载入你的 Python 脚本
  python_file <- system.file("python", "py_functions.py", package = pkgname)
  if (python_file == "") stop("py_functions.py 找不到")
  reticulate::source_python(python_file)
}
