# R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 1. 确保 reticulate 可用
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    packageStartupMessage(
      "⚠️ 需要安装 reticulate：\n",
      "    install.packages('reticulate')"
    )
    return()
  }

  # 2. 检查 Python 环境
  python_ok <- tryCatch({
    reticulate::py_config()
    TRUE
  }, error = function(e) FALSE)

  if (!python_ok) {
    packageStartupMessage(
      "⚠️ 未检测到 Python 环境。\n",
      "  您可以在 R 中运行：\n",
      "    reticulate::install_miniconda()"
    )
    return()
  }

  cfg <- reticulate::py_config()
  packageStartupMessage(
    sprintf("✅ 检测到 Python：%s （版本 %s，%d 位）",
            cfg$python, cfg$version, cfg$bitarchitecture)
  )

  # 3. 检查必需的 Python 包
  required_py_pkgs <- c("numpy", "mpmath", "pandas", "scipy")
  missing <- vapply(required_py_pkgs,
                    function(x) !reticulate::py_module_available(x),
                    logical(1))

  if (any(missing)) {
    pkg_list <- paste(required_py_pkgs[missing], collapse = ", ")
    packageStartupMessage(
      sprintf("⚠️ 缺失 Python 包：%s\n", pkg_list),
      "  您可以在 R 中运行：\n",
      sprintf("    reticulate::py_install(c(%s))",
              paste(sprintf("'%s'", required_py_pkgs[missing]), collapse = ", "))
    )
  }
}

.onAttach <- function(libname, pkgname) {
  # 仅当 Python 环境和 reticulate 都正常时，才去 source
  if (requireNamespace("reticulate", quietly = TRUE) &&
      !inherits(try(reticulate::py_config(), silent = TRUE), "try-error")) {

    python_file <- system.file("python", "py_functions.py", package = pkgname)
    if (python_file == "") {
      packageStartupMessage(
        "❌ 未找到 py_functions.py，跳过加载。请确认它已放在 inst/python/ 目录下并重新安装包。"
      )
    } else {
      reticulate::source_python(python_file)
      packageStartupMessage("✅ 已加载 Python 函数：py_functions.py")
    }
  }
}
