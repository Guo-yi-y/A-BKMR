# R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 1. 检查 reticulate
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    packageStartupMessage(
      "⚠️ 需要安装 reticulate：\n",
      "    install.packages('reticulate')"
    )
    return()
  }

  # 2. 探测 Python 环境
  py_ok <- !inherits(try(reticulate::py_config(), silent = TRUE), "try-error")
  if (!py_ok) {
    packageStartupMessage(
      "⚠️ 未检测到 Python 环境。\n",
      "  请在 R 中运行：\n",
      "    reticulate::install_miniconda() 或 配置好你的系统 Python"
    )
    return()
  }

  cfg <- reticulate::py_config()
  packageStartupMessage(
    sprintf("✅ 检测到 Python：%s （版本 %s，%d 位）",
            cfg$python, cfg$version, cfg$bitarchitecture)
  )

  # 3. 探测必要的 Python 包
  required <- c("numpy", "mpmath", "pandas", "scipy")
  missing <- required[!vapply(required, reticulate::py_module_available, logical(1))]
  if (length(missing)) {
    packageStartupMessage(
      sprintf("⚠️ 缺失 Python 模块：%s\n", paste(missing, collapse = ", ")),
      "  可在 R 中运行：\n",
      sprintf("    reticulate::py_install(c(%s))",
              paste(sprintf("'%s'", missing), collapse = ", "))
    )
  }
}

.onAttach <- function(libname, pkgname) {
  # 只有在 Python 环境和 reticulate 都正常时，才去导入脚本
  if (requireNamespace("reticulate", quietly = TRUE) &&
      !inherits(try(reticulate::py_config(), silent = TRUE), "try-error")) {

    # （可选）再一次确保所有模块都到位
    required <- c("numpy", "mpmath", "pandas", "scipy")
    missing <- required[!vapply(required, reticulate::py_module_available, logical(1))]
    if (length(missing)) {
      packageStartupMessage(
        sprintf("❌ 无法加载 Python 脚本，因为以下模块仍缺失：%s\n", paste(missing, collapse = ", ")),
        "  请先在 R 中运行：\n",
        sprintf("    reticulate::py_install(c(%s))",
                paste(sprintf("'%s'", missing), collapse = ", "))
      )
      return()
    }

    # 真正去加载 py_functions.py
    python_file <- system.file("python", "py_functions.py", package = pkgname)
    if (file.exists(python_file)) {
      reticulate::source_python(python_file)
      packageStartupMessage("✅ 已成功加载 Python 函数 py_functions.py")
    } else {
      packageStartupMessage("❌ 未找到 py_functions.py，跳过 Python 加载。")
    }
  }
}

