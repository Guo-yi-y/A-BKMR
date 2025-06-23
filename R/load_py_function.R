# R/load_py_function.R

.py_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # ... 只做探测和提示，不做 source_python() ...
}

.onAttach <- function(libname, pkgname) {
  if (!interactive()) return()
  if (!requireNamespace("reticulate", quietly = TRUE)) return()
  if (inherits(try(reticulate::py_config(), silent = TRUE), "try-error")) return()

  # 再次检查必要的 Python 包
  req <- c("numpy","mpmath","pandas","scipy")
  miss <- req[!vapply(req, reticulate::py_module_available, logical(1))]
  if (length(miss)) {
    packageStartupMessage(
      sprintf("❌ 缺少 Python 包：%s\n", paste(miss, collapse = ", ")),
      "  请先运行：reticulate::py_install(c(", paste(sprintf("'%s'", miss), collapse = ", "), "))"
    )
    return()
  }

  # import py_functions.py
  python_dir <- system.file("python", package = pkgname)
  fpath <- file.path(python_dir, "py_functions.py")
  if (!file.exists(fpath)) {
    packageStartupMessage("❌ 找不到 py_functions.py，跳过 Python 加载。")
    return()
  }

  mod <- reticulate::import_from_path(
    module  = "py_functions",
    path    = python_dir,
    convert = TRUE
  )
  # 把 mod 存到我们自己创建的环境里
  .py_env$mod <- mod
  packageStartupMessage("✅ 已加载 Python 模块 py_functions 到 .py_env$mod")
}

