# R/load_py_function.R

.onLoad <- function(libname, pkgname) {
  # 1. 先保证 reticulate 包存在
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    packageStartupMessage("[aBKMR] 请先安装 reticulate 包。")
    return()
  }

  # 2. 不触发 Python 会话，探测系统／conda 是否已有可用 Python
  has_py <- reticulate::py_available(initialize = FALSE)

  if (!has_py) {
    # 2.1 如果完全没 Python，就装 Miniconda 并绑定
    packageStartupMessage("[aBKMR] 未检测到任何 Python，正在安装 Miniconda…")
    reticulate::install_miniconda()
    reticulate::use_miniconda(reticulate::miniconda_path(), required = TRUE)

    # 2.2 一次性用 conda 安装所有依赖
    deps <- c("numpy", "pandas", "scipy", "mpmath")
    packageStartupMessage("[aBKMR] 安装 Python 依赖：", paste(deps, collapse = ", "))
    reticulate::conda_install(
      envname  = reticulate::miniconda_path(),
      packages = deps,
      channel  = "conda-forge"
    )
  } else {
    # 3. 已有 Python，真正初始化以拿到配置
    cfg <- reticulate::py_config()
    # 判定是否是在 reticulate 管理的 conda env
    in_conda <- !is.null(cfg$condaenv) && nzchar(cfg$condaenv)

    # 4. 只对缺失的模块做安装
    deps    <- c("numpy", "pandas", "scipy", "mpmath")
    missing <- deps[!vapply(deps, function(pkg) {
      tryCatch(reticulate::py_module_available(pkg), error = function(e) FALSE)
    }, logical(1))]

    if (length(missing)) {
      if (in_conda) {
        packageStartupMessage("[aBKMR] Conda 环境缺少：", paste(missing, collapse = ", "))
        reticulate::conda_install(
          envname  = cfg$condaenv,
          packages = missing,
          channel  = "conda-forge"
        )
      } else {
        packageStartupMessage("[aBKMR] 系统 Python 缺少：", paste(missing, collapse = ", "))
        reticulate::py_install(
          packages = missing,
          pip      = TRUE
        )
      }
    } else {
      packageStartupMessage("[aBKMR] 检测到所有 Python 依赖，无需安装。")
    }
  }

  # 5. 最后，再加载你的 Python 脚本
  pyfile <- system.file("python", "py_functions.py", package = pkgname)
  if (pyfile == "") {
    stop("找不到 Python 脚本 'py_functions.py'。")
  }
  reticulate::source_python(pyfile)
}




utils::globalVariables(c(
  "ComputePostmeanHnew.exact", "K", "Vinv", "est",
  "variable", "variable1", "variable2", "z1", "."
))
