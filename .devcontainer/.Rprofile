Sys.setenv(TERM_PROGRAM = "vscode")
if (interactive()) {
  init <- file.path(
    Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"),
    ".vscode-R", "init.R"
  )
  if (file.exists(init)) try(source(init), silent = TRUE)
}

# Source: https://github.com/REditorSupport/vscode-R/wiki/Plot-viewer#svg-in-httpgd-webpage
if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
  if ("httpgd" %in% .packages(all.available = TRUE)) {
    options(vsc.plot = FALSE)
    options(device = function(...) {
      httpgd::hgd(silent = TRUE)
      .vsc.browser(httpgd::hgd_url(history = FALSE), viewer = "Beside")
    })
  }
}

# Set CRAN mirror with fallback
local_mirror <- Sys.getenv("CRAN_MIRROR", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = local_mirror))