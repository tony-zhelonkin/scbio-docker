# Install and verify bulkiRNA, pinned by commit.
#
# This is a SEPARATE script, run from a SEPARATE late layer, for one reason:
# bulkiRNA changes far more often than anything else in this image. It used to
# be one entry in install_core.R's `github_packages`, which meant a one-line
# version bump invalidated the COPY that gates every CRAN, Bioconductor and
# GitHub R install, plus TinyTeX, the Python venv, samtools, bedtools and every
# runtime layer after them. A release cost a full rebuild to change 40 bytes.
#
# The pin arrives through the environment, set from a Dockerfile ARG, so bumping
# a version does not touch this file either. Only the layer that reads the ARG
# invalidates.
#
# Kept deliberately self-contained rather than sharing install_core.R's
# install_gh_pkg(): that helper exists to install ten heterogeneous repositories
# whose package names do not match their slugs. This installs one pinned, public,
# pure-R package and then proves its identity. Different jobs, so a shared
# abstraction would have to serve both and would drift toward whichever changed
# last.

BULKIRNA_VERSION <- Sys.getenv("BULKIRNA_VERSION", "")
BULKIRNA_SHA <- Sys.getenv("BULKIRNA_SHA", "")
BULKIRNA_REPO <- Sys.getenv("BULKIRNA_REPO", "tony-zhelonkin/bulkiRNA")

if (!nzchar(BULKIRNA_VERSION) || !nzchar(BULKIRNA_SHA)) {
  stop("BULKIRNA_VERSION and BULKIRNA_SHA must both be set. ",
       "They come from Dockerfile ARGs; an unset value means the ARG was not ",
       "declared inside this build stage -- a global ARG has to be redeclared ",
       "after every FROM.", call. = FALSE)
}
if (nchar(BULKIRNA_SHA) < 40L) {
  stop("BULKIRNA_SHA must be a full 40-character commit hash, not a tag or an ",
       "abbreviation: a tag can move and an abbreviation can become ambiguous. ",
       "Got ", nchar(BULKIRNA_SHA), " characters.", call. = FALSE)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  utils::install.packages("remotes", repos = "https://cloud.r-project.org")
}

slug <- paste0(BULKIRNA_REPO, "@", BULKIRNA_SHA)
message("Installing ", slug)

# api.github.com intermittently drops HTTP/2 streams mid-transfer, which loses
# the install outright. Retry, then fall back to a shallow clone.
ok <- FALSE
for (attempt in 1:3) {
  ok <- tryCatch({
    remotes::install_github(slug, quiet = TRUE, upgrade = "never")
    TRUE
  }, error = function(e) {
    message(sprintf("install_github attempt %d/3 failed: %s", attempt,
                    conditionMessage(e)))
    FALSE
  })
  if (ok) break
  Sys.sleep(5 * attempt)
}

if (!ok) {
  message("Falling back to git clone + install_local")
  Sys.setenv(GIT_HTTP_VERSION = "HTTP/1.1", R_LIBCURL_HTTP_VERSION = "1.1")
  tmp <- tempfile()
  dir.create(tmp)
  dest <- file.path(tmp, "bulkiRNA")
  # Clone the default branch, then check out the pinned commit: a bare
  # `--branch <sha>` is not valid for a commit hash, only for a ref name.
  system2("git", c("-c", "http.version=HTTP/1.1", "clone",
                   sprintf("https://github.com/%s.git", BULKIRNA_REPO), dest))
  system2("git", c("-C", dest, "checkout", "--detach", BULKIRNA_SHA))
  remotes::install_local(dest, quiet = TRUE, upgrade = "never",
                         dependencies = TRUE)
  unlink(tmp, recursive = TRUE)
}

# --- identity, verified rather than assumed ---------------------------------
# "0.5.0" once named both a tag and 50 later commits, which is why the version
# alone is not enough. A build that installs the wrong commit must fail here,
# not surface months later as a figure nobody can reproduce.
if (!requireNamespace("bulkiRNA", quietly = TRUE)) {
  stop("bulkiRNA was requested but is not loadable after installation.",
       call. = FALSE)
}
found_version <- as.character(utils::packageVersion("bulkiRNA"))
found_sha <- tryCatch(utils::packageDescription("bulkiRNA")$RemoteSha,
                      error = function(e) NULL)

problems <- character(0)
if (!identical(found_version, BULKIRNA_VERSION)) {
  problems <- c(problems, sprintf("pinned version %s, installed %s",
                                  BULKIRNA_VERSION, found_version))
}
# install_local records no RemoteSha, so the fallback path cannot prove the
# commit. Say so rather than passing silently.
if (is.null(found_sha) || !nzchar(found_sha)) {
  problems <- c(problems,
                "no RemoteSha recorded (installed via the clone fallback), so ",
                "the commit could not be verified")
} else if (!startsWith(BULKIRNA_SHA, found_sha)) {
  problems <- c(problems, sprintf("pinned commit %s, installed %s",
                                  BULKIRNA_SHA, found_sha))
}

message(sprintf("bulkiRNA identity: version %s, commit %s, %d exports",
                found_version,
                if (is.null(found_sha) || !nzchar(found_sha)) "not recorded"
                else found_sha,
                length(getNamespaceExports("bulkiRNA"))))

if (length(problems)) {
  stop("bulkiRNA identity check failed: ", paste(problems, collapse = "; "),
       ".", call. = FALSE)
}

# --- optional-dependency report ---------------------------------------------
# bulkiRNA's Suggests are optional by design, so this keeps its own file rather
# than polluting install_failures.csv, whose meaning is "requested and failed".
# Moved here with the install: a report about a package belongs next to the
# layer that installs it, or it silently reports on the wrong version.
message("\n=== bulkiRNA OPTIONAL DEPENDENCY REPORT ===")
tryCatch({
  deps <- bulkiRNA::bulkirna_check_deps(features = "all", quiet = FALSE,
                                        error = FALSE)
  if (!dir.exists("/opt/settings")) dir.create("/opt/settings", recursive = TRUE)
  out <- "/opt/settings/bulkirna_optional_deps.csv"
  utils::write.csv(as.data.frame(deps), out, row.names = FALSE)
  absent <- deps$package[!deps$installed]
  message(sprintf("bulkiRNA optional dependencies: %d of %d present%s",
                  sum(deps$installed), nrow(deps),
                  if (length(absent))
                    sprintf("; absent: %s", paste(absent, collapse = ", "))
                  else ""))
  message(sprintf("=== report: %s ===", out))
}, error = function(e) {
  # Observational, unlike the identity check above: a missing Suggests makes the
  # image reduced, not wrong.
  message("bulkiRNA optional dependency preflight failed: ",
          conditionMessage(e))
})
