# refdata changelog

Versioned independently of the image. Reference data moves on upstream's clock,
not `scdock-r-dev`'s, so changes here never bump the root `VERSION` and never
appear in `docs/changelog.md`.

## [Unreleased]

### Added
- **`refdata/` — shared reference-data cache tooling.** Generalises the
  snapshot pattern (dated snapshot dir → verify → `MANIFEST.json` → atomic
  `current` flip → prune) into `refdata/refcache.sh` plus pluggable
  `refdata/sources/*.sh`, run by one small `refdata-fetcher` image. Ships two
  sources: `cistarget` (aertslab motif feathers + motif2tf, sha1-verified,
  set-selectable via `CISTARGET_SETS`, `--dry-run` reports upstream sizes
  before you commit to a 73 GB pull) and `coresh` (Synapse `syn66227307`,
  ported off the standalone `coresh-updater`). The repo holds the mechanism
  only — the bytes stay on a host path. See [../refdata/README.md](../refdata/README.md).
- **`init-container.sh --refcache PATH`** binds that cache `:ro` at `/refcache`
  and exports `REFCACHE_ROOT`, so analysis code resolves reference paths from
  an env var instead of hardcoding a host layout or a snapshot tag.

### Fixed
- **Seven R packages were silently absent from the image — including
  `chromVAR` and `motifmatchr`.** Packages that still declare
  `CXX_STD = CXX11` in `src/Makevars` were compiled with `-std=gnu++11`, and
  current RcppArmadillo hard-errors `*** C++14 compiler required` from
  `compiler_check.hpp`. Because `safe_install()` downgrades install failures to
  warnings, the build stayed green and the packages just never appeared.
  Confirmed present in v0.5.10 too, so `chromVAR` has been missing for at least
  two releases while its dependencies (`TFBSTools`, `JASPAR2022`) were installed
  — the reason chromVAR workflows had nothing to run.
  Fixed with a site-wide `/usr/local/lib/R/etc/Makevars.site` raising
  `CXX11STD` / `CXX14STD` to `-std=gnu++17`. Verified to recover `chromVAR`
  and `motifmatchr` (both 1.32.0).

### Known issues
- **`Rfast`, `WGCNA`, `brms`, `mbkmeans` and `rliger` are still absent.** They
  were missing before this release too, and the C++ standard fix did not
  recover them — each installs cleanly at runtime in the finished image, so
  `CXX_STD` was not their blocker. Root cause is undiagnosed, and the build is
  currently not diagnosable: `install_core.R` reports progress via `message()`
  and failures via `warning()`, and only 6 `Installing ...` lines survive in
  the build log for 563 installed packages, with zero `Failed to install`
  warnings recorded. Fixing the reporting is a prerequisite for fixing the
  packages — `safe_install()` should write a machine-readable failure report
  to `/opt/settings/` rather than emitting warnings that get lost.
  Workaround: `install.packages("Rfast")` etc. into the user library.

### Fixed
- **Prune could `rm -rf` unrelated directories.** `refcache.sh` enumerated every
  directory in a source dir with a bare `find -type d`, having dropped the
  `${SYNAPSE_ID}_*` glob guard the original `coresh-updater` had. Reproduced
  destroying a sibling directory holding user data. Prunable now requires both a
  matching tag glob and a `MANIFEST.json`, ordered by `downloaded_at` rather than
  filename — lexical order could otherwise delete the newest snapshot.
- **Prune refuses to run when `current` is not a symlink**, which would otherwise
  have made the live snapshot its own prune candidate.
- **coresh verify gate accepted partial downloads.** A hardcoded floor of 80
  chunks against real counts of 89/85 would pass a run that silently lost 9 human
  and 5 mouse chunks, then flip `current`. Now regression-checks against the
  previous snapshot (`CORESH_ALLOW_SHRINK=1` to override).
- **`flock` per source**, so concurrent runs cannot interleave into one snapshot.
- **`--keep` validated at parse time** instead of failing after the download, in
  the window between the flip and the prune.
- **`REFCACHE_HOST` vs `REFCACHE_ROOT` conflation in the README** — the fetch
  examples used the container path as the host bind source.

### Added (tooling)
- `--verify-only`: re-check a snapshot already on disk against its stored
  `.sha1sum.txt` siblings. The missing half of the recovery story for an
  unbacked-up cache.
- `produced_by` in `MANIFEST.json` (refdata revision + sha256 of driver and
  source script), so a snapshot can be traced to the code that made it.
- Partial-snapshot path reported on failure rather than left silently occupying
  tens of GB.
