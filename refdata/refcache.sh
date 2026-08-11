#!/usr/bin/env bash
# refcache.sh — generic driver for versioned reference-data snapshots.
#
# Generalises the pattern that coresh-updater established: download into a
# dated snapshot dir, verify it, write a MANIFEST.json, flip a `current`
# symlink atomically, prune old snapshots. Every source under sources/ plugs
# into this instead of reimplementing it.
#
# A source script is sourced (not executed) and must define:
#
#   SOURCE_NAME           short id; becomes $CACHE_ROOT/<name>/
#   SOURCE_SNAPSHOT_TAG   snapshot dir name (default: <name>_<UTC date>)
#   src_fetch  <dir>      download everything into $dir
#   src_verify <dir>      exit non-zero if $dir is incomplete/corrupt
#   src_manifest_extra <dir>   OPTIONAL; emit extra JSON fields (no outer braces,
#                              must end with a trailing comma)
#
# Usage:
#   refcache.sh <source> [--dry-run] [--keep N] [--force]
#   refcache.sh --list
#
# Env:
#   CACHE_ROOT       where snapshots live (default /cache; bind-mount target)
#   KEEP_SNAPSHOTS   old snapshots to keep beside current (default 1)
#   DRY_RUN          non-empty = print actions only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCES_DIR="${SCRIPT_DIR}/sources"

: "${CACHE_ROOT:=/cache}"
: "${KEEP_SNAPSHOTS:=1}"
FORCE=0

log() { echo "[refcache] $*"; }
die() { echo "[refcache] ERROR: $*" >&2; exit 1; }

list_sources() {
  echo "Available sources:"
  for f in "${SOURCES_DIR}"/*.sh; do
    [ -e "$f" ] || continue
    printf '  %-14s %s\n' "$(basename "$f" .sh)" \
      "$(sed -n 's/^# *desc: *//p' "$f" | head -1)"
  done
}

[ $# -ge 1 ] || { list_sources; exit 1; }

case "$1" in
  --list|-l) list_sources; exit 0 ;;
  --help|-h)
    sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    list_sources
    exit 0 ;;
esac

SOURCE="$1"; shift
VERIFY_ONLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY_RUN=1; shift ;;
    # Validated here so a bad value fails before the download, not after.
    --keep)
      [ $# -ge 2 ] || die "--keep requires a value"
      case "$2" in (''|*[!0-9]*) die "--keep must be a non-negative integer, got: $2" ;; esac
      KEEP_SNAPSHOTS="$2"; shift 2 ;;
    --force)   FORCE=1; shift ;;
    # Re-verify what is already on disk; no fetch, flip or prune.
    --verify-only) VERIFY_ONLY=1; shift ;;
    *) die "unknown option: $1" ;;
  esac
done

SOURCE_FILE="${SOURCES_DIR}/${SOURCE}.sh"
[ -f "$SOURCE_FILE" ] || { echo "No such source: ${SOURCE}" >&2; list_sources; exit 1; }

# Defaults a source may override.
SOURCE_NAME="$SOURCE"
SOURCE_SNAPSHOT_TAG=""
src_manifest_extra() { :; }

# shellcheck disable=SC1090
. "$SOURCE_FILE"

for fn in src_fetch src_verify; do
  declare -F "$fn" >/dev/null || die "source ${SOURCE} does not define ${fn}()"
done

SOURCE_DIR="${CACHE_ROOT}/${SOURCE_NAME}"
: "${SOURCE_SNAPSHOT_TAG:=${SOURCE_NAME}_$(date -u +%Y%m%d)}"
TARGET_DIR="${SOURCE_DIR}/${SOURCE_SNAPSHOT_TAG}"

# Snapshot names for this source, for pruning; derived from the tag's prefix
# (coresh -> syn66227307_*). Overridable by a source.
: "${SNAPSHOT_GLOB:=${SOURCE_SNAPSHOT_TAG%_*}_*}"

log "source:   ${SOURCE_NAME}"
log "snapshot: ${TARGET_DIR}"

# Checked before the exists-guard: cost estimates must work on an existing tag.
if [ -n "${DRY_RUN:-}" ]; then
  log "DRY_RUN set — would fetch ${SOURCE_NAME} into ${TARGET_DIR}, verify, then flip current."
  [ -d "${TARGET_DIR}" ] && log "NOTE: ${TARGET_DIR} already exists; a real run would need --force."
  declare -F src_dry_run >/dev/null && src_dry_run
  exit 0
fi

mkdir -p "${SOURCE_DIR}"

# Serialise runs of the same source; concurrent runs corrupt the snapshot.
exec 9>"${SOURCE_DIR}/.lock"
if ! flock -n 9; then
  die "another refcache run holds the lock on ${SOURCE_NAME} (${SOURCE_DIR}/.lock). Wait for it or kill it."
fi

if [ "${VERIFY_ONLY}" -eq 1 ]; then
  CHECK_DIR="${TARGET_DIR}"
  [ -d "${CHECK_DIR}" ] || CHECK_DIR="${SOURCE_DIR}/current"
  [ -d "${CHECK_DIR}" ] || die "nothing to verify: neither ${TARGET_DIR} nor ${SOURCE_DIR}/current exists."
  log "verify-only against $(readlink -f "${CHECK_DIR}")"
  src_verify "${CHECK_DIR}" || die "VERIFICATION FAILED for ${CHECK_DIR}"
  log "verification passed."
  exit 0
fi

if [ -d "${TARGET_DIR}" ] && [ "${FORCE}" -eq 0 ]; then
  die "snapshot ${SOURCE_SNAPSHOT_TAG} already exists. Remove it or pass --force to resume into it."
fi

mkdir -p "${TARGET_DIR}"

# Partial snapshots are kept (resumable with --force) but can be tens of GB,
# so report the path rather than leaving it silently.
orphan_notice() {
  local rc=$?
  if [ $rc -ne 0 ] && [ -d "${TARGET_DIR}" ]; then
    echo "[refcache] run failed (exit ${rc}). Partial snapshot left at:" >&2
    echo "[refcache]   ${TARGET_DIR}  ($(du -sh "${TARGET_DIR}" 2>/dev/null | cut -f1))" >&2
    echo "[refcache] Resume with --force, or delete it to reclaim the space." >&2
  fi
}
trap orphan_notice EXIT

log "fetching..."
src_fetch "${TARGET_DIR}"

log "verifying..."
src_verify "${TARGET_DIR}" || die "verification failed for ${TARGET_DIR} — NOT flipping current."

# Excludes MANIFEST.json itself, so a later `du -sb` reads slightly higher.
TOTAL_BYTES=$(du -sb "${TARGET_DIR}" | awk '{print $1}')

{
  echo "{"
  echo "  \"source\": \"${SOURCE_NAME}\","
  echo "  \"snapshot_tag\": \"${SOURCE_SNAPSHOT_TAG}\","
  echo "  \"downloaded_at\": \"$(date -u +%Y-%m-%dT%H:%M:%SZ)\","
  # Lets a snapshot on disk be traced back to the code that produced it.
  echo "  \"produced_by\": {"
  echo "    \"refdata_revision\": \"${REFDATA_REVISION:-unknown}\","
  echo "    \"source_script_sha256\": \"$(sha256sum "${SOURCE_FILE}" | awk '{print $1}')\","
  echo "    \"driver_sha256\": \"$(sha256sum "${BASH_SOURCE[0]}" | awk '{print $1}')\""
  echo "  },"
  src_manifest_extra "${TARGET_DIR}"
  echo "  \"total_size_bytes\": ${TOTAL_BYTES}"
  echo "}"
} > "${TARGET_DIR}/MANIFEST.json"

# Atomic-ish symlink flip. `mv -T` over an existing symlink is a rename(2),
# so readers never observe a missing `current`.
NEW_LINK="${SOURCE_DIR}/current.new"
ln -sfn "${SOURCE_SNAPSHOT_TAG}" "${NEW_LINK}"
mv -T "${NEW_LINK}" "${SOURCE_DIR}/current"
log "current -> ${SOURCE_SNAPSHOT_TAG}"

# Prune old snapshots. Prunable only if the name matches this source's tag glob
# AND the dir holds a MANIFEST.json; ordered by downloaded_at, not filename.
prune_snapshots() {
  local current_target=""
  if [ -e "${SOURCE_DIR}/current" ] && [ ! -L "${SOURCE_DIR}/current" ]; then
    log "WARNING: ${SOURCE_DIR}/current is not a symlink — refusing to prune."
    return 0
  fi
  [ -L "${SOURCE_DIR}/current" ] && current_target="$(readlink "${SOURCE_DIR}/current")"

  local -a candidates=()
  local d name stamp
  while IFS= read -r d; do
    name="$(basename "$d")"
    [ "${name}" = "${current_target}" ] && continue
    [ -f "${d}/MANIFEST.json" ] || { log "skipping ${name} (no MANIFEST.json — not ours)"; continue; }
    stamp="$(sed -n 's/.*"downloaded_at" *: *"\([^"]*\)".*/\1/p' "${d}/MANIFEST.json" | head -1)"
    candidates+=("${stamp:-0000} ${name}")
  done < <(find "${SOURCE_DIR}" -maxdepth 1 -mindepth 1 -type d -name "${SNAPSHOT_GLOB}" -print)

  # current is excluded from candidates, so compare directly.
  local total=${#candidates[@]}
  (( total > KEEP_SNAPSHOTS )) || return 0
  local old
  while IFS= read -r old; do
    old="${old#* }"
    log "pruning ${old}"
    rm -rf "${SOURCE_DIR:?}/${old:?}"
  done < <(printf '%s\n' "${candidates[@]}" | sort | head -n "$((total - KEEP_SNAPSHOTS))")
}
prune_snapshots

trap - EXIT
log "done. $(numfmt --to=iec "${TOTAL_BYTES}" 2>/dev/null || echo "${TOTAL_BYTES} bytes") at ${SOURCE_DIR}/current"
