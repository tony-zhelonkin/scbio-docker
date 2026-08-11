#!/usr/bin/env bash
# desc: CoReSh GEO chunk compendium from Synapse (syn66227307)
#
# Sourced by refcache.sh — see that file for the contract.
# Requires SYNAPSE_AUTH_TOKEN (view + download) at run time; never bake it in:
#   https://accounts.synapse.org/authenticated/personalaccesstokens
#
# Env:
#   SYNAPSE_ID          default syn66227307
#   CORESH_MIN_CHUNKS   per-species minimum for the verify gate (default 80)

SOURCE_NAME=coresh

: "${SYNAPSE_ID:=syn66227307}"
: "${CORESH_MIN_CHUNKS:=80}"

# Historical naming keeps existing snapshots in one series.
: "${SOURCE_SNAPSHOT_TAG:=${SYNAPSE_ID}_$(date -u +%Y%m%d)}"

_counts() {
  local dir=$1 sp=$2
  find "${dir}/preprocessed_chunks/${sp}" -name '*.qs2' 2>/dev/null | wc -l
}

src_dry_run() {
  echo "  would run: synapse get -r ${SYNAPSE_ID}"
  echo "  requires SYNAPSE_AUTH_TOKEN: ${SYNAPSE_AUTH_TOKEN:+set}${SYNAPSE_AUTH_TOKEN:-NOT SET}"
}

src_fetch() {
  local dir=$1
  : "${SYNAPSE_AUTH_TOKEN:?SYNAPSE_AUTH_TOKEN must be set (pass -e at docker run time)}"
  echo "[coresh] authenticating to Synapse..."
  synapse get -r "${SYNAPSE_ID}" --downloadLocation "${dir}"
}

src_verify() {
  local dir=$1
  local hsa mmu rc=0
  hsa=$(_counts "$dir" hsa)
  mmu=$(_counts "$dir" mmu)
  echo "[coresh] hsa chunks: ${hsa}"
  echo "[coresh] mmu chunks: ${mmu}"

  # Refuse to go backwards vs the previous snapshot; the floor is only a
  # first-snapshot backstop and would pass a partial download on its own.
  local prev="${SOURCE_DIR}/current/MANIFEST.json"
  if [ -f "$prev" ]; then
    local phsa pmmu
    phsa=$(sed -n 's/.*"hsa" *: *\([0-9]*\).*/\1/p' "$prev" | head -1)
    pmmu=$(sed -n 's/.*"mmu" *: *\([0-9]*\).*/\1/p' "$prev" | head -1)
    if [ -n "${phsa:-}" ] && [ "$hsa" -lt "$phsa" ]; then
      echo "[coresh] REGRESSION: hsa ${hsa} < previous snapshot's ${phsa}" >&2; rc=1
    fi
    if [ -n "${pmmu:-}" ] && [ "$mmu" -lt "$pmmu" ]; then
      echo "[coresh] REGRESSION: mmu ${mmu} < previous snapshot's ${pmmu}" >&2; rc=1
    fi
    [ $rc -eq 0 ] && echo "[coresh] counts >= previous snapshot (hsa ${phsa}, mmu ${pmmu})"
    [ $rc -ne 0 ] && echo "[coresh] If upstream genuinely shrank, set CORESH_ALLOW_SHRINK=1." >&2
    [ -n "${CORESH_ALLOW_SHRINK:-}" ] && { echo "[coresh] CORESH_ALLOW_SHRINK set — proceeding anyway." >&2; rc=0; }
  else
    echo "[coresh] no previous snapshot; applying floor of ${CORESH_MIN_CHUNKS}"
  fi

  if [ "$hsa" -lt "${CORESH_MIN_CHUNKS}" ] || [ "$mmu" -lt "${CORESH_MIN_CHUNKS}" ]; then
    echo "[coresh] chunk counts below absolute minimum (${CORESH_MIN_CHUNKS})" >&2
    rc=1
  fi
  return $rc
}

src_manifest_extra() {
  local dir=$1
  # Cheap integrity probe. Array-indexed, not `| head -n1`: under pipefail that
  # kills sha256sum with SIGPIPE and aborts the refresh before the flip.
  local hsa_files=( "${dir}"/preprocessed_chunks/hsa/*.qs2 )
  local mmu_files=( "${dir}"/preprocessed_chunks/mmu/*.qs2 )
  echo "  \"synapse_id\": \"${SYNAPSE_ID}\","
  echo "  \"chunk_counts\": { \"hsa\": $(_counts "$dir" hsa), \"mmu\": $(_counts "$dir" mmu) },"
  echo "  \"sha256_sample\": {"
  echo "    \"hsa_first\": \"$(sha256sum "${hsa_files[0]}" | awk '{print $1}')\","
  echo "    \"mmu_first\": \"$(sha256sum "${mmu_files[0]}" | awk '{print $1}')\""
  echo "  },"
}
