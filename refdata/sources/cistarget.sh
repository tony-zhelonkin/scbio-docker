#!/usr/bin/env bash
# desc: aertslab cisTarget motif databases (feathers + motif2tf) for SCENIC+/pySCENIC
#
# Sourced by refcache.sh — see that file for the contract.
#
# Sizes (measured from Content-Length, 2026-08):
#   region-based hg38 screen v10_clust   rankings 32.8G + scores 12.9G = 45.7G
#   region-based mm10 screen v10_clust   rankings 16.6G + scores  7.6G = 24.2G
#   gene-based   hg38 refseq_r80 v10     4 files                       =  1.15G
#   gene-based   mm10 refseq_r80 v10     4 files                       =  0.90G
#   motif2tf     hgnc 94M + mgi 108M                                   =  0.20G
#
# Select with CISTARGET_SETS (comma-separated, default = all of the above):
#   region_hg38, region_mm10, gene_hg38, gene_mm10, motif2tf

SOURCE_NAME=cistarget

BASE=https://resources.aertslab.org/cistarget
: "${CISTARGET_SETS:=region_hg38,region_mm10,gene_hg38,gene_mm10,motif2tf}"

_region_files() {   # $1 = hg38|mm10, $2 = homo_sapiens|mus_musculus
  # Separate `local` statements: in one statement $g/$sp expand before assignment.
  local g=$1
  local sp=$2
  local d="${BASE}/databases/${sp}/${g}/screen/mc_v10_clust/region_based"
  local f
  for f in rankings scores; do
    echo "region_based/${g}_screen_v10_clust.regions_vs_motifs.${f}.feather ${d}/${g}_screen_v10_clust.regions_vs_motifs.${f}.feather"
  done
}

_gene_files() {     # $1 = hg38|mm10, $2 = species dir
  local g=$1
  local sp=$2
  local d="${BASE}/databases/${sp}/${g}/refseq_r80/mc_v10_clust/gene_based"
  local w f
  for w in 10kbp_up_10kbp_down 500bp_up_100bp_down; do
    for f in rankings scores; do
      local n="${g}_${w}_full_tx_v10_clust.genes_vs_motifs.${f}.feather"
      echo "gene_based/${n} ${d}/${n}"
    done
  done
}

# Emits "relative/dest/path<space>url" per line for the selected sets.
_manifest_list() {
  local set
  IFS=, read -ra _sets <<< "${CISTARGET_SETS}"
  for set in "${_sets[@]}"; do
    case "$set" in
      region_hg38) _region_files hg38 homo_sapiens ;;
      region_mm10) _region_files mm10 mus_musculus ;;
      gene_hg38)   _gene_files   hg38 homo_sapiens ;;
      gene_mm10)   _gene_files   mm10 mus_musculus ;;
      motif2tf)
        local t
        for t in hgnc mgi; do
          echo "motif2tf/motifs-v10nr_clust-nr.${t}-m0.001-o0.0.tbl ${BASE}/motif2tf/motifs-v10nr_clust-nr.${t}-m0.001-o0.0.tbl"
        done ;;
      *) echo "[cistarget] unknown set: ${set}" >&2; return 1 ;;
    esac
  done
}

src_dry_run() {
  local total=0 dest url sz
  while read -r dest url; do
    sz=$(curl -sSIL --max-time 30 "$url" | grep -i '^content-length' | tail -1 | tr -dc '0-9')
    total=$(( total + ${sz:-0} ))
    printf '  %10s  %s\n' "$(numfmt --to=iec "${sz:-0}")" "$dest"
  done < <(_manifest_list)
  echo "  ----------"
  printf '  %10s  TOTAL (sets: %s)\n' "$(numfmt --to=iec "$total")" "${CISTARGET_SETS}"
}

src_fetch() {
  local dir=$1 dest url
  while read -r dest url; do
    mkdir -p "${dir}/$(dirname "$dest")"
    echo "[cistarget] -> ${dest}"
    # -C - resumes a partial file so --force does not re-pull 30GB.
    curl -fL --retry 5 --retry-delay 10 --retry-all-errors \
         -C - -o "${dir}/${dest}" "$url"
    # Feathers ship a sibling .sha1sum.txt; grab it for the verify step.
    case "$dest" in
      *.feather)
        curl -fsSL --retry 3 -o "${dir}/${dest}.sha1sum.txt" "${url}.sha1sum.txt" || {
          echo "[cistarget] WARN: no sha1sum published for ${dest}" >&2
          rm -f "${dir}/${dest}.sha1sum.txt"
        } ;;
    esac
  done < <(_manifest_list)
}

src_verify() {
  local dir=$1 dest url rc=0 n=0
  while read -r dest url; do
    n=$((n + 1))
    if [ ! -s "${dir}/${dest}" ]; then
      echo "[cistarget] MISSING or empty: ${dest}" >&2; rc=1; continue
    fi
    # Slow on a 33GB feather, but truncation yields plausible-wrong results
    # rather than an error, so it is checked once per refresh.
    if [ -f "${dir}/${dest}.sha1sum.txt" ]; then
      local want got
      want=$(awk '{print $1}' "${dir}/${dest}.sha1sum.txt")
      got=$(sha1sum "${dir}/${dest}" | awk '{print $1}')
      if [ "$want" != "$got" ]; then
        echo "[cistarget] SHA1 MISMATCH: ${dest} (want ${want}, got ${got})" >&2; rc=1
      else
        echo "[cistarget] ok ${dest}"
      fi
    else
      echo "[cistarget] ok (unverified, no sha1) ${dest}"
    fi
  done < <(_manifest_list)
  [ "$n" -gt 0 ] || { echo "[cistarget] no files selected" >&2; return 1; }
  return $rc
}

src_manifest_extra() {
  local dir=$1
  echo "  \"sets\": \"${CISTARGET_SETS}\","
  echo "  \"motif_collection\": \"v10nr_clust\","
  echo "  \"upstream\": \"${BASE}\","
  echo "  \"file_count\": $(find "$dir" -type f ! -name 'MANIFEST.json' ! -name '*.sha1sum.txt' | wc -l),"
}
