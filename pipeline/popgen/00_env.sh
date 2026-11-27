#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 00_env.sh — provenance for the Code Availability statement.
#
# Writes tool versions, the pinned conda environment, every parameter in
# config.sh and the identity of each input dataset into results/provenance/.
#
# PARAMETERS ARE DISCOVERED, NOT LISTED. The previous revision carried a
# hand-maintained array of parameter names and reported the ones that had gone
# missing. That inverts the problem: a parameter ADDED to config.sh and never
# added here was recorded nowhere and produced no warning, which is how
# MIND_AADR, GENO_AADR and N_KZ_EXPECTED ended up outside the provenance
# record. The name list is now parsed out of config.sh itself, so the two files
# cannot disagree.
#
# Checksums are off by default: md5 over an 82 GB VCF plus the reference panels
# takes the better part of an hour. Size and mtime identify the files well
# enough for day-to-day reruns. Run with CHECKSUM=1 once, for the version that
# goes into the supplement.
#
#   CHECKSUM=1 bash 00_env.sh
# ============================================================================

need "${PLINK1}"; need "${PLINK2}"; need "${ADMIXTURE}"

CHECKSUM="${CHECKSUM:-0}"

# Captured before the provenance block, not inside it. `admixture` with no
# arguments prints its banner and exits non-zero; under `set -o pipefail` that
# makes the whole pipeline non-zero even though grep and head succeeded, so an
# inline `|| echo 'not found'` fires IN ADDITION to the captured banner and the
# substitution returns two lines. The banner is also padded with asterisks and
# runs of spaces, which do not belong in a TSV field.
ADMIX_VER=$(${ADMIXTURE} 2>&1 | grep -i version | head -1 \
            | sed 's/[*]//g; s/^[[:space:]]*//; s/[[:space:]]*$//; s/[[:space:]][[:space:]]*/ /g' || true)
ADMIX_VER="${ADMIX_VER:-version string not found}"
V="${RES_DIR}/provenance"
mkdir -p "${V}"
CFG="$(dirname "$0")/config.sh"

# printf, not echo: bash echo does not expand \t and the file has to be a real TSV.
{
  printf 'key\tvalue\n'
  printf 'run_date\t%s\n'   "$(date -Iseconds)"
  printf 'host\t%s\n'       "$(hostname)"
  printf 'user\t%s\n'       "$(id -un)"
  printf 'kernel\t%s\n'     "$(uname -sr)"
  printf 'project\t%s\n'    "${PROJECT}"
  printf 'plink1.9\t%s\n'   "$(${PLINK1} --version 2>&1 | head -1)"
  printf 'plink2\t%s\n'     "$(${PLINK2} --version 2>&1 | head -1)"
  printf 'admixture\t%s\n'  "${ADMIX_VER}"
  printf 'python\t%s\n'     "$(${PYTHON} --version 2>&1)"
  printf 'numpy\t%s\n'      "$(${PYTHON} -c 'import numpy;print(numpy.__version__)' 2>&1 | tail -1)"
  printf 'scipy\t%s\n'      "$(${PYTHON} -c 'import scipy;print(scipy.__version__)' 2>&1 | tail -1)"
  printf 'bash\t%s\n'       "${BASH_VERSION}"
  printf 'git_commit\t%s\n' "$(git -C "$(dirname "$0")" rev-parse --short HEAD 2>/dev/null || echo 'not a git repo')"
  # Braces + `|| true`: outside a git repo `git status` exits 128, and under
  # `set -o pipefail` that fails the whole `git | wc -l` pipeline, which fires
  # the ERR trap even though nothing is wrong. The trap is there to catch real
  # failures; a false positive on line 60 of the provenance stage teaches you to
  # ignore it.
  printf 'git_dirty\t%s\n'  "$( { git -C "$(dirname "$0")" status --porcelain 2>/dev/null || true; } | wc -l )"
} > "${V}/versions.tsv"

conda env export -n "${CONDA_ENV}" --no-builds > "${V}/environment.yml" 2>/dev/null \
  || log "conda env export failed for '${CONDA_ENV}', record the environment manually"

# A verbatim copy of config.sh goes with the results. The parsed table below is
# for reading; this is the artefact that reproduces the run.
cp "${CFG}" "${V}/config.sh.snapshot"

# ---------------------------------------------------------------------------
# parameters, discovered from config.sh
#
# Matches   export NAME=...   and   NAME=(...)   at the start of a line. Paths,
# binaries and metadata column names are included deliberately: "which HGDP
# file" is as much a parameter of the analysis as MAF is.
# ---------------------------------------------------------------------------
# `=[^=]` and not just `=`: an indented awk line such as `FNR==1 { fn++ }`
# inside a quoted helper otherwise matches, and the run reports FNR as a
# parameter of the study that failed to be set. Awk's own built-ins are
# excluded too, for the same class of false positive.
mapfile -t PARAMS < <(
    grep -oE '^[[:space:]]*(export[[:space:]]+)?[A-Z][A-Z0-9_]*=([^=]|$)' "${CFG}" \
      | sed -E 's/^[[:space:]]*(export[[:space:]]+)?//; s/=.?$//' \
      | grep -vxE 'NR|FNR|NF|FS|OFS|ORS|RS|FILENAME|SUBSEP|RSTART|RLENGTH|CONVFMT|OFMT|ENVIRON' \
      | sort -u || true
)
[[ "${#PARAMS[@]}" -gt 0 ]] || die "could not parse any parameter names out of ${CFG}"

UNSET=()
{
  printf 'parameter\ttype\tvalue\n'
  for p in "${PARAMS[@]}"; do
      if ! declare -p "${p}" >/dev/null 2>&1; then
          printf '%s\t%s\t%s\n' "${p}" "unset" "NOT SET"
          UNSET+=("${p}")
          continue
      fi
      if declare -p "${p}" 2>/dev/null | grep -q '^declare -[aA]'; then
          declare -n _arr="${p}"
          printf '%s\tarray\t%s\n' "${p}" "${_arr[*]}"
          unset -n _arr
      else
          printf '%s\tscalar\t%s\n' "${p}" "${!p}"
      fi
  done
} > "${V}/parameters.tsv"
log "parameters recorded: ${#PARAMS[@]}"

# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
stamp() {   # stamp <file>
    local f=$1
    [[ -f "${f}" ]] || return 0
    if [[ "${CHECKSUM}" -eq 1 ]]; then
        printf '%s\t%s\t%s\t%s\n' "${f}" "$(stat -Lc %s "${f}")" \
               "$(stat -Lc %y "${f}")" "$(md5sum "${f}" | cut -d' ' -f1)"
    else
        printf '%s\t%s\t%s\t%s\n' "${f}" "$(stat -Lc %s "${f}")" \
               "$(stat -Lc %y "${f}")" "not computed (CHECKSUM=0)"
    fi
}

{
  printf 'file\tbytes\tmtime\tmd5\n'
  for f in "${DS_KZ_VCF}" "${DS_1KG}.pgen" "${DS_1KG}.bed" \
           "${DS_SGDP}.bed" "${DS_HGDP}.bed" "${DS_AADR}.bed"; do
      stamp "${f}"
  done
} > "${V}/input_checksums.tsv"

# Hand-edited metadata is part of the method, so it is always checksummed
# regardless of CHECKSUM: these files are small and they are the ones a
# reviewer would want to pin.
{
  printf 'file\tbytes\tmtime\tmd5\n'
  for f in "${MF_SUPERPOP}" "${MF_AADR_KEEP}" "${MF_KZ}" "${MF_TRIBES}" \
           "${MF_PANEL_EXCLUDE}" "${MF_PSEUDOHAPLOID}"; do
      [[ -f "${f}" ]] && printf '%s\t%s\t%s\t%s\n' "${f}" "$(stat -Lc %s "${f}")" \
          "$(stat -Lc %y "${f}")" "$(md5sum "${f}" | cut -d' ' -f1)"
  done
} > "${V}/metadata_checksums.tsv"

log "provenance -> ${V}"
column -t -s $'\t' "${V}/versions.tsv"

if [[ "${#UNSET[@]}" -gt 0 ]]; then
    log "WARNING: ${#UNSET[@]} name(s) parsed from config.sh are not set at runtime:"
    printf '        %s\n' "${UNSET[@]}"
    log "That normally means a conditional assignment; check before publishing."
fi

log "=== 00 done -> 01_prepare.sh ==="
