#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# setup_meta.sh — run once, before check_inputs.sh.
#
# Copies superpop_map.tsv into meta/ (it is edited by hand, so it lives with
# the results rather than with the code) and prints the real header of every
# metadata file so the MC_* column names in config.sh can be set to match.
#
# Nothing else is copied. config.sh points at the metadata where it already
# lives in raw_datasets/, so there is exactly one copy of each file.
# ============================================================================

SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"

# --- superpop_map.tsv -------------------------------------------------------
if [[ -s "${MF_SUPERPOP}" ]]; then
    log "superpop_map.tsv already in meta/, leaving it alone"
elif [[ -s "${SCRIPTDIR}/superpop_map.tsv" ]]; then
    cp "${SCRIPTDIR}/superpop_map.tsv" "${MF_SUPERPOP}"
    log "copied superpop_map.tsv -> ${MF_SUPERPOP}"
else
    log "WARNING: no superpop_map.tsv in ${SCRIPTDIR}, 01 will abort without it"
fi

# --- headers ----------------------------------------------------------------
show() {   # show <label> <file> <expected id> <expected nation> <expected region>
    local LABEL=$1 F=$2 EID=$3 ENAT=$4 EREG=$5
    printf '\n\033[1m=== %s ===\033[0m\n%s\n' "${LABEL}" "${F}"
    if [[ ! -s "${F}" ]]; then
        printf '  \033[31mMISSING\033[0m\n'; return
    fi

    # .janno and .psam start the header with a comment character; some files
    # lead with commentary lines before the real header.
    local HDR
    HDR=$(grep -m1 -P '\t' "${F}" 2>/dev/null || head -1 "${F}")

    printf '  rows: %s\n  columns:\n' "$(( $(wc -l < "${F}") - 1 ))"
    printf '%s\n' "${HDR}" | tr '\t' '\n' | sed 's/^#//' | nl -w4 -s'  ' | sed 's/^/   /'

    for SPEC in "ID:${EID}" "NATION:${ENAT}" "REGION:${EREG}"; do
        local ROLE="${SPEC%%:*}" WANT="${SPEC#*:}"
        if printf '%s\n' "${HDR}" | tr '\t' '\n' | sed 's/^#//' | grep -qx "${WANT}"; then
            printf '  \033[32mOK\033[0m    %-7s %s\n' "${ROLE}" "${WANT}"
        else
            printf '  \033[31mNOT FOUND\033[0m %-7s %s  <- fix the matching MC_* in config.sh\n' \
                   "${ROLE}" "${WANT}"
        fi
    done
}

show "1000 Genomes (psam)" "${MF_1KG}"  "${MC_1KG_ID}"  "${MC_1KG_NATION}"  "${MC_1KG_REGION}"
show "SGDP"                "${MF_SGDP}" "${MC_SGDP_ID}" "${MC_SGDP_NATION}" "${MC_SGDP_REGION}"
show "HGDP"                "${MF_HGDP}" "${MC_HGDP_ID}" "${MC_HGDP_NATION}" "${MC_HGDP_REGION}"
show "AADR (janno)"        "${MF_AADR}" "${MC_AADR_ID}" "${MC_AADR_NATION}" "${MC_AADR_REGION}"

printf '\n\033[1m=== Kazakh cohort ===\033[0m\n'
if [[ -s "${MF_KZ}" ]]; then
    show "Kazakh" "${MF_KZ}" "${MC_KZ_ID}" "${MC_KZ_NATION}" "${MC_KZ_NATION}"
else
    printf '%s\n  absent, every Kazakh will be labelled Nation=Kazakh.\n' "${MF_KZ}"
    printf '  To carry region or zhuz into the figures, create it as:\n'
    printf '    IID<TAB>Nation\n'
fi

printf '\n=== next ===\n'
printf '  1. fix any NOT FOUND above by editing MC_* in config.sh\n'
printf '  2. review %s, especially the HGDP Central_South_Asia split\n' "${MF_SUPERPOP}"
printf '  3. bash check_inputs.sh\n'
