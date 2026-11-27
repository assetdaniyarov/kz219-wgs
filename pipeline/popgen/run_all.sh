#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

# ============================================================================
# run_all.sh — the whole protocol, one command, fully logged.
#
#   ./run_all.sh              run every stage
#   ./run_all.sh 04 05        run only those stages
#
# Each stage can be rerun. Logs go to logs/<stage>_<ts>.log and the exit
# status of every stage to logs/status_<ts>.tsv, so a partial failure is
# visible rather than buried in scrollback.
#
# set -o pipefail is what makes `if bash stage.sh | tee log` test the script's
# exit status and not tee's. Do not remove it.
# ============================================================================

STAGES=(00_env 01_prepare 02_merge 03_relatedness 04_pca 05_admixture 06_fst 07_roh 08_project_aadr 09_kz_tribes)

# Dependency map, for choosing what to rerun after a config change:
#   KZ_GQ_MIN, GENO_PRE, MIND_*, HET_*        -> 01 and everything after
#   KEEP_NATIONS, DROP_SUPERPOPS              -> 02 and everything after
#   KING_*, PRUNE_*, PCA_MAX_PER_NATION,
#   COHORT_CAP, ADMIX_MIN/MAX_SNPS            -> 03 04 05 06 07 08
#   ADMIX_*                                   -> 05
#   FST_*                                     -> 06
#   ROH_*                                     -> 07
#   TRIBE_MIN_N, PERM_N, KZ_K_*               -> 09 only
# 09 reads work/panels/kz and groupA/king_in, so it needs 01 and 03 but not
# 04-08; it can be rerun on its own.

if [[ $# -gt 0 ]]; then
    SELECTED=()
    for want in "$@"; do
        for s in "${STAGES[@]}"; do [[ "${s}" == "${want}"* ]] && SELECTED+=("${s}"); done
    done
    STAGES=("${SELECTED[@]}")
fi

TS=$(date +%Y%m%d_%H%M%S)
STATUS="${LOG_DIR}/status_${TS}.tsv"
printf 'stage\tstart\tend\tseconds\tstatus\n' > "${STATUS}"

for S in "${STAGES[@]}"; do
    LOG="${LOG_DIR}/${S}_${TS}.log"
    T0=$(date +%s)
    echo "==================================================================="
    echo ">>> ${S}   ->   ${LOG}"
    echo "==================================================================="
    if bash "./${S}.sh" 2>&1 | tee "${LOG}"; then
        RC=ok
    else
        RC=FAILED
    fi
    T1=$(date +%s)
    printf '%s\t%s\t%s\t%s\t%s\n' "${S}" "${T0}" "${T1}" "$((T1-T0))" "${RC}" >> "${STATUS}"
    if [[ "${RC}" == "FAILED" ]]; then
        echo "!!! ${S} failed, stopping. See ${LOG}"
        column -t -s $'\t' "${STATUS}"
        exit 1
    fi
done

echo
echo "=== all stages complete ==="
column -t -s $'\t' "${STATUS}"
echo
echo "key tables:"
for F in "${RES_DIR}"/*/pca/pca_eurasia.tsv \
         "${RES_DIR}"/*/admixture/Q_projected.K*.tsv \
         "${RES_DIR}"/*/fst/*/fst_pairs.tsv \
         "${RES_DIR}"/*/roh/froh_by_population.tsv \
         "${RES_DIR}"/fst_platform_concordance/stats.tsv \
         "${RES_DIR}"/kz_tribes/pc_assoc_*.tsv; do
    [[ -s "${F}" ]] && printf '  %s\n' "${F}"
done
echo
echo "results tree:"
find "${RES_DIR}" -maxdepth 3 -type d | sort
