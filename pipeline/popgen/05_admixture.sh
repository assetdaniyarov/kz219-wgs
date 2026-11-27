#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 05_admixture.sh — ADMIXTURE K=K_MIN..K_MAX, N_REPS seeds per K.
#
# THE DESIGN PROBLEM THIS STAGE SOLVES. Reference nations are capped at 25 and
# the cohort is not, so the cohort is ~25% of the Eurasian panel and about nine
# times the size of the largest reference group. ADMIXTURE components follow
# sample size, so a joint fit manufactures a cohort-specific component out of
# the cohort's own weight and then reports it as ancestry. Capping the cohort
# instead throws away the data the paper is about.
#
# The standard resolution is projection:
#
#   1. TRAIN on merged_admix_train — references capped at PCA_MAX_PER_NATION,
#      cohort capped at COHORT_CAP — with N_REPS seeds per K. CV error and the
#      mode structure come from these runs, computed on a balanced sample.
#   2. Take the .P (allele frequencies) of the highest-likelihood run inside the
#      major mode and PROJECT the full panel onto it with `admixture -P`. Every
#      cohort individual gets a Q, and no cohort individual influenced the
#      component definitions.
#   3. As a sensitivity check, fit the uncapped panel jointly at the CV-selected
#      K and correlate the two Q matrices component by component. If capping
#      changed nothing, the correlations are ~1 and the objection is closed with
#      a number rather than an argument.
#
# merged_admix and merged_admix_train share a .bim exactly — 03 verifies it with
# cmp — because ADMIXTURE's -P reads a .P.in whose rows are the target .bed's
# variants in order.
#
# NON-CONVEXITY. One run per K can settle in a local optimum. N_REPS independent
# seeds per K; runs are aligned and grouped into modes by admixture_align.py and
# the major mode is what goes in the figure. CV error is reported as mean +/- SD
# across seeds: a single run gives a point estimate with no spread and "K = x
# minimises CV" does not hold up without one.
#
# INPUT is merged_admix from 03, not merged_pca. The PCA set keeps the African
# outgroup and Oceania so the main PCA can show global position; ADMIXTURE
# cannot afford them. Every continental outgroup left in consumes a component
# before any within-Eurasian structure appears — K=2 separates Africa, K=3
# Oceania, K=4 East from West Eurasia — and Central Asia only starts resolving
# near K=6.
#
# Seeds are deterministic, SEED_BASE + K*100 + REP (+50000 for the joint runs so
# the two sweeps never share a seed), and recorded per run.
#
# TIMING. (K_MAX - K_MIN + 1) * N_REPS runs on the training set, plus one cheap
# projection per K, plus N_REPS joint runs at one K. ADMIX_JOBS processes run
# concurrently with ADMIX_THREADS each rather than one process on all cores:
# ADMIXTURE scales poorly past ~8 threads, so 4x4 beats 2x8 on the same 16
# cores, which matters now that the pruned set is 3-4x larger than it was.
#
# A FAILED RUN NO LONGER KILLS THE STAGE. `wait -n` returns the exit status of
# the job it reaped, and under `set -e` a single ADMIXTURE failure at hour 12
# aborted everything with no summary. Both waits are now guarded and failures
# are counted at the end.
#
# ADMIXTURE writes <bed-basename>.<K>.Q into the working directory, so each run
# executes in its own directory.
# ============================================================================

need "${ADMIXTURE}"
ALIGN="${SCRIPT_DIR}/admixture_align.py"
require_file "${ALIGN}"

# ---------------------------------------------------------------------------
# launch_runs <bfile prefix> <run root> <k_min> <k_max> <n_reps> <seed offset>
# Writes <run root>/runs/K<k>_r<rep>/ and <run root>/cv_raw.tsv.
# ---------------------------------------------------------------------------
launch_runs() {
    local BF=$1 ROOT=$2 KMIN=$3 KMAX=$4 REPS=$5 SOFF=$6
    local STEM; STEM=$(basename "${BF}")
    admix_cache_guard "${BF}" "${ROOT}"
    printf '%s\n' "${STEM}" > "${ROOT}/stem.txt"
    # individuals_run.txt, NOT individuals.txt. run_group writes the FULL
    # panel's list to <R>/individuals.txt and the projected .Q has one row per
    # individual in it; when launch_runs trained on the capped set it used to
    # overwrite that file with the 661 training names, and q-table then failed
    # with "661 individuals but Q has 854 rows" after the six-hour sweep.
    awk '{print $2}' "${BF}.fam" > "${ROOT}/individuals_run.txt"

    local K REP SEED RD
    for K in $(seq "${KMIN}" "${KMAX}"); do
        for REP in $(seq 1 "${REPS}"); do
            SEED=$(( SOFF + K*100 + REP ))
            RD="${ROOT}/runs/K${K}_r${REP}"
            [[ -f "${RD}/${STEM}.${K}.Q" ]] && continue
            mkdir -p "${RD}"
            ln -sf "${BF}.bed" "${RD}/${STEM}.bed"
            ln -sf "${BF}.bim" "${RD}/${STEM}.bim"
            ln -sf "${BF}.fam" "${RD}/${STEM}.fam"
            while [[ "$(jobs -rp | wc -l)" -ge "${ADMIX_JOBS}" ]]; do wait -n || true; done
            log "  launch K=${K} rep=${REP} seed=${SEED}"
            ( cd "${RD}" && ${ADMIXTURE} --cv="${CV_FOLDS}" -j"${ADMIX_THREADS}" \
                -s "${SEED}" "${STEM}.bed" "${K}" > "log.txt" 2>&1 ) &
        done
    done
    wait || true

    local NFAIL=0
    printf 'K\tREP\tSEED\tCV\tLOGLIK\n' > "${ROOT}/cv_raw.tsv"
    for K in $(seq "${KMIN}" "${KMAX}"); do
        for REP in $(seq 1 "${REPS}"); do
            RD="${ROOT}/runs/K${K}_r${REP}"
            SEED=$(( SOFF + K*100 + REP ))
            [[ -f "${RD}/log.txt" ]] || continue
            if [[ ! -f "${RD}/${STEM}.${K}.Q" ]]; then
                log "  FAILED K=${K} rep=${REP}, see ${RD}/log.txt"
                NFAIL=$(( NFAIL + 1 )); continue
            fi
            local CV LL
            CV=$(grep "CV error" "${RD}/log.txt" | tail -1 | awk '{print $NF}')
            # ADMIXTURE prints Loglikelihood on every EM iteration; the final
            # value is the last one, not the first
            LL=$(grep "^Loglikelihood" "${RD}/log.txt" | tail -1 | awk '{print $NF}')
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "${K}" "${REP}" "${SEED}" "${CV:-NA}" "${LL:-NA}" >> "${ROOT}/cv_raw.tsv"
        done
    done

    { printf 'K\tREP\tQ_PATH\n'
      for K in $(seq "${KMIN}" "${KMAX}"); do
        for REP in $(seq 1 "${REPS}"); do
          local Q="${ROOT}/runs/K${K}_r${REP}/${STEM}.${K}.Q"
          [[ -f "${Q}" ]] && printf '%s\t%s\t%s\n' "${K}" "${REP}" "${Q}"
        done
      done
    } > "${ROOT}/runs.tsv"

    awk -F'\t' 'NR>1 && $4!="NA" {n[$1]++; s[$1]+=$4; ss[$1]+=$4*$4}
        END{ printf "K\tn_runs\tCV_mean\tCV_sd\n";
             for (k in n) { m=s[k]/n[k]; v=ss[k]/n[k]-m*m; if (v<0) v=0;
                            printf "%s\t%d\t%.5f\t%.5f\n", k, n[k], m, sqrt(v) } }' \
        "${ROOT}/cv_raw.tsv" | { read -r h; echo "$h"; sort -n; } > "${ROOT}/cv_summary.tsv"

    local NDONE NWANT
    NDONE=$(nbody "${ROOT}/runs.tsv")
    NWANT=$(( (KMAX - KMIN + 1) * REPS ))
    log "  runs completed: ${NDONE} / ${NWANT} (${NFAIL} failed)"
    [[ "${NDONE}" -eq "${NWANT}" ]] \
        || log "  WARNING: ${NWANT} expected; check ${ROOT}/runs/*/log.txt"
    [[ "${NDONE}" -gt 0 ]] || die "no ADMIXTURE run completed under ${ROOT}"
}

run_group() {
    local GROUP=$1
    local G="${WORK_DIR}/${GROUP}"
    local R="${RES_DIR}/${GROUP}/admixture"; mkdir -p "${R}"
    log "=== ${GROUP} ==="

    # A group whose input was never built must not throw away another group's
    # finished work. This used to be require_file -> die, so adding groupA to
    # ADMIX_GROUPS without rebuilding 03 killed the stage at the very end,
    # after groupB had run for six hours.
    if [[ ! -s "${G}/merged_admix.bed" ]]; then
        log "  SKIPPED: ${G}/merged_admix.bed does not exist."
        log "  03 builds merged_admix only for the groups in ADMIX_GROUPS, so add"
        log "  this group there and rerun 03_relatedness.sh before 05."
        return 0
    fi
    local FULL="${G}/merged_admix"
    local STEM_FULL; STEM_FULL=$(basename "${FULL}")

    local TRAIN="${FULL}"
    if [[ "${ADMIX_PROJECT}" -eq 1 ]]; then
        require_file "${G}/merged_admix_train.bed"
        TRAIN="${G}/merged_admix_train"
        cmp -s "${FULL}.bim" "${TRAIN}.bim" \
            || die "merged_admix.bim and merged_admix_train.bim differ; -P projection would be invalid"
    fi
    local STEM_TRAIN; STEM_TRAIN=$(basename "${TRAIN}")

    local NSF NST NV
    NSF=$(nlines "${FULL}.fam"); NST=$(nlines "${TRAIN}.fam"); NV=$(nlines "${FULL}.bim")
    log "  full panel:   ${NSF} samples, ${NV} pruned SNPs"
    log "  training set: ${NST} samples (references <= ${PCA_MAX_PER_NATION}, cohort <= ${COHORT_CAP})"
    [[ "${NV}" -ge "${ADMIX_MIN_SNPS}" ]] \
        || log "  WARNING: ${NV} SNPs is below ADMIX_MIN_SNPS=${ADMIX_MIN_SNPS}"
    log "  $(( (K_MAX - K_MIN + 1) * N_REPS )) training runs, ${ADMIX_JOBS} concurrent x ${ADMIX_THREADS} threads"

    # individual order for the FULL panel, needed to align .Q rows with IDs
    awk '{print $2}' "${FULL}.fam" > "${R}/individuals.txt"
    # Nothing downstream may overwrite this: the projected .Q is one row per
    # individual of the FULL panel, in this order.
    [[ "$(nlines "${R}/individuals.txt")" -eq "${NSF}" ]] \
        || die "individuals.txt does not match merged_admix.fam"
    awk -F'\t' 'NR==FNR{pop[$1]=$3; next} {print ($1 in pop)?pop[$1]:"NA"}' \
        <(tail -n +2 "${META_DIR}/pop_master.tsv") "${R}/individuals.txt" \
        > "${R}/ind2pop.txt"
    [[ "$(grep -c '^NA$' "${R}/ind2pop.txt" || true)" -eq 0 ]] \
        || die "some individuals in merged_admix have no population label"
    awk -F'\t' 'NR==FNR{nat[$1]=$2; next} {print ($1 in nat)?nat[$1]:"NA"}' \
        <(tail -n +2 "${META_DIR}/pop_master.tsv") "${R}/individuals.txt" \
        > "${R}/ind2nation.txt"

    # --- 1. training sweep --------------------------------------------------
    log "[1/4] training runs on $(basename "${TRAIN}")"
    launch_runs "${TRAIN}" "${R}" "${K_MIN}" "${K_MAX}" "${N_REPS}" "${SEED_BASE}"
    log "  CV summary (training set):"; sed 's/^/    /' "${R}/cv_summary.tsv"

    local BEST_K
    BEST_K=$(awk -F'\t' 'NR>1{if(m==""||$3<m){m=$3;k=$1}} END{print k}' "${R}/cv_summary.tsv")
    log "  lowest mean CV error at K=${BEST_K}"
    # Header row, like every other table the pipeline writes. Without one a
    # generic TSV reader takes the single line as the header and finds no data.
    printf 'key\tvalue\nbest_k\t%s\n' "${BEST_K}" > "${R}/best_k.tsv"

    # --- 2. alignment and modes, per K -------------------------------------
    # The .Q rows of the TRAINING runs are in the training .fam order, which is
    # not the full panel's order. Two individual lists, never interchangeable.
    log "[2/4] replicate alignment and mode detection"
    awk '{print $2}' "${TRAIN}.fam" > "${R}/individuals_train.txt"
    ${PYTHON} "${ALIGN}" sweep --admix-dir "${R}" \
        --k-min "${K_MIN}" --k-max "${K_MAX}" --threshold "${ADMIX_SIM_THRESHOLD}" \
        --individuals "${R}/individuals_train.txt" \
        --meta "${META_DIR}/pop_master.tsv" --out-prefix "${R}/train"
    sed 's/^/    /' "${R}/train.modes_by_K.tsv"

    # --- 3. projection of the full panel ------------------------------------
    if [[ "${ADMIX_PROJECT}" -eq 1 ]]; then
        log "[3/4] projecting the full panel with admixture -P"
        mkdir -p "${R}/proj"
        local K
        for K in $(seq "${K_MIN}" "${K_MAX}"); do
            local SUM="${R}/train.K${K}.summary.tsv"
            [[ -s "${SUM}" ]] || { log "  K=${K}: no summary, skipped"; continue; }
            local PSRC; PSRC=$(awk -F'\t' '$1=="representative_P"{print $2}' "${SUM}")
            [[ -s "${PSRC}" ]] || { log "  K=${K}: ${PSRC} missing, skipped"; continue; }

            local PD="${R}/proj/K${K}"; mkdir -p "${PD}"
            ln -sf "${FULL}.bed" "${PD}/${STEM_FULL}.bed"
            ln -sf "${FULL}.bim" "${PD}/${STEM_FULL}.bim"
            ln -sf "${FULL}.fam" "${PD}/${STEM_FULL}.fam"
            cp "${PSRC}" "${PD}/${STEM_FULL}.${K}.P.in"
            # -j and backgrounded. The projection was single-threaded and
            # sequential across K, which on a 64-core box left 63 cores idle
            # for hours while nine fixed-P fits ran one at a time.
            while [[ "$(jobs -rp | wc -l)" -ge "${ADMIX_JOBS}" ]]; do wait -n || true; done
            log "  K=${K}: projecting"
            ( cd "${PD}" && ${ADMIXTURE} -P -j"${ADMIX_THREADS}" \
                "${STEM_FULL}.bed" "${K}" > log.txt 2>&1 ) &
        done
        wait || true

        for K in $(seq "${K_MIN}" "${K_MAX}"); do
            local PD="${R}/proj/K${K}"
            [[ -f "${PD}/${STEM_FULL}.${K}.Q" ]] \
                || { log "  K=${K}: no projected .Q produced, see ${PD}/log.txt"; continue; }
            [[ "$(nlines "${R}/individuals.txt")" -eq "${NSF}" ]] \
                || die "individuals.txt was overwritten (${NSF} expected); the \
projected Q cannot be labelled"
            ${PYTHON} "${ALIGN}" q-table --q "${PD}/${STEM_FULL}.${K}.Q" \
                --individuals "${R}/individuals.txt" --meta "${META_DIR}/pop_master.tsv" \
                --out "${R}/Q_projected.K${K}.tsv"
            log "  K=${K}: projected Q for ${NSF} individuals -> Q_projected.K${K}.tsv"
        done
    else
        log "[3/4] ADMIX_PROJECT=0, no projection; the training sweep IS the result"
        local K
        for K in $(seq "${K_MIN}" "${K_MAX}"); do
            [[ -s "${R}/train.K${K}.Q_major.tsv" ]] && \
                cp "${R}/train.K${K}.Q_major.tsv" "${R}/Q_projected.K${K}.tsv"
        done
    fi

    # --- 4. joint (uncapped) sensitivity at the selected K ------------------
    if [[ "${ADMIX_JOINT_SENSITIVITY}" -eq 1 && "${ADMIX_PROJECT}" -eq 1 ]]; then
        log "[4/4] joint uncapped fit at K=${BEST_K} (sensitivity)"
        local J="${R}/joint"
        launch_runs "${FULL}" "${J}" "${BEST_K}" "${BEST_K}" "${N_REPS}" "$(( SEED_BASE + 50000 ))"
        ${PYTHON} "${ALIGN}" summarise --admix-dir "${J}" --k "${BEST_K}" \
            --threshold "${ADMIX_SIM_THRESHOLD}" \
            --individuals "${J}/individuals_run.txt" --meta "${META_DIR}/pop_master.tsv" \
            --out-prefix "${J}/joint"
        if [[ -s "${R}/Q_projected.K${BEST_K}.tsv" && -s "${J}/joint.Q_major.tsv" ]]; then
            ${PYTHON} "${ALIGN}" compare-q \
                --a "${R}/Q_projected.K${BEST_K}.tsv" --b "${J}/joint.Q_major.tsv" \
                --out "${R}/projected_vs_joint.K${BEST_K}.tsv"
            sed 's/^/    /' "${R}/projected_vs_joint.K${BEST_K}.tsv"
            log "  A per-component correlation near 1 means capping the cohort and"
            log "  projecting reproduces the joint fit: the sample-size objection is"
            log "  answered with a number instead of an argument."
        fi
    else
        log "[4/4] joint sensitivity: disabled"
    fi
}

for GROUP in "${ADMIX_GROUPS[@]}"; do run_group "${GROUP}"; done

log "=== 05 done -> 06_fst.sh ==="
log "Figures come from Q_projected.K<k>.tsv; mode structure from train.modes_by_K.tsv."
log "Report K, the CV mean +/- SD, and the major-mode fraction together."
