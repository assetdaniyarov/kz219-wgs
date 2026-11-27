#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 06_fst.sh — pairwise Fst, Hudson estimator, block-jackknife SE.
#
# Estimator: Hudson, ratio of averages (Bhatia et al. 2013, Genome Res).
# Weir-Cockerham is biased under very unequal sample sizes, and ours are.
#
# Input is merged_fst: unrelated, de-duplicated, missingness-filtered, the full
# uncapped sample set, and deliberately NOT MAF-filtered. A MAF threshold on a
# merged multi-population panel is an ascertainment step that shifts every
# pairwise value, because the sites it removes are not random with respect to
# the populations being compared.
#
# NEGATIVE VALUES ARE REPORTED AS THEY COME. Hudson's estimator is unbiased, so
# a pair with no real divergence produces small negative values about half the
# time. Truncating them at zero makes every such estimate look positive and
# turns a null result into a weak positive one. fst_pairs.tsv carries the raw
# value, the jackknife SE, the 95% interval and a flag for whether that interval
# covers zero, which is the quantity a reader actually needs.
#
# WHICH PANEL CARRIES WHICH CLAIM.
#   groupA is WGS on both sides and is the panel the superpopulation-level Fst
#   should rest on. At NATION level it is nearly useless for this paper: almost
#   every SGDP population is n<10 and is excluded, so what survives is the 1000G
#   populations, and the cohort's actual neighbours are not in it.
#   groupB has those neighbours — Uygur, Hazara, Yakut, Mongola, the Pakistani
#   set — at n=9-51, and is therefore where the nation-level comparisons live.
#   It also mixes 30x WGS cohort genotypes with an array reference, and
#   genotyping differences act in the same direction as real divergence.
#
# So the platform effect is MEASURED. Both panels are scored at superpopulation
# level over the same eight labels, and fst_platform_concordance.tsv reports the
# two values per pair, their difference, and the mean difference split by
# whether the pair involves the cohort. If the array/WGS mixture were inflating
# the cohort's distances specifically, the KAZ pairs would show an offset the
# reference-vs-reference pairs do not. That is a check a reviewer can repeat.
#
# MINIMUM n. Groups below FST_MIN_N_* are removed from the phenotype file rather
# than computed and flagged afterwards. Hudson is unbiased down to n=2, but the
# jackknife interval at n=2-3 is wider than the entire range of the table, and a
# reported point estimate carrying that interval is read as a measurement. The
# excluded groups and their n are written out, so the exclusion is auditable.
#
# report-variants writes one row per variant per pair, tens of gigabytes. Off by
# default; set FST_REPORT_VARIANTS=1 only when a per-locus scan is needed.
# ============================================================================

need "${PLINK2}"
TABLES="${SCRIPT_DIR}/popgen_tables.py"
require_file "${TABLES}"

run_level() {
    local GROUP=$1 LEVEL=$2 COL=$3 MINN=$4
    local G="${WORK_DIR}/${GROUP}"
    local R="${RES_DIR}/${GROUP}/fst/${LEVEL}${FST_LEVEL_SUFFIX}"; mkdir -p "${R}"
    log "--- ${GROUP} / ${LEVEL}${FST_LEVEL_SUFFIX}  (min n = ${MINN}) ---"

    require_file "${G}/merged_fst.bed"
    require_file "${META_DIR}/pop_master.tsv"

    # Phenotype file built FROM THE FILESET'S OWN .fam, carrying FID and IID
    # exactly as they appear there, not from the pop_master IID list.
    #
    # Written the other way round, groupB loaded "1 categorical phenotype
    # (218 values)": only the Kazakh individuals matched and all HGDP
    # individuals came through unlabelled, leaving one population and aborting
    # --fst. The same IIDs matched fine in the awk that builds pop_n.tsv below,
    # so the mismatch was in how plink2 resolves an IID-only phenotype file
    # against a fileset whose FIDs are not blank. Emitting FID and IID straight
    # from the .fam removes the question: there is nothing left to resolve.
    awk -v c="${COL}" '
        NR==FNR { if (FNR>1) { n=split($0,f,"\t"); if (f[c]!="" && f[c]!="NA") p[f[1]]=f[c] }; next }
        ($2 in p) { print $1"\t"$2"\t"p[$2] }' \
        "${META_DIR}/pop_master.tsv" "${G}/merged_fst.fam" > "${R}/pop_all.body"
    { printf '#FID\tIID\tPOP\n'; cat "${R}/pop_all.body"; } > "${R}/pop_all.cat"
    rm -f "${R}/pop_all.body"

    local NFAM NLAB
    NFAM=$(nlines "${G}/merged_fst.fam"); NLAB=$(nbody "${R}/pop_all.cat")
    log "  labelled: ${NLAB} of ${NFAM} individuals in the fileset"
    [[ "${NLAB}" -eq "${NFAM}" ]] || \
        log "  WARNING: $(( NFAM - NLAB )) individuals carry no ${LEVEL} label"

    # sample sizes, counted on the file plink2 will actually read
    awk -F'\t' 'NR>1{n[$3]++} END{for(p in n) printf "%s\t%d\n", p, n[p]}' \
        "${R}/pop_all.cat" | sort -k2,2nr > "${R}/pop_n.body"
    { printf 'POP\tN\n'; cat "${R}/pop_n.body"; } > "${R}/pop_n.tsv"
    rm -f "${R}/pop_n.body"

    awk -F'\t' -v m="${MINN}" 'BEGIN{OFS="\t"; print "POP","N"}
        NR>1 && $2 < m {print $1,$2}' "${R}/pop_n.tsv" > "${R}/pop_excluded.tsv"
    local NKEEP NEXC
    NEXC=$(nbody "${R}/pop_excluded.tsv")
    NKEEP=$(( $(nbody "${R}/pop_n.tsv") - NEXC ))
    log "  groups: ${NKEEP} retained, ${NEXC} below n=${MINN} and excluded"
    [[ "${NKEEP}" -ge 2 ]] || die "${GROUP}/${LEVEL}: only ${NKEEP} groups clear n=${MINN}"

    # phenotype file restricted to the retained groups
    awk -F'\t' -v m="${MINN}" '
        NR==FNR { if (FNR>1 && $2>=m) ok[$1]=1; next }
        FNR==1  { print; next }
        ($3 in ok)' "${R}/pop_n.tsv" "${R}/pop_all.cat" > "${R}/pop.cat"

    # Levels actually present in the file plink2 will read. The count above
    # comes from pop_n.tsv; this one comes from pop.cat itself, so a mismatch
    # between the two is a bug in the filter rather than a property of the data.
    local NLEV
    NLEV=$(awk -F'\t' 'NR>1{v[$3]=1} END{print length(v)}' "${R}/pop.cat")
    log "  levels in pop.cat handed to plink2: ${NLEV}"
    awk -F'\t' 'NR>1{n[$3]++} END{for(k in n) printf "    %-14s %d\n", k, n[k]}' \
        "${R}/pop.cat" | sort -k2,2nr
    [[ "${NLEV}" -ge 2 ]] || die "${GROUP}/${LEVEL}: pop.cat carries ${NLEV} level(s)"

    # --- pairwise Fst -------------------------------------------------------
    local RV=()
    [[ "${FST_REPORT_VARIANTS}" -eq 1 ]] && RV=(report-variants)
    # --no-psam-pheno: the .fam carries a phenotype column inherited from the
    # source panels, and it holds a single category (0 cases, N controls).
    # Loaded beside ours, plink2 reports "2 phenotypes loaded" and --fst can
    # land on that one instead of POP, failing with "requires two or more
    # populations" even though pop.cat has eight.
    ${PLINK2} --bfile "${G}/merged_fst" --no-psam-pheno \
        --pheno "${R}/pop.cat" --fst POP method=hudson \
        blocksize="${FST_BLOCKSIZE}" "${RV[@]+"${RV[@]}"}" \
        --threads "${THREADS}" --out "${R}/fst"
    require_file "${R}/fst.fst.summary"

    # --- publication tables -------------------------------------------------
    ${PYTHON} "${TABLES}" fst-tables --summary "${R}/fst.fst.summary" \
        --pop-n "${R}/pop_n.tsv" \
        --out-pairs "${R}/fst_pairs.tsv" --out-matrix "${R}/fst_matrix.tsv"
    log "  pairs whose 95% interval covers zero: \
$(awk -F'\t' 'NR>1 && $9=="yes"' "${R}/fst_pairs.tsv" | wc -l) of $(nbody "${R}/fst_pairs.tsv")"
    head -4 "${R}/fst_pairs.tsv" | column -t -s $'\t' | sed 's/^/    /'

    # --- the cohort against everything, ranked -----------------------------
    # 595 pairs at groupB/nation is the whole table; this is the one row of it
    # the main text is actually about. SEPARATED_FROM_PREVIOUS says whether
    # consecutive neighbours can be told apart at all, which is what turns
    # "closest to Uygur, then Hazara" from an ordering into a claim.
    # The cohort's label DEPENDS ON THE LEVEL: "KAZ" at superpopulation level,
    # "Kazakh" at nation level. Hardcoding KAZ meant both nation-level runs -
    # the 18-population groupA table and the 35-population groupB one, which is
    # where the cohort's actual neighbours live - silently produced nothing.
    # Derived from pop.cat instead: whichever label the KAZ individuals carry.
    local FOCUS
    FOCUS=$(awk -F'\t' '
        NR==FNR { if (FNR>1 && $3=="KAZ") kz[$1]=1; next }
        FNR>1 && ($2 in kz) { n[$3]++ }
        END { b=""; m=0; for (k in n) if (n[k]>m) { m=n[k]; b=k }; print b }' \
        "${META_DIR}/pop_master.tsv" "${R}/pop.cat")
    if [[ -z "${FOCUS}" ]]; then
        log "  no cohort label survives at this level; ranking and focus skipped"
    else
        log "  cohort label at ${LEVEL} level: ${FOCUS}"
        ${PYTHON} "${TABLES}" fst-focus --summary "${R}/fst.fst.summary" \
            --pop-n "${R}/pop_n.tsv" --focus "${FOCUS}" \
            --out "${R}/fst_vs_cohort.tsv"
        head -8 "${R}/fst_vs_cohort.tsv" | column -t -s $'\t' | sed 's/^/    /'
    fi

    # --- heatmap ------------------------------------------------------------
    # Rows and columns are seriated by average-linkage clustering with optimal
    # leaf ordering, not left alphabetical: an alphabetical 35x35 Fst matrix
    # shows nothing. The PNG needs matplotlib; without it the ordered matrix is
    # still written and can be plotted anywhere.
    ${PYTHON} "${TABLES}" fst-heatmap --matrix "${R}/fst_matrix.tsv" \
        --out-ordered "${R}/fst_matrix_ordered.tsv" \
        --out-png "${R}/fst_heatmap.png" --focus "${FOCUS:-KAZ}" \
        --title "${GROUP} / ${LEVEL}: pairwise Hudson Fst" \
      || log "  heatmap step failed (non-fatal)"

    # --- per-population allele counts, for any downstream f-statistic work ---
    mkdir -p "${R}/counts"
    ${PLINK2} --bfile "${G}/merged_fst" --no-psam-pheno \
        --pheno "${R}/pop.cat" --loop-cats POP --freq counts \
        --threads "${THREADS}" --out "${R}/counts/pop"
    log "  wrote $(ls "${R}"/counts/*.acount 2>/dev/null | wc -l) count files"
}

for GROUP in "${FST_GROUPS[@]}"; do
    run_level "${GROUP}" population 3 "${FST_MIN_N_POP}"     # col 3 = superpopulation
    run_level "${GROUP}" nation     2 "${FST_MIN_N_NATION}"  # col 2 = Nation
done

# ---------------------------------------------------------------------------
# platform concordance
# ---------------------------------------------------------------------------
A_SUM="${RES_DIR}/groupA_kz_1000g_sgdp/fst/population/fst.fst.summary"
B_SUM="${RES_DIR}/groupB_kz_hgdp/fst/population/fst.fst.summary"
if [[ "${FST_PLATFORM_CHECK}" -eq 1 && -s "${A_SUM}" && -s "${B_SUM}" ]]; then
    log "=== platform concordance: groupA (all WGS) vs groupB (WGS cohort + array reference) ==="
    OUT="${RES_DIR}/fst_platform_concordance"; mkdir -p "${OUT}"
    ${PYTHON} "${TABLES}" fst-concordance --a "${A_SUM}" --b "${B_SUM}" \
        --label-a groupA_wgs --label-b groupB_array --focus KAZ \
        --out "${OUT}/pairs.tsv" --out-stats "${OUT}/stats.tsv"
    column -t -s $'\t' "${OUT}/pairs.tsv" | sed 's/^/    /'
    column -t -s $'\t' "${OUT}/stats.tsv" | sed 's/^/    /'
    log "  Read mean_diff_focus against mean_diff_other. A cohort-specific offset"
    log "  is the signature of a platform effect; a shared offset is the panels"
    log "  differing in marker ascertainment, which affects every pair equally."
else
    log "=== platform concordance: skipped ==="
fi

log "=== 06 done -> 07_roh.sh ==="
log "Main text: groupA at superpopulation level (WGS throughout) and groupB at"
log "nation level (the only panel carrying the cohort's neighbours), with the"
log "concordance table cited for the platform caveat."
