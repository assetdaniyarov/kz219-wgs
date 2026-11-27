#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 03_relatedness.sh — KING, reference balancing, final analysis sets.
#
#   merged_raw
#     -> split diploid / pseudo-haploid
#     -> pruning of the diploid part for KING (KING_PRUNE_*, ~100k+ markers)
#     -> KING: cryptic duplicates + relatives
#     -> analysis.keep        unrelated diploid                 -> Fst, ROH
#     -> pca.keep             analysis.keep, reference nations capped
#     -> pca_capped.keep      pca.keep with the COHORT capped too
#     -> merged_pca           MAF+geno+highLD+pruned on pca.keep -> PCA
#     -> merged_pca_capped    merged_pca subset to pca_capped.keep
#     -> merged_admix         merged_pca minus non-Eurasians, re-MAF, re-pruned
#     -> merged_admix_train   merged_admix subset to pca_capped.keep
#
# The order is the point. MAF and geno decide which variants survive, so they
# are computed after the sample set is final. An individual present in two
# reference panels otherwise contributes twice to every allele frequency that
# drives that decision.
#
# KING GETS ITS OWN PRUNING. Kinship precision is set by the number of
# independent markers, and KING wants ~100k. Sharing the PCA threshold left
# groupB with 55,710 markers under the old PRUNE_R2=0.1. KING_PRUNE_R2 is
# separate and looser, and the count is checked rather than assumed.
#
# CUTOFF. KING_CUTOFF=0.0884 removes 2nd degree and closer. The bands are
# KING's own: >0.354 duplicate/MZ, 0.177-0.354 1st, 0.0884-0.177 2nd,
# 0.0442-0.0884 3rd. The previous value, 0.13, sits above the 0.125 expected for
# a genuine 2nd-degree pair and therefore keeps roughly half of them; it was
# chosen so that all 219 Kazakhs survived, and it is not worth defending in
# review. Keep-lists at the alternative cutoffs are written for the Methods
# sensitivity statement.
#
# PSEUDO-HAPLOID SPLIT. Individuals listed in meta/pseudohaploid_iids.txt carry
# one randomly drawn allele per site written as a homozygote. They are removed
# before KING, because kinship is estimated from heterozygote sharing and there
# are none; before the PC computation and ADMIXTURE, because zero within-
# individual heterozygosity is not a property any of those models can absorb;
# and before Fst and ROH, where it corrupts the estimator outright. They survive
# only as projected points, and 08_project_aadr.sh builds that panel from the
# source AADR panel directly, so no analysis fileset here ever contains them.
#
# CAPPING, AND WHAT CHANGED. PCA axes and ADMIXTURE components are driven by the
# populations with the most individuals (McVean 2009). Reference nations are
# capped at PCA_MAX_PER_NATION. The cohort is not capped in the primary set —
# it is the subject of the paper — but that leaves 218 Kazakh against reference
# nations of at most 25, i.e. the cohort is a quarter of the ADMIXTURE panel and
# ~9x the largest reference group. That is the same imbalance the cap exists to
# prevent, pointing the other way, and it manufactures a cohort-specific
# component out of sample size. So a cohort-capped set is built here as well:
# 04 recomputes the PCA on it and correlates the axes, and 05 trains ADMIXTURE
# on it and projects the full cohort with -P. Selection is deterministic —
# lowest genotype missingness, IID as the tie-break.
#
# EURASIAN SET FOR ADMIXTURE. merged_admix drops the African outgroup and
# Oceania. Left in, K=2 splits Africa from everything else, K=3 adds Oceania and
# K=4 splits East from West Eurasia: three components spent before anything
# Central Asian appears. The outgroup stays in merged_fst, where it anchors the
# Fst scale, and nowhere else.
# ============================================================================

need "${PLINK2}"

HIGHLD="${WORK_DIR}/highld_hg19.txt"
require_file "${HIGHLD}"

# apply_cap <cap_input.tsv> <ref cap> <cohort cap> <fam> <out keep>
# cap_input columns: IID  Nation  Population  F_MISS  KEY
# A cap of 0 means "no cap for that side".
apply_cap() {
    local IN=$1 CAP=$2 CCAP=$3 FAM=$4 OUT=$5
    sort -t$'\t' -k5,5 -k4,4g -k1,1 "${IN}" \
      | awk -F'\t' -v cap="${CAP}" -v ccap="${CCAP}" '
            { lim = ($3=="KAZ") ? ccap : cap
              if (lim+0 <= 0) { print $1; next }
              if (++c[$5] <= lim+0) print $1 }' > "${OUT}.ids"
    mk_id_list "${FAM}" "${OUT}" "${OUT}.ids"
    rm -f "${OUT}.ids"
}

run_group() {
    local GROUP=$1
    local G="${WORK_DIR}/${GROUP}"
    local R="${RES_DIR}/${GROUP}"; mkdir -p "${R}"
    log "=== ${GROUP} ==="

    require_file "${G}/merged_raw.bed"
    require_file "${G}/merged_raw.smiss"

    local WANT_ADMIX=0 WANT_FSTROH=0
    in_list "${GROUP}" "${ADMIX_GROUPS[@]}" && WANT_ADMIX=1
    in_list "${GROUP}" "${FST_GROUPS[@]}"   && WANT_FSTROH=1
    [[ "${GROUP}" == "${ROH_GROUP}" ]]      && WANT_FSTROH=1

    # --- 0. diploid / pseudo-haploid split ---------------------------------
    : > "${G}/pseudohaploid.id"
    if [[ -s "${MF_PSEUDOHAPLOID}" ]]; then
        mk_id_list "${G}/merged_raw.fam" "${G}/pseudohaploid.id" "${MF_PSEUDOHAPLOID}"
    fi
    local NPH; NPH=$(nlines "${G}/pseudohaploid.id")
    log "[0/8] pseudo-haploid individuals present: ${NPH}"

    local DIP=("--bfile" "${G}/merged_raw")
    [[ "${NPH}" -gt 0 ]] && DIP+=("--remove" "${G}/pseudohaploid.id")

    # --- 1. pruning for KING ------------------------------------------------
    log "[1/8] pruning for KING (diploid only, r2=${KING_PRUNE_R2}, ${KING_PRUNE_WINDOW})"
    ${PLINK2} "${DIP[@]}" --maf 0.05 --geno 0.05 \
        --exclude range "${HIGHLD}" \
        --indep-pairwise "${KING_PRUNE_WINDOW}" "${KING_PRUNE_STEP}" "${KING_PRUNE_R2}" \
        --threads "${THREADS}" --out "${G}/king_prune"
    ${PLINK2} "${DIP[@]}" --extract "${G}/king_prune.prune.in" \
        --make-bed --threads "${THREADS}" --out "${G}/king_in"
    local NKM; NKM=$(nlines "${G}/king_in.bim")
    log "  KING markers: ${NKM}  samples: $(nlines "${G}/king_in.fam")"
    [[ "${NKM}" -ge "${KING_MIN_MARKERS}" ]] || \
        log "  WARNING: ${NKM} markers is below KING_MIN_MARKERS=${KING_MIN_MARKERS}; \
kinship estimates will be noisy. Raise KING_PRUNE_R2 or widen the intersection."

    # --- 2. pairwise kinship ------------------------------------------------
    # Reported down to KING_REPORT (4th degree) so the whole tail is visible. If
    # drift in an endogamous cohort were inflating kinship, the low bins would be
    # flooded; an empty tail beside a clean 1st-degree mode is the evidence that
    # the detected relationships are real. Note that this holds because the panel
    # is multi-population: KING-robust is centred on the sample as a whole. The
    # same bins applied inside a single endogamous population are NOT meaningful,
    # which is why 09_kz_tribes.sh does not use them.
    log "[2/8] KING pairwise (report >= ${KING_REPORT})"
    ${PLINK2} --bfile "${G}/king_in" --make-king-table \
        --king-table-filter "${KING_REPORT}" \
        --threads "${THREADS}" --out "${R}/king_pairs"
    log "  pairs reported: $(nbody "${R}/king_pairs.kin0")"

    awk -v DUP="${KING_DUP}" 'BEGIN{OFS="\t"}
        NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }
                if (!("IID1" in h) || !("KINSHIP" in h)) {
                    print "unexpected .kin0 header" > "/dev/stderr"; exit 2 }
                print "IID1","IID2","KINSHIP","DEGREE"; next }
        { k=$(h["KINSHIP"]); a=$(h["IID1"]); b=$(h["IID2"]);
          if      (k >  DUP)    d="duplicate/MZ";
          else if (k >  0.177)  d="1st degree";
          else if (k >  0.0884) d="2nd degree";
          else if (k >  0.0442) d="3rd degree";
          else                  d="4th degree";
          print a,b,k,d }' \
        "${R}/king_pairs.kin0" > "${R}/related.tsv"

    awk -F'\t' 'NR>1 && $4=="duplicate/MZ"' "${R}/related.tsv" > "${R}/dup_pairs.tsv"
    log "  cryptic duplicate pairs: $(nlines "${R}/dup_pairs.tsv")"
    awk -F'\t' 'NR>1{c[$4]++} END{for(d in c) printf "  %-14s %d\n", d, c[d]}' \
        "${R}/related.tsv"

    # Kinship histogram for the supplement. Replaces the notebook figure.
    { printf 'bin_lo\tbin_hi\tn_pairs\n'
      awk -F'\t' 'NR>1 { b=int($3/0.005); c[b]++ }
          END { for (b in c) printf "%.3f\t%.3f\t%d\n", b*0.005, (b+1)*0.005, c[b] }' \
          "${R}/related.tsv" | sort -g
    } > "${R}/kinship_hist.tsv"

    # --- 3. unrelated set ---------------------------------------------------
    log "[3/8] KING cutoff ${KING_CUTOFF} (removes 2nd degree and closer)"
    ${PLINK2} --bfile "${G}/king_in" --king-cutoff "${KING_CUTOFF}" \
        --threads "${THREADS}" --out "${G}/rel"
    require_file "${G}/rel.king.cutoff.in.id"
    norm_id_list "${G}/rel.king.cutoff.in.id"  "${G}/merged_raw.fam" "${R}/analysis.keep"
    norm_id_list "${G}/rel.king.cutoff.out.id" "${G}/merged_raw.fam" "${R}/analysis.removed"
    log "  retained $(nlines "${R}/analysis.keep") / removed $(nlines "${R}/analysis.removed")"

    kaz_n() {   # kaz_n <keep file>
        awk -F'\t' 'NR==FNR{k[$2]=1; next} FNR>1 && ($1 in k) && $3=="KAZ"' \
            "$1" "${META_DIR}/pop_master.tsv" | wc -l
    }
    local N_KAZ; N_KAZ=$(kaz_n "${R}/analysis.keep")
    log "  Kazakh individuals in the unrelated set: ${N_KAZ}"

    { printf 'IID1\tIID2\tKINSHIP\tDEGREE\n'
      awk -F'\t' -v C="${KING_CUTOFF}" 'NR>1 && $3 < C' "${R}/related.tsv"; } \
        > "${R}/related_but_retained.tsv"
    log "  reported pairs below the removal cutoff (retained): $(nbody "${R}/related_but_retained.tsv")"

    # --- 3b. cutoff sensitivity, keep-lists only ----------------------------
    # Cheap, and it turns "we used 0.0884" into "here is what every other
    # plausible threshold would have given".
    { printf 'cutoff\tn_retained\tn_removed\tn_kazakh\n'
      printf '%s\t%s\t%s\t%s\n' "${KING_CUTOFF}" \
             "$(nlines "${R}/analysis.keep")" "$(nlines "${R}/analysis.removed")" "${N_KAZ}"
      for C in ${KING_CUTOFF_SENSITIVITY}; do
          ${PLINK2} --bfile "${G}/king_in" --king-cutoff "${C}" \
              --threads "${THREADS}" --out "${G}/rel_${C}" >/dev/null 2>&1 || continue
          norm_id_list "${G}/rel_${C}.king.cutoff.in.id" "${G}/merged_raw.fam" \
                       "${R}/analysis_keep_cutoff_${C}.id"
          printf '%s\t%s\t%s\t%s\n' "${C}" \
                 "$(nlines "${R}/analysis_keep_cutoff_${C}.id")" \
                 "$(( $(nlines "${G}/king_in.fam") - $(nlines "${R}/analysis_keep_cutoff_${C}.id") ))" \
                 "$(kaz_n "${R}/analysis_keep_cutoff_${C}.id")"
      done
    } | sort -g -k1,1 > "${R}/king_cutoff_sensitivity.tsv"
    log "  cutoff sensitivity:"
    sed 's/^/    /' "${R}/king_cutoff_sensitivity.tsv"

    # --- 4. reference balancing --------------------------------------------
    log "[4/8] capping reference nations at ${PCA_MAX_PER_NATION}"
    awk -F'\t' -v OFS='\t' '
        FNR==NR { keep[$2]=1; next }
        FILENAME ~ /pop_master\.tsv$/ {
            if (FNR==1) next
            if ($1 in keep) { nat[$1]=$2; pop[$1]=$3 }
            next }
        {   if (FNR==1) { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; next }
            id=$(h["IID"])
            if (id in nat) print id, nat[id], pop[id], $(h["F_MISS"]),
                                 (pop[id]=="KAZ" ? "KAZ" : nat[id]) }' \
        "${R}/analysis.keep" "${META_DIR}/pop_master.tsv" "${G}/merged_raw.smiss" \
        > "${G}/cap_input.tsv"

    # primary set: references capped, cohort uncapped
    apply_cap "${G}/cap_input.tsv" "${PCA_MAX_PER_NATION}" 0 \
              "${G}/merged_raw.fam" "${R}/pca.keep"
    # sensitivity / ADMIXTURE training set: cohort capped too
    apply_cap "${G}/cap_input.tsv" "${PCA_MAX_PER_NATION}" "${COHORT_CAP}" \
              "${G}/merged_raw.fam" "${R}/pca_capped.keep"

    { printf 'Nation\tPopulation\tn_unrelated\tn_after_cap\tn_after_cohort_cap\n'
      awk -F'\t' -v OFS='\t' '
          FNR==NR                    { sel[$2]=1; next }
          FILENAME ~ /pca_capped/    { sel2[$2]=1; next }
          { tot[$2"\t"$3]++
            if ($1 in sel)  got[$2"\t"$3]++
            if ($1 in sel2) got2[$2"\t"$3]++ }
          END { for (k in tot) print k, tot[k], (k in got?got[k]:0), (k in got2?got2[k]:0) }' \
          "${R}/pca.keep" "${R}/pca_capped.keep" "${G}/cap_input.tsv" \
        | sort -t$'\t' -k2,2 -k1,1
    } > "${R}/pca_cap.tsv"
    log "  PCA set:          $(nlines "${R}/pca.keep") of $(nlines "${R}/analysis.keep") unrelated (cohort uncapped)"
    log "  cohort-capped set: $(nlines "${R}/pca_capped.keep") (KAZ capped at ${COHORT_CAP})"

    # --- 5. PCA set ---------------------------------------------------------
    log "[5/8] PCA set (prune r2=${PRUNE_R2}, ${PRUNE_WINDOW})"
    ${PLINK2} --bfile "${G}/merged_raw" --keep "${R}/pca.keep" \
        --maf "${MAF}" --geno "${GENO}" --exclude range "${HIGHLD}" \
        --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_R2}" \
        --threads "${THREADS}" --out "${G}/pca_prune"
    ${PLINK2} --bfile "${G}/merged_raw" --keep "${R}/pca.keep" \
        --extract "${G}/pca_prune.prune.in" \
        --make-bed --threads "${THREADS}" --out "${G}/merged_pca"
    log "  merged_pca: $(nlines "${G}/merged_pca.fam") samples, $(nlines "${G}/merged_pca.bim") SNPs"

    # Same variants, cohort capped. Kept as a --keep subset rather than
    # recomputed so that the two PCAs rest on an identical marker set and the
    # axis correlation in 04 measures the sample composition and nothing else.
    ${PLINK2} --bfile "${G}/merged_pca" --keep "${R}/pca_capped.keep" \
        --make-bed --threads "${THREADS}" --out "${G}/merged_pca_capped"
    log "  merged_pca_capped: $(nlines "${G}/merged_pca_capped.fam") samples"

    if [[ "${PCA_ALL_SAMPLES}" -eq 1 ]]; then
        ${PLINK2} "${DIP[@]}" --extract "${G}/pca_prune.prune.in" \
            --make-bed --threads "${THREADS}" --out "${G}/merged_pca_all"
    fi

    # --- 6. ADMIXTURE set ---------------------------------------------------
    if [[ "${WANT_ADMIX}" -eq 1 ]]; then
        log "[6/8] ADMIXTURE set (Eurasia only)"
        awk -F'\t' -v drop="${DROP_SUPERPOPS} ${EURASIA_EXCLUDE}" '
            BEGIN { n=split(drop,a," "); for(i=1;i<=n;i++) D[a[i]]=1 }
            NR>1 && ($3 in D) { print $1 }' "${META_DIR}/pop_master.tsv" \
          > "${G}/non_eurasian.ids"
        mk_id_list "${G}/merged_pca.fam" "${G}/non_eurasian.id" "${G}/non_eurasian.ids"
        rm -f "${G}/non_eurasian.ids"
        log "  removed from ADMIXTURE: $(nlines "${G}/non_eurasian.id")"
        # MAF re-applied: variants common across the whole panel can be
        # monomorphic once the outgroup is gone, and ADMIXTURE assumes unlinked
        # polymorphic markers. Re-pruned for the same reason.
        ${PLINK2} --bfile "${G}/merged_pca" --remove "${G}/non_eurasian.id" \
            --maf "${MAF}" \
            --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_R2}" \
            --threads "${THREADS}" --out "${G}/admix_prune"
        ${PLINK2} --bfile "${G}/merged_pca" --remove "${G}/non_eurasian.id" \
            --maf "${MAF}" --extract "${G}/admix_prune.prune.in" \
            --make-bed --threads "${THREADS}" --out "${G}/merged_admix"
        local NAV; NAV=$(nlines "${G}/merged_admix.bim")
        log "  merged_admix: $(nlines "${G}/merged_admix.fam") samples, ${NAV} SNPs"
        [[ "${NAV}" -ge "${ADMIX_MIN_SNPS}" ]] \
            || log "  WARNING: ${NAV} pruned SNPs is below ADMIX_MIN_SNPS=${ADMIX_MIN_SNPS}; \
ADMIXTURE resolution inside Eurasia will be limited. Raise PRUNE_R2."
        [[ "${NAV}" -le "${ADMIX_MAX_SNPS}" ]] \
            || log "  NOTE: ${NAV} SNPs exceeds ADMIX_MAX_SNPS=${ADMIX_MAX_SNPS}; runtime grows \
linearly with no gain in resolution. Consider tightening PRUNE_R2."

        # Training set for the projection design: same variants, cohort capped.
        # Same .bim is mandatory — ADMIXTURE's -P reads a .P.in whose rows are
        # the variants of the target .bed, in order.
        ${PLINK2} --bfile "${G}/merged_admix" --keep "${R}/pca_capped.keep" \
            --make-bed --threads "${THREADS}" --out "${G}/merged_admix_train"
        log "  merged_admix_train: $(nlines "${G}/merged_admix_train.fam") samples, \
$(nlines "${G}/merged_admix_train.bim") SNPs"
        cmp -s "${G}/merged_admix.bim" "${G}/merged_admix_train.bim" \
            || die "merged_admix_train.bim differs from merged_admix.bim; -P projection would be invalid"
    else
        log "[6/8] ADMIXTURE set: not built for this group"
    fi

    # --- 7. Fst / ROH set ---------------------------------------------------
    if [[ "${WANT_FSTROH}" -eq 1 ]]; then
        log "[7/8] Fst / ROH set (no MAF, no pruning, no cap, diploid only)"
        # MAF filtering before Fst is an ascertainment step: it preferentially
        # removes sites rare in the larger panels, which shifts every pairwise
        # value. Bhatia et al. 2013.
        ${PLINK2} --bfile "${G}/merged_raw" --keep "${R}/analysis.keep" \
            --geno "${GENO}" \
            --make-bed --threads "${THREADS}" --out "${G}/merged_fst"
        ln -sf "${G}/merged_fst.bed" "${G}/merged_roh.bed"
        ln -sf "${G}/merged_fst.bim" "${G}/merged_roh.bim"
        ln -sf "${G}/merged_fst.fam" "${G}/merged_roh.fam"
    else
        log "[7/8] Fst / ROH set: not built for this group"
    fi

    # --- 8. sample flow -----------------------------------------------------
    log "[8/8] sample flow"
    { printf 'stage\tsamples\tvariants\n'
      printf 'merged_raw\t%s\t%s\n'        "$(nlines "${G}/merged_raw.fam")" "$(nlines "${G}/merged_raw.bim")"
      printf 'pseudohaploid\t%s\t-\n'      "${NPH}"
      printf 'king_markers\t-\t%s\n'       "${NKM}"
      printf 'unrelated_diploid\t%s\t-\n'  "$(nlines "${R}/analysis.keep")"
      printf 'kazakh_unrelated\t%s\t-\n'   "${N_KAZ}"
      printf 'pca_capped_ref\t%s\t-\n'     "$(nlines "${R}/pca.keep")"
      printf 'pca_capped_both\t%s\t-\n'    "$(nlines "${R}/pca_capped.keep")"
      printf 'merged_pca\t%s\t%s\n'        "$(nlines "${G}/merged_pca.fam")" "$(nlines "${G}/merged_pca.bim")"
      [[ -s "${G}/merged_admix.bim" ]] && \
        printf 'merged_admix\t%s\t%s\n'    "$(nlines "${G}/merged_admix.fam")" "$(nlines "${G}/merged_admix.bim")"
      [[ -s "${G}/merged_admix_train.bim" ]] && \
        printf 'merged_admix_train\t%s\t%s\n' "$(nlines "${G}/merged_admix_train.fam")" "$(nlines "${G}/merged_admix_train.bim")"
      [[ -s "${G}/merged_fst.bim" ]] && \
        printf 'merged_fst\t%s\t%s\n'      "$(nlines "${G}/merged_fst.fam")" "$(nlines "${G}/merged_fst.bim")"
      true
    } > "${R}/sample_flow.tsv"
    column -t -s $'\t' "${R}/sample_flow.tsv"
}

for GROUP in "${ANALYSIS_GROUPS[@]}"; do run_group "${GROUP}"; done

log "=== 03 done -> 04_pca.sh ==="
