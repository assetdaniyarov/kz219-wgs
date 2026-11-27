#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 04_pca.sh — PCA on the diploid, unrelated, reference-capped set.
#
#   pca_main            all retained superpopulations               context
#   pca_eurasia         Eurasians only, recomputed                  main figure
#   pca_eurasia_capped  same markers, cohort capped at COHORT_CAP   sensitivity
#   pca_proj            reference-only PCs, everyone scored onto them
#
# Roles by group:
#   groupA  pca_main is the global WGS context panel (supplementary), and
#           pca_eurasia is the panel-independent replication of the main figure.
#   groupB  pca_eurasia is the main figure.
# 08_project_aadr.sh reuses the groupB axes to place the 46 pseudo-haploid AADR
# individuals; the allele weights written below are what it scores against.
#
# COHORT CAPPING, AND WHY THE SENSITIVITY PANEL EXISTS. Reference nations are
# capped at 25 and the 218-individual cohort is not, so the cohort is a quarter
# of the panel and roughly nine times the largest reference group. Sample size
# drives PC axes (McVean 2009), so the obvious objection is that the cohort's
# apparent position is partly an artefact of its own weight. pca_eurasia_capped
# recomputes the axes with the cohort capped at COHORT_CAP over the SAME marker
# set, the full Eurasian set is then scored onto those axes, and
# pca_eurasia_capped.concordance.tsv reports |r| per axis between the two. That
# is a number a reviewer can check, not a claim.
#
# PROJECTION AND SCALE. plink2's --score reports the mean over scored alleles,
# so .sscore coordinates are NOT on the same scale as the reference .eigenvec —
# they differ by a factor of roughly 1/(2M). Every individual in a given .sscore
# went through the same call, so the geometry inside one file is consistent.
# Never build one figure from .eigenvec for the reference and .sscore for the
# cohort; the clouds separate by a factor of thousands and it looks like signal.
# pca_proj.sscore therefore contains reference AND cohort, and the projected
# panels are drawn entirely from it.
#
# pca_proj is the check that the Kazakh position is not self-induced: a joint
# PCA lets the cohort contribute to the axes it is then plotted on. Here the PCs
# come from the reference panels only. Projected coordinates are shrunk toward
# the origin relative to a joint PCA, so the panel supports the qualitative
# placement, not a numerical comparison of coordinates between figures.
#
# VARIANCE EXPLAINED. plink2 writes only N_PC eigenvalues, so
# eigenvalue/sum(eigenvalues) gives the share of the retained PCs, not of total
# variance. For a standardised GRM the trace equals the sample size, so
#     %variance(PC_k) = eigenvalue_k / N_samples * 100
# That is what pca_*.varexp.tsv contains.
# ============================================================================

need "${PLINK2}"
TABLES="${SCRIPT_DIR}/popgen_tables.py"
require_file "${TABLES}"

# --------------------------------------------------------------------------
# score_cols <eigenvec.allele> -> "<ID col> <A1 col> <PC1 col>"
#
# plink2's .eigenvec.allele carries one row per ALLELE, so a biallelic panel
# gives two rows per variant, and the column layout has changed between builds:
# some write #CHROM ID REF ALT A1 PC1.., others insert POS. Hardcoding "2 5"
# silently pointed --score at the wrong field and every entry was rejected with
# "mismatching allele codes" - 161818 of 161818, which is 2 x 80909 variants.
# --------------------------------------------------------------------------
score_cols() {
    awk 'NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }
                 if (!("ID" in h) || !("A1" in h) || !("PC1" in h)) {
                     print "eigenvec.allele header lacks ID/A1/PC1" > "/dev/stderr"; exit 3 }
                 printf "%d %d %d\n", h["ID"], h["A1"], h["PC1"]; exit }' "$1"
}

# score_onto <bfile> <ref prefix> <out prefix>
score_onto() {
    local BF=$1 REF=$2 OUT=$3
    local IDC A1C PC1C
    read -r IDC A1C PC1C <<< "$(score_cols "${REF}.eigenvec.allele")"
    [[ -n "${PC1C:-}" ]] || die "could not resolve columns in ${REF}.eigenvec.allele"
    ${PLINK2} --bfile "${BF}" \
        --read-freq "${REF}.acount" \
        --score "${REF}.eigenvec.allele" "${IDC}" "${A1C}" header-read \
                no-mean-imputation variance-standardize \
        --score-col-nums "${PC1C}-$(( PC1C + N_PC - 1 ))" \
        --threads "${THREADS}" --out "${OUT}"
    require_file "${OUT}.sscore"
}

# ids_of <acount file> -> variant ID list on stdout
ids_of() {
    awk 'NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }
                 if (!("ID" in h)) { print "acount header lacks ID" > "/dev/stderr"; exit 3 }
                 c=h["ID"]; next }
         { print $c }' "$1"
}

run_group() {
    local GROUP=$1
    local G="${WORK_DIR}/${GROUP}"
    local R="${RES_DIR}/${GROUP}"; mkdir -p "${R}/pca"
    log "=== ${GROUP} ==="

    require_file "${G}/merged_pca.bed"
    require_file "${META_DIR}/pop_master.tsv"

    local NS; NS=$(nlines "${G}/merged_pca.fam")

    # AFR and AMR are already gone from the merge apart from the KEEP_NATIONS
    # outgroup. Both that outgroup and EURASIA_EXCLUDE leave the Eurasian PCA.
    awk -F'\t' -v drop="${DROP_SUPERPOPS} ${EURASIA_EXCLUDE}" '
        BEGIN { n=split(drop,a," "); for(i=1;i<=n;i++) D[a[i]]=1 }
        NR>1 && ($3 in D) { print $1 }' "${META_DIR}/pop_master.tsv" \
      > "${G}/non_eurasian.ids"
    mk_id_list "${G}/merged_pca.fam" "${G}/non_eurasian.id" "${G}/non_eurasian.ids"
    rm -f "${G}/non_eurasian.ids"
    log "  removed from the Eurasian PCA: $(nlines "${G}/non_eurasian.id")"

    # --- 1. main PCA --------------------------------------------------------
    log "[1/4] main PCA (${NS} individuals, $(nlines "${G}/merged_pca.bim") SNPs)"
    ${PLINK2} --bfile "${G}/merged_pca" \
        --freq counts --pca "${N_PC}" allele-wts \
        --threads "${THREADS}" --out "${R}/pca/pca_main"
    ${PYTHON} "${TABLES}" pca-table \
        --eigenvec "${R}/pca/pca_main.eigenvec" --meta "${META_DIR}/pop_master.tsv" \
        --npc "${N_PC}" --eigenval "${R}/pca/pca_main.eigenval" \
        --n-samples "${NS}" --varexp-out "${R}/pca/pca_main.varexp.tsv" \
        --out "${R}/pca/pca_main.tsv"

    # --- 2. Eurasian subset -------------------------------------------------
    # Recomputed, not a zoom. Residual non-Eurasians dominate PC1 without
    # separating anything inside Eurasia. MAF is re-applied because variants
    # common across the whole panel can be monomorphic once they are gone.
    log "[2/4] Eurasian PCA"
    ${PLINK2} --bfile "${G}/merged_pca" --remove "${G}/non_eurasian.id" \
        --maf "${MAF}" --freq counts --pca "${N_PC}" allele-wts \
        --threads "${THREADS}" --out "${R}/pca/pca_eurasia"
    local NSE; NSE=$(nbody "${R}/pca/pca_eurasia.eigenvec")
    ${PYTHON} "${TABLES}" pca-table \
        --eigenvec "${R}/pca/pca_eurasia.eigenvec" --meta "${META_DIR}/pop_master.tsv" \
        --npc "${N_PC}" --eigenval "${R}/pca/pca_eurasia.eigenval" \
        --n-samples "${NSE}" --varexp-out "${R}/pca/pca_eurasia.varexp.tsv" \
        --out "${R}/pca/pca_eurasia.tsv"
    log "  Eurasian PCA: ${NSE} individuals, $(nbody "${R}/pca/pca_eurasia.acount") SNPs"

    # --- 3. cohort-capped sensitivity --------------------------------------
    # Same markers as pca_eurasia, pinned with --extract, so the only thing that
    # differs between the two runs is how many Kazakhs are in the sample.
    log "[3/4] cohort-capped Eurasian PCA (KAZ <= ${COHORT_CAP})"
    if [[ -s "${G}/merged_pca_capped.bed" ]]; then
        ids_of "${R}/pca/pca_eurasia.acount" > "${G}/eurasia_vars.txt"
        ${PLINK2} --bfile "${G}/merged_pca_capped" --remove "${G}/non_eurasian.id" \
            --extract "${G}/eurasia_vars.txt" \
            --freq counts --pca "${N_PC}" allele-wts \
            --threads "${THREADS}" --out "${R}/pca/pca_eurasia_capped"
        ${PYTHON} "${TABLES}" pca-table \
            --eigenvec "${R}/pca/pca_eurasia_capped.eigenvec" \
            --meta "${META_DIR}/pop_master.tsv" --npc "${N_PC}" \
            --eigenval "${R}/pca/pca_eurasia_capped.eigenval" \
            --n-samples "$(nbody "${R}/pca/pca_eurasia_capped.eigenvec")" \
            --varexp-out "${R}/pca/pca_eurasia_capped.varexp.tsv" \
            --out "${R}/pca/pca_eurasia_capped.tsv"

        # |r| per axis over the individuals both runs contain. Eigenvector signs
        # are arbitrary, so the absolute value is the quantity with meaning.
        ${PYTHON} "${TABLES}" pca-compare \
            --a "${R}/pca/pca_eurasia.tsv" --b "${R}/pca/pca_eurasia_capped.tsv" \
            --npc 6 --out "${R}/pca/pca_eurasia_capped.concordance.tsv" \
            --out-subspace "${R}/pca/pca_eurasia_capped.subspace.tsv"
        sed 's/^/    /' "${R}/pca/pca_eurasia_capped.concordance.tsv"
        sed 's/^/    /' "${R}/pca/pca_eurasia_capped.subspace.tsv"

        # The full Eurasian set scored onto the capped axes: the supplementary
        # panel showing all 218 Kazakhs on axes they did not dominate.
        ${PLINK2} --bfile "${G}/merged_pca" --remove "${G}/non_eurasian.id" \
            --extract "${G}/eurasia_vars.txt" \
            --make-bed --threads "${THREADS}" --out "${G}/eurasia_all"
        score_onto "${G}/eurasia_all" "${R}/pca/pca_eurasia_capped" \
                   "${R}/pca/pca_eurasia_on_capped_axes"
        ${PYTHON} "${TABLES}" score-table \
            --sscore "${R}/pca/pca_eurasia_on_capped_axes.sscore" \
            --meta "${META_DIR}/pop_master.tsv" --npc "${N_PC}" \
            --out "${R}/pca/pca_eurasia_on_capped_axes.tsv"
    else
        log "  merged_pca_capped absent, skipping (rerun 03)"
    fi

    # --- 4. cohort projection ----------------------------------------------
    log "[4/4] reference PCA + cohort projection"
    awk -F'\t' 'NR>1 && $3!="KAZ"{print $1}' "${META_DIR}/pop_master.tsv" \
      > "${G}/reference.ids"
    awk -F'\t' 'NR>1 && $3=="KAZ"{print $1}' "${META_DIR}/pop_master.tsv" \
      > "${G}/cohort.ids"
    mk_id_list "${G}/merged_pca.fam" "${G}/reference.id" "${G}/reference.ids"
    mk_id_list "${G}/merged_pca.fam" "${G}/cohort.id"    "${G}/cohort.ids"
    rm -f "${G}/reference.ids" "${G}/cohort.ids"
    log "  reference $(nlines "${G}/reference.id") / cohort $(nlines "${G}/cohort.id")"

    ${PLINK2} --bfile "${G}/merged_pca" --keep "${G}/reference.id" \
        --freq counts --pca "${N_PC}" allele-wts \
        --threads "${THREADS}" --out "${R}/pca/ref"

    # Everyone is scored, reference included, so one file holds the whole figure
    # on one scale.
    score_onto "${G}/merged_pca" "${R}/pca/ref" "${R}/pca/pca_proj"
    ${PYTHON} "${TABLES}" score-table \
        --sscore "${R}/pca/pca_proj.sscore" --meta "${META_DIR}/pop_master.tsv" \
        --npc "${N_PC}" --out "${R}/pca/pca_proj.tsv"
    log "  projected: $(nbody "${R}/pca/pca_proj.sscore") individuals (reference + cohort, one scale)"

    if [[ "${PCA_ALL_SAMPLES}" -eq 1 && -s "${G}/merged_pca_all.bed" ]]; then
        log "[extra] all-sample PCA (supplementary)"
        ${PLINK2} --bfile "${G}/merged_pca_all" --pca "${N_PC}" \
            --threads "${THREADS}" --out "${R}/pca/pca_all"
        ${PYTHON} "${TABLES}" pca-table \
            --eigenvec "${R}/pca/pca_all.eigenvec" --meta "${META_DIR}/pop_master.tsv" \
            --npc "${N_PC}" --eigenval "${R}/pca/pca_all.eigenval" \
            --n-samples "$(nlines "${G}/merged_pca_all.fam")" \
            --varexp-out "${R}/pca/pca_all.varexp.tsv" --out "${R}/pca/pca_all.tsv"
    fi

    { printf 'key\tvalue\n'
      printf 'n_samples\t%s\n'          "${NS}"
      printf 'n_variants\t%s\n'         "$(nlines "${G}/merged_pca.bim")"
      printf 'n_eurasia\t%s\n'          "${NSE}"
      printf 'n_eurasia_variants\t%s\n' "$(nbody "${R}/pca/pca_eurasia.acount")"
      printf 'cap_per_nation\t%s\n'     "${PCA_MAX_PER_NATION}"
      printf 'cohort_cap_sensitivity\t%s\n' "${COHORT_CAP}"
      printf 'pseudohaploid_in_pcs\t0\n'
    } > "${R}/pca/pca_main.n"

    log "  eigenvec rows: $(nbody "${R}/pca/pca_main.eigenvec")"
    head -4 "${R}/pca/pca_eurasia.varexp.tsv" | sed 's/^/    /'
}

for GROUP in "${ANALYSIS_GROUPS[@]}"; do run_group "${GROUP}"; done

log "=== 04 done -> 05_admixture.sh ==="
log "No pseudo-haploid individual contributed to any eigenvector above."
log "Figures are drawn from pca_*.tsv; .eigenvec and .sscore coordinates are on"
log "different scales and must never appear in the same panel."
