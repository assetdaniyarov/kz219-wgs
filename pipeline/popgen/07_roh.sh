#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 07_roh.sh — runs of homozygosity, FROH, and its length distribution.
#
# Six constraints decide the design.
#
#  1. --homozyg exists only in PLINK 1.9. plink2 has no ROH implementation.
#
#  2. ROH must never be called on LD-pruned data. Pruning removes exactly the
#     correlated markers that constitute a run, so FROH collapses toward zero.
#     Input is merged_roh / the Kazakh panel, both unpruned.
#
#  3. Cross-population ROH runs on groupB, not groupA. 1000G Phase 3 is
#     low-coverage plus imputation; imputed haplotypes are smoother than
#     directly called ones and inflate long runs, so a KZ-vs-1000G FROH
#     comparison measures the calling pipeline as much as the demography. SGDP
#     is honest high-coverage but contributes n=1-3 per population, which is not
#     a population estimate. groupB is the Kazakh WGS panel against HGDP 650k,
#     one uniform marker set after the intersection, diploid throughout, and
#     n=9-25 per population. The residual caveat, which belongs in the text: the
#     Kazakh genotypes at those same sites come from 30x WGS, so their error
#     rate is lower and their runs are if anything slightly longer than the
#     reference panels'.
#
#  4. Parameters are density-dependent. The merged groupB set is an order of
#     magnitude sparser than the Kazakh WGS panel, so ROH_ARRAY and ROH_WGS are
#     separate parameter sets and the two FROH values are on different scales
#     and never share an axis.
#
#  5. ROH_MAF restricts to common variants. Rare-variant genotyping error
#     produces isolated heterozygous calls that truncate otherwise real runs.
#     This is a MAF filter, not LD pruning; the two are unrelated here. It is
#     also the parameter most open to challenge, so the sweep includes an arm
#     with no MAF threshold at all.
#
#  6. LENGTH CLASSES, NOT JUST A TOTAL. Total FROH cannot distinguish a genome
#     with many short runs from one with a few very long ones, and those mean
#     different things: short runs are background relatedness in a small
#     historical population, long runs are recent consanguinity. Reporting only
#     a total invites exactly that question in review. FROH is broken into
#     ROH_CLASSES by segment length, per individual and per population.
#
# The sensitivity sweep answers "would this hold under other ROH settings" in
# the supplement instead of in review.
# ============================================================================

need "${PLINK1}"; need "${PLINK2}"
TABLES="${SCRIPT_DIR}/popgen_tables.py"
require_file "${TABLES}"

# homozyg <bfile> <out> <param spec> [keep list] [maf]
# param spec: name:window-snp:window-het:homozyg-snp:homozyg-kb:density
homozyg() {
    # ${5-...} not ${5:-...}: an explicitly empty 5th argument must stay
    # empty, so the nomaf sweep arm can turn the filter off entirely.
    local BFILE=$1 OUT=$2 SPEC=$3 KEEP=${4:-} MAFV=${5-${ROH_MAF}}
    local NAME WS WH HS HK DE
    IFS=: read -r NAME WS WH HS HK DE <<< "${SPEC}"
    log "  ROH '${NAME}': window-snp=${WS} window-het=${WH} snp=${HS} kb=${HK} density=${DE} maf=${MAFV}"
    # shellcheck disable=SC2086
    ${PLINK1} --bfile "${BFILE}" \
        ${KEEP:+--keep "${KEEP}"} \
        ${MAFV:+--maf "${MAFV}"} \
        --homozyg \
        --homozyg-window-snp       "${WS}" \
        --homozyg-window-het       "${WH}" \
        --homozyg-window-missing   "${ROH_WINDOW_MISSING}" \
        --homozyg-window-threshold "${ROH_WINDOW_THRESHOLD}" \
        --homozyg-snp              "${HS}" \
        --homozyg-kb               "${HK}" \
        --homozyg-density          "${DE}" \
        --homozyg-gap              "${ROH_GAP}" \
        --allow-no-sex --threads "${THREADS}" --out "${OUT}"
}

# ---------------------------------------------------------------------------
# 1. cross-population ROH, shared marker set
# ---------------------------------------------------------------------------
GROUP="${ROH_GROUP}"
G="${WORK_DIR}/${GROUP}"
R="${RES_DIR}/${GROUP}/roh"; mkdir -p "${R}"
log "=== cross-population ROH (${GROUP}) ==="
require_file "${G}/merged_roh.bed"
NMARK=$(nlines "${G}/merged_roh.bim")
log "  markers: ${NMARK} (unpruned)"
homozyg "${G}/merged_roh" "${R}/roh_shared" "${ROH_ARRAY}"
log "  segments: $(nbody "${R}/roh_shared.hom")"

${PYTHON} "${TABLES}" roh-summary \
    --hom "${R}/roh_shared.hom" --indiv "${R}/roh_shared.hom.indiv" \
    --meta "${META_DIR}/pop_master.tsv" --autosome-bp "${AUTOSOME_BP}" \
    --classes "${ROH_CLASSES}" --by Population \
    --out "${R}/froh_individual.tsv" --out-by-group "${R}/froh_by_population.tsv"
${PYTHON} "${TABLES}" roh-summary \
    --hom "${R}/roh_shared.hom" --indiv "${R}/roh_shared.hom.indiv" \
    --meta "${META_DIR}/pop_master.tsv" --autosome-bp "${AUTOSOME_BP}" \
    --classes "${ROH_CLASSES}" --by Nation \
    --out "${R}/froh_individual.tsv" --out-by-group "${R}/froh_by_nation.tsv"
log "  FROH by population:"
column -t -s $'\t' "${R}/froh_by_population.tsv" | sed 's/^/    /'
log "  Nation-level FROH (the comparison with the cohort's actual neighbours):"
head -12 "${R}/froh_by_nation.tsv" | column -t -s $'\t' | sed 's/^/    /'

# platform annotation travels with the results so the caveat is enforceable
# rather than remembered
{ echo -e "IID\tPlatform\tSource"
  awk -F'\t' 'NR>1{print $1"\t"$5"\t"$4}' "${META_DIR}/pop_master.tsv"; } \
  > "${R}/platform.tsv"

{ printf 'key\tvalue\n'
  printf 'group\t%s\n'      "${GROUP}"
  printf 'markers\t%s\n'    "${NMARK}"
  printf 'param_set\t%s\n'  "${ROH_ARRAY}"
  printf 'maf\t%s\n'        "${ROH_MAF}"
  printf 'classes\t%s\n'    "${ROH_CLASSES}"
  printf 'froh_denominator_bp\t%s\n' "${AUTOSOME_BP}"
} > "${R}/roh_shared.params.tsv"

# ---------------------------------------------------------------------------
# 2. Kazakh-only ROH on the full WGS marker set
# ---------------------------------------------------------------------------
RK="${RES_DIR}/kz_only/roh"; mkdir -p "${RK}"
KZP="$(panel_path kz)"
log "=== Kazakh-only ROH, full WGS marker set ==="
require_file "${KZP}.bed"
log "  markers: $(nlines "${KZP}.bim")"

# Restrict to the individuals that survived relatedness filtering in groupA,
# which has the most markers and therefore the best-resolved kinship. The FIDs
# in analysis.keep come from the merged fileset, so the IIDs are re-resolved
# against the Kazakh panel's own .fam.
KEEPSRC="${RES_DIR}/groupA_kz_1000g_sgdp/analysis.keep"
require_file "${KEEPSRC}"
awk '{print $2}' "${KEEPSRC}" > "${RK}/kz_unrelated.iids"
mk_id_list "${KZP}.fam" "${RK}/kz_unrelated.id" "${RK}/kz_unrelated.iids"
rm -f "${RK}/kz_unrelated.iids"
log "  unrelated Kazakh individuals: $(nlines "${RK}/kz_unrelated.id")"

homozyg "${KZP}" "${RK}/roh_kz_wgs" "${ROH_WGS}" "${RK}/kz_unrelated.id"
log "  segments: $(nbody "${RK}/roh_kz_wgs.hom")"

${PYTHON} "${TABLES}" roh-summary \
    --hom "${RK}/roh_kz_wgs.hom" --indiv "${RK}/roh_kz_wgs.hom.indiv" \
    --meta "${META_DIR}/pop_master.tsv" --autosome-bp "${AUTOSOME_BP}" \
    --classes "${ROH_CLASSES}" --by Population \
    --out "${RK}/froh_individual.tsv" --out-by-group "${RK}/froh_summary.tsv"
log "  Kazakh WGS FROH:"
column -t -s $'\t' "${RK}/froh_summary.tsv" | sed 's/^/    /'

# ---------------------------------------------------------------------------
# 3. parameter sensitivity sweep, Kazakh WGS set
# ---------------------------------------------------------------------------
log "=== ROH sensitivity sweep ==="
mkdir -p "${RK}/sensitivity"
SWEEP=(
  "lenient:30:2:30:300:50"
  "${ROH_WGS}"
  "strict:100:0:100:1000:50"
  "long_only:50:1:50:1500:50"
)
for S in "${SWEEP[@]}"; do
    NAME="${S%%:*}"
    homozyg "${KZP}" "${RK}/sensitivity/${NAME}" "${S}" "${RK}/kz_unrelated.id" \
      > "${RK}/sensitivity/${NAME}.stdout" 2>&1 \
      || log "    ${NAME} failed, see ${RK}/sensitivity/${NAME}.stdout"
done
# The MAF threshold is the parameter most open to challenge, because filtering
# on frequency before calling ROH is not universal practice. This arm applies
# none at all, on the same default window settings.
homozyg "${KZP}" "${RK}/sensitivity/nomaf" "${ROH_WGS/wgs:/nomaf:}" \
        "${RK}/kz_unrelated.id" "" \
  > "${RK}/sensitivity/nomaf.stdout" 2>&1 \
  || log "    nomaf failed, see ${RK}/sensitivity/nomaf.stdout"

SENS=()
for F in "${RK}"/sensitivity/*.hom.indiv; do [[ -s "${F}" ]] && SENS+=("${F}"); done
if [[ "${#SENS[@]}" -ge 2 ]]; then
    ${PYTHON} "${TABLES}" roh-sensitivity --indiv "${SENS[@]}" \
        --reference "${ROH_WGS%%:*}" --autosome-bp "${AUTOSOME_BP}" \
        --out "${RK}/sensitivity/froh_by_paramset.tsv" \
        --out-stats "${RK}/sensitivity/froh_paramset_stats.tsv"
    log "  sensitivity (correlation with the default parameter set):"
    column -t -s $'\t' "${RK}/sensitivity/froh_paramset_stats.tsv" | sed 's/^/    /'
    log "  A Spearman correlation near 1 across parameter sets means the RANKING"
    log "  of individuals by FROH is parameter-independent even where the absolute"
    log "  values are not. That is the claim the manuscript should make."
else
    log "  fewer than two parameter sets produced output; sweep not summarised"
fi

# ---------------------------------------------------------------------------
# 4. inbreeding coefficients as an independent check on FROH
# ---------------------------------------------------------------------------
log "=== F estimates (method-of-moments and --ibc) ==="
# Same MAF threshold as the ROH run, so FROH and Fhat3 are computed on the same
# marker set and the correlation between them is interpretable. Individuals
# flagged in 01 as heterozygosity-deficit outliers and retained should sit at
# the top of both; that agreement is the evidence they are consanguineous rather
# than technical failures.
${PLINK2} --bfile "${KZP}" --keep "${RK}/kz_unrelated.id" \
    --maf "${ROH_MAF}" --het --threads "${THREADS}" --out "${RK}/kz_het"
# --ibc exists only in PLINK 1.9; plink2 dropped it and rejects the flag.
${PLINK1} --bfile "${KZP}" --keep "${RK}/kz_unrelated.id" \
    --maf "${ROH_MAF}" --ibc --allow-no-sex \
    --threads "${THREADS}" --out "${RK}/kz_ibc"
require_file "${RK}/kz_ibc.ibc"
log "  Fhat estimates: $(nbody "${RK}/kz_ibc.ibc") individuals"

${PYTHON} "${TABLES}" froh-vs-f --froh "${RK}/froh_individual.tsv" \
    --het "${RK}/kz_het.het" --ibc "${RK}/kz_ibc.ibc" \
    --out-join "${RK}/froh_vs_f.tsv" --out-stats "${RK}/froh_vs_f_stats.tsv"
column -t -s $'\t' "${RK}/froh_vs_f_stats.tsv" | sed 's/^/    /'

if [[ -s "${KZP}.het_flag.tsv" ]]; then
    cp "${KZP}.het_flag.tsv" "${RK}/het_deficit_flagged.tsv"
    log "  heterozygosity-deficit individuals flagged in 01 and retained, with their FROH:"
    awk -F'\t' 'NR==FNR{if(FNR>1) f[$1]=$8; next}
        FNR>1 { printf "    %-24s F=%-10s FROH=%s\n", $1, $2, ($1 in f ? f[$1] : "not in ROH set") }' \
        "${RK}/froh_individual.tsv" "${RK}/het_deficit_flagged.tsv"
fi

log "=== 07 done -> 08_project_aadr.sh ==="
log "FROH denominator: ${AUTOSOME_BP} bp (GRCh37 autosomes)"
log "Cross-population FROH (${ROH_ARRAY%%:*} params, ${NMARK} markers) and Kazakh"
log "WGS FROH (${ROH_WGS%%:*} params) are on different marker densities - never one axis."
