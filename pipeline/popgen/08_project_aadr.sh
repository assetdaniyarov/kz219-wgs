#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 08_project_aadr.sh — OPTIONAL, supplementary figure only.
#
# The 21 AADR groups in aadr_keep_groups.txt are the only publicly available
# genome-wide data for Kazakh, Kyrgyz, Uzbek, Turkmen, Karakalpak, Bashkir,
# Tatar, Kalmyk and Teleut. They are also 46 pseudo-haploid individuals at n=1-4
# with 1240K coverage between 228k and 1.23M sites, which is why they are not in
# any analysis panel: a pseudo-haploid genotype is one randomly drawn read
# written as a homozygote, and it breaks heterozygosity, ROH, ADMIXTURE and the
# (n-1) term in Hudson's Fst without raising a warning anywhere.
#
# Projection is the one operation they survive. Principal components are
# computed from diploid individuals only; the pseudo-haploid samples are then
# scored onto those fixed axes and never influence them. This is the standard
# treatment of pseudo-haploid data in a diploid PC space.
#
# HOW MUCH SHRINKAGE, MEASURED TWO WAYS.
#
# The naive comparison - all 46 AADR individuals against all 915 diploids -
# gives a mean radius of 0.073 against 0.117, which looks like heavy shrinkage.
# It is confounded: every AADR group is Central Asian or Siberian and sits near
# the centre of a Eurasian PC space anyway, while the diploid set includes
# Papuans and Africans on the rim. Part of that 62% is geography.
#
# The like-for-like control is the matched pairs in AADR_MATCHED_PAIRS: three
# AADR groups have a diploid counterpart in the same projection, and neither
# side defines the axes, so both carry the same reference-error shrinkage and
# the ratio isolates what pseudo-haploidy and missingness add. Measured:
# centroid norm ratio 0.94 (0.82-1.02) at an angle of 4.5 degrees. Direction
# preserved, magnitude barely touched.
#
# So the caveat has to be scoped, not blanket. The three matched groups are the
# near-complete ones (F_MISS 0.0005-0.25), and for individuals like those the
# projection is close to quantitative. The samples that stay strictly
# qualitative are the ones at F_MISS 0.7-0.8 - Kalmyk, Teleut, Tomsk Tatar,
# Uzbek - which have no diploid counterpart and therefore no direct measurement.
# aadr_groups.tsv gives the reader the missingness per group so the two classes
# are distinguishable on the figure.
#
# A better external Kazakh comparison exists: Serikzhan et al. 2025, n=224
# diploid array genotypes. If that dataset can be merged, this script becomes
# unnecessary and the comparison stops being three pseudo-haploid points.
#
#   bash 08_project_aadr.sh
# ============================================================================

need "${PLINK1}"; need "${PLINK2}"
TABLES="${SCRIPT_DIR}/popgen_tables.py"
require_file "${TABLES}"

GROUP="groupB_kz_hgdp"
G="${WORK_DIR}/${GROUP}"
R="${RES_DIR}/${GROUP}/pca_proj_aadr"; mkdir -p "${R}"
P="${WORK_DIR}/proj_aadr"; mkdir -p "${P}"

require_file "${G}/merged_pca.bed"
require_file "$(panel_path aadr).bed"
require_file "${RES_DIR}/${GROUP}/pca.keep"

log "=== projecting AADR pseudo-haploid samples onto ${GROUP} axes ==="

# --- 1. positions shared by the pruned diploid set and the AADR panel -------
cut -f2 "${G}/merged_pca.bim" | sort -u > "${P}/pca_ids.txt"
cut -d: -f1,2 "$(panel_path aadr).mergeable.snplist" | sort -u > "${P}/aadr_pos.txt"
awk 'NR==FNR{k[$1]=1; next} ($1 in k)' "${P}/aadr_pos.txt" "${P}/pca_ids.txt" \
    > "${P}/common.pos"
NC=$(nlines "${P}/common.pos")
log "  pruned diploid markers: $(nlines "${P}/pca_ids.txt")"
log "  shared with AADR:       ${NC}"
[[ "${NC}" -gt 5000 ]] || die "only ${NC} markers shared with AADR; projection would be meaningless"

# --- 2. AADR subset with CHR:POS IDs, same convention as 02 ----------------
awk -F: 'NR==FNR{k[$1":"$2]=1; next} ($1":"$2) in k' \
    "${P}/common.pos" "$(panel_path aadr).mergeable.snplist" > "${P}/aadr_ids.snplist"
# The panel on disk may predate a change to the inclusion list, so the
# currently selected individuals are re-resolved here. Without this, a group
# dropped in 01 would still be projected and would appear on the figure with
# no label.
awk -F'\t' 'NR>1 && $6=="keep" {print $1}' "${META_DIR}/aadr_selection.tsv" \
    > "${P}/aadr_keep.ids"
mk_id_list "$(panel_path aadr).fam" "${P}/aadr_keep.id" "${P}/aadr_keep.ids"
rm -f "${P}/aadr_keep.ids"
log "  AADR individuals currently selected: $(nlines "${P}/aadr_keep.id") \
of $(nlines "$(panel_path aadr).fam") in the panel"
${PLINK2} --bfile "$(panel_path aadr)" --extract "${P}/aadr_ids.snplist" \
    --keep "${P}/aadr_keep.id" \
    --make-bed --threads "${THREADS}" --out "${P}/aadr_sub"
awk 'BEGIN{OFS="\t"}{$2=$1":"$4; print}' "${P}/aadr_sub.bim" > "${P}/aadr_subn.bim"
ln -sf "${P}/aadr_sub.bed" "${P}/aadr_subn.bed"
cp     "${P}/aadr_sub.fam" "${P}/aadr_subn.fam"

${PLINK2} --bfile "${G}/merged_pca" --extract "${P}/common.pos" \
    --make-bed --threads "${THREADS}" --out "${P}/dip_sub"

# --- 3. merge, resolving strand the same way 02 does -----------------------
# Conflicts are classified rather than all called strand conflicts, using the
# same helper as 02. On a 1240K capture panel against a 650k array most of them
# genuinely are strand, but that has to be shown, not assumed.
rm -f "${P}/both-merge.missnp"
${PLINK1} --bfile "${P}/dip_sub" --bmerge "${P}/aadr_subn" --allow-no-sex \
    --make-bed --threads "${THREADS}" --out "${P}/both" || true
if [[ -f "${P}/both-merge.missnp" ]]; then
    NM=$(nlines "${P}/both-merge.missnp")
    cp "${P}/both-merge.missnp" "${P}/conflicts_initial.missnp"
    log "  merge conflicts: ${NM}"
    classify_missnp "${P}/conflicts_initial.missnp" "${P}/dip_sub.bim" \
                    "${P}/aadr_subn.bim" "${P}/conflicts_initial.tsv" \
      | sort | sed 's/^/    /'
    ${PLINK1} --bfile "${P}/aadr_subn" --flip "${P}/both-merge.missnp" --allow-no-sex \
        --make-bed --threads "${THREADS}" --out "${P}/aadr_flip"
    rm -f "${P}/both-merge.missnp"
    ${PLINK1} --bfile "${P}/dip_sub" --bmerge "${P}/aadr_flip" --allow-no-sex \
        --make-bed --threads "${THREADS}" --out "${P}/both" || true
    if [[ -f "${P}/both-merge.missnp" ]]; then
        NU=$(nlines "${P}/both-merge.missnp")
        cp "${P}/both-merge.missnp" "${P}/conflicts_unresolved.missnp"
        log "  resolved by --flip: $(( NM - NU )); still unresolved: ${NU} (dropped)"
        classify_missnp "${P}/conflicts_unresolved.missnp" "${P}/dip_sub.bim" \
                        "${P}/aadr_subn.bim" "${P}/conflicts_unresolved.tsv" \
          | sort | sed 's/^/    /'
        ${PLINK1} --bfile "${P}/dip_sub"   --exclude "${P}/both-merge.missnp" --allow-no-sex \
            --make-bed --threads "${THREADS}" --out "${P}/dip_cl"
        ${PLINK1} --bfile "${P}/aadr_flip" --exclude "${P}/both-merge.missnp" --allow-no-sex \
            --make-bed --threads "${THREADS}" --out "${P}/aadr_cl"
        ${PLINK1} --bfile "${P}/dip_cl" --bmerge "${P}/aadr_cl" --allow-no-sex \
            --make-bed --threads "${THREADS}" --out "${P}/both"
    else
        log "  all ${NM} resolved by --flip"
    fi
fi
log "  merged: $(nlines "${P}/both.fam") samples, $(nlines "${P}/both.bim") markers"

# --- 4. who defines the axes, who is projected -----------------------------
# Axes come from the diploid reference individuals only: not the AADR samples,
# and not the Kazakh cohort either, so the cohort and the AADR Kazakh samples
# are scored onto the same independent basis and can be compared to each other.
awk -F'\t' 'NR>1 && $3!="KAZ" && $4!="AADR" {print $1}' "${META_DIR}/pop_master.tsv" \
  > "${P}/axis.ids"
mk_id_list "${P}/both.fam" "${P}/axis.id" "${P}/axis.ids"
rm -f "${P}/axis.ids"
log "  axis-defining diploid individuals: $(nlines "${P}/axis.id")"

${PLINK2} --bfile "${P}/both" --keep "${P}/axis.id" \
    --freq counts --pca "${N_PC}" allele-wts \
    --threads "${THREADS}" --out "${R}/ref"

# .eigenvec.allele carries one row per allele and its column layout differs
# between plink2 builds, so the columns are resolved from the header. Passing
# "2 5" blind pointed --score at the wrong field in an earlier revision and
# every entry was rejected as a mismatching allele code.
read -r IDC A1C PC1C <<< "$(awk 'NR==1 {
    for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }
    if (!("ID" in h) || !("A1" in h) || !("PC1" in h)) {
        print "eigenvec.allele header lacks ID/A1/PC1" > "/dev/stderr"; exit 3 }
    printf "%d %d %d\n", h["ID"], h["A1"], h["PC1"]; exit }' "${R}/ref.eigenvec.allele")"
[[ -n "${PC1C:-}" ]] || die "could not resolve columns in ref.eigenvec.allele"
log "  eigenvec.allele columns: ID=${IDC} A1=${A1C} PC1=${PC1C}"

${PLINK2} --bfile "${P}/both" \
    --read-freq "${R}/ref.acount" \
    --score "${R}/ref.eigenvec.allele" "${IDC}" "${A1C}" header-read \
            no-mean-imputation variance-standardize \
    --score-col-nums "${PC1C}-$(( PC1C + N_PC - 1 ))" \
    --threads "${THREADS}" --out "${R}/scores"
require_file "${R}/scores.sscore"

${PYTHON} "${TABLES}" score-table --sscore "${R}/scores.sscore" \
    --meta "${META_DIR}/pop_master.tsv" --npc "${N_PC}" --out "${R}/scores.tsv"

# --- 5. per-sample missingness, because the panel must show it --------------
${PLINK2} --bfile "${P}/both" --missing sample-only \
    --threads "${THREADS}" --out "${R}/scores"

{ printf 'IID\tNation\tPopulation\tSource\tploidy\tF_MISS\n'
  awk -F'\t' -v OFS='\t' '
      FNR==NR { if (FNR>1) { nat[$1]=$2; pop[$1]=$3; src[$1]=$4 }; next }
      FNR==1  { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; next }
      { id=$(h["IID"])
        if (!(id in nat)) next
        print id, nat[id], pop[id], src[id],
              (src[id]=="AADR" ? "pseudo-haploid" : "diploid"), $(h["F_MISS"]) }' \
      "${META_DIR}/pop_master.tsv" "${R}/scores.smiss"
} > "${R}/projected_samples.tsv"

# per-group n and mean missingness: both belong on the figure legend
{ printf 'Nation\tn\tmean_F_MISS\tmax_F_MISS\n'
  awk -F'\t' 'NR>1 && $4=="AADR" { n[$2]++; s[$2]+=$6; if ($6>mx[$2]) mx[$2]=$6 }
      END { for (g in n) printf "%s\t%d\t%.4f\t%.4f\n", g, n[g], s[g]/n[g], mx[g] }' \
      "${R}/projected_samples.tsv" | sort -t$'\t' -k2,2nr
} > "${R}/aadr_groups.tsv"
log "  AADR groups on the panel:"
column -t -s $'\t' "${R}/aadr_groups.tsv" | sed 's/^/    /'

# --- 6. shrinkage, measured ------------------------------------------------
# Projected radius in the PC1-PC2 plane against per-sample missingness. Computed
# separately for the pseudo-haploid and the diploid individuals, so the reader
# can see both that the effect exists and how much of the AADR displacement it
# accounts for.
{ printf 'set\tn\tmean_radius\tpearson_r_radius_vs_FMISS\n'
  for SET in AADR diploid; do
      awk -F'\t' -v set="${SET}" '
          FNR==NR { if (FNR>1) { fm[$1]=$6; sr[$1]=$4 }; next }
          FNR==1  { for(i=1;i<=NF;i++) h[$i]=i; next }
          { id=$(h["IID"])
            if (!(id in fm)) next
            is = (sr[id]=="AADR") ? "AADR" : "diploid"
            if (is != set) next
            p1=$(h["PC1_AVG"]); p2=$(h["PC2_AVG"])
            r=sqrt(p1*p1+p2*p2); m=fm[id]
            n++; sx+=m; sy+=r; sxx+=m*m; syy+=r*r; sxy+=m*r }
          END { if (n<3) { printf "%s\t%d\tNA\tNA\n", set, n+0; exit }
                num=sxy-sx*sy/n; den=sqrt((sxx-sx*sx/n)*(syy-sy*sy/n))
                printf "%s\t%d\t%.6f\t%s\n", set, n, sy/n,
                       (den>0 ? sprintf("%.4f", num/den) : "NA") }' \
          "${R}/projected_samples.tsv" "${R}/scores.sscore"
  done
} > "${R}/shrinkage_stats.tsv"
column -t -s $'\t' "${R}/shrinkage_stats.tsv" | sed 's/^/    /'

# --- 6b. the like-for-like control -----------------------------------------
# Comparing all 46 AADR individuals against all 915 diploids confounds
# shrinkage with geography - every AADR group is Central Asian or Siberian and
# sits near the centre of a Eurasian PC space anyway. These three groups have a
# diploid counterpart in the same projection, which isolates the effect.
${PYTHON} "${TABLES}" proj-shrinkage --scores "${R}/scores.tsv" \
    --pairs "${AADR_MATCHED_PAIRS}" --npc 4 \
    --missingness "${R}/projected_samples.tsv" \
    --out "${R}/matched_shrinkage.tsv" \
  || log "  matched-pair shrinkage skipped (no counterpart resolved)"
[[ -s "${R}/matched_shrinkage.tsv" ]] && \
    column -t -s $'\t' "${R}/matched_shrinkage.tsv" | sed 's/^/    /'

log "  AADR individuals, worst missingness first:"
awk -F'\t' 'NR>1 && $4=="AADR"{print}' "${R}/projected_samples.tsv" \
  | sort -t$'\t' -k6,6gr | head -10 \
  | awk -F'\t' '{printf "    %-26s %-24s F_MISS=%s\n", $1, $2, $6}'

log "=== 08 done -> 09_kz_tribes.sh ==="
log "Figure: ${R}/scores.tsv, metadata ${R}/projected_samples.tsv"
log "Plot AADR as open symbols, label n per group from aadr_groups.tsv, and state"
log "in the legend that these are pseudo-haploid projections that did not"
log "contribute to the axes. Cite matched_shrinkage.tsv, not shrinkage_stats.tsv:"
log "the matched pairs are the like-for-like measurement, the naive AADR-vs-all"
log "comparison confounds shrinkage with the fact that every AADR group is"
log "Central Asian and sits near the centre of the plot regardless."
