#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 09_kz_tribes.sh — within-cohort structure of the Kazakh genomes by tribe.
#
# Runs on the Kazakh panel alone (work/panels/kz, 219 individuals, ~16.6M
# variants). No reference panel, so nothing here is limited by a marker
# intersection and the full WGS density is available.
#
# WHAT THIS CAN AND CANNOT ANSWER. Of 219 individuals, 186 carry a tribe label
# and only three tribes reach n = 10: Argyn 86, Kypshak 31, Kerei 8. Everything
# else is n <= 5, forty-one levels of it. So:
#
#   - Argyn vs Kypshak is the one comparison with real power. Kerei at n = 8 is
#     borderline and reported with its interval, not as a point estimate.
#   - Zhuz level is better balanced: Orta 146, Kishi 20, Uly 10. Tore and Qoja
#     are NOT zhuz - they are status groups outside the three-zhuz system, held
#     separate here rather than merged into one of them.
#   - Every other tribe appears in the PCA as a labelled point and in no test.
#
# EXPECT A NEGATIVE RESULT, AND MAKE IT A RESULT. Kazakh tribes are patrilineal
# and exogamous: marriage within seven generations of the same paternal line is
# prohibited, so wives come from outside the tribe every generation. That rule
# preserves tribal structure on the Y chromosome and actively erases it from the
# autosomes. A near-zero autosomal Fst between Argyn and Kypshak alongside a
# tribe-concordant Y-haplogroup profile is not two contradictory findings, it is
# the direct genetic signature of the marriage system.
#
# A null result is only publishable if it is bounded, so three things are
# reported rather than one:
#   - Fst with its jackknife 95% interval AND the smallest Fst the design could
#     have resolved (1.96 x SE). "Below 1.2e-4" is a measurement; "not
#     significant" is not.
#   - a permutation test of every PC against tribe and zhuz, so "no structure in
#     the PCA" stops being an impression of a scatter plot.
#   - mean pairwise kinship within tribes versus between tribes, permuted. This
#     is the most direct autosomal test of exogamy available and it uses all
#     ~23,000 pairs instead of a two-group summary.
#
# ---------------------------------------------------------------------------
# RELATEDNESS, AND THE ERROR THIS SCRIPT USED TO MAKE
#
# The 2026-08-30 revision ran --make-king-table here and applied KING's standard
# degree bands (>0.354 duplicate, >0.177 1st, >0.0884 2nd, >0.0442 3rd). The
# output was alarming: 23,449 of the 23,871 possible pairs - 98% of the cohort -
# came out "3rd degree or closer", 19,273 of them labelled "3rd degree". Taken
# at face value that is a claim that essentially every pair of Kazakhs in the
# study is a third cousin.
#
# It was an artefact, and NOT the one first suspected. The obvious explanation
# was the estimator: KING-robust builds its denominator from the pair's own
# heterozygote counts, so a population with reduced heterozygosity inflates
# every estimate together. That story is plausible and it was wrong.
#
# The 2026-09-02 run, with KZ_GQ_MIN=20 applied at import in 01 and nothing else
# about the estimator changed, gives:
#
#     kinship median 0.00306, IQR 0.00064-0.00540, max 0.12127
#     23,870 pairs below 2nd degree, exactly ONE pair at 2nd degree
#
# More markers (102,207 -> 175,344) cuts the noise by about a quarter and cannot
# move a median from above 0.044 to 0.003. That is a bias, and the only thing
# that changed the genotypes is genotype-level quality filtering. Low-GQ calls
# in this DRAGEN callset systematically under-call heterozygotes; KING reads
# heterozygote sharing, its denominator is heterozygote counts, and deflating
# those inflates kinship for every pair at once.
#
# So: the cohort has essentially no background relatedness, one 2nd-degree pair,
# and the earlier "98% are third cousins" was a genotype-quality artefact. That
# is worth a sentence in the manuscript - it is direct evidence that the GQ
# filter was necessary rather than decorative.
#
# The safeguards stay in place regardless, because they cost nothing and the
# next callset may not be this clean: nothing below 2nd degree gets a
# relationship name, the whole distribution is reported rather than a count
# above a threshold, PLINK 1.9 --genome PI_HAT cross-checks with a different
# bias structure, and KING is re-run on the groupA marker set restricted to the
# same 219 people so any marker-ascertainment effect is measured.
#
# ADMIXTURE from K = 1. Within a single population the honest question is
# whether K = 1 is beaten at all. Starting at K = 2 assumes the answer.
# ============================================================================

need "${PLINK1}"; need "${PLINK2}"; need "${ADMIXTURE}"
TABLES="${SCRIPT_DIR}/popgen_tables.py"; require_file "${TABLES}"
ALIGN="${SCRIPT_DIR}/admixture_align.py"; require_file "${ALIGN}"

KZP="$(panel_path kz)"
require_file "${KZP}.bed"
require_file "${MF_TRIBES}"

W="${WORK_DIR}/kz_tribes"; mkdir -p "${W}"
R="${RES_DIR}/kz_tribes";  mkdir -p "${R}"

# ---------------------------------------------------------------------------
log "=== [1/7] labels ==="
# ---------------------------------------------------------------------------
# Columns are resolved by NAME, with the historical position as a fallback. The
# previous revision hardcoded 2/4/6, so inserting a column in the hand-curated
# tribe sheet would have silently relabelled the entire analysis.
tribe_col() {   # tribe_col <header name> <fallback index>
    awk -F'\t' -v want="$1" -v fb="$2" 'NR==1{
        for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); gsub(/[ _.-]/,"",g); h[tolower(g)]=i }
        w=tolower(want); gsub(/[ _.-]/,"",w)
        print (w in h) ? h[w] : fb; exit }' "${MF_TRIBES}"
}

awk -F'\t' 'NR==FNR{g[$1]=1; next}
    FNR==1 { for(i=1;i<=NF;i++){ h[$i]=i }
             if (!("IID" in h)) { print "no IID column in the tribe file" > "/dev/stderr"; exit 2 }
             next }
    ($1 in g) { print }' <(fam_iids "${KZP}.fam") "${MF_TRIBES}" > "${W}/tribes_in_panel.tsv"
log "  labelled individuals present in the panel: $(nlines "${W}/tribes_in_panel.tsv") of $(nlines "${KZP}.fam")"

mkcat() {   # mkcat <column index> <name>
    local COL=$1 NAME=$2
    { printf '#IID\t%s\n' "${NAME}"
      awk -F'\t' -v c="${COL}" 'BEGIN{OFS="\t"} $c!="" && $c!="NA" {print $1, $c}' \
          "${W}/tribes_in_panel.tsv"; } > "${W}/${NAME}.cat"
    # same content, IID-keyed, for the Python helpers
    { printf 'IID\t%s\n' "${NAME}"; tail -n +2 "${W}/${NAME}.cat"; } > "${W}/${NAME}.labels.tsv"
    { printf '%s\tN\n' "${NAME}"
      awk -F'\t' 'NR>1{n[$2]++} END{for(k in n) printf "%s\t%d\n", k, n[k]}' \
          "${W}/${NAME}.cat" | sort -k2,2nr; } > "${R}/${NAME}_counts.tsv"
    log "  ${NAME}: $(nbody "${W}/${NAME}.cat") labelled, $(nbody "${R}/${NAME}_counts.tsv") levels (column $(printf %s "${COL}"))"
}
mkcat "$(tribe_col Tribe 2)"        Tribe
mkcat "$(tribe_col Zhuz 4)"         Zhuz
mkcat "$(tribe_col Mother_tribe 6)" MotherTribe

# ---------------------------------------------------------------------------
log "=== [2/7] relatedness within the cohort ==="
# ---------------------------------------------------------------------------
HIGHLD="${WORK_DIR}/highld_hg19.txt"
require_file "${HIGHLD}"

${PLINK2} --bfile "${KZP}" --maf 0.05 --geno 0.02 \
    --exclude range "${HIGHLD}" \
    --indep-pairwise "${KING_PRUNE_WINDOW}" "${KING_PRUNE_STEP}" "${KING_PRUNE_R2}" \
    --threads "${THREADS}" --out "${W}/king_prune"
${PLINK2} --bfile "${KZP}" --extract "${W}/king_prune.prune.in" \
    --make-bed --threads "${THREADS}" --out "${W}/king_in"
NKM=$(nlines "${W}/king_in.bim")
log "  KING markers: ${NKM}"
[[ "${NKM}" -ge "${KING_MIN_MARKERS}" ]] || \
    log "  WARNING: below KING_MIN_MARKERS=${KING_MIN_MARKERS}"

# NO --king-table-filter: every pair is written, because the distribution is the
# result here and a filtered table cannot show it.
${PLINK2} --bfile "${W}/king_in" --make-king-table \
    --threads "${THREADS}" --out "${R}/king_pairs"
NPAIR=$(nbody "${R}/king_pairs.kin0")
log "  pairs: ${NPAIR} (all of them; the tail is the point)"

{ printf 'bin_lo\tbin_hi\tn_pairs\n'
  awk 'NR==1{ for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; next }
       { b=int($(h["KINSHIP"])/0.005); c[b]++ }
       END { for (b in c) printf "%.3f\t%.3f\t%d\n", b*0.005, (b+1)*0.005, c[b] }' \
      "${R}/king_pairs.kin0" | sort -g
} > "${R}/kinship_hist.tsv"

# median/IQR without gawk's asort, which mawk does not have
awk 'NR==1{ for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; c=h["KINSHIP"]; next }
     { print $c }' "${R}/king_pairs.kin0" | sort -g > "${W}/kinship_sorted.txt"
awk 'BEGIN{n=0} {v[++n]=$1}
     END { if (n<4) exit
           printf "  kinship distribution: median %.5f  IQR %.5f-%.5f  max %.5f\n",
                  v[int(n*0.50)+1], v[int(n*0.25)+1], v[int(n*0.75)+1], v[n]
           if (v[int(n*0.50)+1] > 0.02)
               printf "  WARNING: a median this high is background, not genealogy - check that\n  genotype-level QC ran (KZ_GQ_MIN) before reading any degree label\n"
           else
               printf "  median well below the 4th-degree band: the tail is genealogy, not baseline\n" }' \
    "${W}/kinship_sorted.txt"

# Only 2nd degree and closer gets a relationship name. Everything below that is
# reported as a number and called background.
awk -v DUP="${KING_DUP}" 'BEGIN{OFS="\t"}
    NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }
            print "IID1","IID2","KINSHIP","CLASS"; next }
    { k=$(h["KINSHIP"])+0
      if      (k >  DUP)   c="duplicate/MZ"
      else if (k >  0.177) c="1st degree"
      else if (k >  0.0884) c="2nd degree"
      else                  c="below_2nd_degree_cohort_background"
      print $(h["IID1"]), $(h["IID2"]), k, c }' \
    "${R}/king_pairs.kin0" > "${R}/kinship_pairs.tsv"
log "  classes (nothing below 2nd degree is given a relationship name):"
awk -F'\t' 'NR>1{c[$4]++} END{for(k in c) printf "    %-38s %d\n", k, c[k]}' \
    "${R}/kinship_pairs.tsv" | sort -k2,2nr

# --- 2b. independent estimator: method-of-moments PI_HAT -------------------
# --genome decomposes IBD sharing using allele frequencies estimated from the
# sample, so its bias structure differs from KING-robust's. Agreement between
# the two on the top pairs is the evidence that those pairs are genuinely
# related; agreement on the elevated baseline is the evidence that the baseline
# is a property of the cohort and not an artefact of one estimator.
${PLINK1} --bfile "${W}/king_in" --genome --allow-no-sex \
    --threads "${THREADS}" --out "${R}/genome" || \
    log "  WARNING: --genome failed; the cross-check will be missing"

if [[ -s "${R}/genome.genome" ]]; then
    { printf 'IID1\tIID2\tKINSHIP\tPI_HAT\tPI_HAT_over_2\n'
      # `nm` for the header token, `gh` for the .genome header map. The previous
      # revision used `g` for both, and gawk aborts with "attempt to use scalar
      # g as an array" the moment it reaches the second file.
      awk 'FNR==NR { if (FNR==1) { for(i=1;i<=NF;i++){ nm=$i; sub(/^#/,"",nm); h[nm]=i }; next }
                     a=$(h["IID1"]); b=$(h["IID2"])
                     key=(a<b)? a"\t"b : b"\t"a
                     k[key]=$(h["KINSHIP"]); next }
           FNR==1 { for(i=1;i<=NF;i++){ nm=$i; sub(/^#/,"",nm); gh[nm]=i }; next }
           { a=$(gh["IID1"]); b=$(gh["IID2"])
             key=(a<b)? a"\t"b : b"\t"a
             if (key in k) printf "%s\t%s\t%s\t%.6f\n", key, k[key], $(gh["PI_HAT"]), $(gh["PI_HAT"])/2 }' \
          "${R}/king_pairs.kin0" "${R}/genome.genome"
    } > "${R}/kinship_vs_pihat.tsv"
    awk -F'\t' 'NR>1 { n++; sx+=$3; sy+=$5; sxx+=$3*$3; syy+=$5*$5; sxy+=$3*$5 }
        END { if (n<3) exit
              num=sxy-sx*sy/n; den=sqrt((sxx-sx*sx/n)*(syy-sy*sy/n))
              printf "  KING kinship vs PI_HAT/2 over %d pairs: r = %.4f\n", n,
                     (den>0 ? num/den : 0)
              printf "  means: KING %.5f, PI_HAT/2 %.5f\n", sx/n, sy/n }' \
        "${R}/kinship_vs_pihat.tsv"
fi

# --- 2c. same people, marker set ascertained across populations ------------
# The control for the paragraph above. If the elevated within-cohort kinship
# were real genealogy it would not depend on which markers are used; if it is
# the estimator responding to the cohort's own heterozygosity, it moves.
GA="${WORK_DIR}/groupA_kz_1000g_sgdp/king_in"
if [[ -s "${GA}.bed" ]]; then
    awk -F'\t' 'NR>1 && $3=="KAZ"{print $1}' "${META_DIR}/pop_master.tsv" > "${W}/kz.ids"
    mk_id_list "${GA}.fam" "${W}/kz_on_groupA.id" "${W}/kz.ids"
    ${PLINK2} --bfile "${GA}" --keep "${W}/kz_on_groupA.id" --make-king-table \
        --threads "${THREADS}" --out "${R}/king_pairs_groupA_markers" \
      && awk 'NR==1{ for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; next }
              { n++; s+=$(h["KINSHIP"]); if ($(h["KINSHIP"])>0.0442) hi++ }
              END { printf "  same 219 people on the groupA marker set: mean kinship %.5f, %d of %d pairs > 0.0442\n",
                           s/n, hi+0, n }' "${R}/king_pairs_groupA_markers.kin0"
    awk 'NR==1{ for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; next }
        { n++; s+=$(h["KINSHIP"]); if ($(h["KINSHIP"])>0.0442) hi++ }
        END { printf "  cohort-only marker set:                       mean kinship %.5f, %d of %d pairs > 0.0442\n",
                     s/n, hi+0, n }' "${R}/king_pairs.kin0"
else
    log "  groupA king_in not found; the marker-ascertainment control is skipped"
fi

# --- 2d. unrelated set ------------------------------------------------------
${PLINK2} --bfile "${W}/king_in" --king-cutoff "${KING_CUTOFF}" \
    --threads "${THREADS}" --out "${W}/rel"
norm_id_list "${W}/rel.king.cutoff.in.id" "${KZP}.fam" "${R}/unrelated.keep"
log "  retained at kinship < ${KING_CUTOFF}: $(nlines "${R}/unrelated.keep") of $(nlines "${KZP}.fam")"
log "  With a cohort median kinship near 0.003 this cutoff means what it says on"
log "  the tin: genuine 2nd-degree relatives, not a threshold on a noisy baseline."

{ printf 'Tribe\tn_labelled\tn_unrelated\n'
  awk -F'\t' 'FNR==NR{k[$2]=1; next}
      $2=="" || $2=="NA" {next}
      {tot[$2]++; if($1 in k) u[$2]++}
      END{for(t in tot) printf "%s\t%d\t%d\n", t, tot[t], (t in u)?u[t]:0}' \
      "${R}/unrelated.keep" "${W}/tribes_in_panel.tsv" | sort -k2,2nr; } \
  > "${R}/tribe_after_relatedness.tsv"
head -8 "${R}/tribe_after_relatedness.tsv" | sed 's/^/    /'

# ---------------------------------------------------------------------------
log "=== [3/7] marker sets ==="
# ---------------------------------------------------------------------------
# MAF is computed inside the cohort, so 0.05 here means common among Kazakhs.
${PLINK2} --bfile "${KZP}" --keep "${R}/unrelated.keep" \
    --maf "${MAF}" --geno "${GENO}" --exclude range "${HIGHLD}" \
    --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_R2}" \
    --threads "${THREADS}" --out "${W}/pca_prune"
${PLINK2} --bfile "${KZP}" --keep "${R}/unrelated.keep" \
    --extract "${W}/pca_prune.prune.in" \
    --make-bed --threads "${THREADS}" --out "${W}/kz_pruned"
${PLINK2} --bfile "${KZP}" --extract "${W}/pca_prune.prune.in" \
    --make-bed --threads "${THREADS}" --out "${W}/kz_pruned_all"
log "  pruned: $(nlines "${W}/kz_pruned.fam") individuals, $(nlines "${W}/kz_pruned.bim") SNPs"

# Fst set: unrelated, missingness-filtered, deliberately NOT MAF-filtered and
# not pruned.
${PLINK2} --bfile "${KZP}" --keep "${R}/unrelated.keep" --geno "${GENO}" \
    --make-bed --threads "${THREADS}" --out "${W}/kz_fst"
log "  Fst set: $(nlines "${W}/kz_fst.bim") variants"

# ---------------------------------------------------------------------------
log "=== [4/7] PCA ==="
# ---------------------------------------------------------------------------
mkdir -p "${R}/pca"
${PLINK2} --bfile "${W}/kz_pruned" --pca "${N_PC}" \
    --threads "${THREADS}" --out "${R}/pca/pca_unrelated"
${PLINK2} --bfile "${W}/kz_pruned_all" --pca "${N_PC}" \
    --threads "${THREADS}" --out "${R}/pca/pca_all"

# Tribe labels are joined in place of pop_master, which only knows "KAZ".
{ printf 'IID\tNation\tPopulation\tSource\tPlatform\n'
  awk -F'\t' 'NR>1{print $1"\t"$2"\tKAZ\tKZ219\twgs"}' "${W}/Tribe.labels.tsv"
} > "${W}/tribe_meta.tsv"

for SET in pca_unrelated pca_all; do
    NS=$(nbody "${R}/pca/${SET}.eigenvec")
    ${PYTHON} "${TABLES}" pca-table --eigenvec "${R}/pca/${SET}.eigenvec" \
        --meta "${W}/tribe_meta.tsv" --npc "${N_PC}" \
        --eigenval "${R}/pca/${SET}.eigenval" --n-samples "${NS}" \
        --varexp-out "${R}/pca/${SET}.varexp.tsv" --out "${R}/pca/${SET}.tsv"
done
head -6 "${R}/pca/pca_unrelated.varexp.tsv" | sed 's/^/    /'
log "  A flat eigenvalue spectrum here is the expected result, not a failed run."

# ---------------------------------------------------------------------------
log "=== [5/7] is any axis associated with tribe or zhuz? ==="
# ---------------------------------------------------------------------------
# Permutation instead of eyeballing the scatter plot. Nothing is being fitted,
# so this cannot manufacture a result; it can only bound one.
for LEVEL in Tribe Zhuz; do
    ${PYTHON} "${TABLES}" perm-assoc \
        --pc-table "${R}/pca/pca_unrelated.tsv" \
        --labels "${W}/${LEVEL}.labels.tsv" --label-col "${LEVEL}" \
        --npc 10 --min-n "${TRIBE_MIN_N}" --nperm "${PERM_N}" \
        --out "${R}/pc_assoc_${LEVEL}.tsv" \
        --out-levels "${R}/pc_assoc_${LEVEL}_levels.tsv"
    log "  ${LEVEL}:"
    column -t -s $'\t' "${R}/pc_assoc_${LEVEL}.tsv" | sed 's/^/    /'
done

# Mean kinship within a tribe versus between tribes. The most direct autosomal
# test of exogamy available, and it uses every pair rather than a two-group
# summary. Under strict tribal exogamy the difference should be ~0.
for LEVEL in Tribe Zhuz; do
    ${PYTHON} "${TABLES}" kin-by-group --kin "${R}/king_pairs.kin0" \
        --labels "${W}/${LEVEL}.labels.tsv" --label-col "${LEVEL}" \
        --min-n "${TRIBE_MIN_N}" --nperm "${PERM_N}" \
        --out "${R}/kinship_by_${LEVEL}.tsv" \
        --out-levels "${R}/kinship_by_${LEVEL}_levels.tsv" \
      || log "  kin-by-group failed for ${LEVEL}"
    [[ -s "${R}/kinship_by_${LEVEL}.tsv" ]] && \
        { log "  kinship within vs between ${LEVEL}:"
          column -t -s $'\t' "${R}/kinship_by_${LEVEL}.tsv" | sed 's/^/    /'; }
done

# ---------------------------------------------------------------------------
log "=== [6/7] pairwise Fst, levels with n >= ${TRIBE_MIN_N} ==="
# ---------------------------------------------------------------------------
fst_level() {   # fst_level <name>
    local NAME=$1
    local D="${R}/fst/${NAME}"; mkdir -p "${D}"
    awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR>1 && ($1 in keep){n[$2]++}
        END{for(k in n) printf "%s\t%d\n", k, n[k]}' \
        <(fam_iids "${W}/kz_fst.fam") "${W}/${NAME}.cat" | sort -k2,2nr > "${D}/n.body"
    { printf 'POP\tN\n'; cat "${D}/n.body"; } > "${D}/pop_n.tsv"
    awk -F'\t' -v m="${TRIBE_MIN_N}" '$2>=m{print $1}' "${D}/n.body" > "${D}/keep_levels.txt"
    rm -f "${D}/n.body"
    local NK; NK=$(nlines "${D}/keep_levels.txt")
    log "  ${NAME}: ${NK} level(s) reach n>=${TRIBE_MIN_N}"
    { printf 'level\tn\n'
      awk -F'\t' -v m="${TRIBE_MIN_N}" 'NR>1 && $2<m {print $1"\t"$2}' "${D}/pop_n.tsv"
    } > "${D}/excluded.tsv"
    if [[ "${NK}" -lt 2 ]]; then
        log "    fewer than two levels clear the threshold, no Fst computed"
        return
    fi
    awk -F'\t' 'FNR==NR{ok[$1]=1; next} FNR==1{print; next} ($2 in ok)' \
        "${D}/keep_levels.txt" "${W}/${NAME}.cat" > "${D}/pop.cat"
    ${PLINK2} --bfile "${W}/kz_fst" --no-psam-pheno --pheno "${D}/pop.cat" \
        --fst "${NAME}" method=hudson blocksize="${FST_BLOCKSIZE}" \
        --threads "${THREADS}" --out "${D}/fst"
    require_file "${D}/fst.fst.summary"
    ${PYTHON} "${TABLES}" fst-tables --summary "${D}/fst.fst.summary" \
        --pop-n "${D}/pop_n.tsv" \
        --out-pairs "${D}/fst_pairs.tsv" --out-matrix "${D}/fst_matrix.tsv"
    column -t -s $'\t' "${D}/fst_pairs.tsv" | sed 's/^/    /'

    # The detection limit. A null Fst is only interpretable next to the smallest
    # value the design could have resolved, and that is 1.96 x the jackknife SE.
    awk -F'\t' 'NR>1 { printf "    %s vs %s: detectable Fst >= %.2e (1.96 x SE); observed %s, CI covers zero: %s\n",
                              $1, $2, 1.96*$6, $5, $9 }' "${D}/fst_pairs.tsv"
}
fst_level Tribe
fst_level Zhuz

# ---------------------------------------------------------------------------
log "=== [7/7] ADMIXTURE K=${KZ_K_MIN}..${KZ_K_MAX}, ${KZ_N_REPS} seeds ==="
# ---------------------------------------------------------------------------
# K starts at 1 on purpose. Within one endogamous population the question is
# whether any K beats K=1, and a run starting at K=2 cannot answer it. If CV is
# minimised at K=1, that is the finding: no autosomal substructure resolvable at
# this sample size, which is what tribal exogamy predicts.
A="${R}/admixture"; mkdir -p "${A}"
admix_cache_guard "${W}/kz_pruned" "${A}"
STEM=$(basename "${W}/kz_pruned")
awk '{print $2}' "${W}/kz_pruned.fam" > "${A}/individuals.txt"
printf '%s\n' "${STEM}" > "${A}/stem.txt"

for K in $(seq "${KZ_K_MIN}" "${KZ_K_MAX}"); do
    for REP in $(seq 1 "${KZ_N_REPS}"); do
        SEED=$(( SEED_BASE + 90000 + K*100 + REP ))
        RD="${A}/runs/K${K}_r${REP}"
        [[ -f "${RD}/${STEM}.${K}.Q" ]] && continue
        mkdir -p "${RD}"
        for E in bed bim fam; do ln -sf "${W}/kz_pruned.${E}" "${RD}/${STEM}.${E}"; done
        while [[ "$(jobs -rp | wc -l)" -ge "${ADMIX_JOBS}" ]]; do wait -n || true; done
        ( cd "${RD}" && ${ADMIXTURE} --cv="${CV_FOLDS}" -j"${ADMIX_THREADS}" \
            -s "${SEED}" "${STEM}.bed" "${K}" > log.txt 2>&1 ) &
    done
done
wait || true

printf 'K\tREP\tSEED\tCV\tLOGLIK\n' > "${A}/cv_raw.tsv"
for K in $(seq "${KZ_K_MIN}" "${KZ_K_MAX}"); do
    for REP in $(seq 1 "${KZ_N_REPS}"); do
        RD="${A}/runs/K${K}_r${REP}"; SEED=$(( SEED_BASE + 90000 + K*100 + REP ))
        [[ -f "${RD}/log.txt" ]] || continue
        [[ -f "${RD}/${STEM}.${K}.Q" ]] || { log "  FAILED K=${K} rep=${REP}"; continue; }
        CV=$(grep "CV error" "${RD}/log.txt" | tail -1 | awk '{print $NF}')
        LL=$(grep "^Loglikelihood" "${RD}/log.txt" | tail -1 | awk '{print $NF}')
        printf '%s\t%s\t%s\t%s\t%s\n' "${K}" "${REP}" "${SEED}" "${CV:-NA}" "${LL:-NA}" >> "${A}/cv_raw.tsv"
    done
done

{ printf 'K\tREP\tQ_PATH\n'
  for K in $(seq "${KZ_K_MIN}" "${KZ_K_MAX}"); do
    for REP in $(seq 1 "${KZ_N_REPS}"); do
      Q="${A}/runs/K${K}_r${REP}/${STEM}.${K}.Q"
      [[ -f "${Q}" ]] && printf '%s\t%s\t%s\n' "${K}" "${REP}" "${Q}"
    done
  done; } > "${A}/runs.tsv"

awk -F'\t' 'NR>1 && $4!="NA" {n[$1]++; s[$1]+=$4; ss[$1]+=$4*$4}
    END{ printf "K\tn_runs\tCV_mean\tCV_sd\n"
         for (k in n){ m=s[k]/n[k]; v=ss[k]/n[k]-m*m; if(v<0)v=0
                       printf "%s\t%d\t%.5f\t%.5f\n", k, n[k], m, sqrt(v) } }' \
    "${A}/cv_raw.tsv" | { read -r h; echo "$h"; sort -n; } > "${A}/cv_summary.tsv"
log "  CV summary:"; sed 's/^/    /' "${A}/cv_summary.tsv"

BEST=$(awk -F'\t' 'NR>1{if(m==""||$3<m){m=$3;k=$1}} END{print k}' "${A}/cv_summary.tsv")
log "  lowest mean CV error at K=${BEST}"
[[ "${BEST}" == "1" ]] && \
    log "  K=1 wins: no autosomal substructure resolvable at n=$(nlines "${W}/kz_pruned.fam")"
# The alignment now runs WHETHER OR NOT K=1 won. Skipping it when it did was
# exactly backwards: with K=1 selected, the reviewer's next question is "show
# me what K=2 looks like", and the answer - components that line up with no
# tribe - is the strongest form of the null. The runs are already on disk, so
# this costs seconds.
${PYTHON} "${ALIGN}" sweep --admix-dir "${A}" \
    --k-min 2 --k-max "${KZ_K_MAX}" --threshold "${ADMIX_SIM_THRESHOLD}" \
    --individuals "${A}/individuals.txt" --meta "${W}/tribe_meta.tsv" \
    --out-prefix "${A}/kz" || log "  alignment failed"

log "=== 09 done ==="
log "Report together: Fst with its interval AND the detection limit; the"
log "permutation p-values for PC-vs-tribe; the within/between kinship contrast;"
log "and the CV curve from K=1. Four independent statements of the same"
log "negative result are what make it a finding rather than an absence."
