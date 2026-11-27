#!/bin/bash
set -euo pipefail

# ============================================================================
# kz233_relatedness.sh
#
# Reproduces the 233 -> 219 relatedness step of the cohort derivation, using
# parameters identical to those that built the 219-genome panel in
# 01_prepare.sh (harmonise) and 03_relatedness.sh (KING). Nothing here changes
# any analysis file; it writes only to OUT_DIR.
#
# Every flag below is taken from the deployed pipeline:
#   VCF import      KZ_VCF_FILTERS=(--var-filter --vcf-half-call missing)
#                   KZ_GQ_MIN=20            -> --vcf-min-gq 20
#   variant pass    --chr 1-22 --snps-only just-acgt --max-alleles 2
#                   --min-alleles 2 --set-all-var-ids '@:#:$r:$a'
#                   --new-id-max-allele-len 10 truncate --rm-dup force-first
#   position dedup  one record per CHR:POS, bim order (deterministic)
#   lenient pass    GENO_PRE=0.10
#   het outliers    HET_MAF=0.05, HET_SD=3, HET_MODE=excess (drop excess only)
#   sample pass     MIND_KZ=0.05
#   final pass      GENO_PANEL=0.05
#   KING pruning    --maf 0.05 --geno 0.05, high-LD excluded,
#                   --indep-pairwise 200kb 1 0.2
#   KING            --king-table-filter 0.0442, --king-cutoff 0.0884
# ============================================================================

PROJECT=${PROJECT:?set PROJECT to the popgen project root}
VCF233=${VCF233:?set VCF233 to the 233-genome multisample VCF}
HIGHLD="${PROJECT}/results/work/highld_hg19.txt"
PANEL219="${PROJECT}/results/work/panels/kz"

OUT_DIR="${PROJECT}/results/kz233_relatedness"
P="${OUT_DIR}/kz233"
THREADS=16
MEM_MB=131072

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"
exec > >(tee -a "${OUT_DIR}/kz233_relatedness_log.txt") 2>&1
echo "=== started $(date -Is) ==="

for f in "${VCF233}" "${HIGHLD}" "${PANEL219}.fam"; do
    [[ -e "$f" ]] || { echo "missing: $f"; exit 1; }
done

# --- 1. variant pass, canonical IDs, ID-level dedup -------------------------
plink2 --vcf "${VCF233}" --var-filter --vcf-half-call missing --vcf-min-gq 20 \
    --chr 1-22 --allow-extra-chr \
    --snps-only just-acgt --max-alleles 2 --min-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10 truncate \
    --rm-dup force-first list \
    --make-bed --threads "${THREADS}" --memory "${MEM_MB}" --out "${P}_s1"

# --- 2. collapse multi-record positions -------------------------------------
awk '{k=$1":"$4; if (!(k in seen)) {seen[k]=1; print $2}}' "${P}_s1.bim" \
    > "${P}.pos_unique.snplist"

# --- 3. lenient variant pass, before --mind ---------------------------------
plink2 --bfile "${P}_s1" --extract "${P}.pos_unique.snplist" --geno 0.10 \
    --make-bed --threads "${THREADS}" --memory "${MEM_MB}" --out "${P}_s2"

# --- 4. heterozygosity outliers, excess only --------------------------------
plink2 --bfile "${P}_s2" --maf 0.05 --het \
    --threads "${THREADS}" --out "${P}_qc"

awk -v SD=3 'NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; next }
     { n++; F[n]=$(h["F"]); I[n]=$(h["IID"]); s+=F[n]; ss+=F[n]*F[n] }
     END { if (n<2) exit 0
           m=s/n; v=ss/n-m*m; if (v<0) v=0; sd=sqrt(v); if (sd<=0) exit 0
           for(i=1;i<=n;i++) if (F[i] < m-SD*sd) print I[i] }' \
    "${P}_qc.het" > "${P}.het_drop.txt"
echo "het outliers dropped (excess): $(wc -l < "${P}.het_drop.txt")"

RM=()
[[ -s "${P}.het_drop.txt" ]] && {
    awk 'NR==FNR{d[$1]=1; next} ($2 in d){print $1"\t"$2}' \
        "${P}.het_drop.txt" "${P}_s2.fam" > "${P}.drop.id"
    RM=(--remove "${P}.drop.id")
}

# --- 5. sample pass, then final variant pass --------------------------------
plink2 --bfile "${P}_s2" "${RM[@]+"${RM[@]}"}" --mind 0.05 \
    --make-bed --threads "${THREADS}" --memory "${MEM_MB}" --out "${P}_s3"
[[ -f "${P}_s3.mindrem.id" ]] && \
    echo "dropped by --mind: $(awk '!/^#/' "${P}_s3.mindrem.id" | wc -l)"

plink2 --bfile "${P}_s3" --geno 0.05 \
    --make-bed --threads "${THREADS}" --memory "${MEM_MB}" --out "${P}"
rm -f "${P}_s1."{bed,bim,fam} "${P}_s2."{bed,bim,fam} "${P}_s3."{bed,bim,fam}
echo "panel: $(wc -l < "${P}.fam") samples, $(wc -l < "${P}.bim") variants"

# --- 6. KING ----------------------------------------------------------------
plink2 --bfile "${P}" --maf 0.05 --geno 0.05 \
    --exclude range "${HIGHLD}" \
    --indep-pairwise 200kb 1 0.2 \
    --threads "${THREADS}" --out "${P}_prune"
NMARK=$(wc -l < "${P}_prune.prune.in")
echo "KING markers: ${NMARK}"

plink2 --bfile "${P}" --extract "${P}_prune.prune.in" \
    --make-king-table --king-table-filter 0.0442 \
    --threads "${THREADS}" --out "${P}_king"

plink2 --bfile "${P}" --extract "${P}_prune.prune.in" \
    --king-cutoff 0.0884 \
    --threads "${THREADS}" --out "${P}_cut"

# --- 7. pair table, degree bands identical to 03_relatedness.sh -------------
awk -F'\t' 'NR==1{ for(i=1;i<=NF;i++){g=$i;sub(/^#/,"",g);h[g]=i}
                   print "IID1","IID2","KINSHIP","DEGREE"; next }
     { k=$(h["KINSHIP"])
       d = k>0.354 ? "duplicate/MZ" : k>0.177 ? "1st degree" :
           k>0.0884 ? "2nd degree"  : k>0.0442 ? "3rd degree" : "4th degree"
       print $(h["IID1"]), $(h["IID2"]), k, d }' OFS='\t' \
    "${P}_king.kin0" > "${OUT_DIR}/kz233_pairs.tsv"

echo "--- pairs by degree ---"
awk -F'\t' 'NR>1{c[$4]++} END{for(d in c) printf "%-14s %d\n", d, c[d]}' \
    "${OUT_DIR}/kz233_pairs.tsv"

# --- 8. does the KING-retained set equal the analysed 219? ------------------
awk '!/^#/{print $2}' "${P}_cut.king.cutoff.in.id"  | sort > "${OUT_DIR}/keep233.txt"
awk '!/^#/{print $2}' "${P}_cut.king.cutoff.out.id" | sort > "${OUT_DIR}/removed233.txt"
awk '{print $2}' "${PANEL219}.fam"                  | sort > "${OUT_DIR}/actual219.txt"

echo "--- reconciliation ---"
echo "retained by KING on 233 : $(wc -l < "${OUT_DIR}/keep233.txt")"
echo "removed by KING on 233  : $(wc -l < "${OUT_DIR}/removed233.txt")"
echo "analysed cohort         : $(wc -l < "${OUT_DIR}/actual219.txt")"
echo "in KING set, not analysed:"
comm -23 "${OUT_DIR}/keep233.txt" "${OUT_DIR}/actual219.txt" | sed 's/^/    /'
echo "analysed, not in KING set:"
comm -13 "${OUT_DIR}/keep233.txt" "${OUT_DIR}/actual219.txt" | sed 's/^/    /'

# --- 9. threshold sensitivity, for the supplementary table -----------------
{ printf 'cutoff\tn_retained\tn_removed\n'
  for C in 0.0442 0.0884 0.13; do
      plink2 --bfile "${P}" --extract "${P}_prune.prune.in" --king-cutoff "${C}" \
          --threads "${THREADS}" --out "${P}_cut_${C}" >/dev/null 2>&1 || continue
      NK=$(awk '!/^#/' "${P}_cut_${C}.king.cutoff.in.id" | wc -l)
      printf '%s\t%s\t%s\n' "${C}" "${NK}" "$(( $(wc -l < "${P}.fam") - NK ))"
  done
} > "${OUT_DIR}/kz233_cutoff_sensitivity.tsv"
echo "--- cutoff sensitivity ---"
cat "${OUT_DIR}/kz233_cutoff_sensitivity.tsv"

echo "=== done $(date -Is) ==="
echo "send back: kz233_pairs.tsv, kz233_cutoff_sensitivity.tsv, removed233.txt,"
echo "           and the reconciliation block above"
