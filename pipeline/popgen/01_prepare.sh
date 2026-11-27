#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 01_prepare.sh — harmonise each source panel, per-sample QC, population labels.
#
# Per panel: autosomes, biallelic ACGT SNPs, CHR:POS:REF:ALT IDs, duplicate IDs
# removed, one record per position, ambiguous strand SNPs removed, call-rate and
# excess-heterozygosity outliers dropped.
#
# Five choices worth stating in Methods.
#
#   Duplicate variant IDs. --set-all-var-ids writes CHR:POS:REF:ALT, and a panel
#   can still carry two records with the same ID. An awk list of unique IDs does
#   not remove them, because --extract selects by ID and pulls both records back
#   in. --rm-dup force-first resolves it inside plink2, before any ID list is
#   written.
#
#   One record per position. After deduplication by ID, 1:100:A:G and 1:100:A:T
#   both survive, and a position represented twice inside one panel could then
#   satisfy a "present in all panels" intersection on its own. Positions are
#   collapsed to a single record here so the intersection in 02 counts panels.
#   The cost is visible in 02: when two panels keep different representations of
#   the same multi-allelic position, the merge cannot reconcile them and the
#   position is dropped. 02 now counts and classifies those.
#
#   Ambiguous SNPs. A/T and C/G cannot be strand-resolved by allele identity, so
#   they never appear in a merge .missnp list and --flip cannot reach them. Left
#   in, they merge on whichever strand each panel used and put a false
#   between-panel difference into PC1-PC3 and into every Fst value.
#
#   Heterozygosity filter direction. See HET_MODE in config.sh. Only excess
#   heterozygosity is a quality signal; a deficit is consanguinity, and this
#   cohort is the one being measured for it.
#
#   FILTER ORDER. --mind and --geno are applied in three passes, not one:
#   a lenient variant pass (GENO_PRE), then the sample pass (--mind), then the
#   final variant pass. This is what makes genotype-level GQ masking usable on
#   the cohort VCF. Masked calls are concentrated in rare, poorly covered sites;
#   evaluating --mind before removing those sites made per-sample missingness
#   look catastrophic and deleted all 219 samples, which is why the previous
#   revision concluded GQ filtering was impossible and shipped with none.
# ============================================================================

need "${PLINK1}"; need "${PLINK2}"
require_file "${MF_SUPERPOP}"

PANELS="${WORK_DIR}/panels"; mkdir -p "${PANELS}"

# LABELS_ONLY=1 skips every harmonise() call and rebuilds only the metadata:
# pop_master.tsv, the KEEP_NATIONS resolution, the QC summary and the
# validation. The existing panels are reused as they are.
#
# It exists because a metadata-only defect - a mislabelled population, a wrong
# KEEP_NATIONS token - otherwise costs a two-hour rerun of five panel imports
# that would produce byte-identical output. Only 06 and 07 group by Nation, so
# after this the useful rerun is those two, five minutes.
#
#   LABELS_ONLY=1 bash 01_prepare.sh && bash 06_fst.sh && bash 07_roh.sh
#
# Do NOT use it after changing anything that reaches the genotypes: KZ_GQ_MIN,
# GENO_PRE, MIND_*, HET_*, DROP_AMBIGUOUS or any dataset path. It reuses panels
# built under the OLD settings and will silently mix revisions.
LABELS_ONLY="${LABELS_ONLY:-0}"
[[ "${LABELS_ONLY}" -eq 1 ]] && \
    log "LABELS_ONLY=1: reusing the existing panels, rebuilding metadata only"

HARMONISE_EXTRA=()
# Absolute cap on how many samples --mind may remove from the next panel. Empty
# means "use the MIND_MAX_FRAC_REMOVED fraction instead". Set by the caller and
# reset by the caller, like HARMONISE_EXTRA.
HARMONISE_MAX_MIND_REMOVED=""

# ---------------------------------------------------------------------------
# harmonise <in_spec> <out_prefix> <mind> <platform> [geno_panel]
# extra import-time flags go in HARMONISE_EXTRA, reset by the caller
# ---------------------------------------------------------------------------
harmonise() {
    local IN_SPEC=$1 OUT=$2 MIND=$3 PLATFORM=$4 GENO_PANEL=${5:-0.05}

    if [[ "${LABELS_ONLY}" -eq 1 ]]; then
        require_file "${OUT}.bed"
        require_file "${OUT}.mergeable.snplist"
        log "harmonise -> $(basename "${OUT}"): reused ($(nlines "${OUT}.fam") samples, \
$(nlines "${OUT}.bim") variants)"
        return 0
    fi

    # The lenient pass must never be STRICTER than the final one. GENO_PRE=0.10
    # is lenient beside GENO=0.05, but for AADR the final threshold is
    # GENO_AADR=0.90 - because that panel's missingness is coverage, not quality
    # - and a 0.10 pre-pass in front of it took the panel from 1,233,013
    # variants to 532 on the 2026-09-02 run. max(), not the configured value.
    local GENO_PRE_EFF
    GENO_PRE_EFF=$(awk -v a="${GENO_PRE}" -v b="${GENO_PANEL}" 'BEGIN{print (a>b)?a:b}')
    log "harmonise -> $(basename "${OUT}")  (${PLATFORM}, mind=${MIND}, geno=${GENO_PANEL}, geno_pre=${GENO_PRE_EFF})"

    # --- 1. autosomes, biallelic ACGT SNPs, canonical IDs, ID-level dedup ---
    # shellcheck disable=SC2086
    ${PLINK2} ${IN_SPEC} "${HARMONISE_EXTRA[@]+"${HARMONISE_EXTRA[@]}"}" \
        --chr 1-22 --allow-extra-chr \
        --snps-only just-acgt --max-alleles 2 --min-alleles 2 \
        --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10 truncate \
        --rm-dup force-first list \
        --make-bed --threads "${THREADS}" --memory "${MEM_MB}" \
        --out "${OUT}_s1"
    if [[ -s "${OUT}_s1.rmdup.list" ]]; then
        log "  duplicate variant IDs removed: $(nlines "${OUT}_s1.rmdup.list")"
    else
        log "  duplicate variant IDs removed: 0"
    fi

    # --- 2. collapse multi-record positions ---------------------------------
    # bim is position-sorted, so the record kept is deterministic. IDs are
    # unique by now, so --extract is exact.
    awk '{k=$1":"$4; if (!(k in seen)) {seen[k]=1; print $2}}' \
        "${OUT}_s1.bim" > "${OUT}.pos_unique.snplist"
    log "  multi-record positions collapsed: $(( $(nlines "${OUT}_s1.bim") - $(nlines "${OUT}.pos_unique.snplist") ))"

    # --- 3. ambiguous strand SNPs -------------------------------------------
    if [[ "${DROP_AMBIGUOUS}" -eq 1 ]]; then
        awk '{ a=toupper($5); b=toupper($6);
               if ((a=="A"&&b=="T")||(a=="T"&&b=="A")||
                   (a=="C"&&b=="G")||(a=="G"&&b=="C")) print $2 }' \
            "${OUT}_s1.bim" > "${OUT}.ambiguous.snplist"
        local NAMB; NAMB=$(nlines "${OUT}.ambiguous.snplist")
        log "  ambiguous A/T + C/G: ${NAMB} (kept in the panel, excluded from merges)"
        # A source panel legitimately reaching zero here has already been
        # strand-filtered upstream (hgdp650.qc.hg19 has). Saying so in the log
        # is the difference between a documented property and a filter that
        # looks like it silently did nothing.
        [[ "${NAMB}" -eq 0 ]] && \
            log "  NOTE: zero ambiguous sites - this source panel was already strand-filtered upstream"
    else
        : > "${OUT}.ambiguous.snplist"
    fi

    # --- 4. lenient variant pass -------------------------------------------
    # Runs BEFORE --mind. Removes the sites that GQ masking emptied, so the
    # per-sample missingness that --mind sees is a property of the sample and
    # not of the rare tail of the callset.
    # NOTE: the ambiguous list is NOT excluded here. It is applied by 02 and 08
    # when they build a merge, because that is the only place strand ambiguity
    # exists. Removing it at panel-build time propagated a merge-specific filter
    # into the Kazakh-only analyses.
    ${PLINK2} --bfile "${OUT}_s1" \
        --extract "${OUT}.pos_unique.snplist" \
        --geno "${GENO_PRE_EFF}" \
        --make-bed --threads "${THREADS}" --memory "${MEM_MB}" --out "${OUT}_s2"
    rm -f "${OUT}_s1."{bed,bim,fam}
    log "  after geno_pre: $(nlines "${OUT}_s2.bim") variants"

    # sanity: IDs must be unique now, or 02 will merge the wrong records
    local NDUP
    NDUP=$(cut -f2 "${OUT}_s2.bim" | sort | uniq -d | wc -l)
    [[ "${NDUP}" -eq 0 ]] || die "$(basename "${OUT}"): ${NDUP} duplicate variant IDs survived --rm-dup"

    # --- 5. QC metrics, computed before any sample is removed ---------------
    # --het on common variants: F computed over the full spectrum is dominated
    # by rare sites and by their per-panel ascertainment, which makes the panels
    # not comparable to each other.
    ${PLINK2} --bfile "${OUT}_s2" --missing sample-only \
        --threads "${THREADS}" --out "${OUT}_qc"

    : > "${OUT}.drop_samples.id"
    printf 'IID\tF\tpanel_mean\tpanel_sd\tclass\n' > "${OUT}.het_drop.tsv"
    printf 'IID\tF\tpanel_mean\tpanel_sd\tclass\n' > "${OUT}.het_flag.tsv"

    # plink2 refuses --het below 50 samples: it has to impute allele frequencies
    # from the panel itself and will not do so from that few people. The AADR
    # panel is 46 individuals and pseudo-haploid anyway, so F is ~1 for every one
    # of them and there is no outlier to find. Skipping is the correct answer.
    local NSAMP; NSAMP=$(nlines "${OUT}_s2.fam")
    if [[ "${NSAMP}" -lt 50 ]]; then
        log "  --het skipped: ${NSAMP} samples, below plink2's 50-sample minimum"
    else
        ${PLINK2} --bfile "${OUT}_s2" --maf "${HET_MAF}" --het \
            --threads "${THREADS}" --out "${OUT}_qc"

        # F below the mean = excess heterozygosity = contamination -> dropped.
        # F above the mean = heterozygosity deficit = inbreeding -> flagged only.
        # Columns resolved by header name; plink2 has changed .het column order.
        awk -v SD="${HET_SD}" -v MODE="${HET_MODE}" \
            -v DROPF="${OUT}.het_drop.tsv" -v FLAGF="${OUT}.het_flag.tsv" '
            BEGIN { printf "IID\tF\tpanel_mean\tpanel_sd\tclass\n" > DROPF
                    printf "IID\tF\tpanel_mean\tpanel_sd\tclass\n" > FLAGF }
            NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }
                    if (!("F" in h) || !("IID" in h)) { print "HEADER" > "/dev/stderr"; exit 1 }
                    next }
            { n++; F[n]=$(h["F"]); I[n]=$(h["IID"]); s+=F[n]; ss+=F[n]*F[n] }
            END { if (n<2) exit 0
                  m=s/n; v=ss/n-m*m; if (v<0) v=0; sd=sqrt(v)
                  if (sd<=0) exit 0
                  for(i=1;i<=n;i++) {
                      if (F[i] < m-SD*sd)
                          printf "%s\t%.6f\t%.6f\t%.6f\t%s\n", I[i],F[i],m,sd,"excess_heterozygosity" > DROPF
                      else if (F[i] > m+SD*sd) {
                          if (MODE=="both")
                              printf "%s\t%.6f\t%.6f\t%.6f\t%s\n", I[i],F[i],m,sd,"heterozygosity_deficit" > DROPF
                          else
                              printf "%s\t%.6f\t%.6f\t%.6f\t%s\n", I[i],F[i],m,sd,"heterozygosity_deficit_retained" > FLAGF } } }' \
            "${OUT}_qc.het"
        log "  het outliers dropped (excess): $(nbody "${OUT}.het_drop.tsv")"
        log "  het outliers flagged, retained (deficit): $(nbody "${OUT}.het_flag.tsv")"
        [[ $(nbody "${OUT}.het_flag.tsv") -gt 0 ]] && \
            tail -n +2 "${OUT}.het_flag.tsv" | awk '{printf "        %s  F=%s\n", $1, $2}'
    fi

    # --- 6. assemble the removal list --------------------------------------
    # Manual exclusions for this panel only. The file is <panel><TAB><IID>;
    # 23 of the listed IIDs are in both 1000G and SGDP and only the SGDP copy
    # is named.
    #
    # `awk 'NF'` and not `grep -v '^[[:space:]]*$'`: when the block emits nothing
    # at all - which is the NORMAL case for a clean panel with no excess-
    # heterozygosity outliers and no manual exclusions - grep selects no lines
    # and exits 1. Under `set -o pipefail` that fails the whole pipeline and
    # `set -e` then kills 01 silently, right after the het summary and before
    # anything downstream is written. That is exactly how the 2026-09-02 17:47
    # run died: the log simply stopped, 02 then failed on a missing
    # kz.mergeable.snplist, and neither said why. awk exits 0 on empty input.
    #
    # The list is built into a NAMED TEMPORARY FILE and only then converted.
    # The previous revision piped `cut -f2 <drop_samples.id>` into a stage whose
    # own output redirection truncated that same file. All commands in a bash
    # pipeline start concurrently, so the reader raced the truncation and lost.
    # It lost every time on the 1000G panel: 4 excess-heterozygosity samples
    # were identified, silently never removed, and stayed in every downstream
    # analysis. The assertion below makes that class of failure impossible to
    # repeat quietly.
    local PANEL TMPD
    PANEL=$(basename "${OUT}")
    TMPD="${OUT}.drop_samples.list"
    { tail -n +2 "${OUT}.het_drop.tsv" | cut -f1
      [[ -s "${MF_PANEL_EXCLUDE:-/nonexistent}" ]] && \
          awk -F'\t' -v p="${PANEL}" '$0 !~ /^#/ && $1==p {print $2}' "${MF_PANEL_EXCLUDE}"
      true; } | awk 'NF' | sort -u > "${TMPD}"

    mk_id_list "${OUT}_s2.fam" "${OUT}.drop_samples.id" "${TMPD}"

    local NHET NWANT NGOT
    NHET=$(nbody "${OUT}.het_drop.tsv")
    NWANT=$(nlines "${TMPD}")
    NGOT=$(nlines "${OUT}.drop_samples.id")
    log "  removal list: ${NWANT} requested (${NHET} heterozygosity + manual), ${NGOT} matched in this panel"
    [[ "${NGOT}" -ge "${NHET}" ]] \
        || die "${PANEL}: ${NHET} heterozygosity outliers but only ${NGOT} IDs resolved - drop list is being lost"
    rm -f "${TMPD}"

    # --- 7. sample pass, then the final variant pass -----------------------
    rm -f "${OUT}.mindrem.id"
    local RM=()
    [[ -s "${OUT}.drop_samples.id" ]] && RM=(--remove "${OUT}.drop_samples.id")
    ${PLINK2} --bfile "${OUT}_s2" "${RM[@]+"${RM[@]}"}" \
        --mind "${MIND}" \
        --make-bed --threads "${THREADS}" --memory "${MEM_MB}" --out "${OUT}_s3"

    local NMIND=0
    # awk, not `grep -vc ... || echo 0`: grep -vc PRINTS the count and STILL
    # exits 1 when that count is zero, so the `|| echo 0` fires as well and the
    # variable ends up holding two lines. Every later [[ -gt ]] on it is then a
    # syntax error.
    [[ -f "${OUT}.mindrem.id" ]] && NMIND=$(awk '!/^#/' "${OUT}.mindrem.id" | wc -l)
    log "  dropped by --mind: ${NMIND}"

    # A --mind that removes a large share of a panel is a wrong threshold, not a
    # result. Left unguarded this is exactly the failure mode that made the
    # previous revision abandon GQ filtering.
    #
    # Two regimes. Reference panels get a fraction, because losing a handful of
    # 2504 reference individuals is immaterial. The cohort gets an absolute cap,
    # because a fraction of 219 is people the paper is about, and because
    # loosening MIND_KZ to make GQ masking viable also loosens what a fractional
    # brake would catch.
    local NIN; NIN=$(nlines "${OUT}_s2.fam")
    if [[ -n "${HARMONISE_MAX_MIND_REMOVED}" ]]; then
        if [[ "${NMIND}" -gt "${HARMONISE_MAX_MIND_REMOVED}" ]]; then
            log "  --mind would remove these samples:"
            grep -v '^#' "${OUT}.mindrem.id" | sed 's/^/        /'
            die "${PANEL}: --mind ${MIND} removes ${NMIND} sample(s); the limit for this panel is \
${HARMONISE_MAX_MIND_REMOVED}. Expected zero at KZ_GQ_MIN=${KZ_GQ_MIN} - rerun gq_profile.sh and \
compare. If a sample genuinely has to go, add it to ${MF_PANEL_EXCLUDE} so the exclusion is \
recorded in the metadata rather than decided by a threshold."
        fi
    else
        awk -v a="${NMIND}" -v b="${NIN}" -v f="${MIND_MAX_FRAC_REMOVED}" \
            'BEGIN{ exit !(b>0 && a/b > f) }' && \
            die "${PANEL}: --mind ${MIND} removed ${NMIND} of ${NIN} samples (> ${MIND_MAX_FRAC_REMOVED}). \
Either the threshold is wrong for this platform, or this plink2 build applies --geno before --mind \
and GENO_PRE is not doing its job."
    fi

    ${PLINK2} --bfile "${OUT}_s3" --geno "${GENO_PANEL}" \
        --make-bed --threads "${THREADS}" --memory "${MEM_MB}" --out "${OUT}"
    rm -f "${OUT}_s2."{bed,bim,fam} "${OUT}_s3."{bed,bim,fam}

    ${PLINK2} --bfile "${OUT}" --write-snplist --threads "${THREADS}" --out "${OUT}"
    # Two lists. <panel>.snplist is everything in the panel and is what the
    # single-panel analyses use. <panel>.mergeable.snplist drops the ambiguous
    # sites and is what 02 and 08 intersect on.
    # The empty-first-file branch is NOT decoration. With an empty ambiguous
    # list, NR and FNR both start counting on the SECOND file, NR==FNR is true
    # for every line of the snplist, all of them are swallowed into a[] and
    # skipped by `next`, and the output is empty. hgdp is the one panel with
    # zero ambiguous SNPs, so it is the one panel that hit this - and groupB is
    # kz + hgdp, so 02 died on a require_file it could not explain.
    if [[ -s "${OUT}.ambiguous.snplist" ]]; then
        awk 'NR==FNR{a[$1]=1; next} !($1 in a)' \
            "${OUT}.ambiguous.snplist" "${OUT}.snplist" > "${OUT}.mergeable.snplist"
    else
        cp "${OUT}.snplist" "${OUT}.mergeable.snplist"
    fi
    local NSL NML
    NSL=$(nlines "${OUT}.snplist"); NML=$(nlines "${OUT}.mergeable.snplist")
    [[ "${NML}" -gt 0 && "${NML}" -le "${NSL}" ]] \
        || die "$(basename "${OUT}"): mergeable list has ${NML} of ${NSL} variants, which is impossible"
    log "  final: $(nlines "${OUT}.fam") samples, $(nlines "${OUT}.bim") variants"
    log "  mergeable (ambiguous removed): $(nlines "${OUT}.mergeable.snplist")"
}

# ---------------------------------------------------------------------------
log "=== [1/6] Kazakh cohort ==="
HARMONISE_EXTRA=("${KZ_VCF_FILTERS[@]}")
if [[ "${KZ_GQ_MIN}" -gt 0 ]]; then
    HARMONISE_EXTRA+=(--vcf-min-gq "${KZ_GQ_MIN}")
    log "  genotype-level filter: GQ >= ${KZ_GQ_MIN} (calls below are set missing)"
else
    log "  WARNING: KZ_GQ_MIN=0, no genotype-level quality filter on the cohort VCF"
fi
HARMONISE_MAX_MIND_REMOVED="${MIND_MAX_REMOVED_KZ}"
harmonise "--vcf ${DS_KZ_VCF}" "${PANELS}/kz" "${MIND_KZ}" wgs
HARMONISE_EXTRA=()
HARMONISE_MAX_MIND_REMOVED=""
N_KZ=$(nlines "${PANELS}/kz.fam")
[[ "${N_KZ}" -eq "${N_KZ_EXPECTED}" ]] \
    || log "WARNING: Kazakh n=${N_KZ}, expected ${N_KZ_EXPECTED} - see kz.het_drop.tsv and kz.mindrem.id"

log "=== [2/6] 1000 Genomes Phase 3 ==="
if [[ "${DS_1KG_FMT}" == "pfile" ]]; then
    [[ "${DS_1KG_VZS}" -eq 1 ]] && KG_SPEC="--pfile ${DS_1KG} vzs" || KG_SPEC="--pfile ${DS_1KG}"
else
    KG_SPEC="--bfile ${DS_1KG}"
fi
harmonise "${KG_SPEC}" "${PANELS}/kg" "${MIND_WGS}" wgs

log "=== [3/6] SGDP ==="
harmonise "--bfile ${DS_SGDP}" "${PANELS}/sgdp" "${MIND_WGS}" wgs

log "=== [4/6] HGDP 650k ==="
harmonise "--bfile ${DS_HGDP}" "${PANELS}/hgdp" "${MIND_ARRAY}" array

log "=== [5/6] AADR v62 1240K modern ==="
# AADR is filtered at the individual level, not just by group name. Four
# criteria, all recorded in results/meta/aadr_selection.tsv so the exclusions
# are auditable. Run aadr_inspect.sh first to see what these filters will remove.
require_file "${MF_AADR}"

if [[ -s "${MF_AADR_KEEP:-/nonexistent}" ]]; then
    # `|| true`: an inclusion list that is entirely comments is a configuration
    # choice, not a crash. grep exits 1 when it selects nothing.
    { grep -v -e '^#' -e '^[[:space:]]*$' "${MF_AADR_KEEP}" || true; } | sort -u \
        > "${WORK_DIR}/aadr_keep.txt"
    log "  inclusion list: $(nlines "${WORK_DIR}/aadr_keep.txt") groups"
else
    : > "${WORK_DIR}/aadr_keep.txt"
    log "  WARNING: no inclusion list; every group except Ignore_* and *.REF is eligible"
fi

awk -F'\t' -v OFS='\t' \
    -v idc="${MC_AADR_ID}" -v grc="${MC_AADR_NATION}" \
    -v needpass="${AADR_REQUIRE_PASS}" -v needdip="${AADR_REQUIRE_DIPLOID}" \
    -v dropout="${AADR_DROP_OUTLIER}" \
    -v minsnp="${AADR_MIN_SNPS}" \
    -v haskeep="$( [[ -s "${WORK_DIR}/aadr_keep.txt" ]] && echo 1 || echo 0 )" '
    # haskeep guards the NR==FNR idiom: with no inclusion list the first file is
    # empty, NR and FNR both start on the .janno, and every one of its rows gets
    # eaten by this rule. Same trap that emptied hgdp.mergeable.snplist.
    haskeep && NR==FNR { keep[$1]=1; nkeep++; next }
    FNR==1  { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g)
                                  gsub(/^[ \t]+|[ \t]+$/,"",g); h[g]=i }
              if (!(idc in h) || !(grc in h)) {
                  print "AADR: id or group column not found" > "/dev/stderr"; exit 3 }
              hasPL = ("Genotype_Ploidy"  in h)
              hasAS = ("AADR_Assessment"  in h)
              hasNS = ("Nr_SNPs"          in h)
              print "Poseidon_ID","Group_Name","ploidy","assessment","n_snps","verdict"
              next }
    { id=$(h[idc]); grp=$(h[grc])
      pl = hasPL ? $(h["Genotype_Ploidy"]) : "NA"
      as = hasAS ? $(h["AADR_Assessment"]) : "NA"
      ns = hasNS ? $(h["Nr_SNPs"])         : "NA"
      v = "keep"
      if      (grp ~ /^Ignore_/)                       v="drop:ignore_flag"
      else if (grp ~ /\.REF$/)                         v="drop:reference_placeholder"
      else if (dropout && grp ~ /_o(\.[A-Z]+)?$/)       v="drop:aadr_outlier_flag"
      else if (nkeep>0 && !(grp in keep))              v="drop:not_in_inclusion_list"
      else if (needdip && hasPL && pl!="diploid")      v="drop:not_diploid"
      else if (needpass && hasAS && as!="PASS")        v="drop:assessment_" as
      else if (hasNS && ns+0 < minsnp+0)               v="drop:low_coverage"
      print id, grp, pl, as, ns, v }' \
    "${WORK_DIR}/aadr_keep.txt" "${MF_AADR}" \
    > "${META_DIR}/aadr_selection.tsv"

log "  selection summary:"
awk -F'\t' 'NR>1{c[$6]++} END{for(v in c) printf "    %-34s %d\n", v, c[v]}' \
    "${META_DIR}/aadr_selection.tsv" | sort -k2,2nr

awk -F'\t' 'NR>1 && $6 ~ /^drop:/ && $6 !~ /inclusion_list|ignore_flag|reference_placeholder/ {
        printf "    EXCLUDED  %-28s %-24s %s\n", $1, $2, $6 }' \
    "${META_DIR}/aadr_selection.tsv"

awk -F'\t' 'NR>1 && $6=="keep" {print $1}' "${META_DIR}/aadr_selection.tsv" \
    > "${WORK_DIR}/aadr_keep_ids.txt"
mk_id_list "${DS_AADR}.fam" "${WORK_DIR}/aadr_keep.id" "${WORK_DIR}/aadr_keep_ids.txt"
log "  AADR retained: $(nlines "${WORK_DIR}/aadr_keep.id") of $(nlines "${DS_AADR}.fam") genotyped individuals"

# Ploidy is recorded rather than acted on here. Every retained individual whose
# Genotype_Ploidy is not "diploid" carries one randomly drawn allele per site
# written as a homozygote. Downstream that means: excluded from KING (no
# heterozygotes for the estimator to read), from the PC computation and from
# ADMIXTURE (zero within-individual heterozygosity is not something either model
# can absorb), from Fst (within-group heterozygosity is zero and the allele count
# reads 2 per observation) and from ROH (each genome is one run). They enter the
# PCA by projection only. Nothing in plink warns about any of this.
awk -F'\t' 'NR>1 && $6=="keep" && $3!="diploid" {print $1}' \
    "${META_DIR}/aadr_selection.tsv" | sort -u > "${MF_PSEUDOHAPLOID}"
NPH=$(nlines "${MF_PSEUDOHAPLOID}")
log "  pseudo-haploid among the retained: ${NPH}"
if [[ "${NPH}" -gt 0 ]]; then
    awk -F'\t' 'NR>1 && $6=="keep" && $3!="diploid" {n[$2]++}
        END{for(g in n) printf "    %-28s %d\n", g, n[g]}' \
        "${META_DIR}/aadr_selection.tsv" | sort -k2,2nr
fi
[[ -s "${WORK_DIR}/aadr_keep.id" ]] || die "AADR selection matched no genotyped individuals"

harmonise "--bfile ${DS_AADR} --keep ${WORK_DIR}/aadr_keep.id" \
          "${PANELS}/aadr" "${MIND_AADR}" capture "${GENO_AADR}"

# ---------------------------------------------------------------------------
log "=== [6/6] population labels ==="
# ---------------------------------------------------------------------------
# meta/pop_master.tsv is the one label file downstream reads:
#   IID  Nation  Population  Source  Platform
# Nation is the finest label the source panel provides. Population is the
# superpopulation, resolved through meta/superpop_map.tsv rather than copied
# from the source, because each panel names its regions differently.
# ---------------------------------------------------------------------------

extract_meta() {
    local FS_SPEC='\t'
    grep -qP '\t' "$1" || FS_SPEC='[ \t]+'
    awk -F"${FS_SPEC}" -v idc="$2" -v natc="$3" -v regc="$4" '
        NR==1 { for(i=1;i<=NF;i++){ g=$i; gsub(/^[ \t]+|[ \t]+$/,"",g); sub(/^#/,"",g); h[g]=i }
                if (!(idc in h))  { printf "column %s not found\n", idc  > "/dev/stderr"; exit 2 }
                if (!(natc in h)) { printf "column %s not found\n", natc > "/dev/stderr"; exit 2 }
                if (!(regc in h)) { printf "column %s not found\n", regc > "/dev/stderr"; exit 2 }
                next }
        $(h[idc]) != "" { printf "%s\t%s\t%s\n", $(h[idc]), $(h[natc]), $(h[regc]) }' "$1"
}

RAWLAB="${WORK_DIR}/labels_raw.tsv"
: > "${RAWLAB}"

if [[ -s "${MF_KZ}" ]]; then
    extract_meta "${MF_KZ}" \
                 "${MC_KZ_ID}" "${MC_KZ_NATION}" "${MC_KZ_NATION}" \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,"KAZ","KAZ","wgs"}' >> "${RAWLAB}"
else
    log "  ${MF_KZ} absent, labelling every Kazakh as Nation=Kazakh"
    fam_iids "${PANELS}/kz.fam" \
      | awk 'BEGIN{OFS="\t"}{print $1,"Kazakh","KAZ","KAZ","wgs"}' >> "${RAWLAB}"
fi

add_panel() {   # add_panel <meta file> <idc> <natc> <regc> <Source> <Platform>
    local F=$1 SRC=$5 PLT=$6
    require_file "${F}"
    extract_meta "${F}" "$2" "$3" "$4" \
      | awk -F'\t' -v src="${SRC}" -v plt="${PLT}" 'BEGIN{OFS="\t"}
            {print $1,$2,$3,src,plt}' >> "${WORK_DIR}/labels_ref.tsv"
}

: > "${WORK_DIR}/labels_ref.tsv"
add_panel "${MF_1KG}"  "${MC_1KG_ID}"  "${MC_1KG_NATION}"  "${MC_1KG_REGION}"  1KG  wgs
add_panel "${MF_SGDP}" "${MC_SGDP_ID}" "${MC_SGDP_NATION}" "${MC_SGDP_REGION}" SGDP wgs
add_panel "${MF_HGDP}" "${MC_HGDP_ID}" "${MC_HGDP_NATION}" "${MC_HGDP_REGION}" HGDP array

# AADR labels are restricted to the individuals actually selected above. Feeding
# the whole .janno in would put ~3800 rows into pop_master.tsv for people who are
# not genotyped in this study.
awk -F'\t' -v idc="${MC_AADR_ID}" '
    NR==FNR { if (FNR>1 && $6=="keep") sel[$1]=1; next }
    FNR==1  { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g)
                                  gsub(/^[ \t]+|[ \t]+$/,"",g); h[g]=i }
              print; next }
    ($(h[idc]) in sel)' \
    "${META_DIR}/aadr_selection.tsv" "${MF_AADR}" > "${WORK_DIR}/aadr_selected.janno"
add_panel "${WORK_DIR}/aadr_selected.janno" \
          "${MC_AADR_ID}" "${MC_AADR_NATION}" "${MC_AADR_REGION}" AADR capture

awk -F'\t' 'BEGIN{OFS="\t"}
    function base(x) { sub(/\.(DG|SG|AG|HO|REF)$/,"",x); sub(/_o[0-9]*$/,"",x); return x }
    NR==FNR { if ($0 ~ /^#/ || $1=="source_label" || $1=="") next; m[$1]=$2; next }
    { nat=$2; reg=$3
      sp = (nat in m)       ? m[nat] :
           (base(nat) in m) ? m[base(nat)] :
           (reg in m)       ? m[reg] :
           (base(reg) in m) ? m[base(reg)] : "UNMAPPED"
      print $1, nat, sp, $4, $5 }' \
    "${MF_SUPERPOP}" "${WORK_DIR}/labels_ref.tsv" >> "${RAWLAB}"

# ---------------------------------------------------------------------------
# One row per IID.
#
# SGDP's cteam_extended contains HGDP individuals under their HGDP identifiers
# (HGDP01253, HGDP00951, HGDP01344), so the same person is labelled twice. That
# is not a merge problem — groupA uses SGDP, groupB uses HGDP, and no group
# contains both — but pop_master.tsv is keyed on IID and every downstream join
# would silently keep whichever row awk read last.
#
# The superpopulation must agree across copies, because that is what drives the
# figures and the Fst grouping; if it does not, the two metadata files disagree
# about the person and that has to be resolved by hand rather than by sort order.
# ---------------------------------------------------------------------------
MASTER="${META_DIR}/pop_master.tsv"
sort -u "${RAWLAB}" > "${WORK_DIR}/labels_sorted.tsv"

# Nation and Platform are CANONICAL, not joined.
#
# 153 individuals appear in two source panels, and SGDP and HGDP spell the same
# population differently: Basque/French_Basque, Bedouin/BedouinB,
# Bergamo/North_Italian, Cambodian/Cambodians. The previous revision recorded
# both with a "|" join, which is fine as a record and fatal as a label: every
# downstream grouping - the Fst phenotype file, ROH by nation, ADMIXTURE labels
# - read "Basque|French_Basque" as a population of its own. Two of the 22
# French Basques ended up in a group of their own, below the n>=10 threshold,
# and out of the Fst table entirely.
#
# One label wins, by source priority. HGDP first because groupB - the panel that
# carries the Eurasian PCA, ADMIXTURE, ROH and the nation-level Fst - draws
# these individuals from HGDP; the SGDP copies of them are excluded in
# panel_exclude_iids.txt, and no SGDP population reaches n=10 anyway. The
# alternatives are not lost, they go to meta/nation_aliases.tsv.
awk -F'\t' -v OFS='\t' -v CONFLICT="${WORK_DIR}/label_conflicts.tsv" \
            -v COLLAPSED="${META_DIR}/duplicate_labels.tsv" \
            -v ALIASES="${META_DIR}/nation_aliases.tsv" '
    function prio(s) { return (s=="HGDP")?1 : (s=="1KG")?2 : (s=="SGDP")?3 : (s=="AADR")?4 : 5 }
    BEGIN { print "IID","Nation","Population","Source" > CONFLICT
            print "IID","Nations","Populations","Sources" > COLLAPSED
            print "IID","canonical_Nation","all_Nations","Sources" > ALIASES }
    { id=$1
      if (!(id in seen)) { seen[id]=1; all[id]=$2; src[id]=$4; n[id]=1
                           nat[id]=$2; plt[id]=$5; bp[id]=prio($4); pop[id]=$3 }
      else {
          n[id]++
          if (index("|"all[id]"|", "|"$2"|")==0) all[id]=all[id]"|"$2
          if (index("|"src[id]"|", "|"$4"|")==0) src[id]=src[id]"|"$4
          if (prio($4) < bp[id]) { nat[id]=$2; plt[id]=$5; bp[id]=prio($4) }
          if ($3 != pop[id]) { bad[id]=1
                               print id, $2, $3, $4 > CONFLICT } } }
    END { for (id in seen) {
              if (n[id] > 1) { print id, all[id], pop[id], src[id] > COLLAPSED
                               if (all[id] != nat[id])
                                   print id, nat[id], all[id], src[id] > ALIASES }
              print id, nat[id], pop[id], src[id], plt[id] }
          nb=0; for (id in bad) nb++
          if (nb > 0) exit 7 }' \
    "${WORK_DIR}/labels_sorted.tsv" | sort > "${WORK_DIR}/labels_unique.tsv" \
  || die "$(nbody "${WORK_DIR}/label_conflicts.tsv") IIDs get different superpopulations from different panels; see ${WORK_DIR}/label_conflicts.tsv"

{ printf 'IID\tNation\tPopulation\tSource\tPlatform\n'
  cat "${WORK_DIR}/labels_unique.tsv"; } > "${MASTER}"

NCOLL=$(nbody "${META_DIR}/duplicate_labels.tsv")
log "pop_master.tsv: $(nbody "${MASTER}") individuals"
log "  present in more than one panel, collapsed: ${NCOLL}"
log "  of those, spelt differently by their panels: $(nbody "${META_DIR}/nation_aliases.tsv") \
(canonical label chosen by source priority HGDP > 1KG > SGDP > AADR; see meta/nation_aliases.tsv)"
# No Nation may contain a "|" any more: that would mean a join escaped into a
# label and every downstream grouping is about to split a population in two.
awk -F'\t' 'NR>1 && $2 ~ /\|/ {print $1"\t"$2}' "${MASTER}" > "${WORK_DIR}/joined_nations.txt"
[[ ! -s "${WORK_DIR}/joined_nations.txt" ]] \
    || { head -5 "${WORK_DIR}/joined_nations.txt" | sed 's/^/        /'
         die "$(nlines "${WORK_DIR}/joined_nations.txt") Nation labels still contain a pipe"; }
[[ "${NCOLL}" -gt 0 ]] && \
    head -6 "${META_DIR}/duplicate_labels.tsv" | tail -n +2 | sed 's/^/        /'

DUPID=$(tail -n +2 "${MASTER}" | cut -f1 | sort | uniq -d)
[[ -z "${DUPID}" ]] || die "pop_master.tsv still has duplicate IIDs after collapsing"

# --- KEEP_NATIONS resolution ------------------------------------------------
# Verified here rather than discovered in the Fst table. Matching is normalised,
# because HGDP writes "Mbuti Pygmy" where SGDP writes "Mbuti"; the previous
# revision compared literally and dropped the HGDP Mbuti without a word.
{ printf 'token\tmatched_labels\tn_individuals\n'
  for T in ${KEEP_NATIONS}; do
      awk -F'\t' -v t="${T}" '
          function nz(s){ gsub(/[ _.-]/,"",s); return tolower(s) }
          # Nation may hold several pipe-joined labels for one individual, from
          # the duplicate-IID collapse above. Match any component, not the join.
          function hit(s,   n,p,i) { n=split(s,p,"|")
                                     for(i=1;i<=n;i++) if (nz(p[i])==nz(t)) return 1
                                     return 0 }
          NR>1 && hit($2) { n++; lab[$2]=1 }
          END { s=""; for (k in lab) s=(s==""?k:s"|"k)
                printf "%s\t%s\t%d\n", t, (s==""?"NONE":s), n+0 }' "${MASTER}"
  done
} > "${WORK_DIR}/keep_nations_raw.tsv"

# A token that matches nothing because canonicalisation resolved it to another
# spelling is NOT a thinner outgroup - the same individuals are still there
# under the canonical name. Annotate those as alias_of:<canonical> and warn only
# for tokens that are genuinely absent, so the warning keeps meaning something.
awk -F'\t' -v OFS='\t' '
    function nz(s){ gsub(/[ _.-]/,"",s); return tolower(s) }
    NR==FNR { if (FNR>1) { n=split($3,pp,"|")
                           for(i=1;i<=n;i++) alias[nz(pp[i])]=$2 }
              next }
    FNR==1  { print; next }
    { if ($2=="NONE") { k=nz($1); if (k in alias) $2="alias_of:" alias[k] }
      print }' \
    "${META_DIR}/nation_aliases.tsv" "${WORK_DIR}/keep_nations_raw.tsv" \
    > "${META_DIR}/keep_nations_resolution.tsv"
rm -f "${WORK_DIR}/keep_nations_raw.tsv"
log "KEEP_NATIONS resolution:"
sed 's/^/    /' "${META_DIR}/keep_nations_resolution.tsv"
awk -F'\t' 'NR>1 && $2=="NONE"{n++} END{exit !(n>0)}' "${META_DIR}/keep_nations_resolution.tsv" && \
    log "  WARNING: some KEEP_NATIONS tokens match nothing; the Fst outgroup will be thinner than intended"
# Every token failing means the outgroup is gone entirely and every Fst value
# loses its anchor. That is not a warning.
if [[ -n "${KEEP_NATIONS// /}" ]]; then
    awk -F'\t' 'NR>1 && $2!="NONE"{n++} END{exit (n>0)}' \
        "${META_DIR}/keep_nations_resolution.tsv" \
      && die "no KEEP_NATIONS token matched any label in pop_master.tsv; the Fst outgroup would be empty"
fi

# --- validation -------------------------------------------------------------
UNMAPPED=$(awk -F'\t' 'NR>1 && $3=="UNMAPPED"{print $2}' "${MASTER}" | sort -u)
if [[ -n "${UNMAPPED}" ]]; then
    printf '%s\n' "${UNMAPPED}" > "${WORK_DIR}/unmapped_labels.txt"
    log "  labels with no superpopulation: $(nlines "${WORK_DIR}/unmapped_labels.txt")"
fi

FAIL=0
for P in kz kg sgdp hgdp aadr; do
    fam_iids "${PANELS}/${P}.fam" | sort -u > "${PANELS}/${P}.iids"

    awk 'NR==FNR{lab[$1]=1; next} !($1 in lab)' \
        <(tail -n +2 "${MASTER}" | cut -f1) "${PANELS}/${P}.iids" \
        > "${PANELS}/${P}.unlabelled.id"
    if [[ -s "${PANELS}/${P}.unlabelled.id" ]]; then
        log "ERROR: ${P}: $(nlines "${PANELS}/${P}.unlabelled.id") genotyped samples have no label"
        head -10 "${PANELS}/${P}.unlabelled.id" | sed 's/^/        /'
        FAIL=1
    fi

    awk -F'\t' 'NR==FNR{g[$1]=1; next} FNR>1 && ($1 in g) && $3=="UNMAPPED"{print $2}' \
        "${PANELS}/${P}.iids" "${MASTER}" | sort -u > "${PANELS}/${P}.unmapped.txt"
    if [[ -s "${PANELS}/${P}.unmapped.txt" ]]; then
        log "ERROR: ${P}: $(nlines "${PANELS}/${P}.unmapped.txt") labels missing from superpop_map.tsv"
        head -10 "${PANELS}/${P}.unmapped.txt" | sed 's/^/        /'
        FAIL=1
    fi
done
[[ "${FAIL}" -eq 0 ]] || die "run make_map_template.sh, fill in the ??? entries, then rerun 01"

# --- IIDs shared between panels that will meet inside one analysis group -----
# Every downstream join (capping, Fst phenotypes, ROH annotation) is keyed on
# IID, and --keep resolves FID+IID pairs. Two individuals with the same IID
# inside one merged fileset therefore get selected and labelled together. 02
# refuses to proceed if it happens; this is where it is diagnosable.
for SPEC in "groupA:kz kg sgdp" "groupB:kz hgdp"; do
    GNAME="${SPEC%%:*}"
    # shellcheck disable=SC2086
    cat $(for P in ${SPEC#*:}; do printf '%s.iids ' "${PANELS}/${P}"; done) \
      | sort | uniq -d > "${WORK_DIR}/${GNAME}_shared_iids.txt"
    NSH=$(nlines "${WORK_DIR}/${GNAME}_shared_iids.txt")
    if [[ "${NSH}" -gt 0 ]]; then
        log "WARNING: ${GNAME}: ${NSH} IIDs appear in more than one of its panels"
        head -5 "${WORK_DIR}/${GNAME}_shared_iids.txt" | sed 's/^/        /'
        log "         list them in ${MF_PANEL_EXCLUDE} against the panel to drop, then rerun 01"
    fi
done

# --- per-panel QC summary for the supplement --------------------------------
{ printf 'panel\tplatform\tsamples\tvariants\tdup_ids\tambiguous\thet_dropped\thet_flagged\tmind_dropped\n'
  for P in kz kg sgdp hgdp aadr; do
    F="${PANELS}/${P}"
    case "${P}" in
      kz|kg|sgdp) PLT=wgs ;;
      hgdp)       PLT=array ;;
      aadr)       PLT=capture ;;
    esac
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${P}" "${PLT}" \
      "$(nlines "${F}.fam")" \
      "$(nlines "${F}.bim")" \
      "$( [[ -s "${F}_s1.rmdup.list" ]] && nlines "${F}_s1.rmdup.list" || echo 0 )" \
      "$(nlines "${F}.ambiguous.snplist")" \
      "$(nbody "${F}.het_drop.tsv")" \
      "$(nbody "${F}.het_flag.tsv")" \
      "$( [[ -f "${F}.mindrem.id" ]] && awk '!/^#/' "${F}.mindrem.id" | wc -l || echo 0 )"
  done
} > "${META_DIR}/panel_qc.tsv"
column -t -s $'\t' "${META_DIR}/panel_qc.tsv"

awk -F'\t' 'NR>1{c[$3]++} END{for(p in c) printf "  %-6s %d\n", p, c[p]}' "${MASTER}" | sort
log "=== 01 done -> 02_merge.sh ==="
