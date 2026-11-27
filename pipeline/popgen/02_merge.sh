#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# 02_merge.sh — intersect and merge each group. Produces merged_raw only.
#
# No MAF, no pruning, no relatedness filtering and no reference-panel capping at
# this stage. All four are frequency- or sample-dependent, and the sample set is
# not final until KING has run in 03. Filtering here would let individuals who
# are about to be removed decide which variants survive.
#
# Group A  KZ + 1000G + SGDP           all WGS
# Group B  KZ + HGDP650                WGS + array
#
# VARIANT IDs. Inside the merge the ID is CHR:POS and nothing else. An earlier
# revision rewrote it as CHR:POS:<alleles sorted alphabetically>, intending to
# make G:A and A:G match. It did, but it also made A:G and T:C — the same site on
# opposite strands — into two different IDs. plink1.9 matches variants by ID, so
# those sites never appeared in a .missnp list, the --flip branch never fired for
# them, and every genuinely strand-flipped site was carried through the merge as
# two half-empty variants and then deleted by --geno in 03. Silently.
#
# With CHR:POS, plink1.9 sees one variant and compares alleles itself. 01 has
# already collapsed each panel to one record per position, so CHR:POS is unique
# within every panel.
#
# MERGE CONFLICTS ARE CLASSIFIED, NOT ASSUMED. The previous revision logged
# every .missnp entry as a "strand conflict" and flipped. On groupA that
# produced 12885 conflicts before the flip and 12885 after it — the flip
# resolved nothing, because none of them were strand conflicts. They are
# positions where the two panels kept different biallelic representations of the
# same multi-allelic site, which is a direct consequence of the one-record-per-
# position collapse in 01 and which no flip can fix. Both counts went into the
# log under a label that said otherwise. Each conflict is now classified as
# strand_flip / allele_mismatch / same_alleles, the counts go to
# merge_conflicts.tsv, and the log says which is which. groupB, where flipping
# genuinely does work (234 -> 72), reads differently from groupA, as it should.
#
# Superpopulations listed in DROP_SUPERPOPS are removed here, before the merge,
# rather than at plotting time. KEEP_NATIONS holds back a small outgroup so the
# Fst scale keeps an anchor; matching is normalised, because the panels spell
# the same population differently.
# ============================================================================

need "${PLINK1}"; need "${PLINK2}"
require_file "${META_DIR}/pop_master.tsv"

panel_keep() {
    local P=$1 OUT=$2
    awk -F'\t' -v drop="${DROP_SUPERPOPS}" -v keep="${KEEP_NATIONS}" '
        function nz(s){ gsub(/[ _.-]/,"",s); return tolower(s) }
        # Nation may be several pipe-joined labels for an individual present in
        # more than one panel; match any component. Comparing the whole join
        # dropped the SGDP/HGDP Mbuti overlap and left groupB with Yoruba alone.
        function kept(s,   n,p,i) { n=split(s,p,"|")
                                    for(i=1;i<=n;i++) if (nz(p[i]) in K) return 1
                                    return 0 }
        BEGIN { n=split(drop,a," "); for(i=1;i<=n;i++) D[a[i]]=1
                n=split(keep,b," "); for(i=1;i<=n;i++) K[nz(b[i])]=1 }
        NR>1 { if (($3 in D) && !kept($2)) next
               print $1 }' "${META_DIR}/pop_master.tsv" > "${OUT}.ids"
    mk_id_list "${P}.fam" "${OUT}" "${OUT}.ids"
    rm -f "${OUT}.ids"
    local NB NA; NB=$(nlines "${P}.fam"); NA=$(nlines "${OUT}")
    log "  $(basename "${P}"): ${NA} of ${NB} retained after dropping ${DROP_SUPERPOPS}"
}

HIGHLD="${WORK_DIR}/highld_hg19.txt"
cat > "${HIGHLD}" << 'EOF'
1 48000000 52000000 r1
2 86000000 100500000 r2
2 134500000 138000000 r3
2 183000000 190000000 r4
3 47500000 50000000 r5
3 83500000 87000000 r6
3 89000000 97500000 r7
5 44500000 50500000 r8
5 98000000 100500000 r9
5 129000000 132000000 r10
5 135500000 138500000 r11
6 25500000 33500000 r12
6 57000000 64000000 r13
6 140000000 142500000 r14
7 55000000 66000000 r15
8 8000000 12000000 r16
8 43000000 50000000 r17
10 37000000 43000000 r18
11 46000000 57000000 r19
11 87500000 90500000 r20
12 33000000 40000000 r21
12 109500000 112000000 r22
17 40000000 43000000 r23
EOF

# ---------------------------------------------------------------------------
# safe_merge A B OUT — plink1.9 merge, flip B on the .missnp list, then drop
# whatever is still unresolved. Aborts if more than 5% of B stays unresolved:
# at that level the panels are not on the same build.
# ---------------------------------------------------------------------------
safe_merge() {
    local A=$1 B=$2 OUT=$3
    log "  merge $(basename "${A}") + $(basename "${B}")"

    rm -f "${OUT}-merge.missnp"

    ${PLINK1} --bfile "${A}" --bmerge "${B}" --allow-no-sex --make-bed \
              --threads "${THREADS}" --out "${OUT}" || true

    if [[ ! -f "${OUT}-merge.missnp" ]]; then
        log "    no merge conflicts"
        log "    -> $(nlines "${OUT}.fam") samples, $(nlines "${OUT}.bim") variants"
        return 0
    fi

    local NM; NM=$(nlines "${OUT}-merge.missnp")
    cp "${OUT}-merge.missnp" "${OUT}.conflicts_initial.missnp"
    log "    merge conflicts: ${NM}"
    classify_missnp "${OUT}.conflicts_initial.missnp" "${A}.bim" "${B}.bim" \
                    "${OUT}.conflicts_initial.tsv" | sort | sed 's/^/      /'

    ${PLINK1} --bfile "${B}" --flip "${OUT}-merge.missnp" --allow-no-sex \
              --make-bed --threads "${THREADS}" --out "${B}_flip"
    rm -f "${OUT}-merge.missnp"
    ${PLINK1} --bfile "${A}" --bmerge "${B}_flip" --allow-no-sex --make-bed \
              --threads "${THREADS}" --out "${OUT}" || true

    if [[ ! -f "${OUT}-merge.missnp" ]]; then
        log "    all ${NM} resolved by --flip"
        log "    -> $(nlines "${OUT}.fam") samples, $(nlines "${OUT}.bim") variants"
        rm -f "${B}_flip."{bed,bim,fam}
        return 0
    fi

    local NU NB
    NU=$(nlines "${OUT}-merge.missnp"); NB=$(nlines "${B}.bim")
    cp "${OUT}-merge.missnp" "${OUT}.conflicts_unresolved.missnp"
    log "    resolved by --flip: $(( NM - NU )); still unresolved: ${NU} \
($(awk -v a="${NU}" -v b="${NB}" 'BEGIN{printf "%.3f", 100*a/b}')% of $(basename "${B}"))"
    classify_missnp "${OUT}.conflicts_unresolved.missnp" "${A}.bim" "${B}.bim" \
                    "${OUT}.conflicts_unresolved.tsv" | sort | sed 's/^/      /'
    log "    unresolved positions are dropped from BOTH panels; detail in \
$(basename "${OUT}").conflicts_unresolved.tsv"

    awk -v a="${NU}" -v b="${NB}" 'BEGIN{exit !(100*a/b > 5)}' && \
        die "unresolved merge conflicts exceed 5% - panels are probably on different builds"

    ${PLINK1} --bfile "${A}"      --exclude "${OUT}-merge.missnp" --allow-no-sex \
              --make-bed --threads "${THREADS}" --out "${A}_cl"
    ${PLINK1} --bfile "${B}_flip" --exclude "${OUT}-merge.missnp" --allow-no-sex \
              --make-bed --threads "${THREADS}" --out "${B}_cl"
    ${PLINK1} --bfile "${A}_cl" --bmerge "${B}_cl" --allow-no-sex --make-bed \
              --threads "${THREADS}" --out "${OUT}"
    rm -f "${B}_flip."{bed,bim,fam} "${A}_cl."{bed,bim,fam} "${B}_cl."{bed,bim,fam}
    log "    -> $(nlines "${OUT}.fam") samples, $(nlines "${OUT}.bim") variants"
}

# ---------------------------------------------------------------------------
# merge_group GROUP panel [panel ...]
# ---------------------------------------------------------------------------
merge_group() {
    local GROUP=$1; shift
    local NAMES=("$@")
    local NP=${#NAMES[@]}
    local G="${WORK_DIR}/${GROUP}"; mkdir -p "${G}"
    log "=== ${GROUP}  (${NP} panels: ${NAMES[*]}) ==="

    local PREFIX=()
    local N
    for N in "${NAMES[@]}"; do
        PREFIX+=("$(panel_path "${N}")")
        require_file "$(panel_path "${N}").mergeable.snplist"
    done

    # --- 1. intersect on CHR:POS -------------------------------------------
    log "[1/4] intersection"
    rm -f "${G}"/pos_*.txt
    local i=1
    for P in "${PREFIX[@]}"; do
        cut -d: -f1,2 "${P}.mergeable.snplist" | sort -u > "$(printf '%s/pos_%02d.txt' "${G}" "${i}")"
        log "  ${NAMES[$((i-1))]}: $(nlines "$(printf '%s/pos_%02d.txt' "${G}" "${i}")") unique positions"
        i=$((i+1))
    done
    sort -m "${G}"/pos_*.txt | uniq -c | awk -v n="${NP}" '$1==n{print $2}' > "${G}/common.pos"
    local NC; NC=$(nlines "${G}/common.pos")
    log "  common to all ${NP}: ${NC}"
    [[ "${NC}" -gt 10000 ]] || die "${GROUP}: only ${NC} shared positions"

    # --- 2. recover full IDs ------------------------------------------------
    log "[2/4] recovering variant IDs"
    i=1
    for P in "${PREFIX[@]}"; do
        awk -F: 'NR==FNR{k[$1":"$2]=1; next} ($1":"$2) in k' \
            "${G}/common.pos" "${P}.mergeable.snplist" > "${G}/ids${i}.snplist"
        i=$((i+1))
    done

    # --- 3. subset + CHR:POS IDs -------------------------------------------
    log "[3/4] subsetting, ID -> CHR:POS"
    i=1
    for P in "${PREFIX[@]}"; do
        panel_keep "${P}" "${G}/keep${i}.id"
        [[ -s "${G}/keep${i}.id" ]] || die "${GROUP}: ${NAMES[$((i-1))]} has nobody left after DROP_SUPERPOPS"
        ${PLINK2} --bfile "${P}" --extract "${G}/ids${i}.snplist" \
                  --keep "${G}/keep${i}.id" \
                  --make-bed --threads "${THREADS}" --out "${G}/sub${i}"
        awk 'BEGIN{OFS="\t"}{$2=$1":"$4; print}' "${G}/sub${i}.bim" > "${G}/sub${i}n.bim"
        local ND
        ND=$(cut -f2 "${G}/sub${i}n.bim" | sort | uniq -d | wc -l)
        [[ "${ND}" -eq 0 ]] || die "${GROUP}: ${NAMES[$((i-1))]} has ${ND} duplicate CHR:POS after subsetting"
        ln -sf "${G}/sub${i}.bed" "${G}/sub${i}n.bed"
        cp     "${G}/sub${i}.fam" "${G}/sub${i}n.fam"
        i=$((i+1))
    done

    # --- 4. merge -----------------------------------------------------------
    log "[4/4] merge"
    local ACC="${G}/sub1n"
    for ((i=2; i<=NP; i++)); do
        local OUT
        if [[ "${i}" -eq "${NP}" ]]; then OUT="${G}/merged_raw"; else OUT="${G}/m1_${i}"; fi
        safe_merge "${ACC}" "${G}/sub${i}n" "${OUT}"
        ACC="${OUT}"
    done

    local NM; NM=$(nlines "${G}/merged_raw.bim")
    log "  intersection ${NC} -> merged ${NM} ($(awk -v a="${NM}" -v b="${NC}" \
        'BEGIN{printf "%.2f", 100*a/b}')% retained through conflict resolution)"

    # Duplicate IIDs are fatal, not a note. Every downstream join — reference
    # capping, the Fst phenotype file, ROH platform annotation — is keyed on
    # IID, and --keep resolves FID+IID pairs, so two individuals sharing an IID
    # inside one fileset are selected and labelled as one person. 01 reports
    # which panels they come from; MF_PANEL_EXCLUDE is how you choose the copy
    # to drop.
    fam_iids "${G}/merged_raw.fam" | sort | uniq -d > "${G}/dup_iids.txt"
    if [[ -s "${G}/dup_iids.txt" ]]; then
        head -10 "${G}/dup_iids.txt" | sed 's/^/        /'
        die "${GROUP}: $(nlines "${G}/dup_iids.txt") duplicate IIDs in merged_raw. \
Add the copy to drop to ${MF_PANEL_EXCLUDE} as <panel><TAB><IID> and rerun 01."
    fi

    ${PLINK2} --bfile "${G}/merged_raw" --missing --threads "${THREADS}" \
              --out "${G}/merged_raw"

    { printf 'panel\tpositions\n'
      i=1; for N in "${NAMES[@]}"; do
        printf '%s\t%s\n' "${N}" "$(nlines "$(printf '%s/pos_%02d.txt' "${G}" "${i}")")"
        i=$((i+1))
      done
      printf 'intersection\t%s\n' "${NC}"
      printf 'merged\t%s\n'       "${NM}"
    } > "${G}/merge_flow.tsv"

    # One conflict table per group, so the supplement can state exactly what was
    # lost and why rather than "strand conflicts were resolved".
    { printf 'stage\tclass\tn\n'
      for F in "${G}"/*.conflicts_initial.tsv "${G}"/*.conflicts_unresolved.tsv; do
          [[ -s "${F}" ]] || continue
          local ST; ST=$(basename "${F}" .tsv)
          awk -F'\t' -v st="${ST}" 'NR>1{n[$4]++} END{for(k in n) printf "%s\t%s\t%d\n", st, k, n[k]}' "${F}"
      done
    } > "${G}/merge_conflicts.tsv"

    log "${GROUP}: $(nlines "${G}/merged_raw.fam") samples, ${NM} variants"
    rm -f "${G}"/pos_*.txt "${G}"/sub[0-9]*.{bed,bim,fam} "${G}"/ids*.snplist
}

merge_group "groupA_kz_1000g_sgdp" "${GROUPA_PANELS[@]}"
merge_group "groupB_kz_hgdp"       "${GROUPB_PANELS[@]}"

log "=== marker intersections ==="
for G in "${ANALYSIS_GROUPS[@]}"; do
    [[ -s "${WORK_DIR}/${G}/merge_flow.tsv" ]] && \
        { printf '  %s\n' "${G}"; sed 's/^/    /' "${WORK_DIR}/${G}/merge_flow.tsv"; }
    [[ -s "${WORK_DIR}/${G}/merge_conflicts.tsv" ]] && \
        sed 's/^/    /' "${WORK_DIR}/${G}/merge_conflicts.tsv"
done
log "groupB carries the Eurasian PCA, ADMIXTURE and ROH, so its merged count is"
log "the one that has to be comfortable. The AADR panel is not merged here; 08"
log "builds its own two-panel merge for the projection."

log "=== 02 done -> 03_relatedness.sh ==="
