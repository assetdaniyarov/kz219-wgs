#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# check_inputs.sh — run after setup_meta.sh, before run_all.sh.
#
# Confirms binaries, genotypes and metadata are present, reports panel sizes,
# and lists the population labels that superpop_map.tsv does not cover.
#
# The coverage check resolves labels exactly the way 01 does: nation label
# first, region label as a fallback, UNMAPPED only if both miss. A nation label
# with no entry of its own is normal and is not reported — that is what the
# region fallback is for.
# ============================================================================

FAIL=0
ok()   { printf '  \033[32mOK\033[0m    %-56s %s\n' "$1" "${2:-}"; }
bad()  { printf '  \033[31mMISS\033[0m  %-56s %s\n' "$1" "${2:-}"; FAIL=1; }
warn() { printf '  \033[33mWARN\033[0m  %-56s %s\n' "$1" "${2:-}"; }
chk()  { [[ -s "$1" ]] && ok "$(basename "$1")" "$(du -hL "$1" 2>/dev/null | cut -f1)" || bad "$1"; }

echo "=== binaries ==="
for B in "${PLINK1}" "${PLINK2}" "${ADMIXTURE}"; do
    [[ -x "${B}" ]] && ok "${B}" || bad "${B}"
done


echo
echo "=== python and helper scripts ==="
# The pipeline no longer needs a notebook, so these are load-bearing.
if command -v "${PYTHON}" >/dev/null; then
    ok "${PYTHON}" "$(${PYTHON} --version 2>&1)"
    for M in numpy scipy; do
        V=$(${PYTHON} -c "import ${M}; print(${M}.__version__)" 2>/dev/null || true)
        [[ -n "${V}" ]] && ok "python module ${M}" "${V}" \
            || bad "python module ${M}" "pip install ${M} into the ${CONDA_ENV} env"
    done
else
    bad "${PYTHON}" "set PYTHON in config.sh"
fi
for F in popgen_tables.py admixture_align.py; do
    [[ -s "${SCRIPT_DIR}/${F}" ]] && ok "${F}" || bad "${SCRIPT_DIR}/${F}"
done

echo
echo "=== configuration sanity ==="
awk -v c="${KING_CUTOFF}" 'BEGIN{exit !(c>0.125)}' \
    && warn "KING_CUTOFF=${KING_CUTOFF}" "above 0.125: true 2nd-degree pairs will survive" \
    || ok "KING_CUTOFF" "${KING_CUTOFF} (removes 2nd degree and closer)"
awk -v r="${PRUNE_R2}" 'BEGIN{exit !(r<0.15)}' \
    && warn "PRUNE_R2=${PRUNE_R2}" "aggressive; ADMIXTURE and KING may end up short of markers" \
    || ok "PRUNE_R2" "${PRUNE_R2}"
[[ "${KZ_GQ_MIN}" -gt 0 ]] && ok "KZ_GQ_MIN" "${KZ_GQ_MIN}" \
    || warn "KZ_GQ_MIN=0" "no genotype-level filter on the cohort VCF"
[[ "${ADMIX_PROJECT}" -eq 1 ]] && ok "ADMIX_PROJECT" "cohort capped for training, projected with -P" \
    || warn "ADMIX_PROJECT=0" "joint fit only; the cohort is ~25% of the panel"

echo
echo "=== directory layout ==="
for D in "${RAW_DIR}" "${RES_DIR}" "${WORK_DIR}" "${LOG_DIR}" "${META_DIR}"; do
    [[ -d "${D}" ]] && ok "${D#${PROJECT}/}" || bad "${D}"
done

echo
echo "=== genotypes ==="
ok "PROJECT resolved to" "${PROJECT}"
chk "${DS_KZ_VCF}"
if [[ -s "${DS_KZ_VCF}" ]] && command -v bcftools >/dev/null; then
    N=$(bcftools query -l "${DS_KZ_VCF}" | wc -l)
    [[ "${N}" -eq "${N_KZ_EXPECTED}" ]] && ok "KZ sample count" "${N}" \
        || warn "KZ sample count" "${N}, expected ${N_KZ_EXPECTED}"
fi
# KZ_GQ_MIN is applied at import time in 01. Two things are checked here and
# one is deliberately NOT: the fraction of calls below the threshold is NOT
# estimated from a `head` of the stream. That stream starts at the beginning of
# chr1, i.e. the subtelomeric pile-up where every callset looks terrible, and an
# earlier version of this check reported "75.7% of calls below GQ 20" from the
# first 750 variants. That is an artefact of where the reading stopped, not a
# property of the callset, and acting on it either way would be wrong.
# gq_profile.sh samples genome-spread windows and replays the whole 01 filter
# chain at several thresholds instead.
if [[ "${KZ_GQ_MIN}" -gt 0 && -s "${DS_KZ_VCF}" ]] && command -v bcftools >/dev/null; then
    if bcftools view -h "${DS_KZ_VCF}" 2>/dev/null | grep -q '^##FORMAT=<ID=GQ,'; then
        ok "FORMAT/GQ present" "KZ_GQ_MIN=${KZ_GQ_MIN} will mask calls below it"
        if [[ -s "${DS_KZ_VCF}.tbi" || -s "${DS_KZ_VCF}.csi" ]]; then
            if [[ -s "${RES_DIR}/provenance/gq_profile.txt" ]]; then
                ok "  GQ profile present" "$(basename "${RES_DIR}")/provenance/gq_profile.txt"
            else
                warn "  GQ profile not yet run" "bash gq_profile.sh  <- do this before run_all.sh"
            fi
        else
            warn "  cohort VCF is not indexed" \
                 "bcftools index -t '${DS_KZ_VCF}'  (gq_profile.sh needs it to seek)"
        fi
    else
        warn "FORMAT/GQ absent from the cohort VCF" \
             "--vcf-min-gq would mask nothing; set KZ_GQ_MIN=0 or drop the Methods claim"
    fi
fi
chk "${DS_1KG}.pgen"; chk "${DS_1KG}.psam"
if [[ "${DS_1KG_VZS}" -eq 1 ]]; then
    [[ -s "${DS_1KG}.pvar.zst" ]] && ok "all_phase3.pvar.zst" \
        || bad "${DS_1KG}.pvar.zst" "or set DS_1KG_VZS=0"
else
    [[ -s "${DS_1KG}.pvar" ]] && ok "all_phase3.pvar" || bad "${DS_1KG}.pvar"
fi
for SPEC in "SGDP:${DS_SGDP}" "HGDP:${DS_HGDP}" "AADR:${DS_AADR}"; do
    NAME="${SPEC%%:*}"; P="${SPEC#*:}"; M=0
    for E in bed bim fam; do [[ -s "${P}.${E}" ]] || M=1; done
    [[ "${M}" -eq 0 ]] \
        && ok "${NAME}" "$(wc -l < "${P}.fam") samples, $(wc -l < "${P}.bim") variants" \
        || bad "${NAME}  ${P}.{bed,bim,fam}"
done

echo
echo "=== metadata ==="
for F in "${MF_1KG}" "${MF_SGDP}" "${MF_HGDP}" "${MF_AADR}" "${MF_SUPERPOP}"; do chk "${F}"; done
[[ -s "${MF_KZ}" ]] && ok "$(basename "${MF_KZ}")" \
    || warn "kz_219_metadata.tsv" "optional; every Kazakh gets Nation=Kazakh"
[[ -s "${MF_AADR_KEEP:-}" ]] \
    && ok "aadr_keep_groups.txt" "$(grep -vc -e '^#' -e '^$' "${MF_AADR_KEEP}") groups retained" \
    || warn "aadr_keep_groups.txt" "absent: all AADR groups enter the merge"

# --- per-panel metadata vs genotype coverage --------------------------------
# A genotyped sample with no metadata row aborts 01, so it is found here.
cover() {   # cover <label> <meta file> <id col> <fam>
    local NAME=$1 F=$2 IDC=$3 FAM=$4
    [[ -s "${F}" && -s "${FAM}" ]] || return
    local MISS
    MISS=$(awk -F'\t' -v idc="${IDC}" '
        NR==FNR { if (FNR==1) { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g);
                                gsub(/^[ \t]+|[ \t]+$/,"",g); h[g]=i }
                                if (!(idc in h)) exit 3; next }
                  lab[$(h[idc])]=1; next }
        !($2 in lab) { print $2 }' "${F}" "${FAM}")
    local N; N=$(printf '%s\n' "${MISS}" | grep -c . || true)
    if [[ "${N}" -eq 0 ]]; then
        ok "${NAME}: every genotyped sample has a metadata row"
    else
        warn "${NAME}: ${N} genotyped samples with no metadata row" "01 will abort"
        printf '%s\n' "${MISS}" | head -10 | sed 's/^/          /'
    fi
}
echo
echo "=== metadata covers the genotypes ==="
cover SGDP "${MF_SGDP}" "${MC_SGDP_ID}" "${DS_SGDP}.fam"
cover HGDP "${MF_HGDP}" "${MC_HGDP_ID}" "${DS_HGDP}.fam"
cover AADR "${MF_AADR}" "${MC_AADR_ID}" "${DS_AADR}.fam"


echo
echo "=== KEEP_NATIONS resolution ==="
# These tokens hold an African outgroup back from DROP_SUPERPOPS so the Fst
# scale keeps an anchor. Matching is normalised (case, spaces, underscores),
# because HGDP writes "Mbuti Pygmy" where SGDP writes "Mbuti". A token that
# matches nothing used to fail silently and leave groupB with Yoruba alone.
for T in ${KEEP_NATIONS}; do
    HITS=$(for SPEC in "${MF_1KG}:${MC_1KG_NATION}" "${MF_SGDP}:${MC_SGDP_NATION}" \
                       "${MF_HGDP}:${MC_HGDP_NATION}"; do
        F="${SPEC%:*}"; COL="${SPEC##*:}"
        [[ -s "${F}" ]] || continue
        FS_SPEC=$'\t'; grep -qP '\t' "${F}" || FS_SPEC='[ \t]+'
        awk -F"${FS_SPEC}" -v want="${T}" -v col="${COL}" '
            function nz(s){ gsub(/[ _.-]/,"",s); return tolower(s) }
            NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g);
                                        gsub(/^[ \t]+|[ \t]+$/,"",g); h[g]=i }
                    if (!(col in h)) exit; c=h[col]; next }
            nz($c)==nz(want) { n++ } END { if (n>0) printf "%s(%d) ", FILENAME, n }' "${F}"
    done)
    if [[ -n "${HITS}" ]]; then
        ok "KEEP_NATIONS token '${T}'" "${HITS}"
    else
        warn "KEEP_NATIONS token '${T}'" "matches no label in any panel"
        # Matching is exact-after-normalisation on purpose: substring matching
        # would quietly pull in neighbours ("Han" catching "Han-NChina"). So
        # when a token misses, show the labels it ALMOST matches, with the exact
        # spelling to paste into KEEP_NATIONS.
        for SPEC in "${MF_1KG}:${MC_1KG_NATION}" "${MF_SGDP}:${MC_SGDP_NATION}" \
                    "${MF_HGDP}:${MC_HGDP_NATION}"; do
            F="${SPEC%:*}"; COL="${SPEC##*:}"
            [[ -s "${F}" ]] || continue
            FS_SPEC=$'\t'; grep -qP '\t' "${F}" || FS_SPEC='[ \t]+'
            awk -F"${FS_SPEC}" -v want="${T}" -v col="${COL}" -v src="$(basename "${F}")" '
                function nz(s){ gsub(/[ _.-]/,"",s); return tolower(s) }
                function pfx(a,b,  i,m) { m=(length(a)<length(b))?length(a):length(b)
                                          for(i=1;i<=m;i++) if (substr(a,i,1)!=substr(b,i,1)) return i-1
                                          return m }
                NR==1 { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g)
                                            gsub(/^[ \t]+|[ \t]+$/,"",g); h[g]=i }
                        if (!(col in h)) exit; c=h[col]; w=nz(want)
                        lim=(length(w)<5)?length(w):5; next }
                { l=nz($c)
                  if (l!=w && pfx(l,w)>=lim) n[$c]++ }
                END { for (k in n) printf "          did you mean %-24s (%s, n=%d)\n", "\047" k "\047", src, n[k] }' "${F}"
        done
    fi
done

echo
echo "=== superpopulation coverage ==="
[[ -s "${MF_SUPERPOP}" ]] || { bad "${MF_SUPERPOP}"; exit 1; }
ok "entries in superpop_map.tsv" "$(grep -vc -e '^#' -e '^[[:space:]]*$' "${MF_SUPERPOP}")"

# resolve() applies the same rule as 01: nation, then region, then UNMAPPED.
# AADR Group_Name suffixes (.DG/.SG/.AG/.HO/.REF) and outlier markers (_o, _o1)
# are stripped before the lookup, same as in 01.
resolve() {   # resolve <label> <meta file> <nation col> <region col>
    local NAME=$1 F=$2 NATC=$3 REGC=$4
    [[ -s "${F}" ]] || return
    local OUT
    OUT=$(awk -F'\t' -v natc="${NATC}" -v regc="${REGC}" '
        function base(x) { sub(/\.(DG|SG|AG|HO|REF)$/,"",x); sub(/_o[0-9]*$/,"",x); return x }
        NR==FNR { if ($0 ~ /^#/ || $1=="source_label" || $1=="") next; m[$1]=1; next }
        FNR==1  { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g);
                                      gsub(/^[ \t]+|[ \t]+$/,"",g); h[g]=i }
                  if (!(natc in h) || !(regc in h)) exit 3; next }
        { nat=$(h[natc]); reg=$(h[regc])
          if (nat in m || base(nat) in m || reg in m || base(reg) in m) next
          u[nat]++ }
        END { for (k in u) printf "%s\t%d\n", k, u[k] }' \
        "${MF_SUPERPOP}" "${F}" | sort -k2,2nr)
    local N; N=$(printf '%s\n' "${OUT}" | grep -c . || true)
    if [[ "${N}" -eq 0 ]]; then
        ok "${NAME}: every label resolves"
    else
        local NIND; NIND=$(printf '%s\n' "${OUT}" | awk '{s+=$2} END{print s+0}')
        warn "${NAME}: ${N} labels unresolved" "${NIND} individuals; 01 aborts"
        printf '%s\n' "${OUT}" | head -15 | awk -F'\t' '{printf "          %-40s %s\n", $1, $2}'
        [[ "${N}" -gt 15 ]] && printf '          ... and %d more\n' "$(( N - 15 ))"
        FAIL=1
    fi
}
resolve 1KG  "${MF_1KG}"  "${MC_1KG_NATION}"  "${MC_1KG_REGION}"
resolve SGDP "${MF_SGDP}" "${MC_SGDP_NATION}" "${MC_SGDP_REGION}"
resolve HGDP "${MF_HGDP}" "${MC_HGDP_NATION}" "${MC_HGDP_REGION}"

# AADR is checked against the inclusion list only. Labels outside the list
# never enter the merge, so they do not need a superpop_map entry and
# reporting them buries the ones that do.
if [[ -s "${MF_AADR_KEEP:-/nonexistent}" ]]; then
    grep -v -e '^#' -e '^[[:space:]]*$' "${MF_AADR_KEEP}" | sort -u > /tmp/.aadr_keep.$$
    awk -F'\t' -v grc="${MC_AADR_NATION}" '
        NR==FNR { keep[$1]=1; next }
        FNR==1  { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }
                  print; next }
        ($(h[grc]) in keep)' /tmp/.aadr_keep.$$ "${MF_AADR}" > /tmp/.aadr_sub.$$

    NKEEP=$(( $(wc -l < /tmp/.aadr_sub.$$) - 1 ))
    NLIST=$(wc -l < /tmp/.aadr_keep.$$)
    if [[ "${NKEEP}" -eq 0 ]]; then
        bad "AADR inclusion list matches nothing" "${NLIST} group names, 0 individuals"
    else
        ok "AADR inclusion list" "${NLIST} groups -> ${NKEEP} of $(( $(wc -l < "${MF_AADR}") - 1 )) individuals"
        # group names in the list that do not exist in the .janno
        awk -F'\t' -v grc="${MC_AADR_NATION}" '
            FNR==1 && NR!=FNR { for(i=1;i<=NF;i++){ g=$i; sub(/^#/,"",g); h[g]=i }; next }
            NR==FNR { want[$1]=1; next }
            { have[$(h[grc])]=1 }
            END { for (g in want) if (!(g in have)) print g }' \
            /tmp/.aadr_keep.$$ "${MF_AADR}" | sort > /tmp/.aadr_missing.$$
        NMISS=$(wc -l < /tmp/.aadr_missing.$$)
        if [[ "${NMISS}" -gt 0 ]]; then
            warn "AADR list names not present in the .janno" "${NMISS} of ${NLIST}"
            head -15 /tmp/.aadr_missing.$$ | sed 's/^/          /'
        else
            ok "every AADR list name exists in the .janno"
        fi
        resolve AADR /tmp/.aadr_sub.$$ "${MC_AADR_NATION}" "${MC_AADR_REGION}"
    fi
    rm -f /tmp/.aadr_keep.$$ /tmp/.aadr_sub.$$ /tmp/.aadr_missing.$$
else
    warn "no aadr_keep_groups.txt" "all 343 AADR labels need superpop_map entries"
    resolve AADR "${MF_AADR}" "${MC_AADR_NATION}" "${MC_AADR_REGION}"
fi


echo
echo "=== within-cohort analysis (09) ==="
if [[ -s "${MF_TRIBES}" ]]; then
    ok "$(basename "${MF_TRIBES}")" "$(( $(wc -l < "${MF_TRIBES}") - 1 )) rows"
    for COLNAME in IID Tribe Zhuz; do
        head -1 "${MF_TRIBES}" | tr '\t' '\n' | sed 's/^#//' \
          | grep -qix "${COLNAME//_/}" && ok "  column ${COLNAME}" \
          || warn "  column ${COLNAME}" "not found by name; 09 falls back to a fixed position"
    done
else
    warn "$(basename "${MF_TRIBES}")" "absent; 09_kz_tribes.sh cannot run"
fi

echo
echo "=== disk ==="
df -hP "${PROJECT}" 2>/dev/null | tail -1 | awk '{printf "  %s free on %s\n", $4, $6}'
echo "  budget: ~400 GB for results/work"

echo
if [[ "${FAIL}" -eq 0 ]]; then
    echo "Ready. Run ./run_all.sh"
else
    echo "Fix the items above first."
    echo "For unresolved labels: bash make_map_template.sh >> ${MF_SUPERPOP}"
    exit 1
fi
