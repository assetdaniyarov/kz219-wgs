#!/bin/bash
set -uo pipefail
source "$(dirname "$0")/config.sh"

# ============================================================================
# aadr_inspect.sh — read-only report on the AADR .janno. Run before editing
# aadr_keep_groups.txt; writes nothing, changes nothing.
#
# Answers four questions that the .janno can answer and the .fam cannot:
#
#   1. Which groups exist, how large, and what does the inclusion list catch.
#   2. Are any retained individuals PSEUDO-HAPLOID. This is the one that has to
#      be checked before anything else. A pseudo-haploid individual carries one
#      randomly sampled allele per site coded as a homozygote. Mixed into a
#      diploid panel it produces zero heterozygosity, an F of ~1, an Fst against
#      every other group inflated by the missing within-group variance, and a
#      genome-length run of homozygosity. Nothing in plink warns about it.
#   3. AADR_Assessment. AADR flags its own individuals PASS / QUESTIONABLE /
#      QUESTIONABLE_CRITICAL. Only PASS belongs in a reference panel.
#   4. Whether retained AADR individuals are the same people as HGDP, SGDP or
#      1000G individuals already in another panel. AADR .DG entries are
#      re-genotyped copies: Adygei.DG is HGDP01401/HGDP01402, Xibo.DG is
#      HGDP01244-HGDP01251, Yakut.DG is HGDP00951/HGDP00956, and the ACB/YRI/
#      GBR/TSI/FIN/JPT/KHV/STU/ITU/BEB/GIH .DG groups are 1000G individuals.
#      Their Poseidon_ID carries the .DG suffix, so an exact-IID duplicate check
#      never sees them; only KING does, after the merge. The Alternative_IDs and
#      Collection_ID columns carry the original identifier and are checked here.
# ============================================================================

require_file "${MF_AADR}"

hdr() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }

# column index by header name
col() { awk -F'\t' -v want="$1" 'NR==1{for(i=1;i<=NF;i++){g=$i; sub(/^#/,"",g); if(g==want){print i; exit}}}' "${MF_AADR}"; }

C_ID=$(col "${MC_AADR_ID}")
C_GR=$(col "${MC_AADR_NATION}")
C_PL=$(col Genotype_Ploidy)
C_AS=$(col AADR_Assessment)
C_NS=$(col Nr_SNPs)
C_ALT=$(col Alternative_IDs)
C_COL=$(col Collection_ID)

hdr "columns"
for SPEC in "Poseidon_ID:${C_ID}" "Group_Name:${C_GR}" "Genotype_Ploidy:${C_PL}" \
            "AADR_Assessment:${C_AS}" "Nr_SNPs:${C_NS}" \
            "Alternative_IDs:${C_ALT}" "Collection_ID:${C_COL}"; do
    N="${SPEC%%:*}"; I="${SPEC#*:}"
    [[ -n "${I}" ]] && printf '  %-18s col %s\n' "${N}" "${I}" \
                    || printf '  \033[33m%-18s NOT PRESENT\033[0m\n' "${N}"
done

TOT=$(( $(wc -l < "${MF_AADR}") - 1 ))
NGRP=$(awk -F'\t' -v c="${C_GR}" 'NR>1{g[$c]=1} END{print length(g)}' "${MF_AADR}")
printf '\n  %s individuals in %s groups\n' "${TOT}" "${NGRP}"

hdr "genotype ploidy, whole panel"
[[ -n "${C_PL}" ]] && awk -F'\t' -v c="${C_PL}" 'NR>1{n[$c]++} END{for(k in n) printf "  %-14s %d\n", k, n[k]}' \
    "${MF_AADR}" | sort -k2,2nr || echo "  column absent"

hdr "AADR_Assessment, whole panel"
[[ -n "${C_AS}" ]] && awk -F'\t' -v c="${C_AS}" 'NR>1{n[$c]++} END{for(k in n) printf "  %-24s %d\n", k, n[k]}' \
    "${MF_AADR}" | sort -k2,2nr || echo "  column absent"

# ---------------------------------------------------------------------------
if [[ ! -s "${MF_AADR_KEEP:-/nonexistent}" ]]; then
    hdr "inclusion list"
    echo "  ${MF_AADR_KEEP:-unset} absent - nothing further to check"
    exit 0
fi

KEEP=$(mktemp); trap 'rm -f "${KEEP}"' EXIT
grep -v -e '^#' -e '^[[:space:]]*$' "${MF_AADR_KEEP}" | sort -u > "${KEEP}"

hdr "inclusion list: $(wc -l < "${KEEP}") groups"
awk -F'\t' -v c="${C_GR}" -v pl="${C_PL}" -v as="${C_AS}" -v ns="${C_NS}" '
    NR==FNR { want[$1]=1; next }
    FNR==1  { next }
    ($c in want) { n[$c]++
                   if (pl && $pl!="diploid") hap[$c]++
                   if (as && $as!="PASS")    bad[$c]++
                   if (ns) { s[$c]+=$ns; if (mn[$c]==""||$ns<mn[$c]) mn[$c]=$ns } }
    END { printf "  %-32s %5s %6s %6s %11s %11s\n", "group","n","nondip","notPASS","mean_SNPs","min_SNPs"
          for (k in n)
            printf "  %-32s %5d %6d %6d %11d %11d\n", k, n[k],
                   (k in hap?hap[k]:0), (k in bad?bad[k]:0),
                   (ns? s[k]/n[k]:0), (ns? mn[k]:0) }' \
    "${KEEP}" "${MF_AADR}" | { read -r h; echo "$h"; sort -k2,2nr; }

hdr "inclusion-list names absent from the .janno"
awk -F'\t' -v c="${C_GR}" '
    NR==FNR { want[$1]=1; next }
    FNR>1   { have[$c]=1 }
    END     { n=0; for (g in want) if (!(g in have)) { print "  "g; n++ }
              if (n==0) print "  none - every name exists" }' "${KEEP}" "${MF_AADR}"

hdr "retained individuals that are also HGDP / SGDP / 1000G people"
# HGDPnnnnn, LP6005..., NA/HG identifiers appearing in Alternative_IDs or
# Collection_ID mean the same person is already in another panel.
if [[ -n "${C_ALT}${C_COL}" ]]; then
    awk -F'\t' -v c="${C_GR}" -v id="${C_ID}" -v a="${C_ALT}" -v b="${C_COL}" '
        NR==FNR { want[$1]=1; next }
        FNR==1  { next }
        ($c in want) {
            v = (a? $a : "") " " (b? $b : "")
            src=""
            if (v ~ /HGDP[0-9]/)          src="HGDP"
            else if (v ~ /(NA|HG)[0-9]{5}/) src="1000G"
            else if (v ~ /^S_|LP600/)     src="SGDP"
            if (src!="") printf "  %-10s %-28s %-24s %s\n", src, $id, $c, v }
        END {}' "${KEEP}" "${MF_AADR}" | sort -k1,1 -k3,3 | head -40
    echo "  (empty above = the inclusion list contains no re-genotyped copies)"
else
    echo "  Alternative_IDs and Collection_ID both absent, cannot check"
fi

hdr "what the .SG set contributes that HGDP/SGDP/1000G do not"
awk -F'\t' -v c="${C_GR}" '
    NR==FNR { want[$1]=1; next }
    FNR>1 && ($c in want) && $c ~ /\.SG$/ { n[$c]++ }
    END { for (k in n) printf "  %-32s %d\n", k, n[k] }' "${KEEP}" "${MF_AADR}" | sort -k2,2nr

printf '\n'
