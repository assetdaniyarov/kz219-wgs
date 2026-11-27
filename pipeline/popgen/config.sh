#!/bin/bash
# ============================================================================
# config.sh  (rev. 2026-09-02)
#
# Single source of truth. Every script sources this and defines nothing of its
# own — no script re-declares PLINK1/PLINK2/ADMIXTURE locally, or the stages
# stop agreeing on which binary was used.
#
# CHANGES IN THIS REVISION, and why each one exists:
#
#   KING_CUTOFF 0.13 -> 0.0884. 0.13 sits above the expected kinship of a
#   genuine 2nd-degree pair (0.125), so roughly half of all true 2nd-degree
#   relatives survive it by sampling noise. The previous comment claimed "3rd
#   degree and closer", which is 0.0442 and was never the value in use.
#
#   PRUNE_R2 0.1 -> 0.2 and PRUNE_WINDOW 1000kb -> 200kb. The old setting took
#   groupB from 606k merged markers to 52,421, which is thin for ADMIXTURE at
#   K=10 and below the ~100k KING wants.
#
#   Cohort capping. Reference nations were capped at 25 while the 219-sample
#   cohort was not, so the cohort was a quarter of the ADMIXTURE panel and the
#   largest group in it by a factor of ~9. ADMIXTURE now trains on a balanced
#   set and projects the full cohort onto it (-P), which is the standard fix.
#
#   KZ_GQ_MIN. Genotype-level QC on the cohort VCF, made possible by splitting
#   --geno and --mind into separate passes. Reviewers ask what genotype filter
#   was applied and "none" is not an answer.
# ============================================================================

export LC_ALL=C          # deterministic sort collation across all stages

# ---------------------------------------------------------------- paths ----
# Resolution order: an exported PROJECT wins, otherwise data/popgen inside this
# repository (the location notebook 02 reads through KZ_PROJECT). The analysis
# itself was run on the laboratory NAS; only the root differs.
#   PROJECT=/some/other/root ./run_all.sh
export PROJECT="${PROJECT:-}"
if [[ -z "${PROJECT}" ]]; then
    _cand="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/data/popgen"
    [[ -d "${_cand}/raw_datasets" ]] && PROJECT="${_cand}"
    unset _cand
fi
if [[ -z "${PROJECT}" || ! -d "${PROJECT}/raw_datasets" ]]; then
    printf 'config.sh: no project root with a raw_datasets/ directory found.\n' >&2
    printf '  set PROJECT explicitly (default: <repo>/data/popgen). Do NOT let\n' >&2
    printf '  mkdir -p create an empty tree at a path the data is not on.\n' >&2
    # exit, not return: `return 1` from a sourced file only stops the caller
    # under `set -e`, and check_inputs.sh does not use it.
    exit 1
fi

export RAW_DIR="${PROJECT}/raw_datasets"
export RES_DIR="${PROJECT}/results"
export WORK_DIR="${RES_DIR}/work"
export LOG_DIR="${RES_DIR}/logs"
export META_DIR="${RES_DIR}/meta"

# Directory holding the scripts themselves, so the shell stages can call the
# Python helpers without depending on the caller's working directory.
export SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export PYTHON="${PYTHON:-python3}"

# work/ is a few hundred GB of intermediates. When packaging for Zenodo:
#   tar --exclude='work' -czf kz219_results.tar.gz -C "${PROJECT}" results

# -------------------------------------------------------------- binaries ----
# Versions used for the article: PLINK v1.9.0-b.8, PLINK v2.00a5.12, ADMIXTURE v1.3.0.
# Override with the full path of each binary if they are not on PATH.
export PLINK1="${PLINK1:-plink}"          # 1.9: merge, ROH, --ibc, --genome
export PLINK2="${PLINK2:-plink2}"         # everything else
export ADMIXTURE="${ADMIXTURE:-admixture}"
export CONDA_ENV="${CONDA_ENV:-plink2}"

# Concurrency, derived from the host rather than pinned. On the 64-core
# analysis machine this gives THREADS=32, ADMIX_JOBS=16, ADMIX_THREADS=4.
#
# For ADMIXTURE, many concurrent processes with few threads each beats one wide
# process: thread scaling flattens out past ~8, while the runs are completely
# independent. 16x4 is roughly 4.7x the throughput of the 2x8 the first run
# used, on the same cores.
#
# Override any of them from the environment: ADMIX_JOBS=8 ./run_all.sh
_NPROC=$(nproc 2>/dev/null || echo 16)
_MEMMB=$(awk '/MemTotal/{printf "%d", $2/1024}' /proc/meminfo 2>/dev/null || echo 65536)
export THREADS="${THREADS:-$(( _NPROC > 32 ? 32 : _NPROC ))}"
export MEM_MB="${MEM_MB:-$(( _MEMMB/2 > 131072 ? 131072 : _MEMMB/2 ))}"

# ------------------------------------------------------------ genotypes ----
# Sibling of PROJECT on both hosts, so it follows the mount point rather than
# being pinned to one of them.
export DS_KZ_VCF="${DS_KZ_VCF:-$(dirname "${PROJECT}")/kz_219_hg19_multisample_vcf/kz_219_hg19_dragen.vcf.gz}"
export DS_1KG="${RAW_DIR}/1000g_phase3_20160505_37_2504/all_phase3"
export DS_1KG_FMT="pfile"
export DS_1KG_VZS=1
export DS_SGDP="${RAW_DIR}/sgdp_reich_lab_2019/cteam_extended.v4.maf0.1perc"
export DS_HGDP="${RAW_DIR}/hgdp_650k_snps/hgdp650.qc.hg19"
export DS_AADR="${RAW_DIR}/aadr_v62_1240k_modern/AADR_v62_1240K_Modern"

export N_KZ_EXPECTED=219

# --------------------------------------------------------- VCF import ----
# --var-filter keeps PASS records only. --vcf-half-call missing prevents an
# abort if any half-call is present.
export KZ_VCF_FILTERS=(--var-filter --vcf-half-call missing)

# Genotype-level QC on the cohort VCF. The previous revision left this off with
# the note that --vcf-min-gq 20 makes --mind delete all 219 samples. That is
# true only when --geno and --mind run in the same plink2 invocation: --mind is
# evaluated on the raw 20.9M sites, most of which are rare and are exactly the
# ones GQ masking empties. 01_prepare.sh now runs three separate passes —
# lenient --geno, then --mind, then the final --geno — so masking is safe.
#
# Set KZ_GQ_MIN=0 to disable and reproduce the old behaviour.
# --vcf-min-dp is deliberately absent: this DRAGEN VCF carries no FORMAT/DP.
export KZ_GQ_MIN=20
export GENO_PRE=0.10           # lenient variant pass, before --mind
# Abort rather than silently gut a panel if the sample filter removes more than
# this fraction. A --mind that deletes a third of a panel is a broken threshold,
# not a result.
export MIND_MAX_FRAC_REMOVED=0.05

# The cohort gets an ABSOLUTE limit instead of a fraction, and it is zero.
#
# MIND_MAX_FRAC_REMOVED is a disaster brake for the reference panels: 5% of 2504
# 1000G samples is 125 people and losing a few of them changes nothing. 5% of
# the cohort is 10 people out of 219, and n=219 is the headline of the paper.
# Worse, the two knobs are coupled: measured at KZ_GQ_MIN=30, --mind 0.02 would
# have removed 86 samples and tripped the 5% brake, while --mind 0.05 removes 8
# and slips under it. Loosening MIND_KZ therefore also loosened what the brake
# catches, which is not what loosening MIND_KZ was for.
#
# At the configured KZ_GQ_MIN the expected number is zero (gq_profile.sh), so
# any removal at all means an assumption broke and a human should look. If a
# cohort sample genuinely has to go, list it in panel_exclude_iids.txt - then
# the decision is recorded in the metadata instead of happening inside a
# threshold. Raise this only with a reason you would put in Methods.
export MIND_MAX_REMOVED_KZ=0

# ------------------------------------------------------- metadata files ----
export MF_1KG="${DS_1KG}.psam"
export MF_SGDP="${RAW_DIR}/sgdp_reich_lab_2019/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt"
export MF_HGDP="${RAW_DIR}/hgdp_650k_snps/hgdp650_id_pop_coords.txt"
export MF_AADR="${RAW_DIR}/aadr_v62_1240k_modern/AADR_v62_1240K_Modern.janno"
export MF_KZ="${META_DIR}/kz_219_metadata.tsv"
export MF_SUPERPOP="${META_DIR}/superpop_map.tsv"
export MF_AADR_KEEP="${META_DIR}/aadr_keep_groups.txt"
export MF_TRIBES="${META_DIR}/kz_219_tribes.tsv"

# ------------------------------------------------- metadata column names ----
export MC_1KG_ID=IID           MC_1KG_NATION=Population     MC_1KG_REGION=SuperPop
# SGDP: the .fam carries sample IDs (HGDP01253, SAH31, IHW9118), not Illumina
# run IDs. Illumina_ID matches 0 of 345; the alias column matches 343. The
# parentheses have to be quoted or bash reads them as a subshell.
export MC_SGDP_ID='Sample_ID(Aliases)'
export MC_SGDP_NATION=Population_ID MC_SGDP_REGION=Region
export MC_HGDP_ID=id           MC_HGDP_NATION=population    MC_HGDP_REGION=region
export MC_AADR_ID=Poseidon_ID  MC_AADR_NATION=Group_Name    MC_AADR_REGION=Group_Name
export MC_KZ_ID=IID            MC_KZ_NATION=Nation

# ============================================================================
# GROUPS AND THEIR ROLES
#
# groupA  KZ + 1000G + SGDP        diploid WGS, millions of shared sites
#         -> Fst at superpopulation level (main text), global PCA
#            (supplementary). At nation level almost every SGDP population is
#            n<10 and is excluded, so groupA has no Central Asian neighbours;
#            those live in groupB.
#
# groupB  KZ + HGDP650             diploid, one array platform across the whole
#                                  reference side, n = 9-51 per population
#         -> Eurasian PCA (main), ADMIXTURE, cross-population ROH, nation-level
#            Fst against the cohort's actual neighbours (Uygur, Hazara, Yakut,
#            Mongola, the Pakistani set). 06_fst.sh cross-checks the two panels
#            at superpopulation level so the platform difference is measured
#            rather than asserted.
#
# AADR    the 46 retained AADR .SG individuals are PSEUDO-HAPLOID. Verified in
#         the .janno: all 21 groups haploid, n=1-4, 1240K coverage 227,950 to
#         1,231,875. One randomly drawn allele per site written as a homozygote
#         — zero heterozygosity by construction, doubled per-genotype sampling
#         variance, and an allele count reading 2 where the observation is 1.
#         -> In NO analysis group. 08_project_aadr.sh scores them onto the
#            groupB axes as one supplementary panel. Never Fst, never
#            ADMIXTURE, never ROH, never contributing to an axis.
#
# The AADR .SG set is kept despite all of that because no diploid source has
# Uzbek, Turkmen, Karakalpak, Tatar, Bashkir, Kalmyk, Buryat, Teleut or an
# independent Kazakh sample.
# ============================================================================
export ANALYSIS_GROUPS=(groupA_kz_1000g_sgdp groupB_kz_hgdp)
export GROUPA_PANELS=(kz kg sgdp)
export GROUPB_PANELS=(kz hgdp)

# Both panels, so ADMIXTURE gets the same platform replication that Fst and
# PCA already have. Read the groupA result for what it can carry: it is all
# WGS on both sides, so it shows the continental placement is not an artefact
# of mixing a WGS cohort with an array reference - but 1000G has no Central
# Asians and SGDP has n=1-3, so groupA CANNOT reproduce the Uygur/Hazara
# resolution that makes groupB the main figure.
# Cost: 03 must rebuild merged_admix for groupA, then 05 runs the full sweep
# again for it (~6 h on 64 cores). Drop back to (groupB_kz_hgdp) to skip.
export ADMIX_GROUPS=(groupB_kz_hgdp groupA_kz_1000g_sgdp)
export FST_GROUPS=(groupA_kz_1000g_sgdp groupB_kz_hgdp)
export ROH_GROUP="groupB_kz_hgdp"

# The array is NOT called GROUPS. GROUPS is a bash builtin holding the numeric
# group IDs of the current user, assignments to it are silently discarded, and
# `for G in "${GROUPS[@]}"` then iterates over "0" under root instead of over
# the group names — with no error from bash and a confusing missing-file error
# several stages later. The guard below fails loudly if any name here is ever
# shadowed the same way.
for _v in ANALYSIS_GROUPS GROUPA_PANELS GROUPB_PANELS ADMIX_GROUPS FST_GROUPS KZ_VCF_FILTERS; do
    declare -n _chk="${_v}"
    [[ "${#_chk[@]}" -gt 0 && "${_chk[0]}" != "0" ]] \
        || { printf 'config.sh: %s did not survive assignment (shadowed builtin?)\n' "${_v}" >&2; exit 1; }
    unset -n _chk
done
unset _v

# ------------------------------------------------- superpopulation scope ----
# Removed from every panel before the merge. Dropping AFR and AMR removes 1008
# of the 2504 1000G individuals and 121 of the 345 SGDP individuals.
export DROP_SUPERPOPS="AFR AMR"

# Held back from the dropped superpopulations: they anchor the Fst scale and
# outgroup-based statistics need an African outgroup by construction. They
# appear in merged_fst and nowhere else.
#
# Matching is NORMALISED (case-folded, spaces/underscores/hyphens/dots removed)
# because the panels spell the same population differently: HGDP writes
# "Mbuti Pygmy", SGDP writes "Mbuti", 1000G writes "YRI". The previous revision
# split this string on spaces and compared literally, so "Mbuti Pygmy" matched
# neither "Mbuti" nor "Pygmy" and the HGDP Mbuti were silently dropped from
# groupB — leaving that panel's African outgroup as Yoruba alone.
# check_inputs.sh now verifies that every token here matches something.
# HGDP spells it Mbuti_HG (n=15); SGDP writes Mbuti (n=4); 1000G writes YRI.
# MbutiPygmy was a guess and matched nothing.
#
# Matching is normalised AND pipe-aware. 01 collapses individuals present in
# more than one panel into a single pop_master row and joins their labels with
# "|", so the SGDP/HGDP Mbuti overlap produces Nation="Mbuti|Mbuti_HG". Comparing
# the whole string against "Mbuti" fails, which is how the African outgroup for
# groupB silently became Yoruba alone on the 2026-09-02 19:34 run.
export KEEP_NATIONS="YRI Yoruba Mbuti Mbuti_HG"

# Additionally removed from the Eurasian PCA and from ADMIXTURE. AFR and AMR are
# already gone by then, so this is Oceania alone. With the African outgroup and
# Oceania left in, K=2 splits Africa from everything else, K=3 adds Oceania and
# K=4 splits East from West Eurasia — three components spent before anything
# Central Asian appears.
export EURASIA_EXCLUDE="OCE"

# --------------------------------------------------------- sample-level QC ----
export MIND_WGS=0.02
export MIND_ARRAY=0.05

# The cohort gets its own --mind, and it is looser than MIND_WGS on purpose.
#
# MIND_WGS=0.02 was calibrated against UNMASKED data. Genotype masking
# deliberately manufactures missingness, so once KZ_GQ_MIN is on, --mind is
# measuring a different quantity and has to be recalibrated with it. Applying
# the pre-masking threshold to post-masking data is not conservatism, it is a
# unit error, and it silently deletes cohort members - the one thing this study
# cannot afford, since n=219 is the headline.
#
# Measured by gq_profile.sh on 12,958 PASS variants across eight autosomal
# windows: at KZ_GQ_MIN=20 the worst sample sits at 0.0183 per-sample
# missingness, i.e. 92% of the way to 0.02. That is not headroom, it is a
# coin flip on a favourable sample of the genome. 0.05 is ~2.7x the observed
# maximum and still catches a genuinely failed sample; MIND_MAX_FRAC_REMOVED
# catches the disaster case regardless.
#
# Re-run gq_profile.sh and re-read this if KZ_GQ_MIN changes.
export MIND_KZ=0.05

# The reference panels arrive as pfile/bed with no GQ field, so no masking is
# possible there and MIND_WGS keeps its original meaning for them. The
# asymmetry is a consequence of the data, not a double standard.
export HET_SD=3

# AADR .SG individuals carry 227,950 to 1,231,875 of the 1,135,618 panel sites,
# so the sparsest are ~80% missing by construction and MIND_ARRAY removed 16 of
# 46 — preferentially the small-marker groups the panel exists for, among them
# a Kazakh.SG. Missingness here is coverage, not quality, and 08 normalises by
# non-missing genotypes and plots F_MISS beside every point.
export MIND_AADR=0.90
# --geno has to move with it: keeping all 46 means a site would have to be
# called in 44 of them, and the panel collapsed from 1,135,618 variants to 198.
export GENO_AADR=0.90

# Direction of the heterozygosity filter.
#   excess  drop only samples whose F is HET_SD below the panel mean, i.e.
#           excess heterozygosity = contamination or sample mixture.
#   both    drop either tail.
# 'excess' is the correct default. A heterozygosity deficit in a Kazakh WGS
# panel is consanguinity, which is what 07 sets out to measure; removing those
# individuals biases FROH downward and deletes the signal. Deficit outliers go
# to <panel>.het_flag.tsv and are retained.
export HET_MODE=excess
export HET_MAF=0.05

# ------------------------------------------------------ AADR individual QC ----
# REQUIRE_DIPLOID stays 0. Setting it to 1 removes all 46 retained individuals,
# because the entire AADR Central Asian .SG set is pseudo-haploid. They are
# admitted deliberately and confined to the projection role above; the ploidy is
# recorded in meta/pseudohaploid_iids.txt and enforced downstream, not here.
# AADR marks an individual whose profile departs from its own group with an
# "_o" suffix on the group name. That is a curator's verdict about that sample,
# not a spelling of the population, so such groups are dropped rather than
# merged into the parent: folding Tajik_o into Tajik would put a flagged
# outlier into a group of three and move its centroid. Set to 0 to keep them.
export AADR_DROP_OUTLIER=1

export AADR_REQUIRE_PASS=1
export AADR_REQUIRE_DIPLOID=0
# Pseudo-haploid AADR group = its diploid counterpart in the same projection.
# 08 uses these to measure projection shrinkage like-for-like instead of
# comparing Central Asian AADR samples against a diploid set that includes
# Papuans and Africans on the rim of the PC space.
export AADR_MATCHED_PAIRS="Kazakh.SG=Kazakh Uyghur.SG=Uygur Hazara.SG=Hazara"

export AADR_MIN_SNPS=200000    # retained set spans 227,950 to 1,231,875 called
                               # 1240K sites: admits all of them and still
                               # catches anything worse in a future release.

# -------------------------------------------------------- variant-level QC ----
export GENO=0.02
export MAF=0.05                # PCA/ADMIXTURE only. Never applied to the Fst set.
# Ambiguous A/T and C/G SNPs are excluded AT MERGE TIME, not when the panel is
# built. Allele identity cannot resolve strand between two panels, which is why
# they have to go before any merge. Inside a single panel there is no strand
# question: this VCF is called against hg19 forward and every genotype in it
# uses one convention.
#
# Applying the filter at panel-build time cost the cohort 3,565,929 variants -
# 18.9% of the panel - and handed that loss to every Kazakh-only analysis
# (07 ROH on the full WGS set, 09 Fst/PCA/KING) where it buys nothing. 01 now
# writes two variant lists per panel: <panel>.snplist for standalone use and
# <panel>.mergeable.snplist with the ambiguous sites removed, which is what
# 02 and 08 read.
export DROP_AMBIGUOUS=1

# ------------------------------------------------------------- LD pruning ----
# PCA / ADMIXTURE marker set.
export PRUNE_WINDOW=200kb
export PRUNE_STEP=1
export PRUNE_R2=0.2
export ADMIX_MIN_SNPS=100000   # warn below this
export ADMIX_MAX_SNPS=250000   # note above this: runtime grows, resolution does not

# KING is pruned separately and more leniently. Kinship precision is set by the
# number of independent markers, and the PCA threshold left groupB with 55,710.
export KING_PRUNE_WINDOW=200kb
export KING_PRUNE_STEP=1
export KING_PRUNE_R2=0.2
export KING_MIN_MARKERS=100000

# ------------------------------------------------------------- relatedness ----
# 0.0884 = remove 2nd degree and closer. The bands are KING's own:
#   >0.354 duplicate/MZ, 0.177-0.354 1st, 0.0884-0.177 2nd, 0.0442-0.0884 3rd.
# 0.13 (the previous value) lies above the 0.125 expected for a true 2nd-degree
# pair, so it retains roughly half of them; it was chosen to keep all 219
# Kazakhs and costs more in review than the one sample is worth. At 0.0884 one
# Kazakh of a 2nd-degree pair (kinship ~0.118) is removed: 218 of 219.
export KING_CUTOFF=0.0884
export KING_REPORT=0.0442
export KING_DUP=0.354
# Additional cutoffs written as keep-lists only, so the sensitivity statement in
# Methods is backed by files rather than by an assertion.
export KING_CUTOFF_SENSITIVITY="0.0442 0.13"
# KING runs on the diploid individuals only. Kinship estimated from
# pseudo-haploid genotypes is not interpretable: the estimator reads
# heterozygote sharing and there are no heterozygotes.

# ------------------------------------------- reference panel size balancing ----
# PCA axes and ADMIXTURE components are driven by the populations contributing
# the most individuals (McVean 2009). Selection is deterministic: lowest
# genotype missingness first, IID as the tie-break.
export PCA_MAX_PER_NATION=25

# The cohort is NOT capped in the primary PCA/ADMIXTURE set — the cohort is the
# subject of the paper and capping it throws away the data the paper is about.
# Instead both analyses are repeated with the cohort capped, and the two are
# compared:
#   PCA        merged_pca_capped, reference individuals scored onto both sets
#              of axes, correlation of PC1..PC4 reported (04_pca.sh).
#   ADMIXTURE  trained on the capped set, full cohort projected with -P
#              (05_admixture.sh). This is the primary ADMIXTURE result.
export COHORT_CAP=25           # cap applied to KAZ in the sensitivity/training sets
export ADMIX_PROJECT=1         # 1 = train capped + project; 0 = joint run only
export ADMIX_JOINT_SENSITIVITY=1   # also run the uncapped joint fit at the
                                   # CV-selected K and correlate the two Q

# --------------------------------------------------------------- ADMIXTURE ----
export K_MIN=2
export K_MAX=10
export N_REPS=10
export CV_FOLDS=10
export SEED_BASE=20260101
# See the THREADS block above. ADMIX_JOBS defaults to one process per
# ADMIX_THREADS cores.
export ADMIX_THREADS="${ADMIX_THREADS:-4}"
export ADMIX_JOBS="${ADMIX_JOBS:-$(( _NPROC / ADMIX_THREADS < 1 ? 1 : _NPROC / ADMIX_THREADS ))}"
export ADMIX_SIM_THRESHOLD=0.97   # pong's default for "same mode"

# --------------------------------------------------------------------- PCA ----
export N_PC=20
export PCA_ALL_SAMPLES=0       # supplementary "relatedness filter changes
                               # nothing" panel; off by default

# --------------------------------------------------------------------- Fst ----
# Appended to the per-level results directory. Empty by default, so the main
# run keeps writing to fst/nation and fst/population; set it for a second run
# at another threshold and nothing is overwritten.
export FST_LEVEL_SUFFIX="${FST_LEVEL_SUFFIX:-}"

export FST_BLOCKSIZE=1000
export FST_REPORT_VARIANTS=0
export FST_MIN_N_POP="${FST_MIN_N_POP:-10}"
# Overridable from the environment, so a second run at a different threshold
# does not need config.sh edited and reverted:
#     FST_MIN_N_NATION=5 FST_LEVEL_SUFFIX=_minn5 bash 06_fst.sh
export FST_MIN_N_NATION="${FST_MIN_N_NATION:-10}"     # Hudson is unbiased at n=2, but the jackknife
                               # interval there is wider than the whole range of
                               # the table. Groups below this are excluded from
                               # the run and listed, not plotted.
export FST_PLATFORM_CHECK=1    # cross-panel concordance at superpopulation
                               # level: the evidence that groupB's array/WGS
                               # mixture is not driving the KZ-vs-other values

# --------------------------------------------------------------------- ROH ----
export ROH_MAF=0.05
export ROH_WINDOW_MISSING=5
export ROH_WINDOW_THRESHOLD=0.05
export ROH_GAP=1000            # kb
export AUTOSOME_BP=2881033286  # GRCh37 autosomes, FROH denominator

# name : window-snp : window-het : homozyg-snp : homozyg-kb : density (kb/SNP)
# Two sets, because marker density differs by more than an order of magnitude
# between the full Kazakh WGS panel and the merged array set. The two FROH
# values are on different scales and never share an axis.
export ROH_WGS="wgs:50:1:50:500:50"
export ROH_ARRAY="array:50:1:50:1500:50"
# Length classes for the FROH breakdown, in kb. Short runs reflect ancient
# background relatedness, long runs recent consanguinity, and reporting only a
# total FROH hides which of the two the cohort actually has.
export ROH_CLASSES="500-1000 1000-2000 2000-4000 4000-8000 8000-"

# ------------------------------------------------------- within-cohort (09) ----
export TRIBE_MIN_N=10
export KZ_K_MIN=1              # K=1 on purpose: within one population the
                               # question is whether any K beats K=1
export KZ_K_MAX=6
export KZ_N_REPS=10
export PERM_N=10000            # label permutations for the PC-vs-tribe test

# ------------------------------------------------------------------ helpers ----
log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*"; }
die() { printf '[%s] ERROR: %s\n' "$(date +'%F %T')" "$*" >&2; exit 1; }

need() { [[ -x "$1" ]] || die "binary not found or not executable: $1"; }
require_file() { [[ -s "$1" ]] || die "required file missing or empty: $1"; }

fam_iids() { awk '{print $2}' "$1"; }
nlines()   { awk 'END{print NR+0}' "$1"; }
nbody()    { awk 'END{print NR-1}' "$1"; }   # line count minus a header row

# in_list <needle> <haystack...>
in_list() { local n=$1; shift; local x; for x in "$@"; do [[ "${x}" == "${n}" ]] && return 0; done; return 1; }

# mk_id_list <fam> <out> [infile]
# Reads IIDs (one per line) from infile, or from stdin when infile is omitted,
# and writes the matching FID<TAB>IID pairs from <fam>.
#
# NEVER pass the output path as infile: the redirection truncates it before the
# reader runs. 01_prepare.sh lost every heterozygosity outlier in the 1000G
# panel to exactly that, silently, and the panel kept 4 contaminated samples.
mk_id_list() {
    local FAM=$1 OUT=$2 IN=${3:-/dev/stdin}
    [[ "$(readlink -f "${IN}" 2>/dev/null || echo "${IN}")" \
       != "$(readlink -f "${OUT}" 2>/dev/null || echo "${OUT}")" ]] \
        || die "mk_id_list: input and output are the same file (${OUT})"
    awk 'NR==FNR{want[$1]=1; next} ($2 in want){print $1"\t"$2}' \
        "${IN}" "${FAM}" > "${OUT}"
}

norm_id_list() {
    local IN=$1 FAM=$2 OUT=$3
    grep -v '^#' "${IN}" | awk '{print $NF}' > "${OUT}.tmp"
    mk_id_list "${FAM}" "${OUT}" "${OUT}.tmp"
    rm -f "${OUT}.tmp"
}

# ---------------------------------------------------------------------------
# classify_missnp <missnp> <A.bim> <B.bim> <detail out>
#
# Prints "<class>\t<count>" per class on stdout and writes the per-variant
# detail to <detail out>. B.bim must be the ORIGINAL B, never the flipped copy,
# or every strand conflict reads back as an allele mismatch.
# ---------------------------------------------------------------------------
classify_missnp() {
    local ML=$1 AB=$2 BB=$3 OUT=$4
    awk -F'\t' -v OUT="${OUT}" '
        function c(x){ return (x=="A")?"T":(x=="T")?"A":(x=="C")?"G":(x=="G")?"C":"N" }
        function srt(x,y){ return (x<y) ? x y : y x }
        BEGIN { print "ID\tA_alleles\tB_alleles\tclass" > OUT }
        FNR==1 { fn++ }
        fn==1  { m[$1]=1; next }
        fn==2  { if ($2 in m) { A1[$2]=toupper($5); A2[$2]=toupper($6) } next }
        fn==3  { if (!($2 in m)) next
                 b1=toupper($5); b2=toupper($6); a1=A1[$2]; a2=A2[$2]
                 if (a1=="" && a2=="") { k="not_in_A" }
                 else if (srt(a1,a2) == srt(b1,b2))       k="same_alleles"
                 else if (srt(a1,a2) == srt(c(b1),c(b2))) k="strand_flip"
                 else                                     k="allele_mismatch"
                 n[k]++
                 printf "%s\t%s/%s\t%s/%s\t%s\n", $2, a1, a2, b1, b2, k > OUT }
        END    { for (x in n) printf "%s\t%d\n", x, n[x] }
    ' "${ML}" "${AB}" "${BB}"
}

panel_path() { printf '%s/panels/%s' "${WORK_DIR}" "$1"; }

# Label normalisation used wherever a population name from one panel has to be
# matched against a name from another. Inlined into awk as `nz()` where needed.
norm_label() { printf '%s' "$1" | tr -d ' _.-' | tr '[:upper:]' '[:lower:]'; }

# IIDs of pseudo-haploid individuals, written by 01. May legitimately be empty.
export MF_PSEUDOHAPLOID="${META_DIR}/pseudohaploid_iids.txt"

# Optional, <panel><TAB><IID> per line, # for comments. Individuals removed
# from one named panel in 01 for reasons the automated filters cannot see. The
# panel column matters: 23 of these IIDs exist in both 1000G and SGDP and only
# the SGDP copy is dropped.
export MF_PANEL_EXCLUDE="${META_DIR}/panel_exclude_iids.txt"

# ---------------------------------------------------------------------------
# Make a stage that dies under `set -e` say so.
#
# On 2026-09-02 01_prepare.sh was killed by a `grep` that selected no lines
# inside a pipeline, under `set -o pipefail`. The log simply stopped after the
# heterozygosity summary. 02 then failed on a missing kz.mergeable.snplist and
# reported only that, so the actual cause was two stages and one wrong guess
# away from the message. A pipeline that runs for a day cannot afford that.
#
# `set -E` propagates the trap into functions and subshells, which is where all
# the work happens - without it the trap never fires for anything inside
# harmonise(), merge_group(), run_group() and so on.
#
# Applied only when the caller actually has `set -e`. check_inputs.sh runs with
# `set -uo pipefail` on purpose so it can report every problem instead of
# stopping at the first, and an ERR trap there would be noise.
case $- in
    *e*)
        set -E
        trap 'rc=$?; printf "\n[%s] STAGE FAILED\n  file:    %s\n  line:    %s\n  status:  %s\n  command: %s\n\n" \
              "$(date +"%F %T")" "${BASH_SOURCE[0]}" "${LINENO}" "${rc}" "${BASH_COMMAND}" >&2' ERR
        ;;
esac

# Only ever created UNDER an already-verified PROJECT. Never creates PROJECT
# itself: see the note at the top of this file for what that cost once.
mkdir -p "${RES_DIR}" "${WORK_DIR}" "${LOG_DIR}" "${META_DIR}"

# ---------------------------------------------------------------------------
# admix_cache_guard <bfile prefix> <run root>
#
# ADMIXTURE runs are cached by "does the .Q exist", which is only valid while
# the INPUT has not changed. It did change here: 09's runs were left over from
# a pass when unrelated.keep held 203 individuals, the current input holds 218,
# and the stale .Q files were reused. The alignment caught it on a row count.
# Had the counts happened to match, Q would have been silently attached to the
# wrong people and the figure would have looked fine.
#
# So the input is fingerprinted - sample count, variant count, and the .fam and
# .bim checksums - and any run directory built on a different fingerprint is
# removed rather than reused.
# ---------------------------------------------------------------------------
admix_cache_guard() {
    local BF=$1 ROOT=$2
    local FP="${ROOT}/input.fingerprint" NOW
    NOW=$(printf '%s %s %s %s\n' \
            "$(nlines "${BF}.fam")" "$(nlines "${BF}.bim")" \
            "$(md5sum "${BF}.fam" | cut -d" " -f1)" \
            "$(md5sum "${BF}.bim" | cut -d" " -f1)")
    local NDIR NQ REASON=""
    NDIR=$(find "${ROOT}/runs" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | wc -l)

    if [[ "${NDIR}" -gt 0 ]]; then
        if [[ ! -f "${FP}" ]]; then
            # Runs with no fingerprint predate this check, so nothing is known
            # about the input they were built on. The first version treated a
            # missing fingerprint as "cache is fine", which is backwards and is
            # why 60 stale directories survived and the alignment kept failing.
            REASON="no fingerprint: provenance of the cached runs is unknown"
        elif [[ "$(cat "${FP}")" != "${NOW}" ]]; then
            REASON="input changed (was: $(cat "${FP}") / now: ${NOW})"
        else
            # Fingerprints can only vouch for what they recorded. The invariant
            # that actually matters is checked directly: one row per individual.
            # -print -quit, NOT `| head -1`: head closes the pipe after the
            # first line, find dies of SIGPIPE, and under `set -o pipefail`
            # that 141 trips the ERR trap and kills the stage. find stops
            # itself instead, so nothing is left writing to a closed pipe.
            NQ=$(find "${ROOT}/runs" -maxdepth 2 -name '*.Q' -print -quit 2>/dev/null)
            if [[ -n "${NQ}" ]] && [[ "$(nlines "${NQ}")" -ne "$(nlines "${BF}.fam")" ]]; then
                REASON="a cached .Q has $(nlines "${NQ}") rows, the input has $(nlines "${BF}.fam") individuals"
            fi
        fi
    fi

    if [[ -n "${REASON}" ]]; then
        log "  discarding ${NDIR} cached ADMIXTURE run(s): ${REASON}"
        rm -rf "${ROOT}/runs"
    elif [[ "${NDIR}" -gt 0 ]]; then
        log "  reusing ${NDIR} cached ADMIXTURE run(s); input fingerprint matches"
    fi
    mkdir -p "${ROOT}/runs"
    printf '%s\n' "${NOW}" > "${FP}"
}

