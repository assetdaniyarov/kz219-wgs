#!/usr/bin/env python3
"""
popgen_tables.py — the tables that used to be produced in a notebook.

Every number that goes into the manuscript is written here, by the stage that
produced the underlying files, so the figures can be drawn from TSVs and there
is no interactive step between plink's output and the text.

Subcommands
-----------
  pca-table        .eigenvec + pop_master        -> labelled PC table
  score-table      .sscore   + pop_master        -> labelled projected PC table
  pca-compare      two PC tables                 -> per-axis |r| over shared IDs
  fst-tables       fst.fst.summary               -> pairs with CI, square matrix
  fst-focus        fst.fst.summary               -> one population vs all, ranked
  fst-heatmap      fst_matrix.tsv                -> seriated matrix + heatmap
  proj-shrinkage   scores.tsv                    -> matched pseudohaploid/diploid pairs
  fst-concordance  two fst.fst.summary           -> shared pairs, r, slope, bias
  roh-summary      .hom + .hom.indiv             -> FROH, length classes, by pop
  kin-by-group     .kin0 + labels                -> within vs between kinship
  roh-sensitivity  several .hom.indiv            -> FROH per parameter set + rho
  froh-vs-f        FROH + .het + .ibc            -> correlations
  perm-assoc       PC table + labels             -> ANOVA F with permutation p

Requires numpy and nothing else. No pandas, and deliberately no scipy: the
plink2 conda environment on the analysis host has numpy 1.19.5 under Python
3.6.8 and no scipy, and adding a dependency to a working environment for one
function is not worth it. Written to run on Python 3.6.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# small readers. plink2 writes tab-separated files with a '#'-prefixed header;
# plink1.9 writes whitespace-aligned ones. Population labels can contain spaces
# ("Mbuti Pygmy"), so anything TSV must be split on tabs and nothing else.
# ---------------------------------------------------------------------------
def read_tsv(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n").rstrip("\r")
            if ln == "":
                continue
            rows.append(ln.split("\t"))
    if not rows:
        raise ValueError(f"{path}: empty")
    hdr = [c.lstrip("#").strip() for c in rows[0]]
    return hdr, rows[1:]


def read_ws(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            f = ln.split()
            if f:
                rows.append(f)
    if not rows:
        raise ValueError(f"{path}: empty")
    hdr = [c.lstrip("#").strip() for c in rows[0]]
    return hdr, rows[1:]


def idx(hdr, *names):
    """Index of the first header in *names. Raises with the real header shown."""
    for n in names:
        if n in hdr:
            return hdr.index(n)
    raise KeyError(f"none of {names} in header {hdr}")


def write_tsv(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    print(f"  wrote {path} ({len(rows)} rows)", file=sys.stderr)


def load_meta(path):
    """pop_master.tsv -> {IID: {Nation, Population, Source, Platform}}"""
    hdr, rows = read_tsv(path)
    i_id = idx(hdr, "IID")
    out = {}
    for r in rows:
        if len(r) <= i_id:
            continue
        d = {hdr[j]: (r[j] if j < len(r) else "") for j in range(len(hdr))}
        out[r[i_id]] = d
    return out


def pc_columns(hdr):
    """Indices and names of PC columns, in order. Handles PC1 and PC1_AVG."""
    pcs = []
    for j, name in enumerate(hdr):
        base = name[:-4] if name.endswith("_AVG") else name
        if base.startswith("PC") and base[2:].isdigit():
            pcs.append((int(base[2:]), j, base))
    pcs.sort()
    return [(j, name) for _, j, name in pcs]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def pearson(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def spearman(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    return pearson(rankdata(a[ok]), rankdata(b[ok]))


def rankdata(x):
    """Average ranks, ties shared. Avoids a scipy import for one function."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), float)
    ranks[order] = np.arange(1, len(x) + 1)
    # average tied ranks
    i = 0
    xs = x[order]
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return ranks


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------
def cmd_pca_table(a):
    hdr, rows = read_tsv(a.eigenvec)
    i_id = idx(hdr, "IID")
    pcs = pc_columns(hdr)[: a.npc]
    if not pcs:
        raise SystemExit(f"{a.eigenvec}: no PC columns found in {hdr}")
    meta = load_meta(a.meta)
    metacols = ["Nation", "Population", "Source", "Platform"]

    out = []
    unlabelled = 0
    for r in rows:
        iid = r[i_id]
        m = meta.get(iid)
        if m is None:
            unlabelled += 1
            m = {}
        out.append([iid] + [m.get(c, "NA") for c in metacols] +
                   [r[j] for j, _ in pcs])
    if unlabelled:
        print(f"  WARNING: {unlabelled} individuals have no pop_master row",
              file=sys.stderr)
    write_tsv(a.out, ["IID"] + metacols + [n for _, n in pcs], out)

    if a.eigenval and a.varexp_out:
        ev = [fnum(l.split()[0]) for l in open(a.eigenval) if l.strip()]
        n = a.n_samples if a.n_samples else len(rows)
        # For a variance-standardised GRM the trace equals the sample size, so
        # eigenvalue / N is the share of TOTAL genotypic variance. plink2 writes
        # only N_PC eigenvalues, so eigenvalue / sum(eigenvalues) would be the
        # share of the retained PCs and would overstate every axis several-fold.
        write_tsv(a.varexp_out, ["PC", "eigenvalue", "pct_variance"],
                  [[f"PC{k+1}", f"{v:.6g}", f"{100.0 * v / n:.4f}"]
                   for k, v in enumerate(ev)])


def cmd_score_table(a):
    hdr, rows = read_tsv(a.sscore)
    i_id = idx(hdr, "IID")
    pcs = pc_columns(hdr)[: a.npc]
    if not pcs:
        raise SystemExit(f"{a.sscore}: no PC columns found in {hdr}")
    # plink2 --score reports the mean over scored alleles, so these coordinates
    # are NOT on the same scale as the reference .eigenvec. Every individual in
    # this file went through the same --score call, so the geometry within the
    # file is consistent; never plot part of a figure from .eigenvec and the
    # rest from .sscore.
    try:
        i_miss = idx(hdr, "NMISS_ALLELE_CT")
    except KeyError:
        i_miss = None
    meta = load_meta(a.meta)
    metacols = ["Nation", "Population", "Source", "Platform"]
    out = []
    for r in rows:
        m = meta.get(r[i_id], {})
        extra = [r[i_miss]] if i_miss is not None else []
        out.append([r[i_id]] + [m.get(c, "NA") for c in metacols] + extra +
                   [r[j] for j, _ in pcs])
    cols = ["IID"] + metacols + (["NMISS_ALLELE_CT"] if i_miss is not None else []) \
           + [n for _, n in pcs]
    write_tsv(a.out, cols, out)


def cmd_pca_compare(a):
    """|r| per axis between two PC tables, over the individuals they share.

    Eigenvector signs are arbitrary, so the absolute correlation is the
    quantity with meaning. This is what turns "capping the cohort does not
    change the axes" from an assertion into a number.
    """
    def load(path):
        hdr, rows = read_tsv(path)
        i_id = idx(hdr, "IID")
        pcs = pc_columns(hdr)
        d = {}
        for r in rows:
            d[r[i_id]] = [fnum(r[j]) for j, _ in pcs]
        return d, [n for _, n in pcs]

    A, na = load(a.a)
    B, nb = load(a.b)
    shared = sorted(set(A) & set(B))
    if len(shared) < 3:
        raise SystemExit(f"only {len(shared)} shared individuals")
    npc = min(a.npc, len(na), len(nb))
    rows = []
    for k in range(npc):
        x = [A[i][k] for i in shared]
        y = [B[i][k] for i in shared]
        r = pearson(x, y)
        rows.append([f"PC{k+1}", len(shared), f"{r:.6f}", f"{abs(r):.6f}"])
    write_tsv(a.out, ["PC", "n_shared", "r", "abs_r"], rows)
    worst = min(float(r[3]) for r in rows)
    print(f"  lowest |r| over PC1-PC{npc}: {worst:.4f}", file=sys.stderr)

    # ---- subspace agreement -------------------------------------------------
    # Per-axis correlation is the right statistic for PC1 and PC2, whose
    # eigenvalues are well separated. It is the WRONG statistic further down.
    # When eigenvalues are close together - groupB has PC2 at 0.76% and PC3 at
    # 0.57% of variance - the individual eigenvectors are not identified, only
    # the subspace they span is. Any perturbation of the sample rotates them
    # into each other, so a low per-axis r there is expected geometry, not
    # disagreement, and a reviewer asking "did capping move PC4" is asking a
    # question the per-axis number cannot answer.
    #
    # The quantity that IS identified is the overlap between the subspaces
    # spanned by the first k axes: the mean squared canonical correlation
    # (Krzanowski 1979), which is 1 when the two runs span the same space
    # however the axes inside it are labelled or rotated, and 0 when they are
    # orthogonal. Principal angles are reported beside it.
    X = np.array([A[i][:npc] for i in shared], float)
    Y = np.array([B[i][:npc] for i in shared], float)
    X -= X.mean(axis=0)
    Y -= Y.mean(axis=0)
    srows = []
    for k in sorted({2, 4, npc} & set(range(1, npc + 1))):
        Qa, _ = np.linalg.qr(X[:, :k])
        Qb, _ = np.linalg.qr(Y[:, :k])
        sv = np.clip(np.linalg.svd(Qa.T @ Qb, compute_uv=False), -1.0, 1.0)
        overlap = float((sv ** 2).sum() / k)
        angles = np.degrees(np.arccos(sv))
        srows.append([k, f"{overlap:.6f}", f"{angles.max():.2f}",
                      " ".join(f"{v:.1f}" for v in angles)])
    write_tsv(a.out_subspace, ["n_axes", "subspace_overlap",
                               "max_principal_angle_deg", "principal_angles_deg"],
              srows)
    print("  subspace overlap (1 = same space, however the axes inside it are "
          "rotated):", file=sys.stderr)
    for r in srows:
        print(f"    first {r[0]} axes: {r[1]}  (largest principal angle "
              f"{r[2]} deg)", file=sys.stderr)


# ---------------------------------------------------------------------------
# Fst
# ---------------------------------------------------------------------------
def read_fst_summary(path):
    hdr, rows = read_tsv(path)
    i1 = idx(hdr, "POP1")
    i2 = idx(hdr, "POP2")
    ifst = idx(hdr, "HUDSON_FST", "WC_FST", "FST")
    try:
        ise = idx(hdr, "SE")
    except KeyError:
        ise = None
    out = {}
    for r in rows:
        key = tuple(sorted((r[i1], r[i2])))
        out[key] = (fnum(r[ifst]), fnum(r[ise]) if ise is not None else float("nan"))
    return out


def cmd_fst_tables(a):
    fst = read_fst_summary(a.summary)
    ns = {}
    if a.pop_n:
        hdr, rows = read_tsv(a.pop_n)
        i_p, i_n = idx(hdr, "POP"), idx(hdr, "N")
        ns = {r[i_p]: r[i_n] for r in rows}

    z = 1.959963985
    pairs = []
    for (p1, p2), (v, se) in sorted(fst.items()):
        lo = v - z * se if np.isfinite(se) else float("nan")
        hi = v + z * se if np.isfinite(se) else float("nan")
        # Hudson's estimator is unbiased and can take small negative values when
        # the true divergence is at or near zero. Truncating at zero would make
        # every such estimate look positive, so the raw value is reported and a
        # flag says whether the interval covers zero.
        pairs.append([p1, p2, ns.get(p1, "NA"), ns.get(p2, "NA"),
                      f"{v:.8g}", f"{se:.8g}", f"{lo:.8g}", f"{hi:.8g}",
                      "yes" if (np.isfinite(lo) and lo <= 0.0 <= hi) else "no"])
    write_tsv(a.out_pairs,
              ["POP1", "POP2", "N1", "N2", "HUDSON_FST", "SE",
               "CI95_LO", "CI95_HI", "CI_COVERS_ZERO"], pairs)

    pops = sorted({p for k in fst for p in k})
    mat = [[""] * (len(pops) + 1) for _ in pops]
    for i, p in enumerate(pops):
        mat[i][0] = p
        for j, q in enumerate(pops):
            mat[i][j + 1] = "0" if p == q else f"{fst[tuple(sorted((p, q)))][0]:.6g}"
    write_tsv(a.out_matrix, ["POP"] + pops, mat)


def cmd_fst_concordance(a):
    """Same population pairs, two panels. The platform control.

    groupA is WGS on both sides; groupB puts 30x WGS cohort genotypes against an
    array reference. If the array/WGS mixture were inflating the cohort's
    distances, the groupB values involving KAZ would sit above the groupA ones
    while the reference-vs-reference pairs agreed. This table is what lets that
    be checked instead of argued.
    """
    A = read_fst_summary(a.a)
    B = read_fst_summary(a.b)
    shared = sorted(set(A) & set(B))
    if not shared:
        raise SystemExit("the two summaries share no population pair")
    rows, xa, xb = [], [], []
    for k in shared:
        va, vb = A[k][0], B[k][0]
        involves = "yes" if a.focus and a.focus in k else "no"
        rows.append([k[0], k[1], f"{va:.8g}", f"{vb:.8g}",
                     f"{vb - va:.8g}",
                     f"{(vb / va):.4f}" if va not in (0.0,) and np.isfinite(va) else "NA",
                     involves])
        xa.append(va)
        xb.append(vb)
    write_tsv(a.out, ["POP1", "POP2", "FST_" + a.label_a, "FST_" + a.label_b,
                      "DIFF", "RATIO", "INVOLVES_" + (a.focus or "FOCUS")], rows)

    xa, xb = np.array(xa), np.array(xb)
    ok = np.isfinite(xa) & np.isfinite(xb)
    r = pearson(xa, xb)
    slope = float(np.polyfit(xa[ok], xb[ok], 1)[0]) if ok.sum() >= 3 else float("nan")
    stats = [["n_shared_pairs", ok.sum()],
             ["pearson_r", f"{r:.6f}"],
             ["slope_b_on_a", f"{slope:.6f}"],
             ["mean_abs_diff", f"{np.mean(np.abs(xb[ok] - xa[ok])):.8g}"]]
    if a.focus:
        f = np.array([a.focus in k for k in shared])
        for name, sel in (("focus", f & ok), ("other", (~f) & ok)):
            if sel.sum():
                stats.append([f"mean_diff_{name}",
                              f"{np.mean(xb[sel] - xa[sel]):.8g}"])
        stats.append(["interpretation",
                      "a focus-specific offset larger than the other-pair offset "
                      "is the signature of a platform effect on the cohort"])
    write_tsv(a.out_stats, ["key", "value"], stats)


def _read_matrix(path):
    hdr, rows = read_tsv(path)
    pops = hdr[1:]
    M = np.array([[fnum(v) for v in r[1:]] for r in rows], float)
    labs = [r[0] for r in rows]
    return labs, pops, M


def cmd_fst_focus(a):
    """One population against every other, sorted nearest first.

    This is the sentence the main text needs - "the cohort is closest to X,
    then Y" - and a 595-row pair list does not give it. The neighbour ranking
    is only a claim if consecutive entries can be told apart, so each row also
    carries whether its interval overlaps the one above it.
    """
    fst = read_fst_summary(a.summary)
    ns = {}
    if a.pop_n:
        hdr, rows = read_tsv(a.pop_n)
        i_p, i_n = idx(hdr, "POP"), idx(hdr, "N")
        ns = {r[i_p]: r[i_n] for r in rows}

    z = 1.959963985
    items = []
    for (p1, p2), (v, se) in fst.items():
        if a.focus not in (p1, p2):
            continue
        other = p2 if p1 == a.focus else p1
        items.append((v, se, other))
    if not items:
        raise SystemExit(f"'{a.focus}' appears in no pair of {a.summary}")
    items.sort(key=lambda t: t[0])

    out = []
    prev_hi = None
    for rank, (v, se, other) in enumerate(items, 1):
        lo, hi = v - z * se, v + z * se
        sep = "NA" if prev_hi is None else ("yes" if lo > prev_hi else "no")
        out.append([rank, other, ns.get(other, "NA"), f"{v:.8g}", f"{se:.8g}",
                    f"{lo:.8g}", f"{hi:.8g}", sep])
        prev_hi = hi
    write_tsv(a.out, ["rank", "population", "N", f"FST_vs_{a.focus}", "SE",
                      "CI95_LO", "CI95_HI", "SEPARATED_FROM_PREVIOUS"], out)
    near = ", ".join(r[1] for r in out[:5])
    print(f"  {a.focus} is closest to: {near}", file=sys.stderr)


def _pcoa_order(M):
    """1-D seriation by the first principal coordinate of the distance matrix.

    Alphabetical rows and columns make a heatmap that shows nothing. Ordering
    by the leading principal coordinate puts similar populations next to each
    other, which is what makes the block structure visible. Pure numpy, so no
    dependency and no randomness.
    """
    n = M.shape[0]
    D2 = M ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh((B + B.T) / 2.0)
    return np.argsort(V[:, np.argmax(w)])


def _cluster_order(M):
    """Average-linkage seriation with optimal leaf ordering.

    Preferred over _pcoa_order because one distant outgroup takes over the
    first principal coordinate: with an African outgroup in the matrix, PCoA
    axis 1 becomes Africa-versus-everything and the ordering inside Eurasia -
    the part the figure is about - is left to noise. Hierarchical clustering
    keeps the fine structure. Falls back to PCoA when scipy is absent.
    """
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
        from scipy.spatial.distance import squareform
    except ImportError:
        print("  scipy absent, falling back to PCoA seriation", file=sys.stderr)
        return _pcoa_order(M)
    # Hudson Fst is unbiased and goes slightly negative for pairs with no real
    # divergence. Clipped to 0 for the clustering only; the values written out
    # and plotted are untouched.
    D = np.clip((M + M.T) / 2.0, 0.0, None)
    np.fill_diagonal(D, 0.0)
    d = squareform(D, checks=False)
    Z = optimal_leaf_ordering(linkage(d, method="average"), d)
    return np.asarray(leaves_list(Z))


def cmd_fst_heatmap(a):
    labs, pops, M = _read_matrix(a.matrix)
    if labs != pops:
        raise SystemExit("matrix is not square/symmetric in its labels")

    keep = [i for i, l in enumerate(labs)
            if l not in set(a.exclude.split()) if a.exclude]
    if a.exclude:
        labs = [labs[i] for i in keep]
        M = M[np.ix_(keep, keep)]

    if a.order_by == "alpha":
        order = np.argsort(labs)
    elif a.order_by == "focus":
        if a.focus not in labs:
            raise SystemExit(f"focus '{a.focus}' not in the matrix")
        order = np.argsort(M[labs.index(a.focus)])
    elif a.order_by == "pcoa":
        order = _pcoa_order(M)
    else:
        order = _cluster_order(M)

    L = [labs[i] for i in order]
    O = M[np.ix_(order, order)]
    write_tsv(a.out_ordered, ["POP"] + L,
              [[L[i]] + [f"{O[i, j]:.6g}" for j in range(len(L))]
               for i in range(len(L))])

    if not a.out_png:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not installed; the ordered matrix was written and "
              "can be plotted anywhere. To render here: pip install matplotlib",
              file=sys.stderr)
        return

    n = len(L)
    fig, ax = plt.subplots(figsize=(max(6, 0.32 * n + 2), max(5, 0.32 * n + 1.5)))
    Mp = np.array(O, float)
    np.fill_diagonal(Mp, np.nan)
    # One distant outgroup flattens the colour scale: with AFR at 0.11 and
    # every Eurasian pair below 0.02, the figure is one bright row and an
    # otherwise uniform block. --vmax clips the scale; --exclude drops the
    # outgroup from the panel entirely.
    hi = float(np.nanmax(Mp)); p90 = float(np.nanpercentile(Mp, 90))
    if a.vmax is None and p90 > 0 and hi > 5 * p90:
        print(f"  NOTE: the scale is dominated by one distant group "
              f"(max {hi:.4f} vs 90th percentile {p90:.4f}). Consider "
              f"--vmax {p90:.3f} or --exclude <outgroup>.", file=sys.stderr)
    im = ax.imshow(Mp, cmap="viridis", interpolation="nearest", vmax=a.vmax)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(L, rotation=90, fontsize=7)
    ax.set_yticklabels(L, fontsize=7)
    if a.focus in L:
        k = L.index(a.focus)
        for lbl in (ax.get_xticklabels()[k], ax.get_yticklabels()[k]):
            lbl.set_fontweight("bold")
        ax.axhline(k, color="w", lw=0.8, alpha=0.6)
        ax.axvline(k, color="w", lw=0.8, alpha=0.6)
    fig.colorbar(im, ax=ax, shrink=0.75, label="Hudson $F_{ST}$")
    ax.set_title(a.title or "Pairwise $F_{ST}$")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = a.out_png if a.out_png.endswith(ext) else \
            a.out_png.rsplit(".", 1)[0] + "." + ext
        fig.savefig(out, dpi=300)
        print(f"  wrote {out}", file=sys.stderr)
    plt.close(fig)


def cmd_proj_shrinkage(a):
    """Matched pseudo-haploid vs diploid groups on the same projected axes.

    Comparing all 46 AADR individuals against all 915 diploids confounds
    shrinkage with geography: the AADR set is entirely Central Asian and
    Siberian, near the centre of a Eurasian PC space, while the diploid set
    includes populations on the rim. Three AADR groups have a diploid
    counterpart in the same projection - Kazakh, Uyghur, Hazara - and those
    pairs isolate the effect.

    Two numbers per pair. The RATIO of centroid norms is how much the
    pseudo-haploid group is pulled toward the origin. The ANGLE between the two
    centroids is whether it is pulled toward the origin or somewhere else: a
    small angle means the projected position points the right way and is only
    shrunk, which is exactly the qualitative claim the figure needs to support.
    """
    hdr, rows = read_tsv(a.scores)
    i_id = idx(hdr, "IID")
    i_nat = idx(hdr, "Nation")
    i_pop = hdr.index("Population") if "Population" in hdr else None
    pcs = pc_columns(hdr)[: a.npc]
    if not pcs:
        raise SystemExit(f"{a.scores}: no PC columns")

    byname = {}
    for r in rows:
        v = np.array([fnum(r[j]) for j, _ in pcs])
        for key in filter(None, [r[i_nat],
                                 r[i_pop] if i_pop is not None else None]):
            byname.setdefault(key, []).append(v)

    fmiss = {}
    if a.missingness:
        mh, mr = read_tsv(a.missingness)
        j_id, j_fm = idx(mh, "IID"), idx(mh, "F_MISS")
        j_nat = idx(mh, "Nation")
        for r in mr:
            fmiss.setdefault(r[j_nat], []).append(fnum(r[j_fm]))

    out = []
    for spec in a.pairs.split():
        if "=" not in spec:
            continue
        ph, dip = spec.split("=", 1)
        if ph not in byname or dip not in byname:
            print(f"  skipped {spec}: {'both' if ph not in byname and dip not in byname else (ph if ph not in byname else dip)} not found",
                  file=sys.stderr)
            continue
        A = np.vstack(byname[ph])
        B = np.vstack(byname[dip])
        ca, cb = A.mean(axis=0), B.mean(axis=0)
        na, nb = float(np.linalg.norm(ca)), float(np.linalg.norm(cb))
        ratio = na / nb if nb > 0 else float("nan")
        cosang = float(ca @ cb / (na * nb)) if na > 0 and nb > 0 else float("nan")
        ang = float(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
        # Bootstrap over individuals in both groups. With n_ph = 3 the centroid
        # is a 3-point mean and the interval comes out wide - which is the
        # honest answer, not a defect. A ratio of 0.94 with an interval
        # spanning 1.0 says "no shrinkage detected at this sample size", and
        # that is a stronger statement for the figure than an unqualified 0.94.
        rng = np.random.default_rng(a.seed)
        br, bg = [], []
        for _ in range(a.nboot):
            ca_b = A[rng.integers(0, len(A), len(A))].mean(axis=0)
            cb_b = B[rng.integers(0, len(B), len(B))].mean(axis=0)
            na_b, nb_b = np.linalg.norm(ca_b), np.linalg.norm(cb_b)
            if na_b == 0 or nb_b == 0:
                continue
            br.append(na_b / nb_b)
            bg.append(np.degrees(np.arccos(
                np.clip(ca_b @ cb_b / (na_b * nb_b), -1, 1))))
        rlo, rhi = (np.percentile(br, [2.5, 97.5]) if br else (np.nan, np.nan))
        ghi = float(np.percentile(bg, 97.5)) if bg else float("nan")
        mf = (f"{np.mean(fmiss[ph]):.4f}" if ph in fmiss else "NA")
        out.append([ph, len(A), f"{mf}", dip, len(B), f"{na:.6f}", f"{nb:.6f}",
                    f"{ratio:.4f}", f"{rlo:.4f}", f"{rhi:.4f}",
                    f"{ang:.1f}", f"{ghi:.1f}",
                    "yes" if (rlo <= 1.0 <= rhi) else "no"])
    if not out:
        raise SystemExit("no matched pair could be resolved")
    # Mean F_MISS travels with each pair, because the three groups that HAVE a
    # diploid counterpart are the near-complete ones. A small ratio measured on
    # them says pseudo-haploidy alone costs little; it says nothing about the
    # Kalmyk, Teleut and Tomsk Tatar samples at F_MISS 0.7-0.8, which have no
    # counterpart and stay qualitative.
    write_tsv(a.out, ["pseudohaploid_group", "n_ph", "mean_F_MISS_ph",
                      "diploid_group", "n_dip",
                      "centroid_norm_ph", "centroid_norm_dip",
                      "norm_ratio", "ratio_CI_LO", "ratio_CI_HI",
                      "angle_deg", "angle_CI_HI", "RATIO_CI_COVERS_1"], out)
    r = np.array([float(x[7]) for x in out])
    g = np.array([float(x[10]) for x in out])
    ncov = sum(1 for x in out if x[12] == "yes")
    print(f"  matched pairs: norm ratio {r.mean():.3f} (range {r.min():.3f}-"
          f"{r.max():.3f}), angle {g.mean():.1f} deg (max {g.max():.1f})",
          file=sys.stderr)
    print(f"  {ncov} of {len(out)} bootstrap intervals cover a ratio of 1.",
          file=sys.stderr)
    print("  Small angle = the projected position points the right way. A ratio "
          "interval covering 1 means no shrinkage is DETECTABLE at this n, not "
          "that none exists - report the interval, never the point estimate "
          "alone.", file=sys.stderr)


# ---------------------------------------------------------------------------
# ROH
# ---------------------------------------------------------------------------
def parse_classes(spec):
    """'500-1000 1000-2000 8000-' -> [(500,1000),(1000,2000),(8000,inf)]"""
    out = []
    for tok in spec.split():
        lo, _, hi = tok.partition("-")
        out.append((float(lo), float(hi) if hi else float("inf")))
    return out


def cmd_roh_summary(a):
    classes = parse_classes(a.classes)
    meta = load_meta(a.meta) if a.meta else {}

    hdr, rows = read_ws(a.indiv)
    i_id, i_nseg, i_kb = idx(hdr, "IID"), idx(hdr, "NSEG"), idx(hdr, "KB")
    total = {r[i_id]: (int(float(r[i_nseg])), fnum(r[i_kb])) for r in rows}

    per_class = {i: {c: 0.0 for c in range(len(classes))} for i in total}
    nseg_class = {i: {c: 0 for c in range(len(classes))} for i in total}
    long_kb = {i: 0.0 for i in total}
    long_n = {i: 0 for i in total}
    if Path(a.hom).exists():
        hh, hr = read_ws(a.hom)
        j_id, j_kb = idx(hh, "IID"), idx(hh, "KB")
        for r in hr:
            iid = r[j_id]
            kb = fnum(r[j_kb])
            if iid not in per_class:
                continue
            if kb >= a.long_kb:
                long_kb[iid] += kb
                long_n[iid] += 1
            for c, (lo, hi) in enumerate(classes):
                if lo <= kb < hi:
                    per_class[iid][c] += kb
                    nseg_class[iid][c] += 1
                    break

    denom_kb = a.autosome_bp / 1000.0
    cls_names = [f"{int(lo)}kb_plus" if hi == float("inf")
                 else f"{int(lo)}_{int(hi)}kb" for lo, hi in classes]

    out = []
    for iid in sorted(total):
        nseg, kb = total[iid]
        m = meta.get(iid, {})
        row = [iid, m.get("Nation", "NA"), m.get("Population", "NA"),
               m.get("Source", "NA"), m.get("Platform", "NA"),
               nseg, f"{kb:.3f}", f"{kb / denom_kb:.6f}",
               long_n[iid], f"{long_kb[iid] / denom_kb:.6f}"]
        for c in range(len(classes)):
            row += [nseg_class[iid][c], f"{per_class[iid][c] / denom_kb:.6f}"]
        out.append(row)
    # FROH_LONG counts only runs at or above --long-kb. This is the quantity the
    # 1/32 and 1/64 benchmarks refer to (McQuillan 2008, Kirin 2010). Total FROH
    # on a WGS panel with a 500 kb minimum is dominated by ancient background
    # homozygosity that everyone carries, so comparing IT against 1/32 says 59%
    # of a Kazakh cohort are children of first cousins, which is nonsense.
    cols = ["IID", "Nation", "Population", "Source", "Platform",
            "NSEG", "TOTAL_KB", "FROH", "NSEG_LONG", "FROH_LONG"]
    for n in cls_names:
        cols += [f"NSEG_{n}", f"FROH_{n}"]
    write_tsv(a.out, cols, out)

    # by population
    groups = {}
    for r in out:
        key = r[2] if a.by == "Population" else r[1]
        groups.setdefault(key, []).append((float(r[7]), float(r[9])))
    grows = []
    for k in sorted(groups):
        v = np.array([x[0] for x in groups[k]], float)
        w = np.array([x[1] for x in groups[k]], float)
        grows.append([k, len(v), f"{v.mean():.6f}", f"{np.median(v):.6f}",
                      f"{v.std(ddof=1) if len(v) > 1 else 0.0:.6f}",
                      f"{v.min():.6f}", f"{v.max():.6f}",
                      f"{w.mean():.6f}", f"{np.median(w):.6f}",
                      f"{float((w > a.high_froh).mean()):.4f}"])
    write_tsv(a.out_by_group,
              [a.by, "n", "FROH_mean", "FROH_median", "FROH_sd",
               "FROH_min", "FROH_max", "FROH_LONG_mean", "FROH_LONG_median",
               f"frac_FROH_LONG_gt_{a.high_froh}"], grows)


def cmd_roh_sensitivity(a):
    denom_kb = a.autosome_bp / 1000.0
    sets = {}
    for path in a.indiv:
        name = Path(path).name
        for suf in (".hom.indiv", ".indiv"):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        hdr, rows = read_ws(path)
        i_id, i_kb = idx(hdr, "IID"), idx(hdr, "KB")
        sets[name] = {r[i_id]: fnum(r[i_kb]) / denom_kb for r in rows}
    if a.reference not in sets:
        raise SystemExit(f"reference '{a.reference}' not among {sorted(sets)}")

    ids = sorted(set.intersection(*(set(v) for v in sets.values())))
    names = [a.reference] + [n for n in sorted(sets) if n != a.reference]
    write_tsv(a.out, ["IID"] + [f"FROH_{n}" for n in names],
              [[i] + [f"{sets[n][i]:.6f}" for n in names] for i in ids])

    ref = np.array([sets[a.reference][i] for i in ids])
    srows = []
    for n in names:
        v = np.array([sets[n][i] for i in ids])
        srows.append([n, len(ids), f"{v.mean():.6f}", f"{np.median(v):.6f}",
                      f"{pearson(ref, v):.6f}", f"{spearman(ref, v):.6f}"])
    write_tsv(a.out_stats,
              ["param_set", "n", "FROH_mean", "FROH_median",
               f"pearson_vs_{a.reference}", f"spearman_vs_{a.reference}"], srows)


def cmd_froh_vs_f(a):
    hdr, rows = read_tsv(a.froh)
    i_id, i_f = idx(hdr, "IID"), idx(hdr, "FROH")
    froh = {r[i_id]: fnum(r[i_f]) for r in rows}
    # Correlate against FROH_LONG as well when it is present. Total FROH on a
    # WGS panel with a 500 kb floor is mostly ancient homozygosity that every
    # individual carries, so it has little between-individual variance for an F
    # estimator to track; the long runs are where the signal is.
    froh_long = None
    if "FROH_LONG" in hdr:
        j = hdr.index("FROH_LONG")
        froh_long = {r[i_id]: fnum(r[j]) for r in rows}

    series = {}
    if a.het:
        hdr, rows = read_tsv(a.het)
        i_id, i_v = idx(hdr, "IID"), idx(hdr, "F")
        series["F_plink2_het"] = {r[i_id]: fnum(r[i_v]) for r in rows}
    if a.ibc:
        hdr, rows = read_ws(a.ibc)
        i_id = idx(hdr, "IID")
        for name in ("Fhat1", "Fhat2", "Fhat3"):
            if name in hdr:
                j = hdr.index(name)
                series[name] = {r[i_id]: fnum(r[j]) for r in rows}

    ids = sorted(set(froh).intersection(*[set(s) for s in series.values()])) \
        if series else sorted(froh)
    jcols = ["IID", "FROH"] + (["FROH_LONG"] if froh_long else []) + list(series)
    write_tsv(a.out_join, jcols,
              [[i, f"{froh[i]:.6f}"]
               + ([f"{froh_long[i]:.6f}"] if froh_long else [])
               + [f"{series[k][i]:.6f}" for k in series]
               for i in ids])

    x = np.array([froh[i] for i in ids])
    xl = np.array([froh_long[i] for i in ids]) if froh_long else None
    srows = []
    for name, ser in series.items():
        y = np.array([ser[i] for i in ids])
        row = [name, len(ids), f"{pearson(x, y):.4f}", f"{spearman(x, y):.4f}"]
        row += ([f"{pearson(xl, y):.4f}", f"{spearman(xl, y):.4f}"]
                if xl is not None else ["NA", "NA"])
        srows.append(row)
    write_tsv(a.out_stats,
              ["estimator", "n", "pearson_vs_FROH", "spearman_vs_FROH",
               "pearson_vs_FROH_LONG", "spearman_vs_FROH_LONG"], srows)
    print("  FROH and the method-of-moments F are independent estimates of the "
          "same quantity; a high correlation is the evidence that the "
          "heterozygosity-deficit individuals flagged in 01 are consanguineous "
          "rather than technical failures.", file=sys.stderr)


def cmd_kin_by_group(a):
    """Mean pairwise kinship within versus between labelled groups.

    Written for the tribe analysis. Kazakh tribes are patrilineal and
    exogamous, so the autosomal prediction is that two people of the same tribe
    are no more closely related than two people of different tribes. That is a
    testable statement, and it is much stronger evidence about the marriage
    system than a near-zero Fst on its own, because it uses every pair rather
    than a two-group summary.

    The p-value comes from permuting the labels over individuals, which
    preserves the pairing structure and the kinship distribution exactly and
    only breaks the association with the labels.
    """
    hdr, rows = read_tsv(a.kin)
    i1, i2 = idx(hdr, "IID1"), idx(hdr, "IID2")
    ik = idx(hdr, "KINSHIP")

    lhdr, lrows = read_tsv(a.labels)
    li = idx(lhdr, "IID")
    lj = lhdr.index(a.label_col) if a.label_col in lhdr else 1
    lab = {r[li]: r[lj] for r in lrows if len(r) > lj and r[lj] not in ("", "NA")}

    counts = {}
    for v in lab.values():
        counts[v] = counts.get(v, 0) + 1
    levels = sorted(l for l, c in counts.items() if c >= a.min_n)
    if len(levels) < 2:
        write_tsv(a.out, ["key", "value"], [["n_levels", len(levels)],
                                            ["result", "not enough levels"]])
        return
    keep = {i: lab[i] for i in lab if lab[i] in levels}
    ids = sorted(keep)
    pos = {v: k for k, v in enumerate(ids)}

    ii, jj, kk = [], [], []
    for r in rows:
        x, y = r[i1], r[i2]
        if x in pos and y in pos:
            ii.append(pos[x])
            jj.append(pos[y])
            kk.append(fnum(r[ik]))
    if len(kk) < 10:
        raise SystemExit(f"only {len(kk)} usable pairs; was --make-king-table "
                         f"run without --king-table-filter?")
    ii = np.array(ii)
    jj = np.array(jj)
    kk = np.array(kk, float)
    lv = np.array([levels.index(keep[i]) for i in ids])

    same = lv[ii] == lv[jj]
    if same.sum() == 0 or (~same).sum() == 0:
        raise SystemExit("all pairs fall on one side of the within/between split")
    obs = float(kk[same].mean() - kk[~same].mean())

    rng = np.random.default_rng(a.seed)
    ge = 0
    for _ in range(a.nperm):
        p_lv = rng.permutation(lv)
        s = p_lv[ii] == p_lv[jj]
        if s.sum() == 0 or (~s).sum() == 0:
            continue
        if abs(float(kk[s].mean() - kk[~s].mean())) >= abs(obs):
            ge += 1
    p = (ge + 1) / (a.nperm + 1)

    write_tsv(a.out, ["key", "value"], [
        ["label", a.label_col],
        ["n_individuals", len(ids)],
        ["n_levels", len(levels)],
        ["levels", ",".join(levels)],
        ["n_pairs", len(kk)],
        ["n_pairs_within", int(same.sum())],
        ["n_pairs_between", int((~same).sum())],
        ["mean_kinship_within", f"{kk[same].mean():.6f}"],
        ["mean_kinship_between", f"{kk[~same].mean():.6f}"],
        ["median_kinship_all", f"{np.median(kk):.6f}"],
        ["iqr_kinship_all",
         f"{np.percentile(kk, 25):.6f}-{np.percentile(kk, 75):.6f}"],
        ["difference_within_minus_between", f"{obs:.6f}"],
        ["p_two_sided_permutation", f"{p:.5f}"],
        ["n_permutations", a.nperm],
    ])
    write_tsv(a.out_levels, [a.label_col, "n"], [[l, counts[l]] for l in levels])
    print(f"  within {kk[same].mean():.5f} vs between {kk[~same].mean():.5f}, "
          f"difference {obs:+.5f}, p = {p:.4f}", file=sys.stderr)


# ---------------------------------------------------------------------------
# association of PCs with a categorical label, permutation p-value
# ---------------------------------------------------------------------------
def anova_f(values, groups):
    """One-way ANOVA F. groups is an integer array of group indices."""
    k = groups.max() + 1
    n = len(values)
    grand = values.mean()
    ssb = 0.0
    ssw = 0.0
    for g in range(k):
        v = values[groups == g]
        if len(v) == 0:
            continue
        ssb += len(v) * (v.mean() - grand) ** 2
        ssw += ((v - v.mean()) ** 2).sum()
    dfb, dfw = k - 1, n - k
    if dfb <= 0 or dfw <= 0 or ssw <= 0:
        return float("nan")
    return (ssb / dfb) / (ssw / dfw)


def cmd_perm_assoc(a):
    """Is any PC associated with the label, beyond what shuffling would give?

    Reported instead of an eyeballed scatter plot. With a flat eigenvalue
    spectrum and an exogamous marriage system the expected answer is no, and a
    permutation p-value is what makes that a result rather than an absence.
    """
    hdr, rows = read_tsv(a.pc_table)
    i_id = idx(hdr, "IID")
    pcs = pc_columns(hdr)[: a.npc]

    lhdr, lrows = read_tsv(a.labels)
    li_id = idx(lhdr, "IID")
    li_lab = lhdr.index(a.label_col) if a.label_col in lhdr else 1
    lab = {r[li_id]: r[li_lab] for r in lrows if len(r) > li_lab}

    keep_rows = [r for r in rows if lab.get(r[i_id], "") not in ("", "NA")]
    counts = {}
    for r in keep_rows:
        counts[lab[r[i_id]]] = counts.get(lab[r[i_id]], 0) + 1
    levels = sorted(l for l, c in counts.items() if c >= a.min_n)
    if len(levels) < 2:
        write_tsv(a.out, ["PC", "n", "n_levels", "F", "p_perm"],
                  [["NA", len(keep_rows), len(levels), "NA", "NA"]])
        print(f"  fewer than two levels reach n>={a.min_n}; no test run",
              file=sys.stderr)
        return
    lv_index = {l: i for i, l in enumerate(levels)}
    keep_rows = [r for r in keep_rows if lab[r[i_id]] in lv_index]
    groups = np.array([lv_index[lab[r[i_id]]] for r in keep_rows])

    rng = np.random.default_rng(a.seed)
    out = []
    for k, (j, name) in enumerate(pcs):
        v = np.array([fnum(r[j]) for r in keep_rows])
        obs = anova_f(v, groups)
        ge = 0
        for _ in range(a.nperm):
            if anova_f(v, rng.permutation(groups)) >= obs:
                ge += 1
        p = (ge + 1) / (a.nperm + 1)
        out.append([name, len(v), len(levels), f"{obs:.4f}", f"{p:.5f}",
                    f"{min(1.0, p * len(pcs)):.5f}"])
    write_tsv(a.out, ["PC", "n", "n_levels", "F", "p_perm", "p_bonferroni"], out)
    write_tsv(a.out_levels, [a.label_col, "n"],
              [[l, counts[l]] for l in levels])


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd")
    # add_subparsers(required=True) is Python 3.7+; setting the attribute is the
    # 3.6-compatible form and behaves identically.
    sub.required = True

    s = sub.add_parser("pca-table")
    s.add_argument("--eigenvec", required=True)
    s.add_argument("--meta", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--npc", type=int, default=20)
    s.add_argument("--eigenval")
    s.add_argument("--varexp-out")
    s.add_argument("--n-samples", type=int, default=0)
    s.set_defaults(fn=cmd_pca_table)

    s = sub.add_parser("score-table")
    s.add_argument("--sscore", required=True)
    s.add_argument("--meta", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--npc", type=int, default=20)
    s.set_defaults(fn=cmd_score_table)

    s = sub.add_parser("pca-compare")
    s.add_argument("--a", required=True)
    s.add_argument("--b", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--npc", type=int, default=4)
    s.add_argument("--out-subspace", required=True)
    s.set_defaults(fn=cmd_pca_compare)

    s = sub.add_parser("fst-tables")
    s.add_argument("--summary", required=True)
    s.add_argument("--pop-n")
    s.add_argument("--out-pairs", required=True)
    s.add_argument("--out-matrix", required=True)
    s.set_defaults(fn=cmd_fst_tables)

    s = sub.add_parser("fst-concordance")
    s.add_argument("--a", required=True)
    s.add_argument("--b", required=True)
    s.add_argument("--label-a", default="A")
    s.add_argument("--label-b", default="B")
    s.add_argument("--focus", default="KAZ")
    s.add_argument("--out", required=True)
    s.add_argument("--out-stats", required=True)
    s.set_defaults(fn=cmd_fst_concordance)

    s = sub.add_parser("fst-focus")
    s.add_argument("--summary", required=True)
    s.add_argument("--pop-n")
    s.add_argument("--focus", default="KAZ")
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_fst_focus)

    s = sub.add_parser("fst-heatmap")
    s.add_argument("--matrix", required=True)
    s.add_argument("--out-ordered", required=True)
    s.add_argument("--out-png")
    s.add_argument("--focus", default="KAZ")
    s.add_argument("--title")
    s.add_argument("--order-by", default="cluster",
                   choices=["cluster", "pcoa", "focus", "alpha"])
    s.add_argument("--vmax", type=float, default=None)
    s.add_argument("--exclude", default="")
    s.set_defaults(fn=cmd_fst_heatmap)

    s = sub.add_parser("proj-shrinkage")
    s.add_argument("--scores", required=True)
    s.add_argument("--pairs", required=True)
    s.add_argument("--npc", type=int, default=4)
    s.add_argument("--missingness")
    s.add_argument("--nboot", type=int, default=2000)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_proj_shrinkage)

    s = sub.add_parser("roh-summary")
    s.add_argument("--hom", required=True)
    s.add_argument("--indiv", required=True)
    s.add_argument("--meta")
    s.add_argument("--autosome-bp", type=float, required=True)
    s.add_argument("--classes", default="500-1000 1000-2000 2000-4000 4000-8000 8000-")
    s.add_argument("--by", default="Population", choices=["Population", "Nation"])
    s.add_argument("--high-froh", type=float, default=0.03125)
    s.add_argument("--long-kb", type=float, default=1500.0)
    s.add_argument("--out", required=True)
    s.add_argument("--out-by-group", required=True)
    s.set_defaults(fn=cmd_roh_summary)

    s = sub.add_parser("roh-sensitivity")
    s.add_argument("--indiv", nargs="+", required=True)
    s.add_argument("--reference", required=True)
    s.add_argument("--autosome-bp", type=float, required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--out-stats", required=True)
    s.set_defaults(fn=cmd_roh_sensitivity)

    s = sub.add_parser("froh-vs-f")
    s.add_argument("--froh", required=True)
    s.add_argument("--het")
    s.add_argument("--ibc")
    s.add_argument("--out-join", required=True)
    s.add_argument("--out-stats", required=True)
    s.set_defaults(fn=cmd_froh_vs_f)

    s = sub.add_parser("kin-by-group")
    s.add_argument("--kin", required=True)
    s.add_argument("--labels", required=True)
    s.add_argument("--label-col", default="Tribe")
    s.add_argument("--min-n", type=int, default=10)
    s.add_argument("--nperm", type=int, default=10000)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--out", required=True)
    s.add_argument("--out-levels", required=True)
    s.set_defaults(fn=cmd_kin_by_group)

    s = sub.add_parser("perm-assoc")
    s.add_argument("--pc-table", required=True)
    s.add_argument("--labels", required=True)
    s.add_argument("--label-col", default="Tribe")
    s.add_argument("--npc", type=int, default=10)
    s.add_argument("--min-n", type=int, default=10)
    s.add_argument("--nperm", type=int, default=10000)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--out", required=True)
    s.add_argument("--out-levels", required=True)
    s.set_defaults(fn=cmd_perm_assoc)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
