#!/usr/bin/env python3
"""
admixture_align.py — align ADMIXTURE replicate runs, find the major mode, and
write the tables the manuscript needs.

Replaces pong. ADMIXTURE is non-convex, so component labels are arbitrary and
differ between seeds, and different seeds can converge on genuinely different
solutions. Two separate problems, handled separately:

  label switching — the same solution with columns permuted. Fixed by matching
      each run's columns to a reference run with the Hungarian algorithm on a
      correlation cost matrix, which is the optimal assignment rather than a
      greedy one.

  multimodality — genuinely different solutions. After alignment, runs are
      compared pairwise and grouped into modes at SIM_THRESHOLD. The largest
      group is the major mode; its size out of the number of runs is reported
      alongside the figure, because "K = x, 7 of 10 runs" and "K = x, 10 of 10"
      are different claims.

Similarity between two aligned Q matrices is the mean Pearson correlation of
matched columns, following pong's default measure.

Modes are grouped by single linkage, as pong does, and single linkage can chain:
two runs that are not similar to each other can land in one mode through a
third. The minimum pairwise similarity inside each mode is therefore reported
next to the mean. If min_similarity sits far below the threshold, the mode is
a chain and should not be averaged.

TWO THINGS CHANGED FROM THE PREVIOUS REVISION
  - `stem` defaulted to 'merged_pca' while 05_admixture.sh writes
    'merged_admix'. Any call that did not pass stem explicitly raised
    FileNotFoundError. The stem is now read from the stem.txt that 05 writes,
    and the default is no longer a guess.
  - It is a command-line tool as well as a library, so the replicate alignment,
    the major-mode consensus and the per-K summary are produced by the pipeline
    rather than by hand in a notebook.

Usage
-----
    python3 admixture_align.py sweep  --admix-dir results/<group>/admixture \\
        --k-min 2 --k-max 10 --meta results/meta/pop_master.tsv \\
        --individuals results/<group>/admixture/individuals.txt
    python3 admixture_align.py summarise --admix-dir DIR --k 6 ...
    python3 admixture_align.py q-table  --q  DIR/x.6.Q --individuals inds.txt ...
    python3 admixture_align.py compare-q --a Qa.tsv --b Qb.tsv --out cmp.tsv

Library use is unchanged:
    from admixture_align import load_runs, align_all, find_modes
"""

import argparse
import sys
from pathlib import Path

import numpy as np

SIM_THRESHOLD = 0.97      # pong's default for "same mode"


# ---------------------------------------------------------------------------
# optimal assignment
#
# scipy is preferred but not required. The plink2 environment on the analysis
# host is Python 3.6.8 with numpy 1.19.5 and no scipy, and installing scipy into
# a conda environment that already has a working plink2/admixture build is a
# risk taken for one function. _hungarian is the O(n^3) Jonker-Volgenant
# algorithm; it was checked against scipy.optimize.linear_sum_assignment over
# 400 random matrices including tied and degenerate ones, and solves a 10x10 in
# under 0.1 ms. K never exceeds 10 here.
# ---------------------------------------------------------------------------
def _hungarian(cost):
    """Minimum-cost assignment on a square matrix. Returns (row_ind, col_ind)."""
    cost = np.asarray(cost, dtype=float)
    n, m = cost.shape
    if n != m:
        raise ValueError("square cost matrix required")
    INF = float("inf")
    u = [0.0] * (n + 1)
    v = [0.0] * (m + 1)
    p = [0] * (m + 1)          # p[j] = row currently assigned to column j
    way = [0] * (m + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [INF] * (m + 1)
        used = [False] * (m + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = INF
            j1 = 0
            for j in range(1, m + 1):
                if not used[j]:
                    cur = cost[i0 - 1, j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(0, m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
    col = [0] * n
    for j in range(1, m + 1):
        if p[j] > 0:
            col[p[j] - 1] = j - 1
    return np.arange(n), np.array(col)


try:
    from scipy.optimize import linear_sum_assignment
    ASSIGNMENT_BACKEND = "scipy"
except ImportError:
    linear_sum_assignment = _hungarian
    ASSIGNMENT_BACKEND = "builtin"


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def resolve_stem(admix_dir, stem=None):
    """Stem of the ADMIXTURE input fileset. 05 writes it to stem.txt."""
    if stem:
        return stem
    p = Path(admix_dir) / "stem.txt"
    if p.exists():
        s = p.read_text().strip()
        if s:
            return s
    raise ValueError(
        f"stem not given and {p} is missing or empty. The stem is the basename "
        f"of the .bed ADMIXTURE was run on; 05_admixture.sh writes it there.")


def load_runs(admix_dir, K, stem=None, n_reps=None):
    """Return {rep: Q} for one K. Q is (n_individuals, K), rows summing to 1."""
    admix_dir = Path(admix_dir)
    stem = resolve_stem(admix_dir, stem)
    runs = {}
    for d in sorted((admix_dir / "runs").glob(f"K{K}_r*")):
        rep = int(d.name.split("_r")[1])
        if n_reps is not None and rep > n_reps:
            continue
        q = d / f"{stem}.{K}.Q"
        if not q.exists():
            continue
        Q = np.loadtxt(q)
        if Q.ndim == 1:
            Q = Q.reshape(-1, 1)
        if Q.shape[1] != K:
            raise ValueError(f"{q}: {Q.shape[1]} columns, expected {K}")
        runs[rep] = Q
    if not runs:
        raise FileNotFoundError(
            f"no completed runs for K={K} under {admix_dir} (stem '{stem}')")
    return runs


def load_loglik(cv_raw_tsv, K):
    """{rep: loglikelihood} from cv_raw.tsv, for picking the reference run."""
    out = {}
    with open(cv_raw_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ik, ir, il = header.index("K"), header.index("REP"), header.index("LOGLIK")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ik, ir, il) or int(f[ik]) != K or f[il] == "NA":
                continue
            out[int(f[ir])] = float(f[il])
    return out


def load_individuals(path):
    return [l.strip() for l in open(path) if l.strip()]


def load_meta(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    hdr = [c.lstrip("#") for c in rows[0]]
    i = hdr.index("IID")
    return hdr, {r[i]: r for r in rows[1:]}


# ---------------------------------------------------------------------------
# alignment
# ---------------------------------------------------------------------------
def _cost(A, B):
    """Negative correlation between every column of A and every column of B."""
    K = A.shape[1]
    C = np.empty((K, K))
    for i in range(K):
        for j in range(K):
            a, b = A[:, i], B[:, j]
            if a.std() == 0 or b.std() == 0:
                C[i, j] = 0.0
            else:
                C[i, j] = -np.corrcoef(a, b)[0, 1]
    return C


def align_to(Q, ref):
    """Permute Q's columns to best match ref. Returns (Q_aligned, permutation)."""
    _, col = linear_sum_assignment(_cost(ref, Q))
    return Q[:, col], col


def similarity(A, B):
    """Mean correlation of matched columns, after aligning B to A."""
    B2, _ = align_to(B, A)
    r = []
    for i in range(A.shape[1]):
        a, b = A[:, i], B2[:, i]
        if a.std() == 0 or b.std() == 0:
            r.append(1.0 if np.allclose(a, b) else 0.0)
        else:
            r.append(np.corrcoef(a, b)[0, 1])
    return float(np.mean(r))


def align_all(runs, reference=None, loglik=None):
    """
    Align every run to a reference. The reference is the run with the highest
    loglikelihood when one is supplied, otherwise the lowest replicate number,
    so the choice is deterministic either way.

    Returns ({rep: Q_aligned}, reference_rep).
    """
    if reference is None:
        if loglik:
            cand = [r for r in runs if r in loglik]
            reference = max(cand, key=lambda r: loglik[r]) if cand else min(runs)
        else:
            reference = min(runs)
    ref = runs[reference]
    return {rep: (Q if rep == reference else align_to(Q, ref)[0])
            for rep, Q in runs.items()}, reference


# ---------------------------------------------------------------------------
# modes
# ---------------------------------------------------------------------------
class Mode:
    """A group of replicate runs that converged on the same solution.

    A plain class rather than a dataclass: dataclasses are Python 3.7+ and the
    analysis environment is 3.6.8.
    """

    def __init__(self, reps=None, consensus=None,
                 mean_similarity=float("nan"), min_similarity=float("nan")):
        self.reps = list(reps) if reps else []
        self.consensus = consensus
        self.mean_similarity = mean_similarity
        self.min_similarity = min_similarity

    @property
    def size(self):
        return len(self.reps)


def find_modes(aligned, threshold=SIM_THRESHOLD, anchor=None):
    """
    Group aligned runs into modes by single linkage on pairwise similarity.
    Returns a list of Mode, largest first. modes[0] is the major mode and
    modes[0].consensus is the average Q over its members.

    `anchor` is a replicate number. When a mode contains it, the consensus is
    built in THAT run's column order rather than the lowest replicate's. This
    matters because 05 hands the anchor run's .P file to ADMIXTURE's -P
    projection: the projected Q inherits the .P column order, so the consensus
    the figure is drawn from has to be in the same order or the training and
    projected panels would show the same components under different colours.
    """
    reps = sorted(aligned)
    n = len(reps)
    S = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            S[i, j] = S[j, i] = similarity(aligned[reps[i]], aligned[reps[j]])

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if S[i, j] >= threshold:
                a, b = find(i), find(j)
                if a != b:
                    parent[b] = a

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    modes = []
    for members in groups.values():
        # re-align every member to the anchor, so averaging is meaningful
        a_local = members[0]
        if anchor is not None:
            for m in members:
                if reps[m] == anchor:
                    a_local = m
                    break
        A0 = aligned[reps[a_local]]
        stack = [A0] + [align_to(aligned[reps[m]], A0)[0]
                        for m in members if m != a_local]
        pair = [S[a, b] for ai, a in enumerate(members) for b in members[ai + 1:]]
        modes.append(Mode(reps=[reps[m] for m in members],
                          consensus=np.mean(stack, axis=0),
                          mean_similarity=float(np.mean(pair)) if pair else 1.0,
                          min_similarity=float(np.min(pair)) if pair else 1.0))
    modes.sort(key=lambda m: (-m.size, min(m.reps)))
    return modes


def summarise(admix_dir, K, stem=None, threshold=SIM_THRESHOLD):
    """One line of text plus the major-mode Q, for use directly from Python."""
    admix_dir = Path(admix_dir)
    cv = admix_dir / "cv_raw.tsv"
    runs = load_runs(admix_dir, K, stem)
    ll = load_loglik(cv, K) if cv.exists() else None
    aligned, ref = align_all(runs, loglik=ll)
    modes = find_modes(aligned, threshold)
    major = modes[0]
    text = (f"K={K}: {len(runs)} runs, {len(modes)} mode(s); "
            f"major mode {major.size}/{len(runs)} runs "
            f"(mean similarity {major.mean_similarity:.4f}, "
            f"min {major.min_similarity:.4f}), reference run r{ref}")
    return text, major, modes


# ---------------------------------------------------------------------------
# output helpers
# ---------------------------------------------------------------------------
def write_tsv(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(map(str, header)) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"  wrote {path} ({len(rows)} rows)", file=sys.stderr)


def q_rows(Q, inds, meta_hdr, meta):
    K = Q.shape[1]
    if len(inds) != Q.shape[0]:
        raise ValueError(f"{len(inds)} individuals but Q has {Q.shape[0]} rows; "
                         f"the individuals file must match the .fam the run used")
    extra = [c for c in ("Nation", "Population", "Source", "Platform")
             if c in meta_hdr]
    header = ["IID"] + extra + [f"Q{i+1}" for i in range(K)]
    rows = []
    for i, iid in enumerate(inds):
        m = meta.get(iid)
        vals = [m[meta_hdr.index(c)] if m else "NA" for c in extra]
        rows.append([iid] + vals + [f"{v:.6f}" for v in Q[i]])
    return header, rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def cmd_summarise(a):
    admix_dir = Path(a.admix_dir)
    stem = resolve_stem(admix_dir, a.stem)
    runs = load_runs(admix_dir, a.k, stem)
    cv = admix_dir / "cv_raw.tsv"
    ll = load_loglik(cv, a.k) if cv.exists() else {}

    # First pass: group the runs into modes at all.
    aligned, ref = align_all(runs, loglik=ll or None)
    modes = find_modes(aligned, a.threshold)

    # Representative run = highest loglikelihood INSIDE the major mode, not the
    # global best. 05 hands that run's .P to ADMIXTURE's -P projection, so it
    # has to be a member of the mode whose consensus goes in the figure.
    if ll:
        cand = [r for r in modes[0].reps if r in ll]
        rep_run = max(cand, key=lambda r: ll[r]) if cand else min(modes[0].reps)
    else:
        rep_run = min(modes[0].reps)

    # Second pass, anchored on that run, so the consensus columns are in the
    # same order as the .P file the projection will use. Mode membership does
    # not depend on the reference (similarity() re-aligns internally), so this
    # re-labels rather than re-groups.
    aligned, _ = align_all(runs, reference=rep_run)
    modes = find_modes(aligned, a.threshold, anchor=rep_run)
    major = next(m for m in modes if rep_run in m.reps)

    inds = load_individuals(a.individuals)
    meta_hdr, meta = load_meta(a.meta) if a.meta else ([], {})
    hdr, rows = q_rows(major.consensus, inds, meta_hdr, meta)
    write_tsv(f"{a.out_prefix}.Q_major.tsv", hdr, rows)

    write_tsv(f"{a.out_prefix}.modes.tsv",
              ["mode", "n_runs", "reps", "mean_similarity", "min_similarity"],
              [[i + 1, m.size, ",".join(map(str, sorted(m.reps))),
                f"{m.mean_similarity:.6f}", f"{m.min_similarity:.6f}"]
               for i, m in enumerate(modes)])

    chained = "yes" if major.min_similarity < a.threshold else "no"
    write_tsv(f"{a.out_prefix}.summary.tsv", ["key", "value"], [
        ["K", a.k],
        ["n_runs", len(runs)],
        ["n_modes", len(modes)],
        ["major_mode_runs", major.size],
        ["major_mode_fraction", f"{major.size / len(runs):.3f}"],
        ["major_mode_mean_similarity", f"{major.mean_similarity:.6f}"],
        ["major_mode_min_similarity", f"{major.min_similarity:.6f}"],
        ["major_mode_single_linkage_chained", chained],
        ["alignment_reference_rep", ref],
        ["major_mode_representative_rep", rep_run],
        ["representative_P", str(admix_dir / "runs" / f"K{a.k}_r{rep_run}" /
                                 f"{stem}.{a.k}.P")],
        ["column_order", "Q_major columns are in the representative run's own "
                         "order, which is the order of representative_P"],
        ["threshold", a.threshold],
    ])
    print("  " + (f"K={a.k}: {len(runs)} runs, {len(modes)} mode(s); major mode "
                  f"{major.size}/{len(runs)} (mean sim "
                  f"{major.mean_similarity:.4f}, min {major.min_similarity:.4f}); "
                  f"representative r{rep_run}"), file=sys.stderr)
    if chained == "yes":
        print("  WARNING: the major mode is held together by single linkage; "
              "its members are not all mutually similar. Do not average it "
              "without looking at modes.tsv.", file=sys.stderr)


def cmd_sweep(a):
    lines = []
    for K in range(a.k_min, a.k_max + 1):
        try:
            runs = load_runs(a.admix_dir, K, a.stem)
        except FileNotFoundError:
            print(f"  K={K}: no runs, skipped", file=sys.stderr)
            continue
        cv = Path(a.admix_dir) / "cv_raw.tsv"
        ll = load_loglik(cv, K) if cv.exists() else {}
        aligned, ref = align_all(runs, loglik=ll or None)
        modes = find_modes(aligned, a.threshold)
        if ll:
            cand = [r for r in modes[0].reps if r in ll]
            rep_run = max(cand, key=lambda r: ll[r]) if cand else min(modes[0].reps)
        else:
            rep_run = min(modes[0].reps)
        aligned, _ = align_all(runs, reference=rep_run)
        modes = find_modes(aligned, a.threshold, anchor=rep_run)
        major = next(m for m in modes if rep_run in m.reps)
        lines.append([K, len(runs), len(modes), major.size,
                      f"{major.size / len(runs):.3f}",
                      f"{major.mean_similarity:.6f}",
                      f"{major.min_similarity:.6f}", rep_run])

        ns = argparse.Namespace(admix_dir=a.admix_dir, k=K, stem=a.stem,
                                threshold=a.threshold, individuals=a.individuals,
                                meta=a.meta, out_prefix=f"{a.out_prefix}.K{K}")
        cmd_summarise(ns)
    write_tsv(f"{a.out_prefix}.modes_by_K.tsv",
              ["K", "n_runs", "n_modes", "major_mode_runs", "major_mode_fraction",
               "mean_similarity", "min_similarity", "representative_rep"], lines)


def cmd_q_table(a):
    Q = np.loadtxt(a.q)
    if Q.ndim == 1:
        Q = Q.reshape(-1, 1)
    inds = load_individuals(a.individuals)
    meta_hdr, meta = load_meta(a.meta) if a.meta else ([], {})
    hdr, rows = q_rows(Q, inds, meta_hdr, meta)
    write_tsv(a.out, hdr, rows)


def _read_q_table(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    hdr = rows[0]
    i_id = hdr.index("IID")
    qcols = [j for j, c in enumerate(hdr) if c.startswith("Q") and c[1:].isdigit()]
    return {r[i_id]: np.array([float(r[j]) for j in qcols]) for r in rows[1:]}, len(qcols)


def cmd_compare_q(a):
    """Correlate two Q matrices after optimal column matching.

    Used for the projected-versus-joint ADMIXTURE comparison: if capping the
    cohort and projecting reproduces the joint fit, the matched components
    correlate near 1 and the capping objection is answered with a number.
    """
    A, ka = _read_q_table(a.a)
    B, kb = _read_q_table(a.b)
    if ka != kb:
        raise SystemExit(f"K differs: {ka} vs {kb}")
    shared = sorted(set(A) & set(B))
    if len(shared) < 3:
        raise SystemExit(f"only {len(shared)} shared individuals")
    QA = np.array([A[i] for i in shared])
    QB = np.array([B[i] for i in shared])
    QB2, perm = align_to(QB, QA)
    rows = []
    for k in range(ka):
        x, y = QA[:, k], QB2[:, k]
        r = float(np.corrcoef(x, y)[0, 1]) if x.std() and y.std() else float("nan")
        rows.append([f"Q{k+1}", int(perm[k]) + 1, len(shared), f"{r:.6f}",
                     f"{float(np.mean(np.abs(x - y))):.6f}"])
    write_tsv(a.out, ["component_A", "matched_component_B", "n_shared",
                      "pearson_r", "mean_abs_diff"], rows)
    rs = [float(r[3]) for r in rows]
    print(f"  mean |r| over {ka} components: {np.nanmean(np.abs(rs)):.4f}, "
          f"min {np.nanmin(np.abs(rs)):.4f}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd")
    # add_subparsers(required=True) is Python 3.7+; setting the attribute is the
    # 3.6-compatible form and behaves identically.
    sub.required = True

    def common(s):
        s.add_argument("--admix-dir", required=True)
        s.add_argument("--stem", default=None)
        s.add_argument("--threshold", type=float, default=SIM_THRESHOLD)
        s.add_argument("--individuals", required=True)
        s.add_argument("--meta", default=None)
        s.add_argument("--out-prefix", required=True)

    s = sub.add_parser("summarise"); common(s)
    s.add_argument("--k", type=int, required=True)
    s.set_defaults(fn=cmd_summarise)

    s = sub.add_parser("sweep"); common(s)
    s.add_argument("--k-min", type=int, required=True)
    s.add_argument("--k-max", type=int, required=True)
    s.set_defaults(fn=cmd_sweep)

    s = sub.add_parser("q-table")
    s.add_argument("--q", required=True)
    s.add_argument("--individuals", required=True)
    s.add_argument("--meta", default=None)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_q_table)

    s = sub.add_parser("compare-q")
    s.add_argument("--a", required=True)
    s.add_argument("--b", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_compare_q)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
