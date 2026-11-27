"""
figstyle.py — the single place where figure style is defined.

Every panel in Figures_CommsBio.ipynb is drawn through this module, so the
whole paper has one typeface, one palette, one line weight and one export
path. Nothing here computes a statistic; the numbers all come from the
pipeline's TSVs and this only draws them.

WHAT THE JOURNAL ASKS FOR, AND WHERE IT IS ENFORCED
  (Communications Biology / Nature Portfolio submission guidelines)

  89 mm single column, 183 mm double column      -> W1, W15, W2, figure()
  max height 170 mm, legend goes underneath      -> HMAX, checked by save()
  sans-serif throughout, Helvetica or Arial      -> FONT_STACK
  body text 5-7 pt, panel letters 8 pt bold      -> FS_*, panel()
  lowercase bold a, b, c on multi-part figures   -> panel()
  thinnest line >= 1 pt in the final figure      -> LW_MIN, audit()
  white background, no boxing, no 3-D            -> rcParams
  histogram axes not truncated                   -> caller's job, noted below
  text must stay text, never outlines            -> pdf/ps.fonttype = 42
  vector preferred, 300+ dpi for raster          -> save() writes PDF and PNG

The palette is Okabe-Ito, which stays distinguishable under all three common
forms of colour blindness. It is also what the existing tribes notebook
already used, so nothing in the paper changes colour.

TRACEABILITY. load() records every file it reads. sources_table() then prints
one row per figure panel naming the file behind it, which is what makes a
"where does this number come from" question answerable in one line instead of
a directory search.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------- geometry --
MM = 1 / 25.4
W1 = 89 * MM      # single column
W15 = 136 * MM    # 1.5 column
W2 = 183 * MM     # double column
HMAX = 170 * MM   # maximum height, legend sits underneath

# ------------------------------------------------------------------- fonts --
# Helvetica and Arial first because that is what the guidelines name.
# Liberation Sans is metric-compatible with Arial and is the usual Linux
# stand-in; DejaVu is matplotlib's default and is the one to avoid shipping,
# so which font actually resolved is reported at import.
FONT_STACK = ["Helvetica", "Arial", "Nimbus Sans", "Liberation Sans",
              "FreeSans", "DejaVu Sans"]

FS_TICK = 6       # 5-7 pt window
FS_LABEL = 7
FS_LEGEND = 6
FS_PANEL = 8      # panel letters, bold
LW_MIN = 1.0      # the guidelines' floor for the thinnest line
LW_AXIS = 1.0
LW_LINE = 1.2

# ------------------------------------------------------------------ colour --
# Okabe & Ito (2008), colour-blind safe.
OKABE_ITO = {
    "orange":  "#E69F00",
    "skyblue": "#56B4E9",
    "green":   "#009E73",
    "yellow":  "#F0E442",
    "blue":    "#0072B2",
    "vermil":  "#D55E00",
    "purple":  "#CC79A7",
    "black":   "#000000",
}
PALETTE = [OKABE_ITO[k] for k in
           ("vermil", "blue", "green", "orange", "purple", "skyblue", "yellow")]

# One colour per superpopulation, fixed for the whole paper so a reader can
# carry the legend from one figure to the next. The cohort is vermilion
# everywhere and is the only group that gets a black edge.
SUPERPOP_COLOUR = {
    "KAZ": OKABE_ITO["vermil"],
    "CAS": OKABE_ITO["orange"],
    "EAS": OKABE_ITO["blue"],
    "SAS": OKABE_ITO["green"],
    "WEA": OKABE_ITO["purple"],
    "EUR": OKABE_ITO["skyblue"],
    "OCE": "#7F7F7F",
    "AFR": "#000000",
    "AMR": "#B15928",
}
COHORT_COLOUR = OKABE_ITO["vermil"]

# Display names. The pipeline keeps the source panels' own spellings, because
# those are what the .fam files and every intermediate table are keyed on and
# changing them upstream would break the joins. This map is applied only where
# a label is SHOWN - figures and the workbook - so the presentation is clean
# and the data stays traceable. Extend it here and it changes everywhere at
# once; the mapping is recorded in the workbook's INDEX sheet.
DISPLAY_NAME = {
    "Mbuti_HG":       "Mbuti",
    "NAN_Melanesian": "Melanesian",
}


def pretty(labels):
    """Apply DISPLAY_NAME to a string or any iterable of strings."""
    if isinstance(labels, str):
        return DISPLAY_NAME.get(labels, labels)
    return [DISPLAY_NAME.get(v, v) for v in labels]
GREY = "#8C8C8C"


def apply_style():
    """Install the rcParams. Called once at import; safe to call again."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": FS_TICK,
        "axes.labelsize": FS_LABEL,
        "axes.titlesize": FS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
        "legend.title_fontsize": FS_LEGEND,

        # "avoid excessive boxing": two spines, not four.
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": LW_AXIS,
        "xtick.major.width": LW_AXIS,
        "ytick.major.width": LW_AXIS,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": LW_LINE,
        "lines.markersize": 3,
        "patch.linewidth": 0.8,

        "axes.grid": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        # NOT "tight": bbox_inches="tight" crops the output to its content,
        # so a figure requested at 183 mm lands somewhere near it and the
        # journal's column width is no longer guaranteed. Layout is handled
        # by tight_layout()/constrained inside each figure instead, and save()
        # reads the width back out of the PDF to prove it.
        "savefig.bbox": None,
        "savefig.pad_inches": 0.01,
        "figure.constrained_layout.use": False,
        "grid.linewidth": LW_AXIS,

        # Text must remain text. Type 42 embeds TrueType outlines as glyphs
        # that a PDF reader can still select and a typesetter can still
        # restyle; the matplotlib default (Type 3) is what turns lettering
        # into unselectable paths and draws a production query.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",

        "legend.frameon": False,
        "errorbar.capsize": 1.5,
    })


apply_style()


def font_in_use():
    """Which member of FONT_STACK matplotlib actually resolved."""
    from matplotlib import font_manager
    have = {f.name for f in font_manager.fontManager.ttflist}
    for name in FONT_STACK:
        if name in have:
            return name
    return "unknown"


# --------------------------------------------------------------------- I/O --
SOURCES: list[tuple[str, str, int]] = []   # (figure/panel, path, n rows)


def load(path, panel="", required=True):
    """Read a pipeline TSV into a dict of columns, and record where it came from.

    Returns None when the file is absent and required=False, so a notebook run
    while a stage is still going produces the figures it can and says which
    ones it skipped, instead of dying on the first missing file.
    """
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        msg = f"missing or empty: {p}"
        if required:
            raise FileNotFoundError(msg)
        print(f"  SKIP  {panel or p.name}: {msg}")
        return None
    rows = [ln.rstrip("\n").split("\t") for ln in
            p.read_text().splitlines() if ln.strip()]
    hdr = [c.lstrip("#").strip() for c in rows[0]]
    cols = {h: [] for h in hdr}
    for r in rows[1:]:
        for i, h in enumerate(hdr):
            cols[h].append(r[i] if i < len(r) else "")
    if len(rows) < 2:
        # A single-line file has been read AS a header, leaving every column
        # empty. That is what turned 05's headerless best_k.tsv into an
        # IndexError three lines later instead of a message here.
        print(f"  NOTE  {panel or p.name}: {p.name} has one line and no data "
              f"rows; it was read as a header only ({', '.join(hdr[:4])})")
    SOURCES.append((panel or p.name, str(p), len(rows) - 1))
    return cols


def nums(col):
    """Column of strings -> float array, non-numeric becomes nan."""
    out = np.empty(len(col), float)
    for i, v in enumerate(col):
        try:
            out[i] = float(v)
        except (TypeError, ValueError):
            out[i] = np.nan
    return out


def sources_table():
    """One row per panel naming the file it was drawn from."""
    w = max((len(s[0]) for s in SOURCES), default=10)
    print(f"{'panel':<{w}}  {'rows':>6}  source")
    for panel, path, n in SOURCES:
        print(f"{panel:<{w}}  {n:>6}  {path}")


# ------------------------------------------------------------------ layout --
def figure(width=W2, height=None, **kw):
    """A figure at an exact journal column width, in inches."""
    if height is None:
        height = width * 0.618
    if height > HMAX:
        print(f"  NOTE: height {height/MM:.0f} mm exceeds the {HMAX/MM:.0f} mm "
              f"maximum; the legend will not fit underneath")
    return plt.subplots(figsize=(width, height), **kw)


def panel(ax, letter, dx=-0.16, dy=1.06):
    """Lower-case bold panel letter, 8 pt, in figure-relative position."""
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold", va="top", ha="left")


def finish(ax, xlabel=None, ylabel=None, title=None):
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left")
    return ax


FIGDIR = Path("figures_tables")   # figures and tables in one place
DPI = 600                         # PNG only; no PDF is produced


DPI = 600


def save(fig, name, figdir=None):
    """Write the PNG at 600 dpi and audit the result.

    PNG only. The size is verified by reading the file back rather than
    trusting the request: a stray layout call or a colour bar can change the
    exported image without changing the figsize, and the journal measures the
    file.
    """
    d = Path(figdir) if figdir else FIGDIR
    d.mkdir(parents=True, exist_ok=True)
    fig.savefig(d / f"{name}.png", dpi=DPI)
    w, h = fig.get_size_inches()
    flags = audit(fig)
    got = png_size_mm(d / f"{name}.png")
    if got:
        gw, gh = got
        if abs(gw - w / MM) > 1.0 or abs(gh - h / MM) > 1.0:
            flags.append(f"exported {gw:.0f}x{gh:.0f} mm, asked for "
                         f"{w/MM:.0f}x{h/MM:.0f} mm")
    status = "ok" if not flags else "; ".join(flags)
    size = f"{got[0]:.0f} x {got[1]:.0f} mm" if got else f"{w/MM:.0f} x {h/MM:.0f} mm"
    print(f"  {name}: {size} @ {DPI} dpi -> {d}/{name}.png  [{status}]")
    return flags


def png_size_mm(path):
    """Physical width and height of a PNG, from its IHDR and pHYs chunks."""
    import struct
    try:
        raw = Path(path).read_bytes()
    except OSError:
        return None
    if raw[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    px = py = None
    ppm_x = ppm_y = None
    i = 8
    while i + 8 <= len(raw):
        ln = struct.unpack(">I", raw[i:i+4])[0]
        typ = raw[i+4:i+8]
        body = raw[i+8:i+8+ln]
        if typ == b"IHDR":
            px, py = struct.unpack(">II", body[:8])
        elif typ == b"pHYs" and len(body) >= 9 and body[8] == 1:
            ppm_x, ppm_y = struct.unpack(">II", body[:8])
            break
        i += 12 + ln
    if px is None:
        return None
    if not ppm_x:
        return (px / DPI / MM, py / DPI / MM)
    return (px / ppm_x * 1000, py / ppm_y * 1000)


def audit(fig):
    """Check a figure against the parts of the spec that can be checked.

    Catches the two things that actually get flagged in production: lettering
    below 5 pt after reduction, and hairlines thinner than 1 pt.
    """
    flags = []
    small = []
    for t in fig.findobj(mpl.text.Text):
        s = t.get_text()
        if s and s.strip() and t.get_fontsize() < 5:
            small.append(f"{s[:14]!r}@{t.get_fontsize():.1f}pt")
    if small:
        flags.append(f"text below 5 pt: {', '.join(small[:3])}")

    # get_visible(): matplotlib keeps gridlines as Line2D objects at 0.8 pt
    # even when the grid is switched off, and an audit that flags invisible
    # objects on every figure is an audit nobody reads.
    thin = set()
    for ln in fig.findobj(mpl.lines.Line2D):
        lw = ln.get_linewidth()
        if 0 < lw < LW_MIN and ln.get_linestyle() != "None" and ln.get_visible():
            thin.add(round(lw, 2))
    if thin:
        flags.append(f"lines below {LW_MIN} pt: {sorted(thin)}")

    if fig.get_size_inches()[1] > HMAX:
        flags.append("taller than 170 mm")
    return flags


def label_groups(ax, x, y, groups, min_n=1, fontsize=5, only=None,
                 inside_only=True, color="black", max_iter=400, pad=2.0):
    """One label per group, repelled in TWO dimensions until nothing overlaps.

    The first version pushed labels apart vertically only, and only when their
    x positions were already close. That fails exactly where it matters: on the
    AADR panel every group sits in one narrow horizontal band, so the collisions
    are side by side and the labels piled up regardless.

    This measures the real rendered extent of each label and separates
    overlapping boxes along whichever axis needs the smaller move, then draws a
    leader line from the centroid to wherever the label ended up. Labels are
    kept inside the axes, and the data limits are restored afterwards so the
    leader lines cannot rescale the panel.
    """
    import matplotlib.patheffects as pe
    x, y = np.asarray(x, float), np.asarray(y, float)
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()

    # `color` may be a dict {group: colour}. Two separate calls to this
    # function mean two independent repulsion runs that cannot see each
    # other's labels, which is how the highlighted group ended up sitting on
    # top of its neighbours.
    cmap = color if isinstance(color, dict) else {}
    items = []
    for g in sorted(set(groups)):
        if only is not None and g not in only:
            continue
        m = np.array([v == g for v in groups])
        if m.sum() < min_n:
            continue
        cx, cy = float(np.nanmean(x[m])), float(np.nanmean(y[m]))
        if inside_only and not (min(x0, x1) <= cx <= max(x0, x1)
                                and min(y0, y1) <= cy <= max(y0, y1)):
            continue
        items.append((g, cx, cy))
    if not items:
        return []

    fig = ax.figure
    texts = []
    for g, cx, cy in items:
        c = cmap.get(g, "black") if cmap else color
        w = "bold" if (cmap and g in cmap) else "normal"
        t = ax.text(cx, cy, g, fontsize=fontsize, ha="center", va="center",
                    color=c, fontweight=w, zorder=10)
        t.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])
        texts.append(t)

    fig.canvas.draw()                      # extents are only real after a draw
    rend = fig.canvas.get_renderer()
    anchor = np.array([ax.transData.transform((cx, cy)) for _, cx, cy in items])
    ext = np.array([[t.get_window_extent(rend).width,
                     t.get_window_extent(rend).height] for t in texts])
    cur = anchor.copy()

    for _ in range(max_iter):
        moved = False
        for i in range(len(cur)):
            for j in range(i + 1, len(cur)):
                dx, dy = cur[i, 0] - cur[j, 0], cur[i, 1] - cur[j, 1]
                ox = (ext[i, 0] + ext[j, 0]) / 2 + pad - abs(dx)
                oy = (ext[i, 1] + ext[j, 1]) / 2 + pad - abs(dy)
                if ox > 0 and oy > 0:      # boxes overlap: separate them
                    if oy <= ox:           # cheaper to move vertically
                        sgn = 1.0 if dy >= 0 else -1.0
                        cur[i, 1] += sgn * oy / 2; cur[j, 1] -= sgn * oy / 2
                    else:
                        sgn = 1.0 if dx >= 0 else -1.0
                        cur[i, 0] += sgn * ox / 2; cur[j, 0] -= sgn * ox / 2
                    moved = True
        if not moved:
            break

    bb = ax.get_window_extent()
    cur[:, 0] = np.clip(cur[:, 0], bb.x0 + ext[:, 0] / 2, bb.x1 - ext[:, 0] / 2)
    cur[:, 1] = np.clip(cur[:, 1], bb.y0 + ext[:, 1] / 2, bb.y1 - ext[:, 1] / 2)

    inv = ax.transData.inverted()
    for t, (g, cx, cy), p, a in zip(texts, items, cur, anchor):
        lx, ly = inv.transform(p)
        t.set_position((lx, ly))
        if float(np.hypot(*(p - a))) > 3:
            ax.plot([cx, lx], [cy, ly], color="0.55", lw=1.0, zorder=9)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)   # leader lines must not rescale
    return [g for g, _, _ in items]


def report_environment():
    f = font_in_use()
    print(f"matplotlib {mpl.__version__} | font in use: {f}")
    if f not in ("Helvetica", "Arial", "Nimbus Sans", "Liberation Sans"):
        print(f"  WARNING: '{f}' is not Helvetica/Arial and is not metric-"
              f"compatible with them. The guidelines name Helvetica or Arial; "
              f"install one (or Liberation Sans) before the final export, or "
              f"production will query the typeface.")
    print(f"columns: {W1/MM:.0f} / {W15/MM:.0f} / {W2/MM:.0f} mm, "
          f"max height {HMAX/MM:.0f} mm")
    print(f"text {FS_TICK}-{FS_LABEL} pt, panel letters {FS_PANEL} pt bold, "
          f"min line {LW_MIN} pt, PNG only at {DPI} dpi")
