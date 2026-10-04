"""
Figure: where the human/mouse differences sit relative to the cetuximab
epitope, across EGFR domain III.

Produces figures/domain3_epitope_map.png

Usage:  python analysis/05_make_figure.py
"""

from pathlib import Path
from math import dist

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
PDB = ROOT / "data" / "structures" / "6ARU.pdb"
SEQ = ROOT / "data" / "sequences"
OUT = ROOT / "figures" / "domain3_epitope_map.png"

CUTOFF = 5.0
MATURE_START_MOTIF = "LEEKKVCQGTSN"
D3_FIRST, D3_LAST = 310, 481

INK = "#1c1c1c"
GREY = "#c9c9c9"
CONSERVED = "#2f6f9f"
DIVERGENT = "#c4442e"


def load_mature(path):
    seq = "".join(path.read_text().split()).upper()
    i = seq.find(MATURE_START_MOTIF)
    if i < 0:
        raise SystemExit(f"{path.name}: mature start motif not found.")
    return seq[i:]


def main():
    if not PDB.exists():
        raise SystemExit(
            f"Missing {PDB}\n  curl -o data/structures/6ARU.pdb "
            "https://files.rcsb.org/download/6ARU.pdb"
        )

    human = load_mature(SEQ / "human_EGFR_mature.txt")
    mouse = load_mature(SEQ / "mouse_EGFR_precursor.txt")
    h3 = human[D3_FIRST - 1:D3_LAST]
    m3 = mouse[D3_FIRST - 1:D3_LAST]
    variants = {
        p: (h, m)
        for p, (h, m) in enumerate(zip(h3, m3), start=D3_FIRST) if h != m
    }

    atoms = []
    for line in PDB.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[76:78].strip().upper() in {"H", "D"}:
            continue
        atoms.append((line[21].strip(), line[22:26].strip(),
                      (float(line[30:38]), float(line[38:46]),
                       float(line[46:54]))))

    egfr = [a for a in atoms if a[0] == "A"]
    fab = [a for a in atoms if a[0] in ("B", "C")]

    nearest = {}
    for c, rn, xyz in egfr:
        pos = int(rn)
        best = nearest.get(pos, float("inf"))
        for _, _, bxyz in fab:
            d = dist(xyz, bxyz)
            if d < best:
                best = d
        nearest[pos] = best

    epitope = {p for p, d in nearest.items() if d < CUTOFF}

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(13, 5.8), height_ratios=[1.15, 1.5],
        gridspec_kw={"hspace": 0.55},
    )

    # --- Track 1: the domain III ladder -------------------------------
    positions = list(range(D3_FIRST, D3_LAST + 1))
    ax.hlines(0, D3_FIRST - 1, D3_LAST + 1, color=GREY, lw=7,
              zorder=1, capstyle="round")

    ep = sorted(epitope)
    ax.scatter(ep, [0] * len(ep), s=46, color=CONSERVED, zorder=3,
               edgecolor="white", linewidth=0.7)

    both = sorted(epitope & set(variants))
    ax.scatter(both, [0] * len(both), s=46, color=DIVERGENT, zorder=4,
               edgecolor="white", linewidth=0.7)

    other_var = sorted(set(variants) - epitope)
    ax.scatter(other_var, [-0.33] * len(other_var), s=26, color=DIVERGENT,
               alpha=0.45, zorder=3, marker="v")

    # Label divergent epitope residues. Positions closer than MIN_GAP apart
    # would overlap, so they are stacked on alternating rows instead.
    MIN_GAP = 14
    rows, last = [], {}
    for p in both:
        row = 0
        while last.get(row) is not None and p - last[row] < MIN_GAP:
            row += 1
        last[row] = p
        rows.append(row)

    for p, row in zip(both, rows):
        h, m = variants[p]
        dy = 13 + row * 11
        ax.annotate(f"{h}{p}{m}", (p, 0), xytext=(0, dy),
                    textcoords="offset points", ha="center", fontsize=8.5,
                    color=DIVERGENT, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", lw=0.6,
                                    color=DIVERGENT, alpha=0.5,
                                    shrinkA=0, shrinkB=3))

    ax.set_ylim(-0.75, 1.5)
    ax.set_xlim(D3_FIRST - 3, D3_LAST + 3)
    ax.set_yticks([])
    ax.set_xlabel("EGFR domain III — mature numbering", fontsize=9.5,
                  color=INK)
    ax.tick_params(axis="x", labelsize=8.5, colors=INK)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GREY)
    ax.set_title(
        "Seven of the 21 human/mouse differences sit inside the "
        "cetuximab epitope",
        fontsize=12, color=INK, pad=16, loc="left", fontweight="bold",
    )

    ax.legend(handles=[
        Patch(facecolor=CONSERVED, label="epitope residue, conserved"),
        Patch(facecolor=DIVERGENT, label="epitope residue, differs in mouse"),
        Patch(facecolor=DIVERGENT, alpha=0.45,
              label="difference outside the epitope"),
    ], fontsize=8.5, frameon=False, loc="lower left",
        bbox_to_anchor=(0, -0.48), ncol=3)

    # --- Track 2: distance to the Fab ---------------------------------
    xs = [p for p in positions if p in nearest]
    ys = [min(nearest[p], 30) for p in xs]
    ax2.fill_between(xs, ys, 30, color=GREY, alpha=0.3, zorder=1)
    ax2.plot(xs, ys, color=INK, lw=1.1, zorder=2)
    ax2.axhline(CUTOFF, color=CONSERVED, lw=1.1, ls="--", zorder=3)
    ax2.annotate(f"{CUTOFF:.0f} Å contact cutoff", (D3_FIRST + 2, CUTOFF),
                 xytext=(0, -13), textcoords="offset points", ha="left",
                 fontsize=8, color=CONSERVED,
                 bbox=dict(boxstyle="round,pad=0.25", fc="white",
                           ec="none", alpha=0.85))

    vx = [p for p in sorted(variants) if p in nearest]
    vy = [min(nearest[p], 30) for p in vx]
    ax2.scatter(vx, vy, s=30, color=DIVERGENT, zorder=4,
                edgecolor="white", linewidth=0.6)

    ax2.set_ylim(0, 30)
    ax2.set_xlim(D3_FIRST - 3, D3_LAST + 3)
    ax2.invert_yaxis()
    ax2.set_ylabel("min. heavy-atom distance\nto cetuximab Fab (Å)",
                   fontsize=9, color=INK)
    ax2.set_xlabel("EGFR domain III — mature numbering", fontsize=9.5,
                   color=INK)
    ax2.tick_params(labelsize=8.5, colors=INK)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax2.spines[s].set_color(GREY)

    fig.text(0.008, 0.015,
             "Structure: PDB 6ARU (cetuximab Fab mutant + EGFR ECD).  "
             "Sequences: UniProt P00533, Q01279.  "
             "Calibration only — no part of this antibody is used in any "
             "design.",
             fontsize=7.4, color="#666666")

    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"wrote {OUT}")
    print(f"  epitope residues: {len(epitope)}")
    print(f"  of which differ in mouse: {len(both)}")


if __name__ == "__main__":
    main()
