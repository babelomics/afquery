#!/usr/bin/env python3
"""Generate the Gap 2 figure: mixed sequencing technologies inflating the denominator."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---------------------------------------------------------------------------
# Data: schematic gene with 8 exons
# ---------------------------------------------------------------------------
exons = [
    (5, 12),    # Exon 1
    (18, 24),   # Exon 2
    (30, 38),   # Exon 3
    (44, 50),   # Exon 4
    (56, 63),   # Exon 5
    (69, 76),   # Exon 6
    (82, 88),   # Exon 7
    (93, 98),   # Exon 8
]

# Technologies: name, sample count, covered exon indices, color
# WES_KIT_V1 and WES_KIT_V2 are versions of the same kit; V2 improves
# coverage by adding exon 7.
technologies = [
    ("WES_KIT_V1", "(n = 120)", [0, 1, 2, 3, 4, 5],    "#2196F3"),  # exons 1-6
    ("WES_KIT_V2", "(n = 80)",  [0, 1, 2, 3, 4, 5, 6], "#FF9800"),  # exons 1-7
]

n_samples = {"WES_KIT_V1": 120, "WES_KIT_V2": 80}

# ---------------------------------------------------------------------------
# Two query positions
# ---------------------------------------------------------------------------
query_positions = [
    # Label  x-pos  exon_idx
    ("A", 60, 4),   # Exon 5: covered by WGS, V1, V2 (not Panel)
    ("B", 85, 6),   # Exon 7: covered by WGS, V2 only (not V1 nor Panel)
]

query_colors = {
    "A": "#00838F",  # teal
    "B": "#C62828",  # red
}

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10.9, 4.5))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

y_positions = []
bar_height = 0.55
row_spacing = 1.2
n_tech = len(technologies)

# -- Gene backbone -----------------------------------------------------------
gene_y = n_tech * row_spacing + 0.8
ax.plot([0, 100], [gene_y, gene_y], color="#BDBDBD", linewidth=1.5, zorder=1)
for i, (start, end) in enumerate(exons):
    rect = FancyBboxPatch(
        (start, gene_y - 0.20), end - start, 0.40,
        boxstyle="round,pad=0.05", facecolor="#E0E0E0", edgecolor="#9E9E9E",
        linewidth=0.8, zorder=2,
    )
    ax.add_patch(rect)
    ax.text(
        (start + end) / 2, gene_y,
        f"E{i+1}", ha="center", va="center", fontsize=8, color="#616161",
        fontweight="bold", zorder=3,
    )

ax.text(-3, gene_y, "Gene X", ha="right", va="center", fontsize=9,
        fontstyle="italic", color="#616161")

# -- Technology rows ---------------------------------------------------------
for idx, (name, n_label, covered, color) in enumerate(technologies):
    y = (n_tech - 1 - idx) * row_spacing
    y_positions.append(y)

    ax.text(-3, y + 0.08, name, ha="right", va="center", fontsize=9,
            fontweight="bold", color="#333333")
    ax.text(-3, y - 0.28, n_label, ha="right", va="center", fontsize=8.5,
            color="#777777")

    for ei in covered:
        start, end = exons[ei]
        rect = FancyBboxPatch(
            (start, y - bar_height / 2), end - start, bar_height,
            boxstyle="round,pad=0.08", facecolor=color, edgecolor="white",
            alpha=0.75, linewidth=0.5, zorder=2,
        )
        ax.add_patch(rect)

    if covered:
        x_min = exons[covered[0]][0]
        x_max = exons[covered[-1]][1]
        ax.plot([x_min, x_max], [y, y], color=color, linewidth=0.8,
                alpha=0.3, zorder=1)

# -- Query position lines and indicators ------------------------------------
for label, qx, exon_idx in query_positions:
    qcolor = query_colors[label]

    ax.axvline(x=qx, color=qcolor, linewidth=1.4, linestyle="--",
               alpha=0.55, zorder=4)

    # Label above gene
    ax.text(qx, gene_y + 0.65, label, ha="center", va="bottom",
            fontsize=10, fontweight="bold", color=qcolor,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                      edgecolor=qcolor, linewidth=1.0, alpha=0.9),
            zorder=8)

    for idx, (name, _nlabel, covered, _tcolor) in enumerate(technologies):
        y = y_positions[idx]
        if exon_idx in covered:
            ax.plot(qx, y, marker="o", markersize=9, color="#1B5E20",
                    zorder=7, markeredgecolor="white", markeredgewidth=1.0)
            ax.text(qx, y, "\u2713", ha="center", va="center",
                    fontsize=7.5, fontweight="bold", color="white", zorder=8)
        else:
            ax.plot(qx, y, marker="o", markersize=9, color="#D32F2F",
                    zorder=7, markeredgecolor="white", markeredgewidth=1.0)
            ax.text(qx, y, "\u2717", ha="center", va="center",
                    fontsize=7.5, fontweight="bold", color="white", zorder=8)

# -- Footnote ----------------------------------------------------------------
ax.text(
    0.50, -0.06,
    "\u26a0  At position B the correct AN is 160 but na\u00efve counting gives 400 "
    "\u2014 AF is underestimated by more than half",
    transform=ax.transAxes, fontsize=7.5, ha="center", va="top",
    color="#C62828", fontstyle="italic",
)

# -- Axes and layout ---------------------------------------------------------
ax.set_xlim(-18, 105)
ax.set_ylim(-0.8, gene_y + 1.8)
ax.axis("off")

ax.set_title(
    "Mixed Sequencing Technologies: Coverage at Different Query Positions",
    fontsize=12, fontweight="bold", color="#333333", pad=14,
)

plt.tight_layout()

out_path = (Path(__file__).resolve().parent.parent
            / "assets" / "img" / "gap2_mixed_technologies.png")
out_path.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {out_path}")

plt.close()
