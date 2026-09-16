"""Regenerate analytical figures from the archived result tables.

Figures 1 and 3-7 are regenerated into ``figures/``. Figure 2 is retained as a
static contextual composite because its municipal imagery is not a computational
input to the analysis; see ``figures/README.md`` and ``THIRD_PARTY_DATA.md``.
"""
from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
RES = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

blue = "#6FA8C4"
green = "#91B98D"
yellow = "#E8C96A"
orange = "#D9966C"
ink = "#273746"
purple = "#B7A6C8"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
    }
)


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


# Figure 1: study design
fig, ax = plt.subplots(figsize=(10.5, 4.55))
ax.axis("off")
ax.set_xlim(0, 10)
ax.set_ylim(0, 5)
ax.text(5, 4.65, "Study design", ha="center", va="center", fontsize=13, fontweight="bold", color=ink)
boxes = [
    (0.3, 3.05, 1.7, 1.0, blue, "LoD2 roof surfaces\nand ALKIS check"),
    (2.35, 3.05, 1.8, 1.0, green, "Geometry audit and\ncoplanar consolidation"),
    (4.55, 3.05, 1.65, 1.0, yellow, "Complete-module\npacking"),
    (6.6, 3.05, 1.55, 1.0, orange, "Capacity\ncomparison"),
]
for x, y, w, h, c, text in boxes:
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.03", fc=c, ec="#666666", lw=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)
for x1, x2 in [(2.0, 2.35), (4.15, 4.55), (6.2, 6.6)]:
    ax.add_patch(FancyArrowPatch((x1, 3.55), (x2, 3.55), arrowstyle="-|>", mutation_scale=12, lw=1, color="#555555"))
ax.add_patch(FancyBboxPatch((2.0, 0.8), 2.35, 1.0, boxstyle="round,pad=.02", fc="#EEF3F5", ec="#777777", lw=1))
ax.text(3.175, 1.30, "Controlled experiment\nfixed outer roof; increasing subdivision", ha="center", va="center", fontsize=9)
ax.add_patch(FancyBboxPatch((5.3, 0.8), 2.35, 1.0, boxstyle="round,pad=.02", fc="#F4F0E8", ec="#777777", lw=1))
ax.text(6.475, 1.30, "Stuttgart pilot\nobserved LoD2 roof geometry", ha="center", va="center", fontsize=9)
ax.add_patch(FancyArrowPatch((3.25, 3.05), (3.25, 1.82), arrowstyle="-|>", mutation_scale=11, lw=1, color="#555555"))
ax.add_patch(FancyArrowPatch((5.35, 3.05), (6.25, 1.82), arrowstyle="-|>", mutation_scale=11, lw=1, color="#555555"))
ax.text(5, 0.25, "The controlled experiment isolates subdivision; the Stuttgart pilot measures the area-to-module packing shortfall.", ha="center", fontsize=8.6, color="#555555")
save(fig, "Fig1_framework")

# Figure 2 is intentionally retained; fail clearly if it has been removed.
for ext in ("png", "svg"):
    p = OUT / f"Fig2_data_context.{ext}"
    if not p.exists():
        raise FileNotFoundError(f"Missing static contextual figure: {p}")

# Figure 3: controlled experiment
cs = pd.read_csv(RES / "controlled_production_summary.csv")
fig, ax = plt.subplots(figsize=(8.7, 5.1))
x = cs.facet_count.to_numpy()
med = 100 * cs.median_loss.to_numpy()
q1 = 100 * cs.q25.to_numpy()
q3 = 100 * cs.q75.to_numpy()
ax.fill_between(x, q1, q3, color=blue, alpha=.24, label="Interquartile range")
ax.plot(x, med, "-o", color=ink, lw=1.8, ms=4, label="Median incremental loss")
ax.set_xlabel("Number of non-crossable roof facets")
ax.set_ylabel("Incremental packing loss (%)")
ax.set_xticks(x)
ax.set_ylim(0, 55)
ax.grid(True, alpha=.18)
ax.legend(frameon=False, loc="upper left")
ax.set_title("Subdivision reduces feasible module packing when outer roof geometry is fixed", fontweight="bold")
save(fig, "Fig3_controlled_fragmentation")

# Figure 4: production vs denser numerical search
b = pd.read_csv(RES / "packing_benchmark.csv")
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.25), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
ax.scatter(b.stress_count, b.production_count, c=b.usable_area_m2, cmap="viridis", s=48, edgecolor="white", linewidth=.5)
mx = max(b.stress_count.max(), b.production_count.max())
ax.plot([0, mx], [0, mx], "--", color="#888888", lw=1)
ax.set_xlabel("Denser-search module count")
ax.set_ylabel("Production-search module count")
ax.set_title("(a) Deterministic benchmark")
ax.grid(True, alpha=.15)
ax = axs[1]
ax.bar(np.arange(1, len(b) + 1), b.gap_modules, color=orange, edgecolor="#8b6a55", linewidth=.4)
ax.set_xlabel("Benchmark facet")
ax.set_ylabel("Additional modules found")
ax.set_xticks(np.arange(1, len(b) + 1))
ax.set_title("(b) Residual search difference")
ax.grid(axis="y", alpha=.15)
fig.suptitle("Packing-search check across the observed facet-size range", fontweight="bold", y=1.02)
fig.tight_layout()
save(fig, "Fig4_packing_stress_benchmark")

# Figure 5: aggregate and building distribution
base = pd.read_csv(RES / "base_buildings.csv")
cont = base.continuous_kwp.sum() / 1000
pack = base.packed_kwp.sum() / 1000
fig, axs = plt.subplots(1, 2, figsize=(10.7, 4.55))
ax = axs[0]
bars = ax.bar(["Continuous-area\nreference", "Explicit module\npacking"], [cont, pack], color=[yellow, blue], edgecolor="#666666", linewidth=.6, width=.58)
ax.set_ylabel("DC capacity (MWp)")
ax.set_title("(a) Aggregate Stuttgart pilot")
ax.grid(axis="y", alpha=.15)
for bar, value in zip(bars, [cont, pack]):
    ax.text(bar.get_x() + bar.get_width() / 2, value + .15, f"{value:.2f} MWp", ha="center", fontsize=9)
ax = axs[1]
y = base.packing_shortfall_pct.to_numpy()
ax.boxplot([y], positions=[1], widths=.33, patch_artist=True, boxprops=dict(facecolor="#D8C9DF", color="#777777"), medianprops=dict(color=ink, lw=1.5), whiskerprops=dict(color="#777777"), capprops=dict(color="#777777"))
rng = np.random.default_rng(7)
ax.scatter(1 + rng.normal(0, .045, len(y)), y, s=18, color=blue, alpha=.75, edgecolor="white", linewidth=.3)
ax.set_xlim(.65, 1.35)
ax.set_xticks([1])
ax.set_xticklabels(["31 buildings"])
ax.set_ylabel("Building-level packing shortfall (%)")
ax.set_title("(b) Building-level distribution")
ax.grid(axis="y", alpha=.15)
median = np.median(y)
ax.text(1.18, median, f"Median {median:.1f}%", va="center", fontsize=9)
fig.suptitle("Explicit module geometry reduces the continuous-area capacity reference", fontweight="bold", y=1.02)
fig.tight_layout()
save(fig, "Fig5_capacity_distribution")

# Figure 6: facet density association
fig, ax = plt.subplots(figsize=(7.8, 5.45))
x = base.facet_density_per_100m2_gross3d
shortfall = base.packing_shortfall_pct
colour = base.median_facet_area_m2
sc = ax.scatter(x, shortfall, c=colour, cmap="viridis", s=42, edgecolor="white", linewidth=.5)
z = np.polyfit(x, shortfall, 1)
xx = np.linspace(x.min(), x.max(), 100)
ax.plot(xx, np.polyval(z, xx), "--", color="#666666", lw=1.2)
ax.set_xlabel("Roof-facet density (facets per 100 m² gross 3D roof area)")
ax.set_ylabel("Packing shortfall (%)")
ax.set_title("Building-level shortfall increases with roof-facet density", fontweight="bold")
ax.grid(True, alpha=.15)
cb = fig.colorbar(sc, ax=ax)
cb.set_label("Median physical-facet area (m²)")
ax.text(.98, .97, "Spearman ρ = 0.631\nFDR-adjusted p = 0.0010", transform=ax.transAxes, ha="right", va="top", fontsize=9, bbox=dict(boxstyle="round,pad=.3", fc="white", ec="#cccccc"))
fig.tight_layout()
save(fig, "Fig6_empirical_fragmentation")

# Figure 7: sensitivity and influence checks. Tilt subsets come from base clipping.
sens = pd.read_csv(RES / "sensitivity_summary.csv").set_index("scenario")
leave = pd.read_csv(RES / "leave_largest_out.csv")
alignment = pd.read_csv(RES / "alignment_sensitivity.csv")
fig, axs = plt.subplots(2, 2, figsize=(10.7, 7.7))
ax = axs[0, 0]
labels = ["0 m setback", "0.30 m\nbase", "0.50 m", "575 W\nmodule"]
values = [sens.loc["setback0", "aggregate_shortfall_pct"], sens.loc["base", "aggregate_shortfall_pct"], sens.loc["setback05", "aggregate_shortfall_pct"], sens.loc["large575", "aggregate_shortfall_pct"]]
ax.bar(labels, values, color=[blue, blue, green, orange])
ax.set_ylim(0, 20)
ax.set_ylabel("Aggregate shortfall (%)")
ax.set_title("(a) Setback and module size")
ax.grid(axis="y", alpha=.15)
ax = axs[0, 1]
labels = ["0 mm\n(theoretical)", "10 mm\nbase", "20 mm", "40 mm"]
values = [sens.loc["gap0", "aggregate_shortfall_pct"], sens.loc["base", "aggregate_shortfall_pct"], sens.loc["gap02", "aggregate_shortfall_pct"], sens.loc["gap04", "aggregate_shortfall_pct"]]
ax.bar(labels, values, color=[blue, blue, green, orange])
ax.set_ylim(0, 20)
ax.set_title("(b) Inter-module clearance")
ax.grid(axis="y", alpha=.15)
ax = axs[1, 0]
labels = ["All facets", "Pitched ≥5°", "Near-flat <5°"]
values = [sens.loc["base", "aggregate_shortfall_pct"], sens.loc["pitched", "aggregate_shortfall_pct"], sens.loc["flat", "aggregate_shortfall_pct"]]
ax.bar(labels, values, color=[blue, green, purple])
ax.set_ylim(0, 20)
ax.set_ylabel("Packing shortfall (%)")
ax.set_title("(c) Roof-tilt subsets")
ax.grid(axis="y", alpha=.15)
ax = axs[1, 1]
labels = ["Base", "− largest", "− 3 largest", "IoU ≥0.80", "IoU ≥0.90"]
values = [sens.loc["base", "aggregate_shortfall_pct"], leave.loc[leave.removed_largest_n == 1, "aggregate_shortfall_pct"].iloc[0], leave.loc[leave.removed_largest_n == 3, "aggregate_shortfall_pct"].iloc[0], alignment.loc[alignment.iou_threshold == .8, "aggregate_shortfall_pct"].iloc[0], alignment.loc[alignment.iou_threshold == .9, "aggregate_shortfall_pct"].iloc[0]]
ax.plot(labels, values, "-o", color=ink, lw=1.3, ms=4)
ax.set_ylim(0, 20)
ax.set_title("(d) Sample and alignment checks")
ax.grid(axis="y", alpha=.15)
ax.tick_params(axis="x", rotation=18)
fig.suptitle("Sensitivity of the Stuttgart packing shortfall to modelling and sample choices", fontweight="bold", y=.995)
fig.tight_layout(rect=[0, 0, 1, .97])
save(fig, "Fig7_sensitivity_robustness")

print(f"Analytical figures regenerated in {OUT}")
print("Figure 2 retained as static contextual imagery (not a computational input).")
