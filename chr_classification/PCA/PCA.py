#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

mpl.rcParams["svg.fonttype"] = "none"

# ----------------------------
# Helpers: PCA (NumPy-only) + label de-overlap (no extra deps)
# ----------------------------
def pca_svd(X: np.ndarray, n_components: int = 2):
    """
    PCA via SVD on a standardized matrix X (rows=samples, cols=features).
    Returns scores (PC coordinates), explained variance ratio, loadings.
    """
    U, S, Vt = np.linalg.svd(X, full_matrices=False)

    scores = U[:, :n_components] * S[:n_components]

    n = X.shape[0]
    eigs = (S**2) / (n - 1)
    explained = eigs / eigs.sum()

    loadings = Vt[:n_components, :].T  # features x components
    return scores, explained[:n_components], loadings


def repel_texts(ax, texts, max_iter=250, step=2):
    """
    Very lightweight text overlap reducer without external libraries.
    It nudges labels in display (pixel) coords until they stop overlapping.

    - texts: list of matplotlib.text.Text already added to ax
    - step: pixel shift per iteration when overlaps occur
    """
    if not texts:
        return

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    for _ in range(max_iter):
        moved = False
        bboxes = [t.get_window_extent(renderer).expanded(1.05, 1.10) for t in texts]

        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                if bboxes[i].overlaps(bboxes[j]):
                    moved = True

                    ci = bboxes[i].get_points().mean(axis=0)
                    cj = bboxes[j].get_points().mean(axis=0)
                    dx, dy = cj - ci

                    if abs(dx) < 1e-6:
                        dx = 1.0
                    if abs(dy) < 1e-6:
                        dy = 1.0

                    sx = step if dx >= 0 else -step
                    sy = step if dy >= 0 else -step

                    xj, yj = texts[j].get_position()
                    x_disp, y_disp = ax.transData.transform((xj, yj))
                    x_disp += sx
                    y_disp += sy
                    x_new, y_new = ax.transData.inverted().transform((x_disp, y_disp))
                    texts[j].set_position((x_new, y_new))

                    fig.canvas.draw()
                    renderer = fig.canvas.get_renderer()
                    bboxes[j] = texts[j].get_window_extent(renderer).expanded(1.05, 1.10)

        if not moved:
            break


# ----------------------------
# Read data (IMPORTANT: skip the title row)
# ----------------------------
data = pd.read_csv("Supplementary_Table_9.csv", skiprows=1)

# Coerce numeric columns
num_cols = ["Size", "GC%", "Methylation%", "Genes", "Repeats", "Minisatellite"]
for c in num_cols:
    data[c] = pd.to_numeric(data[c], errors="coerce")

# Flags from table
data["IsDot"] = data["Dot"].astype(str).str.upper().eq("TRUE")
data["SizeClass"] = data["Classification_size"].astype(str).str.strip().str.capitalize()  # Macro/Micro

# Sex chromosomes (adjust regex if your naming differs)
data["IsSex"] = data["Name"].astype(str).str.contains(r"chr[ZW]\b|chr[ZW]_", case=False, na=False)

# ----------------------------
# SIZE CORRECTION: use densities per Mb (do NOT include Size in PCA features)
# ----------------------------
mb = data["Size"] / 1e6
data["Genes_per_Mb"] = data["Genes"] / mb
data["Repeats_per_Mb"] = data["Repeats"] / mb
data["Minisatellite_per_Mb"] = data["Minisatellite"] / mb

features = ["GC%", "Methylation%", "Genes_per_Mb", "Repeats_per_Mb", "Minisatellite_per_Mb"]

sub = data[["Name", "Size", "IsDot", "SizeClass", "IsSex"] + features].dropna().copy()

# Standardize features (z-score)
X = sub[features].to_numpy(float)
mu = X.mean(axis=0)
sd = X.std(axis=0, ddof=1)
sd[sd == 0] = 1.0
Xz = (X - mu) / sd

# PCA (NumPy-only)
scores, explained, loadings = pca_svd(Xz, n_components=2)
sub["PC1"] = scores[:, 0]
sub["PC2"] = scores[:, 1]

# Diagnostics: how much does PC1 still correlate with Size?
pc1_size_corr = np.corrcoef(sub["PC1"], np.log10(sub["Size"]))[0, 1]
print("Explained variance (PC1, PC2):", explained)
print("corr(PC1, log10(Size)):", pc1_size_corr)

# ----------------------------
# Plot (updated per your specs)
#   - color: Dot vs Non-dot (ONLY two colors total)
#   - marker: Maternal vs Paternal (o vs s)
#   - macrochromosomes: black outline
#   - labels: strip chr + _mat/_pat and label by chromosome token (incl Z/W)
#   - legend: separate color legend and marker legend
# ----------------------------
fig, ax = plt.subplots(figsize=(10, 7))

# --- Infer maternal/paternal (robust to different column names / encodings)
parent_col = None
for cand in ["Parent", "Haplotype", "Phase", "Origin"]:
    if cand in sub.columns:
        parent_col = cand
        break

def infer_parent(row):
    if parent_col is not None:
        v = str(row[parent_col]).strip().lower()
        if any(k in v for k in ["mat", "maternal", "mother", "m"]):
            return "maternal"
        if any(k in v for k in ["pat", "paternal", "father", "p"]):
            return "paternal"

    n = str(row["Name"]).strip().lower()
    if any(k in n for k in ["_mat", ".mat", "mat_", "maternal"]):
        return "maternal"
    if any(k in n for k in ["_pat", ".pat", "pat_", "paternal"]):
        return "paternal"

    return "maternal"

sub["ParentClass"] = sub.apply(infer_parent, axis=1)

marker_for_parent = {"maternal": "o", "paternal": "s"}

# --- Two colors only: dot vs non-dot (use default cycle entries 0 and 1)
default_colors = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
color_map = {False: default_colors[1], True: default_colors[0]}  # Non-dot, Dot

# --- Scatter with macro black outline
for is_dot in [False, True]:
    for parent in ["maternal", "paternal"]:
        ss = sub[(sub["IsDot"] == is_dot) & (sub["ParentClass"] == parent)]
        if ss.empty:
            continue

        is_macro = ss["SizeClass"].astype(str).eq("Macro")
        ss_macro = ss[is_macro]
        ss_micro = ss[~is_macro]

        if not ss_micro.empty:
            ax.scatter(
                ss_micro["PC1"],
                ss_micro["PC2"],
                c=[color_map[is_dot]],
                marker=marker_for_parent[parent],
                alpha=0.75,
                edgecolors="none",
                linewidths=0.0,
            )

        if not ss_macro.empty:
            ax.scatter(
                ss_macro["PC1"],
                ss_macro["PC2"],
                c=[color_map[is_dot]],
                marker=marker_for_parent[parent],
                alpha=0.75,
                edgecolors="black",
                linewidths=0.9,
            )

# --- Label helper: strip "chr" and trailing hap suffixes
def short_chr_label(name: str) -> str:
    s = str(name).strip()

    # remove leading 'chr' (case-insensitive)
    if s.lower().startswith("chr"):
        s = s[3:]

    # drop common hap/phase suffixes (case-insensitive)
    s_low = s.lower()
    for suf in ["_mat", "_pat", ".mat", ".pat", "_maternal", "_paternal", "_m", "_p"]:
        if s_low.endswith(suf):
            s = s[: -len(suf)]
            s_low = s.lower()
            break

    # if there are still underscores, take the first token (e.g. "Z_something" -> "Z")
    if "_" in s:
        s = s.split("_", 1)[0]

    return s

# --- Label dot chromosomes (as short chr token), then repel
texts = []
dots = sub[sub["IsDot"]]
for _, r in dots.iterrows():
    lbl = short_chr_label(r["Name"])
    texts.append(ax.text(r["PC1"], r["PC2"], lbl, fontsize=9))

# --- Label sex chromosomes directly with Z/W (also short label)
sex = sub[sub["IsSex"]].copy()
for _, r in sex.iterrows():
    lbl = short_chr_label(r["Name"])
    # enforce exact "Z" / "W" if present anywhere in name
    nm = str(r["Name"]).lower()
    if "chrz" in nm or "_z" in nm or " z" in nm:
        lbl = "Z"
    if "chrw" in nm or "_w" in nm or " w" in nm:
        lbl = "W"
    t = ax.text(r["PC1"], r["PC2"], lbl, fontsize=12, fontweight="bold")
    texts.append(t)

repel_texts(ax, texts, max_iter=350, step=3)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title(
    "PCA (size-corrected using densities per Mb)\n"
    "Color = Dot vs Non-dot, Shape = Maternal (o) vs Paternal (s), Macro outlined"
)
ax.grid(True)

# ----------------------------
# Legend: separate color meaning from marker meaning
# ----------------------------
color_handles = [
    Patch(facecolor=color_map[False], edgecolor="none", label="Non-dot"),
    Patch(facecolor=color_map[True], edgecolor="none", label="Dot"),
]
shape_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor="gray",
           markeredgecolor="black", markersize=9, linestyle="None", label="Maternal"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="gray",
           markeredgecolor="black", markersize=9, linestyle="None", label="Paternal"),
]

leg1 = ax.legend(handles=color_handles, title="Color", loc="upper right", frameon=True)
ax.add_artist(leg1)
ax.legend(handles=shape_handles, title="Shape", loc="lower right", frameon=True)

plt.tight_layout()

# Save outputs
out_png = "zebra_finch_chr_pca_density_dot_nondot_parent_macro_outline.png"
out_svg = "zebra_finch_chr_pca_density_dot_nondot_parent_macro_outline.svg"

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_svg, format="svg", bbox_inches="tight")
plt.show()

print("Saved:", out_png)
print("Saved:", out_svg)
