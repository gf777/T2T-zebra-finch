import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# -----------------------------
# Load + clean
# -----------------------------
df = pd.read_csv("Supplementary_Table_9.csv")
headers = df.iloc[0]
data = df[1:].copy()
data.columns = headers
data = data.reset_index(drop=True)

for c in ["Size", "GC%"]:
    data[c] = pd.to_numeric(data[c], errors="coerce")

data["DotLabel"] = data["Dot"].map({"TRUE": "Dot", "FALSE": "Non-dot"})

plot_df = data[["Name", "Size", "GC%", "DotLabel", "Classification_size"]].dropna(
    subset=["Name", "Size", "GC%", "DotLabel"]
).copy()

# haplotype from suffix; label cleanup (strip chr + suffix)
plot_df["Hap"] = np.where(plot_df["Name"].str.endswith("_mat"), "mat",
                   np.where(plot_df["Name"].str.endswith("_pat"), "pat", "unk"))

plot_df["LabelName"] = (
    plot_df["Name"]
    .str.replace("^chr", "", regex=True)
    .str.replace(r"(_mat|_pat)$", "", regex=True)
)

# sex chromosomes: Z/W (they are non-dot macros; no special color/shape)
plot_df["IsSex"] = plot_df["LabelName"].str.match(r"^(Z|W)$", na=False)

# macro chromosomes: driven ONLY by Classification_size
plot_df["IsMacro"] = plot_df["Classification_size"].astype(str).str.contains("macro", case=False, na=False)

# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(9, 6))
ax.set_xscale("log")
ax.set_xlabel("Chromosome size (bp, log scale)")
ax.set_ylabel("GC%")
ax.set_title("GC% vs Chromosome Size")
ax.grid(True)

# two colors: Dot vs Non-dot (use default cycle)
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
dot_color = colors[0]
non_color = colors[1]

marker = {"mat": "o", "pat": "s"}  # maternal circle, paternal square

# scatter all points (including Z/W) with correct dot/non-dot color and hap shape
for dotlab, col, alpha in [("Non-dot", non_color, 0.7), ("Dot", dot_color, 0.9)]:
    for hap in ["mat", "pat"]:
        d = plot_df[(plot_df["DotLabel"] == dotlab) & (plot_df["Hap"] == hap)]

        ax.scatter(
            d["Size"], d["GC%"],
            color=col,
            marker=marker[hap],
            alpha=alpha,
            edgecolor=np.where(d["IsMacro"], "black", "none"),
            linewidth=np.where(d["IsMacro"], 0.8, 0.0)
        )

# -----------------------------
# Add ONLY W/Z labels (no special styling)
# -----------------------------
sex = plot_df[plot_df["IsSex"]]
for _, r in sex.iterrows():
    ax.annotate(
        r["LabelName"],  # "Z" or "W"
        (r["Size"], r["GC%"]),
        textcoords="offset points",
        xytext=(6, 6),
        fontsize=10,
        fontweight="bold"
    )

# -----------------------------
# Legend (semantic)
# -----------------------------
legend_items = [
    Line2D([0], [0], marker="o", color="w", label="Maternal", markerfacecolor="k", markersize=7),
    Line2D([0], [0], marker="s", color="w", label="Paternal", markerfacecolor="k", markersize=7),
    Line2D([0], [0], marker="o", color="w", label="Dot", markerfacecolor=dot_color, markersize=7),
    Line2D([0], [0], marker="o", color="w", label="Non-dot", markerfacecolor=non_color, markersize=7),
    Line2D([0], [0], marker="o", color="w", label="Macro outline",
           markerfacecolor="white", markeredgecolor="black", markersize=7, lw=0),
]
ax.legend(handles=legend_items, frameon=True)

# -----------------------------
# Gentle label repel (for labels < 10.1 Mbp; excludes Z/W)
# -----------------------------
def repel(ax, texts, anchors, iters=200, step=0.4, k_rep=180.0, k_att=0.5, max_dist=80.0):
    if not texts:
        return
    fig = ax.figure
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()

    to_disp = ax.transData.transform
    to_data = ax.transData.inverted().transform
    A = to_disp(np.asarray(anchors, float))

    base = np.array([(10, 8), (10, -8), (-10, 8), (-10, -8),
                     (14, 0), (-14, 0), (0, 12), (0, -12)], float)
    P = A + np.resize(base, (len(A), 2))

    for i, t in enumerate(texts):
        t.set_position(tuple(to_data(P[i])))

    axbox = ax.get_window_extent(ren)

    for _ in range(iters):
        fig.canvas.draw()
        bbs = [t.get_window_extent(ren).expanded(1.02, 1.10) for t in texts]
        F = np.zeros_like(P)
        overlap = False

        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                if bbs[i].overlaps(bbs[j]):
                    overlap = True
                    ox = min(bbs[i].x1, bbs[j].x1) - max(bbs[i].x0, bbs[j].x0)
                    oy = min(bbs[i].y1, bbs[j].y1) - max(bbs[i].y0, bbs[j].y0)
                    v = P[i] - P[j]
                    n = np.hypot(v[0], v[1]) + 1e-9
                    v /= n
                    mag = k_rep * (ox * oy) / (n + 30.0)
                    F[i] += v * mag
                    F[j] -= v * mag

        F += k_att * (A - P)

        n = np.hypot(F[:, 0], F[:, 1]) + 1e-9
        P += (F / n[:, None]) * step

        v = P - A
        d = np.hypot(v[:, 0], v[:, 1]) + 1e-9
        m = d > max_dist
        if np.any(m):
            P[m] = A[m] + v[m] * (max_dist / d[m])[:, None]

        for i, bb in enumerate(bbs):
            cx, cy = P[i]
            w, h = bb.width, bb.height
            cx = np.clip(cx, axbox.x0 + w/2 + 2, axbox.x1 - w/2 - 2)
            cy = np.clip(cy, axbox.y0 + h/2 + 2, axbox.y1 - h/2 - 2)
            P[i] = (cx, cy)

        for i, t in enumerate(texts):
            t.set_position(tuple(to_data(P[i])))

        if not overlap:
            break

# label everything < 10.1 Mbp, but skip Z/W since they have special label already
lab = plot_df[(plot_df["Size"] < 10.1e6) & (~plot_df["IsSex"])]

texts, anchors = [], []
for _, r in lab.iterrows():
    anchors.append((r["Size"], r["GC%"]))
    texts.append(ax.text(r["Size"], r["GC%"], r["LabelName"], fontsize=9, ha="left", va="bottom"))

repel(ax, texts, anchors)

# leader lines
for t, (x, y) in zip(texts, anchors):
    tx, ty = t.get_position()
    ax.annotate("", xy=(x, y), xytext=(tx, ty),
                arrowprops=dict(arrowstyle="-", lw=0.6, alpha=0.6))

plt.tight_layout()
plt.savefig("zebra_finch_gc_vs_size_final.svg", bbox_inches="tight")
plt.show()
