#!/usr/bin/env python3
"""
Main manuscript figures (Figure 1-6) for the DRG pan-cancer study.
Reads   <results>/tables/*.csv   produced by drg_pancancer_pipeline.py (v2)
Writes  <results>/figures_main/Figure{1..6}.{png,pdf,tiff}  (600 dpi)

Usage:  python make_main_figures.py --results "/content/drive/MyDrive/DRG_PanCancer/results"
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

# ── style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
})
INK, INK2, GRID = "#1F1F1F", "#5C5C5C", "#D9D9D9"
C_FULL, C_CORE, C_CYTO = "#7B4FA0", "#D55E00", "#0072B2"      # validated categorical (CVD-safe)
RISK, PROT, NS = "#B2182B", "#2166AC", "#A6A6A6"
DIV = "RdBu_r"

CORE = ["SLC7A11", "SLC3A2", "RPN1", "NCKAP1", "NUBPL", "NDUFA11", "LRPPRC", "OXSM", "NDUFS1", "GYS1"]
CYTO = ["FLNA", "FLNB", "MYH9", "TLN1", "ACTB", "MYL6", "MYH10", "CAPZB", "DSTN", "IQGAP1",
        "ACTN4", "PDLIM1", "CD2AP", "INF2"]
GENES = CORE + CYTO
ARM_COL = {**{g: C_CORE for g in CORE}, **{g: C_CYTO for g in CYTO}}
PW_LABEL = {
    "Pentose_phosphate_NADPH": "Pentose phosphate / NADPH", "Glycogen_metabolism": "Glycogen metabolism",
    "Glycolysis": "Glycolysis", "Hypoxia_core": "Hypoxia", "mTORC1_core": "mTORC1",
    "Fatty_acid_synthesis": "Fatty acid synthesis", "Glutamine_metabolism": "Glutamine metabolism",
    "TCA_cycle": "TCA cycle", "Glutathione_metabolism": "Glutathione metabolism",
    "Serine_one_carbon": "Serine / one-carbon", "Glucose_uptake": "Glucose uptake",
    "NRF2_antioxidant_response": "NRF2 response", "Thioredoxin_system": "Thioredoxin system",
    "Fatty_acid_oxidation": "Fatty acid oxidation", "Oxidative_phosphorylation": "Oxidative phosphorylation",
    "Cysteine_sulfur_metabolism": "Cysteine / sulfur metabolism"}


def stars(q):
    if pd.isna(q):
        return ""
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else ""


def letter(ax, s, dx=-0.02, dy=0.0):
    pos = ax.get_position()
    ax.figure.text(pos.x0 + dx, pos.y1 + dy + 0.008, s, fontsize=10, fontweight="bold",
                   va="bottom", ha="right", color=INK)


def save(fig, out, name):
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    try:
        fig.savefig(out / f"{name}.tiff", bbox_inches="tight", facecolor="white",
                    pil_kwargs={"compression": "tiff_lzw"})
    except Exception:  # noqa: BLE001
        pass
    plt.close(fig)


def heat(ax, M, Q=None, cmap=DIV, vlim=None, center=0, cbar_label="", cax=None, star_size=4.5,
         xrot=90, sequential=False, vmin=None):
    """Heatmap with white cell gaps and optional significance stars; returns mappable."""
    A = M.values.astype(float)
    if sequential:
        norm = plt.Normalize(vmin=0 if vmin is None else vmin, vmax=vlim or np.nanmax(A))
    else:
        lim = vlim or np.nanpercentile(np.abs(A[np.isfinite(A)]), 98)
        norm = TwoSlopeNorm(vcenter=center, vmin=-lim, vmax=lim)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("#F2F2F2")
    im = ax.pcolormesh(np.ma.masked_invalid(A), cmap=cm, norm=norm, edgecolors="white", linewidth=0.4)
    ax.set_xlim(0, A.shape[1])
    ax.set_ylim(A.shape[0], 0)
    ax.set_xticks(np.arange(A.shape[1]) + .5)
    ax.set_xticklabels(M.columns, rotation=xrot)
    ax.set_yticks(np.arange(A.shape[0]) + .5)
    ax.set_yticklabels(M.index)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    if Q is not None:
        Q = Q.reindex(index=M.index, columns=M.columns)
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                s = stars(Q.iat[i, j])
                if s:
                    v = norm(A[i, j]) if np.isfinite(A[i, j]) else .5
                    ax.text(j + .5, i + .62, s, ha="center", va="center", fontsize=star_size,
                            color="white" if (v > .8 or v < .2) else INK)
    if cax is not None:
        cb = plt.colorbar(im, cax=cax)
        cb.set_label(cbar_label, fontsize=6)
        cb.ax.tick_params(labelsize=5.5, length=2)
        cb.outline.set_visible(False)
    return im


def arm_strip(ax, genes, side="left"):
    for i, g in enumerate(genes):
        ax.add_patch(plt.Rectangle((0, i), 1, 1, color=ARM_COL.get(g, NS), lw=0))
    ax.set_xlim(0, 1)
    ax.set_ylim(len(genes), 0)
    ax.axis("off")


def km_panel(ax, d, ep, title):
    d = d[[ep, f"{ep}.time", "score"]].dropna()
    d = d[d[f"{ep}.time"] > 0]
    hi = (d["score"] > d["score"].median()).astype(int)
    _, p = survdiff(d[f"{ep}.time"].values, d[ep].values.astype(int), hi.values)
    for g, lab, col in ((0, "Low", PROT), (1, "High", RISK)):
        m = hi.values == g
        sf = SurvfuncRight(d[f"{ep}.time"].values[m], d[ep].values[m])
        ax.step(np.r_[0, sf.surv_times], np.r_[1, sf.surv_prob], where="post", color=col, lw=1.1,
                label=f"{lab} (n={m.sum()})")
    ax.set_ylim(0, 1.03)
    ax.set_xlim(0, min(d[f"{ep}.time"].quantile(.98), 15))
    ax.set_xlabel("Years")
    ax.set_ylabel("Overall survival")
    ps = f"P = {p:.1e}" if p < 1e-3 else f"P = {p:.3f}"
    ax.set_title(f"{title} (n = {len(d)})\nlog-rank {ps}", fontsize=6.5, loc="left", color=INK, linespacing=1.3)
    ax.grid(axis="y", color=GRID, lw=.4)


def forest(ax, df, hr, lo, hi, p, labels, meta=None, xlabel="HR per SD"):
    df = df.reset_index(drop=True)
    y = np.arange(len(df))
    col = np.where(df[p] < .05, np.where(df[hr] > 1, RISK, PROT), NS)
    ax.hlines(y, df[lo], df[hi], color=col, lw=1)
    ax.scatter(df[hr], y, c=col, s=10, zorder=3, edgecolors="white", linewidths=.4)
    yt, yl = list(y), list(labels)
    if meta is not None:
        ym = len(df) + .8
        ax.fill([meta["HR_low"], meta["HR"], meta["HR_high"], meta["HR"]],
                [ym, ym - .4, ym, ym + .4], color=INK, lw=0)
        yt.append(ym)
        yl.append(f"Random-effects (I² = {meta['I2']:.0f}%)")
        ax.text(meta["HR_high"] * 1.04, ym, f"{meta['HR']:.2f} ({meta['HR_low']:.2f}–{meta['HR_high']:.2f})",
                va="center", fontsize=5.5, color=INK)
    ax.axvline(1, color=INK2, lw=.6, ls="--")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks([0.5, 0.7, 1, 1.5, 2, 3])
    ax.set_yticks(yt)
    ax.set_yticklabels(yl)
    ax.set_ylim(max(yt) + .8, -.8)
    ax.set_xlabel(xlabel)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, lw=.4)


def nudge_labels(ax, xs, ys, labels, min_gap):
    """Simple vertical repulsion so annotation labels do not overlap."""
    order = np.argsort(ys)
    ly = np.array(ys, float)[order].copy()
    for _ in range(200):
        moved = False
        for i in range(1, len(ly)):
            if ly[i] - ly[i - 1] < min_gap:
                shift = (min_gap - (ly[i] - ly[i - 1])) / 2
                ly[i - 1] -= shift
                ly[i] += shift
                moved = True
        if not moved:
            break
    for k, i in enumerate(order):
        ax.annotate(labels[i], (xs[i], ys[i]), xytext=(xs[i] + 0.012, ly[k]), fontsize=5.3, color=INK,
                    va="center", arrowprops=dict(arrowstyle="-", lw=.35, color=INK2, shrinkA=0, shrinkB=1.5))


# ════════════════════════════════════════════════════════════════════════════
def figure1(T, out):
    tn = T["tn"]
    g = tn[tn.gene != "DRG_score"]
    fc = g.pivot(index="gene", columns="cancer", values="log2FC").reindex(GENES)
    q = g.pivot(index="gene", columns="cancer", values="FDR").reindex(GENES)
    sig = g[g.FDR < .05].assign(up=lambda d: d.log2FC > 0)
    up = sig[sig.up].groupby("gene").size().reindex(GENES).fillna(0)
    dn = sig[~sig.up].groupby("gene").size().reindex(GENES).fillna(0)
    co = T["coex"].reindex(index=GENES, columns=GENES)

    fig = plt.figure(figsize=(7.2, 8.6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.05], hspace=0.32)
    top = gridspec.GridSpecFromSubplotSpec(1, 4, gs[0], width_ratios=[.03, 1, .32, .025], wspace=.05)
    axs, axA, axB, cA = (fig.add_subplot(top[i]) for i in range(4))
    arm_strip(axs, GENES)
    heat(axA, fc, q, vlim=2.5, cbar_label="log2 fold change\n(tumour − normal)", cax=cA, star_size=4)
    axA.set_yticklabels([])
    axs.set_ylim(len(GENES), 0)
    for i, gname in enumerate(GENES):
        axs.text(-0.4, i + .5, gname, ha="right", va="center", fontsize=6, fontstyle="italic",
                 transform=axs.transData)
    axA.axhline(len(CORE), color=INK, lw=.8)
    axA.set_title("Tumour vs adjacent normal (BH-FDR: * <0.05, ** <0.01, *** <0.001)", loc="left")
    y = np.arange(len(GENES)) + .5
    axB.barh(y, up.values, height=.7, color=RISK, label="Up in tumour")
    axB.barh(y, -dn.values, height=.7, color=PROT, label="Down in tumour")
    axB.set_ylim(len(GENES), 0)
    axB.set_yticks([])
    axB.axvline(0, color=INK, lw=.5)
    axB.spines["left"].set_visible(False)
    lim = max(up.max(), dn.max()) + 1
    axB.set_xlim(-lim, lim)
    axB.set_xticks([-15, -10, -5, 0, 5, 10, 15])
    axB.set_xticklabels(["15", "10", "5", "0", "5", "10", "15"])
    axB.set_xlabel("No. of cancers\n(FDR < 0.05)")
    axB.legend(loc="lower center", bbox_to_anchor=(.5, 1.0), frameon=False, ncol=1, fontsize=5.5)
    axB.grid(axis="x", color=GRID, lw=.4)

    bot = gridspec.GridSpecFromSubplotSpec(1, 3, gs[1], width_ratios=[1, .03, .6], wspace=.12)
    axC, cC = fig.add_subplot(bot[0]), fig.add_subplot(bot[1])
    heat(axC, co, None, vlim=.6, cbar_label="median Spearman ρ\nacross 33 cancers", cax=cC)
    for lab in axC.get_xticklabels() + axC.get_yticklabels():
        lab.set_color(ARM_COL.get(lab.get_text(), INK))
        lab.set_fontstyle("italic")
    axC.axhline(len(CORE), color=INK, lw=.8)
    axC.axvline(len(CORE), color=INK, lw=.8)
    axC.set_title("Co-expression of DRGs within tumours", loc="left")
    axL = fig.add_subplot(bot[2])
    axL.axis("off")
    axL.legend(handles=[Patch(color=C_CORE, label="Core subset (CRISPR-screen genes)"),
                        Patch(color=C_CYTO, label="Cytoskeleton subset (disulfide-stressed\nactin cytoskeleton proteins)")],
               loc="upper left", frameon=False, fontsize=6, title="Gene colour", title_fontsize=6)
    fig.canvas.draw()
    letter(axs, "A", dx=-.07)
    letter(axB, "B", dx=0)
    letter(axC, "C", dx=-.07)
    save(fig, out, "Figure1")


def figure2(T, out):
    sc = T["sc"]
    os_ = sc[sc.endpoint == "OS"].sort_values("HR_perSD", ascending=False)
    meta = T["meta"].set_index("endpoint").loc["OS"]
    clin = T["clin"]
    gc = T["gcox"]
    g = gc[gc.endpoint == "OS"]
    lhr = np.log2(g.pivot(index="gene", columns="cancer", values="HR")).reindex(GENES)
    q = g.pivot(index="gene", columns="cancer", values="FDR").reindex(GENES)

    fig = plt.figure(figsize=(7.2, 8.8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.15, .85], hspace=.34)
    top = gridspec.GridSpecFromSubplotSpec(1, 2, gs[0], width_ratios=[1, 1.15], wspace=.45)
    axA = fig.add_subplot(top[0])
    forest(axA, os_, "HR_perSD", "HR_low", "HR_high", "p_cox",
           [f"{c} ({int(e)}/{int(n)})" for c, e, n in zip(os_.cancer, os_.events, os_.n)], meta=meta,
           xlabel="HR per SD of DRG score (95% CI)")
    axA.set_title("Overall survival, univariate Cox (events/n)", loc="left")
    axA.legend(handles=[Line2D([], [], marker="o", ls="", color=RISK, ms=3.5, label="Risk (P < 0.05)"),
                        Line2D([], [], marker="o", ls="", color=PROT, ms=3.5, label="Protective (P < 0.05)"),
                        Line2D([], [], marker="o", ls="", color=NS, ms=3.5, label="NS")],
               loc="upper center", bbox_to_anchor=(.45, -.075), ncol=3, frameon=False, fontsize=5.5,
               handletextpad=.1, columnspacing=.8)
    kmg = gridspec.GridSpecFromSubplotSpec(3, 2, top[1], hspace=.95, wspace=.45)
    picks = [c for c in ["LGG", "MESO", "UVM", "HNSC", "LUAD", "KIRC"] if c in set(clin.cancer)]
    km_axes = []
    for k, ca in enumerate(picks[:6]):
        ax = fig.add_subplot(kmg[k // 2, k % 2])
        km_panel(ax, clin[clin.cancer == ca], "OS", ca)
        if k % 2:
            ax.set_ylabel("")
        if k < 4:
            ax.set_xlabel("")
        km_axes.append(ax)
    km_axes[0].legend(handles=[Line2D([], [], color=PROT, lw=1.2, label="Low DRG score"),
                               Line2D([], [], color=RISK, lw=1.2, label="High DRG score")],
                      loc="lower left", bbox_to_anchor=(0, 1.42), ncol=2, frameon=False, fontsize=6,
                      handlelength=1.5)
    bot = gridspec.GridSpecFromSubplotSpec(1, 3, gs[1], width_ratios=[.025, 1, .02], wspace=.03)
    axs, axC, cC = (fig.add_subplot(bot[i]) for i in range(3))
    arm_strip(axs, GENES)
    for i, gname in enumerate(GENES):
        axs.text(-0.4, i + .5, gname, ha="right", va="center", fontsize=5.5, fontstyle="italic")
    heat(axC, lhr, q, vlim=1.2, cbar_label="log2 HR per SD (OS)", cax=cC, star_size=4)
    axC.set_yticklabels([])
    axC.axhline(len(CORE), color=INK, lw=.8)
    axC.set_title("Gene-level univariate Cox regression for overall survival (BH-FDR stars)", loc="left")
    fig.canvas.draw()
    letter(axA, "A", dx=-.09)
    letter(km_axes[0], "B", dx=-.04)
    letter(axs, "C", dx=-.06)
    save(fig, out, "Figure2")


def figure3(T, out):
    sub = T["sub"]
    cor = T["corarm"].sort_values("rho_core_vs_cytoskeleton")
    feats = ["DRG_score", "DRG_core", "DRG_cytoskeleton"]
    flab = {"DRG_score": "Full score", "DRG_core": "Core subset", "DRG_cytoskeleton": "Cytoskeleton subset"}
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = gridspec.GridSpec(1, 3, width_ratios=[.85, 1, 1], wspace=.95)
    axA = fig.add_subplot(gs[0])
    y = np.arange(len(cor))
    v = cor["rho_core_vs_cytoskeleton"].values
    axA.barh(y, v, height=.7, color=np.where(v >= 0, "#8C8C8C", "#BFBFBF"))
    axA.set_yticks(y)
    axA.set_yticklabels(cor.index)
    axA.set_ylim(len(cor) - .4, -.6)
    axA.axvline(0, color=INK, lw=.6)
    med = np.median(v)
    axA.axvline(med, color=C_FULL, lw=1, ls="--")
    axA.set_xlabel("Spearman ρ (core vs cytoskeleton\nsub-score, tumours)")
    axA.legend(handles=[Line2D([], [], color=C_FULL, lw=1, ls="--", label=f"median ρ = {med:.2f}")],
               loc="upper right", frameon=False, fontsize=6)
    axA.set_title("Subset independence", loc="left")
    axA.spines["left"].set_visible(False)
    axA.tick_params(axis="y", length=0)
    axA.grid(axis="x", color=GRID, lw=.4)
    axA.set_xlim(-.6, .6)

    def panel(spec, analysis, tr, lab, ttl, vlim):
        g2 = gridspec.GridSpecFromSubplotSpec(1, 2, spec, width_ratios=[1, .07], wspace=.08)
        ax, cax = fig.add_subplot(g2[0]), fig.add_subplot(g2[1])
        r = sub[sub.analysis == analysis]
        M = tr(r.pivot(index="cancer", columns="feature", values="effect"))[feats]
        Q = r.pivot(index="cancer", columns="feature", values="FDR")[feats]
        M.columns = [flab[c] for c in M.columns]
        Q.columns = M.columns
        heat(ax, M, Q, vlim=vlim, cbar_label=lab, cax=cax, xrot=40, star_size=4.5)
        for t, c in zip(ax.get_xticklabels(), [C_FULL, C_CORE, C_CYTO]):
            t.set_color(c)
            t.set_ha("right")
        ax.set_title(ttl, loc="left")
        return ax
    axB = panel(gs[1], "tumour_vs_normal", lambda x: x, "Δ ssGSEA score (tumour − normal)", "Tumour vs normal", .12)
    axC = panel(gs[2], "cox_OS", np.log2, "log2 HR per SD (OS)", "Overall survival", 1.0)
    fig.canvas.draw()
    letter(axA, "A", dx=-.05)
    letter(axB, "B", dx=-.04)
    letter(axC, "C", dx=-.04)
    save(fig, out, "Figure3")


def figure4(T, out):
    m = T["met"]
    s = m[m.feature == "DRG_score"]
    R = s.pivot(index="pathway", columns="cancer", values="rho_purity_adj")
    Q = s.pivot(index="pathway", columns="cancer", values="FDR_purity_adj")
    order = s.groupby("pathway")["rho_purity_adj"].median().sort_values(ascending=False).index
    R, Q = R.reindex(order), Q.reindex(order)
    R.index = [PW_LABEL.get(p, p) for p in R.index]
    Q.index = R.index
    arms = m[m.feature.isin(["DRG_core", "DRG_cytoskeleton"])].groupby(["pathway", "feature"])["rho_purity_adj"] \
        .median().unstack()
    arms = arms.assign(d=arms["DRG_core"] - arms["DRG_cytoskeleton"]).sort_values("d")
    raw = s.groupby("pathway")[["rho", "rho_purity_adj"]].median()

    fig = plt.figure(figsize=(7.2, 8.4))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, .95], hspace=.42)
    t = gridspec.GridSpecFromSubplotSpec(1, 2, gs[0], width_ratios=[1, .015], wspace=.02)
    axA, cA = fig.add_subplot(t[0]), fig.add_subplot(t[1])
    heat(axA, R, Q, vlim=.6, cbar_label="partial Spearman ρ\n(| ABSOLUTE purity)", cax=cA, star_size=4)
    axA.set_title("DRG score vs metabolic programmes in tumours, adjusted for tumour purity "
                  "(DRGs removed from pathway sets; BH-FDR stars)", loc="left")
    b = gridspec.GridSpecFromSubplotSpec(1, 3, gs[1], width_ratios=[1.25, .78, .42], wspace=.7)
    axB = fig.add_subplot(b[0])
    y = np.arange(len(arms))
    axB.hlines(y, arms["DRG_cytoskeleton"], arms["DRG_core"], color=GRID, lw=2.2, zorder=1)
    axB.scatter(arms["DRG_core"], y, color=C_CORE, s=22, zorder=3, edgecolors="white", linewidths=.6,
                label="Core subset")
    axB.scatter(arms["DRG_cytoskeleton"], y, color=C_CYTO, s=22, zorder=3, edgecolors="white", linewidths=.6,
                label="Cytoskeleton subset")
    axB.set_yticks(y)
    axB.set_yticklabels([PW_LABEL.get(p, p) for p in arms.index])
    axB.axvline(0, color=INK, lw=.6)
    axB.set_xlabel("median partial ρ across cancers")
    axB.set_title("Distinct metabolic coupling of the two subsets", loc="left")
    axB.legend(loc="upper center", bbox_to_anchor=(.45, -.12), ncol=2, frameon=False, fontsize=6)
    axB.spines["left"].set_visible(False)
    axB.tick_params(axis="y", length=0)
    axB.grid(axis="x", color=GRID, lw=.4)
    axC = fig.add_subplot(b[1])
    axC.scatter(raw["rho"], raw["rho_purity_adj"], s=14, color=C_FULL, edgecolors="white", linewidths=.5, zorder=3)
    lo, hi = -.12, .34
    axC.plot([lo, hi], [lo, hi], color=INK2, lw=.6, ls="--")
    axC.axhline(0, color=GRID, lw=.5)
    axC.axvline(0, color=GRID, lw=.5)
    axC.set_xlim(lo, hi + .02)
    axC.set_ylim(lo, hi)
    # label column to the right of the axis, leader lines to points
    xs, ys = raw["rho"].values, raw["rho_purity_adj"].values
    labs = [PW_LABEL.get(p, p) for p in raw.index]
    o = np.argsort(-ys)
    span = np.linspace(hi - .01, lo + .01, len(o))
    for k, i in enumerate(o):
        axC.annotate(labs[i], (xs[i], ys[i]), xytext=(1.04, (span[k] - lo) / (hi - lo)),
                     textcoords="axes fraction", fontsize=5.3, color=INK, va="center",
                     arrowprops=dict(arrowstyle="-", lw=.3, color="#9A9A9A", shrinkA=0, shrinkB=1.5))
    axC.set_xlabel("median ρ, unadjusted")
    axC.set_ylabel("median partial ρ, purity-adjusted")
    axC.set_title("Robustness to tumour purity", loc="left")
    fig.canvas.draw()
    letter(axA, "A", dx=-.06)
    letter(axB, "B", dx=-.17)
    letter(axC, "C", dx=-.05)
    save(fig, out, "Figure4")


FIG5_TABLE = False  # revision: former panel D is Table 2 of the manuscript


def figure5(T, out):
    mv = T["mv"]
    s = mv[(mv.covariate == "DRG_score") & (mv.endpoint == "OS")]
    c1 = s[s.model == "clinical"].set_index("cancer")
    c2 = s[s.model == "clinical+purity"].set_index("cancer")
    keep = sorted(set(c1.index[c1.p < .1]) | set(c2.index[c2.p < .1]), key=lambda c: -c1.HR.get(c, 1))
    perf = T["perf"]
    sp = T["splits"]
    stab = T["stab"]
    rob = perf.loc[perf.robust, "cancer"].tolist()

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, .9], width_ratios=[1, 1.05], hspace=.5, wspace=.55)
    axA = fig.add_subplot(gs[0, 0])
    y = np.arange(len(keep))
    for off, d, mk, lab in ((-.17, c1, "o", "Age + sex + stage"), (.17, c2, "s", "+ ABSOLUTE purity")):
        d = d.reindex(keep)
        col = np.where(d.p < .05, np.where(d.HR > 1, RISK, PROT), NS)
        axA.hlines(y + off, d.HR_low, d.HR_high, color=col, lw=.9)
        axA.scatter(d.HR, y + off, c=col, marker=mk, s=11, zorder=3, edgecolors="white", linewidths=.4)
        axA.scatter([], [], marker=mk, c=INK2, s=11, label=lab)
    axA.axvline(1, color=INK2, lw=.6, ls="--")
    axA.set_xscale("log")
    axA.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
    axA.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axA.set_xticks([0.7, 1, 1.5, 2, 3])
    axA.set_yticks(y)
    axA.set_yticklabels(keep)
    axA.set_ylim(len(keep) - .5, -.6)
    axA.set_xlabel("Adjusted HR per SD of DRG score (OS)")
    axA.set_title("Multivariable Cox regression", loc="left")
    axA.legend(loc="lower right", frameon=False, fontsize=5.5, handletextpad=.2)
    axA.spines["left"].set_visible(False)
    axA.tick_params(axis="y", length=0)
    axA.grid(axis="x", color=GRID, lw=.4)

    axB = fig.add_subplot(gs[0, 1])
    order = perf.sort_values("C_test_median", ascending=False)["cancer"].tolist()
    data = [sp.loc[sp.cancer == c, "C_test"].dropna().values for c in order]
    bp = axB.boxplot(data, patch_artist=True, widths=.6, showfliers=False,
                     medianprops=dict(color=INK, lw=.8), whiskerprops=dict(lw=.6, color=INK2),
                     capprops=dict(lw=.6, color=INK2), boxprops=dict(lw=.5, color=INK2))
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(C_FULL if c in rob else "#E6E6E6")
    axB.set_xticks(np.arange(1, len(order) + 1))
    axB.set_xticklabels(order, rotation=90)
    axB.axhline(.5, color=INK2, lw=.6, ls="--")
    axB.set_ylabel("Test-set C-index\n(100 stratified 70/30 splits)")
    axB.set_title("DRG LASSO-Cox models: out-of-sample performance", loc="left")
    axB.legend(handles=[Patch(color=C_FULL, label="Robust model"), Patch(color="#E6E6E6", label="Not robust")],
               loc="upper right", frameon=False, fontsize=5.5)
    axB.grid(axis="y", color=GRID, lw=.4)

    if FIG5_TABLE:
        c = gridspec.GridSpecFromSubplotSpec(1, 4, gs[1, :], width_ratios=[.55, .02, .2, 1.0], wspace=.08)
        axC, cC = fig.add_subplot(c[0]), fig.add_subplot(c[1])
    else:  # revision: performance summary moved to Table 2 of the manuscript
        c = gridspec.GridSpecFromSubplotSpec(1, 4, gs[1, :], width_ratios=[.42, .55, .02, .55], wspace=.08)
        axC, cC = fig.add_subplot(c[1]), fig.add_subplot(c[2])
    M = stab[stab.cancer.isin(rob)].pivot(index="gene", columns="cancer", values="selection_freq").reindex(GENES)[rob]
    heat(axC, M, None, cmap="Purples", sequential=True, vlim=1, cbar_label="selection frequency\n(100 splits)",
         cax=cC, xrot=0)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iat[i, j]
            if pd.notna(v) and v >= .6:
                axC.text(j + .5, i + .55, f"{v:.2f}", ha="center", va="center", fontsize=5,
                         color="white" if v > .75 else INK)
    for lab in axC.get_yticklabels():
        lab.set_color(ARM_COL.get(lab.get_text(), INK))
        lab.set_fontstyle("italic")
    axC.axhline(len(CORE), color=INK, lw=.8)
    axC.set_title("Gene selection stability (robust models)", loc="left")
    if FIG5_TABLE:
        axT = fig.add_subplot(c[3])
        axT.axis("off")
        pr = perf.set_index("cancer").loc[rob]
        rows = [[ca, f"{r.C_test_median:.2f} ({r.C_test_q025:.2f}–{r.C_test_q975:.2f})",
                 f"{r.AUC_1y_median:.2f} / {r.AUC_3y_median:.2f} / {r.AUC_5y_median:.2f}",
                 f"{r.logrank_oos:.1e}", _wrap((r.stable_genes or "—").replace(";", ", "))] for ca, r in pr.iterrows()]
        tb = axT.table(cellText=rows, colLabels=["Cancer", "Test C (95%)", "AUC 1/3/5 y", "OOS P",
                                                 "Stable genes (≥ 60%)"],
                       loc="center", cellLoc="left", colWidths=[.11, .22, .23, .12, .32])
        tb.auto_set_font_size(False)
        tb.set_fontsize(5.6)
        tb.scale(1, 1.9)
        for (r, cc), cell in tb.get_celld().items():
            cell.set_edgecolor(GRID)
            cell.set_linewidth(.4)
            if r == 0:
                cell.set_text_props(fontweight="bold", color=INK)
                cell.set_facecolor("#F2F2F2")
        axT.set_title("Performance of robust models", loc="left", fontsize=7.5)
    fig.canvas.draw()
    letter(axA, "A", dx=-.06)
    letter(axB, "B", dx=-.07)
    letter(axC, "C", dx=-.07)
    if FIG5_TABLE:
        letter(axT, "D", dx=-.015)
    save(fig, out, "Figure5")


def _wrap(t, n=22):
    out, line = [], ""
    for w in t.split(", "):
        if line and len(line) + len(w) + 2 > n:
            out.append(line + ",")
            line = w
        else:
            line = f"{line}, {w}" if line else w
    return "\n".join(out + [line])


def figure6(T, out):
    cl = T["clin_sub"]
    s = cl[cl.feature == "score"]
    vars_ = ["stage", "grade", "age_group", "gender", "race"]
    vlab = {"stage": "Pathological stage", "grade": "Histological grade", "age_group": "Age (≤60 vs >60)",
            "gender": "Sex", "race": "Race"}
    F = s.pivot(index="variable", columns="cancer", values="FDR").reindex(vars_)
    L = -np.log10(F)
    L.index = [vlab[v] for v in L.index]
    F.index = L.index
    mut = T["mut"].reindex(GENES)
    gain = T["gain"].reindex(GENES)
    loss = T["loss"].reindex(GENES)
    cols = [c for c in gain.columns if c in mut.columns]

    fig = plt.figure(figsize=(7.2, 9.6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[.32, 1, 1], hspace=.55)
    t = gridspec.GridSpecFromSubplotSpec(1, 2, gs[0], width_ratios=[1, .015], wspace=.02)
    axA, cA = fig.add_subplot(t[0]), fig.add_subplot(t[1])
    reds = LinearSegmentedColormap.from_list("r", ["#FFFFFF", "#F4A582", "#B2182B"])
    heat(axA, L, F, cmap=reds, sequential=True, vlim=6, cbar_label="−log10 FDR", cax=cA, star_size=4.5)
    axA.set_title("DRG score across clinical subgroups (UALCAN-style; Kruskal–Wallis / Mann–Whitney)", loc="left")
    m = gridspec.GridSpecFromSubplotSpec(1, 2, gs[1], width_ratios=[1, .015], wspace=.02)
    axB, cB = fig.add_subplot(m[0]), fig.add_subplot(m[1])
    seq = LinearSegmentedColormap.from_list("p", ["#FFFFFF", "#CDB7DF", C_FULL])
    heat(axB, mut[cols], None, cmap=seq, sequential=True, vlim=max(10, np.nanpercentile(mut[cols].values, 99)),
         cbar_label="% mutated\n(non-silent, MC3)", cax=cB)
    for i, gname in enumerate(GENES):
        for j, ca in enumerate(cols):
            v = mut.at[gname, ca]
            if pd.notna(v) and v >= 5:
                axB.text(j + .5, i + .55, f"{v:.0f}", ha="center", va="center", fontsize=4.5,
                         color="white" if v > 10 else INK)
    axB.set_title("Somatic mutation frequency", loc="left")
    c = gridspec.GridSpecFromSubplotSpec(1, 2, gs[2], width_ratios=[1, .015], wspace=.02)
    axC, cC = fig.add_subplot(c[0]), fig.add_subplot(c[1])
    net = gain[cols] - loss[cols]
    heat(axC, net, None, vlim=40, cbar_label="% gain − % loss\n(|log2 ratio| > 0.3)", cax=cC)
    axC.set_title("Copy-number alteration balance (net gain in red, net loss in blue)", loc="left")
    for ax in (axB, axC):
        for lab in ax.get_yticklabels():
            lab.set_color(ARM_COL.get(lab.get_text(), INK))
            lab.set_fontstyle("italic")
        ax.axhline(len(CORE), color=INK, lw=.8)
    fig.canvas.draw()
    letter(axA, "A", dx=-.05)
    letter(axB, "B", dx=-.05)
    letter(axC, "C", dx=-.05)
    save(fig, out, "Figure6")


def load(results: Path):
    t = results / "tables"
    r = lambda f, **k: pd.read_csv(t / f, **k)  # noqa: E731
    T = dict(tn=r("01_tumor_vs_normal.csv"), coex=r("10_DRG_coexpression_median_rho.csv", index_col=0),
             sc=r("04_score_cox_km.csv"), meta=r("04_score_meta_analysis.csv"),
             clin=r("02_tumour_clinical_with_scores.csv", index_col=0),
             gcox=r("03_gene_univariate_cox.csv"), sub=r("11_subscores_tumour_normal_and_survival.csv"),
             corarm=r("11_core_vs_cytoskeleton_correlation.csv", index_col=0),
             met=r("06_metabolic_correlations.csv.gz"), mv=r("07_multivariable_cox.csv"),
             perf=r("08_lasso_repeated_split_performance_OS.csv"), splits=r("08_lasso_all_splits_OS.csv.gz"),
             stab=r("08_lasso_coefficients_and_stability_OS.csv"), clin_sub=r("05_clinical_subgroups.csv"),
             mut=r("09_mutation_frequency_pct.csv", index_col=0), gain=r("09_cnv_gain_pct.csv", index_col=0),
             loss=r("09_cnv_loss_pct.csv", index_col=0))
    return T


def main(argv=None):
    a = argparse.ArgumentParser()
    a.add_argument("--results", required=True)
    args = a.parse_args(argv)
    res = Path(args.results)
    out = res / "figures_main"
    T = load(res)
    for f in (figure1, figure2, figure3, figure4, figure5, figure6):
        f(T, out)
        print("done:", f.__name__)


if __name__ == "__main__":
    main()
