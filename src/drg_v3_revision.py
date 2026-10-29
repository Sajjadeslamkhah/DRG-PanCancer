#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DRG pan-cancer study — v3 REVISION ANALYSES (response to IJMS reviewers)

Reads the outputs of drg_pancancer_pipeline.py (v2) from --base and writes
NEW results to <base>/results_v3/ without touching earlier results.

  R1  Null models with expression-matched random gene sets
        a) subset independence (core vs cytoskeleton) vs random 10/14 pairs
        b) within-subset coherence vs random sets
        c) prognostic signal of DRG / core / cytoskeleton scores vs random sets
           (number of FDR-significant cancers, pooled random-effects HR)
  R2  Metabolic specificity
        strict oxidative-PPP set, adjustment for purity + NRF2-target score,
        gene-level null for SLC7A11 (random expression-matched genes)
  R3  Tumour-microenvironment adjustment
        marker-based stromal (CAF), endothelial, leukocyte, smooth-muscle scores;
        subset Cox and correlations adjusted for purity + microenvironment
  R4  Glioma molecular adjustment
        LGG/GBM: IDH and 1p/19q status (PanCanAtlas subtype calls; fallback
        IDH1/IDH2 mutation from MC3) + grade; IDH-stratified HRs
  R5  Proportional hazards for the reported adjusted models; piecewise
        (time-varying) HRs where PH fails
  R6  LASSO-Cox benchmarks
        clinical-only vs DRG vs clinical+DRG in identical splits, IPCW
        (Uno-type) time-dependent AUC, null of random 24-gene LASSO models

Usage (Colab):
  python drg_v3_revision.py --base "/content/drive/MyDrive/DRG_PanCancer"
"""
from __future__ import annotations

import argparse
import json
import logging
import multiprocessing as mp
import os
import pickle
import re
import shutil
import zlib
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.duration.hazard_regression import PHReg

warnings.filterwarnings("ignore")
LOG = logging.getLogger("drg_v3")

# ── gene sets used only in v3 ─────────────────────────────────────────────
OXPPP_STRICT = ["G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT", "TALDO1"]
NRF2_TARGETS = ["NQO1", "HMOX1", "AKR1C1", "AKR1C2", "AKR1C3", "AKR1B10", "CYP4F11", "ABCC2",
                "ABCC3", "OSGIN1", "GCLM", "TXNRD1", "SRXN1", "GPX2", "CBR1"]
MICROENV = {
    "Stromal_CAF": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A3", "DCN", "LUM", "FAP",
                    "PDGFRB", "THY1", "FBN1", "SPARC", "POSTN", "CTSK", "MMP2"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "ESAM", "EMCN", "KDR", "TIE1", "TEK", "ROBO4"],
    "Leukocyte": ["PTPRC", "CD53", "LCP2", "LAPTM5", "CORO1A", "CD37", "CD48", "IKZF1", "ITGB2",
                  "CD2", "CD3E", "CD247"],
    "Smooth_muscle": ["ACTA2", "MYH11", "CNN1", "DES", "TAGLN", "ACTG2", "LMOD1", "SYNPO2"],
}
SUBTYPE_URLS = [
    "TCGASubtype.20170308.tsv.gz",
    "TCGASubtype.20170308.tsv",
]


# ════════════════════════════════════════════════════════════════════════════
# helpers
# ════════════════════════════════════════════════════════════════════════════
def setup_logging(log_file: Path):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S")
    LOG.setLevel(logging.INFO)
    LOG.handlers.clear()
    for h in (logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, mode="a")):
        h.setFormatter(fmt)
        LOG.addHandler(h)


class RankEngine:
    """Within-sample ranks of the whole transcriptome; ssGSEA (tau = 0.25) for any
    gene set in O(k * n). Identical to drg.ssgsea() before range normalisation."""

    def __init__(self, expr: pd.DataFrame, chunk=500):
        X = expr.values.astype(np.float32)
        X = np.where(np.isnan(X), -np.inf, X)
        self.R = np.empty(X.shape, dtype=np.float32)
        for s0 in range(0, X.shape[1], chunk):
            self.R[:, s0:s0 + chunk] = rankdata(X[:, s0:s0 + chunk], axis=0)
        del X
        self.W = self.R ** 0.25
        self.RW = self.R * self.W
        self.N = self.R.shape[0]
        self.idx = {g: i for i, g in enumerate(expr.index)}
        self.samples = expr.columns

    def score(self, genes):
        ii = np.array([self.idx[g] for g in dict.fromkeys(genes) if g in self.idx])
        if len(ii) < 3:
            return None
        N, k = self.N, len(ii)
        s = self.RW[ii].sum(0) / self.W[ii].sum(0) - (N * (N + 1) / 2 - self.R[ii].sum(0)) / (N - k)
        return pd.Series(s.astype(float), index=self.samples)


class MatchedSampler:
    """Draw random gene sets matched to a target set by mean-expression decile."""

    def __init__(self, expr: pd.DataFrame, exclude, n_bins=20, seed=2026):
        m = expr.mean(axis=1, skipna=True)
        sd = expr.std(axis=1, skipna=True)
        ok = (sd > 0) & m.notna()
        m = m[ok]
        self.bin = pd.Series(pd.qcut(m.rank(method="first"), n_bins, labels=False), index=m.index)
        ex = set(exclude)
        self.pool = {b: np.array([g for g in grp.index if g not in ex])
                     for b, grp in self.bin.groupby(self.bin)}
        self.rng = np.random.default_rng(seed)

    def draw(self, genes):
        out = []
        for g in genes:
            b = self.bin.get(g, self.rng.integers(len(self.pool)))
            pool = self.pool[b]
            for _ in range(20):
                c = pool[self.rng.integers(len(pool))]
                if c not in out:
                    break
            out.append(c)
        return out


def cox_multi(t, e, X, iters=12):
    """Univariate Cox (Breslow) fitted simultaneously for every column of X.
    Returns beta, se, p (Wald)."""
    t, e, X = np.asarray(t, float), np.asarray(e, bool), np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    o = np.argsort(t, kind="mergesort")
    t, e, X = t[o], e[o], X[o]
    first = np.searchsorted(t, t, side="left")
    b = np.zeros(X.shape[1])

    def stats_(b):
        eta = X * b
        eta = eta - eta.max(0)
        w = np.exp(eta)
        rc = lambda a: np.cumsum(a[::-1], 0)[::-1][first]  # noqa: E731
        S0, S1, S2 = rc(w), rc(w * X), rc(w * X * X)
        m = S1 / S0
        U = (X[e] - m[e]).sum(0)
        I = (S2[e] / S0[e] - m[e] ** 2).sum(0)
        return U, I

    for _ in range(iters):
        U, I = stats_(b)
        b = b + np.clip(U / np.maximum(I, 1e-12), -1, 1)
    U, I = stats_(b)
    se = 1 / np.sqrt(np.maximum(I, 1e-12))
    return b, se, 2 * stats.norm.sf(np.abs(b / se))


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        q = p[ok]
        n = len(q)
        o = np.argsort(q)
        r = q[o] * n / np.arange(1, n + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        res = np.empty(n)
        res[o] = np.minimum(r, 1)
        out[ok] = res
    return out


def dl_meta(y, se):
    y, se = np.asarray(y), np.asarray(se)
    w = 1 / se ** 2
    fe = (w * y).sum() / w.sum()
    Q = (w * (y - fe) ** 2).sum()
    k = len(y)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (se ** 2 + tau2)
    mu = (wr * y).sum() / wr.sum()
    s = np.sqrt(1 / wr.sum())
    return mu, s, 2 * stats.norm.sf(abs(mu / s)), tau2


def partial_rank_corr(x, Y, Z=None):
    """Spearman partial correlation of vector x with each column of Y, controlling
    for columns of Z (rank residuals). Returns rho, p arrays."""
    x = rankdata(np.asarray(x, float))
    Y = np.apply_along_axis(rankdata, 0, np.asarray(Y, float).reshape(len(x), -1))
    n = len(x)
    if Z is None or (hasattr(Z, "shape") and Z.size == 0):
        D = np.ones((n, 1))
    else:
        Zr = np.apply_along_axis(rankdata, 0, np.asarray(Z, float).reshape(n, -1))
        D = np.c_[np.ones(n), Zr]
    H = D @ np.linalg.pinv(D)
    rx = x - H @ x
    RY = Y - H @ Y
    r = (rx @ RY) / np.sqrt((rx @ rx) * (RY * RY).sum(0))
    df = n - D.shape[1] - 1
    tt = r * np.sqrt(df / np.maximum(1e-12, 1 - r * r))
    return r, 2 * stats.t.sf(np.abs(tt), df)


def km_surv(times, events):
    """Return step function S(t) (right-continuous) of Kaplan-Meier."""
    times = np.asarray(times, float)
    events = np.asarray(events, int)
    ut = np.unique(times[events == 1])
    s = 1.0
    xs, ys = [0.0], [1.0]
    for u in ut:
        n_risk = (times >= u).sum()
        d = ((times == u) & (events == 1)).sum()
        s *= 1 - d / n_risk
        xs.append(u)
        ys.append(s)
    xs, ys = np.array(xs), np.array(ys)
    return lambda t: ys[np.searchsorted(xs, t, side="right") - 1]


def auc_ipcw(t, e, r, t0, G):
    """Cumulative/dynamic AUC(t0) with inverse-probability-of-censoring weights
    (Uno et al. 2007); G = KM survival function of the censoring distribution."""
    t, e, r = np.asarray(t, float), np.asarray(e, int), np.asarray(r, float)
    cases = (t <= t0) & (e == 1)
    ctrls = t > t0
    if cases.sum() < 3 or ctrls.sum() < 3:
        return np.nan
    w = 1 / np.maximum(G(t[cases] - 1e-9), 1e-3)
    rc, rn = r[cases][:, None], r[ctrls][None, :]
    conc = (rc > rn).sum(1) + 0.5 * (rc == rn).sum(1)
    return float((w * conc).sum() / (w.sum() * ctrls.sum()))


def ph_test(fit, t, e, names):
    """Schoenfeld-residual correlation with ranked event time (per covariate)."""
    sr = fit.schoenfeld_residuals
    out = {}
    tt = np.asarray(t, float)
    ok = ~np.isnan(sr[:, 0])
    rt = rankdata(tt[ok])
    for j, nm in enumerate(names):
        out[nm] = stats.pearsonr(rt, sr[ok, j])[1]
    return out


def fit_cox(d, cols, tcol="OS.time", ecol="OS"):
    return PHReg(d[tcol].values.astype(float), d[cols].values.astype(float),
                 status=d[ecol].values.astype(int), ties="efron").fit(disp=False)


def piecewise_cox(d, var, covs, cut, tcol="OS.time", ecol="OS"):
    """Episode-split Cox with separate HRs for var before/after `cut` (years)."""
    rows = []
    for _, r in d.iterrows():
        T, E = r[tcol], int(r[ecol])
        base = {c: r[c] for c in covs}
        rows.append(dict(base, entry=0.0, time=min(T, cut), event=int(E and T <= cut),
                         early=r[var], late=0.0))
        if T > cut:
            rows.append(dict(base, entry=cut, time=T, event=E, early=0.0, late=r[var]))
    ep = pd.DataFrame(rows)
    ep = ep[ep["time"] > ep["entry"]]
    X = ep[["early", "late"] + covs].values.astype(float)
    m = PHReg(ep["time"].values, X, status=ep["event"].values, entry=ep["entry"].values,
              ties="efron").fit(disp=False)
    ci = np.exp(m.conf_int())
    return dict(HR_early=np.exp(m.params[0]), lo_early=ci[0, 0], hi_early=ci[0, 1], p_early=m.pvalues[0],
                HR_late=np.exp(m.params[1]), lo_late=ci[1, 0], hi_late=ci[1, 1], p_late=m.pvalues[1])


def savefig(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def eligible(clin, ep="OS", min_n=30, min_ev=10):
    out = {}
    for ca, sub in clin.groupby("cancer"):
        d = sub[[ep, f"{ep}.time"]].dropna()
        d = d[d[f"{ep}.time"] > 0]
        if len(d) >= min_n and d[ep].sum() >= min_ev:
            out[ca] = d.index
    return out


# ════════════════════════════════════════════════════════════════════════════
# R1  null models
# ════════════════════════════════════════════════════════════════════════════
def r1_null_models(E, S, expr, clin, core, cyto, cfg, T, F):
    LOG.info(f"[R1] Null models with {cfg.n_null} expression-matched random gene sets")
    tum = clin.index
    groups = {ca: tum[clin.cancer == ca] for ca in clin.cancer.unique() if (clin.cancer == ca).sum() >= 30}

    def med_pair_rho(a, b):
        return np.median([stats.spearmanr(a[ix], b[ix]).correlation for ix in groups.values()])

    sc_core, sc_cyto = E.score(core), E.score(cyto)
    obs_pair = med_pair_rho(sc_core, sc_cyto)

    # coherence: mean pairwise within-tumour Spearman across cancers (median over cancers)
    def coherence(genes):
        vals = []
        for ix in groups.values():
            m = expr.loc[genes, ix].T.rank().corr().values
            vals.append(np.nanmean(m[np.triu_indices(len(genes), 1)]))
        return np.median(vals)

    obs_coh_core, obs_coh_cyto = coherence(core), coherence(cyto)
    null_pair, null_coh_core, null_coh_cyto = [], [], []
    for i in range(cfg.n_null):
        rc, rk = S.draw(core), S.draw(cyto)
        null_pair.append(med_pair_rho(E.score(rc), E.score(rk)))
        if i < cfg.n_null_coherence:
            null_coh_core.append(coherence(rc))
            null_coh_cyto.append(coherence(rk))
        if (i + 1) % 200 == 0:
            LOG.info(f"  structure null {i + 1}/{cfg.n_null}")
    null_pair = np.array(null_pair)
    res_struct = pd.DataFrame([
        dict(test="core-vs-cytoskeleton median rho", observed=obs_pair, null_median=np.median(null_pair),
             null_q025=np.quantile(null_pair, .025), null_q975=np.quantile(null_pair, .975),
             empirical_p_lower=(1 + (null_pair <= obs_pair).sum()) / (1 + len(null_pair)), n_null=len(null_pair)),
        dict(test="core coherence (median pairwise rho)", observed=obs_coh_core, null_median=np.median(null_coh_core),
             null_q025=np.quantile(null_coh_core, .025), null_q975=np.quantile(null_coh_core, .975),
             empirical_p_upper=(1 + (np.array(null_coh_core) >= obs_coh_core).sum()) / (1 + len(null_coh_core)),
             n_null=len(null_coh_core)),
        dict(test="cytoskeleton coherence (median pairwise rho)", observed=obs_coh_cyto,
             null_median=np.median(null_coh_cyto), null_q025=np.quantile(null_coh_cyto, .025),
             null_q975=np.quantile(null_coh_cyto, .975),
             empirical_p_upper=(1 + (np.array(null_coh_cyto) >= obs_coh_cyto).sum()) / (1 + len(null_coh_cyto)),
             n_null=len(null_coh_cyto)),
    ])
    res_struct.to_csv(T / "R1a_structure_null.csv", index=False)
    LOG.info(f"  core-vs-cyto rho observed {obs_pair:.3f}; null median {np.median(null_pair):.3f}")

    # prognostic null (OS)
    el = eligible(clin)
    sets = {"DRG_24": core + cyto, "core_10": core, "cytoskeleton_14": cyto}
    prog_rows, null_store = [], {}
    for name, genes in sets.items():
        obs = E.score(genes)
        mats = [obs.values]
        for _ in range(cfg.n_null):
            mats.append(E.score(S.draw(genes)).values)
        M = np.vstack(mats).T  # samples x (1+n_null)
        M = pd.DataFrame(M, index=E.samples)
        betas, ses = [], []
        for ca, ix in el.items():
            d = clin.loc[ix]
            Xc = M.loc[ix].values
            Xc = (Xc - Xc.mean(0)) / Xc.std(0)
            b, se, p = cox_multi(d["OS.time"].values, d["OS"].values, Xc)
            betas.append(b)
            ses.append(se)
        B, SE = np.array(betas), np.array(ses)  # cancers x sets
        P = 2 * stats.norm.sf(np.abs(B / SE))
        Q = np.apply_along_axis(bh, 0, P)
        nsig = (Q < 0.05).sum(0)
        meta = np.array([dl_meta(B[:, j], SE[:, j]) for j in range(B.shape[1])])
        mean_abs = np.abs(B).mean(0)
        null_store[name] = dict(nsig=nsig[1:], obs_nsig=nsig[0], meta=meta[1:, 0], obs_meta=meta[0, 0],
                                absb=mean_abs[1:], obs_absb=mean_abs[0])
        prog_rows.append(dict(
            score=name, cancers=len(el), obs_n_FDR_sig=int(nsig[0]),
            null_median_n_FDR_sig=float(np.median(nsig[1:])), null_q95_n_FDR_sig=float(np.quantile(nsig[1:], .95)),
            p_n_FDR_sig=(1 + (nsig[1:] >= nsig[0]).sum()) / (1 + cfg.n_null),
            obs_pooled_HR=float(np.exp(meta[0, 0])), null_median_pooled_HR=float(np.exp(np.median(meta[1:, 0]))),
            p_pooled_logHR_abs=(1 + (np.abs(meta[1:, 0]) >= abs(meta[0, 0])).sum()) / (1 + cfg.n_null),
            obs_mean_abs_logHR=float(mean_abs[0]), null_median_mean_abs_logHR=float(np.median(mean_abs[1:])),
            p_mean_abs_logHR=(1 + (mean_abs[1:] >= mean_abs[0]).sum()) / (1 + cfg.n_null)))
        LOG.info(f"  {name}: FDR-sig cancers obs {nsig[0]} vs null median {np.median(nsig[1:]):.0f} "
                 f"(P = {prog_rows[-1]['p_n_FDR_sig']:.3f})")
    prog = pd.DataFrame(prog_rows)
    prog.to_csv(T / "R1b_prognostic_null_OS.csv", index=False)
    return res_struct, prog, dict(null_pair=null_pair, obs_pair=obs_pair, prog=null_store)


# ════════════════════════════════════════════════════════════════════════════
# R2  metabolic specificity
# ════════════════════════════════════════════════════════════════════════════
def r2_metabolic(E, S2, expr, clin, core, cyto, ppp_set, oxphos_set, micro_scores, cfg, T):
    LOG.info("[R2] Metabolic specificity (strict oxPPP, NRF2 adjustment, gene-level null)")
    tum = clin.index
    scores = {
        "PPP_NADPH_full": E.score([g for g in ppp_set if g not in core + cyto]),
        "oxPPP_strict": E.score(OXPPP_STRICT),
        "OXPHOS": E.score([g for g in oxphos_set if g not in core + cyto]),
        "NRF2_targets": E.score([g for g in NRF2_TARGETS if g not in OXPPP_STRICT]),
    }
    feats = {"DRG_score": E.score(core + cyto), "core_subset": E.score(core),
             "cytoskeleton_subset": E.score(cyto)}
    for g in ["SLC7A11", "SLC3A2", "RPN1", "NCKAP1", "LRPPRC", "NDUFS1", "NDUFA11"]:
        if g in expr.index:
            feats[g] = expr.loc[g, tum]
    rows = []
    for ca, sub in clin.groupby("cancer"):
        d = sub.dropna(subset=["purity"])
        if len(d) < 30:
            continue
        ix = d.index
        Z1 = d[["purity"]].values
        Z2 = np.c_[d["purity"].values, scores["NRF2_targets"][ix].values]
        Z3 = np.c_[Z2, np.column_stack([micro_scores[k][ix].values for k in micro_scores])] if micro_scores else Z2
        for fn, fv in feats.items():
            x = fv[ix].values
            for pw in ["PPP_NADPH_full", "oxPPP_strict", "OXPHOS"]:
                y = scores[pw][ix].values[:, None]
                for adj, Z in (("purity", Z1), ("purity+NRF2", Z2), ("purity+NRF2+microenv", Z3)):
                    r, p = partial_rank_corr(x, y, Z)
                    rows.append(dict(cancer=ca, feature=fn, pathway=pw, adjustment=adj, rho=r[0], p=p[0], n=len(ix)))
    res = pd.DataFrame(rows)
    res["FDR"] = res.groupby(["pathway", "adjustment"])["p"].transform(lambda x: bh(x.values))
    res.to_csv(T / "R2a_metabolic_specificity_partial.csv.gz", index=False)
    summ = res.groupby(["feature", "pathway", "adjustment"]).apply(
        lambda g: pd.Series(dict(median_rho=g.rho.median(), n_pos=((g.FDR < .05) & (g.rho > 0)).sum(),
                                 n_neg=((g.FDR < .05) & (g.rho < 0)).sum(), k=len(g)))).reset_index()
    summ.to_csv(T / "R2b_metabolic_specificity_summary.csv", index=False)

    # gene-level null: SLC7A11 vs expression-matched random genes, strict oxPPP | purity + NRF2
    LOG.info("  gene-level null for SLC7A11 ~ strict oxPPP")
    obs_list, null_mat = [], []
    rnd_genes = [S2.draw(["SLC7A11"])[0] for _ in range(cfg.n_null)]
    rnd_genes = [g for g in rnd_genes if g in expr.index]
    for ca, sub in clin.groupby("cancer"):
        d = sub.dropna(subset=["purity"])
        if len(d) < 30:
            continue
        ix = d.index
        Z2 = np.c_[d["purity"].values, scores["NRF2_targets"][ix].values]
        y = scores["oxPPP_strict"][ix].values
        Y = np.c_[expr.loc[["SLC7A11"] + rnd_genes, ix].T.values]
        Yf = np.nan_to_num(Y, nan=np.nanmin(Y))
        r, _ = partial_rank_corr(y, Yf, Z2)
        obs_list.append(r[0])
        null_mat.append(r[1:])
    null_med = np.median(np.array(null_mat), axis=0)
    obs_med = float(np.median(obs_list))
    gnull = pd.DataFrame([dict(gene="SLC7A11", pathway="oxPPP_strict", adjustment="purity+NRF2",
                               observed_median_rho=obs_med, null_median=float(np.median(null_med)),
                               null_q95=float(np.quantile(null_med, .95)),
                               empirical_p=(1 + (null_med >= obs_med).sum()) / (1 + len(null_med)),
                               n_null=len(null_med))])
    gnull.to_csv(T / "R2c_SLC7A11_gene_null.csv", index=False)
    LOG.info(f"  SLC7A11~oxPPP|purity+NRF2 median rho {obs_med:.3f}; null median {np.median(null_med):.3f}")
    return summ, gnull, dict(null_med=null_med, obs_med=obs_med)


# ════════════════════════════════════════════════════════════════════════════
# R3  microenvironment adjustment
# ════════════════════════════════════════════════════════════════════════════
def r3_microenv(E, clin, core, cyto, micro_scores, oxphos_score, cfg, T):
    LOG.info("[R3] Microenvironment-adjusted subset analyses")
    cl = clin.copy()
    cl["core_z_raw"] = E.score(core)[cl.index]
    cl["cyto_z_raw"] = E.score(cyto)[cl.index]
    for k, s in micro_scores.items():
        cl[k] = s[cl.index]
    mcols = list(micro_scores)
    rows, crow = [], []
    for ca, sub in cl.groupby("cancer"):
        d = sub.dropna(subset=["OS", "OS.time", "purity"] + mcols)
        d = d[d["OS.time"] > 0].copy()
        if len(d) < 30 or d.OS.sum() < 10:
            continue
        for c in ["core_z_raw", "cyto_z_raw", "purity"] + mcols:
            d[c + "_z"] = (d[c] - d[c].mean()) / d[c].std()
        adj = ["purity_z"] + [m + "_z" for m in mcols]
        for model, cols in (("cyto+purity+microenv", ["cyto_z_raw_z"] + adj),
                            ("core+purity+microenv", ["core_z_raw_z"] + adj),
                            ("joint+purity+microenv", ["core_z_raw_z", "cyto_z_raw_z"] + adj)):
            try:
                m = fit_cox(d, cols)
            except Exception:  # noqa: BLE001
                continue
            ci = np.exp(m.conf_int())
            for j, cnm in enumerate(cols[:2] if model.startswith("joint") else cols[:1]):
                rows.append(dict(cancer=ca, model=model, term=cnm.replace("_z_raw_z", ""), n=len(d),
                                 events=int(d.OS.sum()), HR=np.exp(m.params[j]), lo=ci[j, 0], hi=ci[j, 1],
                                 p=m.pvalues[j]))
        # partial correlations adjusted for purity + microenvironment
        Z = d[["purity"] + mcols].values
        r1, p1 = partial_rank_corr(d["core_z_raw"].values, d["cyto_z_raw"].values[:, None], Z)
        r2, p2 = partial_rank_corr(d["cyto_z_raw"].values, oxphos_score[d.index].values[:, None], Z)
        r3, _ = partial_rank_corr(d["cyto_z_raw"].values, d["purity"].values[:, None])
        crow.append(dict(cancer=ca, n=len(d), rho_core_cyto_adj=r1[0], p_core_cyto_adj=p1[0],
                         rho_cyto_OXPHOS_adj=r2[0], p_cyto_OXPHOS_adj=p2[0], rho_cyto_purity=r3[0],
                         **{f"rho_cyto_{m}": stats.spearmanr(d["cyto_z_raw"], d[m]).correlation for m in mcols}))
    res = pd.DataFrame(rows)
    if len(res):
        res["FDR"] = res.groupby(["model", "term"])["p"].transform(lambda x: bh(x.values))
    res.to_csv(T / "R3a_subset_cox_purity_microenv.csv", index=False)
    cr = pd.DataFrame(crow)
    if len(cr):
        cr["FDR_core_cyto"] = bh(cr.p_core_cyto_adj)
        cr["FDR_cyto_OXPHOS"] = bh(cr.p_cyto_OXPHOS_adj)
    cr.to_csv(T / "R3b_correlations_purity_microenv.csv", index=False)
    return res, cr


# ════════════════════════════════════════════════════════════════════════════
# R4  glioma molecular adjustment
# ════════════════════════════════════════════════════════════════════════════
def load_glioma_molecular(raw: Path, clin, drg, cfg):
    sub = None
    for name in SUBTYPE_URLS:
        dest = raw / name.replace(".gz", "") if not name.endswith(".gz") else raw / name
        try:
            if not dest.exists():
                drg.download(drg.HUB + name, dest)
            df = pd.read_csv(dest, sep="\t", low_memory=False)
            sub = df
            LOG.info(f"  subtype file loaded: {name} ({df.shape})")
            break
        except Exception as ex:  # noqa: BLE001
            LOG.warning(f"  subtype file {name} unavailable: {ex}")
    mol = pd.DataFrame(index=clin.index)
    if sub is not None:
        idcol = sub.columns[0]
        sub[idcol] = sub[idcol].astype(str).str[:15]
        s = sub.set_index(idcol)
        s = s[~s.index.duplicated()]
        lab = pd.Series(np.nan, index=s.index, dtype=object)
        # (a) explicit IDH/codel labels, if present in any column
        for c in s.columns:
            v = s[c].astype(str)
            ext = v.str.extract(r"(IDHmut-codel|IDHmut-non-codel|IDHwt)", expand=False)
            lab = lab.where(lab.notna(), ext)
        # (b) pan-glioma DNA-methylation clusters (Ceccarelli et al., Cell 2016):
        #     Codel -> IDHmut-codel; G-CIMP-high/-low -> IDHmut-non-codel;
        #     Classic-like, Mesenchymal-like, LGm6-GBM, PA-like -> IDHwt
        meth = {"Codel": "IDHmut-codel", "G-CIMP-high": "IDHmut-non-codel", "G-CIMP-low": "IDHmut-non-codel",
                "Classic-like": "IDHwt", "Mesenchymal-like": "IDHwt", "LGm6-GBM": "IDHwt", "PA-like": "IDHwt"}
        pat = r"GBM_LGG\.(" + "|".join(re.escape(k) for k in meth) + r")\b"
        for c in s.columns:
            ext = s[c].astype(str).str.extract(pat, expand=False).map(meth)
            lab = lab.where(lab.notna(), ext)
        LOG.info(f"  glioma molecular labels: {lab.value_counts(dropna=True).to_dict()}")
        mol["glioma_subtype"] = lab.reindex(mol.index)
        mol["source"] = np.where(mol["glioma_subtype"].notna(), "PanCanAtlas subtype", None)
    else:
        mol["glioma_subtype"] = np.nan
        mol["source"] = None
    # fallback: IDH1/IDH2 non-silent mutation (MC3)
    mc3 = raw / "mc3_nonsilent_gene.xena.gz"
    if mc3.exists():
        try:
            m = drg.load_xena_matrix(mc3, ["IDH1", "IDH2"], list(clin.index))
            idh = (m.fillna(0).sum(0) > 0)
            mol["IDH_mut_MC3"] = idh.reindex(mol.index)
        except Exception as ex:  # noqa: BLE001
            LOG.warning(f"  MC3 IDH fallback failed: {ex}")
    return mol


def r4_glioma(E, clin, core, cyto, mol, drg, cfg, T):
    LOG.info("[R4] Glioma molecular adjustment (IDH, 1p/19q, grade)")
    rows, strat = [], []
    cl = clin.join(mol)
    cl["DRG"] = E.score(core + cyto)[cl.index]
    cl["cyto"] = E.score(cyto)[cl.index]
    cl["core"] = E.score(core)[cl.index]
    for ca in ["LGG", "GBM"]:
        d = cl[(cl.cancer == ca)].dropna(subset=["OS", "OS.time", "age"]).copy()
        d = d[d["OS.time"] > 0]
        if len(d) < 30:
            continue
        if d["glioma_subtype"].notna().mean() > 0.5:
            d["IDHmut"] = d["glioma_subtype"].str.startswith("IDHmut").astype(float).where(d["glioma_subtype"].notna())
            d["codel"] = (d["glioma_subtype"] == "IDHmut-codel").astype(float).where(d["glioma_subtype"].notna())
            src = "PanCanAtlas subtype"
        elif "IDH_mut_MC3" in d and d["IDH_mut_MC3"].notna().mean() > 0.5:
            d["IDHmut"] = d["IDH_mut_MC3"].astype(float)
            d["codel"] = np.nan
            src = "MC3 IDH1/IDH2 mutation (no 1p/19q)"
        else:
            LOG.warning(f"  {ca}: no molecular data, skipped")
            continue
        d["male"] = (d["gender"] == "MALE").astype(float)
        d["grade3"] = (d["grade"] == "G3").astype(float).where(d["grade"].isin(["G2", "G3"]))
        d["age10"] = d["age"] / 10
        for sc in ["DRG", "cyto", "core"]:
            d[sc + "_z"] = (d[sc] - d[sc].mean()) / d[sc].std()
        covsets = {"age+sex+IDH": ["age10", "male", "IDHmut"],
                   "age+sex+IDH+codel": ["age10", "male", "IDHmut", "codel"],
                   "age+sex+IDH+codel+grade": ["age10", "male", "IDHmut", "codel", "grade3"],
                   "age+sex+IDH+codel+grade+purity": ["age10", "male", "IDHmut", "codel", "grade3", "purity"]}
        for sc in ["DRG_z", "cyto_z", "core_z"]:
            for name, cv in covsets.items():
                cv = [c for c in cv if d[c].notna().mean() > 0.5 and d[c].nunique() > 1]
                dd = d.dropna(subset=[sc] + cv)
                if len(dd) < 30 or dd.OS.sum() < 10:
                    continue
                try:
                    m = fit_cox(dd, [sc] + cv)
                except Exception:  # noqa: BLE001
                    continue
                ci = np.exp(m.conf_int())
                rows.append(dict(cancer=ca, score=sc[:-2], adjustment=name, covariates_used="+".join(cv),
                                 source=src, n=len(dd), events=int(dd.OS.sum()), HR=np.exp(m.params[0]),
                                 lo=ci[0, 0], hi=ci[0, 1], p=m.pvalues[0]))
            # IDH-stratified
            for grp, gd in d.dropna(subset=["IDHmut"]).groupby("IDHmut"):
                gd = gd.dropna(subset=[sc])
                if len(gd) < 25 or gd.OS.sum() < 8:
                    continue
                try:
                    m = fit_cox(gd, [sc, "age10"])
                except Exception:  # noqa: BLE001
                    continue
                ci = np.exp(m.conf_int())
                strat.append(dict(cancer=ca, score=sc[:-2], stratum="IDH-mutant" if grp == 1 else "IDH-wildtype",
                                  n=len(gd), events=int(gd.OS.sum()), HR=np.exp(m.params[0]), lo=ci[0, 0],
                                  hi=ci[0, 1], p=m.pvalues[0]))
        # score difference by molecular group
        if d["glioma_subtype"].notna().any():
            grp = d.dropna(subset=["glioma_subtype"]).groupby("glioma_subtype")["DRG"].median()
            LOG.info(f"  {ca} DRG score by subtype: {grp.round(3).to_dict()}")
    res, st = pd.DataFrame(rows), pd.DataFrame(strat)
    res.to_csv(T / "R4a_glioma_molecular_adjusted_cox.csv", index=False)
    st.to_csv(T / "R4b_glioma_IDH_stratified_cox.csv", index=False)
    return res, st


# ════════════════════════════════════════════════════════════════════════════
# R5  proportional hazards for reported adjusted models
# ════════════════════════════════════════════════════════════════════════════
def r5_ph(E, clin, core, cyto, cfg, T):
    LOG.info("[R5] Proportional-hazards checks on adjusted models + piecewise HRs")
    cl = clin.copy()
    cl["DRG"] = E.score(core + cyto)[cl.index]
    cl["core"] = E.score(core)[cl.index]
    cl["cyto"] = E.score(cyto)[cl.index]
    rows, pw = [], []
    for ca, sub in cl.groupby("cancer"):
        d = sub.dropna(subset=["OS", "OS.time", "age"]).copy()
        d = d[d["OS.time"] > 0]
        if len(d) < 30 or d.OS.sum() < 10:
            continue
        d["age10"] = d["age"] / 10
        d["male"] = (d["gender"] == "MALE").astype(float)
        stg = d["stage"].map({"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
        cov = ["age10"] + (["male"] if d["male"].nunique() > 1 else [])
        if stg.notna().mean() > 0.6 and stg.nunique() > 1:
            d["stage_num"] = stg
            cov.append("stage_num")
        for sc in ["DRG", "core", "cyto"]:
            d[sc + "_z"] = (d[sc] - d[sc].mean()) / d[sc].std()
        models = {"DRG+clinical+purity": ["DRG_z"] + cov + ["purity"],
                  "joint core+cyto": ["core_z", "cyto_z"]}
        for mname, cols in models.items():
            dd = d.dropna(subset=cols)
            if len(dd) < 30 or dd.OS.sum() < 10:
                continue
            try:
                m = fit_cox(dd, cols)
                ph = ph_test(m, dd["OS.time"].values, dd["OS"].values, cols)
            except Exception:  # noqa: BLE001
                continue
            for c in cols[:2] if mname.startswith("joint") else cols[:1]:
                rows.append(dict(cancer=ca, model=mname, term=c[:-2], PH_p=ph[c], PH_p_min_any_covariate=min(ph.values())))
            # piecewise where PH violated for the score
            for c in cols[:2] if mname.startswith("joint") else cols[:1]:
                if ph[c] < 0.05:
                    cut = float(np.median(dd.loc[dd.OS == 1, "OS.time"]))
                    try:
                        r = piecewise_cox(dd, c, [x for x in cols if x != c], cut)
                        pw.append(dict(cancer=ca, model=mname, term=c[:-2], cut_years=cut, **r))
                    except Exception as ex:  # noqa: BLE001
                        LOG.warning(f"  piecewise {ca} {c}: {ex}")
    res = pd.DataFrame(rows)
    if len(res):
        res["PH_q"] = res.groupby(["model", "term"])["PH_p"].transform(lambda x: bh(x.values))
    res.to_csv(T / "R5a_PH_tests_adjusted_models.csv", index=False)
    pd.DataFrame(pw).to_csv(T / "R5b_piecewise_HR_PH_violations.csv", index=False)
    return res, pd.DataFrame(pw)


# ════════════════════════════════════════════════════════════════════════════
# R6  LASSO benchmarks
# ════════════════════════════════════════════════════════════════════════════
_W: dict = {}  # per-cancer data shared with forked workers


def _r6_task(task):
    """One LASSO-Cox job: columns of _W['M'], CV seed, split ids -> (Cs, risks|None)."""
    cols, seed, split_ids, want_risk = task
    W = _W
    drg, lc, t, e = W["drg"], W["lc"], W["t"], W["e"]
    Xall = W["M"][:, cols]
    rng = np.random.default_rng(seed)
    Cs, risks = [], []
    for k in split_ids:
        tr, te = W["splits"][k]
        mu, sd = Xall[tr].mean(0), Xall[tr].std(0)
        sd[sd == 0] = 1
        X = (Xall - mu) / sd
        b, _ = drg._lasso_fit_cv(lc, t[tr], e[tr], X[tr], rng, 3, 15, "1se")
        r = X @ b
        Cs.append(drg.cindex(t[te], e[te], r[te]) if np.any(b) else 0.5)
        risks.append(r if want_risk else None)
    return Cs, risks


def _r6_map(tasks, n_jobs):
    if n_jobs <= 1 or len(tasks) == 1:
        return [_r6_task(x) for x in tasks]
    with mp.get_context("fork").Pool(n_jobs) as pool:
        return pool.map(_r6_task, tasks, chunksize=1)


def r6_lasso(expr, clin, genes, S, drg, resdir, cfg, T, ckdir=None):
    LOG.info(f"[R6] LASSO benchmarks ({cfg.n_splits_bench} splits; random-set null in robust cancers)")
    perf_path = resdir / "tables" / "08_lasso_repeated_split_performance_OS.csv"
    robust = []
    if perf_path.exists():
        pf = pd.read_csv(perf_path)
        if "robust" in pf:
            robust = pf.loc[pf.robust.astype(bool), "cancer"].tolist()
    n_jobs = cfg.n_jobs if cfg.n_jobs > 0 else (os.cpu_count() or 1)
    LOG.info(f"  robust cancers: {robust}; parallel workers: {n_jobs}"
             + (f"; checkpoints: {ckdir}" if ckdir else ""))
    if ckdir:
        ckdir.mkdir(parents=True, exist_ok=True)
    lc = drg.LassoCox(max_iter=1500, tol=1e-6)
    rows, null_rows = [], []
    gpos = {g: i for i, g in enumerate(expr.index)}

    eligible = []
    for ca, sub in clin.groupby("cancer"):
        d = sub.dropna(subset=["OS", "OS.time", "age"])
        d = d[d["OS.time"] > 0]
        if len(d) >= 60 and d.OS.sum() >= 30:
            eligible.append(ca)
    t_start = time.time()
    for k_ca, ca in enumerate(eligible, 1):
        ck = ckdir / f"R6_{ca}.pkl" if ckdir else None
        if ck is not None and ck.exists():
            saved = pickle.loads(ck.read_bytes())
            rows += saved["rows"]
            null_rows += saved["null_rows"]
            if saved.get("nulls") is not None:
                np.save(T / f"R6_null_C_{ca}.npy", saved["nulls"])
            LOG.info(f"  R6 [{k_ca}/{len(eligible)}] {ca}: loaded from checkpoint")
            continue
        t_ca = time.time()
        crc = zlib.crc32(ca.encode())
        rng = np.random.default_rng([cfg.seed, crc])
        S.rng = np.random.default_rng([cfg.seed + 7, crc])
        d = clin[clin.cancer == ca].dropna(subset=["OS", "OS.time", "age"]).copy()
        d = d[d["OS.time"] > 0]
        d["age10"] = d["age"] / 10
        d["male"] = (d["gender"] == "MALE").astype(float)
        clin_cols = ["age10"] + (["male"] if d["male"].nunique() > 1 else [])
        stg = d["stage"].map({"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
        if stg.notna().mean() > 0.6 and stg.nunique() > 1:
            d["stage_num"] = stg
            clin_cols.append("stage_num")
        d = d.dropna(subset=clin_cols)
        t, e = d["OS.time"].values.astype(float), d["OS"].values.astype(int)
        M = expr[d.index].T.values.astype(float)
        M = np.where(np.isnan(M), np.nanmean(M, axis=0, keepdims=True), M)
        M = np.nan_to_num(M)
        idx = np.arange(len(e))
        splits = []
        for _ in range(cfg.n_splits_bench):
            tr = np.concatenate([rng.choice(idx[e == c], int(round(0.7 * (e == c).sum())), replace=False)
                                 for c in (0, 1)])
            splits.append((np.sort(tr), np.setdiff1d(idx, tr)))
        _W.clear()
        _W.update(drg=drg, lc=lc, t=t, e=e, M=M, splits=splits)
        gcols = [gpos[g] for g in genes]
        seeds = rng.integers(0, 2**31 - 1, size=cfg.n_splits_bench + cfg.n_null_lasso)
        out = _r6_map([(gcols, int(seeds[k]), [k], True) for k in range(len(splits))], n_jobs)
        Cs_drg = [o[0][0] for o in out]
        risks = [o[1][0] for o in out]
        Xc = d[clin_cols].values.astype(float)
        G = km_surv(t, 1 - e)
        ca_rows = []
        for (tr, te), C_drg, r in zip(splits, Cs_drg, risks):
            row = dict(cancer=ca, n=len(d), events=int(e.sum()), clinical_covariates="+".join(clin_cols),
                       C_DRG=C_drg)
            try:
                m = PHReg(t[tr], Xc[tr], status=e[tr], ties="efron").fit(disp=False)
                rc = Xc @ m.params
                row["C_clinical"] = drg.cindex(t[te], e[te], rc[te])
                m2 = PHReg(t[tr], np.c_[Xc[tr], r[tr]], status=e[tr], ties="efron").fit(disp=False)
                rcd = np.c_[Xc, r] @ m2.params
                row["C_clinical_plus_DRG"] = drg.cindex(t[te], e[te], rcd[te])
            except Exception:  # noqa: BLE001
                row["C_clinical"] = row["C_clinical_plus_DRG"] = np.nan
            for yr in (1, 3, 5):
                row[f"AUC_IPCW_{yr}y"] = auc_ipcw(t[te], e[te], r[te], yr, G)
            ca_rows.append(row)
        LOG.info(f"  R6 [{k_ca}/{len(eligible)}] {ca} (n={len(d)}, events={int(e.sum())}): "
                 f"benchmark median C_DRG {np.median(Cs_drg):.3f} [{time.time() - t_ca:.0f}s]")
        ca_null, nulls = [], None
        if ca in robust and cfg.n_null_lasso > 0:
            sp = list(range(min(cfg.n_splits_null, len(splits))))
            obs = float(np.median([Cs_drg[k] for k in sp]))
            tasks = []
            for i in range(cfg.n_null_lasso):
                rg = [g for g in S.draw(genes) if g in gpos]
                tasks.append(([gpos[g] for g in rg], int(seeds[cfg.n_splits_bench + i]), sp, False))
            nulls, t_n = [], time.time()
            step = max(n_jobs * 4, 10)
            for j in range(0, len(tasks), step):
                nulls += [float(np.median(o[0])) for o in _r6_map(tasks[j:j + step], n_jobs)]
                done = len(nulls)
                eta = (time.time() - t_n) / done * (len(tasks) - done)
                LOG.info(f"    {ca} random-set null {done}/{len(tasks)}  (ETA {eta / 60:.1f} min)")
            nulls = np.array(nulls)
            ca_null.append(dict(cancer=ca, splits=len(sp), observed_median_test_C=obs,
                                null_median=float(np.median(nulls)), null_q95=float(np.quantile(nulls, .95)),
                                empirical_p=(1 + (nulls >= obs).sum()) / (1 + len(nulls)), n_null=len(nulls)))
            LOG.info(f"  {ca}: DRG median C {obs:.3f} vs random 24-gene sets {np.median(nulls):.3f} "
                     f"(P = {ca_null[-1]['empirical_p']:.3f})")
            np.save(T / f"R6_null_C_{ca}.npy", nulls)
        rows += ca_rows
        null_rows += ca_null
        if ck is not None:
            ck.write_bytes(pickle.dumps(dict(rows=ca_rows, null_rows=ca_null, nulls=nulls)))
        el = time.time() - t_start
        LOG.info(f"  R6 progress {k_ca}/{len(eligible)} cancers, elapsed {el / 60:.1f} min")
    _W.clear()
    res = pd.DataFrame(rows)
    res.to_csv(T / "R6a_lasso_benchmark_all_splits.csv.gz", index=False)
    summ = res.groupby("cancer").agg(n=("n", "first"), events=("events", "first"),
                                     clinical_covariates=("clinical_covariates", "first"),
                                     C_DRG=("C_DRG", "median"), C_clinical=("C_clinical", "median"),
                                     C_clinical_plus_DRG=("C_clinical_plus_DRG", "median"),
                                     AUC_IPCW_1y=("AUC_IPCW_1y", "median"), AUC_IPCW_3y=("AUC_IPCW_3y", "median"),
                                     AUC_IPCW_5y=("AUC_IPCW_5y", "median")).reset_index()
    dC = res.assign(dC=res.C_clinical_plus_DRG - res.C_clinical).groupby("cancer")["dC"]
    summ["deltaC_median"] = summ.cancer.map(dC.median())
    summ["deltaC_q025"] = summ.cancer.map(dC.quantile(.025))
    summ["deltaC_q975"] = summ.cancer.map(dC.quantile(.975))
    summ["pct_splits_DRG_adds"] = summ.cancer.map(dC.apply(lambda x: (x > 0).mean() * 100))
    summ["robust_v2"] = summ.cancer.isin(robust)
    summ.sort_values("C_DRG", ascending=False).to_csv(T / "R6b_lasso_benchmark_summary.csv", index=False)
    nr = pd.DataFrame(null_rows)
    nr.to_csv(T / "R6c_lasso_random_geneset_null.csv", index=False)
    return summ, nr


# ════════════════════════════════════════════════════════════════════════════
# figure + report
# ════════════════════════════════════════════════════════════════════════════
def make_figure(R1, R2, R6null, T, F):
    plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    nlas = len(R6null)
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6))
    ax = axes[0, 0]
    ax.hist(R1["null_pair"], bins=40, color="#BDBDBD")
    ax.axvline(R1["obs_pair"], color="#B2182B", lw=1.4)
    ax.set_title("A  Core vs cytoskeleton correlation", loc="left", fontsize=7.5)
    ax.set_xlabel("median within-tumour ρ (random matched pairs)")
    for k, (name, lab) in enumerate((("DRG_24", "B  DRG score (24 genes)"), ("cytoskeleton_14", "C  Cytoskeleton subset"))):
        ax = axes[0, 1 + k]
        st = R1["prog"][name]
        ax.hist(st["nsig"], bins=np.arange(-0.5, max(st["nsig"].max(), st["obs_nsig"]) + 1.5), color="#BDBDBD")
        ax.axvline(st["obs_nsig"], color="#B2182B", lw=1.4)
        ax.set_title(lab, loc="left", fontsize=7.5)
        ax.set_xlabel("cancers with FDR < 0.05 (OS)")
    ax = axes[1, 0]
    ax.hist(R2["null_med"], bins=40, color="#BDBDBD")
    ax.axvline(R2["obs_med"], color="#B2182B", lw=1.4)
    ax.set_title("D  SLC7A11 ~ oxidative PPP | purity, NRF2", loc="left", fontsize=7.5)
    ax.set_xlabel("median partial ρ (random matched genes)")
    for k in range(2):
        ax = axes[1, 1 + k]
        ax.axis("off")
    if nlas:
        ax = axes[1, 1]
        ax.axis("on")
        for i, r in R6null.reset_index(drop=True).iterrows():
            f = T / f"R6_null_C_{r.cancer}.npy"
            if f.exists():
                v = np.load(f)
                ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-.18, .18, len(v)), v, s=3,
                           color="#BDBDBD")
            ax.scatter([i], [r.observed_median_test_C], color="#B2182B", s=22, zorder=3)
        ax.set_xticks(range(nlas))
        ax.set_xticklabels(R6null.cancer)
        ax.axhline(.5, ls="--", color="grey", lw=.6)
        ax.set_ylabel("median test C-index")
        ax.set_title("E  DRG LASSO vs random 24-gene LASSO", loc="left", fontsize=7.5)
    fig.text(0.70, 0.25, "Red: observed DRG statistic\nGrey: expression-matched\nrandom gene sets",
             fontsize=6.5, va="center")
    fig.tight_layout()
    savefig(fig, F / "FigureS9_null_models")


def write_report(out, parts, t0):
    L = ["# DRG pan-cancer — v3 revision analyses\n",
         f"Generated {time.strftime('%Y-%m-%d %H:%M')} | runtime {(time.time() - t0) / 60:.1f} min\n"]
    for title, df in parts:
        if df is None or (hasattr(df, "empty") and df.empty):
            L.append(f"## {title}\n\n(no result)\n")
            continue
        L.append(f"## {title}\n")
        L.append(df.round(4).to_markdown(index=False) if hasattr(df, "to_markdown") else str(df))
        L.append("")
    (out / "REPORT_v3.md").write_text("\n".join(L), encoding="utf-8")


# ════════════════════════════════════════════════════════════════════════════
def main(argv=None):
    a = argparse.ArgumentParser()
    a.add_argument("--base", required=True)
    a.add_argument("--demo", action="store_true")
    a.add_argument("--n-null", type=int, default=1000)
    a.add_argument("--n-null-coherence", type=int, default=200)
    a.add_argument("--n-splits-bench", type=int, default=50)
    a.add_argument("--n-null-lasso", type=int, default=100)
    a.add_argument("--n-splits-null", type=int, default=10)
    a.add_argument("--seed", type=int, default=2026)
    a.add_argument("--n-jobs", type=int, default=0, help="parallel workers for R6 (0 = all CPUs)")
    a.add_argument("--skip", nargs="*", default=[], help="modules to skip, e.g. R6")
    cfg = a.parse_args(argv)
    t0 = time.time()
    base = Path(cfg.base)
    pre = "demo_" if cfg.demo else ""
    out = base / f"{pre}results_v3"
    stage = base / f"{pre}results_v3_running"
    if stage.exists():
        shutil.rmtree(stage)
    T, F = stage / "tables", stage / "figures"
    T.mkdir(parents=True)
    F.mkdir(parents=True)
    setup_logging(base / "logs" / f"{pre}v3_{time.strftime('%Y%m%d_%H%M%S')}.log")
    sys.path.insert(0, str(base / "code"))
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import drg_pancancer_pipeline as drg  # noqa: E402

    LOG.info("=" * 70)
    LOG.info("DRG v3 revision analyses")
    cache = base / f"{pre}data_cache" / "expression.parquet"
    expr = pd.read_parquet(cache) if cache.exists() else pd.read_pickle(cache.with_suffix(".pkl"))
    clin = pd.read_csv(base / f"{pre}results" / "tables" / "02_tumour_clinical_with_scores.csv", index_col=0)
    clin = clin[clin.index.isin(expr.columns)]
    expr = expr[clin.index]
    LOG.info(f"  expression {expr.shape}, tumours {len(clin)}, cancers {clin.cancer.nunique()}")
    core = [g for g in drg.DRG_CORE if g in expr.index]
    cyto = [g for g in drg.DRG_CYTOSKELETON if g in expr.index]
    genes = core + cyto
    LOG.info("  building rank engine (ssGSEA) ...")
    E = RankEngine(expr)
    # sanity check vs pipeline ssGSEA on 50 samples
    chk = drg.ssgsea(expr.iloc[:, :50], {"DRG": genes}, normalize=False).loc["DRG"]
    assert np.allclose(chk.values, E.score(genes).iloc[:50].values, atol=1e-6), "rank engine mismatch"
    S = MatchedSampler(expr, exclude=genes, seed=cfg.seed)
    ppp = drg.METABOLIC_SETS["Pentose_phosphate_NADPH"]
    S2 = MatchedSampler(expr, exclude=set(genes) | set(ppp) | set(OXPPP_STRICT) | set(NRF2_TARGETS), seed=cfg.seed + 1)
    oxphos_set = drg.METABOLIC_SETS["Oxidative_phosphorylation"]
    micro = {k: E.score([g for g in v if g not in genes]) for k, v in MICROENV.items()}
    micro = {k: v for k, v in micro.items() if v is not None}
    LOG.info(f"  microenvironment scores: {list(micro)}")
    pd.DataFrame(micro).to_csv(T / "R3_microenvironment_scores.csv.gz")
    oxphos_score = E.score([g for g in oxphos_set if g not in genes])

    parts, R1d, R2d, R6null = [], None, None, pd.DataFrame()
    steps = [
        ("R1", lambda: r1_null_models(E, S, expr, clin, core, cyto, cfg, T, F)),
        ("R2", lambda: r2_metabolic(E, S2, expr, clin, core, cyto, ppp, oxphos_set, micro, cfg, T)),
        ("R3", lambda: r3_microenv(E, clin, core, cyto, micro, oxphos_score, cfg, T) if micro else (None, None)),
        ("R4", lambda: r4_glioma(E, clin, core, cyto,
                                 load_glioma_molecular(base / f"{pre}data_raw", clin, drg, cfg), drg, cfg, T)),
        ("R5", lambda: r5_ph(E, clin, core, cyto, cfg, T)),
        ("R6", lambda: r6_lasso(expr, clin, genes, S, drg, base / f"{pre}results", cfg, T,
                                 ckdir=base / f"{pre}results_v3_checkpoints")),
    ]
    results = {}
    for name, fn in steps:
        if name in cfg.skip:
            continue
        try:
            results[name] = fn()
        except Exception as ex:  # noqa: BLE001
            LOG.exception(f"  module {name} failed: {ex}")
            results[name] = None
    if results.get("R1"):
        s, p, R1d = results["R1"]
        parts += [("R1a Structure null", s), ("R1b Prognostic null (OS)", p)]
    if results.get("R2"):
        summ, gnull, R2d = results["R2"]
        key = summ[summ.feature.isin(["SLC7A11", "LRPPRC", "NDUFS1", "core_subset", "cytoskeleton_subset", "DRG_score"])]
        parts += [("R2 Metabolic specificity (median partial rho)", key), ("R2c SLC7A11 gene-level null", gnull)]
    if results.get("R3") and results["R3"][0] is not None:
        r3, r3c = results["R3"]
        parts += [("R3a Subset Cox adjusted for purity + microenvironment (FDR<0.1)",
                   r3[r3.FDR < 0.1] if len(r3) else r3),
                  ("R3b Correlations adjusted for purity + microenvironment (median)",
                   r3c.drop(columns=["cancer"]).median(numeric_only=True).to_frame("median").reset_index()
                   if len(r3c) else r3c)]
    if results.get("R4"):
        parts += [("R4a Glioma molecular adjustment", results["R4"][0]), ("R4b IDH-stratified", results["R4"][1])]
    if results.get("R5"):
        r5, r5b = results["R5"]
        parts += [("R5a PH tests (q<0.1 shown)", r5[r5.PH_q < 0.1] if len(r5) else r5),
                  ("R5b Piecewise HRs", r5b)]
    if results.get("R6"):
        summ6, R6null = results["R6"]
        parts += [("R6b LASSO benchmark", summ6), ("R6c Random-gene-set LASSO null", R6null)]
    try:
        if R1d is not None and R2d is not None:
            make_figure(R1d, R2d, R6null, T, F)
    except Exception as ex:  # noqa: BLE001
        LOG.exception(f"  figure failed: {ex}")
    write_report(stage, parts, t0)
    (stage / "run_config_v3.json").write_text(json.dumps(vars(cfg), indent=2))
    if out.exists():
        shutil.rmtree(out)
    stage.rename(out)
    ck = base / f"{pre}results_v3_checkpoints"
    if results.get("R6") is not None and ck.exists():
        shutil.rmtree(ck)
    LOG.info(f"DONE in {(time.time() - t0) / 60:.1f} min → {out}")


if __name__ == "__main__":
    main()
