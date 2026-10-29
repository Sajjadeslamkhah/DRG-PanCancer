#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Disulfidptosis-related gene (DRG) signature — integrated pan-cancer analysis
of prognostic and metabolic effects (TCGA PanCanAtlas; UALCAN-style clinical
subgroup analyses reproduced from TCGA).

Single-run pipeline. All raw data, caches, tables, figures and logs live under
--base (on Colab: /content/drive/MyDrive/DRG_PanCancer).

Modules
  00  Download + cache TCGA PanCanAtlas data (UCSC Xena PanCanAtlas hub)
  01  DRG expression: tumour vs normal (per gene, per cancer)
  02  DRG signature score (ssGSEA, closed form; z-score mean as sensitivity)
  03  Univariate Cox per gene x cancer x endpoint (OS, DSS, DFI, PFI)
  04  Signature score survival: Cox per SD, KM (median / optimal cut),
      random-effects meta-analysis across cancers
  05  UALCAN-style clinical subgroups (stage, grade, gender, age, race)
  06  Metabolic programme analysis (ssGSEA of curated metabolic sets,
      DRGs excluded from sets to avoid circularity)
  07  Multivariable Cox (score + age + gender + stage; + tumour purity)
  08  Per-cancer LASSO-Cox DRG prognostic model, validated over repeated
      stratified 70/30 splits (nested CV; C-index, td-AUC, gene stability)
  09  Genomic alterations of DRGs (mutation, CNV; auto-detects thresholded
      GISTIC vs continuous log2 ratio input)
  10  DRG co-expression, summary report, auto-generated methods text
  11  Sub-scores (core screen vs actin-cytoskeleton arm): tumour/normal,
      survival and metabolic coupling
  v2: tumour-purity (ABSOLUTE) partial correlations for all metabolic analyses

Expression: EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2 (unit log2(norm+1)).
Survival:   TCGA-CDR (Liu et al., Cell 2018) — OS, DSS, DFI, PFI.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import platform
import re
import shutil
import sys
import time
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

warnings.filterwarnings("ignore")
LOG = logging.getLogger("drg")

# ════════════════════════════════════════════════════════════════════════════
# 1. GENE SETS
# ════════════════════════════════════════════════════════════════════════════
# Liu X. et al., Nat Cell Biol 2023;25:404-414 (disulfidptosis discovery):
# CRISPR-screen genes + disulfide-stressed actin cytoskeleton proteins.
DRG_CORE = ["SLC7A11", "SLC3A2", "RPN1", "NCKAP1", "NUBPL",
            "NDUFA11", "LRPPRC", "OXSM", "NDUFS1", "GYS1"]
DRG_CYTOSKELETON = ["FLNA", "FLNB", "MYH9", "TLN1", "ACTB", "MYL6", "MYH10",
                    "CAPZB", "DSTN", "IQGAP1", "ACTN4", "PDLIM1", "CD2AP", "INF2"]
DRG_RAC_WRC = ["RAC1", "CYFIP1", "BRK1", "WASF2", "ABI2"]

SIGNATURES = {
    "core": DRG_CORE,
    "standard": DRG_CORE + DRG_CYTOSKELETON,            # 24-gene set (most used)
    "extended": DRG_CORE + DRG_CYTOSKELETON + DRG_RAC_WRC,
}
GENE_CATEGORY = {**{g: "Core (CRISPR screen)" for g in DRG_CORE},
                 **{g: "Actin cytoskeleton" for g in DRG_CYTOSKELETON},
                 **{g: "Rac1-WRC axis" for g in DRG_RAC_WRC}}

# Curated metabolic programmes (HGNC + legacy symbols present in the 2016
# PanCanAtlas annotation, e.g. ATP5A1). Missing genes are dropped silently.
METABOLIC_SETS = {
    "Glycolysis": ["HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA",
                   "ALDOC", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2",
                   "PKM", "PKM2", "LDHA", "PFKFB3", "PFKFB4", "SLC16A3"],
    "Glucose_uptake": ["SLC2A1", "SLC2A3", "SLC2A4", "SLC2A14", "SLC5A1", "TXNIP", "ARRDC4"],
    "Pentose_phosphate_NADPH": ["G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT",
                                "TALDO1", "H6PD", "ME1", "IDH1", "NNT",
                                "MTHFD1", "MTHFD2", "ALDH1L2"],
    "TCA_cycle": ["CS", "ACO2", "IDH3A", "IDH3B", "IDH3G", "OGDH", "DLST",
                  "SUCLG1", "SUCLA2", "SDHA", "SDHB", "FH", "MDH2", "PDHA1",
                  "PDHB", "DLAT", "DLD"],
    "Oxidative_phosphorylation": [
        "NDUFA1", "NDUFA2", "NDUFA4", "NDUFA9", "NDUFB3", "NDUFB8", "NDUFS2",
        "NDUFS3", "NDUFV1", "NDUFV2", "SDHC", "SDHD", "UQCRC1", "UQCRC2",
        "UQCRFS1", "CYC1", "COX4I1", "COX5A", "COX5B", "COX6C", "COX7B",
        "ATP5A1", "ATP5B", "ATP5C1", "ATP5O", "ATP5G1",
        "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5PO", "ATP5MC1"],
    "Glutathione_metabolism": ["GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2",
                               "GPX4", "GSTP1", "GSTM1", "GGT1", "CHAC1", "IDH2"],
    "Cysteine_sulfur_metabolism": ["CBS", "CTH", "CDO1", "CSAD", "MPST", "GOT1",
                                   "SLC1A4", "SLC7A9", "SLC7A10", "CARS1", "CARS"],
    "Thioredoxin_system": ["TXN", "TXNRD1", "TXNRD2", "TXN2", "PRDX1", "PRDX2",
                           "PRDX3", "PRDX4", "SRXN1"],
    "NRF2_antioxidant_response": ["NFE2L2", "NQO1", "HMOX1", "AKR1C1", "AKR1C2",
                                  "AKR1C3", "AKR1B10", "CYP4F11", "ABCC2",
                                  "ABCC3", "OSGIN1"],
    "Glutamine_metabolism": ["SLC1A5", "SLC38A1", "SLC38A2", "GLS", "GLS2",
                             "GLUD1", "GOT2", "GPT2", "ASNS", "GLUL"],
    "Serine_one_carbon": ["PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2",
                          "MTHFD1L", "MTHFR", "MAT2A"],
    "Fatty_acid_synthesis": ["ACLY", "ACACA", "FASN", "SCD", "ELOVL6",
                             "ACSL1", "ACSS2", "SREBF1", "INSIG1"],
    "Fatty_acid_oxidation": ["CPT1A", "CPT1B", "CPT2", "ACADM", "ACADL",
                             "ACADVL", "HADHA", "HADHB", "ECHS1", "ACAA2", "PPARA"],
    "Glycogen_metabolism": ["GYS2", "GBE1", "PYGL", "PYGM", "PYGB", "UGP2",
                            "GYG1", "AGL", "PGM1", "PPP1R3C"],
    "Hypoxia_core": ["HIF1A", "EPAS1", "CA9", "VEGFA", "ADM", "BNIP3", "PDK1",
                     "NDRG1", "ANKRD37", "P4HA1", "EGLN3", "ANGPTL4"],
    "mTORC1_core": ["MTOR", "RPTOR", "RPS6KB1", "EIF4EBP1", "SLC7A5", "DDIT4",
                    "SESN2", "ATF4"],
}

# TCGA project names (phenotype fallback) -> abbreviations
DISEASE2ABBR = {
    "adrenocortical cancer": "ACC", "bladder urothelial carcinoma": "BLCA",
    "breast invasive carcinoma": "BRCA", "cervical & endocervical cancer": "CESC",
    "cholangiocarcinoma": "CHOL", "colon adenocarcinoma": "COAD",
    "diffuse large B-cell lymphoma": "DLBC", "esophageal carcinoma": "ESCA",
    "glioblastoma multiforme": "GBM", "head & neck squamous cell carcinoma": "HNSC",
    "kidney chromophobe": "KICH", "kidney clear cell carcinoma": "KIRC",
    "kidney papillary cell carcinoma": "KIRP", "acute myeloid leukemia": "LAML",
    "brain lower grade glioma": "LGG", "liver hepatocellular carcinoma": "LIHC",
    "lung adenocarcinoma": "LUAD", "lung squamous cell carcinoma": "LUSC",
    "mesothelioma": "MESO", "ovarian serous cystadenocarcinoma": "OV",
    "pancreatic adenocarcinoma": "PAAD", "pheochromocytoma & paraganglioma": "PCPG",
    "prostate adenocarcinoma": "PRAD", "rectum adenocarcinoma": "READ",
    "sarcoma": "SARC", "skin cutaneous melanoma": "SKCM",
    "stomach adenocarcinoma": "STAD", "testicular germ cell tumor": "TGCT",
    "thyroid carcinoma": "THCA", "thymoma": "THYM",
    "uterine corpus endometrioid carcinoma": "UCEC",
    "uterine carcinosarcoma": "UCS", "uveal melanoma": "UVM",
}

# ════════════════════════════════════════════════════════════════════════════
# 2. DATA SOURCES (UCSC Xena — TCGA PanCanAtlas hub)
# ════════════════════════════════════════════════════════════════════════════
HUB = "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
DATA_FILES = {
    "expression": ("EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
                   "EBpp_AdjustPANCAN_RNASeqV2_geneExp.xena.gz", True),
    "survival": ("Survival_SupplementalTable_S1_20171025_xena_sp",
                 "TCGA_CDR_survival_xena.tsv", True),
    "phenotype": ("TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
                  "TCGA_phenotype_dense.tsv.gz", False),
    "mutation": ("mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz",
                 "mc3_nonsilent_gene.xena.gz", False),
    "cnv": ("broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz",
            "GISTIC2_thresholded_gene.xena.gz", False),
    # ABSOLUTE purity/ploidy (TCGA_mastercalls.abs_tables_JSedit.fixed.txt, GDC PanCanAtlas)
    "purity": ("https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
               "TCGA_ABSOLUTE_purity.txt", False),
}

# Continuous gene-level log2(CN/2) thresholds (used when CNV input is not thresholded)
CNV_GAIN, CNV_LOSS, CNV_AMP, CNV_DEL = 0.3, -0.3, 1.0, -1.0

ENDPOINTS = ["OS", "DSS", "DFI", "PFI"]


# ════════════════════════════════════════════════════════════════════════════
# 3. UTILITIES
# ════════════════════════════════════════════════════════════════════════════
def setup_logging(log_file: Path):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S")
    LOG.setLevel(logging.INFO)
    LOG.handlers.clear()
    for h in (logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, mode="a")):
        h.setFormatter(fmt)
        LOG.addHandler(h)


def download(url: str, dest: Path, retries: int = 3):
    if dest.exists() and dest.stat().st_size > 0:
        LOG.info(f"  cached: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            LOG.info(f"  downloading {dest.name} (attempt {attempt}) ...")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                total = int(r.headers.get("Content-Length", 0))
                done, last = 0, 0
                while True:
                    buf = r.read(1 << 22)
                    if not buf:
                        break
                    f.write(buf)
                    done += len(buf)
                    if total and done - last > total / 10:
                        LOG.info(f"    {done / 1e6:.0f}/{total / 1e6:.0f} MB")
                        last = done
            tmp.rename(dest)
            return dest
        except Exception as e:  # noqa: BLE001
            LOG.warning(f"  download failed: {e}")
            time.sleep(5 * attempt)
    if tmp.exists():
        tmp.unlink()
    raise RuntimeError(f"Could not download {url}")


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def stars(q):
    if pd.isna(q):
        return ""
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def savefig(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    return (x - np.nanmean(x)) / sd if sd > 0 else x * 0


# ════════════════════════════════════════════════════════════════════════════
# 4. SURVIVAL ENGINE (statsmodels + vectorised LASSO-Cox)
# ════════════════════════════════════════════════════════════════════════════
def cox_fit(time_, event, X: pd.DataFrame):
    """Cox PH (Efron ties). Returns DataFrame: HR, CI, p per covariate."""
    m = PHReg(np.asarray(time_, float), np.asarray(X, float),
              status=np.asarray(event, int), ties="efron").fit(disp=False)
    ci = np.exp(m.conf_int())
    return pd.DataFrame({"coef": m.params, "se": m.bse, "HR": np.exp(m.params),
                         "HR_low": ci[:, 0], "HR_high": ci[:, 1],
                         "p": m.pvalues}, index=X.columns)


def logrank(time_, event, group):
    try:
        _, p = survdiff(np.asarray(time_, float), np.asarray(event, int),
                        np.asarray(group))
        return float(p)
    except Exception:  # noqa: BLE001
        return np.nan


def cindex(time_, event, risk):
    """Harrell's C (higher risk -> shorter time)."""
    t, e, r = (np.asarray(a, float) for a in (time_, event, risk))
    ev = e == 1
    T1 = t[ev][:, None]
    R1 = r[ev][:, None]
    comparable = t[None, :] > T1
    conc = (R1 > r[None, :]) & comparable
    ties = (R1 == r[None, :]) & comparable
    n = comparable.sum()
    return float((conc.sum() + 0.5 * ties.sum()) / n) if n else np.nan


def td_auc(time_, event, risk, t0):
    """Cumulative/dynamic AUC at t0 (cases: event<=t0; controls: time>t0)."""
    t, e, r = (np.asarray(a, float) for a in (time_, event, risk))
    cases = (t <= t0) & (e == 1)
    ctrls = t > t0
    if cases.sum() < 3 or ctrls.sum() < 3:
        return np.nan
    rc, rn = r[cases][:, None], r[ctrls][None, :]
    return float(((rc > rn).sum() + 0.5 * (rc == rn).sum()) / (rc.size * rn.size))


class LassoCox:
    """L1-penalised Cox (Breslow) via FISTA; objective -loglik/n + lam*|b|_1."""

    def __init__(self, max_iter=3000, tol=1e-7):
        self.max_iter, self.tol = max_iter, tol

    @staticmethod
    def _prep(t, e, X):
        o = np.argsort(t, kind="mergesort")
        t, e, X = t[o], e[o], X[o]
        _, first = np.unique(t, return_index=True)
        grp = np.searchsorted(t, t, side="left")  # first index of tie group
        last = np.searchsorted(t, t, side="right") - 1
        return t, e, X, grp, last

    def _grad_loss(self, b, e, X, grp, last):
        n = len(e)
        eta = X @ b
        eta -= eta.max()
        w = np.exp(eta)
        S = np.cumsum(w[::-1])[::-1][grp]           # risk set sums
        loss = -(e * (eta - np.log(S))).sum() / n
        cum = np.cumsum(e / S)[last]
        g_eta = -(e - w * cum) / n
        return loss, X.T @ g_eta

    def path(self, t, e, X, lambdas):
        t, e, X, grp, last = self._prep(np.asarray(t, float), np.asarray(e, float),
                                        np.asarray(X, float))
        n, p = X.shape
        L = max(np.linalg.eigvalsh(X.T @ X / n).max(), 1e-8)
        step = 1.0 / L
        b = np.zeros(p)
        coefs = []
        for lam in lambdas:
            y, bk, tk = b.copy(), b.copy(), 1.0
            for _ in range(self.max_iter):
                _, g = self._grad_loss(y, e, X, grp, last)
                z = y - step * g
                bn = np.sign(z) * np.maximum(np.abs(z) - step * lam, 0)
                tn = (1 + np.sqrt(1 + 4 * tk * tk)) / 2
                y = bn + (tk - 1) / tn * (bn - bk)
                if np.max(np.abs(bn - bk)) < self.tol:
                    bk = bn
                    break
                bk, tk = bn, tn
            b = bk
            coefs.append(b.copy())
        return np.array(coefs)

    def lambda_max(self, t, e, X):
        t, e, X, grp, last = self._prep(np.asarray(t, float), np.asarray(e, float),
                                        np.asarray(X, float))
        _, g = self._grad_loss(np.zeros(X.shape[1]), e, X, grp, last)
        return float(np.abs(g).max())


def km_plot(ax, time_, event, group, labels, colors, title, p=None):
    for g, lab, col in zip(sorted(np.unique(group)), labels, colors):
        m = group == g
        if m.sum() < 2:
            continue
        sf = SurvfuncRight(np.asarray(time_)[m], np.asarray(event)[m])
        x = np.r_[0, sf.surv_times]
        y = np.r_[1, sf.surv_prob]
        ax.step(x, y, where="post", color=col, lw=1.4, label=f"{lab} (n={m.sum()})")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Years")
    ax.set_ylabel("Survival probability")
    ax.set_title(title + (f"\nlog-rank p={p:.2g}" if p is not None and not np.isnan(p) else ""),
                 fontsize=8)
    ax.legend(fontsize=6, frameon=False, loc="lower left")
    ax.tick_params(labelsize=7)


# ════════════════════════════════════════════════════════════════════════════
# 5. ssGSEA (closed form — identical to GSVA::ssgsea, tau = 0.25)
# ════════════════════════════════════════════════════════════════════════════
def ssgsea(expr: pd.DataFrame, gene_sets: dict, alpha=0.25, min_size=3,
           chunk=600, normalize=True):
    """
    ES_j = sum_i [P_hit(i) - P_miss(i)] over the ranked list; summing the
    running sums gives the exact closed form
      hit  = sum_{g in S} r_g^(alpha+1) / sum_{g in S} r_g^alpha
      miss = (N(N+1)/2 - sum_{g in S} r_g) / (N - |S|)
    with r = ascending expression rank (highest expressed gene = N).
    Normalisation follows GSVA (divide by range of the whole score matrix).
    """
    idx = {g: i for i, g in enumerate(expr.index)}
    sets = {k: np.array([idx[g] for g in dict.fromkeys(v) if g in idx])
            for k, v in gene_sets.items()}
    sets = {k: v for k, v in sets.items() if len(v) >= min_size}
    N = expr.shape[0]
    total = N * (N + 1) / 2
    X = expr.values
    out = np.zeros((len(sets), X.shape[1]))
    for s0 in range(0, X.shape[1], chunk):
        blk = X[:, s0:s0 + chunk].astype(np.float64)
        blk = np.where(np.isnan(blk), -np.inf, blk)
        R = rankdata(blk, axis=0)
        for k, ii in enumerate(sets.values()):
            r = R[ii, :]
            w = r ** alpha
            out[k, s0:s0 + chunk] = (w * r).sum(0) / w.sum(0) - \
                (total - r.sum(0)) / (N - len(ii))
    df = pd.DataFrame(out, index=list(sets), columns=expr.columns)
    if normalize:
        rng = df.values.max() - df.values.min()
        if rng > 0:
            df = df / rng
    return df


# ════════════════════════════════════════════════════════════════════════════
# 6. DATA LOADING
# ════════════════════════════════════════════════════════════════════════════
def load_expression(path: Path, cache: Path) -> pd.DataFrame:
    if cache.exists():
        LOG.info(f"  loading cached expression matrix: {cache.name}")
        return pd.read_parquet(cache) if cache.suffix == ".parquet" else pd.read_pickle(cache)
    LOG.info("  parsing expression matrix (one-time, ~5-10 min) ...")
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    dtypes = {c: np.float32 for c in header[1:]}
    parts = [ch for ch in pd.read_csv(path, sep="\t", index_col=0, dtype=dtypes,
                                      chunksize=2000, na_values=["NA", "NaN", ""])]
    df = pd.concat(parts)
    genes = df.index.astype(str).str.split("|").str[0]
    df.index = genes
    df = df[(genes != "?") & (genes != "") & (genes != "nan")]
    if df.index.duplicated().any():
        df = df.assign(_m=df.mean(axis=1)).sort_values("_m", ascending=False)
        df = df[~df.index.duplicated(keep="first")].drop(columns="_m")
    df.columns = df.columns.astype(str).str[:15]
    df = df.loc[:, ~df.columns.duplicated(keep="first")]
    df = df.loc[df.isna().mean(axis=1) < 0.5]
    try:
        df.to_parquet(cache)
    except Exception:  # noqa: BLE001
        cache = cache.with_suffix(".pkl")
        df.to_pickle(cache)
    LOG.info(f"  expression: {df.shape[0]} genes x {df.shape[1]} samples (cached)")
    return df


def norm_stage(s):
    if pd.isna(s):
        return np.nan
    m = re.search(r"STAGE\s*(IV|III|II|I)(?![V])", str(s).upper())
    return {"I": "Stage I", "II": "Stage II", "III": "Stage III",
            "IV": "Stage IV"}[m.group(1)] if m else np.nan


def norm_grade(s):
    if pd.isna(s):
        return np.nan
    s = str(s).upper()
    m = re.search(r"G([1-4])", s)
    if m:
        return f"G{m.group(1)}"
    if "HIGH" in s:
        return "High grade"
    if "LOW" in s:
        return "Low grade"
    return np.nan


def build_annotation(expr_cols, surv: pd.DataFrame, pheno: pd.DataFrame | None):
    ann = pd.DataFrame(index=pd.Index(expr_cols, name="sample"))
    ann["patient"] = ann.index.str[:12]
    ann["code"] = ann.index.str[13:15]
    ctcol = next((c for c in surv.columns if "cancer type" in c.lower()), None)
    pat_ct = surv.dropna(subset=[ctcol]).groupby("_PATIENT")[ctcol].first() if ctcol else pd.Series(dtype=str)
    ann["cancer"] = ann["patient"].map(pat_ct)
    if pheno is not None and ann["cancer"].isna().any():
        dcol = next((c for c in pheno.columns if "primary_disease" in c), None)
        if dcol:
            ph = pheno.set_index(pheno.columns[0])[dcol].str.lower().map(DISEASE2ABBR)
            ph.index = ph.index.astype(str).str[:15]
            ph = ph[~ph.index.duplicated()]
            fill = ann.index.to_series().map(ph)
            ann["cancer"] = ann["cancer"].fillna(fill)
    tumour = (ann["code"] == "01") | ((ann["cancer"] == "LAML") & (ann["code"] == "03")) | \
             ((ann["cancer"] == "SKCM") & (ann["code"] == "06"))
    ann["group"] = np.where(tumour, "Tumor", np.where(ann["code"] == "11", "Normal", "Other"))
    ann = ann.dropna(subset=["cancer"])
    return ann


def build_clinical(ann: pd.DataFrame, surv: pd.DataFrame):
    """One primary tumour per patient + clinical/survival variables."""
    tum = ann[ann["group"] == "Tumor"].copy()
    tum["pref"] = (tum["code"] != "01").astype(int)
    tum = tum.sort_values(["patient", "pref"]).drop_duplicates("patient")
    s = surv.drop_duplicates("_PATIENT").set_index("_PATIENT")
    clin = tum[["patient", "cancer"]].join(s, on="patient", rsuffix="_s")
    out = pd.DataFrame(index=clin.index)
    out["patient"] = clin["patient"]
    out["cancer"] = clin["cancer"]
    out["age"] = pd.to_numeric(clin.get("age_at_initial_pathologic_diagnosis"), errors="coerce")
    out["gender"] = clin.get("gender", pd.Series(index=clin.index, dtype=object)).str.upper()
    out["race"] = clin.get("race", pd.Series(index=clin.index, dtype=object))
    st = clin.get("ajcc_pathologic_tumor_stage", pd.Series(index=clin.index, dtype=object)).map(norm_stage)
    st2 = clin.get("clinical_stage", pd.Series(index=clin.index, dtype=object)).map(norm_stage)
    out["stage"] = st.fillna(st2)
    out["grade"] = clin.get("histological_grade", pd.Series(index=clin.index, dtype=object)).map(norm_grade)
    for ep in ENDPOINTS:
        out[ep] = pd.to_numeric(clin.get(ep), errors="coerce")
        out[f"{ep}.time"] = pd.to_numeric(clin.get(f"{ep}.time"), errors="coerce") / 365.25
    out["age_group"] = pd.Series(np.where(out["age"] > 60, ">60", "<=60"),
                                 index=out.index, dtype=object).where(out["age"].notna())
    out["race"] = out["race"].where(~out["race"].astype(str).str.contains(
        r"\[|NOT|UNKNOWN|nan", case=False, regex=True), np.nan)
    return out


def load_purity(path: Path) -> pd.Series:
    """ABSOLUTE purity per 15-char sample barcode (robust to column naming / gzip)."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    df = pd.read_csv(path, sep="\t", compression="gzip" if gz else None, low_memory=False)
    pcol = next(c for c in df.columns if c.strip().lower() == "purity")
    idcol = next(c for c in df.columns
                 if df[c].astype(str).str.match(r"^TCGA-\w{2}-\w{4}-\d{2}").mean() > 0.9)
    s = pd.to_numeric(df[pcol], errors="coerce")
    s.index = df[idcol].astype(str).str[:15]
    return s[~s.index.duplicated()].dropna()


def partial_spearman(x, y, z):
    """Spearman correlation of x and y controlling for z (rank residual method)."""
    x, y, z = (rankdata(np.asarray(a, float)) for a in (x, y, z))
    Z = np.c_[np.ones_like(z), z]
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    r = np.corrcoef(rx, ry)[0, 1]
    n = len(x)
    t = r * np.sqrt((n - 3) / max(1e-12, 1 - r * r))
    return r, 2 * stats.t.sf(abs(t), n - 3)


SCORE_LABELS = {"score": "DRG_score", "score_core": "DRG_core", "score_cyto": "DRG_cytoskeleton"}


# ════════════════════════════════════════════════════════════════════════════
# 7. ANALYSIS MODULES
# ════════════════════════════════════════════════════════════════════════════
def mod01_tumor_normal(expr, ann, genes, score_all, P, cfg):
    LOG.info("[01] Tumour vs normal expression")
    rows = []
    for ca, sub in ann[ann["group"].isin(["Tumor", "Normal"])].groupby("cancer"):
        t = sub.index[sub["group"] == "Tumor"]
        n = sub.index[sub["group"] == "Normal"]
        if len(n) < cfg.min_normals or len(t) < cfg.min_normals:
            continue
        feats = {g: expr.loc[g] for g in genes}
        feats["DRG_score"] = score_all
        for g, v in feats.items():
            a, b = v[t].dropna(), v[n].dropna()
            if len(a) < 3 or len(b) < 3:
                continue
            p = stats.mannwhitneyu(a, b).pvalue
            rows.append(dict(cancer=ca, gene=g, n_tumor=len(a), n_normal=len(b),
                             mean_tumor=a.mean(), mean_normal=b.mean(),
                             log2FC=a.mean() - b.mean(), p=p))
    res = pd.DataFrame(rows)
    if res.empty:
        LOG.warning("  no cancer with enough normals")
        return res
    res["FDR"] = np.nan
    gm = res["gene"] != "DRG_score"
    res.loc[gm, "FDR"] = bh(res.loc[gm, "p"])
    res.loc[~gm, "FDR"] = bh(res.loc[~gm, "p"])
    res.to_csv(P["tab"] / "01_tumor_vs_normal.csv", index=False)

    g = res[gm]
    fc = g.pivot(index="gene", columns="cancer", values="log2FC").reindex(genes)
    q = g.pivot(index="gene", columns="cancer", values="FDR").reindex(genes)
    fig, ax = plt.subplots(figsize=(0.45 * fc.shape[1] + 3, 0.32 * fc.shape[0] + 1.5))
    lim = np.nanpercentile(np.abs(fc.values), 98) if np.isfinite(fc.values).any() else 1
    sns.heatmap(fc, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, ax=ax,
                annot=q.map(stars) if hasattr(q, "map") else q.applymap(stars),
                fmt="", annot_kws={"size": 6}, linewidths=.3, linecolor="white",
                cbar_kws={"label": "log2 fold change (tumour − normal)"})
    ax.set_title("DRG expression: tumour vs normal (BH-FDR: * <.05, ** <.01, *** <.001)", fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")
    savefig(fig, P["fig"] / "01_tumor_normal_log2FC_heatmap")

    # Signature score boxplot
    sub = ann[ann["group"].isin(["Tumor", "Normal"])]
    keep = res.loc[~gm, "cancer"].unique()
    df = pd.DataFrame({"cancer": sub["cancer"], "group": sub["group"],
                       "score": score_all.reindex(sub.index)})
    df = df[df["cancer"].isin(keep)]
    if len(df):
        fig, ax = plt.subplots(figsize=(0.55 * len(keep) + 2, 3.6))
        sns.boxplot(data=df, x="cancer", y="score", hue="group", order=sorted(keep),
                    hue_order=["Normal", "Tumor"], palette={"Normal": "#4C78A8", "Tumor": "#E45756"},
                    fliersize=0, linewidth=.6, ax=ax)
        sc = res[~gm].set_index("cancer")
        ymax = df["score"].max()
        for i, ca in enumerate(sorted(keep)):
            ax.text(i, ymax, stars(sc.loc[ca, "FDR"]), ha="center", fontsize=8)
        ax.set_ylabel("DRG ssGSEA score")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=90, labelsize=7)
        ax.legend(frameon=False, fontsize=7)
        savefig(fig, P["fig"] / "01_DRG_score_tumor_vs_normal_boxplot")
    return res


def mod03_gene_cox(expr, clin, genes, P, cfg):
    LOG.info("[03] Univariate Cox: gene x cancer x endpoint")
    rows = []
    for ca, sub in clin.groupby("cancer"):
        for ep in cfg.endpoints:
            d = sub[[ep, f"{ep}.time"]].dropna()
            d = d[d[f"{ep}.time"] > 0]
            if len(d) < cfg.min_samples or d[ep].sum() < cfg.min_events:
                continue
            for g in genes:
                x = zscore(expr.loc[g, d.index].values)
                if np.isnan(x).any() or np.nanstd(x) == 0:
                    continue
                try:
                    r = cox_fit(d[f"{ep}.time"], d[ep], pd.DataFrame({g: x}, index=d.index)).iloc[0]
                    rows.append(dict(cancer=ca, endpoint=ep, gene=g, n=len(d),
                                     events=int(d[ep].sum()), HR=r.HR, HR_low=r.HR_low,
                                     HR_high=r.HR_high, p=r.p))
                except Exception:  # noqa: BLE001
                    pass
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["FDR"] = res.groupby("endpoint")["p"].transform(lambda s: bh(s.values))
    res["effect"] = np.where(res["p"] >= 0.05, "NS", np.where(res["HR"] > 1, "Risk", "Protective"))
    res.to_csv(P["tab"] / "03_gene_univariate_cox.csv", index=False)
    for ep in res["endpoint"].unique():
        r = res[res["endpoint"] == ep]
        lhr = np.log2(r.pivot(index="gene", columns="cancer", values="HR")).reindex(genes)
        pv = r.pivot(index="gene", columns="cancer", values="p").reindex(genes)
        fig, ax = plt.subplots(figsize=(0.42 * lhr.shape[1] + 3, 0.32 * lhr.shape[0] + 1.5))
        lim = np.nanpercentile(np.abs(lhr.values), 98)
        ann_ = pv.apply(lambda c: c.map(lambda v: "" if pd.isna(v) else
                                        ("**" if v < 0.01 else "*" if v < 0.05 else "")))
        sns.heatmap(lhr, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=ann_, fmt="",
                    annot_kws={"size": 6}, linewidths=.3, linecolor="white", ax=ax,
                    cbar_kws={"label": "log2 HR per SD (red = risk)"})
        ax.set_title(f"Univariate Cox — {ep} (* p<.05, ** p<.01)", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("")
        savefig(fig, P["fig"] / f"03_gene_cox_heatmap_{ep}")
    # summary count plot
    cnt = res[res["effect"] != "NS"].groupby(["gene", "effect"]).size().unstack(fill_value=0)
    if len(cnt):
        cnt = cnt.reindex(genes).fillna(0)
        fig, ax = plt.subplots(figsize=(7, 0.28 * len(genes) + 1))
        cnt.plot.barh(stacked=True, ax=ax, color={"Risk": "#E45756", "Protective": "#4C78A8"},
                      width=.8)
        ax.invert_yaxis()
        ax.set_xlabel("No. of cancer × endpoint combinations (p<0.05)")
        ax.legend(frameon=False, fontsize=7)
        savefig(fig, P["fig"] / "03_gene_cox_risk_protective_counts")
    return res


def optimal_cut(t, e, x, lo=0.2, hi=0.8):
    best, bp = np.median(x), 1.0
    for c in np.unique(np.quantile(x, np.linspace(lo, hi, 25))):
        g = (x > c).astype(int)
        p = logrank(t, e, g)
        if not np.isnan(p) and p < bp:
            best, bp = c, p
    return best


def random_effects_meta(logHR, se):
    logHR, se = np.asarray(logHR), np.asarray(se)
    w = 1 / se ** 2
    fe = (w * logHR).sum() / w.sum()
    Q = (w * (logHR - fe) ** 2).sum()
    df = len(logHR) - 1
    tau2 = max(0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())) if df > 0 else 0
    wr = 1 / (se ** 2 + tau2)
    re = (wr * logHR).sum() / wr.sum()
    se_re = np.sqrt(1 / wr.sum())
    I2 = max(0, (Q - df) / Q) * 100 if Q > 0 else 0
    p = 2 * stats.norm.sf(abs(re / se_re))
    return dict(HR=np.exp(re), HR_low=np.exp(re - 1.96 * se_re),
                HR_high=np.exp(re + 1.96 * se_re), p=p, I2=I2, tau2=tau2, Q_p=stats.chi2.sf(Q, df) if df > 0 else np.nan)


def mod04_score_survival(clin, P, cfg):
    LOG.info("[04] Signature score survival (Cox, KM, meta-analysis)")
    rows, km_data = [], {}
    for ca, sub in clin.groupby("cancer"):
        for ep in cfg.endpoints:
            d = sub[["score", ep, f"{ep}.time"]].dropna()
            d = d[d[f"{ep}.time"] > 0]
            if len(d) < cfg.min_samples or d[ep].sum() < cfg.min_events:
                continue
            x = zscore(d["score"].values)
            try:
                r = cox_fit(d[f"{ep}.time"], d[ep], pd.DataFrame({"score": x}, index=d.index)).iloc[0]
            except Exception:  # noqa: BLE001
                continue
            cut = np.median(d["score"]) if cfg.cutoff == "median" else \
                optimal_cut(d[f"{ep}.time"].values, d[ep].values, d["score"].values)
            grp = (d["score"] > cut).astype(int).values
            p_lr = logrank(d[f"{ep}.time"], d[ep], grp)
            hl = cox_fit(d[f"{ep}.time"], d[ep], pd.DataFrame({"high": grp}, index=d.index)).iloc[0]
            rows.append(dict(cancer=ca, endpoint=ep, n=len(d), events=int(d[ep].sum()),
                             HR_perSD=r.HR, HR_low=r.HR_low, HR_high=r.HR_high, coef=r.coef,
                             se=r.se, p_cox=r.p, cutoff=cut, HR_high_vs_low=hl.HR,
                             p_logrank=p_lr, C_index=cindex(d[f"{ep}.time"], d[ep], x)))
            km_data[(ca, ep)] = (d[f"{ep}.time"].values, d[ep].values, grp, p_lr)
    res = pd.DataFrame(rows)
    if res.empty:
        return res, pd.DataFrame()
    res["FDR_cox"] = res.groupby("endpoint")["p_cox"].transform(lambda s: bh(s.values))
    res.to_csv(P["tab"] / "04_score_cox_km.csv", index=False)

    metas = []
    for ep in res["endpoint"].unique():
        r = res[res["endpoint"] == ep].sort_values("HR_perSD")
        m = random_effects_meta(r["coef"], r["se"])
        metas.append(dict(endpoint=ep, k=len(r), **m))
        # forest
        fig, ax = plt.subplots(figsize=(5.2, 0.24 * len(r) + 1.6))
        y = np.arange(len(r))
        col = np.where(r["p_cox"] < 0.05, np.where(r["HR_perSD"] > 1, "#E45756", "#4C78A8"), "grey")
        ax.errorbar(r["HR_perSD"], y, xerr=[r["HR_perSD"] - r["HR_low"], r["HR_high"] - r["HR_perSD"]],
                    fmt="none", ecolor=col, elinewidth=1)
        ax.scatter(r["HR_perSD"], y, c=col, s=14, zorder=3)
        ax.scatter([m["HR"]], [-1.5], marker="D", c="black", s=30)
        ax.hlines(-1.5, m["HR_low"], m["HR_high"], color="black")
        ax.set_yticks(list(y) + [-1.5])
        ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(r["cancer"], r["n"])] +
                           [f"RE meta (I²={m['I2']:.0f}%)"], fontsize=6.5)
        ax.axvline(1, ls="--", c="k", lw=.6)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("HR per SD of DRG score (95% CI)")
        ax.set_title(f"DRG signature score — {ep}", fontsize=9)
        savefig(fig, P["fig"] / f"04_score_forest_{ep}")
        # KM grid
        keys = [k for k in km_data if k[1] == ep]
        nc = 6
        nr = int(np.ceil(len(keys) / nc))
        fig, axes = plt.subplots(nr, nc, figsize=(nc * 2.6, nr * 2.4), squeeze=False)
        for a, k in zip(axes.flat, sorted(keys)):
            t, e, g, p = km_data[k]
            km_plot(a, t, e, g, ["Low", "High"], ["#4C78A8", "#E45756"], f"{k[0]} — {ep}", p)
        for a in axes.flat[len(keys):]:
            a.axis("off")
        fig.tight_layout()
        savefig(fig, P["fig"] / f"04_KM_grid_{ep}")
    meta = pd.DataFrame(metas)
    meta.to_csv(P["tab"] / "04_score_meta_analysis.csv", index=False)
    return res, meta


def mod05_clinical(expr, clin, genes, P, cfg):
    LOG.info("[05] UALCAN-style clinical subgroup analyses")
    vars_ = {"stage": ["Stage I", "Stage II", "Stage III", "Stage IV"],
             "grade": ["G1", "G2", "G3", "G4", "Low grade", "High grade"],
             "gender": None, "age_group": ["<=60", ">60"], "race": None}
    feats = ["score"] + list(genes)
    rows = []
    for ca, sub in clin.groupby("cancer"):
        for var in vars_:
            v = sub[var].dropna()
            vc = v.value_counts()
            lv = vc[vc >= cfg.min_group].index
            if len(lv) < 2:
                continue
            idx = v[v.isin(lv)].index
            for f in feats:
                x = sub.loc[idx, "score"] if f == "score" else expr.loc[f, idx]
                groups = [x[v[idx] == l].dropna().values for l in lv]
                try:
                    p = stats.kruskal(*groups).pvalue if len(lv) > 2 else \
                        stats.mannwhitneyu(*groups).pvalue
                except ValueError:
                    continue
                rho_trend = np.nan
                if var in ("stage", "grade") and vars_[var]:
                    order = {l: i for i, l in enumerate(vars_[var])}
                    ok = v[idx].map(order)
                    if ok.nunique() > 2:
                        rho_trend = stats.spearmanr(ok, x).correlation
                rows.append(dict(cancer=ca, variable=var, feature=f, levels="|".join(map(str, lv)),
                                 n=len(idx), p=p, spearman_trend=rho_trend,
                                 medians="|".join(f"{l}:{np.median(gv):.3f}" for l, gv in zip(lv, groups))))
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["FDR"] = res.groupby(["variable", res["feature"] == "score"])["p"].transform(lambda s: bh(s.values))
    res.to_csv(P["tab"] / "05_clinical_subgroups.csv", index=False)
    # heatmap of -log10 FDR for score
    s = res[res["feature"] == "score"].pivot(index="variable", columns="cancer", values="FDR")
    if len(s):
        fig, ax = plt.subplots(figsize=(0.4 * s.shape[1] + 2.5, 2.4))
        sns.heatmap(-np.log10(s), cmap="Reds", ax=ax, linewidths=.3, linecolor="white",
                    cbar_kws={"label": "−log10 FDR"}, annot=s.apply(lambda c: c.map(stars)),
                    fmt="", annot_kws={"size": 7})
        ax.set_title("DRG score vs clinical subgroups (UALCAN-style)", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("")
        savefig(fig, P["fig"] / "05_score_clinical_heatmap")
    # boxplots: stage for cancers with FDR<0.05
    for var in ("stage", "grade"):
        sig = res[(res["feature"] == "score") & (res["variable"] == var) & (res["FDR"] < 0.05)]["cancer"].tolist()
        if not sig:
            continue
        nc = min(4, len(sig))
        nr = int(np.ceil(len(sig) / nc))
        fig, axes = plt.subplots(nr, nc, figsize=(nc * 2.8, nr * 2.6), squeeze=False)
        for a, ca in zip(axes.flat, sig):
            d = clin[(clin["cancer"] == ca)].dropna(subset=[var])
            order = [l for l in vars_[var] if (d[var] == l).sum() >= cfg.min_group]
            sns.boxplot(data=d[d[var].isin(order)], x=var, y="score", order=order, ax=a,
                        color="#F2B701", fliersize=1, linewidth=.6)
            q = res[(res.cancer == ca) & (res.variable == var) & (res.feature == "score")]["FDR"].iloc[0]
            a.set_title(f"{ca}  FDR={q:.2g}", fontsize=8)
            a.set_xlabel("")
            a.tick_params(axis="x", rotation=45, labelsize=7)
        for a in axes.flat[len(sig):]:
            a.axis("off")
        fig.tight_layout()
        savefig(fig, P["fig"] / f"05_score_by_{var}_significant")
    return res


def mod06_metabolism(expr, clin, ann, genes, P, cfg, extra_sets):
    LOG.info("[06] Metabolic programme analysis")
    drg = set(genes)
    sets = {k: [g for g in v if g not in drg] for k, v in {**METABOLIC_SETS, **extra_sets}.items()}
    tum = clin.index
    met = ssgsea(expr[tum], sets, min_size=4).T          # samples x pathways
    met.to_csv(P["tab"] / "06_metabolic_ssgsea_scores_tumours.csv.gz")
    sizes = {k: len([g for g in v if g in expr.index]) for k, v in sets.items()}
    pd.Series(sizes, name="n_genes_found").to_csv(P["tab"] / "06_metabolic_set_sizes.csv")

    score_cols = [c for c in SCORE_LABELS if c in clin.columns]
    has_pur = "purity" in clin.columns and clin["purity"].notna().sum() > 100
    rows = []
    for ca, sub in clin.groupby("cancer"):
        if len(sub) < cfg.min_samples:
            continue
        pur_ok = has_pur and sub["purity"].notna().sum() >= cfg.min_samples
        sp = sub.index[sub["purity"].notna()] if pur_ok else None
        feats = {SCORE_LABELS[c]: sub[c] for c in score_cols}
        feats.update({g: expr.loc[g, sub.index] for g in genes})
        for pw in met.columns:
            mv = met.loc[sub.index, pw]
            for f, v in feats.items():
                rho, p = stats.spearmanr(v, mv)
                row = dict(cancer=ca, feature=f, pathway=pw, rho=rho, p=p, n=len(sub))
                if pur_ok:
                    row["rho_purity_adj"], row["p_purity_adj"] = partial_spearman(
                        v[sp], mv[sp], sub.loc[sp, "purity"])
                    row["n_purity"] = len(sp)
                rows.append(row)
    res = pd.DataFrame(rows)
    res["FDR"] = bh(res["p"])
    if "p_purity_adj" in res:
        res["FDR_purity_adj"] = bh(res["p_purity_adj"])
    res.to_csv(P["tab"] / "06_metabolic_correlations.csv.gz", index=False)
    if "rho_purity_adj" in res:
        s = res[res["feature"] == "DRG_score"].dropna(subset=["rho_purity_adj"])
        R = s.pivot(index="pathway", columns="cancer", values="rho_purity_adj")
        Q = s.pivot(index="pathway", columns="cancer", values="FDR_purity_adj")
        fig, ax = plt.subplots(figsize=(0.42 * R.shape[1] + 3, 0.34 * R.shape[0] + 1.5))
        sns.heatmap(R, cmap="RdBu_r", center=0, vmin=-.8, vmax=.8, ax=ax, linewidths=.3,
                    linecolor="white", annot=Q.apply(lambda c: c.map(stars)), fmt="",
                    annot_kws={"size": 6}, cbar_kws={"label": "partial Spearman ρ | purity"})
        ax.set_title("DRG score vs metabolic programmes — adjusted for tumour purity (ABSOLUTE)", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("")
        savefig(fig, P["fig"] / "06_score_metabolic_correlation_heatmap_purity_adjusted")
        both = s.groupby("pathway")[["rho", "rho_purity_adj"]].median()
        fig, ax = plt.subplots(figsize=(4.2, 4))
        ax.scatter(both["rho"], both["rho_purity_adj"], s=18, c="#E45756")
        for pw, r in both.iterrows():
            ax.annotate(pw.replace("_", " "), (r["rho"], r["rho_purity_adj"]), fontsize=5.5,
                        xytext=(3, 2), textcoords="offset points")
        lo = min(both.min().min(), 0) - .05
        hi = both.max().max() + .05
        ax.plot([lo, hi], [lo, hi], "k--", lw=.6)
        ax.axhline(0, c="grey", lw=.4)
        ax.axvline(0, c="grey", lw=.4)
        ax.set_xlabel("median ρ (unadjusted)")
        ax.set_ylabel("median partial ρ (purity-adjusted)")
        ax.set_title("Robustness of metabolic coupling to tumour purity", fontsize=8)
        savefig(fig, P["fig"] / "06_purity_adjustment_robustness")

    s = res[res["feature"] == "DRG_score"]
    R = s.pivot(index="pathway", columns="cancer", values="rho")
    Q = s.pivot(index="pathway", columns="cancer", values="FDR")
    fig, ax = plt.subplots(figsize=(0.42 * R.shape[1] + 3, 0.34 * R.shape[0] + 1.5))
    sns.heatmap(R, cmap="RdBu_r", center=0, vmin=-.8, vmax=.8, ax=ax, linewidths=.3,
                linecolor="white", annot=Q.apply(lambda c: c.map(stars)), fmt="",
                annot_kws={"size": 6}, cbar_kws={"label": "Spearman ρ"})
    ax.set_title("DRG score vs metabolic programmes (tumours; DRGs excluded from sets)", fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")
    savefig(fig, P["fig"] / "06_score_metabolic_correlation_heatmap")

    gsum = res[~res["feature"].isin(list(SCORE_LABELS.values()))].groupby(["feature", "pathway"])["rho"].median().unstack()
    gsum = gsum.reindex(genes)
    fig, ax = plt.subplots(figsize=(0.4 * gsum.shape[1] + 3, 0.3 * gsum.shape[0] + 2))
    sns.heatmap(gsum, cmap="RdBu_r", center=0, vmin=-.6, vmax=.6, ax=ax, linewidths=.3,
                linecolor="white", cbar_kws={"label": "median Spearman ρ across cancers"})
    ax.set_title("DRG genes vs metabolic programmes (pan-cancer median)", fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")
    savefig(fig, P["fig"] / "06_gene_metabolic_median_rho")

    # High vs low score: pathway shift (Cliff's-delta-like via median difference)
    hl = []
    for ca, sub in clin.groupby("cancer"):
        if len(sub) < cfg.min_samples:
            continue
        hi = sub.index[sub["score"] > sub["score"].median()]
        lo = sub.index[sub["score"] <= sub["score"].median()]
        for pw in met.columns:
            a, b = met.loc[hi, pw], met.loc[lo, pw]
            p = stats.mannwhitneyu(a, b).pvalue
            d = (np.sign(a.values[:, None] - b.values[None, :])).mean()
            hl.append(dict(cancer=ca, pathway=pw, cliffs_delta=d, p=p))
    hl = pd.DataFrame(hl)
    hl["FDR"] = bh(hl["p"])
    hl.to_csv(P["tab"] / "06_high_vs_low_score_metabolic.csv", index=False)

    # Metabolic pathways prognostic Cox (OS) + interaction with DRG score
    cx = []
    for ca, sub in clin.groupby("cancer"):
        d = sub[["score", "OS", "OS.time"]].dropna()
        d = d[d["OS.time"] > 0]
        if len(d) < cfg.min_samples or d["OS"].sum() < cfg.min_events:
            continue
        for pw in met.columns:
            X = pd.DataFrame({"score": zscore(d["score"]), "pathway": zscore(met.loc[d.index, pw])},
                             index=d.index)
            X["interaction"] = X["score"] * X["pathway"]
            try:
                r = cox_fit(d["OS.time"], d["OS"], X)
                cx.append(dict(cancer=ca, pathway=pw, HR_pathway=r.loc["pathway", "HR"],
                               p_pathway=r.loc["pathway", "p"], HR_interaction=r.loc["interaction", "HR"],
                               p_interaction=r.loc["interaction", "p"]))
            except Exception:  # noqa: BLE001
                pass
    cx = pd.DataFrame(cx)
    if len(cx):
        cx["FDR_interaction"] = bh(cx["p_interaction"])
        cx.to_csv(P["tab"] / "06_score_x_metabolism_interaction_cox_OS.csv", index=False)
    return res, hl


def mod07_multivariable(clin, P, cfg):
    LOG.info("[07] Multivariable Cox (score + age + gender + stage [+ purity])")
    rows = []
    models = ["clinical"] + (["clinical+purity"] if "purity" in clin.columns
                             and clin["purity"].notna().sum() > 100 else [])
    for ca, sub in clin.groupby("cancer"):
      for model in models:
        for ep in cfg.endpoints:
            cols = ["score", "age", "gender", "stage", ep, f"{ep}.time"] + \
                (["purity"] if model == "clinical+purity" else [])
            d = sub[cols].copy()
            d = d[d[f"{ep}.time"] > 0]
            X = pd.DataFrame(index=d.index)
            X["DRG_score"] = zscore(d["score"])
            X["age_per10y"] = d["age"] / 10
            if d["gender"].nunique() > 1:
                X["male"] = (d["gender"] == "MALE").astype(float)
            stg = d["stage"].map({"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
            if stg.notna().mean() > 0.6 and stg.nunique() > 1:
                X["stage_num"] = stg
            if model == "clinical+purity":
                X["purity_per10pct"] = d["purity"] * 10
            dd = pd.concat([X, d[[ep, f"{ep}.time"]]], axis=1).dropna()
            if len(dd) < cfg.min_samples or dd[ep].sum() < cfg.min_events:
                continue
            dd = dd.loc[:, dd.nunique() > 1]
            covs = [c for c in X.columns if c in dd.columns]
            if "DRG_score" not in covs:
                continue
            try:
                r = cox_fit(dd[f"{ep}.time"], dd[ep], dd[covs])
            except Exception:  # noqa: BLE001
                continue
            for cov, rr in r.iterrows():
                rows.append(dict(cancer=ca, model=model, endpoint=ep, covariate=cov, n=len(dd),
                                 events=int(dd[ep].sum()), HR=rr.HR, HR_low=rr.HR_low,
                                 HR_high=rr.HR_high, p=rr.p))
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res.to_csv(P["tab"] / "07_multivariable_cox.csv", index=False)
    s = res[(res["covariate"] == "DRG_score") & (res["endpoint"] == "OS") &
            (res["model"] == "clinical")].sort_values("HR")
    if len(s):
        fig, ax = plt.subplots(figsize=(5, 0.24 * len(s) + 1.2))
        y = np.arange(len(s))
        col = np.where(s["p"] < 0.05, np.where(s["HR"] > 1, "#E45756", "#4C78A8"), "grey")
        ax.errorbar(s["HR"], y, xerr=[s["HR"] - s["HR_low"], s["HR_high"] - s["HR"]], fmt="none", ecolor=col)
        ax.scatter(s["HR"], y, c=col, s=14, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(s["cancer"], fontsize=7)
        ax.axvline(1, ls="--", c="k", lw=.6)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("Adjusted HR per SD (age, gender, stage)")
        ax.set_title("DRG score as independent OS factor", fontsize=9)
        savefig(fig, P["fig"] / "07_multivariable_forest_OS")
    return res


def _lasso_fit_cv(lc, t, e, X, rng, nfold, nlam, rule):
    """Choose lambda by k-fold CV C-index on (t, e, X); return coefficients."""
    lmax = lc.lambda_max(t, e, X)
    lams = lmax * np.logspace(0, -2, nlam)
    folds = np.array_split(rng.permutation(len(t)), nfold)
    cv = np.full((nfold, nlam), np.nan)
    for k, f in enumerate(folds):
        fit = np.setdiff1d(np.arange(len(t)), f)
        if e[f].sum() < 2:
            continue
        for j, b in enumerate(lc.path(t[fit], e[fit], X[fit], lams)):
            cv[k, j] = cindex(t[f], e[f], X[f] @ b) if np.any(b) else 0.5
    m = np.nanmean(cv, 0)
    sd = np.nanstd(cv, 0) / np.sqrt(np.sum(~np.isnan(cv), 0).clip(1))
    jb = int(np.nanargmax(m))
    j1 = int(np.where(m >= m[jb] - sd[jb])[0].min())
    path = lc.path(t, e, X, lams)
    b = path[j1] if rule == "1se" else path[jb]
    if not np.any(b):
        b = path[jb]
    return b, float(m[jb])


def mod08_lasso_model(expr, clin, genes, P, cfg):
    LOG.info(f"[08] LASSO-Cox DRG models — {cfg.n_splits} repeated stratified 70/30 splits")
    rng = np.random.default_rng(cfg.seed)
    ep = cfg.model_endpoint
    lc = LassoCox(max_iter=1500, tol=1e-6)
    perf, splits, coefs, km = [], [], [], {}
    for ca, sub in clin.groupby("cancer"):
        d = sub[[ep, f"{ep}.time"]].dropna()
        d = d[d[f"{ep}.time"] > 0]
        if len(d) < max(60, cfg.min_samples) or d[ep].sum() < cfg.min_events_model:
            continue
        Xall = expr.loc[genes, d.index].T.astype(float)
        Xall = Xall.loc[:, Xall.std() > 0].dropna(axis=1)
        gnames = list(Xall.columns)
        t, e = d[f"{ep}.time"].values, d[ep].values.astype(int)
        idx = np.arange(len(d))
        sel = np.zeros(len(gnames))
        oos_risk = [[] for _ in idx]                      # out-of-sample risk percentiles
        for rep_ in range(cfg.n_splits):
            tr = np.concatenate([rng.choice(idx[e == k], int(round(0.7 * (e == k).sum())), replace=False)
                                 for k in (0, 1)])
            te = np.setdiff1d(idx, tr)
            mu, sdv = Xall.iloc[tr].mean(), Xall.iloc[tr].std()
            X = ((Xall - mu) / sdv).values
            b, _ = _lasso_fit_cv(lc, t[tr], e[tr], X[tr], rng, 3, 15, cfg.lambda_rule)
            sel += b != 0
            if not np.any(b):
                splits.append(dict(cancer=ca, split=rep_, n_genes=0, C_train=np.nan, C_test=0.5))
                continue
            rtr, rte = X[tr] @ b, X[te] @ b
            cut = np.median(rtr)
            row = dict(cancer=ca, split=rep_, n_genes=int((b != 0).sum()),
                       C_train=cindex(t[tr], e[tr], rtr), C_test=cindex(t[te], e[te], rte),
                       logrank_test=logrank(t[te], e[te], (rte > cut).astype(int)))
            for yr in (1, 3, 5):
                row[f"AUC_{yr}y_test"] = td_auc(t[te], e[te], rte, yr)
            splits.append(row)
            for i, r in zip(te, rte > cut):
                oos_risk[i].append(int(r))
        sp = pd.DataFrame([r for r in splits if r["cancer"] == ca])
        ct = sp["C_test"].dropna()
        # final model on the full cohort (for coefficients / nomogram use)
        mu, sdv = Xall.mean(), Xall.std()
        Xf = ((Xall - mu) / sdv).values
        bf, cvC = _lasso_fit_cv(lc, t, e, Xf, rng, 5, 25, cfg.lambda_rule)
        for g, c, f in zip(gnames, bf, sel / cfg.n_splits):
            coefs.append(dict(cancer=ca, gene=g, coef_full_model=c, selection_freq=f))
        # out-of-sample risk group = majority vote across splits where the patient was in test
        grp = np.array([np.mean(v) > 0.5 if v else np.nan for v in oos_risk], dtype=float)
        ok = ~np.isnan(grp)
        p_oos = logrank(t[ok], e[ok], grp[ok].astype(int))
        km[ca] = (t[ok], e[ok], grp[ok].astype(int), p_oos)
        perf.append(dict(cancer=ca, endpoint=ep, n=len(d), events=int(e.sum()),
                         n_splits=cfg.n_splits, C_test_median=ct.median(),
                         C_test_q025=ct.quantile(.025), C_test_q975=ct.quantile(.975),
                         C_train_median=sp["C_train"].median(),
                         optimism=sp["C_train"].median() - ct.median(),
                         pct_splits_C_gt_0_6=(ct > 0.6).mean() * 100,
                         AUC_1y_median=sp.get("AUC_1y_test", pd.Series(dtype=float)).median(),
                         AUC_3y_median=sp.get("AUC_3y_test", pd.Series(dtype=float)).median(),
                         AUC_5y_median=sp.get("AUC_5y_test", pd.Series(dtype=float)).median(),
                         logrank_oos=p_oos, n_genes_full=int((bf != 0).sum()), cv_C_full=cvC,
                         stable_genes=";".join(g for g, f in zip(gnames, sel / cfg.n_splits) if f >= 0.6)))
        LOG.info(f"  {ca}: median test C = {ct.median():.3f} "
                 f"[{ct.quantile(.025):.2f}-{ct.quantile(.975):.2f}]")
    perf = pd.DataFrame(perf)
    if perf.empty:
        return perf, pd.DataFrame()
    perf["robust"] = (perf["C_test_q025"] > 0.5) & (perf["C_test_median"] >= 0.6) & (perf["logrank_oos"] < 0.05)
    perf = perf.sort_values("C_test_median", ascending=False)
    coefs = pd.DataFrame(coefs)
    perf.to_csv(P["tab"] / f"08_lasso_repeated_split_performance_{ep}.csv", index=False)
    pd.DataFrame(splits).to_csv(P["tab"] / f"08_lasso_all_splits_{ep}.csv.gz", index=False)
    coefs.to_csv(P["tab"] / f"08_lasso_coefficients_and_stability_{ep}.csv", index=False)

    # Figure: test C-index distribution per cancer
    sp = pd.DataFrame(splits)
    order = perf["cancer"].tolist()
    fig, ax = plt.subplots(figsize=(0.42 * len(order) + 2, 3.4))
    pal = {c: ("#E45756" if r else "#BBBBBB") for c, r in zip(perf["cancer"], perf["robust"])}
    sns.boxplot(data=sp, x="cancer", y="C_test", order=order, palette=pal, fliersize=1, linewidth=.6, ax=ax)
    ax.axhline(0.5, ls="--", c="k", lw=.6)
    ax.set_ylabel(f"Test-set C-index ({cfg.n_splits} splits)")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.set_title("DRG LASSO-Cox models: out-of-sample discrimination (red = robust)", fontsize=9)
    savefig(fig, P["fig"] / f"08_lasso_test_Cindex_distribution_{ep}")
    # Figure: selection stability (robust cancers)
    rob = perf.loc[perf["robust"], "cancer"].tolist()
    if rob:
        sm = coefs[coefs["cancer"].isin(rob)].pivot(index="gene", columns="cancer", values="selection_freq")
        sm = sm.reindex([g for g in genes if g in sm.index])
        fig, ax = plt.subplots(figsize=(0.5 * sm.shape[1] + 3, 0.3 * sm.shape[0] + 1.5))
        sns.heatmap(sm, cmap="Greens", vmin=0, vmax=1, ax=ax, annot=True, fmt=".2f",
                    annot_kws={"size": 5.5}, linewidths=.3, linecolor="white",
                    cbar_kws={"label": "selection frequency"})
        ax.set_title("Gene selection stability across splits (robust models)", fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("")
        savefig(fig, P["fig"] / f"08_lasso_gene_stability_{ep}")
        nc = min(5, len(rob))
        nr = int(np.ceil(len(rob) / nc))
        fig, axes = plt.subplots(nr, nc, figsize=(nc * 2.7, nr * 2.5), squeeze=False)
        for a, k in zip(axes.flat, rob):
            tt, ee, g, p = km[k]
            km_plot(a, tt, ee, g, ["Low risk", "High risk"], ["#4C78A8", "#E45756"],
                    f"{k} out-of-sample", p)
        for a in axes.flat[len(rob):]:
            a.axis("off")
        fig.tight_layout()
        savefig(fig, P["fig"] / f"08_lasso_KM_out_of_sample_{ep}")
    return perf, coefs


def load_xena_matrix(path, genes, samples):
    keep = set(genes)
    parts = []
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=3000):
        ch.index = ch.index.astype(str).str.split("|").str[0]
        parts.append(ch[ch.index.isin(keep)])
    m = pd.concat(parts)
    m = m[~m.index.duplicated()]
    m.columns = m.columns.astype(str).str[:15]
    m = m.loc[:, ~m.columns.duplicated()]
    return m.reindex(columns=[s for s in samples if s in m.columns])


def mod09_genomics(paths, expr, clin, genes, P, cfg):
    LOG.info("[09] Genomic alterations (mutation / CNV)")
    out = {}
    tum = clin.index.tolist()
    if paths.get("mutation") and paths["mutation"].exists():
        try:
            mut = load_xena_matrix(paths["mutation"], genes, tum)
            cc = clin.loc[mut.columns, "cancer"]
            freq = mut.T.groupby(cc).mean().T * 100
            freq.to_csv(P["tab"] / "09_mutation_frequency_pct.csv")
            fig, ax = plt.subplots(figsize=(0.4 * freq.shape[1] + 3, 0.3 * freq.shape[0] + 1.5))
            sns.heatmap(freq.reindex([g for g in genes if g in freq.index]), cmap="Reds", ax=ax,
                        annot=True, fmt=".0f", annot_kws={"size": 5}, cbar_kws={"label": "% mutated"})
            ax.set_title("Non-silent somatic mutation frequency (MC3)", fontsize=9)
            ax.set_xlabel("")
            ax.set_ylabel("")
            savefig(fig, P["fig"] / "09_mutation_frequency")
            out["mutation"] = freq
        except Exception as ex:  # noqa: BLE001
            LOG.warning(f"  mutation module skipped: {ex}")
    if paths.get("cnv") and paths["cnv"].exists():
        try:
            cnv = load_xena_matrix(paths["cnv"], genes, tum)
            cc = clin.loc[cnv.columns, "cancer"]
            vals = pd.Series(cnv.values.ravel()).dropna()
            thresholded = vals.isin([-2, -1, 0, 1, 2]).mean() > 0.99
            if thresholded:
                LOG.info("  CNV input: GISTIC2 thresholded (-2..2)")
                g_, l_, a_, d_ = 1, -1, 2, -2
                cnv_kind = "GISTIC2 thresholded calls (gain >=1, loss <=-1; high-level ±2)"
            else:
                LOG.info("  CNV input: continuous gene-level log2(CN/2) — applying thresholds")
                g_, l_, a_, d_ = CNV_GAIN, CNV_LOSS, CNV_AMP, CNV_DEL
                cnv_kind = (f"continuous gene-level log2 copy-number ratio (gain > {CNV_GAIN}, "
                            f"loss < {CNV_LOSS}; high-level amplification > {CNV_AMP}, deep deletion < {CNV_DEL})")
            pd.Series({"cnv_input": cnv_kind}).to_csv(P["tab"] / "09_cnv_input_type.csv")
            amp = (cnv >= g_ if thresholded else cnv > g_).T.groupby(cc).mean().T * 100
            dele = (cnv <= l_ if thresholded else cnv < l_).T.groupby(cc).mean().T * 100
            hamp = (cnv >= a_ if thresholded else cnv > a_).T.groupby(cc).mean().T * 100
            hdel = (cnv <= d_ if thresholded else cnv < d_).T.groupby(cc).mean().T * 100
            amp.to_csv(P["tab"] / "09_cnv_gain_pct.csv")
            dele.to_csv(P["tab"] / "09_cnv_loss_pct.csv")
            hamp.to_csv(P["tab"] / "09_cnv_highlevel_amplification_pct.csv")
            hdel.to_csv(P["tab"] / "09_cnv_deep_deletion_pct.csv")
            rows = []
            for ca in cc.unique():
                s = cc.index[cc == ca]
                for g in cnv.index:
                    if g in expr.index:
                        r, p = stats.spearmanr(cnv.loc[g, s], expr.loc[g, s], nan_policy="omit")
                        rows.append(dict(cancer=ca, gene=g, rho=r, p=p))
            ce = pd.DataFrame(rows)
            ce["FDR"] = bh(ce["p"])
            ce.to_csv(P["tab"] / "09_cnv_expression_correlation.csv", index=False)
            fig, axes = plt.subplots(1, 2, figsize=(2 * (0.35 * amp.shape[1] + 2), 0.3 * amp.shape[0] + 1.5))
            order = [g for g in genes if g in amp.index]
            for a, mtx, ttl, cm in ((axes[0], amp, "CNV gain (%)", "Reds"), (axes[1], dele, "CNV loss (%)", "Blues")):
                sns.heatmap(mtx.reindex(order), cmap=cm, ax=a, cbar_kws={"label": "%"})
                a.set_title(ttl, fontsize=9)
                a.set_xlabel("")
                a.set_ylabel("")
            fig.suptitle(cnv_kind, fontsize=7, y=1.0)
            savefig(fig, P["fig"] / "09_cnv_frequency")
            out["cnv"] = (amp, dele)
        except Exception as ex:  # noqa: BLE001
            LOG.warning(f"  CNV module skipped: {ex}")
    return out


def mod10_coexpression(expr, clin, genes, P):
    LOG.info("[10] DRG co-expression")
    mats = []
    for ca, sub in clin.groupby("cancer"):
        if len(sub) >= 30:
            mats.append(expr.loc[genes, sub.index].T.corr(method="spearman").values)
    med = pd.DataFrame(np.nanmedian(np.stack(mats), axis=0), index=genes, columns=genes)
    med.to_csv(P["tab"] / "10_DRG_coexpression_median_rho.csv")
    cats = pd.Series(genes).map(GENE_CATEGORY)
    pal = dict(zip(cats.unique(), ["#E45756", "#72B7B2", "#F2B701"]))
    g = sns.clustermap(med.fillna(0), cmap="RdBu_r", center=0, vmin=-.8, vmax=.8,
                       row_colors=cats.map(pal).values, figsize=(8, 8),
                       cbar_kws={"label": "median Spearman ρ"})
    g.fig.suptitle("DRG co-expression (median across cancers)", y=1.01, fontsize=10)
    g.savefig(P["fig"] / "10_DRG_coexpression_clustermap.png", dpi=300, bbox_inches="tight")
    g.savefig(P["fig"] / "10_DRG_coexpression_clustermap.pdf", bbox_inches="tight")
    plt.close("all")
    return med


def mod11_subscores(clin, ann, scores_all, P, cfg):
    LOG.info("[11] Sub-scores: core screen arm vs actin-cytoskeleton arm")
    feats = [("DRG_score", "score"), ("DRG_core", "score_core"), ("DRG_cytoskeleton", "score_cyto")]
    feats = [(f, c) for f, c in feats if c in scores_all.columns]
    rows = []
    tn = ann[ann["group"].isin(["Tumor", "Normal"])]
    for f, col in feats:
        for ca, sub in tn.groupby("cancer"):
            ti, ni = sub.index[sub.group == "Tumor"], sub.index[sub.group == "Normal"]
            if len(ni) < cfg.min_normals or len(ti) < cfg.min_normals:
                continue
            a, b = scores_all.loc[ti, col], scores_all.loc[ni, col]
            rows.append(dict(feature=f, analysis="tumour_vs_normal", cancer=ca, n=len(ti) + len(ni),
                             effect=a.mean() - b.mean(), p=stats.mannwhitneyu(a, b).pvalue))
        for ca, sub in clin.groupby("cancer"):
            for ep in [e for e in ("OS", "DSS", "PFI") if e in cfg.endpoints]:
                d = sub[[col, ep, f"{ep}.time"]].dropna()
                d = d[d[f"{ep}.time"] > 0]
                if len(d) < cfg.min_samples or d[ep].sum() < cfg.min_events:
                    continue
                try:
                    r = cox_fit(d[f"{ep}.time"], d[ep], pd.DataFrame({"x": zscore(d[col])}, index=d.index)).iloc[0]
                except Exception:  # noqa: BLE001
                    continue
                rows.append(dict(feature=f, analysis=f"cox_{ep}", cancer=ca, n=len(d), effect=r.HR,
                                 low=r.HR_low, high=r.HR_high, p=r.p))
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["FDR"] = res.groupby(["feature", "analysis"])["p"].transform(lambda x: bh(x.values))
    res.to_csv(P["tab"] / "11_subscores_tumour_normal_and_survival.csv", index=False)
    cor = clin.groupby("cancer").apply(lambda d: stats.spearmanr(d["score_core"], d["score_cyto"]).correlation
                                       if len(d) > 10 else np.nan).rename("rho_core_vs_cytoskeleton")
    cor.to_csv(P["tab"] / "11_core_vs_cytoskeleton_correlation.csv")
    LOG.info(f"  core vs cytoskeleton sub-score: median rho across cancers = {cor.median():.2f}")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 0.24 * clin.cancer.nunique() + 1.8),
                             gridspec_kw={"width_ratios": [1, 1]})
    for a, an, lab, tr in ((axes[0], "tumour_vs_normal", "Δ score (tumour − normal)", lambda v: v),
                           (axes[1], "cox_OS", "log2 HR per SD (OS)", np.log2)):
        r = res[res.analysis == an]
        if r.empty:
            a.axis("off")
            continue
        M = tr(r.pivot(index="cancer", columns="feature", values="effect"))[[f for f, _ in feats]]
        Q = r.pivot(index="cancer", columns="feature", values="FDR")[[f for f, _ in feats]]
        lim = np.nanpercentile(np.abs(M.values), 98)
        sns.heatmap(M, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, ax=a, linewidths=.3, linecolor="white",
                    annot=Q.apply(lambda c: c.map(stars)), fmt="", annot_kws={"size": 6},
                    cbar_kws={"label": lab, "shrink": .6})
        a.set_xlabel("")
        a.set_ylabel("")
        a.tick_params(axis="x", rotation=30, labelsize=7)
        a.tick_params(axis="y", labelsize=7)
        a.set_title("Tumour vs normal" if an == "tumour_vs_normal" else "Overall survival", fontsize=9)
    fig.suptitle("Full DRG score vs its two biological arms (BH-FDR * <.05, ** <.01, *** <.001)", fontsize=9)
    fig.tight_layout()
    savefig(fig, P["fig"] / "11_subscores_core_vs_cytoskeleton")
    return res


# ════════════════════════════════════════════════════════════════════════════
# 8. REPORT
# ════════════════════════════════════════════════════════════════════════════
def write_report(P, cfg, genes, missing, ann, clin, tn, sc, meta, mv, met, perf, t0, subres=None):
    L = []
    L.append("# Disulfidptosis-related gene signature — pan-cancer analysis report\n")
    L.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M')} | runtime {(time.time() - t0) / 60:.1f} min\n")
    L.append("## Data\n")
    L.append(f"- Signature `{cfg.signature}`: {len(genes)} genes found "
             f"({', '.join(genes)})" + (f"; missing: {', '.join(missing)}" if missing else ""))
    L.append(f"- Tumour samples (one per patient): {len(clin)} across {clin['cancer'].nunique()} cancers")
    L.append(f"- Normal samples: {(ann['group'] == 'Normal').sum()}\n")
    if tn is not None and len(tn):
        s = tn[tn["gene"] == "DRG_score"]
        up = s[(s.FDR < .05) & (s.log2FC > 0)]["cancer"].tolist()
        dn = s[(s.FDR < .05) & (s.log2FC < 0)]["cancer"].tolist()
        L.append("## Tumour vs normal\n")
        L.append(f"- DRG score higher in tumour: {', '.join(up) or 'none'}")
        L.append(f"- DRG score lower in tumour: {', '.join(dn) or 'none'}")
        g = tn[(tn.gene != "DRG_score") & (tn.FDR < .05)]
        top = g.assign(sgn=np.sign(g.log2FC)).groupby("gene")["sgn"].agg(
            up=lambda x: (x > 0).sum(), down=lambda x: (x < 0).sum()).sort_values("up", ascending=False)
        L.append("\n| gene | cancers up | cancers down |\n|---|---|---|")
        L += [f"| {i} | {r.up} | {r.down} |" for i, r in top.iterrows()]
        L.append("")
    if sc is not None and len(sc):
        L.append("## Prognosis (DRG score, HR per SD)\n")
        for ep in cfg.endpoints:
            s = sc[sc.endpoint == ep]
            if s.empty:
                continue
            risk = s[(s.p_cox < .05) & (s.HR_perSD > 1)]["cancer"].tolist()
            prot = s[(s.p_cox < .05) & (s.HR_perSD < 1)]["cancer"].tolist()
            L.append(f"- **{ep}** — risk: {', '.join(risk) or 'none'}; protective: {', '.join(prot) or 'none'}")
        if meta is not None and len(meta):
            L.append("\n| endpoint | k | pooled HR (95% CI) | p | I² |\n|---|---|---|---|---|")
            L += [f"| {r.endpoint} | {r.k} | {r.HR:.2f} ({r.HR_low:.2f}–{r.HR_high:.2f}) | {r.p:.2g} | {r.I2:.0f}% |"
                  for r in meta.itertuples()]
        L.append("")
    if mv is not None and len(mv):
        for model, ttl in (("clinical", "age/gender/stage"), ("clinical+purity", "age/gender/stage + ABSOLUTE purity")):
            s = mv[(mv.covariate == "DRG_score") & (mv.endpoint == "OS") & (mv.p < .05) & (mv.model == model)]
            if model not in set(mv.model):
                continue
            L.append(f"## Independent prognostic value (OS, adjusted for {ttl})\n")
            L.append(", ".join(f"{r.cancer} (HR {r.HR:.2f}, p={r.p:.2g})" for r in s.itertuples()) or "none")
            L.append("")
    if met is not None and len(met):
        s = met[met.feature == "DRG_score"]
        agg = s.groupby("pathway").agg(median_rho=("rho", "median"),
                                       pos=("rho", lambda x: ((x > 0) & (s.loc[x.index, "FDR"] < .05)).sum()),
                                       neg=("rho", lambda x: ((x < 0) & (s.loc[x.index, "FDR"] < .05)).sum()))
        agg = agg.sort_values("median_rho", ascending=False)
        L.append("## Metabolic coupling of the DRG score (tumours)\n")
        L.append("| pathway | median ρ | cancers + (FDR<.05) | cancers − |\n|---|---|---|---|")
        L += [f"| {i} | {r.median_rho:.2f} | {int(r.pos)} | {int(r.neg)} |" for i, r in agg.iterrows()]
        L.append("")
        if "rho_purity_adj" in met:
            s = met[met.feature == "DRG_score"].dropna(subset=["rho_purity_adj"])
            agg = s.groupby("pathway").agg(
                rho=("rho", "median"), rho_adj=("rho_purity_adj", "median"),
                pos=("rho_purity_adj", lambda x: ((x > 0) & (s.loc[x.index, "FDR_purity_adj"] < .05)).sum()),
                neg=("rho_purity_adj", lambda x: ((x < 0) & (s.loc[x.index, "FDR_purity_adj"] < .05)).sum()))
            agg = agg.sort_values("rho_adj", ascending=False)
            L.append("### Purity-adjusted (partial Spearman | ABSOLUTE purity)\n")
            L.append("| pathway | median ρ raw | median ρ adjusted | cancers + | cancers − |\n|---|---|---|---|---|")
            L += [f"| {i} | {r.rho:.2f} | {r.rho_adj:.2f} | {int(r.pos)} | {int(r.neg)} |" for i, r in agg.iterrows()]
            L.append("")
        sub = met[met.feature.isin(["DRG_core", "DRG_cytoskeleton"])]
        if len(sub):
            t = sub.groupby(["pathway", "feature"])["rho"].median().unstack()
            L.append("### Sub-score coupling (median ρ)\n")
            L.append("| pathway | core screen arm | cytoskeleton arm |\n|---|---|---|")
            L += [f"| {i} | {r.get('DRG_core', np.nan):.2f} | {r.get('DRG_cytoskeleton', np.nan):.2f} |"
                  for i, r in t.iterrows()]
            L.append("")
    if perf is not None and len(perf) and "C_test_median" in perf:
        L.append(f"## LASSO-Cox models ({cfg.model_endpoint}) — {cfg.n_splits} repeated 70/30 splits\n")
        L.append("Robust = 2.5th percentile of test C > 0.5, median test C ≥ 0.6 and out-of-sample log-rank p < 0.05.\n")
        L.append("| cancer | robust | median test C [95% range] | optimism | AUC 1y | AUC 3y | AUC 5y | OOS log-rank p | stable genes (≥60%) |\n|---|---|---|---|---|---|---|---|---|")
        L += [f"| {r.cancer} | {'yes' if r.robust else ''} | {r.C_test_median:.2f} [{r.C_test_q025:.2f}–{r.C_test_q975:.2f}] | "
              f"{r.optimism:.2f} | {r.AUC_1y_median:.2f} | {r.AUC_3y_median:.2f} | {r.AUC_5y_median:.2f} | "
              f"{r.logrank_oos:.2g} | {r.stable_genes or '—'} |" for r in perf.itertuples()]
        L.append("")
    if subres is not None and len(subres):
        L.append("## Sub-scores: core screen arm vs actin-cytoskeleton arm (OS, HR per SD, p<0.05)\n")
        for f in ("DRG_core", "DRG_cytoskeleton"):
            s = subres[(subres.feature == f) & (subres.analysis == "cox_OS") & (subres.p < .05)]
            L.append(f"- **{f}** — risk: {', '.join(s[s.effect > 1].cancer) or 'none'}; "
                     f"protective: {', '.join(s[s.effect < 1].cancer) or 'none'}")
            s = subres[(subres.feature == f) & (subres.analysis == "tumour_vs_normal") & (subres.FDR < .05)]
            L.append(f"  - tumour > normal: {', '.join(s[s.effect > 0].cancer) or 'none'}; "
                     f"tumour < normal: {', '.join(s[s.effect < 0].cancer) or 'none'}")
        L.append("")
    L.append("## Notes / limitations\n")
    L.append("- UALCAN has no programmatic API; its TCGA-based analyses (tumour vs normal, stage, grade, "
             "gender, age, race) are reproduced here from the same TCGA data. CPTAC protein and promoter "
             "methylation panels should be exported manually from ualcan.path.uab.edu for selected genes.")
    L.append("- Optimal-cutpoint KM inflates type-I error; median split is the primary analysis.")
    L.append("- Time-dependent AUC is the naive cumulative/dynamic estimator (no IPCW).")
    if "purity" in clin.columns:
        L.append(f"- ABSOLUTE purity available for {clin['purity'].notna().sum()} / {len(clin)} tumours.")
    else:
        L.append("- Tumour purity file unavailable: purity-adjusted analyses were skipped.")
    (P["base_out"] / "REPORT.md").write_text("\n".join(L), encoding="utf-8")


def write_methods(P, cfg, genes, clin):
    txt = f"""## Methods (auto-generated draft — verify before submission)

**Data acquisition.** Batch-corrected RNA-seq expression (EB++ adjusted, log2[normalised count + 1]),
TCGA Clinical Data Resource survival endpoints and sample annotations of the TCGA Pan-Cancer Atlas
were obtained from the UCSC Xena PanCanAtlas hub. One primary tumour sample per patient
(sample type 01; 03 for LAML; 06 included for SKCM) was retained, yielding {len(clin)} tumours from
{clin['cancer'].nunique()} cancer types; solid-tissue normal samples (type 11) were used for
tumour–normal comparisons in cancers with ≥{cfg.min_normals} normals.

**Disulfidptosis-related genes.** A {len(genes)}-gene disulfidptosis-related gene (DRG) set was compiled from
the original description of disulfidptosis (Liu et al., Nat Cell Biol 2023), comprising genes identified in the
glucose-starvation CRISPR screen of SLC7A11-high cells and actin-cytoskeleton proteins undergoing aberrant
disulfide bonding.

**Signature scoring.** Per-sample DRG activity was quantified by single-sample gene set enrichment analysis
(ssGSEA; τ = 0.25, GSVA normalisation) computed over the whole transcriptome; the mean of within-sample
gene z-scores was used as a sensitivity score.

**Differential expression.** Tumour–normal differences were tested with two-sided Mann–Whitney U tests; log2
fold-change was the difference of mean log2 expression. P values were adjusted with the Benjamini–Hochberg
procedure.

**Survival analyses.** Overall survival (OS), disease-specific survival (DSS), disease-free interval (DFI) and
progression-free interval (PFI) were analysed. Univariate Cox proportional hazards models (Efron ties) were
fitted per cancer for each z-standardised gene and for the DRG score (HR per SD), in cancers with
≥{cfg.min_samples} patients and ≥{cfg.min_events} events. Kaplan–Meier curves were compared between high and low
score groups ({cfg.cutoff} cut-off) using log-rank tests. Cancer-specific log-HRs were pooled with a
DerSimonian–Laird random-effects model (I² for heterogeneity). Multivariable Cox models adjusted for age,
gender and AJCC stage (where available).

**Clinical subgroups.** Following the UALCAN framework, DRG score and gene expression were compared across
pathological stage, histological grade, gender, age (≤60/>60 years) and race using Kruskal–Wallis or
Mann–Whitney U tests (groups with ≥{cfg.min_group} samples), with Spearman correlation for ordinal trend.

**Sub-scores.** Because the DRG set combines genes from the glucose-starvation CRISPR screen and actin-cytoskeleton
proteins, ssGSEA sub-scores were additionally computed for the core screen arm and the cytoskeleton arm and analysed
for tumour–normal differences, survival and metabolic coupling.

**Metabolic analyses.** ssGSEA scores were computed for curated metabolic programmes (glycolysis, glucose
uptake, pentose phosphate/NADPH, TCA cycle, oxidative phosphorylation, glutathione, cysteine/sulfur,
thioredoxin, NRF2, glutamine, serine/one-carbon, fatty acid synthesis/oxidation, glycogen, hypoxia, mTORC1);
DRGs were removed from these sets to avoid circularity. Associations with the DRG score and individual DRGs were
assessed by Spearman correlation; metabolic shifts between score-high and score-low tumours by Mann–Whitney U
(Cliff's δ), and score × pathway interactions by Cox regression. To exclude confounding by tumour
cellularity, all score–pathway correlations were repeated as partial Spearman correlations controlling for
ABSOLUTE tumour purity, and multivariable Cox models were refitted with purity as an additional covariate.

**Prognostic modelling.** For each cancer with ≥{cfg.min_events_model} {cfg.model_endpoint} events, LASSO-penalised
Cox models of the DRGs were developed and validated over {cfg.n_splits} repeated stratified 70/30 train–test splits
(seed {cfg.seed}). Within each training set, genes were standardised and the penalty was selected by 3-fold
cross-validated C-index ({'1-SE' if cfg.lambda_rule == '1se' else 'minimum'} rule; nested, so test patients never informed
model selection). Performance was summarised as the median and 2.5–97.5th percentile of test-set Harrell's C-index,
time-dependent AUC at 1, 3 and 5 years and optimism (train − test C). Gene selection frequency across splits
quantified stability, and each patient's out-of-sample risk group (majority vote across the splits in which the
patient was in the test set) was compared by log-rank test. Models were considered robust when the 2.5th percentile
of test C exceeded 0.5, the median test C was ≥0.6 and the out-of-sample log-rank P was <0.05.

**Genomic alterations.** Non-silent somatic mutation frequencies were derived from the MC3 call set; gene-level
copy number was summarised as gain/loss frequencies (see 09_cnv_input_type.csv for the thresholds applied to the
input format) and correlated with expression by Spearman correlation.

**Software.** Python {platform.python_version()} (pandas, NumPy, SciPy, statsmodels, seaborn). Analyses
were performed with a custom, fully reproducible pipeline.
"""
    (P["base_out"] / "METHODS_draft.md").write_text(txt, encoding="utf-8")


# ════════════════════════════════════════════════════════════════════════════
# 9. DEMO DATA (synthetic — for smoke-testing only)
# ════════════════════════════════════════════════════════════════════════════
def make_demo_data(raw: Path, seed=1):
    rng = np.random.default_rng(seed)
    allg = sorted(set(SIGNATURES["extended"]) | {g for v in METABOLIC_SETS.values() for g in v})
    allg += [f"GENE{i}" for i in range(1500)]
    cancers = ["BRCA", "LUAD", "KIRC", "LIHC", "COAD", "GBM"]
    samples, rows = [], []
    k = 0
    for ca in cancers:
        n_t = 160
        n_n = 20 if ca != "GBM" else 3
        for i in range(n_t + n_n):
            k += 1
            pid = f"TCGA-{ca[:2]}-{k:04d}"
            code = "01" if i < n_t else "11"
            samples.append((f"{pid}-{code}", pid, ca, code))
    S = pd.DataFrame(samples, columns=["sample", "patient", "cancer", "code"])
    n = len(S)
    latent = rng.normal(size=n)
    X = rng.normal(6, 1.2, size=(len(allg), n))
    drg_i = [allg.index(g) for g in SIGNATURES["extended"]]
    X[drg_i] += 0.8 * latent + (S["code"] == "01").values * 0.7
    gly = [allg.index(g) for g in METABOLIC_SETS["Glycolysis"]]
    X[gly] += 0.5 * latent
    expr = pd.DataFrame(X.astype(np.float32), index=allg, columns=S["sample"])
    expr.index.name = "sample"
    with gzip.open(raw / DATA_FILES["expression"][1], "wt") as f:
        expr.to_csv(f, sep="\t", float_format="%.4f")
    tum = S[S.code == "01"].copy()
    lt = latent[S.code == "01"]
    haz = np.exp(0.5 * lt) * 0.15
    T = rng.exponential(1 / haz)
    C = rng.uniform(0.5, 10, size=len(T))
    surv = pd.DataFrame({"sample": tum["sample"].values, "_PATIENT": tum["patient"].values,
                         "cancer type abbreviation": tum["cancer"].values,
                         "age_at_initial_pathologic_diagnosis": rng.integers(30, 85, len(tum)),
                         "gender": rng.choice(["MALE", "FEMALE"], len(tum)),
                         "race": rng.choice(["WHITE", "ASIAN", "BLACK OR AFRICAN AMERICAN"], len(tum)),
                         "ajcc_pathologic_tumor_stage": np.where(lt > 0.5, "Stage IIIA",
                                                                 rng.choice(["Stage I", "Stage IIB", "Stage IV"], len(tum))),
                         "clinical_stage": np.nan,
                         "histological_grade": rng.choice(["G1", "G2", "G3"], len(tum))})
    for ep in ENDPOINTS:
        surv[ep] = (T <= C).astype(int)
        surv[f"{ep}.time"] = (np.minimum(T, C) * 365.25).round()
    norm = S[S.code == "11"]
    surv_n = pd.DataFrame({"sample": norm["sample"], "_PATIENT": norm["patient"],
                           "cancer type abbreviation": norm["cancer"]})
    pd.concat([surv, surv_n]).to_csv(raw / DATA_FILES["survival"][1], sep="\t", index=False)
    pd.DataFrame({"array": tum["sample"].values, "call status": "called",
                  "purity": rng.uniform(0.2, 0.95, len(tum)).round(2)}).to_csv(
        raw / DATA_FILES["purity"][1], sep="\t", index=False)
    drg = SIGNATURES["extended"]
    cn = pd.DataFrame(rng.normal(0, 0.25, (len(drg), len(tum))).round(4), index=drg, columns=tum["sample"])
    cn.index.name = "sample"
    with gzip.open(raw / DATA_FILES["cnv"][1], "wt") as f:
        cn.to_csv(f, sep="\t")
    mu = pd.DataFrame((rng.random((len(drg), len(tum))) < 0.03).astype(int), index=drg, columns=tum["sample"])
    mu.index.name = "sample"
    with gzip.open(raw / DATA_FILES["mutation"][1], "wt") as f:
        mu.to_csv(f, sep="\t")
    LOG.info(f"  demo data written: {expr.shape}")


# ════════════════════════════════════════════════════════════════════════════
# 10. MAIN
# ════════════════════════════════════════════════════════════════════════════
def read_gmt(path):
    sets = {}
    for line in Path(path).read_text().splitlines():
        p = line.strip().split("\t")
        if len(p) > 2:
            sets[p[0]] = p[2:]
    return sets


def parse_args(argv=None):
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--base", default=".", help="project folder (raw data, cache, results)")
    a.add_argument("--signature", default="standard", choices=list(SIGNATURES))
    a.add_argument("--custom-genes", default=None, help="text file, one gene per line (overrides --signature)")
    a.add_argument("--gmt", default=None, help="optional extra metabolic gene sets (GMT)")
    a.add_argument("--cancers", nargs="*", default=None, help="restrict to these TCGA codes")
    a.add_argument("--endpoints", nargs="*", default=ENDPOINTS)
    a.add_argument("--min-normals", type=int, default=5)
    a.add_argument("--min-samples", type=int, default=30)
    a.add_argument("--min-events", type=int, default=10)
    a.add_argument("--min-events-model", type=int, default=30)
    a.add_argument("--min-group", type=int, default=5)
    a.add_argument("--cutoff", default="median", choices=["median", "optimal"])
    a.add_argument("--model-endpoint", default="OS")
    a.add_argument("--lambda-rule", default="1se", choices=["1se", "min"])
    a.add_argument("--n-splits", type=int, default=100, help="repeated train/test splits for LASSO models")
    a.add_argument("--no-genomics", action="store_true")
    a.add_argument("--skip-model", action="store_true")
    a.add_argument("--seed", type=int, default=2026)
    a.add_argument("--demo", action="store_true", help="run on synthetic data (smoke test)")
    return a.parse_args(argv)


def main(argv=None):
    cfg = parse_args(argv)
    t0 = time.time()
    base = Path(cfg.base).expanduser()
    sub = "demo_" if cfg.demo else ""
    final_out = base / f"{sub}results"
    staging = base / f"{sub}results_running"          # new run is written here, swapped in at the end
    if staging.exists():
        shutil.rmtree(staging)
    P = {"raw": base / f"{sub}data_raw", "cache": base / f"{sub}data_cache", "base_out": staging}
    P["tab"], P["fig"] = P["base_out"] / "tables", P["base_out"] / "figures"
    for p in P.values():
        p.mkdir(parents=True, exist_ok=True)
    setup_logging(base / "logs" / f"{sub}pipeline_{time.strftime('%Y%m%d_%H%M%S')}.log")
    LOG.info("=" * 70)
    LOG.info("DRG pan-cancer pipeline")
    LOG.info(f"base: {base}")
    sns.set_theme(style="white", font_scale=0.8)
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "DejaVu Sans"})

    # ---- 00 data
    LOG.info("[00] Data")
    paths = {}
    if cfg.demo:
        if not (P["raw"] / DATA_FILES["expression"][1]).exists():
            make_demo_data(P["raw"], cfg.seed)
        paths = {k: P["raw"] / v[1] for k, v in DATA_FILES.items()}
    else:
        for key, (remote, local, required) in DATA_FILES.items():
            if key in ("mutation", "cnv") and cfg.no_genomics:
                continue
            try:
                url = remote if remote.startswith("http") else HUB + remote
                paths[key] = download(url, P["raw"] / local)
            except Exception as e:  # noqa: BLE001
                if required:
                    raise
                LOG.warning(f"  optional file '{key}' unavailable: {e}")
    expr = load_expression(paths["expression"], P["cache"] / "expression.parquet")
    surv = pd.read_csv(paths["survival"], sep="\t", low_memory=False)
    if "_PATIENT" not in surv.columns:
        surv["_PATIENT"] = surv["sample"].astype(str).str[:12]
    pheno = None
    if paths.get("phenotype") and paths["phenotype"].exists():
        pheno = pd.read_csv(paths["phenotype"], sep="\t", low_memory=False)

    ann = build_annotation(expr.columns, surv, pheno)
    if cfg.cancers:
        ann = ann[ann["cancer"].isin(cfg.cancers)]
    expr = expr[ann.index]
    clin = build_clinical(ann, surv)
    if paths.get("purity") and paths["purity"].exists():
        try:
            pur = load_purity(paths["purity"])
            clin["purity"] = pur.reindex(clin.index)
            LOG.info(f"  ABSOLUTE purity matched for {clin['purity'].notna().sum()} / {len(clin)} tumours")
        except Exception as ex:  # noqa: BLE001
            LOG.warning(f"  purity file could not be parsed ({ex}); purity adjustment skipped")
    ann.to_csv(P["tab"] / "00_sample_annotation.csv")
    LOG.info(f"  tumours: {len(clin)} | normals: {(ann.group == 'Normal').sum()} | "
             f"cancers: {clin.cancer.nunique()}")
    pd.crosstab(ann["cancer"], ann["group"]).to_csv(P["tab"] / "00_samples_per_cancer.csv")

    genes_req = [l.strip() for l in Path(cfg.custom_genes).read_text().splitlines() if l.strip()] \
        if cfg.custom_genes else SIGNATURES[cfg.signature]
    genes = [g for g in genes_req if g in expr.index]
    missing = [g for g in genes_req if g not in expr.index]
    if missing:
        LOG.warning(f"  genes not in matrix: {missing}")
    pd.DataFrame({"gene": genes_req, "category": [GENE_CATEGORY.get(g, "custom") for g in genes_req],
                  "found": [g in expr.index for g in genes_req]}).to_csv(P["tab"] / "00_DRG_gene_list.csv", index=False)

    # ---- 02 score
    LOG.info("[02] DRG signature scoring (ssGSEA + z-score)")
    use = ann.index[ann.group.isin(["Tumor", "Normal"])]
    arms = {"DRG": genes,
            "core": [g for g in genes if g in DRG_CORE],
            "cyto": [g for g in genes if g in DRG_CYTOSKELETON]}
    arms = {k: v for k, v in arms.items() if len(v) >= 3}
    sc_mat = ssgsea(expr[use], arms)
    score_all = sc_mat.loc["DRG"]
    scores_all = pd.DataFrame({"score": score_all})
    for k, col in (("core", "score_core"), ("cyto", "score_cyto")):
        if k in sc_mat.index:
            scores_all[col] = sc_mat.loc[k]
    z = expr.loc[genes, use]
    zmean = ((z.T - z.T.mean()) / z.T.std()).mean(axis=1)
    scores_all.assign(zscore_mean=zmean).to_csv(P["tab"] / "02_DRG_scores.csv")
    for col in scores_all.columns:
        clin[col] = scores_all[col].reindex(clin.index)
    clin["score_z"] = zmean.reindex(clin.index)
    rho = stats.spearmanr(clin["score"], clin["score_z"]).correlation
    LOG.info(f"  ssGSEA vs z-score Spearman rho = {rho:.3f}")
    clin.to_csv(P["tab"] / "02_tumour_clinical_with_scores.csv")

    extra = read_gmt(cfg.gmt) if cfg.gmt else {}
    res = {}
    steps = [
        ("tn", lambda: mod01_tumor_normal(expr, ann, genes, score_all, P, cfg)),
        ("gcox", lambda: mod03_gene_cox(expr, clin, genes, P, cfg)),
        ("sc", lambda: mod04_score_survival(clin, P, cfg)),
        ("clin", lambda: mod05_clinical(expr, clin, genes, P, cfg)),
        ("met", lambda: mod06_metabolism(expr, clin, ann, genes, P, cfg, extra)),
        ("mv", lambda: mod07_multivariable(clin, P, cfg)),
        ("model", (lambda: (pd.DataFrame(), pd.DataFrame())) if cfg.skip_model
         else (lambda: mod08_lasso_model(expr, clin, genes, P, cfg))),
        ("gen", (lambda: {}) if cfg.no_genomics
         else (lambda: mod09_genomics(paths, expr, clin, genes, P, cfg))),
        ("coex", lambda: mod10_coexpression(expr, clin, genes, P)),
        ("sub", (lambda: mod11_subscores(clin, ann, scores_all, P, cfg))
         if {"score_core", "score_cyto"} <= set(clin.columns) else (lambda: pd.DataFrame())),
    ]
    for name, fn in steps:
        try:
            res[name] = fn()
        except Exception as e:  # noqa: BLE001 — one failing module must not kill the run
            LOG.exception(f"  module '{name}' failed: {e}")
            res[name] = None
    sc, meta = res["sc"] if res.get("sc") is not None else (None, None)
    met = res["met"][0] if res.get("met") is not None else None
    perf = res["model"][0] if res.get("model") is not None else None
    write_report(P, cfg, genes, missing, ann, clin, res.get("tn"), sc, meta, res.get("mv"), met, perf, t0,
                 subres=res.get("sub"))
    write_methods(P, cfg, genes, clin)
    (P["base_out"] / "run_config.json").write_text(json.dumps(
        {**vars(cfg), "python": platform.python_version(), "pandas": pd.__version__,
         "numpy": np.__version__, "n_genes": len(genes)}, indent=2, default=str))
    swap_results(staging, final_out)
    LOG.info(f"DONE in {(time.time() - t0) / 60:.1f} min → {final_out}")


def swap_results(staging: Path, final_out: Path):
    """Replace previous pipeline outputs with the new run; keep user-added folders (e.g. ualcan_manual)."""
    final_out.mkdir(parents=True, exist_ok=True)
    managed = ["tables", "figures", "REPORT.md", "METHODS_draft.md", "run_config.json"]
    new_files = {p.relative_to(staging) for p in staging.rglob("*") if p.is_file()}
    old_files = {p.relative_to(final_out) for m in managed for p in
                 ([final_out / m] if (final_out / m).is_file() else (final_out / m).rglob("*"))
                 if p.is_file()}
    obsolete = sorted(old_files - new_files)
    for m in managed:
        tgt = final_out / m
        if tgt.is_dir():
            shutil.rmtree(tgt)
        elif tgt.exists():
            tgt.unlink()
        if (staging / m).exists():
            shutil.move(str(staging / m), str(tgt))
    shutil.rmtree(staging, ignore_errors=True)
    LOG.info(f"  previous outputs replaced; {len(obsolete)} obsolete file(s) removed")
    for f in obsolete:
        LOG.info(f"    removed: {f}")


if __name__ == "__main__":
    main()
