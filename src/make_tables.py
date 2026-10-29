#!/usr/bin/env python3
"""
Manuscript Tables 1 and 2 for the DRG pan-cancer study.
Reads   <results>/tables/*.csv   produced by drg_pancancer_pipeline.py
Writes  <results>/tables_main/Table1_samples_survival.csv
        <results>/tables_main/Table2_robust_lasso_models.csv

Usage:  python make_tables.py --results "/content/drive/MyDrive/DRG_PanCancer/results"
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def reverse_km_median(time, event):
    """Median follow-up by the reverse Kaplan-Meier method (censoring treated as the event)."""
    t = np.asarray(time, float)
    c = 1 - np.asarray(event, int)
    o = np.argsort(t)
    t, c = t[o], c[o]
    s = 1.0
    for u in np.unique(t[c == 1]):
        s *= 1 - ((t == u) & (c == 1)).sum() / (t >= u).sum()
        if s <= 0.5:
            return float(u)
    return np.nan


def table1(tab: Path) -> pd.DataFrame:
    clin = pd.read_csv(tab / "02_tumour_clinical_with_scores.csv", index_col=0)
    spc = pd.read_csv(tab / "00_samples_per_cancer.csv", index_col=0)
    tn_cancers = set(pd.read_csv(tab / "01_tumor_vs_normal.csv")["cancer"])
    rows = []
    for ca, g in clin.groupby("cancer"):
        os_ = g.dropna(subset=["OS", "OS.time"])
        os_ = os_[os_["OS.time"] > 0]
        pfi = g.dropna(subset=["PFI", "PFI.time"])
        pfi = pfi[pfi["PFI.time"] > 0]
        rows.append({
            "Cancer type": ca,
            "Used in tumour-normal comparison": ca in tn_cancers,
            "Tumours": len(g),
            "Normal": int(spc.loc[ca, "Normal"]) if "Normal" in spc.columns else 0,
            "OS (n)": len(os_),
            "OS events": int(os_["OS"].sum()),
            "PFI events": int(pfi["PFI"].sum()) if len(pfi) else np.nan,
            "Median follow-up (years)": round(reverse_km_median(os_["OS.time"], os_["OS"]), 1),
            "Purity (n)": int(g["purity"].notna().sum()) if "purity" in g else np.nan,
        })
    t = pd.DataFrame(rows)
    all_os = clin.dropna(subset=["OS", "OS.time"])
    all_os = all_os[all_os["OS.time"] > 0]
    total = {"Cancer type": "Total", "Used in tumour-normal comparison": "",
             **{c: t[c].sum() for c in ["Tumours", "Normal", "OS (n)", "OS events", "PFI events", "Purity (n)"]},
             "Median follow-up (years)": round(reverse_km_median(all_os["OS.time"], all_os["OS"]), 1)}
    return pd.concat([t, pd.DataFrame([total])], ignore_index=True)


def table2(tab: Path) -> pd.DataFrame:
    p = pd.read_csv(tab / "08_lasso_repeated_split_performance_OS.csv")
    p = p[p["robust"].astype(bool)]
    return pd.DataFrame({
        "Cancer type": p["cancer"],
        "Patients/events": p["n"].astype(str) + "/" + p["events"].astype(str),
        "Median test C-index (95% range)": [f"{a:.2f} ({b:.2f}-{c:.2f})" for a, b, c in
                                            zip(p.C_test_median, p.C_test_q025, p.C_test_q975)],
        "Optimism": p["optimism"].round(2),
        "AUC at 1/3/5 years": [f"{a:.2f}/{b:.2f}/{c:.2f}" for a, b, c in
                               zip(p.AUC_1y_median, p.AUC_3y_median, p.AUC_5y_median)],
        "Out-of-sample log-rank P": p["logrank_oos"].map(lambda v: f"{v:.2g}"),
        "Stable genes (>=60% of splits)": p["stable_genes"].str.replace(";", ", "),
    })


def main(argv=None):
    a = argparse.ArgumentParser()
    a.add_argument("--results", required=True)
    args = a.parse_args(argv)
    res = Path(args.results)
    out = res / "tables_main"
    out.mkdir(parents=True, exist_ok=True)
    table1(res / "tables").to_csv(out / "Table1_samples_survival.csv", index=False)
    table2(res / "tables").to_csv(out / "Table2_robust_lasso_models.csv", index=False)
    print("wrote", out)


if __name__ == "__main__":
    main()
