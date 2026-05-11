#!/usr/bin/env python3
"""Summarise temporal patterns among local isolates."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Integrated per-genome TSV table.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--genes", nargs="*", default=None, help="Optional binary gene columns to summarise by year.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    required = {"source_label", "collection_year", "ST_label", "T3SS_type", "total_vf_genes"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    local = df[df["source_label"] == "local"].copy()
    local["collection_year"] = pd.to_numeric(local["collection_year"], errors="coerce")
    local = local.dropna(subset=["collection_year"])
    local["collection_year"] = local["collection_year"].astype(int)

    year_counts = local.groupby("collection_year").size().reset_index(name="n_isolates")
    year_counts.to_csv(outdir / "year_counts.tsv", sep="\t", index=False)

    t3ss = (
        local.groupby(["collection_year", "T3SS_type"])
        .size()
        .reset_index(name="n")
        .merge(year_counts, on="collection_year", how="left")
    )
    t3ss["prevalence_pct"] = (t3ss["n"] / t3ss["n_isolates"] * 100).round(1)
    t3ss.to_csv(outdir / "t3ss_by_year.tsv", sep="\t", index=False)

    st_counts = (
        local.groupby(["collection_year", "ST_label"])
        .size()
        .reset_index(name="n")
        .sort_values(["collection_year", "n"], ascending=[True, False])
    )
    st_counts.to_csv(outdir / "st_counts_by_year.tsv", sep="\t", index=False)

    vf_summary = local.groupby("collection_year")["total_vf_genes"].mean().reset_index()
    vf_summary = vf_summary.rename(columns={"total_vf_genes": "mean_total_vf_genes"})
    if len(vf_summary) >= 3:
        rho, p_value = spearmanr(vf_summary["collection_year"], vf_summary["mean_total_vf_genes"])
    else:
        rho, p_value = float("nan"), float("nan")
    vf_summary["spearman_rho_all_years"] = rho
    vf_summary["spearman_p_value_all_years"] = p_value
    vf_summary.to_csv(outdir / "virulence_burden_by_year.tsv", sep="\t", index=False)

    genes = [g for g in (args.genes or []) if g in local.columns]
    gene_rows = []
    for gene in genes:
        tmp = local.groupby("collection_year")[gene].mean().reset_index()
        for _, row in tmp.iterrows():
            gene_rows.append(
                {
                    "collection_year": int(row["collection_year"]),
                    "gene": gene,
                    "prevalence_pct": round(float(row[gene]) * 100, 1),
                }
            )
    if gene_rows:
        pd.DataFrame(gene_rows).to_csv(outdir / "selected_gene_prevalence_by_year.tsv", sep="\t", index=False)

    print(f"Wrote temporal summaries to {outdir}")


if __name__ == "__main__":
    main()

