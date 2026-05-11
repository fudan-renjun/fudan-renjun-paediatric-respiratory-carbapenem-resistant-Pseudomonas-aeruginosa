#!/usr/bin/env python3
"""Compare local and public genomes within matched STs at module level."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_MODULES = [
    "T3SS_effectors",
    "T3SS_apparatus",
    "Alginate",
    "Biofilm_Fap",
    "Flagella",
    "Pili_T4",
    "Pyoverdine",
    "Pyochelin",
    "Protease",
    "T6SS",
    "QS_Las",
    "QS_Rhl",
    "Toxins",
    "Phenazine",
    "Carbapenem",
    "Beta_lactam",
    "Aminoglycoside",
    "Fluoroquinolone",
    "Efflux_MexAB",
    "Efflux_MexCD",
    "Efflux_MexEF",
    "Efflux_MexXY",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Integrated per-genome TSV table.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument(
        "--modules",
        nargs="*",
        default=None,
        help="Module columns to compare. Defaults to common virulence/resistance modules present in the input.",
    )
    return parser.parse_args()


def prevalence(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").fillna(0)
    return float(values.mean() * 100) if len(values) else 0.0


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t")
    required = {"source_label", "ST_label", "T3SS_type", "total_vf_genes", "total_resist_genes"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    modules = args.modules or [m for m in DEFAULT_MODULES if m in df.columns]
    if not modules:
        raise ValueError("No module columns found. Provide --modules or add module columns to the input table.")

    rows = []
    for st, st_df in sorted(df.groupby("ST_label"), key=lambda x: x[0]):
        local = st_df[st_df["source_label"] == "local"]
        public = st_df[st_df["source_label"] == "public"]
        if local.empty or public.empty:
            continue

        for module in modules:
            local_prev = prevalence(local[module])
            public_prev = prevalence(public[module])
            rows.append(
                {
                    "ST_label": st,
                    "module": module,
                    "local_n": len(local),
                    "public_n": len(public),
                    "local_prevalence_pct": round(local_prev, 1),
                    "public_prevalence_pct": round(public_prev, 1),
                    "delta_local_minus_public": round(local_prev - public_prev, 1),
                }
            )

    module_table = pd.DataFrame(rows)
    module_table.to_csv(outdir / "matched_st_module_comparison.tsv", sep="\t", index=False)

    summary_rows = []
    for st, st_df in sorted(df.groupby("ST_label"), key=lambda x: x[0]):
        local = st_df[st_df["source_label"] == "local"]
        public = st_df[st_df["source_label"] == "public"]
        if local.empty or public.empty:
            continue
        summary_rows.append(
            {
                "ST_label": st,
                "local_n": len(local),
                "public_n": len(public),
                "local_T3SS": ";".join(sorted(local["T3SS_type"].dropna().unique())),
                "public_T3SS": ";".join(sorted(public["T3SS_type"].dropna().unique())),
                "local_mean_vf_genes": round(float(local["total_vf_genes"].mean()), 2),
                "public_mean_vf_genes": round(float(public["total_vf_genes"].mean()), 2),
                "local_mean_resist_genes": round(float(local["total_resist_genes"].mean()), 2),
                "public_mean_resist_genes": round(float(public["total_resist_genes"].mean()), 2),
            }
        )

    pd.DataFrame(summary_rows).to_csv(outdir / "matched_st_summary.tsv", sep="\t", index=False)
    print(f"Wrote module comparison outputs to {outdir}")


if __name__ == "__main__":
    main()
