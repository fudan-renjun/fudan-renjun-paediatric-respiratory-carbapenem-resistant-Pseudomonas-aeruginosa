#!/usr/bin/env python3
"""Run gene-level local-versus-public Fisher exact tests with FDR correction."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Long gene-presence TSV table.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--min-local-n", type=int, default=3, help="Minimum local genomes per ST.")
    parser.add_argument("--min-public-n", type=int, default=3, help="Minimum public genomes per ST.")
    return parser.parse_args()


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    pvalues = pvalues.astype(float)
    ranked = pvalues.sort_values()
    n = len(ranked)
    adjusted = ranked * n / pd.Series(range(1, n + 1), index=ranked.index)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    return adjusted.reindex(pvalues.index)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    required = {"genome_id", "source_label", "ST_label", "gene", "present"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if "kind" not in df.columns:
        df["kind"] = "gene"

    df["present"] = pd.to_numeric(df["present"], errors="coerce").fillna(0).astype(int)
    rows = []
    for (st, kind, gene), sub in df.groupby(["ST_label", "kind", "gene"]):
        local = sub[sub["source_label"] == "local"]
        public = sub[sub["source_label"] == "public"]
        local_n = local["genome_id"].nunique()
        public_n = public["genome_id"].nunique()
        if local_n < args.min_local_n or public_n < args.min_public_n:
            continue

        local_pos = int(local.groupby("genome_id")["present"].max().sum())
        public_pos = int(public.groupby("genome_id")["present"].max().sum())
        table = [[local_pos, local_n - local_pos], [public_pos, public_n - public_pos]]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        local_prev = local_pos / local_n * 100
        public_prev = public_pos / public_n * 100
        rows.append(
            {
                "ST_label": st,
                "kind": kind,
                "gene": gene,
                "local_n": local_n,
                "public_n": public_n,
                "local_positive": local_pos,
                "public_positive": public_pos,
                "local_prevalence_pct": round(local_prev, 1),
                "public_prevalence_pct": round(public_prev, 1),
                "delta_local_minus_public": round(local_prev - public_prev, 1),
                "odds_ratio": odds_ratio,
                "p_value": p_value,
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No tests were run after filtering.")
    out["fdr_global"] = benjamini_hochberg(out["p_value"])
    out["fdr_within_st"] = out.groupby("ST_label")["p_value"].transform(benjamini_hochberg)
    out["fdr_within_st_kind"] = out.groupby(["ST_label", "kind"])["p_value"].transform(benjamini_hochberg)
    out = out.sort_values(["ST_label", "p_value", "gene"])
    out.to_csv(outdir / "gene_level_fisher_fdr.tsv", sep="\t", index=False)

    summary = (
        out.groupby("ST_label")
        .agg(
            n_tests=("gene", "size"),
            n_fdr_global_005=("fdr_global", lambda x: int((x < 0.05).sum())),
            n_fdr_within_st_005=("fdr_within_st", lambda x: int((x < 0.05).sum())),
            n_fdr_within_st_kind_005=("fdr_within_st_kind", lambda x: int((x < 0.05).sum())),
        )
        .reset_index()
    )
    summary.to_csv(outdir / "gene_level_fisher_fdr_summary.tsv", sep="\t", index=False)
    print(f"Wrote gene-level statistics to {outdir}")


if __name__ == "__main__":
    main()

