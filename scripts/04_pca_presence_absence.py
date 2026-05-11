#!/usr/bin/env python3
"""Run PCA on binary gene or module presence-absence features."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Integrated per-genome TSV table.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument(
        "--feature-prefix",
        default="gene__",
        help="Use columns beginning with this prefix as PCA features.",
    )
    parser.add_argument("--n-components", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    id_cols = [c for c in ["genome_id", "source_label", "ST_label", "T3SS_type", "collection_year"] if c in df.columns]
    feature_cols = [c for c in df.columns if c.startswith(args.feature_prefix)]
    if len(feature_cols) < 2:
        raise ValueError(f"Need at least two feature columns with prefix {args.feature_prefix!r}.")

    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    variable_cols = X.columns[X.nunique() > 1].tolist()
    X = X[variable_cols]
    if len(variable_cols) < 2:
        raise ValueError("Fewer than two variable features remain after filtering.")

    n_components = min(args.n_components, len(variable_cols), len(df))
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=n_components, random_state=42)
    coords = pca.fit_transform(X_scaled)

    coord_df = df[id_cols].copy()
    for i in range(n_components):
        coord_df[f"PC{i + 1}"] = coords[:, i]
    coord_df.to_csv(outdir / "pca_coordinates.tsv", sep="\t", index=False)

    var_df = pd.DataFrame(
        {
            "PC": [f"PC{i + 1}" for i in range(n_components)],
            "explained_variance_ratio": pca.explained_variance_ratio_,
        }
    )
    var_df.to_csv(outdir / "pca_variance.tsv", sep="\t", index=False)

    loadings = pd.DataFrame(pca.components_.T, index=variable_cols, columns=[f"PC{i + 1}" for i in range(n_components)])
    loadings.reset_index(names="feature").to_csv(outdir / "pca_loadings.tsv", sep="\t", index=False)
    print(f"Wrote PCA outputs to {outdir}")


if __name__ == "__main__":
    main()

