#!/usr/bin/env python3
"""Exploratory local-versus-public classification using interpretable features."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, permutation_test_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


DEFAULT_FEATURES = [
    "exoS",
    "exoU",
    "total_vf_genes",
    "total_resist_genes",
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
    parser.add_argument("--features", nargs="*", default=None, help="Feature columns to use.")
    parser.add_argument("--focus-sts", nargs="*", default=None, help="Optional ST labels for within-ST models.")
    parser.add_argument("--permutations", type=int, default=30)
    return parser.parse_args()


def make_pipeline(features: list[str], include_st: bool) -> Pipeline:
    transformers = [
        (
            "num",
            Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]),
            features,
        )
    ]
    if include_st:
        transformers.append(("st", OneHotEncoder(handle_unknown="ignore"), ["ST_label"]))
    preprocessor = ColumnTransformer(transformers)
    classifier = LogisticRegression(
        penalty="l1",
        solver="liblinear",
        class_weight="balanced",
        max_iter=5000,
        random_state=42,
    )
    return Pipeline([("preprocessor", preprocessor), ("classifier", classifier)])


def evaluate(df: pd.DataFrame, features: list[str], include_st: bool, label: str, permutations: int) -> tuple[dict, pd.DataFrame]:
    X_cols = features + (["ST_label"] if include_st else [])
    X = df[X_cols].copy()
    y = (df["source_label"] == "local").astype(int).to_numpy()
    min_class = np.bincount(y).min()
    if min_class < 2:
        raise ValueError(f"Dataset {label} has fewer than two samples in one class.")

    cv = RepeatedStratifiedKFold(n_splits=min(4, min_class), n_repeats=3, random_state=42)
    pipe = make_pipeline(features, include_st=include_st)
    bal_scores = []
    auc_scores = []
    coef_tables = []

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        warnings.filterwarnings("ignore", category=UserWarning)
        for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
            pipe.fit(X.iloc[train_idx], y[train_idx])
            prob = pipe.predict_proba(X.iloc[test_idx])[:, 1]
            pred = (prob >= 0.5).astype(int)
            bal_scores.append(balanced_accuracy_score(y[test_idx], pred))
            if len(np.unique(y[test_idx])) == 2:
                auc_scores.append(roc_auc_score(y[test_idx], prob))
            names = pipe.named_steps["preprocessor"].get_feature_names_out()
            coef = pipe.named_steps["classifier"].coef_[0]
            coef_tables.append(pd.DataFrame({"dataset": label, "fold": fold, "feature": names, "coef": coef}))

        perm_score, perm_scores, perm_p = permutation_test_score(
            pipe,
            X,
            y,
            scoring="balanced_accuracy",
            cv=cv,
            n_permutations=permutations,
            random_state=42,
            n_jobs=1,
        )

    metrics = {
        "dataset": label,
        "include_st": include_st,
        "n_samples": len(df),
        "n_local": int(y.sum()),
        "n_public": int((1 - y).sum()),
        "mean_balanced_accuracy": float(np.mean(bal_scores)),
        "std_balanced_accuracy": float(np.std(bal_scores)),
        "mean_auc": float(np.mean(auc_scores)) if auc_scores else np.nan,
        "std_auc": float(np.std(auc_scores)) if auc_scores else np.nan,
        "permutation_balanced_accuracy": float(perm_score),
        "permutation_p_value": float(perm_p),
    }
    coef_df = pd.concat(coef_tables, ignore_index=True)
    coef_summary = (
        coef_df.assign(abs_coef=lambda x: x["coef"].abs())
        .groupby(["dataset", "feature"])
        .agg(mean_coef=("coef", "mean"), mean_abs_coef=("abs_coef", "mean"))
        .reset_index()
        .sort_values(["dataset", "mean_abs_coef"], ascending=[True, False])
    )
    return metrics, coef_summary


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    required = {"source_label", "ST_label"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    features = args.features or [f for f in DEFAULT_FEATURES if f in df.columns]
    if not features:
        raise ValueError("No usable feature columns found.")
    for feature in features:
        df[feature] = pd.to_numeric(df[feature], errors="coerce")

    metrics = []
    feature_tables = []
    for include_st, label in [(False, "pooled_modules_only"), (True, "pooled_modules_plus_st")]:
        res, feats = evaluate(df, features, include_st=include_st, label=label, permutations=args.permutations)
        metrics.append(res)
        feature_tables.append(feats)

    focus_sts = args.focus_sts or sorted(df["ST_label"].dropna().unique())
    for st in focus_sts:
        sub = df[df["ST_label"] == st].copy()
        if sub["source_label"].nunique() < 2:
            continue
        try:
            res, feats = evaluate(sub, features, include_st=False, label=f"{st}_within_st", permutations=args.permutations)
        except ValueError:
            continue
        metrics.append(res)
        feature_tables.append(feats)

    pd.DataFrame(metrics).to_csv(outdir / "ml_performance.tsv", sep="\t", index=False)
    pd.concat(feature_tables, ignore_index=True).to_csv(outdir / "ml_feature_importance.tsv", sep="\t", index=False)
    print(f"Wrote machine-learning outputs to {outdir}")


if __name__ == "__main__":
    main()
