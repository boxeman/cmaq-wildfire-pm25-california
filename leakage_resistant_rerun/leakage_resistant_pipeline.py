"""Leakage-resistant CMAQ error modelling workflow for Bonilla et al.

Primary estimand: daily CMAQ total PM2.5 error (modelled - observed).
Primary QC: paired fire-season rows with valid observed PM2.5 in [0, 557]
and CMAQ total PM2.5 <= 557 ug/m3. Values above 557 are quarantined, never
capped, and retained in the raw-data sensitivity analysis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SEED = 20260905
N_FOLDS = 5
TEST_GROUP_FRACTION = 0.20
CMAQ_QUARANTINE_THRESHOLD = 557.0
RF_PARAMS = {
    "n_estimators": 160,
    "max_depth": 24,
    "min_samples_leaf": 10,
    "max_features": 0.7,
    "n_jobs": -1,
    "random_state": SEED,
}
PRIMARY_FEATURES = [
    "cmaq_fire_pm25",
    "latitude",
    "longitude",
    "year",
    "month",
    "dayofyear",
    "distance_to_fire_km",
    "distance_to_fire_missing",
    "n_obs_records",
]
PROHIBITED_PRIMARY_FEATURES = [
    "observed_pm25",
    "cmaq_total_pm25",
    "cmaq_nofire_pm25",
    "pm25_regime",
    "fire_name",
    "fire_window_tag",
]
REGIME_ORDER = ["Background (<12)", "Moderate (12-35)", "High (>35)"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def prepare_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df.insert(0, "source_row_id", np.arange(len(df), dtype=np.int64))
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month
    df["dayofyear"] = df["date"].dt.dayofyear
    required = [
        "date", "site_id", "observed_pm25", "cmaq_total_pm25",
        "cmaq_fire_pm25", "cmaq_nofire_pm25", "latitude", "longitude",
        "year", "month", "dayofyear", "n_obs_records",
    ]
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    df = df.loc[df["month"].between(6, 10)].dropna(subset=required).copy()
    df = df.loc[df["observed_pm25"].between(0, 557)].copy()
    df["distance_to_fire_missing"] = df["distance_to_fire_km"].isna().astype("int8")
    df["error"] = df["cmaq_total_pm25"] - df["observed_pm25"]
    df["group_id"] = df["site_id"].astype(str) + "_" + df["year"].astype(int).astype(str)
    df["pm25_regime"] = pd.cut(
        df["observed_pm25"], [-np.inf, 12, 35, np.inf], right=False,
        labels=REGIME_ORDER,
    ).astype(str)
    df["qc_status"] = np.where(
        df["cmaq_total_pm25"] > CMAQ_QUARANTINE_THRESHOLD,
        "quarantined_extreme_cmaq", "primary_included",
    )
    return df


def assign_splits(df: pd.DataFrame) -> pd.DataFrame:
    groups = df["group_id"].to_numpy()
    gss = GroupShuffleSplit(n_splits=1, test_size=TEST_GROUP_FRACTION, random_state=SEED)
    train_idx, test_idx = next(gss.split(df, groups=groups))
    split = pd.DataFrame(index=df.index)
    split["split_role"] = "train"
    split.loc[df.index[test_idx], "split_role"] = "test"
    split["cv_fold"] = -1
    train_df = df.iloc[train_idx]
    gkf = GroupKFold(n_splits=N_FOLDS)
    for fold, (_, val_pos) in enumerate(gkf.split(train_df, groups=train_df["group_id"])):
        split.loc[train_df.index[val_pos], "cv_fold"] = fold
    if (split.loc[split["split_role"] == "train", "cv_fold"] < 0).any():
        raise RuntimeError("At least one training row lacks a CV fold.")
    train_groups = set(df.loc[split["split_role"] == "train", "group_id"])
    test_groups = set(df.loc[split["split_role"] == "test", "group_id"])
    if train_groups & test_groups:
        raise RuntimeError("Group leakage detected between train and test.")
    return split


def rf_pipeline(features: list[str] | None = None) -> Pipeline:
    features = features or PRIMARY_FEATURES
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
        ("model", RandomForestRegressor(**RF_PARAMS)),
    ])


def linear_pipeline(features: list[str]) -> Pipeline:
    return Pipeline([
        ("prep", ColumnTransformer([
            ("numeric", Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
            ]), features),
        ], remainder="drop")),
        ("model", LinearRegression()),
    ])


def metrics(y_true, y_pred) -> dict[str, float | int]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return {
        "n": int(len(y_true)),
        "r2_coefficient_of_determination": float(r2_score(y_true, y_pred)),
        "pearson_r": float(np.corrcoef(y_true, y_pred)[0, 1]) if len(y_true) > 1 else np.nan,
        "rmse_ug_m3": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "mae_ug_m3": float(mean_absolute_error(y_true, y_pred)),
        "mean_residual_pred_minus_actual_ug_m3": float(np.mean(y_pred - y_true)),
    }


def fit_cv(primary_train: pd.DataFrame, outdir: Path):
    rows, predictions = [], []
    for fold in range(N_FOLDS):
        fit = primary_train["cv_fold"] != fold
        val = primary_train["cv_fold"] == fold
        model = rf_pipeline()
        model.fit(primary_train.loc[fit, PRIMARY_FEATURES], primary_train.loc[fit, "error"])
        pred = model.predict(primary_train.loc[val, PRIMARY_FEATURES])
        row = {"analysis": "primary_grouped_cv", "model": "random_forest", "fold": fold}
        row.update(metrics(primary_train.loc[val, "error"], pred))
        row["n_training_groups"] = primary_train.loc[fit, "group_id"].nunique()
        row["n_validation_groups"] = primary_train.loc[val, "group_id"].nunique()
        rows.append(row)
        predictions.append(pd.DataFrame({
            "source_row_id": primary_train.loc[val, "source_row_id"].to_numpy(),
            "cv_fold": fold,
            "actual_error_ug_m3": primary_train.loc[val, "error"].to_numpy(),
            "predicted_error_ug_m3": pred,
        }))
    cv = pd.DataFrame(rows)
    cv.to_csv(outdir / "grouped_cv_metrics.csv", index=False)
    summary_metrics = [
        "r2_coefficient_of_determination", "pearson_r", "rmse_ug_m3",
        "mae_ug_m3", "mean_residual_pred_minus_actual_ug_m3",
    ]
    pd.DataFrame([
        {"statistic": stat, **{metric: float(getattr(cv[metric], stat)()) for metric in summary_metrics}}
        for stat in ["mean", "std"]
    ]).to_csv(outdir / "grouped_cv_summary.csv", index=False)
    pd.concat(predictions, ignore_index=True).sort_values("source_row_id").to_csv(
        outdir / "grouped_cv_oof_predictions.csv", index=False
    )
    return cv


def evaluate_model(name, model, train, test, features, analysis="primary"):
    model.fit(train[features], train["error"])
    pred = model.predict(test[features])
    row = {"analysis": analysis, "model": name, "feature_set": " | ".join(features)}
    row.update(metrics(test["error"], pred))
    return model, pred, row


def compute_ale(model: Pipeline, X: pd.DataFrame, feature: str, bins=20, max_rows=20000):
    sample = X.sample(min(max_rows, len(X)), random_state=SEED).copy()
    vals = sample[feature].dropna()
    edges = np.unique(np.quantile(vals, np.linspace(0, 1, bins + 1)))
    if len(edges) < 3:
        raise ValueError(f"Not enough unique values for ALE: {feature}")
    bin_id = np.clip(np.digitize(sample[feature].fillna(vals.median()), edges[1:-1]), 0, len(edges) - 2)
    effects, counts, mids = [], [], []
    for i in range(len(edges) - 1):
        mask = bin_id == i
        n = int(mask.sum())
        counts.append(n)
        mids.append((edges[i] + edges[i + 1]) / 2)
        if n == 0:
            effects.append(0.0)
            continue
        lo, hi = sample.loc[mask].copy(), sample.loc[mask].copy()
        lo[feature], hi[feature] = edges[i], edges[i + 1]
        effects.append(float(np.mean(model.predict(hi) - model.predict(lo))))
    accumulated = np.cumsum(effects)
    weights = np.asarray(counts) / np.sum(counts)
    centered = accumulated - np.sum(accumulated * weights)
    return pd.DataFrame({
        "feature": feature, "bin": np.arange(len(mids)), "bin_lower": edges[:-1],
        "bin_upper": edges[1:], "bin_midpoint": mids, "n": counts,
        "ale_effect_ug_m3": centered,
    })


def label(feature):
    return {
        "cmaq_fire_pm25": "CMAQ fire-only PM$_{2.5}$ (µg m$^{-3}$)",
        "cmaq_nofire_pm25": "CMAQ no-fire PM$_{2.5}$ (µg m$^{-3}$)",
        "latitude": "Latitude", "longitude": "Longitude", "year": "Year",
        "month": "Month", "dayofyear": "Day of year",
        "distance_to_fire_km": "Distance to fire (km)",
        "distance_to_fire_missing": "Distance-to-fire missing indicator",
        "n_obs_records": "Daily observation count",
    }.get(feature, feature)


def save_figures(primary, test, pred, importance, ale_frames, outdir):
    plt.rcParams.update({"font.family": "Arial", "font.size": 10})
    imp = importance.sort_values("importance_mean", ascending=True)
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    ax.barh([label(x) for x in imp["feature"]], imp["importance_mean"],
            xerr=imp["importance_sd"], color="#2F6B8A", alpha=.9)
    ax.axvline(0, color="black", lw=.8)
    ax.set_xlabel("Decrease in test-set R² after permutation")
    ax.set_title("Figure S10. Leakage-resistant permutation importance")
    fig.tight_layout(); fig.savefig(outdir / "Figure_S10_permutation_importance.png", dpi=300); plt.close(fig)

    for number, ale in zip([11, 12], ale_frames):
        feature = ale["feature"].iloc[0]
        fig, ax = plt.subplots(figsize=(7.5, 5.2))
        ax.plot(ale["bin_midpoint"], ale["ale_effect_ug_m3"], marker="o", color="#B24A3B")
        ax.axhline(0, color="black", lw=.8)
        if feature == "cmaq_fire_pm25":
            ax.set_xscale("symlog", linthresh=1)
        ax.set_xlabel(label(feature)); ax.set_ylabel("Accumulated local effect on predicted error (µg m$^{-3}$)")
        ax.set_title(f"Figure S{number}. Leakage-resistant ALE: {label(feature)}")
        fig.tight_layout(); fig.savefig(outdir / f"Figure_S{number}_ALE_{feature}.png", dpi=300); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    hb = ax.hexbin(test["error"], pred, gridsize=55, bins="log", mincnt=1, cmap="viridis")
    lim = np.nanpercentile(np.r_[test["error"].to_numpy(), pred], [0.5, 99.5])
    ax.plot(lim, lim, "--", color="white", lw=1.3); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Actual CMAQ error (µg m$^{-3}$)"); ax.set_ylabel("Predicted CMAQ error (µg m$^{-3}$)")
    ax.set_title("Figure 6a. Untouched site-year holdout")
    fig.colorbar(hb, ax=ax, label="log10(count)"); fig.tight_layout()
    fig.savefig(outdir / "Figure_6a_predicted_vs_actual_error.png", dpi=300); plt.close(fig)

    plot = pd.DataFrame({"actual": test["error"].to_numpy(), "predicted": pred,
                         "regime": test["pm25_regime"].to_numpy()})
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharex=False, sharey=False)
    for ax, regime in zip(axes, REGIME_ORDER):
        sub = plot[plot["regime"] == regime]
        hb = ax.hexbin(sub["actual"], sub["predicted"], gridsize=35, bins="log", mincnt=1, cmap="magma")
        lo, hi = np.nanpercentile(np.r_[sub["actual"], sub["predicted"]], [0.5, 99.5])
        ax.plot([lo, hi], [lo, hi], "--", color="white", lw=1); ax.set_title(f"{regime}\nn={len(sub):,}")
        ax.set_xlabel("Actual error"); ax.set_ylabel("Predicted error")
    fig.suptitle("Figure 6b. Leakage-resistant test predictions by observed PM$_{2.5}$ regime")
    fig.tight_layout(); fig.savefig(outdir / "Figure_6b_predictions_by_regime.png", dpi=300); plt.close(fig)

    vals = [primary.loc[primary["pm25_regime"] == r, "error"].to_numpy() for r in REGIME_ORDER]
    means = [np.mean(x) for x in vals]
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.boxplot(vals, tick_labels=REGIME_ORDER, showfliers=False, widths=.55)
    ax.scatter(range(1, 4), means, marker="D", color="#B24A3B", label="Mean")
    ax.axhline(0, color="black", lw=.8); ax.set_ylabel("CMAQ total error (µg m$^{-3}$)")
    ax.set_title("Figure S13. Error distribution by observed concentration regime")
    ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(outdir / "Figure_S13_error_by_observed_regime.png", dpi=300); plt.close(fig)


def run(data_path: Path, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    full = prepare_data(data_path)
    split = assign_splits(full)
    full = full.join(split)
    primary = full[full["qc_status"] == "primary_included"].copy()
    train, test = primary[primary["split_role"] == "train"], primary[primary["split_role"] == "test"]

    assignments = full[["source_row_id", "date", "site_id", "year", "group_id", "split_role", "cv_fold", "qc_status"]].copy()
    assignments["date"] = assignments["date"].dt.strftime("%Y-%m-%d")
    assignments.to_csv(outdir / "fold_assignments.csv", index=False)
    fit_cv(train, outdir)

    comparison = []
    final_rf, pred, row = evaluate_model("leakage_resistant_random_forest", rf_pipeline(), train, test, PRIMARY_FEATURES)
    comparison.append(row)
    mean_pred = np.repeat(train["error"].mean(), len(test))
    row = {"analysis": "primary", "model": "training_mean_baseline", "feature_set": "none"}; row.update(metrics(test["error"], mean_pred)); comparison.append(row)
    _, linear_pred, row = evaluate_model("linear_regression", linear_pipeline(PRIMARY_FEATURES), train, test, PRIMARY_FEATURES); comparison.append(row)
    _, fire_pred, row = evaluate_model("fire_only_random_forest", rf_pipeline(["cmaq_fire_pm25"]), train, test, ["cmaq_fire_pm25"]); comparison.append(row)
    nofire_features = ["cmaq_nofire_pm25" if x == "cmaq_fire_pm25" else x for x in PRIMARY_FEATURES]
    _, nofire_pred, row = evaluate_model("no_fire_component_sensitivity_rf", rf_pipeline(nofire_features), train, test, nofire_features, "component_sensitivity"); comparison.append(row)

    raw_train, raw_test = full[full["split_role"] == "train"], full[full["split_role"] == "test"]
    _, raw_pred, row = evaluate_model("random_forest_including_quarantined_extremes", rf_pipeline(), raw_train, raw_test, PRIMARY_FEATURES, "raw_extreme_sensitivity"); comparison.append(row)
    pd.DataFrame(comparison).to_csv(outdir / "model_comparison_test_metrics.csv", index=False)

    pred_df = test[["source_row_id", "date", "site_id", "group_id", "pm25_regime", "observed_pm25", "cmaq_total_pm25", "error"]].copy()
    pred_df["predicted_error_ug_m3"] = pred
    pred_df["residual_pred_minus_actual_ug_m3"] = pred - pred_df["error"]
    pred_df["date"] = pred_df["date"].dt.strftime("%Y-%m-%d")
    pred_df.to_csv(outdir / "untouched_test_predictions.csv", index=False)

    regime_rows = []
    for regime in ["All"] + REGIME_ORDER:
        mask = np.ones(len(test), dtype=bool) if regime == "All" else test["pm25_regime"].to_numpy() == regime
        row = {"analysis": "primary_untouched_test", "regime": regime}; row.update(metrics(test.loc[mask, "error"], pred[mask])); regime_rows.append(row)
    pd.DataFrame(regime_rows).to_csv(outdir / "primary_test_metrics_by_regime.csv", index=False)

    descriptive_rows = []
    for analysis_name, frame in [("primary_qc", primary), ("raw_extreme_sensitivity", full)]:
        for regime in ["All"] + REGIME_ORDER:
            sub = frame if regime == "All" else frame[frame["pm25_regime"] == regime]
            descriptive_rows.append({
                "analysis": analysis_name, "regime": regime, "n": int(len(sub)),
                "observed_mean_ug_m3": float(sub["observed_pm25"].mean()),
                "cmaq_total_mean_ug_m3": float(sub["cmaq_total_pm25"].mean()),
                "mean_error_cmaq_minus_observed_ug_m3": float(sub["error"].mean()),
                "median_error_ug_m3": float(sub["error"].median()),
                "rmse_ug_m3": float(np.sqrt(np.mean(np.square(sub["error"])))),
                "mae_ug_m3": float(sub["error"].abs().mean()),
                "error_p05_ug_m3": float(sub["error"].quantile(.05)),
                "error_p95_ug_m3": float(sub["error"].quantile(.95)),
            })
    pd.DataFrame(descriptive_rows).to_csv(outdir / "Table_S3_regime_error_summary.csv", index=False)

    imp_sample = test.sample(min(20000, len(test)), random_state=SEED)
    perm = permutation_importance(final_rf, imp_sample[PRIMARY_FEATURES], imp_sample["error"], scoring="r2", n_repeats=10, random_state=SEED, n_jobs=-1)
    importance = pd.DataFrame({"feature": PRIMARY_FEATURES, "importance_mean": perm.importances_mean, "importance_sd": perm.importances_std}).sort_values("importance_mean", ascending=False)
    importance.to_csv(outdir / "permutation_importance_test_set.csv", index=False)
    contextual = importance.loc[~importance["feature"].isin(["cmaq_fire_pm25", "distance_to_fire_missing"]), "feature"].tolist()
    second_feature = contextual[0]
    ale_frames = [compute_ale(final_rf, test[PRIMARY_FEATURES], "cmaq_fire_pm25"), compute_ale(final_rf, test[PRIMARY_FEATURES], second_feature)]
    pd.concat(ale_frames, ignore_index=True).to_csv(outdir / "ale_values.csv", index=False)
    save_figures(primary, test, pred, importance, ale_frames, outdir)

    manifest = {
        "workflow": "Bonilla et al leakage-resistant CMAQ error rerun",
        "random_seed": SEED, "grouping_unit": "site_id-year", "test_group_fraction": TEST_GROUP_FRACTION,
        "cv_folds": N_FOLDS, "target": "cmaq_total_pm25 - observed_pm25",
        "primary_features": PRIMARY_FEATURES, "prohibited_primary_features": PROHIBITED_PRIMARY_FEATURES,
        "primary_qc_rule": "observed_pm25 in [0,557]; cmaq_total_pm25 <= 557; no capping/winsorization",
        "sensitivity_rule": "retain all otherwise-valid rows, including cmaq_total_pm25 > 557",
        "n_full_hard_valid": int(len(full)), "n_primary": int(len(primary)),
        "n_quarantined": int((full["qc_status"] != "primary_included").sum()),
        "n_groups": int(full["group_id"].nunique()), "n_test_groups": int(test["group_id"].nunique()),
        "rf_parameters": RF_PARAMS, "second_ale_feature": second_feature,
        "data_file_name": data_path.name, "data_sha256": sha256(data_path),
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__, "matplotlib": matplotlib.__version__},
        "metric_note": "R2 is sklearn's coefficient of determination; Pearson r is reported separately and is not squared or relabelled.",
    }
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"status": "complete", "outdir": str(outdir), "primary_rows": len(primary), "test_rows": len(test), "quarantined": manifest["n_quarantined"], "second_ale_feature": second_feature}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    run(args.data_csv.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
