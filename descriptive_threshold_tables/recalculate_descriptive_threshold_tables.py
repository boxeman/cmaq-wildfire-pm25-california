"""Recalculate Bonilla et al. descriptive and threshold-detection tables."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/paired_cmaq_aqs_fireseason_2008_2018.csv"
OUT = Path(__file__).resolve().parent
THRESHOLDS = (12.0, 35.0)
REGIMES = ["Background (<12)", "Moderate (12-<35)", "High (>=35)"]


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(DATA, low_memory=False)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month
    required = [
        "date", "year", "month", "site_id", "observed_pm25", "cmaq_total_pm25",
        "cmaq_fire_pm25", "cmaq_nofire_pm25", "latitude", "longitude", "n_obs_records",
    ]
    hard = df.loc[df["month"].between(6, 10)].dropna(subset=required).copy()
    hard = hard.loc[hard["observed_pm25"].between(0, 557)].copy()
    hard["error"] = hard["cmaq_total_pm25"] - hard["observed_pm25"]
    hard["regime"] = pd.cut(
        hard["observed_pm25"], [-np.inf, 12, 35, np.inf], right=False, labels=REGIMES
    ).astype(str)
    primary = hard.loc[hard["cmaq_total_pm25"] <= 557].copy()
    return hard, primary


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else math.nan


def descriptive_record(analysis: str, stratum_type: str, stratum: str, frame: pd.DataFrame) -> dict:
    obs = frame["observed_pm25"]
    mod = frame["cmaq_total_pm25"]
    err = frame["error"]
    pearson = obs.corr(mod) if len(frame) > 1 else math.nan
    return {
        "analysis": analysis,
        "stratum_type": stratum_type,
        "stratum": stratum,
        "n_pairs": int(len(frame)),
        "n_sites": int(frame["site_id"].nunique()),
        "observed_mean_ug_m3": float(obs.mean()),
        "observed_median_ug_m3": float(obs.median()),
        "cmaq_total_mean_ug_m3": float(mod.mean()),
        "cmaq_total_median_ug_m3": float(mod.median()),
        "cmaq_fire_mean_ug_m3": float(frame["cmaq_fire_pm25"].mean()),
        "cmaq_nofire_mean_ug_m3": float(frame["cmaq_nofire_pm25"].mean()),
        "pearson_r": float(pearson),
        "r_squared": float(pearson ** 2),
        "mean_bias_ug_m3": float(err.mean()),
        "median_error_ug_m3": float(err.median()),
        "mae_ug_m3": float(err.abs().mean()),
        "rmse_ug_m3": float(np.sqrt(np.mean(np.square(err)))),
        "nmb_percent": 100 * safe_ratio(err.sum(), obs.sum()),
        "nme_percent": 100 * safe_ratio(err.abs().sum(), obs.sum()),
        "overestimated_percent": 100 * float((err > 0).mean()),
        "underestimated_percent": 100 * float((err < 0).mean()),
        "exact_match_percent": 100 * float((err == 0).mean()),
        "error_p05_ug_m3": float(err.quantile(0.05)),
        "error_p95_ug_m3": float(err.quantile(0.95)),
    }


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def detection_record(analysis: str, period: str, threshold: float, frame: pd.DataFrame) -> dict:
    observed_event = frame["observed_pm25"] >= threshold
    modeled_event = frame["cmaq_total_pm25"] >= threshold
    tp = int((observed_event & modeled_event).sum())
    fn = int((observed_event & ~modeled_event).sum())
    fp = int((~observed_event & modeled_event).sum())
    tn = int((~observed_event & ~modeled_event).sum())
    sensitivity = safe_ratio(tp, tp + fn)
    specificity = safe_ratio(tn, tn + fp)
    precision = safe_ratio(tp, tp + fp)
    npv = safe_ratio(tn, tn + fn)
    accuracy = safe_ratio(tp + tn, len(frame))
    sensitivity_ci = wilson(tp, tp + fn)
    specificity_ci = wilson(tn, tn + fp)
    precision_ci = wilson(tp, tp + fp)
    return {
        "analysis": analysis,
        "period": period,
        "threshold_ug_m3": threshold,
        "n_pairs": int(len(frame)),
        "observed_exceedances": int(observed_event.sum()),
        "modeled_exceedances": int(modeled_event.sum()),
        "true_positive": tp,
        "false_negative": fn,
        "false_positive": fp,
        "true_negative": tn,
        "observed_prevalence_percent": 100 * float(observed_event.mean()),
        "modeled_prevalence_percent": 100 * float(modeled_event.mean()),
        "sensitivity": sensitivity,
        "sensitivity_ci95_low": sensitivity_ci[0],
        "sensitivity_ci95_high": sensitivity_ci[1],
        "specificity": specificity,
        "specificity_ci95_low": specificity_ci[0],
        "specificity_ci95_high": specificity_ci[1],
        "precision_ppv": precision,
        "precision_ci95_low": precision_ci[0],
        "precision_ci95_high": precision_ci[1],
        "negative_predictive_value": npv,
        "accuracy": accuracy,
        "balanced_accuracy": (sensitivity + specificity) / 2,
        "f1_score": safe_ratio(2 * tp, 2 * tp + fp + fn),
        "critical_success_index": safe_ratio(tp, tp + fp + fn),
        "false_negative_rate": safe_ratio(fn, tp + fn),
        "false_positive_rate": safe_ratio(fp, fp + tn),
        "frequency_bias": safe_ratio(tp + fp, tp + fn),
    }


def write_table(name: str, frame: pd.DataFrame, payload: dict) -> None:
    frame.to_csv(OUT / f"{name}.csv", index=False)
    payload[name] = {"columns": list(frame.columns), "rows": frame.where(pd.notna(frame), None).values.tolist()}


def main() -> None:
    hard, primary = prepare()
    analyses = [("primary_qc", primary), ("raw_extreme_sensitivity", hard)]
    payload: dict = {}

    annual_rows = []
    for analysis, frame in analyses:
        annual_rows.append(descriptive_record(analysis, "year", "All years", frame))
        for year, group in frame.groupby("year", sort=True):
            annual_rows.append(descriptive_record(analysis, "year", str(int(year)), group))
    write_table("annual_descriptive_statistics", pd.DataFrame(annual_rows), payload)

    regime_rows = []
    for analysis, frame in analyses:
        regime_rows.append(descriptive_record(analysis, "observed_pm25_regime", "All regimes", frame))
        for regime in REGIMES:
            regime_rows.append(descriptive_record(analysis, "observed_pm25_regime", regime, frame[frame["regime"] == regime]))
    write_table("regime_descriptive_statistics", pd.DataFrame(regime_rows), payload)

    case_rows = []
    for analysis, frame in analyses:
        cases = frame.loc[frame["fire_window_tag"].fillna(False).astype(bool)]
        for fire_name in ["Klamath", "Rim"]:
            case_rows.append(descriptive_record(analysis, "case_study_fire_window", fire_name, cases[cases["fire_name"] == fire_name]))
    write_table("case_study_descriptive_statistics", pd.DataFrame(case_rows), payload)

    overall_detection = []
    annual_detection = []
    for analysis, frame in analyses:
        for threshold in THRESHOLDS:
            overall_detection.append(detection_record(analysis, "All years", threshold, frame))
        for year, group in frame.groupby("year", sort=True):
            for threshold in THRESHOLDS:
                annual_detection.append(detection_record(analysis, str(int(year)), threshold, group))
    write_table("threshold_detection_overall", pd.DataFrame(overall_detection), payload)
    write_table("threshold_detection_by_year", pd.DataFrame(annual_detection), payload)

    qc = pd.DataFrame([
        {"item": "Source paired rows passing hard validity", "value": len(hard), "unit": "rows", "definition": "June-October; complete core variables; observed PM2.5 in [0,557]"},
        {"item": "Primary-QC rows", "value": len(primary), "unit": "rows", "definition": "Hard-valid rows with CMAQ total PM2.5 <=557 ug/m3"},
        {"item": "Quarantined extreme CMAQ rows", "value": len(hard) - len(primary), "unit": "rows", "definition": "CMAQ total PM2.5 >557 ug/m3; retained only in raw sensitivity"},
        {"item": "Source SHA-256", "value": file_hash(DATA), "unit": "hash", "definition": DATA.name},
        {"item": "Threshold 12", "value": 12, "unit": "ug/m3", "definition": "Observed/modelled daily paired classification; event is value >= threshold"},
        {"item": "Threshold 35", "value": 35, "unit": "ug/m3", "definition": "Observed/modelled daily paired classification; event is value >= threshold"},
    ])
    write_table("qc_and_threshold_definitions", qc, payload)

    definitions = pd.DataFrame([
        ["Error", "CMAQ total PM2.5 minus observed PM2.5"],
        ["Pearson r", "Linear correlation between paired observed and CMAQ total PM2.5"],
        ["R-squared", "Square of Pearson r for descriptive tables; not the ML coefficient of determination"],
        ["NMB", "100 * sum(CMAQ-observed) / sum(observed)"],
        ["NME", "100 * sum(abs(CMAQ-observed)) / sum(observed)"],
        ["Sensitivity", "TP / (TP + FN): fraction of observed exceedances detected by CMAQ"],
        ["Specificity", "TN / (TN + FP): fraction of observed non-exceedances correctly classified"],
        ["Precision (PPV)", "TP / (TP + FP): fraction of modeled exceedances that were observed exceedances"],
        ["Balanced accuracy", "Mean of sensitivity and specificity"],
        ["F1 score", "Harmonic mean of sensitivity and precision"],
        ["Critical success index", "TP / (TP + FP + FN)"],
        ["Frequency bias", "Modeled exceedance count / observed exceedance count; 1 is unbiased frequency"],
        ["Wilson 95% CI", "Binomial Wilson score interval for sensitivity, specificity, and precision"],
    ], columns=["metric", "definition"])
    write_table("metric_definitions", definitions, payload)

    klamath = pd.DataFrame(case_rows).query("analysis == 'primary_qc' and stratum == 'Klamath'").iloc[0]
    rim = pd.DataFrame(case_rows).query("analysis == 'primary_qc' and stratum == 'Rim'").iloc[0]
    primary_high = pd.DataFrame(regime_rows).query("analysis == 'primary_qc' and stratum == 'High (>=35)'").iloc[0]
    raw_high = pd.DataFrame(regime_rows).query("analysis == 'raw_extreme_sensitivity' and stratum == 'High (>=35)'").iloc[0]
    reconciliation = pd.DataFrame([
        ["Table S3 high regime", "n_pairs", 2160, int(primary_high.n_pairs), int(raw_high.n_pairs), "Primary excludes 9 high-regime CMAQ extremes; raw sensitivity reproduces legacy count."],
        ["Table S3 high regime", "mean_error_ug_m3", -3.34, primary_high.mean_bias_ug_m3, raw_high.mean_bias_ug_m3, "Raw sensitivity reproduces legacy rounded value; primary result changes materially."],
        ["Table S3 high regime", "mae_ug_m3", 38.22, primary_high.mae_ug_m3, raw_high.mae_ug_m3, "Raw sensitivity reproduces legacy rounded value; primary result changes materially."],
        ["Table S3 high regime", "rmse_ug_m3", 150.37, primary_high.rmse_ug_m3, raw_high.rmse_ug_m3, "Raw sensitivity reproduces legacy rounded value; primary result changes materially."],
        ["Table 2S Klamath", "r_squared", 0.66, klamath.r_squared, klamath.r_squared, "Reproduces after rounding."],
        ["Table 2S Klamath", "mae_ug_m3", 15.75, klamath.mae_ug_m3, klamath.mae_ug_m3, "Does not reproduce from the current 958 paired fire-window rows."],
        ["Table 2S Klamath", "rmse_ug_m3", 29.37, klamath.rmse_ug_m3, klamath.rmse_ug_m3, "Does not reproduce from the current 958 paired fire-window rows."],
        ["Table 2S Klamath", "mean_bias_ug_m3", 13.04, klamath.mean_bias_ug_m3, klamath.mean_bias_ug_m3, "Does not reproduce; no quarantined extremes occur in 2008."],
        ["Table 2S Rim", "r_squared", 0.42, rim.r_squared, rim.r_squared, "Reproduces after rounding."],
        ["Table 2S Rim", "mae_ug_m3", 3.88, rim.mae_ug_m3, rim.mae_ug_m3, "Reproduces after rounding."],
        ["Table 2S Rim", "rmse_ug_m3", 6.57, rim.rmse_ug_m3, rim.rmse_ug_m3, "Reproduces after rounding."],
        ["Table 2S Rim", "mean_bias_ug_m3", -0.93, rim.mean_bias_ug_m3, rim.mean_bias_ug_m3, "Reproduces after rounding."],
    ], columns=["manuscript_table", "metric", "legacy_reported", "recalculated_primary", "recalculated_raw_sensitivity", "interpretation"])
    write_table("legacy_reconciliation", reconciliation, payload)

    (OUT / "tables_payload.json").write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    (OUT / "calculation_manifest.json").write_text(json.dumps({
        "source_file": DATA.name,
        "source_sha256": file_hash(DATA),
        "thresholds_ug_m3": THRESHOLDS,
        "primary_qc": "Observed PM2.5 in [0,557]; CMAQ total PM2.5 <=557; no capping",
        "sensitivity": "All otherwise-valid paired rows, including CMAQ total PM2.5 >557",
        "n_hard_valid": len(hard),
        "n_primary": len(primary),
        "n_quarantined": len(hard) - len(primary),
    }, indent=2), encoding="utf-8")
    print(json.dumps({"status": "complete", "tables": len(payload), "primary_rows": len(primary), "quarantined": len(hard)-len(primary)}))


if __name__ == "__main__":
    main()
