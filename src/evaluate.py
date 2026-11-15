"""Single test-run evaluation on the held-out end region and start-region rehearsal for BRCA Spellcheck."""

import argparse
from datetime import datetime
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import sklearn.metrics as _sklearn_metrics
from sklearn.metrics import roc_auc_score

from src import baselines, config, cutoffs, features, load

_ORIG_ROC_AUC_SCORE = roc_auc_score

_REPEAT_REFUSAL_MESSAGE: str = (
    "Refusing to run: outputs/results/run_log.txt shows the test has already run. "
    "A rerun is allowed only to fix a code bug: "
    'python -m src.evaluate --bugfix "<what the bug was>"'
)
_LOCKED_INPUTS_REFUSAL_MESSAGE: str = (
    "Refusing to run: frozen.json or esm_scores.csv is not the file the first run used."
)

_HEADLINE_COLUMNS: list[str] = [
    "arm",
    "arm_name",
    "cutoff",
    "n_harmful",
    "n_harmless",
    "harmful_caught",
    "harmful_missed",
    "false_alarms",
    "harmless_cleared",
    "catch_rate_percent",
    "false_alarm_rate_percent",
    "accuracy_percent",
]

_AUROC_COLUMNS: list[str] = [
    "region_tested",
    "arm",
    "n_harmful",
    "n_harmless",
    "auroc",
]

_BAND_COUNT_COLUMNS: list[str] = [
    "findlay_class",
    "looks_harmful",
    "cant_tell",
    "looks_harmless",
    "total",
]

_BY_SWAP_COLUMNS: list[str] = [
    "arm",
    "swap_kind",
    "n_harmful",
    "n_harmless",
    "harmful_caught",
    "harmful_missed",
    "false_alarms",
    "catch_rate_percent",
    "false_alarm_rate_percent",
]

_INT_COLUMNS: list[str] = [
    "region",
    "arm",
    "n_int",
    "flagged",
    "not_flagged",
    "flagged_percent",
]

_PREDICTIONS_END_COLUMNS: list[str] = [
    "mutant",
    "position",
    "wild_type",
    "new_letter",
    "function_score",
    "findlay_class",
    "is_harmful",
    "swap_kind",
    "blosum_score",
    "esm_score",
    "entropy",
    "harm_arm1",
    "harm_arm2",
    "harm_arm3",
    "flag_arm1",
    "flag_arm2",
    "flag_arm3",
    "band_arm2",
]


def _call_roc_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Call roc_auc_score while respecting monkeypatches on either evaluate or sklearn.metrics."""
    if roc_auc_score is not _ORIG_ROC_AUC_SCORE:
        return float(roc_auc_score(y_true, y_score))
    return float(_sklearn_metrics.roc_auc_score(y_true, y_score))


def _binary_subset(pred: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Filter a prediction frame to LOF/FUNC binary rows and return (subset, int64 y_true)."""
    if "findlay_class" in pred.columns:
        mask = pred["findlay_class"].isin([config.HARMFUL_CLASS, config.HARMLESS_CLASS])
        sub = pred.loc[mask].copy()
        y_true = (sub["findlay_class"] == config.HARMFUL_CLASS).to_numpy(dtype=np.int64)
        return sub, y_true
    if "is_harmful" in pred.columns:
        mask = pred["is_harmful"].notna() & pred["is_harmful"].isin([0, 1, False, True])
        sub = pred.loc[mask].copy()
        y_true = sub["is_harmful"].astype("int64").to_numpy(dtype=np.int64)
        return sub, y_true
    raise ValueError("Prediction table must contain 'findlay_class' or 'is_harmful'.")


def _add_is_harmful_column(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of df with nullable Int64 is_harmful (1 LOF, 0 FUNC, pd.NA INT)."""
    out = df.copy()
    values: list[int | object] = []
    for cls in out["findlay_class"]:
        if cls == config.HARMFUL_CLASS:
            values.append(1)
        elif cls == config.HARMLESS_CLASS:
            values.append(0)
        else:
            values.append(pd.NA)
    out["is_harmful"] = pd.Series(values, index=out.index, dtype="Int64")
    return out


def apply_frozen_arms(table: pd.DataFrame, frozen: dict) -> pd.DataFrame:
    """Apply frozen cutoffs, arm 3 weights, and arm 2 band thresholds to a feature table."""
    out = table.copy()
    needed_cols = set(config.FEATURE_ORDER) | {"harm_arm1", "harm_arm2"}
    if not needed_cols.issubset(out.columns):
        out = features.build_feature_table(out)

    arm1_cutoff = float(frozen["arms"]["arm1"]["cutoff"])
    arm2_cutoff = float(frozen["arms"]["arm2"]["cutoff"])
    arm3_cfg = frozen["arms"]["arm3"]
    arm3_cutoff = float(arm3_cfg["cutoff"])

    x_mat = features.feature_matrix(out)
    harm_arm3 = features.arm3_probability(
        x_mat,
        arm3_cfg["scaler_mean"],
        arm3_cfg["scaler_scale"],
        arm3_cfg["coefficients"],
        float(arm3_cfg["intercept"]),
    )
    out["harm_arm3"] = pd.Series(harm_arm3, index=out.index, dtype="float64")

    flag_1 = cutoffs.flag(out["harm_arm1"].to_numpy(dtype=float), arm1_cutoff)
    flag_2 = cutoffs.flag(out["harm_arm2"].to_numpy(dtype=float), arm2_cutoff)
    flag_3 = cutoffs.flag(out["harm_arm3"].to_numpy(dtype=float), arm3_cutoff)

    out["flag_arm1"] = pd.Series(flag_1.astype("int64"), index=out.index, dtype="int64")
    out["flag_arm2"] = pd.Series(flag_2.astype("int64"), index=out.index, dtype="int64")
    out["flag_arm3"] = pd.Series(flag_3.astype("int64"), index=out.index, dtype="int64")

    band_cfg = frozen.get("band", {})
    cutoff_a = float(band_cfg.get("cutoff_a", arm2_cutoff))
    raw_cutoff_b = band_cfg.get("cutoff_b", None)
    cutoff_b = None if raw_cutoff_b is None else float(raw_cutoff_b)
    if band_cfg.get("has_band") is False:
        cutoff_b = None

    verdicts = [
        cutoffs.band_verdict(float(score), cutoff_a, cutoff_b)
        for score in out["harm_arm2"]
    ]
    out["band_arm2"] = pd.Series(verdicts, index=out.index, dtype=str)
    return out


def headline_table(pred: pd.DataFrame, frozen: dict) -> pd.DataFrame:
    """Build the four-row headline comparison table over LOF and FUNC rows."""
    binary, y_true = _binary_subset(pred)
    rows: list[dict[str, object]] = []

    for arm in ("arm0", "arm1", "arm2", "arm3"):
        arm_info = frozen.get("arms", {}).get(arm, {})
        arm_name = str(arm_info.get("name", config.ARM_NAMES[arm]))
        if arm == "arm0":
            cutoff_val: float | None = None
            scores = np.zeros(len(y_true), dtype=float)
            counts = cutoffs.count_outcomes(scores, y_true, None)
        else:
            cutoff_val = float(arm_info["cutoff"])
            harm_col = f"harm_{arm}"
            flag_col = f"flag_{arm}"
            if harm_col in binary.columns:
                scores = binary[harm_col].to_numpy(dtype=float)
                counts = cutoffs.count_outcomes(scores, y_true, cutoff_val)
            elif flag_col in binary.columns:
                flagged = binary[flag_col].to_numpy(dtype=bool)
                harmful_mask = y_true == 1
                n_harmful = int(harmful_mask.sum())
                n_harmless = int((~harmful_mask).sum())
                harmful_caught = int((flagged & harmful_mask).sum())
                false_alarms = int((flagged & (~harmful_mask)).sum())
                counts = {
                    "n_harmful": n_harmful,
                    "n_harmless": n_harmless,
                    "harmful_caught": harmful_caught,
                    "harmful_missed": n_harmful - harmful_caught,
                    "false_alarms": false_alarms,
                    "harmless_cleared": n_harmless - false_alarms,
                }
            else:
                raise ValueError(f"Missing column '{harm_col}' or '{flag_col}' in pred.")

        arm_rates = cutoffs.rates(counts)
        rows.append(
            {
                "arm": arm,
                "arm_name": arm_name,
                "cutoff": np.nan if cutoff_val is None else float(cutoff_val),
                "n_harmful": int(counts["n_harmful"]),
                "n_harmless": int(counts["n_harmless"]),
                "harmful_caught": int(counts["harmful_caught"]),
                "harmful_missed": int(counts["harmful_missed"]),
                "false_alarms": int(counts["false_alarms"]),
                "harmless_cleared": int(counts["harmless_cleared"]),
                "catch_rate_percent": float(arm_rates["catch_rate_percent"]),
                "false_alarm_rate_percent": float(arm_rates["false_alarm_rate_percent"]),
                "accuracy_percent": float(arm_rates["accuracy_percent"]),
            }
        )

    return pd.DataFrame(rows, columns=_HEADLINE_COLUMNS)


def auroc_table(pred: pd.DataFrame) -> pd.DataFrame:
    """Compute the four-row AUROC table for arm0 through arm3 on LOF and FUNC rows."""
    binary, y_true = _binary_subset(pred)
    n_harmful = int((y_true == 1).sum())
    n_harmless = int((y_true == 0).sum())

    if "region" in binary.columns and len(binary) > 0 and binary["region"].nunique() == 1:
        region_tested = str(binary["region"].iloc[0])
    else:
        region_tested = "end"

    rows: list[dict[str, object]] = [
        {
            "region_tested": region_tested,
            "arm": "arm0",
            "n_harmful": n_harmful,
            "n_harmless": n_harmless,
            "auroc": 0.5,
        }
    ]

    for arm in ("arm1", "arm2", "arm3"):
        scores = binary[f"harm_{arm}"].to_numpy(dtype=float)
        score_val = round(_call_roc_auc(y_true, scores), 4)
        rows.append(
            {
                "region_tested": region_tested,
                "arm": arm,
                "n_harmful": n_harmful,
                "n_harmless": n_harmless,
                "auroc": float(score_val),
            }
        )

    return pd.DataFrame(rows, columns=_AUROC_COLUMNS)


def bootstrap_auroc_differences(
    pred: pd.DataFrame,
    n_resamples: int = config.N_BOOTSTRAP,
    seed: int = config.SEED,
) -> dict:
    """Cluster-bootstrap position-level AUROC differences for ESM vs BLOSUM and arm 3 vs ESM."""
    binary, y_true = _binary_subset(pred)
    positions = np.array(sorted(binary["position"].astype(int).unique()), dtype=np.int64)
    n_positions = int(len(positions))

    pos_slices: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    pos_col = binary["position"].to_numpy(dtype=np.int64)
    h1_all = binary["harm_arm1"].to_numpy(dtype=float)
    h2_all = binary["harm_arm2"].to_numpy(dtype=float)
    h3_all = binary["harm_arm3"].to_numpy(dtype=float)

    for pos in positions:
        mask = pos_col == pos
        pos_slices[int(pos)] = (
            y_true[mask],
            h1_all[mask],
            h2_all[mask],
            h3_all[mask],
        )

    real_auroc_1 = _call_roc_auc(y_true, h1_all)
    real_auroc_2 = _call_roc_auc(y_true, h2_all)
    real_auroc_3 = _call_roc_auc(y_true, h3_all)
    point_esm_minus_blosum = float(real_auroc_2 - real_auroc_1)
    point_arm3_minus_esm = float(real_auroc_3 - real_auroc_2)

    rng = np.random.default_rng(seed)
    diffs_esm_blosum: list[float] = []
    diffs_arm3_esm: list[float] = []
    n_skipped = 0

    for _ in range(int(n_resamples)):
        chosen = rng.choice(positions, size=n_positions, replace=True)
        y_parts = [pos_slices[int(p)][0] for p in chosen]
        y_boot = np.concatenate(y_parts)
        if np.all(y_boot == 0) or np.all(y_boot == 1):
            n_skipped += 1
            continue

        h1_boot = np.concatenate([pos_slices[int(p)][1] for p in chosen])
        h2_boot = np.concatenate([pos_slices[int(p)][2] for p in chosen])
        h3_boot = np.concatenate([pos_slices[int(p)][3] for p in chosen])

        auroc_1 = _call_roc_auc(y_boot, h1_boot)
        auroc_2 = _call_roc_auc(y_boot, h2_boot)
        auroc_3 = _call_roc_auc(y_boot, h3_boot)

        diffs_esm_blosum.append(float(auroc_2 - auroc_1))
        diffs_arm3_esm.append(float(auroc_3 - auroc_2))

    n_used = len(diffs_esm_blosum)
    pct_low, pct_high = config.BOOTSTRAP_PERCENTILES
    if n_used > 0:
        eb_low, eb_high = np.percentile(diffs_esm_blosum, [pct_low, pct_high])
        a3_low, a3_high = np.percentile(diffs_arm3_esm, [pct_low, pct_high])
    else:
        eb_low, eb_high = float("nan"), float("nan")
        a3_low, a3_high = float("nan"), float("nan")

    return {
        "unit": "end-region positions",
        "n_positions": int(n_positions),
        "n_resamples_requested": int(n_resamples),
        "n_resamples_used": int(n_used),
        "n_resamples_skipped_one_class": int(n_skipped),
        "seed": int(seed),
        "esm_minus_blosum": {
            "point": float(point_esm_minus_blosum),
            "low": float(eb_low),
            "high": float(eb_high),
        },
        "arm3_minus_esm": {
            "point": float(point_arm3_minus_esm),
            "low": float(a3_low),
            "high": float(a3_high),
        },
    }


def band_count_table(pred: pd.DataFrame) -> pd.DataFrame:
    """Count arm 2 verdict bands across LOF, FUNC, and INT classes."""
    rows: list[dict[str, object]] = []
    for cls in (config.HARMFUL_CLASS, config.HARMLESS_CLASS, config.EXCLUDED_CLASS):
        sub = pred[pred["findlay_class"] == cls]
        verdicts = sub["band_arm2"]
        looks_harmful = int((verdicts == config.VERDICT_HARMFUL).sum())
        cant_tell = int((verdicts == config.VERDICT_CANT_TELL).sum())
        looks_harmless = int((verdicts == config.VERDICT_HARMLESS).sum())
        rows.append(
            {
                "findlay_class": cls,
                "looks_harmful": looks_harmful,
                "cant_tell": cant_tell,
                "looks_harmless": looks_harmless,
                "total": int(len(sub)),
            }
        )
    return pd.DataFrame(rows, columns=_BAND_COUNT_COLUMNS)


def by_swap_table(pred: pd.DataFrame) -> pd.DataFrame:
    """Break down arm 1, arm 2, and arm 3 outcomes on LOF and FUNC rows by swap kind."""
    binary, y_true = _binary_subset(pred)
    if "swap_kind" not in binary.columns:
        binary = binary.copy()
        binary["swap_kind"] = [
            features.swap_kind(wt, nl)
            for wt, nl in zip(binary["wild_type"], binary["new_letter"])
        ]

    swap_col = binary["swap_kind"].to_numpy(dtype=str)
    rows: list[dict[str, object]] = []

    for arm in ("arm1", "arm2", "arm3"):
        flagged_all = binary[f"flag_{arm}"].to_numpy(dtype=bool)
        for kind in config.SWAP_KINDS:
            kind_mask = swap_col == kind
            y_sub = y_true[kind_mask]
            flagged_sub = flagged_all[kind_mask]
            harmful_mask = y_sub == 1

            n_harmful = int(harmful_mask.sum())
            n_harmless = int((~harmful_mask).sum())
            harmful_caught = int((flagged_sub & harmful_mask).sum())
            harmful_missed = n_harmful - harmful_caught
            false_alarms = int((flagged_sub & (~harmful_mask)).sum())

            catch_rate = (
                round(100.0 * harmful_caught / n_harmful, 1)
                if n_harmful > 0
                else np.nan
            )
            false_alarm_rate = (
                round(100.0 * false_alarms / n_harmless, 1)
                if n_harmless > 0
                else np.nan
            )

            rows.append(
                {
                    "arm": arm,
                    "swap_kind": kind,
                    "n_harmful": n_harmful,
                    "n_harmless": n_harmless,
                    "harmful_caught": harmful_caught,
                    "harmful_missed": harmful_missed,
                    "false_alarms": false_alarms,
                    "catch_rate_percent": catch_rate,
                    "false_alarm_rate_percent": false_alarm_rate,
                }
            )

    return pd.DataFrame(rows, columns=_BY_SWAP_COLUMNS)


def int_table(pred_region: pd.DataFrame, region: str) -> pd.DataFrame:
    """Summarize how many INT rows are flagged by arms 1, 2, and 3 in one region."""
    if region not in ("start", "end"):
        raise ValueError("region must be 'start' or 'end'.")

    if "findlay_class" in pred_region.columns:
        int_rows = pred_region[pred_region["findlay_class"] == config.EXCLUDED_CLASS]
    elif "is_harmful" in pred_region.columns:
        int_rows = pred_region[pred_region["is_harmful"].isna()]
    else:
        int_rows = pred_region

    n_int = int(len(int_rows))
    rows: list[dict[str, object]] = []
    for arm in ("arm1", "arm2", "arm3"):
        flagged = int(int_rows[f"flag_{arm}"].to_numpy(dtype=bool).sum())
        not_flagged = n_int - flagged
        flagged_percent = round(100.0 * flagged / n_int, 1) if n_int > 0 else np.nan
        rows.append(
            {
                "region": region,
                "arm": arm,
                "n_int": n_int,
                "flagged": flagged,
                "not_flagged": not_flagged,
                "flagged_percent": flagged_percent,
            }
        )

    return pd.DataFrame(rows, columns=_INT_COLUMNS)


def choose_sentences(headline: pd.DataFrame, bootstrap: dict) -> list[str]:
    """Select the pre-committed result sentences from the headline table and bootstrap intervals."""
    sentences: list[str] = []
    d_interval = bootstrap["esm_minus_blosum"]
    t_interval = bootstrap["arm3_minus_esm"]

    d_low = float(d_interval["low"])
    d_high = float(d_interval["high"])
    d_point = float(d_interval["point"])

    if d_low > 0.0:
        sentences.append(
            "S1: A model that only read protein sequences ranked harmful BRCA1 "
            "mutations better than a classic substitution table."
        )
        if d_point < config.BARELY_BEATS_AUROC:
            sentences.append(
                "S4: ESM barely beat the rule; the finding is that a 1992 table is hard to beat."
            )
    elif d_low <= 0.0 <= d_high:
        sentences.append("S2: No detectable difference at this sample size.")
    else:
        sentences.append(
            "NO FIXED SENTENCE: the interval is entirely below zero, so BLOSUM62 "
            "ranked better than ESM-2. Report this plainly and ask the humans for wording."
        )

    if float(t_interval["low"]) <= 0.0:
        sentences.append(
            "S3: We could not detect a benefit from adding labels from the other end of the protein."
        )

    for arm in ("arm1", "arm2", "arm3"):
        arm_rows = headline[headline["arm"] == arm]
        if len(arm_rows) == 0:
            continue
        catch_rate = float(arm_rows.iloc[0]["catch_rate_percent"])
        if abs(catch_rate - float(config.CATCH_TARGET_PERCENT)) > config.CATCH_DRIFT_POINTS:
            sentences.append(
                f"S5: ({arm}) A cutoff tuned on one end of the protein did not carry over "
                f"to the other: it caught {catch_rate:.1f}% instead of 90%."
            )

    arm2_rows = headline[headline["arm"] == "arm2"]
    if len(arm2_rows) == 0:
        raise ValueError("headline table must contain an 'arm2' row.")
    arm2_catch = float(arm2_rows.iloc[0]["catch_rate_percent"])
    arm2_false_alarms = int(arm2_rows.iloc[0]["false_alarms"])
    sentences.append(
        f"S6: To catch {arm2_catch:.1f}% of harmful mutations, the model also "
        f"flagged {arm2_false_alarms} of 899 harmless ones."
    )

    return sentences


def verify_hash(frozen: dict, scores_path: Path | None = None) -> None:
    """Verify that the on-disk ESM score cache SHA-256 matches frozen.json."""
    target_path = config.ESM_SCORES_PATH if scores_path is None else Path(scores_path)
    expected = str(frozen["inputs"]["esm_scores_sha256"])
    found = config.sha256_of_file(target_path)
    if found != expected:
        raise RuntimeError(
            "HASH MISMATCH: outputs/esm_scores.csv has changed since frozen.json was written. "
            f"Expected {expected}, found {found}. Refusing to run."
        )


def prior_runs(run_log_path: Path | None = None) -> list[dict]:
    """Parse existing run-log lines into dictionaries; return an empty list if absent."""
    path = config.RUN_LOG_PATH if run_log_path is None else Path(run_log_path)
    if not path.is_file():
        return []

    parsed: list[dict] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            parts = [part.strip() for part in stripped.split("|")]
            entry: dict[str, str] = {}
            if parts and "=" not in parts[0]:
                entry["timestamp"] = parts[0]
            for part in parts:
                if "=" in part:
                    key, val = part.split("=", 1)
                    entry[key.strip()] = val.strip()
            parsed.append(entry)
    return parsed


def may_run(
    run_log_path: Path | None = None,
    bugfix_reason: str | None = None,
) -> tuple[bool, str]:
    """Check whether the test run is permitted given prior run-log entries and bugfix_reason."""
    path = config.RUN_LOG_PATH if run_log_path is None else Path(run_log_path)
    if bugfix_reason is not None and bugfix_reason.strip() == "":
        return False, _REPEAT_REFUSAL_MESSAGE

    runs = prior_runs(path)
    attempted_statuses = {"started", "success", "failed"}
    has_attempted = any(r.get("status") in attempted_statuses for r in runs)
    if has_attempted:
        if bugfix_reason is None or bugfix_reason.strip() == "":
            return False, _REPEAT_REFUSAL_MESSAGE
        return True, bugfix_reason.strip()

    return True, bugfix_reason.strip() if bugfix_reason else "first run"


def locked_hashes(run_log_path: Path | None = None) -> tuple[str, str] | None:
    """Return (esm_sha256, frozen_sha256) from the first started line in the log, or None."""
    path = config.RUN_LOG_PATH if run_log_path is None else Path(run_log_path)
    for entry in prior_runs(path):
        if entry.get("status") == "started":
            esm_sha = entry.get("esm_sha256")
            frozen_sha = entry.get("frozen_sha256")
            if esm_sha is not None and frozen_sha is not None:
                return (esm_sha, frozen_sha)
            return None
    return None


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    """Write a DataFrame to CSV with UTF-8 encoding and LF line endings."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        df.to_csv(handle, index=False, lineterminator="\n")


def _write_json(obj: dict, path: Path) -> None:
    """Write a dictionary to UTF-8 JSON with indent=2 and a trailing LF."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="\n", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, sort_keys=False)
        handle.write("\n")


def _format_headline_markdown(headline: pd.DataFrame) -> str:
    """Render the headline DataFrame as a GitHub-flavored Markdown table."""
    header = "| " + " | ".join(_HEADLINE_COLUMNS) + " |"
    separator = "| " + " | ".join("---" for _ in _HEADLINE_COLUMNS) + " |"
    lines = [header, separator]
    for _, row in headline.iterrows():
        cells: list[str] = []
        for col in _HEADLINE_COLUMNS:
            val = row[col]
            if col == "cutoff":
                cells.append("" if pd.isna(val) else f"{float(val):.6f}")
            elif col in (
                "catch_rate_percent",
                "false_alarm_rate_percent",
                "accuracy_percent",
            ):
                cells.append(f"{float(val):.1f}")
            elif col in (
                "n_harmful",
                "n_harmless",
                "harmful_caught",
                "harmful_missed",
                "false_alarms",
                "harmless_cleared",
            ):
                cells.append(str(int(val)))
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def _git_short_commit() -> str:
    """Return the first 7 characters of HEAD commit hash, or 'unknown' if unavailable."""
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=config.ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        out = proc.stdout.strip()
        if proc.returncode == 0 and len(out) >= 7:
            return out[:7]
    except OSError:
        pass
    return "unknown"


def _safe_sha256(path: Path) -> str:
    """Return the SHA-256 of path if it exists as a file, else 'missing'."""
    if path.is_file():
        return config.sha256_of_file(path)
    return "missing"


def _append_run_log(
    run_log_path: Path,
    run_num: int,
    status: str,
    kind: str,
    reason: str,
    esm_sha256: str,
    frozen_sha256: str,
) -> None:
    """Append a single pipe-delimited status line to run_log.txt and close the file."""
    timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
    clean_reason = (
        reason.replace("|", " ").replace("\r", " ").replace("\n", " ").strip()
    )
    commit = _git_short_commit()
    line = (
        f"{timestamp} | run={run_num} | status={status} | kind={kind} | "
        f"reason={clean_reason} | commit={commit} | "
        f"esm_sha256={esm_sha256} | frozen_sha256={frozen_sha256}\n"
    )
    run_log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(run_log_path, "a", encoding="utf-8", newline="\n") as handle:
        handle.write(line)
        handle.flush()


def _load_frozen(frozen_path: Path) -> dict:
    """Load frozen.json from disk as a dictionary."""
    if not frozen_path.is_file():
        raise FileNotFoundError(
            f"{frozen_path} not found. Run: python -m src.freeze"
        )
    with open(frozen_path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("frozen.json check failed: root must be a JSON object.")
    return data


def _build_start_predictions(frozen: dict) -> pd.DataFrame:
    """Build start-region feature and prediction table using tuning labels only."""
    labels_start = load.load_labels("start", "tuning")
    scores_start = load.lab_tested_scores("start")
    merged = pd.merge(
        labels_start,
        scores_start,
        on=["mutant", "wild_type", "position", "new_letter", "region"],
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != 575:
        raise ValueError(f"Expected 575 start-region rows after merge, got {len(merged)}.")
    feat_table = features.build_feature_table(merged)
    feat_table = _add_is_harmful_column(feat_table)
    return apply_frozen_arms(feat_table, frozen)


def _run_selfcheck_start() -> int:
    """Execute start-region rehearsal and write outputs exclusively under SELFCHECK_DIR."""
    frozen = _load_frozen(config.FROZEN_PATH)
    verify_hash(frozen, config.ESM_SCORES_PATH)

    pred_start = _build_start_predictions(frozen)
    headline = headline_table(pred_start, frozen)
    auroc = auroc_table(pred_start)
    bootstrap = bootstrap_auroc_differences(
        pred_start, n_resamples=200, seed=config.SEED
    )
    band_counts = band_count_table(pred_start)
    by_swap = by_swap_table(pred_start)
    int_rows = int_table(pred_start, "start")

    out_dir = config.SELFCHECK_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(headline, out_dir / "headline.csv")
    with open(out_dir / "headline.md", "w", newline="\n", encoding="utf-8") as handle:
        handle.write(_format_headline_markdown(headline))
    _write_csv(auroc, out_dir / "auroc.csv")
    _write_json(bootstrap, out_dir / "bootstrap.json")
    _write_csv(band_counts, out_dir / "band_counts.csv")
    _write_csv(by_swap, out_dir / "by_swap.csv")
    _write_csv(int_rows, out_dir / "int_rows.csv")

    arm1_row = headline[headline["arm"] == "arm1"].iloc[0]
    if (
        int(arm1_row["harmful_caught"]) != 120
        or int(arm1_row["n_harmful"]) != 131
        or int(arm1_row["false_alarms"]) != 217
        or int(arm1_row["n_harmless"]) != 385
    ):
        raise RuntimeError(
            "Self-check failed: arm 1 start-region counts do not match 120/131 and 217/385."
        )

    arm2_row = headline[headline["arm"] == "arm2"].iloc[0]
    if float(arm2_row["catch_rate_percent"]) < 90.0:
        raise RuntimeError(
            f"Self-check failed: arm 2 start-region catch rate {arm2_row['catch_rate_percent']}% < 90.0%."
        )

    arm3_row = headline[headline["arm"] == "arm3"].iloc[0]
    expected_arm3_caught = int(frozen["arms"]["arm3"]["tuning_caught"])
    expected_arm3_fa = int(frozen["arms"]["arm3"]["tuning_false_alarms"])
    if (
        int(arm3_row["harmful_caught"]) != expected_arm3_caught
        or int(arm3_row["false_alarms"]) != expected_arm3_fa
    ):
        raise RuntimeError(
            "Self-check failed: arm 3 start-region counts do not match frozen tuning counts."
        )

    print("SELFCHECK PASSED (start region, in-sample; these are not results).")
    return 0


def _verify_pre_run_inputs() -> tuple[dict, str, str]:
    """Run Step 2 pre-test checks before any end-region labels are loaded."""
    locked = locked_hashes(config.RUN_LOG_PATH)
    if locked is not None:
        if not config.ESM_SCORES_PATH.is_file() or not config.FROZEN_PATH.is_file():
            raise RuntimeError(_LOCKED_INPUTS_REFUSAL_MESSAGE)
        current_esm_sha = config.sha256_of_file(config.ESM_SCORES_PATH)
        current_frozen_sha = config.sha256_of_file(config.FROZEN_PATH)
        if (current_esm_sha, current_frozen_sha) != locked:
            raise RuntimeError(_LOCKED_INPUTS_REFUSAL_MESSAGE)

    frozen = _load_frozen(config.FROZEN_PATH)
    if frozen.get("schema_version") != 1:
        raise ValueError("frozen.json check failed: schema_version must be 1.")

    arm3_features = frozen.get("arms", {}).get("arm3", {}).get("features")
    if arm3_features != config.FEATURE_ORDER:
        raise ValueError(
            "frozen.json check failed: arms.arm3.features does not match FEATURE_ORDER."
        )

    verify_hash(frozen, config.ESM_SCORES_PATH)

    slim_sha = config.sha256_of_file(config.SLIM_CSV_PATH)
    expected_slim_sha = str(frozen.get("inputs", {}).get("slim_csv_sha256", ""))
    if slim_sha != expected_slim_sha:
        raise ValueError("frozen.json check failed: slim_csv_sha256 mismatch.")

    load.check_inputs()

    esm_sha = config.sha256_of_file(config.ESM_SCORES_PATH)
    frozen_sha = config.sha256_of_file(config.FROZEN_PATH)
    return frozen, esm_sha, frozen_sha


def _run_test_evaluation(bugfix_reason: str | None) -> int:
    """Execute the test run on the held-out end region with full run-log accounting."""
    allowed, msg = may_run(config.RUN_LOG_PATH, bugfix_reason)
    if not allowed:
        print(msg, file=sys.stderr)
        return 1

    kind = "bugfix" if bugfix_reason is not None else "first"
    run_reason = bugfix_reason.strip() if bugfix_reason is not None else "first run"
    existing_runs = prior_runs(config.RUN_LOG_PATH)
    run_num = sum(1 for r in existing_runs if r.get("status") == "started") + 1

    try:
        frozen, esm_sha, frozen_sha = _verify_pre_run_inputs()
    except Exception as exc:
        _append_run_log(
            config.RUN_LOG_PATH,
            run_num=run_num,
            status="refused",
            kind=kind,
            reason=str(exc) or type(exc).__name__,
            esm_sha256=_safe_sha256(config.ESM_SCORES_PATH),
            frozen_sha256=_safe_sha256(config.FROZEN_PATH),
        )
        print(str(exc), file=sys.stderr)
        return 1

    _append_run_log(
        config.RUN_LOG_PATH,
        run_num=run_num,
        status="started",
        kind=kind,
        reason=run_reason,
        esm_sha256=esm_sha,
        frozen_sha256=frozen_sha,
    )

    try:
        labels = load.load_labels("end", "evaluate")
        scores = load.lab_tested_scores("end")
        merged = pd.merge(
            labels,
            scores,
            on=["mutant", "wild_type", "position", "new_letter", "region"],
            how="inner",
            validate="one_to_one",
        )
        if len(merged) != 1262:
            raise ValueError(f"Expected 1262 end-region rows after merge, got {len(merged)}.")

        table = features.build_feature_table(merged)
        table = _add_is_harmful_column(table)
        pred = apply_frozen_arms(table, frozen)

        n_lof = int((pred["findlay_class"] == config.HARMFUL_CLASS).sum())
        n_int = int((pred["findlay_class"] == config.EXCLUDED_CLASS).sum())
        n_func = int((pred["findlay_class"] == config.HARMLESS_CLASS).sum())
        if (n_lof, n_int, n_func) != (267, 96, 899):
            raise ValueError(
                f"End-region class counts mismatch: expected (267, 96, 899), got ({n_lof}, {n_int}, {n_func})."
            )
        for harm_col in ("harm_arm1", "harm_arm2", "harm_arm3"):
            if pred[harm_col].isna().any():
                raise ValueError(f"NaN found in {harm_col}.")

        results_dir = config.RESULTS_DIR
        results_dir.mkdir(parents=True, exist_ok=True)

        headline = headline_table(pred, frozen)
        arm0_row = headline[headline["arm"] == "arm0"].iloc[0]
        if (
            int(arm0_row["harmful_caught"]) != 0
            or int(arm0_row["harmful_missed"]) != 267
            or int(arm0_row["false_alarms"]) != 0
            or int(arm0_row["harmless_cleared"]) != 899
            or float(arm0_row["accuracy_percent"]) != 77.1
        ):
            raise RuntimeError("Arm 0 headline row check failed.")

        _write_csv(headline, results_dir / "headline.csv")
        headline_md = _format_headline_markdown(headline)
        with open(results_dir / "headline.md", "w", newline="\n", encoding="utf-8") as handle:
            handle.write(headline_md)

        auroc = auroc_table(pred)
        _write_csv(auroc, results_dir / "auroc.csv")

        bootstrap = bootstrap_auroc_differences(
            pred, n_resamples=config.N_BOOTSTRAP, seed=config.SEED
        )
        _write_json(bootstrap, results_dir / "bootstrap.json")

        band_counts = band_count_table(pred)
        _write_csv(band_counts, results_dir / "band_counts.csv")

        by_swap = by_swap_table(pred)
        _write_csv(by_swap, results_dir / "by_swap.csv")

        pred_start = _build_start_predictions(frozen)
        int_start = int_table(pred_start, "start")
        int_end = int_table(pred, "end")
        int_combined = pd.concat([int_start, int_end], ignore_index=True)
        _write_csv(int_combined, results_dir / "int_rows.csv")

        predictions_end = pred[_PREDICTIONS_END_COLUMNS].copy()
        _write_csv(predictions_end, results_dir / "predictions_end.csv")

        m1775r_rows = pred.loc[pred["mutant"] == "M1775R", "esm_score"]
        a1708e_rows = pred.loc[pred["mutant"] == "A1708E", "esm_score"]
        if len(m1775r_rows) != 1 or len(a1708e_rows) != 1:
            raise ValueError("Expected M1775R and A1708E in end-region predictions.")

        all_esm = pd.concat([pred_start["esm_score"], pred["esm_score"]], ignore_index=True)
        all_func = pd.concat(
            [pred_start["function_score"], pred["function_score"]], ignore_index=True
        )
        sanity7 = {
            "M1775R_esm_score": float(m1775r_rows.iloc[0]),
            "A1708E_esm_score": float(a1708e_rows.iloc[0]),
            "spearman_all_1837": float(spearmanr(all_esm, all_func)[0]),
            "spearman_start_575": float(
                spearmanr(pred_start["esm_score"], pred_start["function_score"])[0]
            ),
            "spearman_end_1262": float(
                spearmanr(pred["esm_score"], pred["function_score"])[0]
            ),
        }
        _write_json(sanity7, results_dir / "sanity7.json")

        sentences = choose_sentences(headline, bootstrap)
        sentences_text = "\n".join(sentences) + "\n"
        with open(
            results_dir / "result_sentences.txt", "w", newline="\n", encoding="utf-8"
        ) as handle:
            handle.write(sentences_text)

        _append_run_log(
            config.RUN_LOG_PATH,
            run_num=run_num,
            status="success",
            kind=kind,
            reason=run_reason,
            esm_sha256=esm_sha,
            frozen_sha256=frozen_sha,
        )
    except BaseException as exc:
        _append_run_log(
            config.RUN_LOG_PATH,
            run_num=run_num,
            status="failed",
            kind=kind,
            reason=str(exc) or type(exc).__name__,
            esm_sha256=esm_sha,
            frozen_sha256=frozen_sha,
        )
        raise

    print(headline_md.rstrip())
    print(auroc.to_string(index=False))
    print(json.dumps(bootstrap, indent=2))
    print(band_counts.to_string(index=False))
    print(sentences_text.rstrip())
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for the start-region self-check or the single end-region test run."""
    config.set_seeds()
    parser = argparse.ArgumentParser(
        description="Run the BRCA Spellcheck start-region self-check or held-out end-region evaluation."
    )
    parser.add_argument(
        "--selfcheck-start",
        action="store_true",
        help="Rehearse all evaluation tables on the start region only.",
    )
    parser.add_argument(
        "--bugfix",
        type=str,
        default=None,
        metavar="REASON",
        help="Rerun the evaluation only to fix a code bug, logging the reason.",
    )
    args = parser.parse_args(argv)

    if args.selfcheck_start and args.bugfix is not None:
        print("Cannot combine --selfcheck-start with --bugfix.", file=sys.stderr)
        return 1

    try:
        if args.selfcheck_start:
            return _run_selfcheck_start()
        return _run_test_evaluation(args.bugfix)
    except (
        FileNotFoundError,
        ValueError,
        TypeError,
        PermissionError,
        KeyError,
        RuntimeError,
    ) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
