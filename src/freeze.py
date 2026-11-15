"""Start-region tuning and pre-test parameter freezing into outputs/frozen.json."""

import argparse
from datetime import datetime
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd

from src import baselines, config, cutoffs, features, load

_TEST_STARTED_REFUSAL_MSG: str = (
    "Refusing to freeze: the test run has started. Cutoffs may not change."
)
_OVERWRITE_REFUSAL_MSG: str = (
    "Refusing to overwrite outputs/frozen.json. "
    "Pass --force only if a human asked for a re-freeze before the test run."
)


def _start_tuning_table() -> pd.DataFrame:
    """Build the 516-row start-region LOF/FUNC tuning feature table with INT rows excluded."""
    labels = load.load_labels("start", "tuning")
    expected_start = config.EXPECTED_COUNTS["start"]
    expected_total = (
        expected_start[config.HARMFUL_CLASS]
        + expected_start[config.EXCLUDED_CLASS]
        + expected_start[config.HARMLESS_CLASS]
    )
    if len(labels) != expected_total:
        raise ValueError(
            f"Expected {expected_total} start-region label rows, got {len(labels)}."
        )

    scores = load.lab_tested_scores("start")
    if len(scores) != expected_total:
        raise ValueError(
            f"Expected {expected_total} start-region score rows, got {len(scores)}."
        )

    score_cols_to_merge = [
        col
        for col in scores.columns
        if col not in ("wild_type", "position", "new_letter", "region")
    ]
    merged = pd.merge(
        labels,
        scores[score_cols_to_merge],
        on="mutant",
        validate="one_to_one",
    )
    if len(merged) != expected_total:
        raise ValueError(
            f"Expected {expected_total} merged start-region rows, got {len(merged)}."
        )

    table = features.build_feature_table(merged)
    binary = table[
        table["findlay_class"].isin([config.HARMFUL_CLASS, config.HARMLESS_CLASS])
    ].copy()

    n_lof = int((binary["findlay_class"] == config.HARMFUL_CLASS).sum())
    n_func = int((binary["findlay_class"] == config.HARMLESS_CLASS).sum())
    expected_binary = (
        expected_start[config.HARMFUL_CLASS] + expected_start[config.HARMLESS_CLASS]
    )
    if (
        len(binary) != expected_binary
        or n_lof != expected_start[config.HARMFUL_CLASS]
        or n_func != expected_start[config.HARMLESS_CLASS]
    ):
        raise ValueError(
            f"Start tuning table counts mismatch: got {len(binary)} rows "
            f"({n_lof} LOF, {n_func} FUNC)."
        )

    binary["is_harmful"] = (
        binary["findlay_class"] == config.HARMFUL_CLASS
    ).astype("int64")
    return binary.reset_index(drop=True)


def _fit_arm3(table: pd.DataFrame) -> dict[str, Any]:
    """Extract features, binary labels, and position groups from table and fit Arm 3."""
    x = features.feature_matrix(table)
    y = table["is_harmful"].to_numpy()
    groups = table["position"].to_numpy()
    return features.fit_arm3(x, y, groups)


def _git_metadata() -> tuple[str, bool]:
    """Return (HEAD commit SHA, code_dirty flag ignoring outputs/frozen.json)."""
    commit_proc = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=config.ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    code_commit = commit_proc.stdout.strip()

    status_proc = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=config.ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    code_dirty = False
    for line in status_proc.stdout.splitlines():
        if not line.strip():
            continue
        path_part = line[3:].strip().replace("\\", "/")
        if path_part != "outputs/frozen.json":
            code_dirty = True
            break
    return code_commit, code_dirty


def _print_summary(
    frozen_obj: dict[str, Any],
    n_tuning_rows: int,
) -> None:
    """Print a human-readable summary of all frozen cutoffs, tuning counts, band, and Arm 3 weights."""
    n_harmful = int(frozen_obj["tuning_counts"]["harmful"])
    n_harmless = int(frozen_obj["tuning_counts"]["harmless"])

    print("Frozen tuning summary (start region, 516 LOF + FUNC rows):")
    for arm_key in ("arm1", "arm2"):
        arm_info = frozen_obj["arms"][arm_key]
        caught = int(arm_info["tuning_caught"])
        fa = int(arm_info["tuning_false_alarms"])
        flagged = caught + fa
        catch_pct = round(100.0 * caught / n_harmful, 1)
        fa_pct = round(100.0 * fa / n_harmless, 1)
        flag_pct = round(100.0 * flagged / n_tuning_rows, 1)
        print(
            f"  {arm_key} ({arm_info['name']}): cutoff={arm_info['cutoff']}, "
            f"caught={caught}/{n_harmful} ({catch_pct:.1f}%), "
            f"false_alarms={fa}/{n_harmless} ({fa_pct:.1f}%), "
            f"flagged={flagged}/{n_tuning_rows} ({flag_pct:.1f}%)"
        )

    arm3_info = frozen_obj["arms"]["arm3"]
    oof_caught = int(arm3_info["out_of_fold_caught"])
    oof_fa = int(arm3_info["out_of_fold_false_alarms"])
    oof_flagged = oof_caught + oof_fa
    oof_catch_pct = round(100.0 * oof_caught / n_harmful, 1)
    oof_fa_pct = round(100.0 * oof_fa / n_harmless, 1)
    oof_flag_pct = round(100.0 * oof_flagged / n_tuning_rows, 1)

    tun_caught = int(arm3_info["tuning_caught"])
    tun_fa = int(arm3_info["tuning_false_alarms"])
    tun_flagged = tun_caught + tun_fa
    tun_catch_pct = round(100.0 * tun_caught / n_harmful, 1)
    tun_fa_pct = round(100.0 * tun_fa / n_harmless, 1)
    tun_flag_pct = round(100.0 * tun_flagged / n_tuning_rows, 1)

    print(
        f"  arm3 ({arm3_info['name']}): cutoff={arm3_info['cutoff']}, "
        f"out_of_fold caught={oof_caught}/{n_harmful} ({oof_catch_pct:.1f}%), "
        f"out_of_fold false_alarms={oof_fa}/{n_harmless} ({oof_fa_pct:.1f}%), "
        f"out_of_fold flagged={oof_flagged}/{n_tuning_rows} ({oof_flag_pct:.1f}%), "
        f"in_sample caught={tun_caught}/{n_harmful} ({tun_catch_pct:.1f}%), "
        f"in_sample false_alarms={tun_fa}/{n_harmless} ({tun_fa_pct:.1f}%), "
        f"in_sample flagged={tun_flagged}/{n_tuning_rows} ({tun_flag_pct:.1f}%)"
    )

    band_info = frozen_obj["band"]
    print(
        f"  band (arm2): has_band={band_info['has_band']}, "
        f"cutoff_a={band_info['cutoff_a']}, cutoff_b={band_info['cutoff_b']}"
    )
    print(f"  arm3 converged: {arm3_info['converged']}")
    print("  arm3 coefficients (standardized features):")
    for feat_name, coef_val in zip(arm3_info["features"], arm3_info["coefficients"]):
        print(f"    {feat_name}: {coef_val}")
    print(f"    intercept: {arm3_info['intercept']}")


def main(argv: list[str] | None = None) -> int:
    """Tune start-region cutoffs and Arm 3 model, then atomically write outputs/frozen.json."""
    import sklearn

    config.set_seeds()
    parser = argparse.ArgumentParser(
        description="Tune on the start region and write outputs/frozen.json."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing outputs/frozen.json if the test run has not started.",
    )
    args = parser.parse_args(argv)

    if load.test_attempted():
        print(_TEST_STARTED_REFUSAL_MSG, file=sys.stderr)
        return 1

    if config.FROZEN_PATH.exists() and not args.force:
        print(_OVERWRITE_REFUSAL_MSG, file=sys.stderr)
        return 1

    try:
        load.check_inputs()
        all_scores = load.load_scores()
        start_counts = load.class_counts()["start"]
        expected_start = config.EXPECTED_COUNTS["start"]
        if start_counts != expected_start:
            raise ValueError(
                f"Start-region class counts {start_counts} do not match {expected_start}."
            )

        table = _start_tuning_table()
        is_harmful = table["is_harmful"].to_numpy()

        # Verify Arm 1 baseline helper consistency and compute cutoff.
        _ = baselines.arm0_flags(len(table))
        cutoff_1 = (
            float(
                cutoffs.choose_catch_cutoff(
                    table["harm_arm1"], is_harmful, config.CATCH_TARGET_PERCENT
                )
            )
            + 0.0
        )
        counts_1 = cutoffs.count_outcomes(table["harm_arm1"], is_harmful, cutoff_1)
        if (
            cutoff_1 != 0.0
            or counts_1["harmful_caught"] != 120
            or counts_1["false_alarms"] != 217
        ):
            raise ValueError(
                f"Arm 1 start-region check failed: cutoff={cutoff_1}, "
                f"caught={counts_1['harmful_caught']}, false_alarms={counts_1['false_alarms']}."
            )

        cutoff_a = (
            float(
                cutoffs.choose_catch_cutoff(
                    table["harm_arm2"], is_harmful, config.CATCH_TARGET_PERCENT
                )
            )
            + 0.0
        )
        counts_2 = cutoffs.count_outcomes(table["harm_arm2"], is_harmful, cutoff_a)

        raw_cutoff_b = cutoffs.choose_false_alarm_cutoff(
            table["harm_arm2"], is_harmful, config.FALSE_ALARM_TARGET_PERCENT
        )
        cutoff_b = None if raw_cutoff_b is None else float(raw_cutoff_b) + 0.0
        band_exists = bool(cutoffs.has_band(cutoff_a, cutoff_b))

        arm3_res = _fit_arm3(table)
        cutoff_3 = float(arm3_res["cutoff"]) + 0.0
        oof_counts_3 = cutoffs.count_outcomes(
            arm3_res["out_of_fold_probabilities"], is_harmful, cutoff_3
        )
        in_sample_probs_3 = features.arm3_probability(
            features.feature_matrix(table),
            arm3_res["scaler_mean"],
            arm3_res["scaler_scale"],
            arm3_res["coefficients"],
            arm3_res["intercept"],
        )
        tuning_counts_3 = cutoffs.count_outcomes(
            in_sample_probs_3, is_harmful, cutoff_3
        )

        esm_sha256 = config.sha256_of_file(config.ESM_SCORES_PATH)
        slim_sha256 = config.sha256_of_file(config.SLIM_CSV_PATH)
        code_commit, code_dirty = _git_metadata()

        frozen_obj: dict[str, Any] = {
            "schema_version": 1,
            "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "code_commit": str(code_commit),
            "code_dirty": bool(code_dirty),
            "seed": int(config.SEED),
            "model_id": str(all_scores["model_id"].iloc[0]),
            "inputs": {
                "esm_scores_path": "outputs/esm_scores.csv",
                "esm_scores_sha256": str(esm_sha256),
                "esm_scores_rows": int(len(all_scores)),
                "slim_csv_sha256": str(slim_sha256),
            },
            "tuning_region": "start",
            "test_region": "end",
            "labels": {
                "harmful_class": str(config.HARMFUL_CLASS),
                "harmless_class": str(config.HARMLESS_CLASS),
                "excluded_class": str(config.EXCLUDED_CLASS),
                "lof_below": float(config.LOF_BELOW),
                "func_above": float(config.FUNC_ABOVE),
            },
            "tuning_counts": {
                "harmful": int(start_counts[config.HARMFUL_CLASS]),
                "harmless": int(start_counts[config.HARMLESS_CLASS]),
                "int_excluded": int(start_counts[config.EXCLUDED_CLASS]),
            },
            "catch_target_percent": int(config.CATCH_TARGET_PERCENT),
            "false_alarm_target_percent": int(config.FALSE_ALARM_TARGET_PERCENT),
            "arms": {
                "arm0": {
                    "name": str(config.ARM_NAMES["arm0"]),
                    "cutoff": None,
                },
                "arm1": {
                    "name": str(config.ARM_NAMES["arm1"]),
                    "harm_score": "0 - blosum62(wild_type, new_letter)",
                    "cutoff": float(cutoff_1),
                    "tuning_caught": int(counts_1["harmful_caught"]),
                    "tuning_false_alarms": int(counts_1["false_alarms"]),
                },
                "arm2": {
                    "name": str(config.ARM_NAMES["arm2"]),
                    "harm_score": "0 - esm_score",
                    "cutoff": float(cutoff_a),
                    "tuning_caught": int(counts_2["harmful_caught"]),
                    "tuning_false_alarms": int(counts_2["false_alarms"]),
                },
                "arm3": {
                    "name": str(config.ARM_NAMES["arm3"]),
                    "harm_score": "predicted probability of harmful",
                    "features": [str(f) for f in config.FEATURE_ORDER],
                    "scaler_mean": [float(v) for v in arm3_res["scaler_mean"]],
                    "scaler_scale": [float(v) for v in arm3_res["scaler_scale"]],
                    "coefficients": [float(v) for v in arm3_res["coefficients"]],
                    "intercept": float(arm3_res["intercept"]),
                    "cutoff": float(cutoff_3),
                    "cutoff_source": (
                        "out-of-fold probabilities, "
                        "StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42), "
                        "groups = position"
                    ),
                    "model_settings": {
                        "C": 1.0,
                        "solver": "lbfgs",
                        "max_iter": 100,
                        "random_state": int(config.SEED),
                    },
                    "converged": bool(arm3_res["converged"]),
                    "out_of_fold_caught": int(oof_counts_3["harmful_caught"]),
                    "out_of_fold_false_alarms": int(oof_counts_3["false_alarms"]),
                    "tuning_caught": int(tuning_counts_3["harmful_caught"]),
                    "tuning_false_alarms": int(tuning_counts_3["false_alarms"]),
                },
            },
            "band": {
                "arm": "arm2",
                "cutoff_a": float(cutoff_a),
                "cutoff_b": None if cutoff_b is None else float(cutoff_b),
                "has_band": bool(band_exists),
            },
            "versions": {
                "python": platform.python_version(),
                "numpy": str(np.__version__),
                "pandas": str(pd.__version__),
                "scikit_learn": str(sklearn.__version__),
                "biopython": str(importlib.metadata.version("biopython")),
            },
        }

        frozen_path = config.FROZEN_PATH
        frozen_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = frozen_path.parent / f"{frozen_path.name}.tmp"
        with open(tmp_path, "w", newline="\n", encoding="utf-8") as handle:
            json.dump(frozen_obj, handle, indent=2, sort_keys=False)
            handle.write("\n")
        os.replace(tmp_path, frozen_path)

        with open(frozen_path, "r", encoding="utf-8") as handle:
            round_trip = json.load(handle)
        if round_trip != frozen_obj:
            raise RuntimeError("outputs/frozen.json round-trip verification failed.")

        _print_summary(frozen_obj, len(table))
        return 0
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
