"""Export self-contained data and configuration files for the BRCA Spellcheck Gradio demo."""

import argparse
import json
import os
from pathlib import Path
import sys

import pandas as pd

from src import config, load

_APP_DATA_COLUMNS: list[str] = [
    "position",
    "wild_type",
    "new_letter",
    "esm_score",
    "lab_class",
    "lab_score",
]


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


def _format_lab_score(val: object) -> str:
    """Format a lab function score rounded to 3 decimals, or empty string when missing."""
    if pd.isna(val):
        return ""
    rounded = round(float(val), 3)
    if rounded == 0.0:
        rounded = 0.0
    return f"{rounded:.3f}"


def _write_csv_atomic(df: pd.DataFrame, path: Path) -> None:
    """Atomically write a DataFrame to CSV with UTF-8 encoding and LF line endings."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        df.to_csv(handle, index=False, lineterminator="\n")
    os.replace(tmp_path, path)


def _write_json_atomic(obj: dict, path: Path) -> None:
    """Atomically write a dictionary to UTF-8 JSON with indent=2 and a trailing LF."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with open(tmp_path, "w", newline="\n", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, sort_keys=False)
        handle.write("\n")
    os.replace(tmp_path, path)


def export_app_files() -> tuple[int, int, bool]:
    """Build and write app/app_data.csv and app/app_config.json and return summary counts."""
    if not config.FROZEN_PATH.is_file():
        raise FileNotFoundError(
            "outputs/frozen.json not found. Run: python -m src.freeze"
        )
    with open(config.FROZEN_PATH, "r", encoding="utf-8") as handle:
        frozen = json.load(handle)

    verify_hash(frozen, config.ESM_SCORES_PATH)
    scores = load.load_scores()
    done = load.test_run_done()

    start_labels = load.load_labels("start", "tuning")[
        ["mutant", "findlay_class", "function_score"]
    ]
    if done:
        end_labels = load.load_labels("end", "post_test")[
            ["mutant", "findlay_class", "function_score"]
        ]
        labels = pd.concat([start_labels, end_labels], ignore_index=True)
    else:
        labels = start_labels.copy()

    merged = scores.merge(
        labels,
        on="mutant",
        how="left",
        validate="one_to_one",
        sort=False,
    )
    if not (merged["mutant"].to_numpy() == scores["mutant"].to_numpy()).all():
        raise RuntimeError("Row order changed during merge with lab labels.")

    lab_class_series = merged["findlay_class"].fillna("").astype(str)
    lab_score_series = pd.Series(
        [_format_lab_score(v) for v in merged["function_score"]],
        index=merged.index,
        dtype=str,
    )
    esm_score_series = pd.Series(
        [f"{float(v):.6f}" for v in merged["esm_score"]],
        index=merged.index,
        dtype=str,
    )

    app_df = pd.DataFrame(
        {
            "position": merged["position"].astype("int64"),
            "wild_type": merged["wild_type"].astype(str),
            "new_letter": merged["new_letter"].astype(str),
            "esm_score": esm_score_series,
            "lab_class": lab_class_series,
            "lab_score": lab_score_series,
        },
        columns=_APP_DATA_COLUMNS,
    )

    band_cfg = frozen["band"]
    raw_cutoff_b = band_cfg.get("cutoff_b")
    app_cfg = {
        "model_id": str(frozen["model_id"]),
        "cutoff_a": float(band_cfg["cutoff_a"]),
        "cutoff_b": None if raw_cutoff_b is None else float(raw_cutoff_b),
        "has_band": bool(band_cfg["has_band"]),
        "end_labels_included": bool(done),
        "protein_length": int(config.PROTEIN_LENGTH),
        "n_tested_positions": int(config.N_POSITIONS),
    }

    _write_csv_atomic(app_df, config.APP_DATA_PATH)
    _write_json_atomic(app_cfg, config.APP_CONFIG_PATH)

    n_rows = int(len(app_df))
    n_labeled = int((app_df["lab_class"] != "").sum())
    return n_rows, n_labeled, bool(done)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point to export app/app_data.csv and app/app_config.json."""
    config.set_seeds()
    parser = argparse.ArgumentParser(
        description="Build app/app_data.csv and app/app_config.json for the Gradio demo."
    )
    parser.parse_args(argv)

    try:
        n_rows, n_labeled, done = export_app_files()
    except (FileNotFoundError, ValueError, PermissionError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    print(
        f"Wrote app/app_data.csv ({n_rows} rows, {n_labeled} rows with a lab_class); "
        f"end_labels_included: {done}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
