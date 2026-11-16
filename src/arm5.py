"""Exploratory arm 5 for FoldFlux: arm 2's scores with a balanced cutoff, added after the test run.

Arm 2 uses the "catch at least 90% of harmful rows" cutoff, which accepts many false alarms.
Arm 5 keeps the same ESM-2 scores and changes only the cutoff rule: the cutoff that maximises
catch rate minus false-alarm rate (Youden's J) on the start region. The ranking, and so the
AUROC, is identical to arm 2; only the yes/no calls change. It was chosen after the test run,
with knowledge of arm 2's end-region result, so it is exploratory and never replaces the
headline. The cutoff is frozen on the start region before one run on the end region.
Two steps, in order: `freeze`, `evaluate` (once).
"""

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from src import config, cutoffs, load

# Fixed in advance. The interval is the 95% position bootstrap of arm 5's accuracy minus
# arm 2's on the end region, in percentage points.
SENTENCE_BETTER: str = "On the held-out region the balanced cutoff gave more correct calls than `arm2`."
SENTENCE_NO_DIFFERENCE: str = "No detectable difference in correct calls from `arm2` at this sample size."
SENTENCE_WORSE: str = "On the held-out region the balanced cutoff gave fewer correct calls than `arm2`."

_PARAGRAPH_MARKER: str = "<!-- arm5.json -->"


def choose_balanced_cutoff(harm_scores: np.ndarray, is_harmful: np.ndarray) -> float:
    """Observed harm score maximising catch rate minus false-alarm rate; a tie goes to the higher cutoff."""
    harm = np.asarray(harm_scores, dtype=float)
    y = np.asarray(is_harmful, dtype=np.int64)
    n_harmful = int(y.sum())
    n_harmless = int(len(y) - n_harmful)
    if n_harmful == 0 or n_harmless == 0:
        raise ValueError("Need both harmful and harmless rows to choose a cutoff.")
    best_cut, best_j = None, None
    for cut in np.sort(np.unique(harm))[::-1]:
        flagged = harm >= cut
        j = (flagged & (y == 1)).sum() / n_harmful - (flagged & (y == 0)).sum() / n_harmless
        if best_j is None or j > best_j + 1e-12:
            best_cut, best_j = float(cut), float(j)
    return best_cut


def choose_sentence(difference: dict) -> str:
    """Map the bootstrap interval of arm 5 minus arm 2 accuracy to one fixed sentence."""
    if difference["low"] > 0:
        return SENTENCE_BETTER
    if difference["high"] < 0:
        return SENTENCE_WORSE
    return SENTENCE_NO_DIFFERENCE


def accuracy_bootstrap(
    positions: np.ndarray, correct_a: np.ndarray, correct_b: np.ndarray, n_resamples: int, seed: int
) -> dict:
    """Accuracy of a minus accuracy of b in points, with a 95% interval from resampling positions."""
    positions = np.asarray(positions, dtype=np.int64)
    a = np.asarray(correct_a, dtype=float)
    b = np.asarray(correct_b, dtype=float)
    unique = np.unique(positions)
    index = {p: np.flatnonzero(positions == p) for p in unique}
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(n_resamples):
        picks = rng.choice(unique, size=len(unique), replace=True)
        rows = np.concatenate([index[p] for p in picks])
        diffs.append(100.0 * (a[rows].mean() - b[rows].mean()))
    return {
        "point": float(100.0 * (a.mean() - b.mean())),
        "low": float(np.percentile(diffs, 2.5)),
        "high": float(np.percentile(diffs, 97.5)),
        "n_positions": int(len(unique)),
        "n_resamples_used": int(n_resamples),
    }


def log_statuses(run_log_path: Path) -> list[str]:
    """Return every status value in the arm 5 run log, in order; a missing file gives none."""
    if not run_log_path.is_file():
        return []
    statuses: list[str] = []
    for line in run_log_path.read_text(encoding="utf-8").splitlines():
        for part in (piece.strip() for piece in line.split("|")):
            if part.startswith("status="):
                statuses.append(part[len("status="):])
    return statuses


def _append_log(status: str, reason: str) -> None:
    """Append one line to the arm 5 run log. Only this module writes that file."""
    stamp = datetime.now().astimezone().isoformat(timespec="seconds")
    clean = reason.replace("|", " ").replace("\n", " ")
    line = (
        f"{stamp} | run=1 | status={status} | reason={clean} | "
        f"esm_scores_sha256={config.sha256_of_file(config.ESM_SCORES_PATH)} | "
        f"frozen_arm5_sha256={config.sha256_of_file(config.ARM5_FROZEN_PATH)}\n"
    )
    config.ARM5_RUN_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(config.ARM5_RUN_LOG_PATH, "a", newline="\n", encoding="utf-8") as handle:
        handle.write(line)


def _write_text(text: str, path: Path) -> None:
    """Write text through a temporary file so a reader never sees half a file."""
    tmp_path = path.parent / f"{path.name}.tmp"
    with open(tmp_path, "w", newline="\n", encoding="utf-8") as handle:
        handle.write(text)
    os.replace(tmp_path, path)


def _write_json(obj: dict, path: Path) -> None:
    """Write indented JSON with a final newline."""
    _write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n", path)


def _scored(region: str, labels: pd.DataFrame) -> pd.DataFrame:
    """Join one region's labels with the arm 2 ESM-2 scores, one row per lab-tested swap."""
    esm2 = load.lab_tested_scores(region)[["mutant", "esm_score"]]
    merged = labels.merge(esm2, on="mutant", validate="one_to_one")
    if len(merged) != len(labels):
        raise ValueError("Arm 5 check failed: a lab-tested mutation has no score.")
    merged["harm"] = 0.0 - merged["esm_score"]
    return merged


def _binary_rows(table: pd.DataFrame) -> pd.DataFrame:
    """Keep LOF and FUNC rows and add is_harmful; INT rows are left out, as everywhere."""
    keep = table["findlay_class"].isin([config.HARMFUL_CLASS, config.HARMLESS_CLASS])
    out = table.loc[keep].copy()
    out["is_harmful"] = (out["findlay_class"] == config.HARMFUL_CLASS).astype("int64")
    return out


def run_freeze() -> int:
    """Step 1: choose the balanced cutoff on the start region only, then freeze it."""
    if config.ARM5_FROZEN_PATH.exists():
        print("Refusing to overwrite outputs/frozen_arm5.json.", file=sys.stderr)
        return 1
    if log_statuses(config.ARM5_RUN_LOG_PATH):
        print("Refusing to freeze: the arm 5 run has started.", file=sys.stderr)
        return 1
    if not load.test_run_done():
        print("Refusing to freeze: arm 5 is an addition after the one test run.", file=sys.stderr)
        return 1

    start = _binary_rows(_scored("start", load.load_labels("start", "tuning")))
    if len(start) != 516:
        raise ValueError("Arm 5 check failed: expected 516 start-region LOF and FUNC rows.")
    y = start["is_harmful"].to_numpy()
    harm = start["harm"].to_numpy()
    cutoff = choose_balanced_cutoff(harm, y)
    counts = cutoffs.count_outcomes(harm, y, cutoff)
    arm2_cutoff = cutoffs.choose_catch_cutoff(harm, y, config.CATCH_TARGET_PERCENT)
    arm2_counts = cutoffs.count_outcomes(harm, y, arm2_cutoff)
    with open(config.FROZEN_PATH, "r", encoding="utf-8") as handle:
        arm2_frozen_cutoff = float(json.load(handle)["arms"]["arm2"]["cutoff"])
    if abs(arm2_cutoff - arm2_frozen_cutoff) > 1e-6:
        raise ValueError("Arm 5 check failed: recomputed arm 2 cutoff differs from frozen.json.")

    frozen = {
        "what": (
            "EXPLORATORY arm 5, added after the one test run. Same ESM-2 scores as arm 2; "
            "only the cutoff rule differs, chosen on the start region only. Nothing here "
            "changes outputs/frozen.json."
        ),
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "seed": config.SEED,
        "tuning_region": "start",
        "test_region": "end",
        "scores": "ESM-2 650M zero-shot (the arm 2 scores, outputs/esm_scores.csv)",
        "cutoff_rule": config.ARM5_CUTOFF_RULE,
        "cutoff": float(cutoff),
        "tuning": {**counts, **cutoffs.rates(counts)},
        "arm2_cutoff_same_rows": float(arm2_cutoff),
        "arm2_tuning_same_rows": {**arm2_counts, **cutoffs.rates(arm2_counts)},
        "inputs": {
            "esm_scores_sha256": config.sha256_of_file(config.ESM_SCORES_PATH),
            "frozen_sha256": config.sha256_of_file(config.FROZEN_PATH),
        },
    }
    _write_json(frozen, config.ARM5_FROZEN_PATH)
    print(json.dumps(frozen, indent=2))
    print("Start-region numbers are in-sample; they are not results.")
    return 0


def paragraph_text(result: dict) -> str:
    """The one write-up paragraph for arm 5, built only from the arm 5 result file."""
    e5 = result["end"]
    e2 = result["arm2_end_same_rows"]
    diff = result["arm5_minus_arm2_accuracy_points"]
    return (
        "**Second extra arm, added after the test run (exploratory):** `arm5` keeps `arm2`'s "
        "ESM-2 scores and changes only the cutoff: instead of catching at least 90% of harmful "
        "mutations, it takes the cutoff with the best catch rate minus false-alarm rate on "
        "`1–101`, frozen in `outputs/frozen_arm5.json` before one run on the end region "
        "(`outputs/results/arm5_run_log.txt`). It caught "
        f"{e5['harmful_caught']} of {e5['n_harmful']} harmful mutations "
        f"({e5['harmful_missed']} missed, vs {e2['harmful_missed']} for `arm2`) with "
        f"{e5['false_alarms']} of {e5['n_harmless']} false alarms (vs {e2['false_alarms']}); "
        f"{e5['accuracy_percent']:.1f}% of calls were right vs {e2['accuracy_percent']:.1f}% "
        f"(difference `{diff['point']:+.1f}` points, 95% interval `{diff['low']:+.1f}` to "
        f"`{diff['high']:+.1f}`). {result['sentence']} Most of that gain comes from the end region "
        f"having {e5['n_harmless']} harmless vs {e5['n_harmful']} harmful rows: always calling "
        f"'harmless' would already be right {100.0 * e5['n_harmless'] / (e5['n_harmless'] + e5['n_harmful']):.1f}% "
        f"of the time, and the catch rate fell from {e2['catch_rate_percent']:.1f}% to "
        f"{e5['catch_rate_percent']:.1f}%. The ranking, and so the AUROC, is the "
        "same as `arm2`'s: this is a different trade between missed harmful mutations and "
        "false alarms, not a better model. We picked it after seeing `arm2`'s end-region "
        f"results, so it is exploratory and the headline table is unchanged. {_PARAGRAPH_MARKER}\n"
    )


def run_evaluate() -> int:
    """Step 2: apply the frozen arm 5 cutoff to the end region, once."""
    if not load.test_run_done():
        print("Refusing to run: the one test run is not done.", file=sys.stderr)
        return 1
    if not config.ARM5_FROZEN_PATH.is_file():
        print("Refusing to run: outputs/frozen_arm5.json not found. Run freeze first.", file=sys.stderr)
        return 1
    if log_statuses(config.ARM5_RUN_LOG_PATH):
        print(
            "Refusing to run: outputs/results/arm5_run_log.txt shows arm 5 has already run. "
            "Ask a human.",
            file=sys.stderr,
        )
        return 1
    with open(config.ARM5_FROZEN_PATH, "r", encoding="utf-8") as handle:
        frozen = json.load(handle)
    if config.sha256_of_file(config.ESM_SCORES_PATH) != frozen["inputs"]["esm_scores_sha256"]:
        print("HASH MISMATCH: esm_scores.csv changed since frozen_arm5.json was written.", file=sys.stderr)
        return 1
    if config.sha256_of_file(config.FROZEN_PATH) != frozen["inputs"]["frozen_sha256"]:
        print("HASH MISMATCH: frozen.json changed since frozen_arm5.json was written.", file=sys.stderr)
        return 1
    with open(config.FROZEN_PATH, "r", encoding="utf-8") as handle:
        arm2_cutoff = float(json.load(handle)["arms"]["arm2"]["cutoff"])

    _append_log("started", "arm 5 run")
    try:
        end_all = _scored("end", load.load_labels("end", "post_test"))
        end = _binary_rows(end_all)
        if len(end_all) != 1262 or len(end) != 1166:
            raise ValueError("Arm 5 check failed: unexpected end-region row counts.")
        y = end["is_harmful"].to_numpy()
        harm = end["harm"].to_numpy()
        counts5 = cutoffs.count_outcomes(harm, y, frozen["cutoff"])
        counts2 = cutoffs.count_outcomes(harm, y, arm2_cutoff)
        correct5 = cutoffs.flag(harm, frozen["cutoff"]) == (y == 1)
        correct2 = cutoffs.flag(harm, arm2_cutoff) == (y == 1)
        difference = accuracy_bootstrap(
            end["position"].to_numpy(), correct5, correct2, config.N_BOOTSTRAP, config.SEED
        )
        int_rows = end_all.loc[end_all["findlay_class"] == config.EXCLUDED_CLASS]
        result = {
            "what": (
                "EXPLORATORY arm 5 on the held-out end region, run once after the headline "
                "test run. Same 1,166 rows and same ESM-2 scores as arm 2; only the cutoff "
                "differs. Not a replacement for arm 2."
            ),
            "cutoff_rule": frozen["cutoff_rule"],
            "cutoff": float(frozen["cutoff"]),
            "arm2_cutoff": arm2_cutoff,
            "end": {**counts5, **cutoffs.rates(counts5)},
            "arm2_end_same_rows": {**counts2, **cutoffs.rates(counts2)},
            "arm5_minus_arm2_accuracy_points": difference,
            "int_rows_end": {
                "n_int": int(len(int_rows)),
                "flagged": int(cutoffs.flag(int_rows["harm"].to_numpy(), frozen["cutoff"]).sum()),
            },
            "sentence": choose_sentence(difference),
        }
        _write_json(result, config.ARM5_RESULT_PATH)
        _write_text(paragraph_text(result), config.ARM5_PARAGRAPH_PATH)
    except BaseException as exc:
        _append_log("failed", str(exc) or type(exc).__name__)
        raise
    _append_log("success", "arm 5 run")
    print(json.dumps(result, indent=2))
    print(config.ARM5_PARAGRAPH_PATH.read_text(encoding="utf-8"))
    return 0


def main(argv: list[str] | None = None) -> int:
    """Command-line entry: python -m src.arm5 {freeze,evaluate}."""
    parser = argparse.ArgumentParser(description="Exploratory arm 5 (balanced cutoff on arm 2 scores).")
    parser.add_argument("step", choices=["freeze", "evaluate"])
    args = parser.parse_args(argv)
    if args.step == "freeze":
        return run_freeze()
    return run_evaluate()


if __name__ == "__main__":
    sys.exit(main())
