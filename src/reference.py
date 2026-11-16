"""Post-test context for BRCA Spellcheck: how our numbers sit next to public ones.

This module runs only after the one test run has succeeded. It adds no arm and
tunes nothing we report as ours. It does two things. First, it reads ProteinGym's
published leaderboard row for this lab table and records where ESM-2 650M ranks.
Second, it scores AlphaMissense's published BRCA1 predictions (a stronger,
third-party tool) on the same held-out end-region rows, with a cutoff chosen on
the start region by the same 90%-catch rule. It never writes frozen.json, the
score cache, the run log or any stage 5 result file.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from src import config, cutoffs, load

_ALPHAMISSENSE_COLUMNS: list[str] = ["protein_variant", "am_pathogenicity", "am_class"]
_LEADERBOARD_COLUMNS: list[str] = ["model", "spearman", "auc"]


def leaderboard_summary(board: pd.DataFrame, our_model: str) -> dict:
    """Rank one model on a leaderboard row for both Spearman and AUC (ties share a rank)."""
    if list(board.columns) != _LEADERBOARD_COLUMNS:
        raise ValueError("Leaderboard check failed: unexpected columns.")
    if board["model"].duplicated().any() or board[["spearman", "auc"]].isna().any().any():
        raise ValueError("Leaderboard check failed: duplicate model or missing value.")
    if our_model not in set(board["model"]):
        raise ValueError(f"Leaderboard check failed: {our_model} not found.")

    summary: dict = {"n_models": int(len(board)), "our_model": our_model}
    for metric in ("spearman", "auc"):
        values = board[metric].to_numpy(dtype=float)
        ours = float(board.loc[board["model"] == our_model, metric].iloc[0])
        best_index = int(np.argmax(values))
        summary[metric] = {
            "our_model_value": ours,
            "our_model_rank": int((values > ours).sum()) + 1,
            "best_model": str(board["model"].iloc[best_index]),
            "best_value": float(values[best_index]),
            "median_value": float(np.median(values)),
        }
    return summary


def parse_alphamissense(raw: pd.DataFrame, wild_type_seq: str) -> pd.DataFrame:
    """Check the AlphaMissense table against our wild type and return mutant, harm_am."""
    if list(raw.columns) != _ALPHAMISSENSE_COLUMNS:
        raise ValueError("AlphaMissense check failed: unexpected columns.")
    if raw["protein_variant"].duplicated().any() or raw.isna().any().any():
        raise ValueError("AlphaMissense check failed: duplicate or missing value.")
    scores = raw["am_pathogenicity"].to_numpy(dtype=float)
    if not np.all((scores >= 0.0) & (scores <= 1.0)):
        raise ValueError("AlphaMissense check failed: score outside 0 to 1.")
    n_wrong = 0
    for mutant in raw["protein_variant"]:
        wild_type, position, _new_letter = load.parse_mutant(mutant)
        if position > len(wild_type_seq) or wild_type_seq[position - 1] != wild_type:
            n_wrong += 1
    if n_wrong > 0:
        raise ValueError(
            f"AlphaMissense check failed: {n_wrong} rows disagree with our wild-type sequence."
        )
    # Higher AlphaMissense pathogenicity already means more likely harmful.
    return pd.DataFrame({"mutant": raw["protein_variant"], "harm_am": scores})


def paired_position_bootstrap(
    positions: np.ndarray,
    is_harmful: np.ndarray,
    harm_a: np.ndarray,
    harm_b: np.ndarray,
    n_resamples: int,
    seed: int,
) -> dict:
    """AUROC of a minus AUROC of b, with a percentile interval from resampling positions."""
    positions = np.asarray(positions, dtype=np.int64)
    is_harmful = np.asarray(is_harmful, dtype=np.int64)
    harm_a = np.asarray(harm_a, dtype=float)
    harm_b = np.asarray(harm_b, dtype=float)
    distinct = np.array(sorted(set(positions.tolist())), dtype=np.int64)
    rows_of = {int(p): np.flatnonzero(positions == p) for p in distinct}

    point = float(roc_auc_score(is_harmful, harm_a) - roc_auc_score(is_harmful, harm_b))
    rng = np.random.default_rng(seed)
    differences: list[float] = []
    n_skipped = 0
    for _ in range(int(n_resamples)):
        chosen = rng.choice(distinct, size=len(distinct), replace=True)
        index = np.concatenate([rows_of[int(p)] for p in chosen])
        y_boot = is_harmful[index]
        if np.all(y_boot == 0) or np.all(y_boot == 1):
            n_skipped += 1
            continue
        differences.append(
            float(roc_auc_score(y_boot, harm_a[index]) - roc_auc_score(y_boot, harm_b[index]))
        )
    low, high = np.percentile(differences, list(config.BOOTSTRAP_PERCENTILES))
    return {
        "point": point,
        "low": float(low),
        "high": float(high),
        "n_positions": int(len(distinct)),
        "n_resamples_used": int(len(differences)),
        "n_resamples_skipped_one_class": int(n_skipped),
    }


def _binary_rows(table: pd.DataFrame) -> pd.DataFrame:
    """Keep LOF and FUNC rows and add is_harmful (INT rows are left out, as in the headline)."""
    keep = table["findlay_class"].isin([config.HARMFUL_CLASS, config.HARMLESS_CLASS])
    out = table.loc[keep].copy()
    out["is_harmful"] = (out["findlay_class"] == config.HARMFUL_CLASS).astype("int64")
    return out


def alphamissense_reference() -> dict:
    """Score AlphaMissense on the held-out end region next to arm 2, same rows, same cutoff rule."""
    if config.sha256_of_file(config.ALPHAMISSENSE_PATH) != config.ALPHAMISSENSE_SHA256:
        raise ValueError("AlphaMissense check failed: file hash differs from the recorded one.")
    raw = pd.read_csv(
        config.ALPHAMISSENSE_PATH,
        dtype={"protein_variant": str, "am_pathogenicity": "float64", "am_class": str},
        keep_default_na=False,
    )
    am = parse_alphamissense(raw, load.load_wild_type())

    start = load.load_labels("start", "tuning").merge(am, on="mutant", validate="one_to_one")
    end_labels = load.load_labels("end", "post_test")
    end = end_labels.merge(am, on="mutant", validate="one_to_one")
    end = end.merge(
        load.lab_tested_scores("end")[["mutant", "esm_score"]], on="mutant", validate="one_to_one"
    )
    if len(start) != 575 or len(end) != 1262:
        raise ValueError("AlphaMissense check failed: a lab-tested mutation has no score.")
    end["harm_arm2"] = 0.0 - end["esm_score"]

    start_binary = _binary_rows(start)
    end_binary = _binary_rows(end)
    y_end = end_binary["is_harmful"].to_numpy()

    # Same rule as every arm: highest cutoff catching at least 90% of start-region harmful rows.
    cutoff = cutoffs.choose_catch_cutoff(
        start_binary["harm_am"].to_numpy(),
        start_binary["is_harmful"].to_numpy(),
        config.CATCH_TARGET_PERCENT,
    )
    start_counts = cutoffs.count_outcomes(
        start_binary["harm_am"].to_numpy(), start_binary["is_harmful"].to_numpy(), cutoff
    )
    end_counts = cutoffs.count_outcomes(end_binary["harm_am"].to_numpy(), y_end, cutoff)

    auroc_am = float(roc_auc_score(y_end, end_binary["harm_am"].to_numpy()))
    auroc_arm2 = float(roc_auc_score(y_end, end_binary["harm_arm2"].to_numpy()))
    difference = paired_position_bootstrap(
        end_binary["position"].to_numpy(),
        y_end,
        end_binary["harm_am"].to_numpy(),
        end_binary["harm_arm2"].to_numpy(),
        config.N_BOOTSTRAP,
        config.SEED,
    )
    both = pd.concat([start, end], ignore_index=True)

    return {
        "what": (
            "REFERENCE ONLY, added after the one test run. AlphaMissense is a third-party "
            "tool, not one of our arms. Same held-out end-region rows as the headline."
        ),
        "source_url": config.ALPHAMISSENSE_URL,
        "source_sha256": config.ALPHAMISSENSE_SHA256,
        "retrieved": config.REFERENCE_RETRIEVED,
        "caveat": (
            "AlphaMissense also uses predicted protein structure and population-frequency "
            "data, so this is context from a different kind of tool, not a like-for-like contest."
        ),
        "cutoff_rule": "highest cutoff catching at least 90% of start-region harmful rows",
        "cutoff": float(cutoff),
        "start_tuning": {k: int(v) for k, v in start_counts.items()},
        "end": {
            **{k: int(v) for k, v in end_counts.items()},
            **cutoffs.rates(end_counts),
            "auroc": round(auroc_am, 4),
            "auroc_arm2_same_rows": round(auroc_arm2, 4),
        },
        "alphamissense_minus_esm": difference,
        "spearman_with_lab_score": {
            # Negated so that, like the ESM score, higher means the protein still works.
            "all_1837": float(spearmanr(-both["harm_am"], both["function_score"])[0]),
            "end_1262": float(spearmanr(-end["harm_am"], end["function_score"])[0]),
        },
    }


def proteingym_reference() -> dict:
    """Put our overall Spearman next to ProteinGym's published row for this lab table."""
    board = pd.read_csv(
        config.PROTEINGYM_ROW_PATH,
        dtype={"model": str, "spearman": "float64", "auc": "float64"},
        keep_default_na=False,
    )
    summary = leaderboard_summary(board, config.PROTEINGYM_OUR_MODEL)
    with open(config.RESULTS_DIR / "sanity7.json", "r", encoding="utf-8") as handle:
        ours = json.load(handle)
    summary["what"] = (
        "REFERENCE ONLY. ProteinGym's published zero-shot scores on all 1,837 rows of this "
        "table. Its AUC uses ProteinGym's own yes/no labels, so it is not comparable with "
        "our end-region AUROC; the Spearman values are comparable with spearman_all_1837."
    )
    summary["source_url"] = config.PROTEINGYM_URL
    summary["retrieved"] = config.REFERENCE_RETRIEVED
    summary["our_spearman_all_1837"] = float(ours["spearman_all_1837"])
    return summary


def _write_json(obj: dict, path: os.PathLike) -> None:
    """Write JSON through a temporary file so a reader never sees half a file."""
    tmp_path = f"{path}.tmp"
    with open(tmp_path, "w", newline="\n", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, sort_keys=False)
        handle.write("\n")
    os.replace(tmp_path, path)


def main(argv: list[str] | None = None) -> int:
    """Write the two reference files. Refuses before the one test run has succeeded."""
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    config.set_seeds()
    if not load.test_run_done():
        print("Refusing to run: the test run is not done.", file=sys.stderr)
        return 1
    try:
        proteingym = proteingym_reference()
        alphamissense = alphamissense_reference()
    except (ValueError, FileNotFoundError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    _write_json(proteingym, config.REFERENCE_PROTEINGYM_RESULT_PATH)
    _write_json(alphamissense, config.REFERENCE_ALPHAMISSENSE_RESULT_PATH)
    print(json.dumps(proteingym, indent=2))
    print(json.dumps(alphamissense, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
