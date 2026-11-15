"""Baseline arm 0 (always harmless) and arm 1 (BLOSUM62 rule) scoring and start-region CLI."""

import argparse
import sys

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices

from src import config

_BLOSUM62_CACHE: dict[tuple[str, str], int] | None = None


def load_blosum62() -> dict[tuple[str, str], int]:
    """Load, validate, and cache the 20x20 BLOSUM62 substitution dictionary."""
    global _BLOSUM62_CACHE
    if _BLOSUM62_CACHE is not None:
        return dict(_BLOSUM62_CACHE)

    matrix = substitution_matrices.load("BLOSUM62")
    table: dict[tuple[str, str], int] = {}
    for aa_a in config.AA_ORDER:
        for aa_b in config.AA_ORDER:
            if aa_a not in matrix.alphabet or aa_b not in matrix.alphabet:
                raise ValueError(
                    f"BLOSUM62 check failed: missing letter for pair ({aa_a}, {aa_b})."
                )
            raw_val = float(matrix[aa_a, aa_b])
            if not raw_val.is_integer():
                raise ValueError(
                    f"BLOSUM62 check failed: non-integer value {raw_val} for ({aa_a}, {aa_b})."
                )
            table[(aa_a, aa_b)] = int(raw_val)

    if len(table) != 400:
        raise ValueError(f"BLOSUM62 check failed: expected 400 pairs, got {len(table)}.")

    for aa_a in config.AA_ORDER:
        for aa_b in config.AA_ORDER:
            if table[(aa_a, aa_b)] != table[(aa_b, aa_a)]:
                raise ValueError(
                    f"BLOSUM62 check failed: asymmetric entry for ({aa_a}, {aa_b})."
                )

    spot_expected = {("C", "G"): -3, ("D", "E"): 2, ("M", "I"): 1, ("W", "W"): 11}
    for pair, expected in spot_expected.items():
        if table[pair] != expected:
            raise ValueError(
                f"BLOSUM62 check failed: spot value {pair} expected {expected}, got {table[pair]}."
            )

    off_diag = [val for (aa_a, aa_b), val in table.items() if aa_a != aa_b]
    if len(off_diag) != 380 or min(off_diag) != -4 or max(off_diag) != 3:
        raise ValueError(
            f"BLOSUM62 check failed: off-diagonal min/max expected (-4, 3), "
            f"got ({min(off_diag)}, {max(off_diag)})."
        )

    _BLOSUM62_CACHE = table
    return dict(_BLOSUM62_CACHE)


def blosum_score(wild_type: str, new_letter: str) -> int:
    """Return the integer BLOSUM62 score for replacing wild_type with new_letter."""
    return load_blosum62()[(wild_type, new_letter)]


def harm_arm1(wild_type: str, new_letter: str) -> float:
    """Return the negated BLOSUM62 harm score, avoiding negative zero."""
    return 0.0 - blosum_score(wild_type, new_letter)


def add_arm1(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of df with integer blosum_score and float harm_arm1 columns added."""
    missing = [col for col in ("wild_type", "new_letter") if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}.")

    result = df.copy()
    table = load_blosum62()
    scores = [table[(wt, new)] for wt, new in zip(result["wild_type"], result["new_letter"])]
    harms = [0.0 - score for score in scores]
    result["blosum_score"] = pd.Series(scores, index=result.index, dtype="int64")
    result["harm_arm1"] = pd.Series(harms, index=result.index, dtype="float64")
    return result


def arm0_flags(n_rows: int) -> np.ndarray:
    """Return an all-False boolean array representing arm 0's predictions."""
    if isinstance(n_rows, bool) or not isinstance(n_rows, (int, np.integer)) or n_rows < 0:
        raise ValueError("n_rows must be a non-negative integer.")
    return np.zeros(int(n_rows), dtype=bool)


def main(argv: list[str] | None = None) -> int:
    """Evaluate arm 1 on start-region tuning rows and report arm 0 end-region counts."""
    from src import cutoffs, load

    config.set_seeds()
    parser = argparse.ArgumentParser(
        description="Compute arm 1 baseline cutoff on start-region tuning data."
    )
    parser.parse_args(argv)

    try:
        labels = load.load_labels("start", "tuning")
        binary = labels[
            labels["findlay_class"].isin([config.HARMFUL_CLASS, config.HARMLESS_CLASS])
        ].copy()
        table = add_arm1(binary)
        is_harmful = (table["findlay_class"] == config.HARMFUL_CLASS).to_numpy()
        cutoff = (
            cutoffs.choose_catch_cutoff(
                table["harm_arm1"],
                is_harmful,
                config.CATCH_TARGET_PERCENT,
            )
            + 0.0
        )
        counts = cutoffs.count_outcomes(table["harm_arm1"], is_harmful, cutoff)
        arm1_rates = cutoffs.rates(counts)
        total_rows = counts["n_harmful"] + counts["n_harmless"]
        total_flagged = counts["harmful_caught"] + counts["false_alarms"]
        flagged_percent = round(100.0 * total_flagged / total_rows, 1)

        end_counts = load.class_counts()["end"]
        end_harmful = end_counts[config.HARMFUL_CLASS]
        end_harmless = end_counts[config.HARMLESS_CLASS]
        end_total = end_harmful + end_harmless
        arm0_end_counts = {
            "n_harmful": end_harmful,
            "n_harmless": end_harmless,
            "harmful_caught": 0,
            "harmful_missed": end_harmful,
            "false_alarms": 0,
            "harmless_cleared": end_harmless,
        }
        arm0_end_rates = cutoffs.rates(arm0_end_counts)
    except (FileNotFoundError, ValueError, TypeError, PermissionError, KeyError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    print(f"Arm 1 ({config.ARM_NAMES['arm1']}) on start region (tuning):")
    print(f"  cutoff (harm scale): {cutoff}")
    print(
        f"  caught: {counts['harmful_caught']} of {counts['n_harmful']} "
        f"({arm1_rates['catch_rate_percent']:.1f}%)"
    )
    print(
        f"  false alarms: {counts['false_alarms']} of {counts['n_harmless']} "
        f"({arm1_rates['false_alarm_rate_percent']:.1f}%)"
    )
    print(f"  accuracy: {arm1_rates['accuracy_percent']:.1f}%")
    print(f"  flagged: {total_flagged} of {total_rows} ({flagged_percent:.1f}%)")
    print(f"Arm 0 ({config.ARM_NAMES['arm0']}) on end region (from class_counts):")
    print(
        f"  caught: {arm0_end_counts['harmful_caught']} of {arm0_end_counts['n_harmful']} "
        f"({arm0_end_rates['catch_rate_percent']:.1f}%)"
    )
    print(
        f"  right on: {arm0_end_counts['harmless_cleared']} of {end_total} "
        f"({arm0_end_rates['accuracy_percent']:.1f}%)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
