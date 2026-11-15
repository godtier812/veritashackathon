"""Stage 6 failure analysis: select Arm 2's five worst misses and five worst false alarms on the end region."""

import argparse
from collections.abc import Mapping
from pathlib import Path
import sys

import pandas as pd

from src import config, load

_NOT_DONE_REFUSAL_MESSAGE: str = "Refusing to run: the test run is not done."
_STOP_SPLICING_MESSAGE: str = (
    "STOP: human decision needed. Tell the humans the column names above. "
    "Do not guess which column is the RNA score."
)
_NO_FULL_TABLE_MESSAGE: str = (
    "No full study table found. Splicing is stated as a limit we could not test."
)

_KNOWN_DATA_FILES: frozenset[str] = frozenset(
    {
        "BRCA1_HUMAN_Findlay_2018.csv",
        "brca1_missense_scores_slim.csv",
        "BRCA1_HUMAN_Findlay_2018_wildtype.fasta",
        ".DS_Store",
    }
)

FAILURE_COLUMNS: list[str] = [
    "kind",
    "rank",
    "mutant",
    "position",
    "wild_type",
    "new_letter",
    "function_score",
    "findlay_class",
    "esm_score",
    "harm_arm2",
    "band_arm2",
    "blosum_score",
    "entropy",
    "swap_kind",
    "flag_arm1",
    "flag_arm3",
    "guess_why",
]

CURATED_GUESSES: dict[str, str] = {
    "A1823S": (
        "Guess: Position 1823 has high entropy (2.572233) and the A-to-S swap "
        "(swap_kind=other) has a positive BLOSUM62 score (1) and a positive ESM-2 "
        "score (0.746405, harm_arm2=-0.746405, band_arm2=Looks harmless), so all "
        "three arms left it unflagged (flag_arm1=0, flag_arm3=0) despite the lab "
        "LOF score (-3.642831), possibly because the DNA change disrupted RNA "
        "splicing that a protein-sequence model cannot see."
    ),
    "A1823T": (
        "Guess: At high-entropy position 1823 (entropy=2.572233), the A-to-T swap "
        "(swap_kind=other, blosum_score=0) received a positive ESM-2 score "
        "(0.259899, harm_arm2=-0.259899, band_arm2=Looks harmless) and was missed "
        "by Arm 2 and Arm 3 (flag_arm3=0, flag_arm1=1 only at the 0.0 tie cutoff), "
        "possibly because its lab LOF score (-3.186519) reflects an untested DNA "
        "splicing effect rather than a protein-sequence penalty."
    ),
    "I1824M": (
        "Guess: The I-to-M swap (swap_kind=other) at high-entropy position 1824 "
        "(entropy=2.362860) has a positive BLOSUM62 score (1) and a mild ESM-2 "
        "score (-0.979514, harm_arm2=0.979514, band_arm2=Looks harmless), causing "
        "all three arms to leave it unflagged (flag_arm1=0, flag_arm3=0) despite "
        "the lab LOF score (-2.078020), possibly due to a DNA splicing effect or "
        "local constraint not reflected in sequence statistics."
    ),
    "R1670K": (
        "Guess: Although position 1670 has low entropy (0.733420), the R-to-K "
        "swap preserves charge (swap_kind=other) with a high BLOSUM62 score (2) "
        "and a mild ESM-2 score (-1.005743, harm_arm2=1.005743, band_arm2=Looks "
        "harmless), leading all three arms to leave it unflagged (flag_arm1=0, "
        "flag_arm3=0) despite the lab LOF score (-1.552698), possibly because the "
        "assay required the exact R letter or intact DNA splicing."
    ),
    "A1693D": (
        "Guess: Position 1693 has moderate-to-high entropy (1.779203) so ESM-2 "
        "gave this A-to-D charge_change swap a mild score (-1.014801, "
        "harm_arm2=1.014801, band_arm2=Looks harmless) below Arm 2's cutoff, "
        "whereas Arm 1 and Arm 3 both flagged it (flag_arm1=1, flag_arm3=1) from "
        "its negative BLOSUM62 score (-2) and chemistry features, matching the "
        "lab LOF score (-1.539991)."
    ),
    "L1844P": (
        "Guess: The L-to-P swap (swap_kind=to_proline) carries a low "
        "BLOSUM62 score (-3) and an extreme ESM-2 penalty "
        "(esm_score=-12.934240, harm_arm2=12.934240, band_arm2=Looks harmful) "
        "that caused all three arms to flag it (flag_arm1=1, flag_arm3=1), but "
        "high-entropy position 1844 (entropy=2.041959) tolerated the swap in the "
        "lab assay (function_score=-0.232682, FUNC)."
    ),
    "E1661V": (
        "Guess: Low entropy at position 1661 (entropy=0.392916) and the E-to-V "
        "charge_change swap (blosum_score=-2) produced a large ESM-2 penalty "
        "(esm_score=-12.026711, harm_arm2=12.026711, band_arm2=Looks harmful) and "
        "flags from all three arms (flag_arm1=1, flag_arm3=1), yet the lab "
        "function score (-0.584913) remained in the FUNC range."
    ),
    "C1787F": (
        "Guess: Near-zero entropy at position 1787 (entropy=0.044551) and the "
        "C-to-F loss_of_cysteine swap (blosum_score=-2) led to a steep ESM-2 "
        "penalty (esm_score=-11.638639, harm_arm2=11.638639, band_arm2=Looks "
        "harmful) and flags across all three arms (flag_arm1=1, flag_arm3=1), "
        "even though this swap scored as FUNC (-0.299372) in the lab assay."
    ),
    "Q1779P": (
        "Guess: Low entropy at position 1779 (entropy=0.328521) combined with a "
        "Q-to-P to_proline swap (blosum_score=-1) produced a heavy ESM-2 penalty "
        "(esm_score=-11.623756, harm_arm2=11.623756, band_arm2=Looks harmful) and "
        "flags from all three arms (flag_arm1=1, flag_arm3=1), but the variant "
        "scored as FUNC (-0.597504) in the lab assay."
    ),
    "V1653G": (
        "Guess: Very low entropy at position 1653 (entropy=0.153824) and a "
        "negative BLOSUM62 score (-3) for the V-to-G swap (swap_kind=other) drove "
        "a steep ESM-2 penalty (esm_score=-11.412864, harm_arm2=11.412864, "
        "band_arm2=Looks harmful) and flags from all three arms (flag_arm1=1, "
        "flag_arm3=1), although the lab assay scored it as FUNC (-0.263024)."
    ),
}


def _find_optional_full_table(explicit_path: str | Path | None = None) -> Path | None:
    """Return the explicit --full-table path or any unknown regular file in config.DATA_DIR."""
    if explicit_path is not None:
        return Path(explicit_path)
    if not config.DATA_DIR.is_dir():
        return None
    for entry in sorted(config.DATA_DIR.iterdir(), key=lambda p: p.name):
        if entry.is_file() and entry.name not in _KNOWN_DATA_FILES:
            return entry
    return None


def _inspect_full_table(path: Path) -> None:
    """Read only the header and first three rows of a candidate full table and print the stop notice."""
    preview = pd.read_csv(path, nrows=3, encoding="utf-8", keep_default_na=False)
    print(f"Full study table candidate: {path.name}")
    print(f"Columns ({len(preview.columns)}): {list(preview.columns)}")
    print(preview.to_string(index=False))
    print(_STOP_SPLICING_MESSAGE)


def select_failures(
    pred: pd.DataFrame,
    guesses: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Select Arm 2's 5 worst misses and 5 worst false alarms and attach curated guesses."""
    guess_map: Mapping[str, str] = CURATED_GUESSES if guesses is None else guesses

    work = pred.copy()
    work["position"] = work["position"].astype("int64")
    work["function_score"] = work["function_score"].astype("float64")
    work["esm_score"] = work["esm_score"].astype("float64")
    work["harm_arm2"] = work["harm_arm2"].astype("float64")
    work["blosum_score"] = work["blosum_score"].astype("int64")
    work["entropy"] = work["entropy"].astype("float64")
    work["flag_arm1"] = work["flag_arm1"].astype("int64")
    work["flag_arm2"] = work["flag_arm2"].astype("int64")
    work["flag_arm3"] = work["flag_arm3"].astype("int64")

    misses = (
        work[
            (work["findlay_class"] == config.HARMFUL_CLASS)
            & (work["flag_arm2"] == 0)
        ]
        .sort_values(
            by=["harm_arm2", "function_score", "mutant"],
            ascending=[True, True, True],
            kind="mergesort",
        )
        .head(5)
        .copy()
    )
    misses["kind"] = "miss"
    misses["rank"] = list(range(1, len(misses) + 1))

    false_alarms = (
        work[
            (work["findlay_class"] == config.HARMLESS_CLASS)
            & (work["flag_arm2"] == 1)
        ]
        .sort_values(
            by=["harm_arm2", "function_score", "mutant"],
            ascending=[False, False, True],
            kind="mergesort",
        )
        .head(5)
        .copy()
    )
    false_alarms["kind"] = "false_alarm"
    false_alarms["rank"] = list(range(1, len(false_alarms) + 1))

    combined = pd.concat([misses, false_alarms], ignore_index=True)
    if len(combined) == 0:
        return pd.DataFrame(columns=FAILURE_COLUMNS)

    combined["rank"] = combined["rank"].astype("int64")
    combined["guess_why"] = [
        str(guess_map.get(str(m), "")) for m in combined["mutant"]
    ]
    return combined[FAILURE_COLUMNS].reset_index(drop=True)


def format_failures_markdown(failures_df: pd.DataFrame) -> str:
    """Render the failures table, the two-sentence explained wrong case, and the splicing note."""
    header = "| " + " | ".join(FAILURE_COLUMNS) + " |"
    separator = "| " + " | ".join("---" for _ in FAILURE_COLUMNS) + " |"
    table_lines = [header, separator]

    int_cols = {"rank", "position", "blosum_score", "flag_arm1", "flag_arm3"}
    float_cols = {"function_score", "esm_score", "harm_arm2", "entropy"}

    for _, row in failures_df.iterrows():
        cells: list[str] = []
        for col in FAILURE_COLUMNS:
            val = row[col]
            if col in int_cols:
                cells.append(str(int(val)))
            elif col in float_cols:
                cells.append(f"{float(val):.6f}")
            else:
                cells.append(str(val))
        table_lines.append("| " + " | ".join(cells) + " |")

    lines: list[str] = [
        "# Stage 6: Arm 2 Worst Misses and False Alarms (End Region)",
        "",
        *table_lines,
        "",
        "## Explained wrong case",
        "",
        (
            "**Miss highlighted: `A1823S` (rank 1 miss).** "
            "Guess: At position 1823, sequence entropy is high (`2.572233`) and "
            "swapping `A` to `S` (`swap_kind=other`) has a positive BLOSUM62 score "
            "(`1`) and a positive ESM-2 score (`0.746405`, `harm_arm2=-0.746405`, "
            "`band_arm2=Looks harmless`), so Arm 2 and both other arms left it "
            "unflagged (`flag_arm1=0`, `flag_arm3=0`) even though the lab measured "
            "severe loss of function (`function_score=-3.642831`, `LOF`). "
            "Guess: Because the lab edited DNA rather than protein directly, the "
            "`A1823S` DNA change may have disrupted RNA splicing or a local requirement "
            "that a protein-only sequence model cannot observe."
        ),
        "",
        f"**Splicing note:** {_NO_FULL_TABLE_MESSAGE}",
        "",
    ]
    return "\n".join(lines)


def _write_failures_outputs(failures_df: pd.DataFrame) -> str:
    """Write outputs/results/failures.csv and outputs/results/failures.md and return the Markdown text."""
    results_dir = config.RESULTS_DIR
    results_dir.mkdir(parents=True, exist_ok=True)

    csv_path = results_dir / "failures.csv"
    with open(csv_path, "w", encoding="utf-8", newline="") as handle:
        failures_df.to_csv(handle, index=False, lineterminator="\n")

    md_text = format_failures_markdown(failures_df)
    md_path = results_dir / "failures.md"
    with open(md_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(md_text)

    return md_text


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for Stage 6 worst misses and worst false alarms analysis."""
    config.set_seeds()
    parser = argparse.ArgumentParser(
        description="List Arm 2's five worst misses and five worst false alarms on the end region."
    )
    parser.add_argument(
        "--full-table",
        type=str,
        default=None,
        metavar="PATH",
        help="Optional path to the full Findlay 2018 DNA/RNA table for splicing inspection.",
    )
    args = parser.parse_args(argv)

    if not load.test_run_done():
        print(_NOT_DONE_REFUSAL_MESSAGE, file=sys.stderr)
        return 1

    try:
        full_table_path = _find_optional_full_table(args.full_table)
        if full_table_path is not None:
            _inspect_full_table(full_table_path)
            return 0

        print(_NO_FULL_TABLE_MESSAGE)

        pred_path = config.RESULTS_DIR / "predictions_end.csv"
        if not pred_path.is_file():
            raise FileNotFoundError(
                f"{pred_path} not found. Run the Stage 5 evaluation first."
            )

        pred_df = pd.read_csv(pred_path, encoding="utf-8", keep_default_na=False)
        failures_df = select_failures(pred_df)
        _write_failures_outputs(failures_df)
        print(failures_df.to_string(index=False))
        return 0
    except (FileNotFoundError, ValueError, TypeError, KeyError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
