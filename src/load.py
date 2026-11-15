"""Single data-access gate for reading, validating, and controlling access to BRCA1 inputs and ESM score caches."""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src import config

_MUTANT_RE = re.compile(
    rf"^([{config.AA_ORDER}])([1-9][0-9]{{0,3}})([{config.AA_ORDER}])$"
)

_MUTATION_COLUMNS: list[str] = [
    "mutant",
    "wild_type",
    "position",
    "new_letter",
    "region",
]
_LABEL_COLUMNS: list[str] = _MUTATION_COLUMNS + ["function_score", "findlay_class"]
_SLIM_COLUMNS: list[str] = _LABEL_COLUMNS + ["proteingym_bin"]
_SLIM_DTYPES: dict[str, object] = {
    "mutant": str,
    "wild_type": str,
    "position": "int64",
    "new_letter": str,
    "region": str,
    "function_score": "float64",
    "findlay_class": str,
    "proteingym_bin": "int64",
}

_RAW_COLUMNS: list[str] = ["mutant", "mutated_sequence", "DMS_score", "DMS_score_bin"]
_RAW_DTYPES: dict[str, object] = {
    "mutant": str,
    "mutated_sequence": str,
    "DMS_score": "float64",
    "DMS_score_bin": "int64",
}

_SCORE_COLUMNS: list[str] = [
    "position",
    "region",
    "wild_type",
    "new_letter",
    "mutant",
    "window",
    "token_index",
    "token_at_index",
    "esm_score",
    "entropy",
    "is_lab_tested",
    "model_id",
    "transformers_version",
    "precision",
]
_SCORE_DTYPES: dict[str, object] = {
    "position": "int64",
    "region": str,
    "wild_type": str,
    "new_letter": str,
    "mutant": str,
    "window": str,
    "token_index": "int64",
    "token_at_index": str,
    "esm_score": "float64",
    "entropy": "float64",
    "is_lab_tested": "int64",
    "model_id": str,
    "transformers_version": str,
    "precision": str,
}


def parse_mutant(text: str) -> tuple[str, int, str]:
    """Parse strict mutation notation (e.g., C61G) into (wild_type, position, new_letter)."""
    if not isinstance(text, str):
        raise TypeError("mutant must be a string")
    match = _MUTANT_RE.fullmatch(text)
    if match is None:
        raise ValueError(
            f"Cannot parse mutant '{text}': expected one letter, a position, "
            "one letter, for example C61G."
        )
    wild_type, pos_str, new_letter = match.groups()
    if wild_type == new_letter:
        raise ValueError(f"Mutant '{text}' does not change the letter.")
    return wild_type, int(pos_str), new_letter


def region_of(position: int) -> str:
    """Return 'start' or 'end' for a 1-based position inside a lab-tested region."""
    start_lo, start_hi = config.START_REGION
    end_lo, end_hi = config.END_REGION
    if not isinstance(position, bool):
        if start_lo <= position <= start_hi:
            return "start"
        if end_lo <= position <= end_hi:
            return "end"
    raise ValueError(
        f"Position {position} is outside the two lab-tested regions "
        f"({start_lo}-{start_hi} and {end_lo}-{end_hi})."
    )


def class_of(function_score: float) -> str:
    """Map a numeric function score to 'LOF', 'FUNC', or 'INT'."""
    if not np.isfinite(function_score):
        raise ValueError("function_score must be finite.")
    if function_score < config.LOF_BELOW:
        return config.HARMFUL_CLASS
    if function_score > config.FUNC_ABOVE:
        return config.HARMLESS_CLASS
    return config.EXCLUDED_CLASS


def load_wild_type() -> str:
    """Read and validate the 1,863-letter BRCA1 wild-type sequence from FASTA."""
    with open(config.FASTA_PATH, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if not line.startswith(">")]
    seq = "".join(lines)
    if len(seq) != config.PROTEIN_LENGTH:
        raise ValueError(
            f"Wild-type sequence check failed: expected length {config.PROTEIN_LENGTH}, got {len(seq)}."
        )
    if not seq.startswith(config.SEQUENCE_START):
        raise ValueError(
            f"Wild-type sequence check failed: does not start with {config.SEQUENCE_START}."
        )
    invalid = sorted(set(seq) - set(config.AA_ORDER))
    if invalid:
        raise ValueError(
            f"Wild-type sequence check failed: unexpected letters {''.join(invalid)}."
        )
    return seq


def _read_slim_csv() -> pd.DataFrame:
    """Read the slim CSV with explicit dtypes and preserved string values."""
    return pd.read_csv(
        config.SLIM_CSV_PATH,
        dtype=_SLIM_DTYPES,
        keep_default_na=False,
    )


def load_mutations() -> pd.DataFrame:
    """Return a label-free copy of the 1,837 mutation rows from the slim CSV."""
    df = _read_slim_csv()
    return df[_MUTATION_COLUMNS].copy()


def class_counts() -> dict[str, dict[str, int]]:
    """Return LOF, INT, and FUNC counts per region without printing or returning rows."""
    df = _read_slim_csv()
    classes = (config.HARMFUL_CLASS, config.EXCLUDED_CLASS, config.HARMLESS_CLASS)
    counts: dict[str, dict[str, int]] = {}
    for region in ("start", "end"):
        sub = df[df["region"] == region]
        counts[region] = {
            cls: int((sub["findlay_class"] == cls).sum()) for cls in classes
        }
    return counts


def _read_log_parts(run_log_path: Path | None = None) -> list[str]:
    """Read stripped pipe-separated tokens from all non-empty lines in the run log."""
    path = config.RUN_LOG_PATH if run_log_path is None else Path(run_log_path)
    if not path.is_file():
        return []
    parts: list[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            parts.extend(part.strip() for part in stripped.split("|"))
    return parts


def test_run_done(run_log_path: Path | None = None) -> bool:
    """Return True if the run log contains a part exactly equal to 'status=success'."""
    return "status=success" in _read_log_parts(run_log_path)


def test_attempted(run_log_path: Path | None = None) -> bool:
    """Return True if the run log contains started, success, or failed status."""
    attempted = {"status=started", "status=success", "status=failed"}
    return any(part in attempted for part in _read_log_parts(run_log_path))


def load_labels(region: str, purpose: str) -> pd.DataFrame:
    """Return labeled rows for one region, enforcing end-region access locks."""
    if region not in ("start", "end"):
        raise ValueError("region must be 'start' or 'end'.")
    if purpose not in ("tuning", "evaluate", "post_test"):
        raise ValueError("purpose must be 'tuning', 'evaluate' or 'post_test'.")

    if region == "end":
        if purpose == "tuning":
            raise PermissionError("End-region lab labels may not be used for tuning.")
        if purpose == "evaluate" and not test_attempted():
            raise PermissionError(
                "End-region lab labels may be loaded only inside the test run, "
                "after its started line is written."
            )
        if purpose == "post_test" and not test_run_done():
            raise PermissionError(
                "End-region lab labels are locked until the one test run has succeeded "
                "(outputs/results/run_log.txt has no status=success line)."
            )

    df = _read_slim_csv()
    return df.loc[df["region"] == region, _LABEL_COLUMNS].reset_index(drop=True)


def validate_scores(df: pd.DataFrame) -> None:
    """Validate structure, indexing, wild-type consistency, and metadata of ESM scores."""
    if list(df.columns) != _SCORE_COLUMNS:
        raise ValueError("esm_scores.csv check failed: unexpected columns.")
    if len(df) != config.N_SCORE_ROWS:
        raise ValueError(
            f"esm_scores.csv check failed: expected {config.N_SCORE_ROWS} rows, got {len(df)}."
        )
    if df["mutant"].nunique() != config.N_SCORE_ROWS:
        raise ValueError("esm_scores.csv check failed: mutant values are not unique.")
    if df["position"].nunique() != config.N_POSITIONS:
        raise ValueError(
            f"esm_scores.csv check failed: expected {config.N_POSITIONS} distinct positions."
        )
    if not df.groupby("position").size().eq(19).all():
        raise ValueError("esm_scores.csv check failed: each position must have 19 rows.")
    if df.isna().any().any():
        raise ValueError("esm_scores.csv check failed: nulls found.")
    if not np.isfinite(df["esm_score"]).all() or not np.isfinite(df["entropy"]).all():
        raise ValueError("esm_scores.csv check failed: non-finite esm_score or entropy.")
    if not (df["token_at_index"] == df["wild_type"]).all():
        raise ValueError("esm_scores.csv check failed: token_at_index != wild_type.")

    seq = load_wild_type()
    if not df["position"].between(1, len(seq)).all():
        raise ValueError("esm_scores.csv check failed: position out of sequence bounds.")
    expected_wt = df["position"].map(lambda pos: seq[int(pos) - 1])
    if not (df["wild_type"] == expected_wt).all():
        raise ValueError("esm_scores.csv check failed: wild_type does not match sequence.")

    start_lo, start_hi = config.START_REGION
    end_lo, end_hi = config.END_REGION
    in_start = df["position"].between(start_lo, start_hi)
    in_end = df["position"].between(end_lo, end_hi)
    if not (in_start | in_end).all():
        raise ValueError("esm_scores.csv check failed: position outside scored regions.")

    expected_region = np.where(in_start, "start", "end")
    if not (df["region"].to_numpy() == expected_region).all():
        raise ValueError("esm_scores.csv check failed: region mismatch.")

    expected_window = np.where(in_start, "A", "B")
    if not (df["window"].to_numpy() == expected_window).all():
        raise ValueError("esm_scores.csv check failed: window does not match formula.")

    expected_token_index = np.where(
        in_start,
        df["position"] - config.WINDOW_START["A"] + 1,
        df["position"] - config.WINDOW_START["B"] + 1,
    )
    if not (df["token_index"].to_numpy() == expected_token_index).all():
        raise ValueError("esm_scores.csv check failed: token_index does not match formula.")

    if not df["is_lab_tested"].isin([0, 1]).all():
        raise ValueError("esm_scores.csv check failed: is_lab_tested must be 0 or 1.")
    if int(df["is_lab_tested"].sum()) != config.N_ROWS:
        raise ValueError(
            f"esm_scores.csv check failed: is_lab_tested sum is not {config.N_ROWS}."
        )
    slim_mutants = set(load_mutations()["mutant"])
    tested_mutants = set(df.loc[df["is_lab_tested"] == 1, "mutant"])
    if tested_mutants != slim_mutants:
        raise ValueError("esm_scores.csv check failed: is_lab_tested mutants do not match slim CSV.")

    if not ((df["entropy"] >= 0.0) & (df["entropy"] <= 2.9958)).all():
        raise ValueError("esm_scores.csv check failed: entropy out of [0, 2.9958].")
    if not df.groupby("position")["entropy"].nunique().eq(1).all():
        raise ValueError("esm_scores.csv check failed: entropy is not constant within position.")

    for col in ("model_id", "transformers_version", "precision"):
        if df[col].nunique() != 1 or df[col].iloc[0] == "":
            raise ValueError(
                f"esm_scores.csv check failed: {col} must have one non-empty distinct value."
            )


def load_scores() -> pd.DataFrame:
    """Load and validate the cached ESM score table without any lab labels."""
    if not config.ESM_SCORES_PATH.is_file():
        raise FileNotFoundError(
            "outputs/esm_scores.csv not found. Run: python -m src.esm_score"
        )
    df = pd.read_csv(
        config.ESM_SCORES_PATH,
        dtype=_SCORE_DTYPES,
        keep_default_na=False,
    )
    validate_scores(df)
    return df


def lab_tested_scores(region: str | None = None) -> pd.DataFrame:
    """Return validated ESM scores filtered to lab-tested rows and optional region."""
    if region is not None and region not in ("start", "end"):
        raise ValueError("region must be 'start' or 'end'.")
    df = load_scores()
    tested = df[df["is_lab_tested"] == 1]
    if region is not None:
        tested = tested[tested["region"] == region]
    return tested.reset_index(drop=True)


def check_inputs() -> None:
    """Verify SHA-256 hashes, slim CSV schema, row consistency, counts, and boundaries."""
    hash_targets = (
        ("raw CSV", config.RAW_CSV_PATH, config.RAW_CSV_SHA256),
        ("slim CSV", config.SLIM_CSV_PATH, config.SLIM_CSV_SHA256),
        ("FASTA", config.FASTA_PATH, config.FASTA_SHA256),
    )
    for label, path, expected_sha in hash_targets:
        actual_sha = config.sha256_of_file(path)
        if actual_sha != expected_sha:
            raise ValueError(f"Input check failed: {label} SHA-256 mismatch.")
    print("Check 1/6 passed: SHA-256 hashes match for raw CSV, slim CSV, and FASTA.")

    try:
        df = _read_slim_csv()
    except (ValueError, pd.errors.ParserError) as exc:
        raise ValueError(f"Input check failed: slim CSV cannot be parsed ({exc}).") from exc

    if len(df) != config.N_ROWS:
        raise ValueError(
            f"Input check failed: slim CSV must have {config.N_ROWS} rows, got {len(df)}."
        )
    if list(df.columns) != _SLIM_COLUMNS:
        raise ValueError("Input check failed: slim CSV columns do not match Section 6.2.")
    if df.isna().any().any():
        raise ValueError("Input check failed: slim CSV contains nulls.")
    if df["mutant"].nunique() != config.N_ROWS:
        raise ValueError("Input check failed: slim CSV mutant values are not unique.")
    print("Check 2/6 passed: slim CSV has 1837 rows, 8 columns, no nulls, unique mutants.")

    seq = load_wild_type()
    bad_parse = 0
    bad_region = 0
    bad_class = 0
    bad_wt = 0
    for mutant, wt, pos, new_aa, region, score, fclass in zip(
        df["mutant"],
        df["wild_type"],
        df["position"],
        df["new_letter"],
        df["region"],
        df["function_score"],
        df["findlay_class"],
    ):
        pos_int = int(pos)
        try:
            if parse_mutant(mutant) != (wt, pos_int, new_aa):
                bad_parse += 1
        except (ValueError, TypeError):
            bad_parse += 1
        try:
            if region_of(pos_int) != region:
                bad_region += 1
        except ValueError:
            bad_region += 1
        try:
            if class_of(float(score)) != fclass:
                bad_class += 1
        except ValueError:
            bad_class += 1
        if not (1 <= pos_int <= len(seq)) or wt != seq[pos_int - 1]:
            bad_wt += 1

    if bad_parse:
        raise ValueError(
            f"Input check failed: parse_mutant mismatch in {bad_parse} rows."
        )
    if bad_region:
        raise ValueError(
            f"Input check failed: region mismatch in {bad_region} rows."
        )
    if bad_class:
        raise ValueError(
            f"Input check failed: findlay_class mismatch in {bad_class} rows."
        )
    if bad_wt:
        raise ValueError(
            f"Input check failed: wild_type sequence mismatch in {bad_wt} rows."
        )
    print("Check 3/6 passed: all 1837 rows match mutant parsing, region, class, and wild-type sequence.")

    if class_counts() != config.EXPECTED_COUNTS:
        raise ValueError("Input check failed: class_counts() does not match EXPECTED_COUNTS.")
    print("Check 4/6 passed: class counts match expected counts.")

    if df["position"].nunique() != config.N_POSITIONS:
        raise ValueError(
            f"Input check failed: expected {config.N_POSITIONS} distinct positions, "
            f"got {df['position'].nunique()}."
        )
    print("Check 5/6 passed: 326 distinct positions found.")

    at_cutoff = (df["function_score"] == config.LOF_BELOW) | (
        df["function_score"] == config.FUNC_ABOVE
    )
    if at_cutoff.any():
        raise ValueError(
            f"Input check failed: {int(at_cutoff.sum())} rows equal a class cutoff."
        )
    print("Check 6/6 passed: no function_score equals -1.328 or -0.748.")


def rebuild_and_compare() -> None:
    """Rebuild wild-type sequence and slim table from raw CSV and compare with disk files."""
    raw_df = pd.read_csv(
        config.RAW_CSV_PATH,
        dtype=_RAW_DTYPES,
        keep_default_na=False,
    )
    if list(raw_df.columns) != _RAW_COLUMNS:
        raise ValueError("Rebuild check failed: raw CSV columns do not match Section 6.1.")
    if len(raw_df) != config.N_ROWS:
        raise ValueError(
            f"Rebuild check failed: raw CSV expected {config.N_ROWS} rows, got {len(raw_df)}."
        )
    if raw_df.isna().any().any():
        raise ValueError("Rebuild check failed: raw CSV contains nulls.")
    if not raw_df["mutant"].is_unique:
        raise ValueError("Rebuild check failed: raw CSV mutant values are not unique.")

    parsed = [parse_mutant(m) for m in raw_df["mutant"]]
    first_seq = raw_df["mutated_sequence"].iloc[0]
    first_wt, first_pos, _ = parsed[0]
    if len(first_seq) != config.PROTEIN_LENGTH or not (1 <= first_pos <= config.PROTEIN_LENGTH):
        raise ValueError(
            f"Row {raw_df['mutant'].iloc[0]} disagrees with the rebuilt wild-type sequence."
        )
    wt_chars = list(first_seq)
    wt_chars[first_pos - 1] = first_wt
    rebuilt_seq = "".join(wt_chars)

    for mutant, mutated_seq, (wt, pos, new_aa) in zip(
        raw_df["mutant"], raw_df["mutated_sequence"], parsed
    ):
        if (
            len(mutated_seq) != len(rebuilt_seq)
            or not (1 <= pos <= len(rebuilt_seq))
            or rebuilt_seq[pos - 1] != wt
            or mutated_seq != rebuilt_seq[: pos - 1] + new_aa + rebuilt_seq[pos:]
        ):
            raise ValueError(f"Row {mutant} disagrees with the rebuilt wild-type sequence.")

    rebuilt_slim = pd.DataFrame(
        {
            "mutant": raw_df["mutant"],
            "wild_type": [wt for wt, _, _ in parsed],
            "position": pd.Series([pos for _, pos, _ in parsed], dtype="int64"),
            "new_letter": [new_aa for _, _, new_aa in parsed],
            "region": [region_of(pos) for _, pos, _ in parsed],
            "function_score": raw_df["DMS_score"].astype("float64"),
            "findlay_class": [class_of(float(s)) for s in raw_df["DMS_score"]],
            "proteingym_bin": raw_df["DMS_score_bin"].astype("int64"),
        },
        columns=_SLIM_COLUMNS,
    )

    disk_seq = load_wild_type()
    if rebuilt_seq != disk_seq:
        raise ValueError(
            "Rebuild check failed: rebuilt wild-type sequence does not match FASTA."
        )

    disk_slim = _read_slim_csv()
    if list(disk_slim.columns) != _SLIM_COLUMNS or len(disk_slim) != len(rebuilt_slim):
        raise ValueError("Rebuild check failed: slim CSV shape or columns do not match.")
    if disk_slim.isna().any().any():
        raise ValueError("Rebuild check failed: slim CSV contains nulls.")

    for col in _SLIM_COLUMNS:
        n_diff = int((rebuilt_slim[col].to_numpy() != disk_slim[col].to_numpy()).sum())
        if n_diff > 0:
            raise ValueError(
                f"Rebuild check failed: column {col} differs in {n_diff} rows."
            )

    print("REBUILD MATCHES: slim CSV and FASTA agree with the raw ProteinGym file.")


def _format_counts_table(counts: dict[str, dict[str, int]]) -> str:
    """Format per-region class counts as an aligned plain-text table."""
    lines = ["region  LOF  INT  FUNC  total"]
    for region in ("start", "end"):
        row = counts[region]
        lof = row[config.HARMFUL_CLASS]
        inter = row[config.EXCLUDED_CLASS]
        func = row[config.HARMLESS_CLASS]
        total = lof + inter + func
        lines.append(f"{region:<6}  {lof:>3}  {inter:>3}  {func:>4}  {total:>5}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """Run input integrity checks or raw-to-slim rebuild verification from the CLI."""
    config.set_seeds()
    parser = argparse.ArgumentParser(description="Validate BRCA1 input datasets.")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run input checks and print per-region class counts.",
    )
    parser.add_argument(
        "--rebuild-check",
        action="store_true",
        help="Rebuild slim CSV and FASTA from raw ProteinGym data and compare.",
    )
    args = parser.parse_args(argv)
    if not args.check and not args.rebuild_check:
        parser.print_help(sys.stderr)
        return 1
    try:
        if args.check:
            check_inputs()
            print(_format_counts_table(class_counts()))
        if args.rebuild_check:
            rebuild_and_compare()
    except (FileNotFoundError, ValueError, TypeError, PermissionError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
