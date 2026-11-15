"""Unit tests across Stages 1–3 using synthetic data, temporary paths, and monkeypatching."""

import math
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import baselines, config, cutoffs, esm_score, load


def _write_log_lines(tmp_path: Path, filename: str, entries: list[tuple[str, str]]) -> Path:
    """Write synthetic pipe-delimited run-log lines into a temporary file."""
    path = tmp_path / filename
    lines = [
        f"2026-11-15T09:05:1{i}-06:00 | run=1 | status={status} | kind=first | reason={reason}\n"
        for i, (status, reason) in enumerate(entries)
    ]
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _assert_log_cases(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    predicate: Callable[..., bool],
    cases: list[tuple[list[tuple[str, str]] | None, bool]],
) -> None:
    """Check a run-log predicate via direct path argument and monkeypatched RUN_LOG_PATH."""
    for idx, (entries, expected) in enumerate(cases):
        if entries is None:
            log_path = tmp_path / f"missing_{idx}.txt"
        else:
            log_path = _write_log_lines(tmp_path, f"log_{idx}.txt", entries)
        assert predicate(log_path) is expected
        monkeypatch.setattr(config, "RUN_LOG_PATH", log_path)
        assert predicate() is expected


def test_parse_mutant_good() -> None:
    """Valid mutation strings parse into (wild_type, position, new_letter)."""
    assert load.parse_mutant("C61G") == ("C", 61, "G")
    assert load.parse_mutant("S1631A") == ("S", 1631, "A")
    assert load.parse_mutant("M1I") == ("M", 1, "I")


def test_parse_mutant_bad() -> None:
    """Malformed mutation strings raise ValueError and None raises TypeError."""
    bad_inputs = [
        "",
        "c61g",
        "C61",
        "61G",
        "C61GG",
        "X61G",
        "C0G",
        "C61G ",
        "C61C",
        "p.C61G",
    ]
    for bad in bad_inputs:
        with pytest.raises(ValueError):
            load.parse_mutant(bad)

    with pytest.raises(TypeError):
        load.parse_mutant(None)  # type: ignore[arg-type]


def test_region_of() -> None:
    """Boundary positions map to start/end and out-of-range positions raise ValueError."""
    assert load.region_of(1) == "start"
    assert load.region_of(101) == "start"
    assert load.region_of(1631) == "end"
    assert load.region_of(1855) == "end"
    for pos in (0, 102, 1630, 1856):
        with pytest.raises(ValueError):
            load.region_of(pos)


def test_class_of_boundaries() -> None:
    """Strict boundary cutoffs classify scores into LOF, INT, and FUNC."""
    assert load.class_of(-1.329) == "LOF"
    assert load.class_of(-1.328) == "INT"
    assert load.class_of(-1.0) == "INT"
    assert load.class_of(-0.748) == "INT"
    assert load.class_of(-0.747) == "FUNC"


def test_sha256_of_file(tmp_path: Path) -> None:
    """sha256_of_file returns the expected digest for a file containing b'abc'."""
    sample = tmp_path / "abc.bin"
    sample.write_bytes(b"abc")
    assert (
        config.sha256_of_file(sample)
        == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_test_attempted(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """test_attempted returns True for started, failed, or success logs, False otherwise."""
    cases: list[tuple[list[tuple[str, str]] | None, bool]] = [
        (None, False),
        ([("refused", "hash check"), ("refused", "dirty tree")], False),
        ([("started", "first run")], True),
        ([("failed", "runtime error")], True),
        ([("success", "first run")], True),
    ]
    _assert_log_cases(tmp_path, monkeypatch, load.test_attempted, cases)


def test_test_run_done(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """test_run_done returns True only when a log field is exactly status=success."""
    cases: list[tuple[list[tuple[str, str]] | None, bool]] = [
        (None, False),
        ([("refused", "hash check")], False),
        ([("started", "first run")], False),
        ([("failed", "runtime error")], False),
        ([("success", "first run")], True),
        ([("failed", "still failed after success")], False),
    ]
    _assert_log_cases(tmp_path, monkeypatch, load.test_run_done, cases)


def test_load_labels_end_locked(monkeypatch: pytest.MonkeyPatch) -> None:
    """End-region labels raise PermissionError when locked; start tuning returns 575 rows."""
    with pytest.raises(PermissionError):
        load.load_labels("end", "tuning")

    monkeypatch.setattr(load, "test_run_done", lambda *a, **kw: False)
    with pytest.raises(PermissionError):
        load.load_labels("end", "post_test")

    monkeypatch.setattr(load, "test_attempted", lambda *a, **kw: False)
    with pytest.raises(PermissionError):
        load.load_labels("end", "evaluate")

    start_df = load.load_labels("start", "tuning")
    assert len(start_df) == 575
    assert "proteingym_bin" not in start_df.columns


def test_catch_cutoff_worked_example() -> None:
    """Section 8.3 worked example 1 matches cutoffs and caught/false-alarm counts."""
    harm_scores = [0.9, 0.7, 0.5, 0.5, 0.5, 0.1]
    is_harmful = [1, 0, 1, 0, 1, 0]

    assert cutoffs.choose_catch_cutoff(harm_scores, is_harmful, 90) == 0.5
    assert cutoffs.choose_catch_cutoff(harm_scores, is_harmful, 60) == 0.5
    assert cutoffs.choose_catch_cutoff(harm_scores, is_harmful, 30) == 0.9

    counts_05 = cutoffs.count_outcomes(harm_scores, is_harmful, 0.5)
    assert counts_05 == {
        "n_harmful": 3,
        "n_harmless": 3,
        "harmful_caught": 3,
        "harmful_missed": 0,
        "false_alarms": 2,
        "harmless_cleared": 1,
    }

    counts_09 = cutoffs.count_outcomes(harm_scores, is_harmful, 0.9)
    assert counts_09 == {
        "n_harmful": 3,
        "n_harmless": 3,
        "harmful_caught": 1,
        "harmful_missed": 2,
        "false_alarms": 0,
        "harmless_cleared": 3,
    }


def test_catch_cutoff_matches_brute_force() -> None:
    """choose_catch_cutoff equals max observed score catching >= ceil(target% * n_harmful)."""
    rng = np.random.default_rng(42)
    for _ in range(200):
        n = int(rng.integers(5, 31))
        scores = np.round(rng.uniform(-2.0, 2.0, size=n), 1)
        y = rng.integers(0, 2, size=n)
        if y.sum() == 0:
            y[0] = 1
        elif y.sum() == n:
            y[0] = 0

        n_harmful = int(y.sum())
        observed = {float(s) for s in scores}
        for target in (50, 80, 90):
            needed = -(-(target * n_harmful) // 100)
            expected = max(
                c
                for c in observed
                if int(y[cutoffs.flag(scores, c)].sum()) >= needed
            )
            got = cutoffs.choose_catch_cutoff(scores, y, target)
            assert isinstance(got, float)
            assert got == expected


def test_catch_cutoff_needed_counts() -> None:
    """Distinct harmful scores at target 90 catch 118 of 131 and 241 of 267."""
    for n_harmful, expected_caught in ((131, 118), (267, 241)):
        harmful_scores = np.arange(1, n_harmful + 1, dtype=float)
        harmless_scores = np.arange(-50, 0, dtype=float)
        scores = np.concatenate([harmful_scores, harmless_scores])
        is_harmful = np.concatenate(
            [np.ones(n_harmful, dtype=int), np.zeros(50, dtype=int)]
        )

        cutoff = cutoffs.choose_catch_cutoff(scores, is_harmful, 90)
        counts = cutoffs.count_outcomes(scores, is_harmful, cutoff)
        assert counts["harmful_caught"] == expected_caught
        assert counts["harmful_missed"] == n_harmful - expected_caught


def test_false_alarm_cutoff_worked_example() -> None:
    """Section 8.3 worked example 2 returns 6.5 (a harmful row's score)."""
    harmless = [8, 6, 5, 5, 4, 3, 3, 2, 1, 0]
    harmful = [9, 7, 6.5, 2]
    scores = harmless + harmful
    is_harmful = [0] * len(harmless) + [1] * len(harmful)

    cutoff_b = cutoffs.choose_false_alarm_cutoff(scores, is_harmful, 10)
    assert isinstance(cutoff_b, float)
    assert cutoff_b == 6.5
    assert cutoffs.count_outcomes(scores, is_harmful, cutoff_b)["false_alarms"] == 1


def test_false_alarm_cutoff_none() -> None:
    """Second case in worked example 2 returns None when no score meets target."""
    harmless = [5, 5, 5, 1]
    harmful = [4, 2]
    scores = harmless + harmful
    is_harmful = [0] * len(harmless) + [1] * len(harmful)

    assert cutoffs.choose_false_alarm_cutoff(scores, is_harmful, 10) is None


def test_band_verdict() -> None:
    """Section 8.3 band_verdict and has_band worked examples with and without a band."""
    assert cutoffs.has_band(2.0, 5.0) is True
    for score, expected in (
        (5.0, "Looks harmful"),
        (7.1, "Looks harmful"),
        (4.99, "Can't tell"),
        (2.0, "Can't tell"),
        (1.99, "Looks harmless"),
    ):
        assert cutoffs.band_verdict(score, 2.0, 5.0) == expected

    for cutoff_b in (None, 2.0, 1.0):
        assert cutoffs.has_band(2.0, cutoff_b) is False
        assert cutoffs.band_verdict(2.0, 2.0, cutoff_b) == "Looks harmful"
        assert cutoffs.band_verdict(1.99, 2.0, cutoff_b) == "Looks harmless"


def test_flag_is_inclusive() -> None:
    """flag uses >= so scores equal to cutoff are flagged."""
    flags = cutoffs.flag([1.0, 2.0, 3.0], 2.0)
    assert isinstance(flags, np.ndarray)
    assert flags.dtype == bool
    assert flags.tolist() == [False, True, True]


def test_count_outcomes_arm0() -> None:
    """With cutoff=None, nothing is flagged and all harmful rows are missed."""
    harm_scores = [4.0, 1.5, 0.0, -1.0, 3.2]
    is_harmful = [1, 1, 0, 0, 0]

    counts = cutoffs.count_outcomes(harm_scores, is_harmful, None)
    assert counts == {
        "n_harmful": 2,
        "n_harmless": 3,
        "harmful_caught": 0,
        "harmful_missed": 2,
        "false_alarms": 0,
        "harmless_cleared": 3,
    }
    assert counts["harmful_missed"] == counts["n_harmful"]
    assert counts["false_alarms"] == 0
    assert counts["harmless_cleared"] == counts["n_harmless"]

    rate_summary = cutoffs.rates(counts)
    assert rate_summary == {
        "catch_rate_percent": 0.0,
        "false_alarm_rate_percent": 0.0,
        "accuracy_percent": 60.0,
    }


def test_blosum_loads() -> None:
    """BLOSUM62 table, blosum_score, add_arm1, and arm0_flags work on synthetic inputs."""
    table = baselines.load_blosum62()
    assert len(table) == 400
    assert table[("C", "G")] == -3
    assert table[("D", "E")] == 2
    assert table[("M", "I")] == 1
    assert table[("W", "W")] == 11

    for a in config.AA_ORDER:
        for b in config.AA_ORDER:
            assert (a, b) in table
            assert type(table[(a, b)]) is int
            assert table[(a, b)] == table[(b, a)]

    assert baselines.blosum_score("C", "G") == -3
    assert baselines.blosum_score("D", "E") == 2
    assert baselines.blosum_score("M", "I") == 1
    assert baselines.blosum_score("W", "W") == 11
    with pytest.raises(KeyError):
        baselines.blosum_score("X", "A")

    df = pd.DataFrame({"wild_type": ["C", "D", "A"], "new_letter": ["G", "E", "T"]})
    enriched = baselines.add_arm1(df)
    assert list(df.columns) == ["wild_type", "new_letter"]
    assert enriched["blosum_score"].tolist() == [-3, 2, 0]
    assert pd.api.types.is_integer_dtype(enriched["blosum_score"])
    assert enriched["harm_arm1"].tolist() == [3.0, -2.0, 0.0]
    assert pd.api.types.is_float_dtype(enriched["harm_arm1"])

    flags = baselines.arm0_flags(4)
    assert isinstance(flags, np.ndarray)
    assert flags.dtype == bool
    assert flags.shape == (4,)
    assert not flags.any()


def test_harm_arm1_no_negative_zero() -> None:
    """harm_arm1 returns +0.0 (never -0.0) when BLOSUM62 score is 0."""
    zero_harm = baselines.harm_arm1("A", "T")
    assert isinstance(zero_harm, float)
    assert zero_harm == 0.0
    assert str(zero_harm) == "0.0"
    assert not math.copysign(1.0, zero_harm) < 0

    assert baselines.harm_arm1("C", "G") == 3.0


def test_token_index_worked_examples() -> None:
    """Section 8.5 worked and boundary positions map to expected windows and token indices."""
    cases = [
        (61, "A", 61),
        (1775, "B", 412),
        (1, "A", 1),
        (101, "A", 101),
        (1631, "B", 268),
        (1855, "B", 492),
    ]
    for pos, expected_window, expected_index in cases:
        assert esm_score.window_of(pos) == expected_window
        assert esm_score.token_index_of(pos) == expected_index

    for bad_pos in (102, 1630):
        with pytest.raises(ValueError):
            esm_score.window_of(bad_pos)
        with pytest.raises(ValueError):
            esm_score.token_index_of(bad_pos)


def test_window_sequence() -> None:
    """Windows A and B slice the wild-type sequence to expected lengths, ends, and residues."""
    seq = load.load_wild_type()
    win_a = esm_score.window_sequence(seq, "A")
    win_b = esm_score.window_sequence(seq, "B")

    assert len(win_a) == 500
    assert win_a.startswith("MDLSA")
    assert win_a.endswith("RPLTN")
    assert win_a[60] == "C"

    assert len(win_b) == 500
    assert win_b.startswith("NLGEA")
    assert win_b.endswith("PHSHY")
    assert win_b[411] == "M"


def test_only_refuses_end_region(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """--only refuses end-region mutants before the test run without loading model or device."""
    monkeypatch.setattr(load, "test_run_done", lambda *a, **kw: False)

    def _fail(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("should not load model or choose device")

    monkeypatch.setattr(esm_score, "load_model", _fail)
    monkeypatch.setattr(esm_score, "choose_device", _fail)

    exit_code = esm_score.main(["--only", "S1631A"])
    captured = capsys.readouterr()

    assert exit_code == 1
    assert (
        "Refusing to score a single end-region mutation before the test run. "
        "Use --letter-test-only to check end-region indexing."
    ) in captured.err

