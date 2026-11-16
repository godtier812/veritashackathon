"""Unit tests across Stages 1–3 using synthetic data, temporary paths, and monkeypatching."""

import json
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


def _patch_stage4_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Path]:
    """Point all Stage 4 output and log paths in config at a temporary directory."""
    outputs_dir = tmp_path / "outputs"
    results_dir = outputs_dir / "results"
    selfcheck_dir = outputs_dir / "selfcheck_start"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    selfcheck_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "OUTPUTS_DIR": outputs_dir,
        "ESM_SCORES_PATH": outputs_dir / "esm_scores.csv",
        "FROZEN_PATH": outputs_dir / "frozen.json",
        "RESULTS_DIR": results_dir,
        "RUN_LOG_PATH": results_dir / "run_log.txt",
        "SELFCHECK_DIR": selfcheck_dir,
    }
    for attr, path in paths.items():
        monkeypatch.setattr(config, attr, path)
    return paths


def _synthetic_frozen_dict(
    esm_sha256: str = "a" * 64,
    cutoff_a: float = 2.0,
    cutoff_b: float | None = 5.0,
    has_band_flag: bool = True,
) -> dict[str, object]:
    """Build a synthetic frozen configuration dictionary matching Section 6.5."""
    return {
        "schema_version": 1,
        "created_at": "2026-11-15T03:10:42-06:00",
        "code_commit": "0" * 40,
        "code_dirty": False,
        "seed": 42,
        "model_id": config.MODEL_ID,
        "inputs": {
            "esm_scores_path": "outputs/esm_scores.csv",
            "esm_scores_sha256": esm_sha256,
            "esm_scores_rows": config.N_SCORE_ROWS,
            "slim_csv_sha256": config.SLIM_CSV_SHA256,
        },
        "tuning_region": "start",
        "test_region": "end",
        "labels": {
            "harmful_class": config.HARMFUL_CLASS,
            "harmless_class": config.HARMLESS_CLASS,
            "excluded_class": config.EXCLUDED_CLASS,
            "lof_below": config.LOF_BELOW,
            "func_above": config.FUNC_ABOVE,
        },
        "tuning_counts": {"harmful": 131, "harmless": 385, "int_excluded": 59},
        "catch_target_percent": config.CATCH_TARGET_PERCENT,
        "false_alarm_target_percent": config.FALSE_ALARM_TARGET_PERCENT,
        "arms": {
            "arm0": {
                "name": config.ARM_NAMES["arm0"],
                "cutoff": None,
            },
            "arm1": {
                "name": config.ARM_NAMES["arm1"],
                "harm_score": "0 - blosum62(wild_type, new_letter)",
                "cutoff": 0.0,
                "tuning_caught": 120,
                "tuning_false_alarms": 217,
            },
            "arm2": {
                "name": config.ARM_NAMES["arm2"],
                "harm_score": "0 - esm_score",
                "cutoff": cutoff_a,
                "tuning_caught": 118,
                "tuning_false_alarms": 99,
            },
            "arm3": {
                "name": config.ARM_NAMES["arm3"],
                "harm_score": "predicted probability of harmful",
                "features": list(config.FEATURE_ORDER),
                "scaler_mean": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                "scaler_scale": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                "coefficients": [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                "intercept": 0.0,
                "cutoff": 0.9,
                "cutoff_source": (
                    "out-of-fold probabilities, StratifiedGroupKFold(n_splits=5, "
                    "shuffle=True, random_state=42), groups = position"
                ),
                "model_settings": {
                    "C": 1.0,
                    "solver": "lbfgs",
                    "max_iter": 100,
                    "random_state": 42,
                },
                "converged": True,
                "out_of_fold_caught": 118,
                "out_of_fold_false_alarms": 99,
                "tuning_caught": 120,
                "tuning_false_alarms": 90,
            },
        },
        "band": {
            "arm": "arm2",
            "cutoff_a": cutoff_a,
            "cutoff_b": cutoff_b,
            "has_band": has_band_flag,
        },
        "versions": {
            "python": "3.12.10",
            "numpy": "2.4.4",
            "pandas": "2.3.3",
            "scikit_learn": "1.8.0",
            "biopython": "1.88",
        },
    }


def _synthetic_six_row_feature_table() -> pd.DataFrame:
    """Build a 6-row synthetic feature table with known BLOSUM and ESM values."""
    from src import features

    raw_scores = pd.DataFrame(
        {
            "position": [1631, 1632, 1633, 1634, 1635, 1636],
            "region": ["end"] * 6,
            "wild_type": ["C", "C", "A", "D", "M", "S"],
            "new_letter": ["G", "G", "T", "E", "I", "A"],
            "mutant": ["C1631G", "C1632G", "A1633T", "D1634E", "M1635I", "S1636A"],
            "esm_score": [-6.0, -5.0, -3.0, -2.0, -0.5, 1.0],
            "entropy": [0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            "function_score": [-2.1, -1.8, 0.1, 0.2, -1.9, 0.3],
            "findlay_class": ["LOF", "LOF", "FUNC", "FUNC", "LOF", "FUNC"],
            "is_harmful": pd.Series([1, 1, 0, 0, 1, 0], dtype="Int64"),
        }
    )
    return features.build_feature_table(raw_scores)


def test_swap_kind_each_kind() -> None:
    """Section 7.9 swap_kind priority order and charge helper match all worked pairs."""
    from src import features

    assert features.charge("K") == 1
    assert features.charge("R") == 1
    assert features.charge("D") == -1
    assert features.charge("E") == -1
    assert features.charge("H") == 0
    assert features.charge("A") == 0

    cases = [
        ("L", "P", "to_proline"),
        ("C", "G", "loss_of_cysteine"),
        ("C", "R", "loss_of_cysteine"),
        ("C", "P", "to_proline"),
        ("S", "D", "charge_change"),
        ("L", "K", "charge_change"),
        ("H", "R", "charge_change"),
        ("K", "R", "other"),
        ("L", "S", "other"),
    ]
    for wild_type, new_letter, expected in cases:
        assert features.swap_kind(wild_type, new_letter) == expected


def test_features_c61g() -> None:
    """Section 8.6 worked example for C61G matches all biochemical features and raw matrix."""
    from src import features

    scores = pd.DataFrame(
        {
            "position": [61],
            "region": ["start"],
            "wild_type": ["C"],
            "new_letter": ["G"],
            "mutant": ["C61G"],
            "esm_score": [-9.123456],
            "entropy": [0.412345],
        }
    )
    table = features.build_feature_table(scores)
    row = table.iloc[0]

    assert int(row["blosum_score"]) == -3
    assert float(row["harm_arm1"]) == 3.0
    assert float(row["harm_arm2"]) == pytest.approx(9.123456)
    assert float(row["size_change"]) == pytest.approx(-48.4)
    assert float(row["hydrophobicity_change"]) == pytest.approx(-2.9)
    assert int(row["to_proline"]) == 0
    assert int(row["from_cysteine"]) == 1
    assert int(row["from_glycine"]) == 0
    assert row["swap_kind"] == "loss_of_cysteine"

    mat = features.feature_matrix(table)
    assert mat.shape == (1, 8)
    assert mat.dtype == np.float64
    assert np.allclose(
        mat[0],
        [-9.123456, -3.0, 0.412345, -48.4, -2.9, 0.0, 1.0, 0.0],
    )


def test_arm3_probability_matches_sklearn() -> None:
    """arm3_probability matches StandardScaler + LogisticRegression predict_proba[:, 1] within 1e-9."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    from src import features

    rng = np.random.default_rng(42)
    x = rng.normal(loc=1.2, scale=2.3, size=(80, 8))
    y = rng.integers(0, 2, size=80)
    y[0], y[1] = 0, 1

    pipeline = make_pipeline(StandardScaler(), LogisticRegression(random_state=42))
    pipeline.fit(x, y)
    scaler: StandardScaler = pipeline.named_steps["standardscaler"]
    model: LogisticRegression = pipeline.named_steps["logisticregression"]

    manual = features.arm3_probability(
        x,
        scaler.mean_,
        scaler.scale_,
        model.coef_[0],
        float(model.intercept_[0]),
    )
    expected = pipeline.predict_proba(x)[:, 1]
    assert manual.shape == (80,)
    assert np.max(np.abs(manual - expected)) < 1e-9

    with pytest.raises(ValueError):
        features.arm3_probability(
            x[:, :7],
            scaler.mean_,
            scaler.scale_,
            model.coef_[0],
            float(model.intercept_[0]),
        )


def test_fit_arm3_groups_do_not_leak(monkeypatch: pytest.MonkeyPatch) -> None:
    """fit_arm3 on 40 groups yields finite out-of-fold probabilities with disjoint group folds."""
    import sklearn.model_selection as model_selection

    from src import features

    real_sgkf = model_selection.StratifiedGroupKFold
    folds_checked: list[tuple[set[int], set[int]]] = []

    class _CheckingSGKFold:
        def __init__(self, *args: object, **kwargs: object) -> None:
            self._inner = real_sgkf(*args, **kwargs)

        def split(
            self,
            x_arr: np.ndarray,
            y_arr: np.ndarray,
            groups_arr: np.ndarray,
        ):
            for train_idx, test_idx in self._inner.split(x_arr, y_arr, groups_arr):
                train_groups = {int(g) for g in np.asarray(groups_arr)[train_idx]}
                test_groups = {int(g) for g in np.asarray(groups_arr)[test_idx]}
                assert train_groups.isdisjoint(test_groups)
                folds_checked.append((train_groups, test_groups))
                yield train_idx, test_idx

    monkeypatch.setattr(model_selection, "StratifiedGroupKFold", _CheckingSGKFold)
    if hasattr(features, "StratifiedGroupKFold"):
        monkeypatch.setattr(features, "StratifiedGroupKFold", _CheckingSGKFold)

    rng = np.random.default_rng(42)
    groups = np.repeat(np.arange(1, 41, dtype=int), 3)
    x = rng.normal(size=(len(groups), 8))
    y = np.array([idx % 2 for idx in range(len(groups))], dtype=int)

    fitted = features.fit_arm3(x, y, groups)
    oof = np.asarray(fitted["out_of_fold_probabilities"], dtype=np.float64)

    assert len(folds_checked) == config.N_FOLDS
    assert oof.shape == (len(groups),)
    assert not np.isnan(oof).any()
    assert ((oof >= 0.0) & (oof <= 1.0)).all()
    assert len(fitted["scaler_mean"]) == 8
    assert len(fitted["scaler_scale"]) == 8
    assert len(fitted["coefficients"]) == 8
    assert isinstance(fitted["cutoff"], float)
    assert isinstance(fitted["converged"], bool)


def test_apply_frozen_arms_synthetic() -> None:
    """apply_frozen_arms produces hand-computed flags, probabilities, and band verdicts on 6 rows."""
    from src import evaluate

    table = _synthetic_six_row_feature_table()
    frozen = _synthetic_frozen_dict(cutoff_a=2.0, cutoff_b=5.0, has_band_flag=True)

    pred = evaluate.apply_frozen_arms(table, frozen)
    expected_harm3 = [1.0 / (1.0 + math.exp(-h)) for h in [6.0, 5.0, 3.0, 2.0, 0.5, -1.0]]

    assert np.allclose(pred["harm_arm3"].to_numpy(dtype=float), expected_harm3, atol=1e-9)
    assert [int(v) for v in pred["flag_arm1"].tolist()] == [1, 1, 1, 0, 0, 0]
    assert [int(v) for v in pred["flag_arm2"].tolist()] == [1, 1, 1, 1, 0, 0]
    assert [int(v) for v in pred["flag_arm3"].tolist()] == [1, 1, 1, 0, 0, 0]
    assert pred["band_arm2"].tolist() == [
        "Looks harmful",
        "Looks harmful",
        "Can't tell",
        "Can't tell",
        "Looks harmless",
        "Looks harmless",
    ]


def test_headline_table_synthetic() -> None:
    """headline_table matches hand-computed outcome counts and rates across all four arms."""
    from src import evaluate, features

    table = _synthetic_six_row_feature_table()
    int_row = features.build_feature_table(
        pd.DataFrame(
            {
                "position": [1637],
                "region": ["end"],
                "wild_type": ["C"],
                "new_letter": ["G"],
                "mutant": ["C1637G"],
                "esm_score": [-8.0],
                "entropy": [0.3],
                "function_score": [-1.0],
                "findlay_class": ["INT"],
                "is_harmful": pd.Series([pd.NA], dtype="Int64"),
            }
        )
    )
    combined = pd.concat([table, int_row], ignore_index=True)
    frozen = _synthetic_frozen_dict(cutoff_a=2.0, cutoff_b=5.0, has_band_flag=True)
    pred = evaluate.apply_frozen_arms(combined, frozen)
    headline = evaluate.headline_table(pred, frozen)

    assert headline["arm"].tolist() == ["arm0", "arm1", "arm2", "arm3"]
    expected_rows = {
        "arm0": (3, 3, 0, 3, 0, 3, 0.0, 0.0, 50.0),
        "arm1": (3, 3, 2, 1, 1, 2, 66.7, 33.3, 66.7),
        "arm2": (3, 3, 2, 1, 2, 1, 66.7, 66.7, 50.0),
        "arm3": (3, 3, 2, 1, 1, 2, 66.7, 33.3, 66.7),
    }
    for _, row in headline.iterrows():
        arm = str(row["arm"])
        exp = expected_rows[arm]
        assert int(row["n_harmful"]) == exp[0]
        assert int(row["n_harmless"]) == exp[1]
        assert int(row["harmful_caught"]) == exp[2]
        assert int(row["harmful_missed"]) == exp[3]
        assert int(row["false_alarms"]) == exp[4]
        assert int(row["harmless_cleared"]) == exp[5]
        assert int(row["harmful_caught"]) + int(row["harmful_missed"]) == int(row["n_harmful"])
        assert int(row["false_alarms"]) + int(row["harmless_cleared"]) == int(row["n_harmless"])
        assert float(row["catch_rate_percent"]) == pytest.approx(exp[6])
        assert float(row["false_alarm_rate_percent"]) == pytest.approx(exp[7])
        assert float(row["accuracy_percent"]) == pytest.approx(exp[8])


def test_auroc_direction() -> None:
    """auroc_table returns 1.0 when harmful rows rank above harmless rows and 0.0 when reversed."""
    from src import evaluate

    pred_perfect = pd.DataFrame(
        {
            "findlay_class": ["LOF", "LOF", "LOF", "FUNC", "FUNC", "FUNC"],
            "is_harmful": pd.Series([1, 1, 1, 0, 0, 0], dtype="Int64"),
            "harm_arm1": [3.0, 2.0, 1.0, 0.0, -1.0, -2.0],
            "harm_arm2": [9.0, 8.0, 7.0, 2.0, 1.0, 0.5],
            "harm_arm3": [0.95, 0.90, 0.85, 0.25, 0.15, 0.05],
        }
    )
    auroc_high = evaluate.auroc_table(pred_perfect).set_index("arm")["auroc"].to_dict()
    assert auroc_high["arm0"] == pytest.approx(0.5)
    assert auroc_high["arm1"] == pytest.approx(1.0)
    assert auroc_high["arm2"] == pytest.approx(1.0)
    assert auroc_high["arm3"] == pytest.approx(1.0)

    pred_reversed = pred_perfect.copy()
    pred_reversed["harm_arm1"] = -pred_perfect["harm_arm1"]
    pred_reversed["harm_arm2"] = -pred_perfect["harm_arm2"]
    pred_reversed["harm_arm3"] = 1.0 - pred_perfect["harm_arm3"]

    auroc_low = evaluate.auroc_table(pred_reversed).set_index("arm")["auroc"].to_dict()
    assert auroc_low["arm0"] == pytest.approx(0.5)
    assert auroc_low["arm1"] == pytest.approx(0.0)
    assert auroc_low["arm2"] == pytest.approx(0.0)
    assert auroc_low["arm3"] == pytest.approx(0.0)


def test_bootstrap_resamples_positions(monkeypatch: pytest.MonkeyPatch) -> None:
    """bootstrap_auroc_differences is deterministic, resamples positions with duplicates, and skips 1-class draws."""
    from sklearn import metrics as sk_metrics

    from src import evaluate

    rows: list[dict[str, object]] = []
    counts_per_pos = [1, 2, 3, 4, 1, 2, 3, 4, 2, 5]
    rng = np.random.default_rng(42)
    for idx, n_rows in enumerate(counts_per_pos):
        pos = 1631 + idx
        for r_idx in range(n_rows):
            harmful = 1 if (idx + r_idx) % 2 == 0 else 0
            rows.append(
                {
                    "position": pos,
                    "findlay_class": "LOF" if harmful else "FUNC",
                    "is_harmful": harmful,
                    "harm_arm1": float(harmful * 2.0 + rng.normal(0.0, 0.5)),
                    "harm_arm2": float(harmful * 3.0 + rng.normal(0.0, 0.5)),
                    "harm_arm3": float(np.clip(0.2 + 0.6 * harmful + rng.normal(0.0, 0.1), 0.01, 0.99)),
                }
            )
    pred = pd.DataFrame(rows)

    res_a = evaluate.bootstrap_auroc_differences(pred, n_resamples=60, seed=42)
    res_b = evaluate.bootstrap_auroc_differences(pred, n_resamples=60, seed=42)
    assert res_a == res_b
    assert res_a["n_positions"] == 10
    assert res_a["n_resamples_requested"] == 60
    assert res_a["n_resamples_used"] + res_a["n_resamples_skipped_one_class"] == 60

    real_roc_auc = sk_metrics.roc_auc_score
    seen_sizes: list[int] = []

    def _recording_auc(y_true: object, y_score: object) -> float:
        seen_sizes.append(len(np.asarray(y_true)))
        return float(real_roc_auc(y_true, y_score))

    monkeypatch.setattr(sk_metrics, "roc_auc_score", _recording_auc)
    if hasattr(evaluate, "roc_auc_score"):
        monkeypatch.setattr(evaluate, "roc_auc_score", _recording_auc)

    evaluate.bootstrap_auroc_differences(pred, n_resamples=60, seed=42)
    assert len(set(seen_sizes)) > 1

    one_class_pred = pd.DataFrame(
        {
            "position": [1631, 1631, 1632, 1632, 1633],
            "findlay_class": ["FUNC", "FUNC", "FUNC", "FUNC", "LOF"],
            "is_harmful": [0, 0, 0, 0, 1],
            "harm_arm1": [0.0, 1.0, -1.0, 0.0, 3.0],
            "harm_arm2": [1.0, 2.0, 0.5, 1.5, 8.0],
            "harm_arm3": [0.1, 0.2, 0.15, 0.25, 0.9],
        }
    )
    res_skip = evaluate.bootstrap_auroc_differences(one_class_pred, n_resamples=80, seed=42)
    assert res_skip["n_resamples_skipped_one_class"] > 0
    assert res_skip["n_resamples_used"] + res_skip["n_resamples_skipped_one_class"] == 80


def test_choose_sentences_rules() -> None:
    """choose_sentences selects S1/S2/S3/S4/S5/S6 and NO FIXED SENTENCE according to Section 8.9 rules."""
    from src import evaluate

    def _make_headline(arm2_catch: float, arm2_fa: int = 120) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "arm": "arm0",
                    "arm_name": config.ARM_NAMES["arm0"],
                    "cutoff": np.nan,
                    "n_harmful": 267,
                    "n_harmless": 899,
                    "harmful_caught": 0,
                    "harmful_missed": 267,
                    "false_alarms": 0,
                    "harmless_cleared": 899,
                    "catch_rate_percent": 0.0,
                    "false_alarm_rate_percent": 0.0,
                    "accuracy_percent": 77.1,
                },
                {
                    "arm": "arm1",
                    "arm_name": config.ARM_NAMES["arm1"],
                    "cutoff": 0.0,
                    "n_harmful": 267,
                    "n_harmless": 899,
                    "harmful_caught": 243,
                    "harmful_missed": 24,
                    "false_alarms": 400,
                    "harmless_cleared": 499,
                    "catch_rate_percent": 91.0,
                    "false_alarm_rate_percent": 44.5,
                    "accuracy_percent": 63.6,
                },
                {
                    "arm": "arm2",
                    "arm_name": config.ARM_NAMES["arm2"],
                    "cutoff": 5.0,
                    "n_harmful": 267,
                    "n_harmless": 899,
                    "harmful_caught": 243,
                    "harmful_missed": 24,
                    "false_alarms": arm2_fa,
                    "harmless_cleared": 899 - arm2_fa,
                    "catch_rate_percent": arm2_catch,
                    "false_alarm_rate_percent": round(100.0 * arm2_fa / 899, 1),
                    "accuracy_percent": 80.0,
                },
                {
                    "arm": "arm3",
                    "arm_name": config.ARM_NAMES["arm3"],
                    "cutoff": 0.5,
                    "n_harmful": 267,
                    "n_harmless": 899,
                    "harmful_caught": 243,
                    "harmful_missed": 24,
                    "false_alarms": 110,
                    "harmless_cleared": 789,
                    "catch_rate_percent": 91.0,
                    "false_alarm_rate_percent": 12.2,
                    "accuracy_percent": 88.5,
                },
            ]
        )

    s1_only = evaluate.choose_sentences(
        _make_headline(91.0, 120),
        {
            "esm_minus_blosum": {"point": 0.06, "low": 0.02, "high": 0.10},
            "arm3_minus_esm": {"point": 0.00, "low": -0.01, "high": 0.02},
        },
    )
    assert (
        "S1: A model that only read protein sequences ranked harmful BRCA1 mutations "
        "better than a classic substitution table."
    ) in s1_only
    assert not any(line.startswith("S4") for line in s1_only)
    assert (
        "S3: We could not detect a benefit from adding labels from the other end of the protein."
        in s1_only
    )
    assert not any(line.startswith("S5") for line in s1_only)
    assert (
        "S6: To catch 91.0% of harmful mutations, the model also flagged 120 of 899 harmless ones."
        in s1_only
    )

    s1_and_s4 = evaluate.choose_sentences(
        _make_headline(91.0, 120),
        {
            "esm_minus_blosum": {"point": 0.02, "low": 0.01, "high": 0.04},
            "arm3_minus_esm": {"point": 0.03, "low": 0.01, "high": 0.05},
        },
    )
    assert any(line.startswith("S1:") for line in s1_and_s4)
    assert any(line.startswith("S4:") for line in s1_and_s4)
    assert not any(line.startswith("S3") for line in s1_and_s4)
    assert any(line.startswith("S6:") for line in s1_and_s4)

    s2_case = evaluate.choose_sentences(
        _make_headline(91.0, 120),
        {
            "esm_minus_blosum": {"point": 0.01, "low": -0.02, "high": 0.05},
            "arm3_minus_esm": {"point": 0.00, "low": -0.01, "high": 0.01},
        },
    )
    assert "S2: No detectable difference at this sample size." in s2_case
    assert not any(line.startswith("S1") for line in s2_case)
    assert not any(line.startswith("S4") for line in s2_case)

    below_zero_and_s5 = evaluate.choose_sentences(
        _make_headline(83.0, 95),
        {
            "esm_minus_blosum": {"point": -0.04, "low": -0.08, "high": -0.01},
            "arm3_minus_esm": {"point": -0.01, "low": -0.03, "high": 0.00},
        },
    )
    assert (
        "NO FIXED SENTENCE: the interval is entirely below zero, so BLOSUM62 ranked "
        "better than ESM-2. Report this plainly and ask the humans for wording."
    ) in below_zero_and_s5
    assert any(
        line.startswith("S5")
        and "A cutoff tuned on one end of the protein did not carry over to the other: "
        "it caught 83.0% instead of 90%." in line
        for line in below_zero_and_s5
    )
    assert (
        "S6: To catch 83.0% of harmful mutations, the model also flagged 95 of 899 harmless ones."
        in below_zero_and_s5
    )


def test_verify_hash_on_temp_copy(tmp_path: Path) -> None:
    """verify_hash passes on a matching temporary CSV and raises RuntimeError on a 1-byte change."""
    from src import evaluate

    temp_csv = tmp_path / "esm_scores.csv"
    temp_csv.write_bytes(b"position,mutant,esm_score\n61,C61G,-9.123456\n")
    expected_sha = config.sha256_of_file(temp_csv)
    frozen = {"inputs": {"esm_scores_sha256": expected_sha}}

    evaluate.verify_hash(frozen, temp_csv)

    with open(temp_csv, "ab") as handle:
        handle.write(b" ")

    with pytest.raises(RuntimeError, match="HASH MISMATCH"):
        evaluate.verify_hash(frozen, temp_csv)


def test_may_run_refusal(tmp_path: Path) -> None:
    """may_run allows first run or refused-only logs and requires non-empty bugfix reason after any attempt."""
    from src import evaluate

    missing_log = tmp_path / "missing_log.txt"
    assert evaluate.prior_runs(missing_log) == []
    assert evaluate.may_run(missing_log, None)[0] is True

    refused_log = _write_log_lines(
        tmp_path, "refused.txt", [("refused", "hash check"), ("refused", "dirty")]
    )
    assert evaluate.may_run(refused_log, None)[0] is True

    started_log = _write_log_lines(tmp_path, "started.txt", [("started", "first run")])
    allowed_started, msg_started = evaluate.may_run(started_log, None)
    assert allowed_started is False
    assert "Refusing to run:" in msg_started

    success_log = _write_log_lines(
        tmp_path, "success.txt", [("started", "first run"), ("success", "first run")]
    )
    allowed_success, msg_success = evaluate.may_run(success_log, None)
    assert allowed_success is False
    assert "Refusing to run:" in msg_success

    failed_log = _write_log_lines(
        tmp_path, "failed.txt", [("started", "first run"), ("failed", "boom")]
    )
    allowed_failed, msg_failed = evaluate.may_run(failed_log, None)
    assert allowed_failed is False
    assert "Refusing to run:" in msg_failed

    assert evaluate.may_run(success_log, "fixed by-swap bug")[0] is True
    assert evaluate.may_run(success_log, "")[0] is False
    assert evaluate.may_run(success_log, "   ")[0] is False


def test_freeze_refuses_overwrite(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """freeze.main refuses to overwrite frozen.json without --force and refuses after test_attempted even with --force."""
    from src import freeze

    paths = _patch_stage4_paths(tmp_path, monkeypatch)
    paths["FROZEN_PATH"].write_text("{}\n", encoding="utf-8")

    code_existing = freeze.main([])
    captured_existing = capsys.readouterr()
    assert code_existing != 0
    assert "Refusing to overwrite outputs/frozen.json" in captured_existing.err

    for status in ("success", "started", "failed"):
        _write_log_lines(tmp_path, "temp_log.txt", [(status, "run attempt")])
        paths["RUN_LOG_PATH"].write_text(
            (tmp_path / "temp_log.txt").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        code_no_force = freeze.main([])
        code_force = freeze.main(["--force"])
        captured = capsys.readouterr()
        assert code_no_force != 0
        assert code_force != 0
        assert (
            "Refusing to freeze: the test run has started. Cutoffs may not change."
            in captured.err
        )


def test_locked_hashes_block_changed_inputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """locked_hashes returns (esm_sha256, frozen_sha256) and evaluate.main refuses when frozen.json hash changed."""
    import json

    from src import evaluate

    paths = _patch_stage4_paths(tmp_path, monkeypatch)

    assert evaluate.locked_hashes(paths["RUN_LOG_PATH"]) is None
    paths["RUN_LOG_PATH"].write_text(
        "2026-11-15T09:05:10-06:00 | run=1 | status=refused | kind=first | "
        "reason=check | commit=0123abc | esm_sha256=" + ("1" * 64) + " | frozen_sha256=" + ("2" * 64) + "\n",
        encoding="utf-8",
    )
    assert evaluate.locked_hashes(paths["RUN_LOG_PATH"]) is None

    paths["ESM_SCORES_PATH"].write_bytes(b"position,mutant,esm_score\n61,C61G,-9.0\n")
    esm_sha = config.sha256_of_file(paths["ESM_SCORES_PATH"])
    frozen_obj = _synthetic_frozen_dict(esm_sha256=esm_sha)
    paths["FROZEN_PATH"].write_text(json.dumps(frozen_obj, indent=2) + "\n", encoding="utf-8")
    logged_frozen_sha = "f" * 64

    paths["RUN_LOG_PATH"].write_text(
        f"2026-11-15T09:05:10-06:00 | run=1 | status=started | kind=first | "
        f"reason=first run | commit=0123abc | esm_sha256={esm_sha} | frozen_sha256={logged_frozen_sha}\n",
        encoding="utf-8",
    )
    assert evaluate.locked_hashes(paths["RUN_LOG_PATH"]) == (esm_sha, logged_frozen_sha)

    def _forbid_load_labels(*_args: object, **_kwargs: object) -> pd.DataFrame:
        raise AssertionError("load_labels must not be called when locked hashes mismatch")

    monkeypatch.setattr(load, "load_labels", _forbid_load_labels)

    exit_code = evaluate.main(["--bugfix", "x"])
    captured = capsys.readouterr()

    assert exit_code != 0
    assert (
        "Refusing to run: frozen.json or esm_scores.csv is not the file the first run used."
        in captured.err
    )
    log_lines = [
        line.strip()
        for line in paths["RUN_LOG_PATH"].read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(log_lines) == 2
    assert "status=refused" in log_lines[1]


_APP_MOD_CACHE: object | None = None


def _load_app_module():
    """Load app/app.py by path or skip if Stage 8 app data has not been exported yet."""
    global _APP_MOD_CACHE
    app_py = config.APP_DIR / "app.py"
    if not config.APP_DATA_PATH.exists() or not app_py.exists():
        pytest.skip("app data not exported yet (stage 8)")
    if _APP_MOD_CACHE is not None:
        return _APP_MOD_CACHE

    import importlib.util

    spec = importlib.util.spec_from_file_location("app_module", app_py)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _APP_MOD_CACHE = module
    return module


def test_app_parser_good() -> None:
    """One-letter, p.-prefixed, whitespace-padded, and three-letter C61G forms parse to ('C', 61, 'G')."""
    app_mod = _load_app_module()
    good_inputs = [
        "C61G",
        "c61g",
        "C61g",
        " C61G ",
        "p.C61G",
        "P.c61g",
        "Cys61Gly",
        "p.Cys61Gly",
        "CYS61GLY",
    ]
    for raw in good_inputs:
        assert app_mod.parse_input(raw) == ("C", 61, "G")


def test_app_refusals() -> None:
    """App inputs failing R1-R6 rules produce the exact Section 8.13 refusal strings."""
    app_mod = _load_app_module()

    r1_expected = (
        "Type a mutation such as C61G: the normal letter, the position, then the new letter."
    )
    for raw_r1 in ("", None):
        out_r1 = app_mod.check_mutation(raw_r1)
        assert out_r1.splitlines()[0] == r1_expected

    r2_cases = [
        ("hello", "hello"),
        ("X61G", "X61G"),
        ("c.181T>G", "c.181TG"),
        ("C61*", "C61"),
        ("Cys61Ter", "Cys61Ter"),
    ]
    for raw_r2, sanitized in r2_cases:
        expected_r2 = (
            f'I could not read "{sanitized}". Use the form C61G: the normal letter, '
            "the position, then the new letter. Three-letter names such as Cys61Gly also work."
        )
        out_r2 = app_mod.check_mutation(raw_r2)
        assert out_r2.splitlines()[0] == expected_r2

    out_r3 = app_mod.check_mutation("C61C")
    assert out_r3.splitlines()[0] == '"C61C" is not a change: both letters are C.'

    for raw_r4, pos_num in (("C0G", 0), ("C5000G", 5000)):
        out_r4 = app_mod.check_mutation(raw_r4)
        assert (
            out_r4.splitlines()[0]
            == f"Position {pos_num} is outside BRCA1, which has 1,863 positions."
        )

    out_r5 = app_mod.check_mutation("C999G")
    assert out_r5.splitlines()[0] == (
        "Position 999 is not one of the 326 lab-tested positions (1 to 101 and 1631 to 1855). "
        "This demo only answers where its answers can be checked against the lab, so it will not guess here."
    )

    out_r6 = app_mod.check_mutation("A61G")
    assert (
        out_r6.splitlines()[0]
        == "BRCA1 has C at position 61, not A. Did you mean C61G?"
    )


def test_app_notes_always_present() -> None:
    """Every valid and refused check_mutation call includes both fixed notes word for word."""
    app_mod = _load_app_module()
    for raw in ("C61G", "D2E", "C61A", "", "hello", "X61G", "C999G"):
        output = app_mod.check_mutation(raw)
        assert isinstance(output, str)
        assert config.DEMO_NOTE in output
        assert config.DISCLAIMER in output


def test_app_lab_display() -> None:
    """C61G displays Loss of function (harmful) and C61A displays not lab-tested."""
    app_mod = _load_app_module()
    out_c61g = app_mod.check_mutation("C61G")
    assert "Loss of function (harmful)" in out_c61g

    out_c61a = app_mod.check_mutation("C61A")
    assert "not lab-tested" in out_c61a


def test_app_end_region_locked(monkeypatch: pytest.MonkeyPatch) -> None:
    """When end_labels_included is False, S1631A returns refusal R8 with no verdict/score and both notes."""
    app_mod = _load_app_module()
    monkeypatch.setitem(app_mod.APP_CONFIG, "end_labels_included", False)

    output = app_mod.check_mutation("S1631A")
    expected_r8 = (
        "End-region answers (positions 1631 to 1855) are hidden until the one test run is done."
    )
    assert output.splitlines()[0] == expected_r8
    for forbidden in (
        "Looks harmful",
        "Can't tell",
        "Looks harmless",
        "**Model verdict:**",
        "**Model score:**",
    ):
        assert forbidden not in output
    assert config.DEMO_NOTE in output
    assert config.DISCLAIMER in output


def test_app_imports_no_model() -> None:
    """Importing app/app.py in a fresh subprocess does not import torch or transformers."""
    import subprocess
    import sys

    app_py = config.APP_DIR / "app.py"
    if not config.APP_DATA_PATH.exists() or not app_py.exists():
        pytest.skip("app data not exported yet (stage 8)")

    code = (
        "import importlib.util, sys; "
        f"spec = importlib.util.spec_from_file_location('app_module', {str(app_py)!r}); "
        "mod = importlib.util.module_from_spec(spec); "
        "spec.loader.exec_module(mod); "
        "assert 'torch' not in sys.modules, 'torch imported'; "
        "assert 'transformers' not in sys.modules, 'transformers imported'"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr


def test_app_verdict_end_region_track_record(monkeypatch: pytest.MonkeyPatch) -> None:
    """check_mutation includes the held-out end-region track record line per verdict when end labels are enabled."""
    app_mod = _load_app_module()

    out_harmful = app_mod.check_mutation("C61G")
    for snippet in (
        "End-region track record:",
        "196 of 288",
        "68.1%",
        "92 were FUNC harmless",
        "29 INT",
    ):
        assert snippet in out_harmful

    out_harmless = app_mod.check_mutation("D2E")
    for snippet in (
        "End-region track record:",
        "462 of 477",
        "96.9%",
        "15 were LOF harmful",
        "20 INT",
    ):
        assert snippet in out_harmless

    cutoff_a = float(app_mod.APP_CONFIG["cutoff_a"])
    raw_cutoff_b = app_mod.APP_CONFIG.get("cutoff_b")
    cutoff_b = None if raw_cutoff_b is None else float(raw_cutoff_b)
    cant_tell_mutant: str | None = None
    for row in app_mod.APP_DATA.to_dict("records"):
        verdict = app_mod.band_verdict(0.0 - float(row["esm_score"]), cutoff_a, cutoff_b)
        if verdict == "Can't tell":
            cant_tell_mutant = f"{row['wild_type']}{int(row['position'])}{row['new_letter']}"
            break
    assert cant_tell_mutant is not None

    out_cant_tell = app_mod.check_mutation(cant_tell_mutant)
    for snippet in (
        "End-region track record:",
        "56 LOF harmful, 345 FUNC harmless, and 47 INT",
        "448 total",
    ):
        assert snippet in out_cant_tell

    monkeypatch.setitem(app_mod.APP_CONFIG, "end_labels_included", False)
    out_locked = app_mod.check_mutation("C61G")
    assert "End-region track record:" not in out_locked
    assert config.DEMO_NOTE in out_locked
    assert config.DISCLAIMER in out_locked




# --- Post-test reference context (src/reference.py) ---


def test_reference_leaderboard_summary() -> None:
    """Ranks count strictly better models, so ties share a rank."""
    from src import reference

    board = pd.DataFrame(
        {
            "model": ["a", "b", "c", "d"],
            "spearman": [0.6, 0.5, 0.5, 0.2],
            "auc": [0.7, 0.9, 0.8, 0.6],
        }
    )
    summary = reference.leaderboard_summary(board, "c")
    assert summary["n_models"] == 4
    assert summary["spearman"]["our_model_rank"] == 2
    assert summary["spearman"]["best_model"] == "a"
    assert summary["spearman"]["best_value"] == 0.6
    assert summary["auc"]["our_model_rank"] == 2
    assert summary["auc"]["best_model"] == "b"
    with pytest.raises(ValueError):
        reference.leaderboard_summary(board, "missing")
    with pytest.raises(ValueError):
        reference.leaderboard_summary(board.rename(columns={"auc": "x"}), "c")


def test_reference_parse_alphamissense() -> None:
    """The third-party table must agree with our wild type and keep scores in 0 to 1."""
    from src import reference

    good = pd.DataFrame(
        {
            "protein_variant": ["M1A", "D2E"],
            "am_pathogenicity": [0.9, 0.1],
            "am_class": ["LPath", "LBen"],
        }
    )
    parsed = reference.parse_alphamissense(good, "MDLS")
    assert list(parsed.columns) == ["mutant", "harm_am"]
    assert parsed["harm_am"].tolist() == [0.9, 0.1]

    wrong_letter = good.assign(protein_variant=["M1A", "L2E"])
    with pytest.raises(ValueError, match="1 rows disagree"):
        reference.parse_alphamissense(wrong_letter, "MDLS")
    with pytest.raises(ValueError, match="outside 0 to 1"):
        reference.parse_alphamissense(good.assign(am_pathogenicity=[1.5, 0.1]), "MDLS")
    with pytest.raises(ValueError, match="unexpected columns"):
        reference.parse_alphamissense(good.drop(columns=["am_class"]), "MDLS")


def test_reference_bootstrap_resamples_positions() -> None:
    """The paired difference is seeded, signed a minus b, and resamples whole positions."""
    from src import reference

    positions = np.repeat(np.arange(10), 3)
    is_harmful = np.tile([1, 0, 0], 10)
    harm_a = is_harmful.astype(float)
    harm_b = 1.0 - harm_a

    first = reference.paired_position_bootstrap(positions, is_harmful, harm_a, harm_b, 50, 42)
    second = reference.paired_position_bootstrap(positions, is_harmful, harm_a, harm_b, 50, 42)
    assert first == second
    assert first["point"] == 1.0
    assert first["low"] == 1.0 and first["high"] == 1.0
    assert first["n_positions"] == 10
    assert first["n_resamples_used"] == 50

    same = reference.paired_position_bootstrap(positions, is_harmful, harm_a, harm_a, 50, 42)
    assert same["point"] == 0.0 and same["low"] == 0.0 and same["high"] == 0.0

    # One position holds the only harmful row, so some resamples have one class and are skipped.
    lone = np.zeros(30, dtype=int)
    lone[0] = 1
    skipped = reference.paired_position_bootstrap(positions, lone, harm_a, harm_b, 200, 42)
    assert skipped["n_resamples_skipped_one_class"] > 0
    assert skipped["n_resamples_used"] + skipped["n_resamples_skipped_one_class"] == 200


def test_reference_refuses_before_test(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Before the one test run succeeds, the reference module writes nothing."""
    from src import reference

    monkeypatch.setattr(config, "RUN_LOG_PATH", tmp_path / "run_log.txt")
    monkeypatch.setattr(
        config, "REFERENCE_ALPHAMISSENSE_RESULT_PATH", tmp_path / "reference_alphamissense.json"
    )
    monkeypatch.setattr(
        config, "REFERENCE_PROTEINGYM_RESULT_PATH", tmp_path / "reference_proteingym.json"
    )

    def _must_not_load(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("labels must not be loaded before the test run")

    monkeypatch.setattr(load, "load_labels", _must_not_load)
    assert reference.main([]) == 1
    assert list(tmp_path.iterdir()) == []


# --- Exploratory arm 4 (src/arm4.py) ---


def _arm4_synthetic_region(n_lof: int, n_int: int, n_func: int, first_position: int) -> pd.DataFrame:
    """Made-up merged table: harmful rows get more negative scores, three rows per position."""
    rng = np.random.default_rng(7)
    classes = ["LOF"] * n_lof + ["INT"] * n_int + ["FUNC"] * n_func
    rng.shuffle(classes)
    n = len(classes)
    shift = np.array([{"LOF": -6.0, "INT": -3.0, "FUNC": -1.0}[c] for c in classes])
    return pd.DataFrame(
        {
            "mutant": [f"A{i + 1}C" for i in range(n)],
            "position": first_position + np.arange(n) // 3,
            "findlay_class": classes,
            "esm_score": shift + rng.normal(0.0, 2.0, n),
            "esmc_score": shift + rng.normal(0.0, 1.0, n),
        }
    )


def test_arm4_pure_rules(tmp_path: Path) -> None:
    """Candidate harm scores, the tie rule, the fixed sentences and the log reader."""
    from src import arm4

    esm2 = np.array([-4.0, 1.0])
    esmc = np.array([-2.0, -1.0])
    assert arm4.candidate_harm(esm2, esmc, "esmc").tolist() == [2.0, 1.0]
    assert arm4.candidate_harm(esm2, esmc, "esm2_esmc_mean").tolist() == [3.0, 0.0]
    with pytest.raises(ValueError):
        arm4.candidate_harm(esm2, esmc, "other")

    assert arm4.choose_candidate({"esmc": 0.80, "esm2_esmc_mean": 0.85}) == "esm2_esmc_mean"
    assert arm4.choose_candidate({"esmc": 0.85, "esm2_esmc_mean": 0.80}) == "esmc"
    assert arm4.choose_candidate({"esmc": 0.85, "esm2_esmc_mean": 0.85}) == "esmc"

    assert arm4.choose_sentence({"low": 0.01, "high": 0.05}) == arm4.SENTENCE_BETTER
    assert arm4.choose_sentence({"low": -0.01, "high": 0.05}) == arm4.SENTENCE_NO_DIFFERENCE
    assert arm4.choose_sentence({"low": 0.0, "high": 0.05}) == arm4.SENTENCE_NO_DIFFERENCE
    assert arm4.choose_sentence({"low": -0.05, "high": -0.01}) == arm4.SENTENCE_WORSE

    assert arm4.log_statuses(tmp_path / "missing.txt") == []
    log = tmp_path / "log.txt"
    log.write_text(
        "t | run=1 | status=started | reason=status=success in text\n"
        "t | run=1 | status=failed | reason=x\n",
        encoding="utf-8",
    )
    assert arm4.log_statuses(log) == ["started", "failed"]


def test_arm4_freeze_then_one_evaluate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """On made-up data: freeze uses the start region only, evaluate runs once and then refuses."""
    from src import arm4

    for name in (
        "ARM4_SCORES_PATH", "ARM4_FROZEN_PATH", "ARM4_RUN_LOG_PATH", "ARM4_RESULT_PATH",
        "ARM4_PARAGRAPH_PATH", "ESM_SCORES_PATH", "FROZEN_PATH",
    ):
        monkeypatch.setattr(config, name, tmp_path / getattr(config, name).name)
    for name in ("ARM4_SCORES_PATH", "ESM_SCORES_PATH", "FROZEN_PATH"):
        getattr(config, name).write_text(name, encoding="utf-8")

    tables = {
        "start": _arm4_synthetic_region(131, 59, 385, 1),
        "end": _arm4_synthetic_region(267, 96, 899, 1631),
    }
    regions_asked: list[tuple[str, str]] = []

    def _fake_labels(region: str, purpose: str) -> pd.DataFrame:
        regions_asked.append((region, purpose))
        return tables[region]

    monkeypatch.setattr(load, "load_labels", _fake_labels)
    monkeypatch.setattr(load, "test_run_done", lambda *_a, **_k: True)
    monkeypatch.setattr(arm4, "_merged_scores", lambda region, labels: labels.copy())
    monkeypatch.setattr(
        arm4, "read_arm4_scores", lambda: pd.DataFrame({"model_id": [config.ARM4_MODEL_ID]})
    )

    assert arm4.main(["evaluate"]) == 1  # nothing frozen yet
    assert regions_asked == []

    assert arm4.main(["freeze"]) == 0
    assert regions_asked == [("start", "tuning")]
    frozen = json.loads(config.ARM4_FROZEN_PATH.read_text(encoding="utf-8"))
    assert frozen["chosen"] in config.ARM4_CANDIDATES
    assert frozen["tuning_caught"] >= 118
    assert set(frozen["candidates_start_auroc"]) == set(config.ARM4_CANDIDATES)
    assert arm4.main(["freeze"]) == 1  # never overwritten

    assert arm4.main(["writeup"]) == 1  # no successful run yet
    assert arm4.main(["evaluate"]) == 0
    assert regions_asked[-1] == ("end", "post_test")
    assert arm4.log_statuses(config.ARM4_RUN_LOG_PATH) == ["started", "success"]
    result = json.loads(config.ARM4_RESULT_PATH.read_text(encoding="utf-8"))
    end = result["end"]
    assert end["n_harmful"] == 267 and end["n_harmless"] == 899
    assert end["harmful_caught"] + end["harmful_missed"] == 267
    assert result["int_rows_end"]["n_int"] == 96
    assert result["sentence"] in (
        arm4.SENTENCE_BETTER, arm4.SENTENCE_NO_DIFFERENCE, arm4.SENTENCE_WORSE
    )
    paragraph = config.ARM4_PARAGRAPH_PATH.read_text(encoding="utf-8")
    assert "exploratory" in paragraph and result["sentence"] in paragraph
    assert f"{end['harmful_caught']} of 267" in paragraph

    n_asked = len(regions_asked)
    assert arm4.main(["evaluate"]) == 1  # one run only
    assert len(regions_asked) == n_asked
    assert arm4.log_statuses(config.ARM4_RUN_LOG_PATH) == ["started", "success"]


def test_arm4_evaluate_refuses_changed_scores(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A score file that changed after the arm 4 freeze stops the run before any label is read."""
    from src import arm4

    for name in ("ARM4_SCORES_PATH", "ARM4_FROZEN_PATH", "ARM4_RUN_LOG_PATH", "ESM_SCORES_PATH"):
        monkeypatch.setattr(config, name, tmp_path / getattr(config, name).name)
    config.ARM4_SCORES_PATH.write_text("scores", encoding="utf-8")
    config.ESM_SCORES_PATH.write_text("esm", encoding="utf-8")
    frozen = {
        "inputs": {
            "arm4_scores_sha256": config.sha256_of_file(config.ARM4_SCORES_PATH),
            "esm_scores_sha256": config.sha256_of_file(config.ESM_SCORES_PATH),
        }
    }
    config.ARM4_FROZEN_PATH.write_text(json.dumps(frozen), encoding="utf-8")
    config.ARM4_SCORES_PATH.write_text("scores changed", encoding="utf-8")

    def _must_not_load(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("labels must not be loaded")

    monkeypatch.setattr(load, "load_labels", _must_not_load)
    monkeypatch.setattr(load, "test_run_done", lambda *_a, **_k: True)
    assert arm4.main(["evaluate"]) == 1
    assert not config.ARM4_RUN_LOG_PATH.exists()
