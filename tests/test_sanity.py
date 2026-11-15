"""Sanity tests against the verified real input files and score cache."""

import numpy as np
import pytest
import scipy.stats

from src import config, esm_score, load


def test_sanity_1_counts() -> None:
    """Input checks pass and counts, row count, uniqueness, and position count match spec."""
    load.check_inputs()
    assert load.class_counts() == {
        "start": {"LOF": 131, "INT": 59, "FUNC": 385},
        "end": {"LOF": 267, "INT": 96, "FUNC": 899},
    }
    mutations = load.load_mutations()
    assert len(mutations) == 1837
    assert mutations["mutant"].is_unique
    assert mutations["position"].nunique() == 326


def test_sanity_2_wild_type() -> None:
    """Wild-type FASTA has expected length/prefix and agrees with every mutation row."""
    seq = load.load_wild_type()
    assert len(seq) == 1863
    assert seq.startswith("MDLSALRVEEVQNVINAMQK")
    mutations = load.load_mutations()
    for row in mutations.to_dict("records"):
        assert row["wild_type"] == seq[row["position"] - 1]


def test_sanity_6_arm0() -> None:
    """Arm 0 end-region baseline counts and accuracy match spec using counts only."""
    counts = load.class_counts()
    end_counts = counts["end"]
    n_harmful = end_counts["LOF"]
    n_harmless = end_counts["FUNC"]
    assert n_harmful == 267
    assert n_harmless == 899

    total_no_int = n_harmful + n_harmless
    assert total_no_int == 1166

    caught = 0
    correct = n_harmless
    assert caught == 0
    assert correct == 899
    assert round(100 * correct / total_no_int, 1) == 77.1


def _require_score_cache() -> None:
    """Skip Stage 3 cached-score sanity checks until outputs/esm_scores.csv exists."""
    if not config.ESM_SCORES_PATH.exists():
        pytest.skip("scores not computed yet (stage 3)")


def test_sanity_3a_letter_test_cached() -> None:
    """Cached ESM score table matches wild-type residues, window names, and token indices."""
    _require_score_cache()
    df = load.load_scores()
    assert len(df) == config.N_SCORE_ROWS
    assert bool((df["token_at_index"] == df["wild_type"]).all())
    positions = [int(pos) for pos in df["position"]]
    assert df["window"].tolist() == [esm_score.window_of(pos) for pos in positions]
    assert df["token_index"].tolist() == [esm_score.token_index_of(pos) for pos in positions]


def test_sanity_3b_letter_test_live() -> None:
    """Live tokenizer letter test passes across all 326 scored positions."""
    try:
        import transformers

        tokenizer = transformers.AutoTokenizer.from_pretrained(config.MODEL_ID)
    except Exception as exc:
        pytest.skip(str(exc) or type(exc).__name__)

    res = esm_score.run_letter_test(
        tokenizer, load.load_wild_type(), esm_score.scored_positions()
    )
    assert len(res) == config.N_POSITIONS
    assert bool(res["ok"].all())


def test_sanity_4_known_harmful_start() -> None:
    """Known harmful start-region variants C61G and C64G score negative and in lowest quartile."""
    _require_score_cache()
    scores = load.lab_tested_scores("start")
    assert len(scores) == 575

    by_mutant = scores.set_index("mutant", verify_integrity=True)["esm_score"]
    q25 = float(np.percentile(scores["esm_score"].to_numpy(), 25))
    for mutant in ("C61G", "C64G"):
        val = float(by_mutant.loc[mutant])
        assert val < 0
        assert val <= q25


def test_sanity_5_start_correlation() -> None:
    """Start-region ESM scores correlate positively above tripwire with lab function scores."""
    _require_score_cache()
    labels = load.load_labels("start", "tuning")
    scores = load.lab_tested_scores("start")
    merged = scores.merge(labels, on="mutant", validate="one_to_one")
    assert len(merged) == 575

    res = scipy.stats.spearmanr(merged["esm_score"], merged["function_score"])
    assert res.statistic > config.SPEARMAN_MIN_START


def test_sanity_7_after_test() -> None:
    """Post-test sanity check verifies negative M1775R/A1708E scores and positive overall Spearman."""
    import json

    if not load.test_run_done():
        pytest.skip("runs only after the one test run")

    sanity7_path = config.RESULTS_DIR / "sanity7.json"
    assert sanity7_path.exists()
    with open(sanity7_path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    assert data["M1775R_esm_score"] < 0
    assert data["A1708E_esm_score"] < 0
    assert data["spearman_all_1837"] > 0

