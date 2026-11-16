"""Exploratory arm 4 for FoldFlux: a stronger zero-shot model, added after the test run.

Arm 4 scores every swap with ESM-C 600M, using the same hide-a-letter code as arm 2.
It was added after the one test run, with knowledge of arm 2's end-region result, so it
is exploratory and never replaces the headline. It keeps the project's discipline all the
same: the choice between two candidates (ESM-C alone, or the mean of the ESM-2 and ESM-C
scores) and the cutoff are made on the start region only, written to its own frozen file,
and then applied once to the end region with its own run log. Four steps, in order:
`score` (needs a GPU), `freeze`, `evaluate` (once), `writeup`.
"""

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from src import config, cutoffs, esm_score, load, reference

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

# Fixed in advance, like the sentences in BUILD_SPEC.md section 3.4. The interval is the
# 95% position bootstrap of arm 4's AUROC minus arm 2's on the end region.
SENTENCE_BETTER: str = "On the held-out region the added arm ranked harmful mutations better than `arm2`."
SENTENCE_NO_DIFFERENCE: str = "No detectable difference from `arm2` at this sample size."
SENTENCE_WORSE: str = "On the held-out region the added arm ranked harmful mutations worse than `arm2`."

_PARAGRAPH_MARKER: str = "<!-- arm4.json -->"
_LIMITS_ANCHOR: str = "**Limits:** (1) Blind to splicing"


def read_arm4_scores() -> pd.DataFrame:
    """Read and validate the arm 4 score file (same layout and checks as esm_scores.csv)."""
    if not config.ARM4_SCORES_PATH.is_file():
        raise FileNotFoundError(
            "outputs/arm4_esmc_scores.csv not found. Run: python -m src.arm4 score"
        )
    df = pd.read_csv(config.ARM4_SCORES_PATH, dtype=_SCORE_DTYPES, keep_default_na=False)
    load.validate_scores(df)
    return df


def candidate_harm(esm2_score: np.ndarray, esmc_score: np.ndarray, candidate: str) -> np.ndarray:
    """Harm score of one candidate; higher means more likely harmful, as for every arm."""
    esm2 = np.asarray(esm2_score, dtype=float)
    esmc = np.asarray(esmc_score, dtype=float)
    if candidate == "esmc":
        return 0.0 - esmc
    if candidate == "esm2_esmc_mean":
        return 0.0 - (esm2 + esmc) / 2.0
    raise ValueError(f"Unknown arm 4 candidate '{candidate}'.")


def choose_candidate(start_aurocs: dict[str, float]) -> str:
    """Pick the candidate with the higher start-region AUROC; a tie goes to the simpler one."""
    first, second = config.ARM4_CANDIDATES
    return second if start_aurocs[second] > start_aurocs[first] else first


def choose_sentence(difference: dict) -> str:
    """Map the bootstrap interval of arm 4 minus arm 2 to one fixed sentence."""
    if difference["low"] > 0:
        return SENTENCE_BETTER
    if difference["high"] < 0:
        return SENTENCE_WORSE
    return SENTENCE_NO_DIFFERENCE


def log_statuses(run_log_path: Path) -> list[str]:
    """Return every status value in the arm 4 run log, in order; a missing file gives none."""
    if not run_log_path.is_file():
        return []
    statuses: list[str] = []
    for line in run_log_path.read_text(encoding="utf-8").splitlines():
        for part in (piece.strip() for piece in line.split("|")):
            if part.startswith("status="):
                statuses.append(part[len("status="):])
    return statuses


def _append_log(status: str, reason: str) -> None:
    """Append one line to the arm 4 run log. Only this module writes that file."""
    stamp = datetime.now().astimezone().isoformat(timespec="seconds")
    clean = reason.replace("|", " ").replace("\n", " ")
    line = (
        f"{stamp} | run=1 | status={status} | reason={clean} | "
        f"arm4_scores_sha256={config.sha256_of_file(config.ARM4_SCORES_PATH)} | "
        f"frozen_arm4_sha256={config.sha256_of_file(config.ARM4_FROZEN_PATH)}\n"
    )
    config.ARM4_RUN_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(config.ARM4_RUN_LOG_PATH, "a", newline="\n", encoding="utf-8") as handle:
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


def _merged_scores(region: str, labels: pd.DataFrame) -> pd.DataFrame:
    """Join one region's labels with the arm 2 and arm 4 scores, one row per lab-tested swap."""
    esm2 = load.lab_tested_scores(region)[["mutant", "esm_score"]]
    esmc = read_arm4_scores()[["mutant", "esm_score"]].rename(columns={"esm_score": "esmc_score"})
    merged = labels.merge(esm2, on="mutant", validate="one_to_one")
    merged = merged.merge(esmc, on="mutant", validate="one_to_one")
    if len(merged) != len(labels):
        raise ValueError("Arm 4 check failed: a lab-tested mutation has no score.")
    return merged


def _binary_rows(table: pd.DataFrame) -> pd.DataFrame:
    """Keep LOF and FUNC rows and add is_harmful; INT rows are left out, as everywhere."""
    keep = table["findlay_class"].isin([config.HARMFUL_CLASS, config.HARMLESS_CLASS])
    out = table.loc[keep].copy()
    out["is_harmful"] = (out["findlay_class"] == config.HARMFUL_CLASS).astype("int64")
    return out


def run_score(allow_cpu: bool, batch_size: int) -> int:
    """Step 1 (GPU): score all 326 positions with ESM-C 600M. Uses no labels."""
    if config.ARM4_SCORES_PATH.exists():
        print("Refusing to score: outputs/arm4_esmc_scores.csv already exists.", file=sys.stderr)
        return 1
    try:
        import torch
        import transformers
        from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer
    except ImportError as exc:
        print(
            f"Arm 4 needs the ESM-C package ({exc}). Install it with: pip install esm",
            file=sys.stderr,
        )
        return 1

    device = esm_score.choose_device(allow_cpu)
    tokenizer = EsmcTokenizer()
    model = EsmcForMaskedLM.from_pretrained(
        config.ARM4_MODEL_ID, device=device, dtype=torch.float32
    )
    model.eval()

    seq = load.load_wild_type()
    # score_positions runs the letter test first, so a tokenizer that splits or pads
    # differently from ESM-2's stops here instead of producing shifted scores.
    scores = esm_score.score_positions(
        model, tokenizer, seq, esm_score.scored_positions(), device, batch_size
    )
    tested = set(load.load_mutations()["mutant"])
    scores["is_lab_tested"] = scores["mutant"].isin(tested).astype("int64")
    scores["model_id"] = config.ARM4_MODEL_ID
    scores["transformers_version"] = str(transformers.__version__)
    scores["precision"] = "float32"
    scores = scores[list(_SCORE_DTYPES)]

    tmp_path = config.ARM4_SCORES_PATH.parent / f"{config.ARM4_SCORES_PATH.name}.tmp"
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        scores.to_csv(handle, index=False, lineterminator="\n", float_format="%.6f")
    os.replace(tmp_path, config.ARM4_SCORES_PATH)

    checked = read_arm4_scores()
    c61g = float(checked.loc[checked["mutant"] == "C61G", "esm_score"].iloc[0])
    print(f"Wrote outputs/arm4_esmc_scores.csv: {len(checked)} rows, model {config.ARM4_MODEL_ID}.")
    print(f"C61G (start region): esmc_score = {c61g:.6f} (must be negative)")
    print(f"sha256 = {config.sha256_of_file(config.ARM4_SCORES_PATH)}")
    if c61g >= 0:
        print("STOP: C61G is not negative. Do not continue; report this.", file=sys.stderr)
        return 1
    return 0


def run_freeze() -> int:
    """Step 2: choose the candidate and its cutoff on the start region only, then freeze."""
    if config.ARM4_FROZEN_PATH.exists():
        print("Refusing to overwrite outputs/frozen_arm4.json.", file=sys.stderr)
        return 1
    if log_statuses(config.ARM4_RUN_LOG_PATH):
        print("Refusing to freeze: the arm 4 run has started.", file=sys.stderr)
        return 1
    if not load.test_run_done():
        print("Refusing to freeze: arm 4 is an addition after the one test run.", file=sys.stderr)
        return 1

    start = _binary_rows(_merged_scores("start", load.load_labels("start", "tuning")))
    if len(start) != 516:
        raise ValueError("Arm 4 check failed: expected 516 start-region LOF and FUNC rows.")
    y = start["is_harmful"].to_numpy()
    esm2 = start["esm_score"].to_numpy()
    esmc = start["esmc_score"].to_numpy()

    aurocs = {
        name: float(roc_auc_score(y, candidate_harm(esm2, esmc, name)))
        for name in config.ARM4_CANDIDATES
    }
    chosen = choose_candidate(aurocs)
    harm = candidate_harm(esm2, esmc, chosen)
    cutoff = cutoffs.choose_catch_cutoff(harm, y, config.CATCH_TARGET_PERCENT)
    counts = cutoffs.count_outcomes(harm, y, cutoff)

    frozen = {
        "what": (
            "EXPLORATORY arm 4, added after the one test run. Candidate and cutoff were "
            "chosen on the start region only. Nothing here changes outputs/frozen.json."
        ),
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "seed": config.SEED,
        "model_id": str(read_arm4_scores()["model_id"].iloc[0]),
        "tuning_region": "start",
        "test_region": "end",
        "selection_rule": "higher start-region AUROC (LOF vs FUNC, 516 rows); tie goes to esmc",
        "candidates_start_auroc": aurocs,
        "arm2_start_auroc_in_sample": float(roc_auc_score(y, 0.0 - esm2)),
        "chosen": chosen,
        "chosen_name": config.ARM4_CANDIDATE_NAMES[chosen],
        "cutoff_rule": "highest cutoff catching at least 90% of start-region harmful rows",
        "cutoff": float(cutoff),
        "tuning_caught": int(counts["harmful_caught"]),
        "tuning_false_alarms": int(counts["false_alarms"]),
        "inputs": {
            "arm4_scores_sha256": config.sha256_of_file(config.ARM4_SCORES_PATH),
            "esm_scores_sha256": config.sha256_of_file(config.ESM_SCORES_PATH),
            "frozen_sha256": config.sha256_of_file(config.FROZEN_PATH),
        },
    }
    _write_json(frozen, config.ARM4_FROZEN_PATH)
    print(json.dumps(frozen, indent=2))
    print("Start-region numbers are in-sample; they are not results.")
    return 0


def paragraph_text(frozen: dict, result: dict) -> str:
    """The one write-up paragraph for arm 4, built only from the two arm 4 files."""
    other = [name for name in config.ARM4_CANDIDATES if name != frozen["chosen"]][0]
    aurocs = frozen["candidates_start_auroc"]
    end = result["end"]
    diff = result["arm4_minus_arm2"]
    return (
        "**Extra arm, added after the test run (exploratory):** We then added `arm4`, "
        f"{frozen['chosen_name']}, scored with the same hide-a-letter method. It was chosen over "
        f"{config.ARM4_CANDIDATE_NAMES[other]} by start-region AUROC alone "
        f"({aurocs[frozen['chosen']]:.4f} vs {aurocs[other]:.4f}), with its cutoff set by the "
        "same 90% rule on `1–101` and frozen in `outputs/frozen_arm4.json` before one run on the "
        "end region (`outputs/results/arm4_run_log.txt`). It caught "
        f"{end['harmful_caught']} of {end['n_harmful']} harmful mutations "
        f"({end['harmful_missed']} missed) with {end['false_alarms']} of {end['n_harmless']} "
        f"false alarms; AUROC {end['auroc']:.4f} vs {end['auroc_arm2_same_rows']:.4f} for `arm2` "
        f"(difference `{diff['point']:+.3f}`, 95% interval `{diff['low']:+.3f}` to "
        f"`{diff['high']:+.3f}`). {result['sentence']} We picked this arm after seeing `arm2`'s "
        "end-region results and ProteinGym's leaderboard, so it is exploratory and the headline "
        f"table is unchanged. {_PARAGRAPH_MARKER}\n"
    )


def run_evaluate() -> int:
    """Step 3: apply the frozen arm 4 to the end region, once."""
    if not load.test_run_done():
        print("Refusing to run: the one test run is not done.", file=sys.stderr)
        return 1
    if not config.ARM4_FROZEN_PATH.is_file():
        print("Refusing to run: outputs/frozen_arm4.json not found. Run freeze first.", file=sys.stderr)
        return 1
    if log_statuses(config.ARM4_RUN_LOG_PATH):
        print(
            "Refusing to run: outputs/results/arm4_run_log.txt shows arm 4 has already run. "
            "Ask a human.",
            file=sys.stderr,
        )
        return 1
    with open(config.ARM4_FROZEN_PATH, "r", encoding="utf-8") as handle:
        frozen = json.load(handle)
    if (
        config.sha256_of_file(config.ARM4_SCORES_PATH) != frozen["inputs"]["arm4_scores_sha256"]
        or config.sha256_of_file(config.ESM_SCORES_PATH) != frozen["inputs"]["esm_scores_sha256"]
    ):
        print("HASH MISMATCH: a score file changed since frozen_arm4.json was written.", file=sys.stderr)
        return 1

    _append_log("started", "arm 4 run")
    try:
        end_all = _merged_scores("end", load.load_labels("end", "post_test"))
        end_all["harm_arm4"] = candidate_harm(
            end_all["esm_score"].to_numpy(), end_all["esmc_score"].to_numpy(), frozen["chosen"]
        )
        end_all["harm_arm2"] = 0.0 - end_all["esm_score"]
        end = _binary_rows(end_all)
        if len(end_all) != 1262 or len(end) != 1166:
            raise ValueError("Arm 4 check failed: unexpected end-region row counts.")
        y = end["is_harmful"].to_numpy()
        harm4 = end["harm_arm4"].to_numpy()
        harm2 = end["harm_arm2"].to_numpy()

        counts = cutoffs.count_outcomes(harm4, y, frozen["cutoff"])
        difference = reference.paired_position_bootstrap(
            end["position"].to_numpy(), y, harm4, harm2, config.N_BOOTSTRAP, config.SEED
        )
        int_rows = end_all.loc[end_all["findlay_class"] == config.EXCLUDED_CLASS]
        result = {
            "what": (
                "EXPLORATORY arm 4 on the held-out end region, run once after the headline "
                "test run. Same 1,166 rows as the headline. Not a replacement for arm 2."
            ),
            "chosen": frozen["chosen"],
            "chosen_name": frozen["chosen_name"],
            "model_id": frozen["model_id"],
            "cutoff": float(frozen["cutoff"]),
            "end": {
                **{k: int(v) for k, v in counts.items()},
                **cutoffs.rates(counts),
                "auroc": round(float(roc_auc_score(y, harm4)), 4),
                "auroc_arm2_same_rows": round(float(roc_auc_score(y, harm2)), 4),
            },
            "arm4_minus_arm2": difference,
            "int_rows_end": {
                "n_int": int(len(int_rows)),
                "flagged": int(cutoffs.flag(int_rows["harm_arm4"].to_numpy(), frozen["cutoff"]).sum()),
            },
            "sentence": choose_sentence(difference),
        }
        _write_json(result, config.ARM4_RESULT_PATH)
        _write_text(paragraph_text(frozen, result), config.ARM4_PARAGRAPH_PATH)
    except BaseException as exc:
        _append_log("failed", str(exc) or type(exc).__name__)
        raise
    _append_log("success", "arm 4 run")
    print(json.dumps(result, indent=2))
    print(config.ARM4_PARAGRAPH_PATH.read_text(encoding="utf-8"))
    return 0


def _insert_once(path: Path, anchor: str, addition: str, marker: str, after: bool) -> str:
    """Insert text next to an anchor unless the marker is already in the file."""
    text = path.read_text(encoding="utf-8")
    if marker in text:
        return f"{path.name}: already done"
    if text.count(anchor) != 1:
        raise ValueError(f"Write-up check failed: anchor not found exactly once in {path.name}.")
    new = anchor + addition if after else addition + anchor
    _write_text(text.replace(anchor, new), path)
    return f"{path.name}: updated"


def run_writeup() -> int:
    """Step 4: put the arm 4 paragraph and credits into the write-up files, once."""
    if "success" not in log_statuses(config.ARM4_RUN_LOG_PATH):
        print("Refusing to edit the write-up: arm 4 has no successful run.", file=sys.stderr)
        return 1
    paragraph = config.ARM4_PARAGRAPH_PATH.read_text(encoding="utf-8")
    writeup = config.ROOT / "writeup"
    source = (
        " ESM-C 600M `biohub/ESMC-600M` (ESM Cambrian, EvolutionaryScale/Biohub; "
        "https://huggingface.co/biohub/ESMC-600M, https://github.com/evolutionaryscale/esm), "
        "used only for the exploratory `arm4`;"
    )
    credit = (
        "\n- ESM-C 600M (ESM Cambrian; EvolutionaryScale, now Biohub), model `biohub/ESMC-600M`, "
        "run through the `esm` Python package: https://huggingface.co/biohub/ESMC-600M and "
        "https://github.com/evolutionaryscale/esm — used only for the exploratory `arm4`, "
        "added after the one test run."
    )
    decision = (
        "\n## Exploratory Arm 4 (Human Decision, After the Test Run)\n\n"
        "- **Arm 4 (`src/arm4.py`, written by Claude Code):** At a human's request, added an "
        "exploratory arm scored with ESM-C 600M (`biohub/ESMC-600M`) through the same "
        "hide-a-letter code as `arm2`. The candidate (ESM-C alone, or the mean of the ESM-2 and "
        "ESM-C scores) and the 90%-catch cutoff were chosen on the start region only, frozen in "
        "`outputs/frozen_arm4.json`, and applied once to the end region "
        "(`outputs/results/arm4_run_log.txt`, `outputs/results/arm4.json`). It was picked after "
        "seeing `arm2`'s end-region results and ProteinGym's leaderboard, so it is reported as "
        "exploratory; `outputs/frozen.json`, `outputs/esm_scores.csv`, the original run log and "
        "the headline table are unchanged.\n"
    )
    anchor_models = "https://huggingface.co/facebook/esm2_t33_650M_UR50D, https://github.com/facebookresearch/esm);"
    for name in ("README.md", "writeup/writeup.md"):
        path = config.ROOT / name
        print(_insert_once(path, _LIMITS_ANCHOR, paragraph + "\n", _PARAGRAPH_MARKER, after=False))
        print(_insert_once(path, anchor_models, source, "biohub/ESMC-600M", after=True))
    print(
        _insert_once(
            writeup / "credits.md",
            "https://huggingface.co/facebook/esm2_t33_650M_UR50D and https://github.com/facebookresearch/esm",
            credit,
            "biohub/ESMC-600M",
            after=True,
        )
    )
    log_path = writeup / "decision_log.md"
    log_text = log_path.read_text(encoding="utf-8")
    if "Exploratory Arm 4" not in log_text:
        _write_text(log_text.rstrip("\n") + "\n" + decision, log_path)
    if (config.ROOT / "README.md").read_bytes() != (writeup / "writeup.md").read_bytes():
        print("STOP: README.md and writeup/writeup.md differ. Fix by hand.", file=sys.stderr)
        return 1
    print("Write-up updated. README.md and writeup/writeup.md are identical.")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Command line for the four arm 4 steps: score, freeze, evaluate, writeup."""
    parser = argparse.ArgumentParser(description="Exploratory arm 4 (ESM-C 600M).")
    parser.add_argument("step", choices=["score", "freeze", "evaluate", "writeup"])
    parser.add_argument("--allow-cpu", action="store_true", help="Score on CPU (slow).")
    parser.add_argument("--batch-size", type=int, default=config.ESM_BATCH_SIZE)
    args = parser.parse_args(argv)
    config.set_seeds()
    try:
        if args.step == "score":
            return run_score(args.allow_cpu, args.batch_size)
        if args.step == "freeze":
            return run_freeze()
        if args.step == "evaluate":
            return run_evaluate()
        return run_writeup()
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
