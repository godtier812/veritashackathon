"""Exploratory arm 6 for FoldFlux: an ensemble built from ideas the leading models use.

On ProteinGym's leaderboard for this lab table, 15 zero-shot models rank above ESM-2 650M.
Most of them add one of three things to a single protein language model:

- Several language models together (ESM-C 600M ranks 5th; ProtSSN and ESM-1v are ensembles).
- 3D structure (S3F and S3F-MSA rank 1st and 3rd; SaProt, ProtSSN and VenusREM also use it).
  Buried residues in a folded domain tolerate fewer swaps than exposed or disordered ones.
- A substitution matrix with conservation (VESPA and VespaG).

Arm 6 borrows these ideas in the simplest form we can run without new training: the
standardised ESM-2 and ESM-C scores, a burial score from the AlphaFold model of BRCA1
(weighted contact number, scaled by AlphaFold's own confidence), and BLOSUM62. Parts are
combined by an equal-weight average of z-scores, the way ProteinGym builds its ensembles, so
no weight is fitted to labels. Alignment-based models (EVE, TranceptEVE, PoET) are not copied:
they need a multiple sequence alignment we do not have. AlphaMissense is not used as a part:
it was calibrated on ClinVar, which already draws on this lab's data.

Arm 6 was added after the one test run, so it is exploratory. Only the choice between three
fixed versions and the cutoff use labels, both on the start region, frozen in their own file
before one run on the end region. Steps, in order: `features` (no labels), `freeze`,
`evaluate` (once), `writeup`.
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

from src import arm4, baselines, config, cutoffs, load, reference

SENTENCE_BETTER: str = "On the held-out region the added arm ranked harmful mutations better than `arm2`."
SENTENCE_NO_DIFFERENCE: str = "No detectable difference from `arm2` at this sample size."
SENTENCE_WORSE: str = "On the held-out region the added arm ranked harmful mutations worse than `arm2`."

_PARAGRAPH_MARKER: str = "<!-- arm6.json -->"
_LIMITS_ANCHOR: str = "**Limits:** (1) Blind to splicing"
_PARTS: dict[str, tuple[str, ...]] = {
    "plm": ("z_esm2", "z_esmc"),
    "plm_structure": ("z_esm2", "z_esmc", "z_burial"),
    "plm_structure_blosum": ("z_esm2", "z_esmc", "z_burial", "z_blosum"),
}


# ---------------------------------------------------------------- pure helpers (tested)

def read_structure(text: str) -> pd.DataFrame:
    """Per-residue position, letter, pLDDT and the atom used for contacts (CB; CA for Gly)."""
    three_to_one = config.THREE_TO_ONE
    rows: dict[int, dict] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        res = line[17:20].strip()
        pos = int(line[22:26])
        if res not in three_to_one:
            continue
        entry = rows.setdefault(
            pos, {"position": pos, "letter": three_to_one[res], "plddt": float(line[60:66])}
        )
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if atom == "CA":
            entry["ca"] = xyz
        elif atom == "CB":
            entry["cb"] = xyz
    out = []
    for pos in sorted(rows):
        entry = rows[pos]
        point = entry.get("cb", entry.get("ca"))
        if point is None:
            raise ValueError(f"Structure check failed: residue {pos} has no CA or CB atom.")
        out.append({"position": pos, "letter": entry["letter"], "plddt": entry["plddt"],
                    "x": point[0], "y": point[1], "z": point[2]})
    return pd.DataFrame(out)


def weighted_contact_number(coords: np.ndarray) -> np.ndarray:
    """Sum of 1/d^2 to every other residue: high for buried residues, low for exposed ones."""
    coords = np.asarray(coords, dtype=float)
    diff = coords[:, None, :] - coords[None, :, :]
    d2 = (diff ** 2).sum(axis=-1)
    np.fill_diagonal(d2, np.inf)
    return (1.0 / d2).sum(axis=1)


def burial_score(structure: pd.DataFrame) -> pd.Series:
    """Burial scaled by AlphaFold confidence, so disordered stretches count as exposed."""
    wcn = weighted_contact_number(structure[["x", "y", "z"]].to_numpy())
    return pd.Series(wcn * structure["plddt"].to_numpy() / 100.0, index=structure["position"].to_numpy())


def zscore(values: np.ndarray) -> np.ndarray:
    """Standardise to mean 0 and SD 1 over the given rows (no labels involved)."""
    values = np.asarray(values, dtype=float)
    sd = values.std()
    if not np.isfinite(sd) or sd == 0:
        raise ValueError("Arm 6 check failed: a part has no spread.")
    return (values - values.mean()) / sd


def candidate_harm(features: pd.DataFrame, candidate: str) -> np.ndarray:
    """Equal-weight mean of the candidate's z-scored parts; higher means more likely harmful."""
    if candidate not in _PARTS:
        raise ValueError(f"Unknown arm 6 candidate '{candidate}'.")
    return features[list(_PARTS[candidate])].to_numpy(dtype=float).mean(axis=1)


def choose_candidate(start_aurocs: dict[str, float]) -> str:
    """Highest start-region AUROC; a tie goes to the simpler (earlier) version."""
    best = config.ARM6_CANDIDATES[0]
    for name in config.ARM6_CANDIDATES[1:]:
        if start_aurocs[name] > start_aurocs[best]:
            best = name
    return best


def choose_sentence(difference: dict) -> str:
    """Map the bootstrap interval of arm 6 minus arm 2 to one fixed sentence."""
    if difference["low"] > 0:
        return SENTENCE_BETTER
    if difference["high"] < 0:
        return SENTENCE_WORSE
    return SENTENCE_NO_DIFFERENCE


# ---------------------------------------------------------------- file helpers

def _append_log(status: str, reason: str) -> None:
    """Append one line to the arm 6 run log. Only this module writes that file."""
    stamp = datetime.now().astimezone().isoformat(timespec="seconds")
    clean = reason.replace("|", " ").replace("\n", " ")
    line = (
        f"{stamp} | run=1 | status={status} | reason={clean} | "
        f"arm6_features_sha256={config.sha256_of_file(config.ARM6_FEATURES_PATH)} | "
        f"frozen_arm6_sha256={config.sha256_of_file(config.ARM6_FROZEN_PATH)}\n"
    )
    config.ARM6_RUN_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(config.ARM6_RUN_LOG_PATH, "a", newline="\n", encoding="utf-8") as handle:
        handle.write(line)


def read_features() -> pd.DataFrame:
    """Read the arm 6 feature file (one row per scored swap, no labels)."""
    if not config.ARM6_FEATURES_PATH.is_file():
        raise FileNotFoundError("outputs/arm6_features.csv not found. Run: python -m src.arm6 features")
    return pd.read_csv(config.ARM6_FEATURES_PATH, keep_default_na=False)


def _merged(region: str, labels: pd.DataFrame) -> pd.DataFrame:
    """Join one region's labels with the arm 6 features and the arm 2 score."""
    feats = read_features()
    merged = labels.merge(feats, on="mutant", validate="one_to_one", suffixes=("", "_f"))
    if len(merged) != len(labels):
        raise ValueError("Arm 6 check failed: a lab-tested mutation has no features.")
    return merged


# ---------------------------------------------------------------- steps

def run_features() -> int:
    """Step 1: build every part for all 6,194 scored swaps. Uses no labels."""
    if config.ARM6_FEATURES_PATH.exists():
        print("Refusing to rebuild: outputs/arm6_features.csv already exists.", file=sys.stderr)
        return 1
    if config.sha256_of_file(config.ARM6_STRUCTURE_PATH) != config.ARM6_STRUCTURE_SHA256:
        raise ValueError("Structure check failed: AlphaFold file hash differs from the recorded one.")
    structure = read_structure(config.ARM6_STRUCTURE_PATH.read_text(encoding="utf-8"))
    seq = load.load_wild_type()
    if len(structure) != len(seq) or "".join(structure["letter"]) != seq:
        raise ValueError("Structure check failed: AlphaFold sequence differs from our wild type.")
    burial = burial_score(structure)
    plddt = pd.Series(structure["plddt"].to_numpy(), index=structure["position"].to_numpy())

    esm2 = load.load_scores()[["mutant", "position", "wild_type", "new_letter", "esm_score"]]
    esmc = arm4.read_arm4_scores()[["mutant", "esm_score"]].rename(columns={"esm_score": "esmc_score"})
    feats = esm2.merge(esmc, on="mutant", validate="one_to_one")
    if len(feats) != config.N_SCORE_ROWS:
        raise ValueError("Arm 6 check failed: ESM-2 and ESM-C score files do not line up.")
    feats["burial"] = feats["position"].map(burial)
    feats["plddt"] = feats["position"].map(plddt)
    feats["blosum"] = [
        baselines.blosum_score(w, n) for w, n in zip(feats["wild_type"], feats["new_letter"])
    ]
    # Every part is turned so that higher means more likely harmful, then standardised
    # over all scored swaps.
    feats["z_esm2"] = zscore(0.0 - feats["esm_score"].to_numpy())
    feats["z_esmc"] = zscore(0.0 - feats["esmc_score"].to_numpy())
    feats["z_burial"] = zscore(feats["burial"].to_numpy())
    feats["z_blosum"] = zscore(0.0 - feats["blosum"].to_numpy(dtype=float))
    feats = feats.drop(columns=["position", "wild_type", "new_letter"])
    tmp = config.ARM6_FEATURES_PATH.parent / f"{config.ARM6_FEATURES_PATH.name}.tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        feats.to_csv(handle, index=False, lineterminator="\n", float_format="%.6f")
    os.replace(tmp, config.ARM6_FEATURES_PATH)
    print(f"Wrote outputs/arm6_features.csv: {len(feats)} rows.")
    print(f"sha256 = {config.sha256_of_file(config.ARM6_FEATURES_PATH)}")
    return 0


def run_freeze() -> int:
    """Step 2: choose one of three fixed versions and its cutoff on the start region only."""
    if config.ARM6_FROZEN_PATH.exists():
        print("Refusing to overwrite outputs/frozen_arm6.json.", file=sys.stderr)
        return 1
    if arm4.log_statuses(config.ARM6_RUN_LOG_PATH):
        print("Refusing to freeze: the arm 6 run has started.", file=sys.stderr)
        return 1
    if not load.test_run_done():
        print("Refusing to freeze: arm 6 is an addition after the one test run.", file=sys.stderr)
        return 1
    start = arm4._binary_rows(_merged("start", load.load_labels("start", "tuning")))
    if len(start) != 516:
        raise ValueError("Arm 6 check failed: expected 516 start-region LOF and FUNC rows.")
    y = start["is_harmful"].to_numpy()
    aurocs = {name: float(roc_auc_score(y, candidate_harm(start, name))) for name in config.ARM6_CANDIDATES}
    chosen = choose_candidate(aurocs)
    harm = candidate_harm(start, chosen)
    cutoff = cutoffs.choose_catch_cutoff(harm, y, config.CATCH_TARGET_PERCENT)
    counts = cutoffs.count_outcomes(harm, y, cutoff)
    frozen = {
        "what": (
            "EXPLORATORY arm 6, added after the one test run: an equal-weight ensemble of ideas "
            "used by models above ESM-2 on ProteinGym. Version and cutoff chosen on the start "
            "region only. Nothing here changes outputs/frozen.json."
        ),
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "seed": config.SEED,
        "tuning_region": "start",
        "test_region": "end",
        "selection_rule": "highest start-region AUROC (LOF vs FUNC, 516 rows); a tie goes to the simpler version",
        "candidates_start_auroc": aurocs,
        "arm2_start_auroc_in_sample": float(roc_auc_score(y, 0.0 - start["esm_score"].to_numpy())),
        "chosen": chosen,
        "chosen_name": config.ARM6_CANDIDATE_NAMES[chosen],
        "parts": list(_PARTS[chosen]),
        "cutoff_rule": "highest cutoff catching at least 90% of start-region harmful rows",
        "cutoff": float(cutoff),
        "tuning_caught": int(counts["harmful_caught"]),
        "tuning_false_alarms": int(counts["false_alarms"]),
        "inputs": {
            "arm6_features_sha256": config.sha256_of_file(config.ARM6_FEATURES_PATH),
            "structure_sha256": config.ARM6_STRUCTURE_SHA256,
            "arm4_scores_sha256": config.sha256_of_file(config.ARM4_SCORES_PATH),
            "esm_scores_sha256": config.sha256_of_file(config.ESM_SCORES_PATH),
            "frozen_sha256": config.sha256_of_file(config.FROZEN_PATH),
        },
    }
    arm4._write_json(frozen, config.ARM6_FROZEN_PATH)
    print(json.dumps(frozen, indent=2))
    print("Start-region numbers are in-sample; they are not results.")
    return 0


def paragraph_text(frozen: dict, result: dict) -> str:
    """The one write-up paragraph for arm 6, built only from the two arm 6 files."""
    aurocs = frozen["candidates_start_auroc"]
    end = result["end"]
    diff = result["arm6_minus_arm2"]
    others = ", ".join(f"{aurocs[n]:.4f}" for n in config.ARM6_CANDIDATES if n != frozen["chosen"])
    return (
        "**Extra arm, added after the test run (exploratory):** `arm6` borrows ideas from models "
        "that rank above ESM-2 on ProteinGym for this table: several language models together, "
        "3D structure (burial in the AlphaFold model of BRCA1, as S3F and SaProt use structure) "
        "and a substitution matrix (as VESPA does), combined by an equal-weight average of "
        f"standardised scores. The version used, {frozen['chosen_name']}, was chosen from three by "
        f"start-region AUROC alone ({aurocs[frozen['chosen']]:.4f} vs {others}), with its cutoff "
        "set by the same 90% rule on `1–101` and frozen in `outputs/frozen_arm6.json` before one "
        "run on the end region. It caught "
        f"{end['harmful_caught']} of {end['n_harmful']} harmful mutations ({end['harmful_missed']} "
        f"missed) with {end['false_alarms']} of {end['n_harmless']} false alarms; AUROC "
        f"{end['auroc']:.4f} vs {end['auroc_arm2_same_rows']:.4f} for `arm2` (difference "
        f"`{diff['point']:+.3f}`, 95% interval `{diff['low']:+.3f}` to `{diff['high']:+.3f}`). "
        f"{result['sentence']} It was designed after seeing the end-region results, so it is "
        f"exploratory and the headline table is unchanged. {_PARAGRAPH_MARKER}\n"
    )


def run_evaluate() -> int:
    """Step 3: apply the frozen arm 6 to the end region, once."""
    if not load.test_run_done():
        print("Refusing to run: the one test run is not done.", file=sys.stderr)
        return 1
    if not config.ARM6_FROZEN_PATH.is_file():
        print("Refusing to run: outputs/frozen_arm6.json not found. Run freeze first.", file=sys.stderr)
        return 1
    if arm4.log_statuses(config.ARM6_RUN_LOG_PATH):
        print("Refusing to run: arm 6 has already run. Ask a human.", file=sys.stderr)
        return 1
    with open(config.ARM6_FROZEN_PATH, "r", encoding="utf-8") as handle:
        frozen = json.load(handle)
    if config.sha256_of_file(config.ARM6_FEATURES_PATH) != frozen["inputs"]["arm6_features_sha256"]:
        print("HASH MISMATCH: outputs/arm6_features.csv changed since freezing.", file=sys.stderr)
        return 1

    _append_log("started", "arm 6 run")
    try:
        end_all = _merged("end", load.load_labels("end", "post_test"))
        end_all["harm_arm6"] = candidate_harm(end_all, frozen["chosen"])
        end_all["harm_arm2"] = 0.0 - end_all["esm_score"]
        end = arm4._binary_rows(end_all)
        if len(end_all) != 1262 or len(end) != 1166:
            raise ValueError("Arm 6 check failed: unexpected end-region row counts.")
        y = end["is_harmful"].to_numpy()
        harm6 = end["harm_arm6"].to_numpy()
        harm2 = end["harm_arm2"].to_numpy()
        counts = cutoffs.count_outcomes(harm6, y, frozen["cutoff"])
        difference = reference.paired_position_bootstrap(
            end["position"].to_numpy(), y, harm6, harm2, config.N_BOOTSTRAP, config.SEED
        )
        int_rows = end_all.loc[end_all["findlay_class"] == config.EXCLUDED_CLASS]
        result = {
            "what": (
                "EXPLORATORY arm 6 on the held-out end region, run once after the headline test "
                "run. Same 1,166 rows as the headline. Not a replacement for arm 2."
            ),
            "chosen": frozen["chosen"],
            "chosen_name": frozen["chosen_name"],
            "cutoff": float(frozen["cutoff"]),
            "end": {
                **{k: int(v) for k, v in counts.items()},
                **cutoffs.rates(counts),
                "auroc": round(float(roc_auc_score(y, harm6)), 4),
                "auroc_arm2_same_rows": round(float(roc_auc_score(y, harm2)), 4),
            },
            "arm6_minus_arm2": difference,
            "int_rows_end": {
                "n_int": int(len(int_rows)),
                "flagged": int(cutoffs.flag(int_rows["harm_arm6"].to_numpy(), frozen["cutoff"]).sum()),
            },
            "sentence": choose_sentence(difference),
        }
        arm4._write_json(result, config.ARM6_RESULT_PATH)
        arm4._write_text(paragraph_text(frozen, result), config.ARM6_PARAGRAPH_PATH)
    except BaseException as exc:
        _append_log("failed", str(exc) or type(exc).__name__)
        raise
    _append_log("success", "arm 6 run")
    print(json.dumps(result, indent=2))
    print(config.ARM6_PARAGRAPH_PATH.read_text(encoding="utf-8"))
    return 0


def run_writeup() -> int:
    """Step 4: put the arm 6 paragraph and the AlphaFold credit into the write-up, once."""
    if "success" not in arm4.log_statuses(config.ARM6_RUN_LOG_PATH):
        print("Refusing to edit the write-up: arm 6 has no successful run.", file=sys.stderr)
        return 1
    paragraph = config.ARM6_PARAGRAPH_PATH.read_text(encoding="utf-8")
    credit = (
        "\n- AlphaFold Protein Structure Database model of BRCA1 (AF-P38398-F1, v6; Jumper J et al., "
        f"2021, *Nature* 596:583–589; Varadi M et al., 2024): {config.ARM6_STRUCTURE_URL} — used "
        "only for the burial part of the exploratory `arm6`."
    )
    for name in ("README.md", "writeup/writeup.md"):
        path = config.ROOT / name
        if path.is_file() and _LIMITS_ANCHOR in path.read_text(encoding="utf-8"):
            print(arm4._insert_once(path, _LIMITS_ANCHOR, paragraph + "\n", _PARAGRAPH_MARKER, after=False))
    credits = config.ROOT / "writeup" / "credits.md"
    if credits.is_file() and "AF-P38398-F1-model_v6" not in credits.read_text(encoding="utf-8"):
        arm4._write_text(credits.read_text(encoding="utf-8").rstrip("\n") + credit + "\n", credits)
        print("credits.md: updated")
    print("Write-up updated.")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Command line for the four arm 6 steps: features, freeze, evaluate, writeup."""
    parser = argparse.ArgumentParser(description="Exploratory arm 6 (ensemble).")
    parser.add_argument("step", choices=["features", "freeze", "evaluate", "writeup"])
    args = parser.parse_args(argv)
    config.set_seeds()
    try:
        if args.step == "features":
            return run_features()
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
