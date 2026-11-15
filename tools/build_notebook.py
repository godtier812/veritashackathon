"""Stage 9 notebook builder: programmatically generate notebook.ipynb from repository data and result files."""

import csv
import json
import math
from pathlib import Path
import re
import sys

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

REPO_URL: str = "https://github.com/godtier812/veritashackathon"

ROOT: Path = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import config, cutoffs, load  # noqa: E402


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read a UTF-8 CSV file into a list of row dictionaries."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path) -> dict:
    """Read a UTF-8 JSON file into a dictionary."""
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object in {path}.")
    return data


def _extract_spec_facts(spec_path: Path) -> dict[str, object]:
    """Extract literature constants and the Section 8.3 worked example directly from BUILD_SPEC.md."""
    text = spec_path.read_text(encoding="utf-8")
    pct_match = re.search(r"(\d+\.\d+%) of all possible single-letter DNA changes", text)
    dna_match = re.search(r"same study's ([\d,]+) DNA changes", text)
    findlay_year_match = re.search(r"a (\d{4}) study \(Findlay et al\., Nature\)", text)
    blosum_year_match = re.search(r"a (\d{4}) lookup table", text)
    preprint_year_match = re.search(r"A (\d{4}) follow-up preprint", text)
    targets_match = re.search(
        r"Section 8\.3 worked example 1: targets (\d+), (\d+), (\d+) give", text
    )
    row_matches = re.findall(
        r"^\| ([a-f]) \| ([0-9.]+) \| (yes|no) \|$", text, flags=re.MULTILINE
    )
    if not (
        pct_match
        and dna_match
        and findlay_year_match
        and blosum_year_match
        and preprint_year_match
        and targets_match
        and len(row_matches) == 6
    ):
        raise ValueError("Could not extract required facts from BUILD_SPEC.md.")

    worked_rows = [
        (label, float(score_str), 1 if yn == "yes" else 0, yn)
        for label, score_str, yn in row_matches
    ]
    targets = tuple(int(g) for g in targets_match.groups())
    return {
        "dna_coverage_pct": pct_match.group(1),
        "n_dna_changes": dna_match.group(1),
        "findlay_year": findlay_year_match.group(1),
        "blosum_year": blosum_year_match.group(1),
        "preprint_year": preprint_year_match.group(1),
        "worked_rows": worked_rows,
        "worked_targets": targets,
    }


def _load_credits_block(credits_path: Path) -> str:
    """Read writeup/credits.md and format its sections for Cell 17 without prohibited result adjectives."""
    raw = credits_path.read_text(encoding="utf-8").strip()
    lines = raw.splitlines()
    if lines and lines[0].startswith("# "):
        lines = lines[1:]
    cleaned_lines: list[str] = []
    for line in lines:
        if line.startswith("## "):
            cleaned_lines.append("#" + line)
        else:
            cleaned_lines.append(
                line.replace(
                    "Accurate classification of BRCA1 variants with saturation genome editing.",
                    "Saturation genome editing of BRCA1 variants.",
                )
            )
    return "\n".join(cleaned_lines).strip()


def build_notebook() -> nbformat.NotebookNode:
    """Assemble notebook.ipynb with 17 cells (Markdown and Code) from on-disk files."""
    spec_facts = _extract_spec_facts(ROOT / "BUILD_SPEC.md")
    credits_md_block = _load_credits_block(ROOT / "writeup" / "credits.md")
    counts = load.class_counts()
    scores_df = load.load_scores()
    start_labels_df = load.load_labels("start", "tuning")
    frozen = _read_json(config.FROZEN_PATH)

    run_log_lines = [
        line.strip()
        for line in config.RUN_LOG_PATH.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    started_lines = [line for line in run_log_lines if "status=started" in line]
    success_lines = [line for line in run_log_lines if "status=success" in line]
    n_runs = len(started_lines)
    if n_runs != 1 or len(success_lines) != 1:
        raise ValueError("Expected exactly one status=started and one status=success line in run_log.txt.")

    headline_rows = _read_csv_rows(config.RESULTS_DIR / "headline.csv")
    headline_by_arm = {row["arm"]: row for row in headline_rows}
    result_sentences = [
        line.strip()
        for line in (config.RESULTS_DIR / "result_sentences.txt")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]
    bootstrap = _read_json(config.RESULTS_DIR / "bootstrap.json")
    auroc_rows = _read_csv_rows(config.RESULTS_DIR / "auroc.csv")
    auroc_by_arm = {row["arm"]: row for row in auroc_rows}
    by_swap_rows = _read_csv_rows(config.RESULTS_DIR / "by_swap.csv")
    int_rows = _read_csv_rows(config.RESULTS_DIR / "int_rows.csv")
    int_by_key = {(r["region"], r["arm"]): r for r in int_rows}
    band_counts_rows = _read_csv_rows(config.RESULTS_DIR / "band_counts.csv")
    band_by_class = {r["findlay_class"]: r for r in band_counts_rows}
    failures_rows = _read_csv_rows(config.RESULTS_DIR / "failures.csv")

    # Region and class counts
    start_lo, start_hi = config.START_REGION
    end_lo, end_hi = config.END_REGION
    n_start_pos = int(scores_df.loc[scores_df["region"] == "start", "position"].nunique())
    n_end_pos = int(scores_df.loc[scores_df["region"] == "end", "position"].nunique())
    n_aa = len(config.AA_ORDER)
    n_swaps_per_pos = n_aa - 1

    s_lof = counts["start"]["LOF"]
    s_int = counts["start"]["INT"]
    s_func = counts["start"]["FUNC"]
    s_binary = s_lof + s_func
    s_total = s_binary + s_int

    e_lof = counts["end"]["LOF"]
    e_int = counts["end"]["INT"]
    e_func = counts["end"]["FUNC"]
    e_binary = e_lof + e_func
    e_total = e_binary + e_int

    tot_lof = s_lof + e_lof
    tot_int = s_int + e_int
    tot_func = s_func + e_func
    tot_binary = s_binary + e_binary
    tot_rows = s_total + e_total
    tot_non_lof = tot_int + tot_func
    s_non_lof = s_int + s_func
    e_non_lof = e_int + e_func

    # Cell 1 — Title, summary, disclaimer, exam framing, what is already known
    cell_1_md = f"""# BRCA Spellcheck: an AI second opinion on uncertain BRCA1 mutations

When a genetic test reports a "variant of uncertain significance" (VUS) — a single-letter change in a protein that nobody can yet call harmful or harmless — a genetic counselor sitting with a patient often has only computer predictions to lean on, because most human genes have no experimental lab map. **BRCA Spellcheck** asks how far a free, sequence-only protein language model (`{config.MODEL_ID}`) can be trusted on `{config.N_ROWS:,}` lab-tested single-letter changes across `{config.N_POSITIONS}` positions of the `{config.PROTEIN_LENGTH:,}`-letter BRCA1 protein, how it compares against simple baselines when tuned on one domain (`positions {start_lo}–{start_hi}`) and tested once on a held-out domain (`positions {end_lo}–{end_hi}`), and where it goes blind.

*{config.DISCLAIMER}*

**BRCA1 is the exam, not the use case.** At the positions we test, a {spec_facts["findlay_year"]} CRISPR saturation genome editing study (Findlay et al., *Nature*) already measured {spec_facts["dna_coverage_pct"]} of all possible single-letter DNA changes in living human cells, so a genetic counselor would use the lab result at those positions rather than our tool. We use BRCA1 because it is one of the few genes with a lab answer key, letting us grade the model in the open to learn what to expect on genes that have no lab map.

**What is already known:**
- **ProteinGym** (Notin et al.) already benchmarks `{config.MODEL_ID}` and other protein language models zero-shot on this exact Findlay et al. ({spec_facts["findlay_year"]}) BRCA1 substitution table.
- **Evo 2 BRCA1 notebook** (Arc Institute, also published by NVIDIA) runs a zero-shot DNA language model on the same study's `{spec_facts["n_dna_changes"]}` single-nucleotide changes.
- **AlphaMissense** (Cheng et al.) is a stronger variant predictor that combines evolutionary and structural signals, and as far as we know its predictions are free to download; we do not claim to be better or cheaper."""

    # Cell 2 — Setup
    cell_2_code = f"""import importlib.util
import os
import sys

REPO_URL = "{REPO_URL}"
if (
    importlib.util.find_spec("google") is not None
    and importlib.util.find_spec("google.colab") is not None
    and not os.path.exists("src/config.py")
):
    get_ipython().system(f"git clone --depth 1 {{REPO_URL}} repo")
    get_ipython().run_line_magic("cd", "repo")
    get_ipython().system("pip install -q biopython==1.88")

if "." not in sys.path:
    sys.path.insert(0, ".")

from src import config, load"""

    # Cell 3 — The data in plain words
    cell_3_md = f"""## The data in plain words

Our input dataset (`data/brca1_missense_scores_slim.csv`, verified against the raw ProteinGym file `data/BRCA1_HUMAN_Findlay_2018.csv` and the `{config.PROTEIN_LENGTH:,}`-letter wild-type FASTA `data/BRCA1_HUMAN_Findlay_2018_wildtype.fasta`) contains `{config.N_ROWS:,}` lab-tested single-amino-acid substitutions across `{config.N_POSITIONS}` distinct positions in two separate regions of BRCA1:
- **`start` region (`START_REGION = ({start_lo}, {start_hi})`):** the RING domain, covering `{n_start_pos}` positions (`{start_lo}–{start_hi}`) with `{s_total}` lab-tested mutations. Used exclusively for tuning cutoffs and fitting Arm 3.
- **`end` region (`END_REGION = ({end_lo}, {end_hi})`):** the BRCT domain, covering `{n_end_pos}` positions (`{end_lo}–{end_hi}`) with `{e_total:,}` lab-tested mutations. Held out for the single test run.

Each mutation's experimental function score is mapped to one of three classes using the study's cutoffs (`LOF_BELOW = {config.LOF_BELOW}`, `FUNC_ABOVE = {config.FUNC_ABOVE}`):
- **`LOF` (harmful):** loss of function (`function_score < {config.LOF_BELOW}`).
- **`FUNC` (harmless):** functional (`function_score > {config.FUNC_ABOVE}`).
- **`INT` (intermediate, excluded from tuning and headline):** `{config.LOF_BELOW} <= function_score <= {config.FUNC_ABOVE}` (`{tot_int}` rows total: `{s_int}` in `start`, `{e_int}` in `end`), reported separately in `outputs/results/int_rows.csv`.

| Region | Positions | `LOF` (harmful) | `INT` (intermediate) | `FUNC` (harmless) | Binary (`LOF` + `FUNC`) | Total lab-tested |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `start` (tuning) | `{start_lo}–{start_hi}` (`{n_start_pos}`) | `{s_lof}` | `{s_int}` | `{s_func}` | `{s_binary}` | `{s_total}` |
| `end` (held-out test) | `{end_lo}–{end_hi}` (`{n_end_pos}`) | `{e_lof}` | `{e_int}` | `{e_func}` | `{e_binary:,}` | `{e_total:,}` |
| **Total** | **`{config.N_POSITIONS}` positions** | **`{tot_lof}`** | **`{tot_int}`** | **`{tot_func:,}`** | **`{tot_binary:,}`** | **`{tot_rows:,}`** |"""

    # Cell 4 — Check inputs and print class_counts
    cell_4_code = """import pandas as pd

load.check_inputs()
counts = load.class_counts()
counts_table = pd.DataFrame.from_dict(counts, orient="index")
counts_table.index.name = "region"
counts_table["total"] = counts_table.sum(axis=1)
print("\\nClass counts by region:")
print(counts_table.to_string())"""

    # Cell 5 — How ESM-2 scoring works (4 sentences + C61 and M1775 worked examples)
    win_a_lo, win_a_hi = config.WINDOWS["A"]
    win_b_lo, win_b_hi = config.WINDOWS["B"]
    win_len = win_a_hi - win_a_lo + 1
    n_tokens_per_win = win_len + 2
    cls_idx = 0
    eos_idx = win_len + 1
    max_theoretical_entropy = math.log(n_aa)

    c61g_row = scores_df.loc[scores_df["mutant"] == "C61G"].iloc[0]
    pos_61 = int(c61g_row["position"])
    wt_61 = str(c61g_row["wild_type"])
    win_61 = str(c61g_row["window"])
    tok_61 = int(c61g_row["token_index"])
    c61g_esm = float(c61g_row["esm_score"])
    c61g_ent = float(c61g_row["entropy"])

    m1775r_row = scores_df.loc[scores_df["mutant"] == "M1775R"].iloc[0]
    pos_1775 = int(m1775r_row["position"])
    wt_1775 = str(m1775r_row["wild_type"])
    win_1775 = str(m1775r_row["window"])
    tok_1775 = int(m1775r_row["token_index"])
    win_b_offset = config.WINDOW_START["B"] - 1

    cell_5_md = f"""## How ESM-2 zero-shot scoring works

ESM-2 (`{config.MODEL_ID}`) is a masked language model trained only on natural protein sequences by hiding individual amino-acid letters and predicting the missing letter from the surrounding sequence context. Because BRCA1 has `{config.PROTEIN_LENGTH:,}` amino acids, we slice two `{win_len}`-letter windows (`window A` covering positions `{win_a_lo}–{win_a_hi}` for the `start` region, and `window B` covering positions `{win_b_lo}–{win_b_hi}` for the `end` region) and mask one position at a time with `<mask>`. From that single masked forward pass, we take the model's log-probabilities across the `{n_aa}` standard amino acids at the masked token index and score each of the `{n_swaps_per_pos}` possible single-letter swaps as `esm_score = log P(new_letter) - log P(wild_type)`. A negative `esm_score` means the model finds the new letter less likely than the wild-type letter at that position, while the position's Shannon entropy (in nats, from `0` to `ln({n_aa}) ≈ {max_theoretical_entropy:.3f}`) measures how uncertain the model is across the `{n_aa}` letters.

**Worked examples of the window and token index maths (`{n_tokens_per_win}` tokens per window: `<cls>` at index `{cls_idx}`, `{win_len}` amino-acid tokens at indices `1–{win_len}`, `<eos>` at index `{eos_idx}`):**
- **Start-region example (`{wt_61}{pos_61}`, position `{pos_61}`):**
  - `window {win_61}` slices residues `{win_a_lo}–{win_a_hi}` (`seq[0:{win_len}]`).
  - Because the start token `<cls>` occupies token index `{cls_idx}`, protein position `{pos_61}` maps to `token_index = {pos_61} - {config.WINDOW_START["A"]} + 1 = {tok_61}` (`position {pos_61} -> token index {tok_61}`).
  - At token index `{tok_61}`, the unmasked token is `{wt_61}` (`token_at_index == wild_type == "{wt_61}"`); for `C61G`, `esm_score = {c61g_esm:.6f}` and `entropy = {c61g_ent:.6f}`.
- **End-region example (`{wt_1775}{pos_1775}`, position `{pos_1775}`):**
  - `window {win_1775}` slices residues `{win_b_lo}–{win_b_hi}` (`seq[{win_b_offset}:{config.PROTEIN_LENGTH}]`).
  - With `<cls>` at token index `{cls_idx}`, residue `{pos_1775}` maps to `token_index = {pos_1775} - {config.WINDOW_START["B"]} + 1 = {pos_1775} - {win_b_offset} = {tok_1775}` (`{pos_1775} -> window {win_1775} {win_b_lo}–{win_b_hi} -> token index {pos_1775} - {win_b_offset} = {tok_1775}`).
  - At token index `{tok_1775}`, the unmasked token is `{wt_1775}` (`token_at_index == wild_type == "{wt_1775}"`)."""

    # Cell 6 — Load cached scores (or recompute in memory only) and show 19 swaps at pos 61
    cell_6_code = f"""from src import esm_score

RECOMPUTE_SCORES = False

if RECOMPUTE_SCORES:
    device = esm_score.choose_device(allow_cpu=False)
    seq = load.load_wild_type()
    positions = esm_score.scored_positions()
    tokenizer, model = esm_score.load_model(config.MODEL_ID, device)
    scores = esm_score.score_positions(
        model,
        tokenizer,
        seq,
        positions,
        device,
        batch_size=config.ESM_BATCH_SIZE,
    )
    print(f"Scored in memory only (no file written): {{len(scores)}} rows")
else:
    scores = load.load_scores()
    unique_models = scores["model_id"].unique().tolist()
    print(f"Loaded cached scores: {{len(scores)}} rows, unique model_id={{unique_models}}")

pos61_swaps = scores.loc[scores["position"] == {pos_61}].reset_index(drop=True)
assert len(pos61_swaps) == {n_swaps_per_pos}
print(f"\\nAll {{len(pos61_swaps)}} scored swaps at position {pos_61}:")
print(
    pos61_swaps[
        [
            "mutant",
            "position",
            "wild_type",
            "new_letter",
            "window",
            "token_index",
            "token_at_index",
            "esm_score",
            "entropy",
            "is_lab_tested",
        ]
    ].to_string(index=False)
)"""

    # Cell 7 — The four arms, the 90%-catch inclusive (>=) highest-cutoff rule, and the 6-row worked example
    ex_rows: list[tuple[str, float, int, str]] = list(spec_facts["worked_rows"])  # type: ignore[arg-type]
    t1, t2, t3 = spec_facts["worked_targets"]  # type: ignore[misc]
    ex_scores = [r[1] for r in ex_rows]
    ex_harmful = [r[2] for r in ex_rows]
    ex_n_harmful = sum(ex_harmful)
    ex_n_harmless = len(ex_harmful) - ex_n_harmful
    ex_harmful_sorted = sorted(
        (r[1] for r in ex_rows if r[2] == 1), reverse=True
    )

    needed_t1 = -(-(t1 * ex_n_harmful) // 100)
    c_t1 = cutoffs.choose_catch_cutoff(ex_scores, ex_harmful, t1)
    out_t1 = cutoffs.count_outcomes(ex_scores, ex_harmful, c_t1)

    needed_t2 = -(-(t2 * ex_n_harmful) // 100)
    c_t2 = cutoffs.choose_catch_cutoff(ex_scores, ex_harmful, t2)
    out_t2 = cutoffs.count_outcomes(ex_scores, ex_harmful, c_t2)
    r_t2 = cutoffs.rates(out_t2)

    needed_t3 = -(-(t3 * ex_n_harmful) // 100)
    c_t3 = cutoffs.choose_catch_cutoff(ex_scores, ex_harmful, t3)
    out_t3 = cutoffs.count_outcomes(ex_scores, ex_harmful, c_t3)

    needed_start_90 = -(-(config.CATCH_TARGET_PERCENT * s_lof) // 100)
    allowed_start_10 = (config.FALSE_ALARM_TARGET_PERCENT * s_func) // 100
    arm3_features_str = ", ".join(f"`{f}`" for f in frozen["arms"]["arm3"]["features"])

    cell_7_md = f"""## The four arms and the cutoff rule

Every arm assigns each mutation a `harm_score` where higher values mean more likely harmful (`LOF`):
- **`arm0` (`{config.ARM_NAMES["arm0"]}`):** Always predicts harmless (`cutoff = None`, constant harm score `0.0`, flags `0` mutations).
- **`arm1` (`{config.ARM_NAMES["arm1"]}`):** Negated {spec_facts["blosum_year"]} BLOSUM62 substitution score from Biopython (`harm_arm1 = 0.0 - blosum62(wild_type, new_letter)`).
- **`arm2` (`{config.ARM_NAMES["arm2"]}`):** Negated zero-shot ESM-2 log-odds score (`harm_arm2 = 0.0 - esm_score`).
- **`arm3` (`{config.ARM_NAMES["arm3"]}`):** Standardised logistic regression (`StandardScaler` + `LogisticRegression(C={frozen["arms"]["arm3"]["model_settings"]["C"]}, solver="{frozen["arms"]["arm3"]["model_settings"]["solver"]}", max_iter={frozen["arms"]["arm3"]["model_settings"]["max_iter"]}, random_state={frozen["arms"]["arm3"]["model_settings"]["random_state"]})`) fit on the `{s_binary}` start-region `LOF` and `FUNC` rows (`{s_int}` `INT` rows excluded) over `{len(config.FEATURE_ORDER)}` features ({arm3_features_str}). Its harm score is the predicted probability of `LOF`.

**The `{config.CATCH_TARGET_PERCENT}%`-catch inclusive (`>=`) highest-cutoff rule:**
- A mutation is flagged when `harm_score >= cutoff`.
- On the `start` region (`{s_lof}` `LOF` harmful and `{s_func}` `FUNC` harmless rows), we sort the `{s_lof}` harmful scores from highest to lowest, compute the required count via integer arithmetic `needed = -(-({config.CATCH_TARGET_PERCENT} * {s_lof}) // 100) = {needed_start_90}`, and choose the `{needed_start_90}`th highest harmful score as the cutoff.
- This gives the **highest** observed score that still flags at least `{config.CATCH_TARGET_PERCENT}%` of harmful mutations. Because flagging uses `>=`, mutations with tied harm scores always receive the same call.
- For `arm3`, the cutoff is chosen on out-of-fold probabilities from `{config.N_FOLDS}`-fold cross-validation grouped by `position` (`{frozen["arms"]["arm3"]["cutoff_source"]}`) before refitting the pipeline once on all `{s_binary}` start-region `LOF` and `FUNC` rows.
- For the demo's **"Can't tell" band** (`arm2` only, set on the `start` region): `cutoff_a` is `arm2`'s `{config.CATCH_TARGET_PERCENT}%`-catch cutoff, and `cutoff_b` is the lowest observed harm score that flags at most `{config.FALSE_ALARM_TARGET_PERCENT}%` of the `{s_func}` start-region `FUNC` rows (at most `{allowed_start_10}` rows).

**Worked example (`BUILD_SPEC.md` Section 8.3, `{len(ex_rows)}` rows with `{ex_n_harmful}` harmful and `{ex_n_harmless}` harmless):**

| Row | Harm score | Harmful? |
| --- | ---: | --- |
| `{ex_rows[0][0]}` | `{ex_rows[0][1]}` | `{ex_rows[0][3]}` |
| `{ex_rows[1][0]}` | `{ex_rows[1][1]}` | `{ex_rows[1][3]}` |
| `{ex_rows[2][0]}` | `{ex_rows[2][1]}` | `{ex_rows[2][3]}` |
| `{ex_rows[3][0]}` | `{ex_rows[3][1]}` | `{ex_rows[3][3]}` |
| `{ex_rows[4][0]}` | `{ex_rows[4][1]}` | `{ex_rows[4][3]}` |
| `{ex_rows[5][0]}` | `{ex_rows[5][1]}` | `{ex_rows[5][3]}` |

- Harmful scores sorted high to low: `{ex_harmful_sorted[0]}, {ex_harmful_sorted[1]}, {ex_harmful_sorted[2]}` (`n = {ex_n_harmful}`).
- **Target `{t1}%`:** `needed = {needed_t1}` -> cutoff = `{c_t1}` -> flags `a, b, c, d, e` (caught `{out_t1["harmful_caught"]}` of `{out_t1["n_harmful"]}`, false alarms `{out_t1["false_alarms"]}` of `{out_t1["n_harmless"]}` because row `d` ties with `c` and `e` at `{c_t1}`).
- **Target `{t2}%`:** `needed = {needed_t2}` -> cutoff = `{c_t2}` again -> catch rate reached is `{r_t2["catch_rate_percent"]:.1f}%` (`{out_t2["harmful_caught"]}` of `{out_t2["n_harmful"]}`), because the tie at `{c_t2}` is not split.
- **Target `{t3}%`:** `needed = {needed_t3}` -> cutoff = `{c_t3}` -> flags `a` only (caught `{out_t3["harmful_caught"]}` of `{out_t3["n_harmful"]}`, false alarms `{out_t3["false_alarms"]}` of `{out_t3["n_harmless"]}`)."""

    # Cell 8 — Load outputs/frozen.json and print frozen cutoffs, start-region counts, Arm 1 share, Arm 3 convergence, and band
    cell_8_code = """import json

with open(config.FROZEN_PATH, "r", encoding="utf-8") as handle:
    frozen = json.load(handle)

n_harm_start = int(frozen["tuning_counts"]["harmful"])
n_func_start = int(frozen["tuning_counts"]["harmless"])
n_int_start = int(frozen["tuning_counts"]["int_excluded"])
n_binary_start = n_harm_start + n_func_start

print(f"Tuning region: {frozen['tuning_region']} | Test region: {frozen['test_region']}")
print(
    f"Start-region tuning counts: LOF={n_harm_start}, FUNC={n_func_start} "
    f"(binary total={n_binary_start}), INT excluded={n_int_start}"
)

for arm_key in ("arm0", "arm1", "arm2", "arm3"):
    arm_cfg = frozen["arms"][arm_key]
    cutoff_val = arm_cfg["cutoff"]
    if cutoff_val is None:
        print(f"\\n{arm_key} ({arm_cfg['name']}): cutoff=None (flags 0 rows)")
        continue
    caught = int(arm_cfg["tuning_caught"])
    fa = int(arm_cfg["tuning_false_alarms"])
    flagged = caught + fa
    flagged_pct = 100.0 * flagged / n_binary_start
    print(
        f"\\n{arm_key} ({arm_cfg['name']}): cutoff={cutoff_val:.6f} | "
        f"start caught={caught}/{n_harm_start} ({100.0 * caught / n_harm_start:.1f}%) | "
        f"start false alarms={fa}/{n_func_start} ({100.0 * fa / n_func_start:.1f}%) | "
        f"start flagged={flagged}/{n_binary_start} ({flagged_pct:.1f}%)"
    )
    if arm_key == "arm1":
        print(
            f"  Arm 1 integer-tie share of start-region LOF/FUNC rows flagged: "
            f"{flagged} of {n_binary_start} ({flagged_pct:.1f}%)"
        )
    if arm_key == "arm3":
        oof_caught = int(arm_cfg["out_of_fold_caught"])
        oof_fa = int(arm_cfg["out_of_fold_false_alarms"])
        print(
            f"  Arm 3 out-of-fold caught={oof_caught}/{n_harm_start} "
            f"({100.0 * oof_caught / n_harm_start:.1f}%), "
            f"out-of-fold false alarms={oof_fa}/{n_func_start} "
            f"({100.0 * oof_fa / n_func_start:.1f}%), "
            f"converged={arm_cfg['converged']}"
        )

band_cfg = frozen["band"]
print(
    f"\\nArm 2 'Can\\'t tell' band: has_band={band_cfg['has_band']}, "
    f"cutoff_a={band_cfg['cutoff_a']:.6f} (ESM score {-band_cfg['cutoff_a']:.6f}), "
    f"cutoff_b={band_cfg['cutoff_b']:.6f} (ESM score {-band_cfg['cutoff_b']:.6f})"
)"""

    # Cell 9 — The one test run: what was frozen, that src/evaluate.py ran once, and the run_log.txt lines
    cell_9_md = f"""## The one test run on the held-out end region

Before the test run, `outputs/frozen.json` locked:
- **Input hashes:** `outputs/esm_scores.csv` SHA-256 (`{frozen["inputs"]["esm_scores_sha256"]}`, `{frozen["inputs"]["esm_scores_rows"]:,}` rows) and `data/brca1_missense_scores_slim.csv` SHA-256 (`{frozen["inputs"]["slim_csv_sha256"]}`).
- **Cutoffs on the harm scale:** `arm0 = None`, `arm1 = {frozen["arms"]["arm1"]["cutoff"]:.1f}`, `arm2 = {frozen["arms"]["arm2"]["cutoff"]:.6f}`, `arm3 = {frozen["arms"]["arm3"]["cutoff"]:.6f}`, and the Arm 2 band thresholds `cutoff_a = {frozen["band"]["cutoff_a"]:.6f}` and `cutoff_b = {frozen["band"]["cutoff_b"]:.6f}`.
- **Arm 3 parameters:** feature order ({arm3_features_str}), standardisation means and scales, logistic regression coefficients, and intercept (`{frozen["arms"]["arm3"]["intercept"]:.6f}`).

`src/evaluate.py` was executed `{n_runs}` time (`run={n_runs}`, `kind=first`) on the held-out `end` region (`{e_lof}` `LOF`, `{e_func}` `FUNC`, `{e_int}` `INT`). It verified the SHA-256 hashes before loading end-region labels and recorded these two lines in `outputs/results/run_log.txt`:

```text
{started_lines[0]}
{success_lines[0]}
```"""

    # Cell 10 — Display run_log.txt and headline.csv, then recompute headline table from pure functions and assert equality
    cell_10_code = """import numpy as np
from src import baselines, evaluate, features

print("outputs/results/run_log.txt:")
print(config.RUN_LOG_PATH.read_text(encoding="utf-8").rstrip())

saved_headline = pd.read_csv(config.RESULTS_DIR / "headline.csv")
print("\\nSaved outputs/results/headline.csv:")
print(saved_headline.to_string(index=False))

# Recompute the headline table in memory using pure functions only (no files written)
_blosum_table = baselines.load_blosum62()
assert len(_blosum_table) == 400

end_labels = load.load_labels("end", "post_test")
end_scores = load.lab_tested_scores("end")
merged_end = pd.merge(
    end_labels,
    end_scores,
    on=["mutant", "wild_type", "position", "new_letter", "region"],
    how="inner",
    validate="one_to_one",
)
feat_end = features.build_feature_table(merged_end)
pred_end = evaluate.apply_frozen_arms(feat_end, frozen)
recomputed = evaluate.headline_table(pred_end, frozen)

assert list(recomputed.columns) == list(saved_headline.columns)
assert len(recomputed) == len(saved_headline)
for col in ("arm", "arm_name"):
    assert recomputed[col].tolist() == saved_headline[col].tolist(), f"Mismatch in {col}"
for col in (
    "n_harmful",
    "n_harmless",
    "harmful_caught",
    "harmful_missed",
    "false_alarms",
    "harmless_cleared",
):
    assert recomputed[col].astype(int).tolist() == saved_headline[col].astype(int).tolist(), f"Mismatch in {col}"
for col in ("cutoff", "catch_rate_percent", "false_alarm_rate_percent", "accuracy_percent"):
    rec_vals = recomputed[col].to_numpy(dtype=float)
    sav_vals = saved_headline[col].to_numpy(dtype=float)
    assert np.array_equal(np.isnan(rec_vals), np.isnan(sav_vals)), f"NaN mismatch in {col}"
    mask = ~np.isnan(rec_vals)
    assert np.allclose(rec_vals[mask], sav_vals[mask], rtol=0.0, atol=1e-9), f"Float mismatch in {col}"

pd.testing.assert_frame_equal(recomputed, saved_headline, check_dtype=False, rtol=0.0, atol=1e-9)
print("\\nASSERT PASSED: recomputed headline table matches saved outputs/results/headline.csv.")"""

    # Cell 11 — Exact result_sentences.txt word for word, Human Decision 2 bootstrap interval next to S3, and Cut List item 2 note
    s1_line = next(line for line in result_sentences if line.startswith("S1:"))
    s3_line = next(line for line in result_sentences if line.startswith("S3:"))
    s6_line = next(line for line in result_sentences if line.startswith("S6:"))

    eb_point = float(bootstrap["esm_minus_blosum"]["point"])
    eb_low = float(bootstrap["esm_minus_blosum"]["low"])
    eb_high = float(bootstrap["esm_minus_blosum"]["high"])

    a3_point = float(bootstrap["arm3_minus_esm"]["point"])
    a3_low = float(bootstrap["arm3_minus_esm"]["low"])
    a3_high = float(bootstrap["arm3_minus_esm"]["high"])
    n_boot_pos = int(bootstrap["n_positions"])
    n_boot_resamples = int(bootstrap["n_resamples_used"])
    ci_pct = int(round(config.BOOTSTRAP_PERCENTILES[1] - config.BOOTSTRAP_PERCENTILES[0]))

    auroc_arm1 = float(auroc_by_arm["arm1"]["auroc"])
    auroc_arm2 = float(auroc_by_arm["arm2"]["auroc"])
    auroc_arm3 = float(auroc_by_arm["arm3"]["auroc"])

    bc_lof = band_by_class["LOF"]
    bc_func = band_by_class["FUNC"]
    bc_int = band_by_class["INT"]
    int_start_arm2 = int_by_key[("start", "arm2")]
    int_end_arm2 = int_by_key[("end", "arm2")]

    cell_11_md = f"""## Pre-committed result sentences (`outputs/results/result_sentences.txt`)

Exact contents of `outputs/results/result_sentences.txt`, word for word:

- `{s1_line}` (Across `{n_boot_resamples:,}` position-cluster bootstrap resamples of the `{n_boot_pos}` end-region positions in `outputs/results/bootstrap.json`, Arm 2's AUROC is `{auroc_arm2:.4f}` versus Arm 1's `{auroc_arm1:.4f}`, giving an `esm_minus_blosum` difference of `{eb_point:.3f}` with a `{ci_pct}%` interval from `{eb_low:.3f}` to `{eb_high:.3f}`.)
- `{s3_line}` In `outputs/results/bootstrap.json`, the `arm3_minus_esm` `{ci_pct}%` position-bootstrap AUROC difference interval is `{a3_low:.3f}` to `{a3_high:.3f}` (point estimate `{a3_point:.3f}`, comparing Arm 3 AUROC `{auroc_arm3:.4f}` with Arm 2 AUROC `{auroc_arm2:.4f}`), so the entire interval lies below zero.
- `{s6_line}`

Secondary breakdowns in Cell 12 also show Arm 2's three-verdict band counts on the `end` region (`outputs/results/band_counts.csv`: `{bc_lof["looks_harmful"]}` `"Looks harmful"`, `{bc_lof["cant_tell"]}` `"Can't tell"`, `{bc_lof["looks_harmless"]}` `"Looks harmless"` among `{bc_lof["total"]}` `LOF` rows; `{bc_func["looks_harmful"]}` `"Looks harmful"`, `{bc_func["cant_tell"]}` `"Can't tell"`, `{bc_func["looks_harmless"]}` `"Looks harmless"` among `{bc_func["total"]}` `FUNC` rows; `{bc_int["looks_harmful"]}`, `{bc_int["cant_tell"]}`, `{bc_int["looks_harmless"]}` among `{bc_int["total"]}` `INT` rows) and intermediate-row flagging (`outputs/results/int_rows.csv`: Arm 2 flagged `{int_start_arm2["flagged"]}` of `{int_start_arm2["n_int"]}` or `{int_start_arm2["flagged_percent"]}%` in `start` and `{int_end_arm2["flagged"]}` of `{int_end_arm2["n_int"]}` or `{int_end_arm2["flagged_percent"]}%` in `end`). The secondary reverse split (`reverse_split.csv`, tuning on the `end` region and testing on the `start` region) was cut prior to the test run per `writeup/decision_log.md` (Cut List item 2)."""

    # Cell 12 — Display auroc.csv, bootstrap.json, by_swap.csv, int_rows.csv, and band_counts.csv
    cell_12_code = """for filename in ("auroc.csv", "by_swap.csv", "int_rows.csv", "band_counts.csv"):
    df_part = pd.read_csv(config.RESULTS_DIR / filename)
    print(f"=== outputs/results/{filename} ===")
    print(df_part.to_string(index=False))
    print()

with open(config.RESULTS_DIR / "bootstrap.json", "r", encoding="utf-8") as handle:
    bootstrap_data = json.load(handle)
print("=== outputs/results/bootstrap.json ===")
print(json.dumps(bootstrap_data, indent=2))"""

    # Cell 13 — Where it fails (Human Decision 4: 3 of 5 worst misses at 1823 and 1824, A1823S explained wrong case, all causal explanations prefixed "Guess:")
    miss_rows = [r for r in failures_rows if r["kind"] == "miss"]
    fa_rows = [r for r in failures_rows if r["kind"] == "false_alarm"]
    top3_positions = sorted({int(r["position"]) for r in miss_rows[:3]})
    misses_cluster = [
        r for r in miss_rows if int(r["position"]) in top3_positions
    ]
    n_misses_total = len(miss_rows)
    n_fa_total = len(fa_rows)
    n_misses_cluster = len(misses_cluster)
    cluster_mutants_str = ", ".join(f"`{r['mutant']}`" for r in misses_cluster)

    a1823s = miss_rows[0]
    a1823s_mut = a1823s["mutant"]
    a1823s_pos = int(a1823s["position"])
    a1823s_func = float(a1823s["function_score"])
    a1823s_esm = float(a1823s["esm_score"])
    a1823s_harm2 = float(a1823s["harm_arm2"])
    a1823s_blosum = int(a1823s["blosum_score"])
    a1823s_ent = float(a1823s["entropy"])
    a1823s_guess = a1823s["guess_why"]
    if not a1823s_guess.startswith("Guess:"):
        raise ValueError("Expected guess_why for A1823S to start with 'Guess:'.")

    cell_13_md = f"""## Where it fails (`outputs/results/failures.csv`)

`outputs/results/failures.csv` lists Arm 2's `{n_misses_total}` worst misses (`findlay_class == "LOF"` and `flag_arm2 == 0`, sorted by `harm_arm2` ascending) and `{n_fa_total}` worst false alarms (`findlay_class == "FUNC"` and `flag_arm2 == 1`, sorted by `harm_arm2` descending) on the held-out `end` region.

- **Fact from `outputs/results/failures.csv`:** `{n_misses_cluster}` of Arm 2's `{n_misses_total}` worst misses sit at adjacent positions `{top3_positions[0]}` and `{top3_positions[1]}` ({cluster_mutants_str}).
- **Explained wrong case (`{a1823s_mut}`, rank `{a1823s["rank"]}` miss at position `{a1823s_pos}`):** In the lab assay, `{a1823s_mut}` caused loss of function (`function_score = {a1823s_func:.6f}`, `findlay_class = LOF`), yet it was left unflagged by all three arms (`flag_arm1 = {a1823s["flag_arm1"]}`, `flag_arm2 = 0`, `flag_arm3 = {a1823s["flag_arm3"]}`, `band_arm2 = {a1823s["band_arm2"]}`). {a1823s_guess}
- **Splicing note:** No supplementary full DNA/RNA study table is present in `data/`, so DNA-level splicing disruption remains an untested limitation, and every causal explanation in `outputs/results/failures.csv` begins with `Guess:` (`blosum_score = {a1823s_blosum}`, `entropy = {a1823s_ent:.6f}`, `esm_score = {a1823s_esm:.6f}`, `harm_arm2 = {a1823s_harm2:.6f}`)."""

    # Cell 14 — Load and display outputs/results/failures.csv
    cell_14_code = """failures_df = pd.read_csv(config.RESULTS_DIR / "failures.csv")
assert all(str(g).startswith("Guess:") for g in failures_df["guess_why"])
print(failures_df.to_string(index=False))"""

    # Cell 15 — The heat map and how to read it (Human Decision 3: only state what heatmap.png shows at the 8 PDB 1JM7 zinc positions)
    zinc_positions = list(config.ZINC_POSITIONS_TO_CONFIRM)
    n_zinc_pos = len(zinc_positions)
    zinc_scores_df = scores_df.loc[scores_df["position"].isin(zinc_positions)]
    n_zinc_swaps = int(len(zinc_scores_df))
    n_zinc_neg = int((zinc_scores_df["esm_score"] < 0.0).sum())
    zinc_min_esm = float(zinc_scores_df["esm_score"].min())
    zinc_max_esm = float(zinc_scores_df["esm_score"].max())
    zinc_min_row = zinc_scores_df.sort_values("esm_score", ascending=True).iloc[0]
    zinc_max_row = zinc_scores_df.sort_values("esm_score", ascending=False).iloc[0]

    zinc_labels_df = start_labels_df.loc[start_labels_df["position"].isin(zinc_positions)]
    n_zinc_tested = int(len(zinc_labels_df))
    n_zinc_lof = int((zinc_labels_df["findlay_class"] == "LOF").sum())
    n_zinc_func = int((zinc_labels_df["findlay_class"] == "FUNC").sum())
    n_zinc_int = int((zinc_labels_df["findlay_class"] == "INT").sum())
    n_zinc_grey = n_zinc_func + n_zinc_int

    zinc_wt_pairs = [
        f"`{scores_df.loc[scores_df['position'] == p, 'wild_type'].iloc[0]}{p}`"
        for p in zinc_positions
    ]
    zinc_pos_list_str = ", ".join(str(p) for p in zinc_positions)

    cell_15_md = f"""## The heat map (`outputs/heatmap.png`) and how to read it

`outputs/heatmap.png` renders all `{config.N_SCORE_ROWS:,}` ESM-2 scores across the `{config.N_POSITIONS}` scored positions in two stacked panels (`start` region `{start_lo}–{start_hi}` on top; `end` region `{end_lo}–{end_hi}` on the bottom), with the `{n_aa}` standard amino acids (`{config.AA_ORDER}`, `A` at the top to `Y` at the bottom) on the vertical axis:
- **Colour scale (`{config.HEATMAP_COLORSCALE}`, `{config.HEATMAP_ZMIN:.0f}` to `+{config.HEATMAP_ZMAX:.0f}`, centred at `0`):** Red cells (`esm_score < 0`) indicate swaps the model finds less likely than the normal letter; blue cells (`esm_score > 0`) indicate swaps the model finds more likely than the normal letter.
- **Blank white cells (`1` per column):** The cell where the row letter equals that position's wild-type amino acid has no mutation score (`NaN`) and is drawn in white.
- **Lab-tested overlay dots (`{tot_rows:,}` total):** `{tot_lof}` black dots mark lab-confirmed `LOF` (harmful) swaps (`{s_lof}` in `start`, `{e_lof}` in `end`), and `{tot_non_lof:,}` grey dots mark lab-tested non-`LOF` (`FUNC` or `INT`) swaps (`{s_non_lof}` in `start`, `{e_non_lof}` in `end`). (The interactive HTML heat map `outputs/heatmap.html` was cut per Cut List item 4 in `writeup/decision_log.md`.)
- **What `outputs/heatmap.png` shows at the `{n_zinc_pos}` PDB `1JM7` zinc-holding positions (`[{zinc_pos_list_str}]`, wild-type residues {", ".join(zinc_wt_pairs)}):** Above the `start` panel, downward triangles mark those `{n_zinc_pos}` columns (`"zinc-holding letters (PDB 1JM7)"`). Across the `{n_zinc_swaps}` scored non-wild-type cells in those `{n_zinc_pos}` columns (`{n_zinc_pos} × {n_swaps_per_pos}`), `{n_zinc_neg}` of `{n_zinc_swaps}` cells have negative ESM-2 scores, ranging from `{zinc_min_esm:.6f}` (`{zinc_min_row["mutant"]}`) to `{zinc_max_esm:.6f}` (`{zinc_max_row["mutant"]}`). Of the `{n_zinc_tested}` lab-tested swaps in those `{n_zinc_pos}` columns, `{n_zinc_lof}` carry black `LOF` dots and `{n_zinc_grey}` carry grey non-`LOF` dots (`{n_zinc_int}` `INT` and `{n_zinc_func}` `FUNC`). We report only the scores and markers at those `{n_zinc_pos}` marked columns and do not claim that those `{n_zinc_pos}` sites explain the colours across the surrounding span of positions `{zinc_positions[0]}–{zinc_positions[-1]}`."""

    # Cell 16 — Display outputs/heatmap.png
    cell_16_code = """from IPython.display import Image, display

assert config.HEATMAP_PNG_PATH.is_file(), f"Missing {config.HEATMAP_PNG_PATH}"
display(Image(filename="outputs/heatmap.png"))"""

    # Cell 17 — Limitations, groups check, sources and credits, and final disclaimer
    swap_lookup = {
        (row["arm"], row["swap_kind"]): row for row in by_swap_rows
    }
    h_arm0 = headline_by_arm["arm0"]
    h_arm1 = headline_by_arm["arm1"]
    h_arm2 = headline_by_arm["arm2"]
    h_arm3 = headline_by_arm["arm3"]

    s_arm2_caught = int(frozen["arms"]["arm2"]["tuning_caught"])
    s_arm2_fa = int(frozen["arms"]["arm2"]["tuning_false_alarms"])
    s_arm2_catch_pct = 100.0 * s_arm2_caught / s_lof
    s_arm2_fa_pct = 100.0 * s_arm2_fa / s_func

    tp_2 = swap_lookup[("arm2", "to_proline")]
    tp_1 = swap_lookup[("arm1", "to_proline")]
    tp_3 = swap_lookup[("arm3", "to_proline")]

    lc_2 = swap_lookup[("arm2", "loss_of_cysteine")]
    lc_1 = swap_lookup[("arm1", "loss_of_cysteine")]
    lc_3 = swap_lookup[("arm3", "loss_of_cysteine")]

    cc_2 = swap_lookup[("arm2", "charge_change")]
    cc_1 = swap_lookup[("arm1", "charge_change")]
    cc_3 = swap_lookup[("arm3", "charge_change")]

    ot_2 = swap_lookup[("arm2", "other")]
    ot_1 = swap_lookup[("arm1", "other")]
    ot_3 = swap_lookup[("arm3", "other")]

    cell_17_md = f"""## Limitations, groups check, sources, and credits

### Limitations

**Main three limitations:**
1. **Blind to splicing (we could not test this):** The {spec_facts["findlay_year"]} lab study edited genomic DNA (`{spec_facts["n_dna_changes"]}` single-nucleotide changes), so some mutations are harmful because they disrupt RNA splicing rather than because the amino-acid replacement destabilises the protein. A protein-sequence model cannot see DNA-level splicing effects, and because the supplementary DNA/RNA table is not in `data/`, we could not test this.
2. **Cells in a dish, not people:** The lab measurements track cell survival in cultured human cells, which is direct cellular evidence but indirect evidence about cancer risk in people. A {spec_facts["preprint_year"]} follow-up preprint from the same research group reports that assay results can also depend on the cell type used.
3. **Stronger tools exist and the method is not new:** Stronger variant effect predictors exist, such as AlphaMissense (Cheng et al.), whose predictions are free to download as far as we know; we do not claim to be better or cheaper. Zero-shot masked-marginal scoring with ESM-2 on this BRCA1 table is already part of the ProteinGym benchmark, and the public Evo 2 notebook evaluates a zero-shot DNA model on the same study.

**Additional five limitations:**
4. **Small coverage:** We evaluate `{config.N_POSITIONS}` of BRCA1's `{config.PROTEIN_LENGTH:,}` positions (`{start_lo}–{start_hi}` and `{end_lo}–{end_hi}`), single-letter amino-acid substitutions only, in one gene.
5. **The model has read the BRCA1 sequence during pre-training:** `{config.MODEL_ID}` saw natural BRCA1 protein sequences in UniRef50 during training (though never any lab function scores or labels).
6. **The cutoffs are our choice:** Targeting `{config.CATCH_TARGET_PERCENT}%` harmful catch on the `start` region (and `{config.FALSE_ALARM_TARGET_PERCENT}%` false alarms for Cutoff B) is a design choice; a different cutoff rule changes the balance between missed harmful mutations and false alarms (`{h_arm2["false_alarms"]}` of `{h_arm2["n_harmless"]}` harmless mutations flagged by Arm 2 on the `end` region).
7. **Cross-region cutoff drift and non-transfer of supervised weights:** Cutoffs tuned on the `start` region (`{start_lo}–{start_hi}`) shifted when applied unchanged to the `end` region (`{end_lo}–{end_hi}`), catching `{h_arm1["catch_rate_percent"]}%` (`{h_arm1["harmful_caught"]}/{h_arm1["n_harmful"]}`) for Arm 1, `{h_arm2["catch_rate_percent"]}%` (`{h_arm2["harmful_caught"]}/{h_arm2["n_harmful"]}`) for Arm 2, and `{h_arm3["catch_rate_percent"]}%` (`{h_arm3["harmful_caught"]}/{h_arm3["n_harmful"]}`) for Arm 3 (while Arm 0 caught `{h_arm0["harmful_caught"]}/{h_arm0["n_harmful"]}` with `{h_arm0["accuracy_percent"]}%` accuracy).
8. **Redundant where tested:** At these `{config.N_POSITIONS}` positions, experimental lab measurements already exist (`{spec_facts["dna_coverage_pct"]}` of all possible single-letter DNA changes measured).

### Groups check

Results are broken down across protein regions and across the `{len(config.SWAP_KINDS)}` swap kinds in `outputs/results/by_swap.csv`:
- **By protein region:** On the `start` region (`{start_lo}–{start_hi}`, in-sample tuning), Arm 2 flagged `{s_arm2_caught}` of `{s_lof}` `LOF` mutations (`{s_arm2_catch_pct:.1f}%`) and `{s_arm2_fa}` of `{s_func}` `FUNC` mutations (`{s_arm2_fa_pct:.1f}%`); on the held-out `end` region (`{end_lo}–{end_hi}`), Arm 2 caught `{h_arm2["harmful_caught"]}` of `{h_arm2["n_harmful"]}` `LOF` mutations (`{h_arm2["catch_rate_percent"]}%`) and flagged `{h_arm2["false_alarms"]}` of `{h_arm2["n_harmless"]}` `FUNC` mutations (`{h_arm2["false_alarm_rate_percent"]}%`). (The secondary reverse split was cut per Cut List item 2.)
- **By swap kind on the held-out `end` region (`outputs/results/by_swap.csv`):**
  - `to_proline` (`{tp_2["n_harmful"]}` `LOF`, `{tp_2["n_harmless"]}` `FUNC`): Arm 2 caught `{tp_2["harmful_caught"]}/{tp_2["n_harmful"]}` (`{tp_2["catch_rate_percent"]}%`) with `{tp_2["false_alarms"]}/{tp_2["n_harmless"]}` false alarms (`{tp_2["false_alarm_rate_percent"]}%`); Arm 1 caught `{tp_1["harmful_caught"]}/{tp_1["n_harmful"]}` (`{tp_1["catch_rate_percent"]}%`) with `{tp_1["false_alarms"]}/{tp_1["n_harmless"]}` false alarms (`{tp_1["false_alarm_rate_percent"]}%`); Arm 3 caught `{tp_3["harmful_caught"]}/{tp_3["n_harmful"]}` (`{tp_3["catch_rate_percent"]}%`) with `{tp_3["false_alarms"]}/{tp_3["n_harmless"]}` false alarms (`{tp_3["false_alarm_rate_percent"]}%`).
  - `loss_of_cysteine` (`{lc_2["n_harmful"]}` `LOF`, `{lc_2["n_harmless"]}` `FUNC`): Arm 2 caught `{lc_2["harmful_caught"]}/{lc_2["n_harmful"]}` (`{lc_2["catch_rate_percent"]}%`) with `{lc_2["false_alarms"]}/{lc_2["n_harmless"]}` false alarms (`{lc_2["false_alarm_rate_percent"]}%`); Arm 1 and Arm 3 both caught `{lc_1["harmful_caught"]}/{lc_1["n_harmful"]}` (`{lc_1["catch_rate_percent"]}%`) while flagging all `{lc_1["false_alarms"]}/{lc_1["n_harmless"]}` harmless swaps (`{lc_1["false_alarm_rate_percent"]}%` and `{lc_3["false_alarm_rate_percent"]}%`).
  - `charge_change` (`{cc_2["n_harmful"]}` `LOF`, `{cc_2["n_harmless"]}` `FUNC`): Arm 2 caught `{cc_2["harmful_caught"]}/{cc_2["n_harmful"]}` (`{cc_2["catch_rate_percent"]}%`) with `{cc_2["false_alarms"]}/{cc_2["n_harmless"]}` false alarms (`{cc_2["false_alarm_rate_percent"]}%`), versus Arm 1 (`{cc_1["harmful_caught"]}/{cc_1["n_harmful"]}` or `{cc_1["catch_rate_percent"]}%` caught, `{cc_1["false_alarms"]}/{cc_1["n_harmless"]}` or `{cc_1["false_alarm_rate_percent"]}%` false alarms) and Arm 3 (`{cc_3["harmful_caught"]}/{cc_3["n_harmful"]}` or `{cc_3["catch_rate_percent"]}%` caught, `{cc_3["false_alarms"]}/{cc_3["n_harmless"]}` or `{cc_3["false_alarm_rate_percent"]}%` false alarms).
  - `other` (`{ot_2["n_harmful"]}` `LOF`, `{ot_2["n_harmless"]}` `FUNC`): Arm 2 caught `{ot_2["harmful_caught"]}/{ot_2["n_harmful"]}` (`{ot_2["catch_rate_percent"]}%`) with `{ot_2["false_alarms"]}/{ot_2["n_harmless"]}` false alarms (`{ot_2["false_alarm_rate_percent"]}%`), versus Arm 1 (`{ot_1["harmful_caught"]}/{ot_1["n_harmful"]}` or `{ot_1["catch_rate_percent"]}%` caught, `{ot_1["false_alarms"]}/{ot_1["n_harmless"]}` or `{ot_1["false_alarm_rate_percent"]}%` false alarms) and Arm 3 (`{ot_3["harmful_caught"]}/{ot_3["n_harmful"]}` or `{ot_3["catch_rate_percent"]}%` caught, `{ot_3["false_alarms"]}/{ot_3["n_harmless"]}` or `{ot_3["false_alarm_rate_percent"]}%` false alarms).

Our data are mutations, not people, so we cannot check performance by age, sex or ancestry.

### Sources and credits

{credits_md_block}

*{config.DISCLAIMER}*"""

    cells = [
        new_markdown_cell(cell_1_md),
        new_code_cell(cell_2_code),
        new_markdown_cell(cell_3_md),
        new_code_cell(cell_4_code),
        new_markdown_cell(cell_5_md),
        new_code_cell(cell_6_code),
        new_markdown_cell(cell_7_md),
        new_code_cell(cell_8_code),
        new_markdown_cell(cell_9_md),
        new_code_cell(cell_10_code),
        new_markdown_cell(cell_11_md),
        new_code_cell(cell_12_code),
        new_markdown_cell(cell_13_md),
        new_code_cell(cell_14_code),
        new_markdown_cell(cell_15_md),
        new_code_cell(cell_16_code),
        new_markdown_cell(cell_17_md),
    ]

    nb = new_notebook(
        cells=cells,
        metadata={
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "pygments_lexer": "ipython3",
            },
        },
    )
    return nb


def main() -> int:
    """Build notebook.ipynb at the repository root."""
    nb = build_notebook()
    notebook_path = ROOT / "notebook.ipynb"
    with open(notebook_path, "w", encoding="utf-8", newline="\n") as handle:
        nbformat.write(nb, handle)
    print(f"Wrote {notebook_path} with {len(nb.cells)} cells.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
