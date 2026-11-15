# Decision Log

## Stage 1 — Setup and Data Foundations

- **Label rule:** Classify `DMS_score < -1.328` as harmful (`LOF`), `DMS_score > -0.748` as harmless (`FUNC`), and exclude intermediate (`INT`, `-1.328 <= DMS_score <= -0.748`) variants from tuning and headline evaluation, matching Findlay et al. (2018).
- **Region split:** Tune exclusively on the RING-domain start region (positions 1–101) and test once on the BRCT-domain end region (positions 1631–1855), never splitting rows at random across positions.
- **Reproducibility seed:** Fix seed `42` across `random`, `numpy`, `torch`, `scikit-learn`, and bootstrap resampling (`numpy.random.default_rng(42)`).
- **Hydrophobicity scale:** Use Kyte & Doolittle (1982, *J. Mol. Biol.* 157:105–132) hydropathy values for all amino acids.
- **Residue volume scale:** Use Zamyatnin (1972, *Prog. Biophys. Mol. Biol.* 24:107–123) amino acid residue volumes in Å³.
- **Histidine charge:** Count histidine (`H`) as neutral (`0`) at intracellular pH when computing charge changes.
- **Signed biochemical deltas:** Compute `size_change` and `hydrophobicity_change` as signed `new - wild_type` differences rather than unsigned magnitudes so directionality is preserved.
- **Cross-validation scheme:** Use `StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)` grouped by `position` so mutations at the same amino acid site never leak across folds.
- **Region derivation:** Derive `region` from `position` (`1–101` as `start`, `1631–1855` as `end`) because the raw ProteinGym table (`BRCA1_HUMAN_Findlay_2018.csv`) has no region column of its own.

## Stage 2 — Baselines and Cutoff Logic

- **Inclusive highest-cutoff rule (`>=`):** Flag mutations at or above the highest observed harm score achieving at least 90% harmful catch rate (using exact integer arithmetic `-(-(target * n) // 100)`), ensuring identical harm scores always receive identical calls.
- **BLOSUM62 integer ties (Rule 14):** Because BLOSUM62 scores are discrete integers, 8 tied harmful rows at `harm_arm1 = 0.0` push Arm 1 start-region catch from 112 to 120 of 131 (91.6%) and flag 337 of 516 rows (65.3%) with 217 of 385 false alarms (56.4%); keep exact ties with no noise or tie-breakers added.
- **Sign formatting (`0.0 - x`):** Compute negated harm scores via `0.0 - score` so zero values serialize as positive `0.0` rather than `-0.0`.

## Stage 3 — ESM-2 Zero-Shot Scoring and Schedule Cuts

- **ESM-2 scoring configuration:** Scored all 326 positions (6,194 single-amino-acid substitutions) with `facebook/esm2_t33_650M_UR50D` in `float32` via masked-marginal log-odds (`log P(new) - log P(wild_type)`) over two 500-residue windows (`A`: 1–500, `B`: 1364–1863), running on CPU (`--allow-cpu`) with human approval (419 seconds, float32) because this build machine has no NVIDIA GPU.
- **Cut list item 2 (human decision):** Cut the secondary reverse split (tuned on `end`, tested on `start`; `BUILD_SPEC.md` Section 14, item 2) per human decision.
- **Cut list item 4 (human decision):** Cut the interactive Plotly HTML heat map (`outputs/heatmap.html`) and keep the static Matplotlib PNG heat map (`outputs/heatmap.png`; `BUILD_SPEC.md` Section 14, item 4) per human decision.

## Checkpoint A — Human Decisions for Stages 4–12

1. **CPU Stage 3 approval:** Full Stage 3 ESM-2 650M run on CPU with `--allow-cpu` approved by a human (`419` seconds, `float32`).
2. **Spec anonymization pause:** Replaced every team personal name in `BUILD_SPEC.md` with `"the operator"` or `"the reviewer"` for Round 1 blind judging, staged with `git add -A` without committing, and paused for `"continue"`.
3. **Uninterrupted Stage 4:** After `"continue"`, run Stage 4 in full without further stops through Checkpoint B.
4. **Local test runner substitution:** Do not wait for external `test-runner` after any stage; local `pytest` runs stand in for it.
5. **Checkpoint B local audit substitution:** At Checkpoint B, run a fresh read-only audit of `src/` and `tests/` against all 15 locked rules in Section 3 and paste its full report in place of external `leak-auditor`.
6. **Stage 5b version-zero README:** When `"run the test"` is given, complete Stage 5, then immediately write and commit a version-zero `README.md` (`stage 5b: version zero README`) before stopping at Checkpoint C: Section 12.1 headings, one to three sentences each, exact disclaimer line, headline table copied from `outputs/results/headline.md`, and `outputs/results/result_sentences.txt` lines word for word.
7. **Chained Stages 6 and 7:** When `"go C"` is given, run Stages 6 and 7 back-to-back without stopping between them.
8. **Zinc-coordinating sites confirmed (PDB 1JM7):** Confirmed Section 7.10 zinc-binding RING positions `[24, 27, 39, 41, 44, 47, 61, 64]` against BRCA1-BARD1 RING structure PDB entry `1JM7`; draw heat map markers with the label `"zinc-holding letters (PDB 1JM7)"`.
9. **Updated fixed sentence S3 wording:** Sentence `S3` is now: `"We could not detect a benefit from adding labels from the other end of the protein."`
10. **Write-up takeaway line placement:** Put one takeaway line directly under the disclaimer, made of the main fixed sentence followed by `S6`.

## Stage 4 — Feature Engineering, Start-Region Freezing, and Evaluation Guards

- **Arm 3 feature set and fixed order:** Combine `esm_score`, `blosum_score`, `entropy`, `size_change` (Zamyatnin 1972 residue volume in Å³, `new - wild_type`), `hydrophobicity_change` (Kyte & Doolittle 1982 hydropathy, `new - wild_type`), `to_proline`, `from_cysteine`, and `from_glycine` in fixed `FEATURE_ORDER`.
- **Swap-kind priority hierarchy:** Assign each mutation to one swap kind in strict priority order: `to_proline` (`new_letter == "P"`), `loss_of_cysteine` (`wild_type == "C"`), `charge_change` (net charge change with `K, R = +1`, `D, E = -1`, and `H = 0`), then `other`.
- **Arm 3 out-of-fold cutoff and refit (Rules 4, 7, 13):** Choose the Arm 3 90%-catch probability cutoff on out-of-fold probabilities from `StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)` grouped by `position` over the 516 start-region `LOF` and `FUNC` rows (`INT` excluded) using a `StandardScaler` + default `LogisticRegression(random_state=42)` pipeline, then refit on all 516 start-region rows and freeze the fitted parameters in `outputs/frozen.json`.
- **Arm 2 "Can't tell" band (Rule 10):** Set Cutoff A to Arm 2's 90%-catch start-region cutoff and Cutoff B to the lowest observed harm score flagging at most 10% of the 385 start-region `FUNC` rows (at most 38 rows), enabling a `"Can't tell"` band only if `cutoff_b > cutoff_a`.
- **Pre-test freeze and run-log locks (Rules 8, 9, 15):** Record the SHA-256 hashes of `outputs/esm_scores.csv` and `data/brca1_missense_scores_slim.csv` inside `outputs/frozen.json`, block overwriting or score recomputation once frozen or once `outputs/results/run_log.txt` records any attempt, and require `status=started` to be written to `run_log.txt` immediately before end-region labels are loaded in `src/evaluate.py`.

## Stage 5 — End-Region Test Run

- **Single frozen end-region evaluation (Rules 2, 5, 6, 8, 13, 15):** Executed `python -m src.evaluate` once (`run=1`, `kind=first`, `status=started` followed by `status=success` in `outputs/results/run_log.txt`) on the held-out BRCT end region (`267` `LOF`, `899` `FUNC`, `96` `INT`) using only the pre-frozen parameters from `outputs/frozen.json`, omitting `reverse_split.csv` per Cut List item 2 and writing `headline.csv`, `headline.md`, `auroc.csv`, `bootstrap.json`, `band_counts.csv`, `by_swap.csv`, `int_rows.csv`, `predictions_end.csv`, `sanity7.json`, and `result_sentences.txt`.

## Stage 6 — Failure Analysis

- **Worst misses and false alarms selection:** Extracted Arm 2's five worst misses (`findlay_class == "LOF"`, `flag_arm2 == 0`, sorted by `harm_arm2` ascending, then `function_score` ascending, then `mutant` ascending) and five worst false alarms (`findlay_class == "FUNC"`, `flag_arm2 == 1`, sorted by `harm_arm2` descending, then `function_score` descending, then `mutant` ascending) from `outputs/results/predictions_end.csv` into `outputs/results/failures.csv` and `outputs/results/failures.md`.
- **Labelled guesses and explained wrong case:** Prefixed every `guess_why` entry with `Guess:`, grounded each hypothesis strictly in observed table columns (`wild_type`, `new_letter`, `swap_kind`, `blosum_score`, `entropy`, `esm_score`, `harm_arm2`, `band_arm2`, `flag_arm1`, `flag_arm3`) with no biological mechanism stated as fact and no zinc references on BRCT positions, and selected rank-1 miss `A1823S` as the explained wrong case for the write-up.
- **Splicing table check:** Confirmed no supplementary full study table is present in `data/`, recording DNA-level RNA splicing disruption as an untested limitation of protein-sequence-only scoring.

## Stage 7 — Static ESM-2 Heat Map

- **Static Matplotlib heat map (`outputs/heatmap.png`):** Rendered both the start region (`1–101`) and end region (`1631–1855`) on a fixed diverging `RdBu` scale (`vmin = -12.0`, `vmax = 12.0` centred at `0.0`) with white blank cells at each wild-type residue, overlaying `398` black `LOF` dots and `1,439` grey non-`LOF` (`FUNC`/`INT`) dots while omitting `outputs/heatmap.html` per Cut List item 4.
- **RING zinc-coordinating site markers (Human Decision 8):** Marked RING positions `[24, 27, 39, 41, 44, 47, 61, 64]` above the start-region panel with the legend label `"zinc-holding letters (PDB 1JM7)"`.

## Stage 8 — Gradio Prototype Demo

- **Self-contained app export (`src/export_app.py`):** Exported `app/app_data.csv` (`6,194` rows, `1,837` rows with `lab_class` and 3-decimal `lab_score`) and `app/app_config.json` (`cutoff_a = 3.034498`, `cutoff_b = 7.751908`, `has_band = true`, `end_labels_included = true`) using a self-contained SHA-256 check on `outputs/esm_scores.csv` without importing `src.evaluate`, `Bio`, `torch`, or `transformers`.
- **Standalone Gradio app (`app/app.py`):** Built the Hugging Face Space demo importing only `gradio`, `pandas`, `json`, `re`, and `pathlib` (no `src`, `torch`, or `transformers` imports), reading only `app/app_data.csv` and `app/app_config.json`, enforcing exact R1–R8 refusal messages, appending both fixed notes (`DEMO_NOTE` and `DISCLAIMER`) to every response and refusal, and calling `demo.launch()` with no arguments.
- **Space packaging and Stage 8 unit tests:** Wrote `app/requirements.txt` (`pandas==2.3.3`) and `app/README.md` (`sdk: gradio`, `sdk_version: 6.32.0`, exam-not-use-case framing, and exact disclaimer), and verified all six Stage 8 unit tests (`test_app_parser_good`, `test_app_refusals`, `test_app_notes_always_present`, `test_app_lab_display`, `test_app_end_region_locked`, `test_app_imports_no_model`).

