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
