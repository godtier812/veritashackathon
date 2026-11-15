# BRCA Spellcheck: an AI second opinion on uncertain BRCA1 mutations

**Track:** Health
**Team:** [names — the humans fill this in]

*This project is a prototype and is not medical advice.*

A model that only read protein sequences ranked harmful BRCA1 mutations better than a classic substitution table. To catch 94.4% of harmful mutations, the model also flagged 437 of 899 harmless ones.

## The problem

When a patient's genetic test returns a "variant of uncertain significance" (VUS), a genetic counselor has to sit with that patient while no one can yet call the mutation harmful or harmless. Most genes have no experimental lab map of which mutations break the protein, so counselors often have to lean on computer predictions and need to know how far to trust them and where they go blind.

## What we built

BRCA Spellcheck takes a single-letter BRCA1 protein substitution such as `C61G` and returns a sequence-model harm score and verdict alongside the experimental CRISPR lab result. In this project, BRCA1 is the exam rather than the use case: at the positions we test, Findlay et al. (2018) already measured 96.5% of all possible single-letter DNA changes, so a counselor would use the lab map for BRCA1 while using this exam to judge how far a sequence-only model can be trusted on genes without a lab map. Stage 5 held-out evaluation is complete; the failure-analysis, demo app, and notebook stages follow next.

## What we added

Our own idea. Grading a language model zero-shot against these BRCA1 lab measurements has been done before in public benchmarks (ProteinGym and the Arc Institute / NVIDIA Evo 2 BRCA1 notebook), and stronger predictors such as AlphaMissense also exist. On top of prior zero-shot scoring, we hold out an entire protein region with decision cutoffs frozen in advance, test whether supervised labels from the start region transfer to the held-out end region, count missed harmful mutations and false alarms separately at a fixed 90% catch target, and inspect specific failure cases where the sequence model goes wrong.

## How it works

We use the Findlay et al. (2018) saturation genome editing dataset via ProteinGym, covering 1,837 lab-tested single-letter missense mutations across the start region (positions 1–101) and end region (positions 1631–1855), with 155 intermediate (`INT`) rows excluded from tuning and headline rates (`131` harmful `LOF` and `385` harmless `FUNC` in the start region; `267` `LOF` and `899` `FUNC` in the end region). Each mutation is scored zero-shot using masked-marginal ("hide-a-letter") log-probability differences from ESM-2 (`facebook/esm2_t33_650M_UR50D`), which saw the natural BRCA1 protein sequence during pretraining but never saw any patient records or lab labels. We compare ESM-2 zero-shot (`arm2`) against an always-harmless baseline (`arm0`), the 1992 `BLOSUM62` substitution table (`arm1`), and a start-region logistic regression combining ESM-2, `BLOSUM62`, and amino-acid chemistry features (`arm3`), freezing all 90%-catch cutoffs on positions 1–101 (`outputs/frozen.json`) before running a single evaluation (`run=1`) on held-out positions 1631–1855.

## Results

The held-out evaluation was executed once (`run=1` in `outputs/results/run_log.txt`) on the 1,166 end-region `LOF` and `FUNC` variants (267 harmful, 899 harmless), reporting missed harmful cases and false alarms separately alongside cutoff-free AUROC (`0.5000` for `arm0`, `0.6681` for `arm1`, `0.8995` for `arm2`, and `0.8662` for `arm3`). What didn't work: adding supervised start-region labels in `arm3` reduced held-out end-region AUROC from `0.8995` (zero-shot `arm2`) to `0.8662` (position-bootstrap difference `-0.0333`, 95% CI `[-0.0535, -0.0144]`) and missed 22 harmful variants compared with 15 for `arm2`, while the tied `BLOSUM62` rule (`arm1`) flagged 571 of 899 harmless variants (63.5%) as false alarms.

| arm | arm_name | cutoff | n_harmful | n_harmless | harmful_caught | harmful_missed | false_alarms | harmless_cleared | catch_rate_percent | false_alarm_rate_percent | accuracy_percent |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| arm0 | Always harmless |  | 267 | 899 | 0 | 267 | 0 | 899 | 0.0 | 0.0 | 77.1 |
| arm1 | BLOSUM62 rule | 0.000000 | 267 | 899 | 232 | 35 | 571 | 328 | 86.9 | 63.5 | 48.0 |
| arm2 | ESM-2 zero-shot | 3.034498 | 267 | 899 | 252 | 15 | 437 | 462 | 94.4 | 48.6 | 61.2 |
| arm3 | Logistic regression on ESM-2 + BLOSUM62 + chemistry | 0.099752 | 267 | 899 | 245 | 22 | 405 | 494 | 91.8 | 45.1 | 63.4 |

```text
S1: A model that only read protein sequences ranked harmful BRCA1 mutations better than a classic substitution table.
S3: We could not detect a benefit from adding labels from the other end of the protein.
S6: To catch 94.4% of harmful mutations, the model also flagged 437 of 899 harmless ones.
```

## What we'd do next

First, we would pair the protein sequence model with a DNA-level splicing predictor to cover mutations that disrupt RNA splicing outside protein-level visibility. Second, we would repeat the same frozen-cutoff, held-out-region exam on a second gene with an experimental lab map to test how these miss and false-alarm rates behave on another protein.

## Sources

Experimental BRCA1 saturation genome editing data and class cutoffs come from Findlay et al. (2018, *Nature*, University of Washington SGE resource) as distributed by ProteinGym (Notin et al., 2023), and mutation scoring uses ESM-2 (`facebook/esm2_t33_650M_UR50D`, Lin et al., 2023) with masked-marginal scoring (Meier et al., 2021), the BLOSUM62 matrix (Henikoff & Henikoff, 1992) via Biopython, hydrophobicity values from Kyte & Doolittle (1982), and residue volumes from Zamyatnin (1972). We also cite prior zero-shot and variant-effect work on this dataset and gene—ProteinGym, the Arc Institute / NVIDIA Evo 2 BRCA1 notebook, AlphaMissense (Cheng et al., 2023), and the 2025 follow-up preprint on cell-type context—alongside the Python software stack (PyTorch, Hugging Face Transformers, NumPy, pandas, SciPy, scikit-learn, Biopython, Plotly, Matplotlib, Gradio, pytest, and Jupyter). AI tools used: Claude Code (Anthropic) for planning, writing the build specification, and review, and Gemini in Google Antigravity for writing and running the project code; full links and citations are listed in [`writeup/credits.md`](writeup/credits.md).
