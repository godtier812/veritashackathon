# BRCA Spellcheck: an AI second opinion on uncertain BRCA1 mutations

**Track:** Health
**Team:** [names — the humans fill this in]

*This project is a prototype and is not medical advice.*

A model that only read protein sequences ranked harmful BRCA1 mutations better than a classic substitution table. To catch 94.4% of harmful mutations, the model also flagged 437 of 899 harmless ones. <!-- result_sentences.txt S1 + S6 -->

## The problem

A genetic counselor sitting with a patient whose test returned a "variant of uncertain significance" (VUS) faces a mutation nobody can yet call harmful or harmless. Most genes have no lab map, so counselors lean on computer predictions and need to know how far to trust them.

## What we built

BRCA Spellcheck takes a single-letter BRCA1 mutation (`C61G`) and outputs a verdict (`Looks harmful`, `Can't tell`, `Looks harmless`), score, start-region cutoffs (`-7.75`, `-3.03` <!-- app/app_config.json cutoff_b cutoff_a -->), the CRISPR lab answer, and the verdict's held-out track record. For `C61G`, the demo returns `Looks harmful`, score `-10.50` <!-- app/app_data.csv C61G esm_score -->, `Loss of function (harmful). Lab function score -1.738.` <!-- app/app_data.csv C61G lab_score -->, and held-out track record `196 of 288` clear-class mutations (`68.1%`; `92` `FUNC` harmless, `29` `INT` in between <!-- band_counts.csv looks_harmful -->). BRCA1 is the exam, not the use case: Findlay et al. (2018) measured 96.5% <!-- BUILD_SPEC.md section 1 --> of possible single-letter DNA changes at these positions, so a counselor would use the lab result there while using this exam to judge sequence-only models on unmapped genes.

## What we added

Our own idea. Grading a language model on this lab data has been done before (ProteinGym; the Evo 2 BRCA1 notebook), and stronger tools such as AlphaMissense exist. We added four things: holding out a whole protein region with cutoffs frozen in advance; testing whether start-region labels help at the other end; counting misses and false alarms separately; and explaining a failure case.

## How it works

Findlay et al. (2018) via ProteinGym has 1,837 <!-- brca1_missense_scores_slim.csv rows --> single-letter changes in two regions (`1–101` start, `1631–1855` end) with study cutoffs (`< -1.328` `LOF`, `> -0.748` `FUNC` <!-- frozen.json labels -->) and 155 excluded `INT` rows (59 start, 96 end <!-- int_rows.csv n_int -->); we rebuild the 1,863-letter wild-type sequence to verify every row. Zero-shot `ESM-2` (`facebook/esm2_t33_650M_UR50D` <!-- frozen.json model_id -->)—which read the BRCA1 sequence in pretraining but never any lab label—hides each position and scores `log P(new) - log P(normal)` (`arm2`); in `outputs/heatmap.png`, at the 8 `zinc-holding letters (PDB 1JM7)` (`24, 27, 39, 41, 44, 47, 61, 64`), all 152/152 swaps are negative (`-15.61` to `-6.02`) with 47 black `LOF` and 2 grey `INT` dots. Against `arm0` (always harmless), `arm1` (`BLOSUM62`, flagging 337/516 = 65.3% start rows from integer ties at `0.0` <!-- frozen.json arm1 -->), and `arm3` (logistic regression on 8 features: ESM, BLOSUM62, entropy, Zamyatnin size, Kyte-Doolittle hydrophobicity, to-proline, from-cysteine, from-glycine), we freeze 90%-catch cutoffs on `1–101` (`outputs/frozen.json`) and test once on `1631–1855`, never splitting rows at random because same-position swaps behave alike.

## Results

Tested 1 time (`run=1` in `outputs/results/run_log.txt` <!-- run_log.txt run=1 -->) on 1,166 held-out end-region rows (267 `LOF`, 899 `FUNC` <!-- headline.csv n_harmful n_harmless -->):

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

Held-out AUROC was 0.5000 (`arm0`), 0.6681 (`arm1`), 0.8995 (`arm2`), and 0.8662 (`arm3`) <!-- auroc.csv auroc -->, with `esm_minus_blosum` `+0.231` (95% CI `+0.188` to `+0.272` <!-- bootstrap.json esm_minus_blosum -->; Spearman 0.506 overall, 0.430 start, 0.538 end <!-- sanity7.json spearman -->). Next to `S3`, the `arm3_minus_esm` 95% interval is `-0.054` to `-0.014` (point `-0.033` <!-- bootstrap.json arm3_minus_esm -->), entirely below zero; **what didn't work** was that start-region labels (`arm3`) lowered AUROC from 0.8995 to 0.8662 and missed 22 vs 15 harmful cases <!-- headline.csv harmful_missed -->, while `BLOSUM62` ties (`arm1`) flagged 571/899 (63.5%) harmless cases <!-- headline.csv arm1 false_alarms --> (and 72, 76, 67 of 96 end `INT` rows across `arm1`–`arm3` <!-- int_rows.csv end -->). Three of Arm 2's five worst misses sit at adjacent positions 1823–1824 (`A1823S`, `A1823T`, `I1824M` <!-- failures.csv rank 1..3 -->); all three arms missed our explained wrong case `A1823S` (`LOF` `-3.643`, `ESM-2` `+0.75`, `BLOSUM62` `+1` <!-- failures.csv A1823S -->). Guess: High entropy (2.572 <!-- failures.csv A1823S entropy -->) and a mild `A`-to-`S` swap suggest DNA-level splicing disruption invisible to a protein model.

**Groups check:** By region (`start` tuning vs held-out `end`; the secondary reverse split was cut before the test run) and swap kind (`by_swap.csv` <!-- by_swap.csv -->), `arm2` catches 93.5%–100.0% across swap kinds (`29/31` `to_proline`, `9/9` `loss_of_cysteine`, `91/96` `charge_change`, `123/131` `other`), while false alarms are 75.0% (`15/20`) on `loss_of_cysteine` and 57.4% (`178/310`) on `charge_change` vs 43.8% (`14/32`) on `to_proline` and 42.8% (`230/537`) on `other` <!-- by_swap.csv arm2 -->. Our data are mutations, not people, so we cannot check performance by age, sex or ancestry.

**Context, added after the test run (reference only, not our arms):** On ProteinGym's published leaderboard for this table, `ESM-2` 650M ranks 16th of 97 zero-shot models by Spearman (published 0.515, ours 0.506; best 0.560, median 0.456 <!-- reference_proteingym.json spearman, our_spearman_all_1837 -->), and the larger `ESM-2` models score lower there (0.497 for 3B, 0.394 for 15B <!-- outputs/reference/proteingym_brca1_leaderboard.csv spearman -->). AlphaMissense, which also uses protein structure and population data, was scored from its published predictions on the same 1,166 held-out rows with the same start-region cutoff rule: AUROC 0.9155 vs 0.8995 for `arm2` (difference `+0.016`, 95% interval `-0.003` to `+0.038`), catching 224 of 267 harmful mutations (43 missed) with 124 of 899 false alarms <!-- reference_alphamissense.json end, alphamissense_minus_esm -->.

**Extra arm, added after the test run (exploratory):** We then added `arm4`, ESM-C 600M zero-shot, scored with the same hide-a-letter method. It was chosen over the mean of the ESM-2 650M and ESM-C 600M scores by start-region AUROC alone (0.8655 vs 0.8572), with its cutoff set by the same 90% rule on `1–101` and frozen in `outputs/frozen_arm4.json` before one run on the end region (`outputs/results/arm4_run_log.txt`). It caught 250 of 267 harmful mutations (17 missed) with 423 of 899 false alarms; AUROC 0.9024 vs 0.8995 for `arm2` (difference `+0.003`, 95% interval `-0.015` to `+0.023`). No detectable difference from `arm2` at this sample size. We picked this arm after seeing `arm2`'s end-region results and ProteinGym's leaderboard, so it is exploratory and the headline table is unchanged. <!-- arm4.json -->

**Second extra arm, added after the test run (exploratory):** `arm5` keeps `arm2`'s ESM-2 scores and changes only the cutoff: instead of catching at least 90% of harmful mutations, it takes the cutoff with the best catch rate minus false-alarm rate on `1–101`, frozen in `outputs/frozen_arm5.json` before one run on the end region (`outputs/results/arm5_run_log.txt`). It caught 220 of 267 harmful mutations (47 missed, vs 15 for `arm2`) with 136 of 899 false alarms (vs 437); 84.3% of calls were right vs 61.2% (difference `+23.1` points, 95% interval `+19.7` to `+27.1`). On the held-out region the balanced cutoff gave more correct calls than `arm2`. Most of that gain comes from the end region having 899 harmless vs 267 harmful rows: always calling 'harmless' would already be right 77.1% of the time, and the catch rate fell from 94.4% to 82.4%. The ranking, and so the AUROC, is the same as `arm2`'s: this is a different trade between missed harmful mutations and false alarms, not a better model. We picked it after seeing `arm2`'s end-region results, so it is exploratory and the headline table is unchanged. <!-- arm5.json -->

**Limits:** (1) Blind to splicing—we could not test this. (2) Cells in a dish, not people (a 2025 preprint notes cell-type effects). (3) Stronger tools exist (AlphaMissense, free to download), and zero-shot scoring here is not a new method.

## What we'd do next

Add a DNA splicing predictor for the splicing blind spot, and repeat the frozen-cutoff held-out-region exam on a second gene with a lab map.

## Sources

- **Data:** Findlay GM et al. (2018), *Nature* 562:217–222 (https://doi.org/10.1038/s41586-018-0461-z); UW BRCA1 SGE (https://sge.gs.washington.edu/BRCA1/); ProteinGym (Notin P et al., 2023, https://github.com/OATML-Markslab/ProteinGym, also prior work).
- **Models & methods:** ESM-2 `facebook/esm2_t33_650M_UR50D` (Lin Z et al., 2023, *Science* 379:1123–1130, https://huggingface.co/facebook/esm2_t33_650M_UR50D, https://github.com/facebookresearch/esm); ESM-C 600M `biohub/ESMC-600M` (ESM Cambrian, EvolutionaryScale/Biohub; https://huggingface.co/biohub/ESMC-600M, https://github.com/evolutionaryscale/esm), used only for the exploratory `arm4`; Meier J et al. (2021, *NeurIPS*, https://doi.org/10.1101/2021.07.09.450648); BLOSUM62 (Henikoff & Henikoff, 1992) via Biopython (https://biopython.org); Kyte & Doolittle (1982); Zamyatnin (1972); PDB `1JM7` (https://www.rcsb.org/structure/1JM7).
- **Prior work & limits:** Evo 2 BRCA1 notebook (https://github.com/ArcInstitute/evo2, https://docs.nvidia.com/bionemo-recipes/2.7.1/main/examples/bionemo-evo2/zeroshot_brca1/); AlphaMissense (Cheng J et al., 2023, *Science* 381:eadg7492; published BRCA1 predictions from the AlphaFold Database, used only as a reference row, https://alphafold.ebi.ac.uk/files/AF-P38398-F1-aa-substitutions.csv); 2025 cell-type preprint (https://www.medrxiv.org/content/10.1101/2025.08.11.25333423.full.pdf).
- **Libraries & AI tools:** PyTorch, Transformers, NumPy, pandas, SciPy, scikit-learn, Biopython, Plotly, Matplotlib, Gradio, pytest, Jupyter (`nbformat`, `nbconvert`, `ipykernel`); **Claude Code (Anthropic)** (planning, build spec, review/audit/mock judge, and the post-test reference comparison in `src/reference.py`) and **Gemini in Google Antigravity** (writing and running project code).
