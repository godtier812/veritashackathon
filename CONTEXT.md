# BRCA Spellcheck: full project context

Written Sunday Nov 15 2026, about 5:15 PM US Central. Submission deadline: tonight, 10:59 PM Central (11:59 PM Eastern).

This file is a handover. It should let a person or an AI assistant with no history understand what the project is, what has been built, what the results are, and what is still open. The judging rules are in `RUBRIC_AND_RULES.md`. The organisers' original documents are in `docs/`.

**How current this is:** everything up to the stage 10 report (about 1 PM) comes from the build agent's own reports. Anything after that, including whether the repository was pushed and the form submitted, is not known here. Check before relying on it.

---

## 1. The project in one paragraph

BRCA Spellcheck is an entry in the Veritas AI International High School Hackathon 2026, AI + Health track, own idea. A user types a single-letter protein mutation in the BRCA1 gene, such as `C61G`, and the tool says whether it looks harmful. The verdict comes from ESM-2, a protein language model that learned only from natural protein sequences and was never shown patients or labels. The project grades that model against lab truth (a 2018 study that used CRISPR to test BRCA1 mutations in living cells), on a region of the protein that was kept out of all tuning.

**One-sentence pitch:** BRCA Spellcheck helps a genetic counselor whose patient got an "uncertain" BRCA1 test result by taking the mutation and saying whether it looks harmful, using an AI that learned the language of proteins and was graded in the open against CRISPR lab results.

**Honest framing used throughout:** BRCA1 is the exam, not the use case. At the positions tested here the lab answer already exists, so a counselor would use that. The point is to measure how far a sequence-only model can be trusted before someone leans on it where no lab map exists.

---

## 2. The team and how the work was done

- Two high school students. One led the biology and machine-learning thinking; one ran the coding agent and the front end.
- The idea was chosen after reviewing and discarding an earlier plan (a diabetes screener on NHANES data, judged not distinctive enough).
- Planning, the build specification and audits were done with Claude (Anthropic).
- The code was written and run by Gemini 4 Argon (High) as a multi-agent "swarm" in Google Antigravity, working from `BUILD_SPEC.md`, on a Mac, with CPU only.
- All code was written during the hackathon weekend, as the rules require.
- Every AI tool must be credited with its real role in the write-up and on the form.

---

## 3. Data

| File | What it is |
| --- | --- |
| `data/BRCA1_HUMAN_Findlay_2018.csv` | The ProteinGym copy of the lab results, unedited. 1,837 rows. |
| `data/brca1_missense_scores_slim.csv` | The same rows in a compact form, with region and class columns added. |
| `data/BRCA1_HUMAN_Findlay_2018_wildtype.fasta` | The normal BRCA1 protein, 1,863 letters. |

Verified facts:

- 1,837 single-letter protein changes at 326 positions, in two regions.
- **Start region (RING), positions 1–101:** 575 mutations (131 harmful, 59 in-between, 385 fine).
- **End region (BRCT), positions 1631–1855:** 1,262 mutations (267 harmful, 96 in-between, 899 fine).
- Classes use the study's own cutoffs: functional above −0.748, loss of function below −1.328, intermediate between.
- The raw file is byte-identical to the official ProteinGym download.
- Source study: Findlay et al. 2018, Nature, "Accurate classification of BRCA1 variants with saturation genome editing". The lab tested 3,893 DNA changes; this table holds the ones that change exactly one protein letter.
- The data are lab measurements on cells. There are no patient records.

---

## 4. Method

### The four arms

| Arm | What it does |
| --- | --- |
| 0. Always "harmless" | The simplest-guess baseline |
| 1. BLOSUM62 rule | A 1992 lookup table of how unusual each letter swap is |
| 2. ESM-2 zero-shot | Hide the letter, ask the model what belongs there, compare the mutant letter with the normal one |
| 3. ESM-2 plus labels | Logistic regression on the ESM score, the BLOSUM score and six simple chemistry features |

### Locked design

- **Tune on the start region. Test on the end region, once.** Rows are never split at random, because mutations at one position behave alike and a random split would leak.
- **Label:** harmful = loss of function; harmless = functional. The in-between rows are left out of the headline and reported separately.
- **Cutoff rule:** for each arm, on the start region, the loosest cutoff that still catches at least 90% of harmful mutations. Frozen, then applied unchanged to the end region.
- **Freeze before testing:** every cutoff and fitted number was written to `outputs/frozen.json` and hash-locked before the test run.
- **One test run,** recorded in a run log.
- **Model:** `facebook/esm2_t33_650M_UR50D`, two 500-letter windows (residues 1–500 and 1364–1863), masked scoring.
- **Bootstrap:** 1,000 resamples of end-region positions (not rows).

### Decisions made during the build

- The full scoring run was done on CPU (419 seconds), approved by a human.
- **Cut:** the reverse split (tune on end, test on start) and the interactive heat map. A static PNG heat map was kept.
- **Sentence S3** was reworded to: "We could not detect a benefit from adding labels from the other end of the protein."
- The eight zinc-holding letters of the RING region (C24, C27, C39, H41, C44, C47, C61, C64) were confirmed against the solved structure, Protein Data Bank entry 1JM7, and are marked on the heat map.
- A version-zero README was written right after the test run so the team could submit early.
- A track-record line was added to the demo (see section 6).

---

## 5. Results

All numbers are from the one test run on the held-out end region (267 harmful, 899 harmless), as reported by the build agent.

### Headline table

| Method | Harmful caught | Harmful missed | False alarms | Catch rate | AUROC |
| --- | --- | --- | --- | --- | --- |
| Always harmless | 0 | 267 | 0 | 0.0% | 0.500 |
| BLOSUM62 rule | 232 | 35 | 571 | 86.9% | 0.668 |
| ESM-2 zero-shot | 252 | 15 | 437 | 94.4% | 0.900 |
| ESM-2 plus labels | 245 | 22 | 405 | 91.8% | 0.866 |

The "always harmless" baseline is 77.1% accurate and catches nothing.

### Fixed result sentences (use word for word)

- **S1:** "A model that only read protein sequences ranked harmful BRCA1 mutations better than a classic substitution table."
- **S3:** "We could not detect a benefit from adding labels from the other end of the protein."
- **S6:** "To catch 94.4% of harmful mutations, the model also flagged 437 of 899 harmless ones."

### Intervals and correlations

- **ESM-2 minus BLOSUM62 (AUROC):** +0.231, 95% interval 0.188 to 0.272.
- **ESM-2 plus labels minus ESM-2 (AUROC):** −0.033, 95% interval −0.054 to −0.014. This lies entirely below zero: adding labels from the other end made ranking measurably worse.
- **Correlation with the lab score (Spearman):** 0.430 on the start region, 0.538 on the end region, 0.506 overall.

### What the results mean

- The sequence-only model clearly beats the 1992 table at ranking harmful mutations.
- The cost is large and should be said plainly: catching 94% of harmful mutations flags about half the harmless ones.
- It held up on a region it never saw.
- Labels from one end of the protein did not help at the other end, and slightly hurt.

### Tuning-region numbers (for reference, not results)

- BLOSUM62: cutoff 0.0, 120 of 131 caught, 217 of 385 false alarms. It flags 65.3% of tuning rows.
- ESM-2: cutoff 3.034498, 118 of 131 caught, 196 of 385 false alarms.
- "Can't tell" band for the demo: between 3.034498 and 7.751908.
- Checks passed before the test: the label-free letter test (326 of 326 positions), and `C61G` scoring −10.50, in the lowest quarter.

---

## 6. The demo

A Gradio app that needs no GPU. It reads pre-computed scores.

- **Input:** a mutation such as `C61G` (also accepts `c61g` and `p.Cys64Gly`).
- **Output:** "Looks harmful", "Can't tell" or "Looks harmless"; the model score; the lab result if that exact swap was tested; and a track-record line.
- **Refuses** unreadable input, positions outside the 326 tested ones, and a wrong first letter, each with a clear message.
- **Always shows** the disclaimer and the note that a counselor would use the lab result where one exists.

### Track record on the held-out region

| Verdict | Lab said harmful | Lab said harmless | Lab said in-between |
| --- | --- | --- | --- |
| Looks harmful | 196 | 92 | 29 |
| Can't tell | 56 | 345 | 47 |
| Looks harmless | 15 | 462 | 20 |

"Looks harmless" matched the lab on 462 of 477 clear cases. "Looks harmful" matched on 196 of 288. The tool is much more trustworthy when it clears a mutation than when it flags one.

---

## 7. Where it fails

- **Five worst misses** (harmful, called harmless): A1823S, A1823T, I1824M, R1670K, A1693D.
- **Five worst false alarms** (harmless, flagged): L1844P, E1661V, C1787F, Q1779P, V1653G.
- **A fact from the data:** three of the five worst misses sit at adjacent positions 1823 and 1824.
- **A guess, not a finding:** a cluster like that fits a splicing effect, which a protein model cannot see. This was not verified. The study's full table (with silent and splice-region changes) was never obtained, so splicing is stated as a limit that could not be tested.
- **Explained wrong case in the write-up:** A1823S. The model found the swap unsurprising, and the lab found loss of function.

---

## 8. Limits

- Lab truth comes from cells in a dish, which is strong but indirect evidence about cancer risk.
- A protein model is blind to changes that break how the gene is read (splicing).
- Coverage is small: 326 of 1,863 positions, one gene, single-letter protein changes only.
- ESM-2 has read the BRCA1 sequence during its training, though never the lab labels.
- The 90% catch rule is a choice; a different rule changes the trade-off.
- Stronger tools exist, such as AlphaMissense. This project does not claim to match them.

## 9. Never claim

- That the method is new. ProteinGym already grades protein language models on this table, and an NVIDIA Evo 2 notebook does similar BRCA1 scoring.
- That it beats or is cheaper than another tool.
- That it diagnoses, or classifies variants for clinical use.
- That scores at particular positions are explained by biology, beyond what the heat map visibly shows.

---

## 10. Build status at the last report (about 1 PM)

| Stage | Status |
| --- | --- |
| 1–3: setup, baselines, ESM-2 scores | Done |
| 4: features, freeze | Done, audit reported clean |
| 5: the one test run | Done, run 1, success |
| 5b: version-zero README | Done |
| 6: failure analysis | Done |
| 7: heat map (static PNG) | Done |
| 8 and 8b: demo, track-record line | Done |
| 9: notebook, executed with outputs | Done |
| 10: write-up, README, form answers | Done |
| 11: video script | Not confirmed |
| 12: judge Q&A sheet | Not confirmed |

- 49 tests passed at the last report.
- Repository address used in the notebook: `https://github.com/godtier812/veritashackathon`.
- The reviews were run by the same agent swarm that wrote the code, so they are not an independent check.

## 11. Open items

- [ ] Confirm the repository is pushed and opens in a private window without signing in.
- [ ] Confirm the notebook runs in Colab with Run all. The clone step could not be tested before the push.
- [ ] Upload the five `app/` files to a Hugging Face Space and test `C61G`, `D2E` and `hello`.
- [ ] Fill in the `Team:` line and the form's team details.
- [ ] Submit the form, then resubmit after any improvement. The organisers keep the latest version.
- [ ] Record the two-minute video (optional but recommended).
- [ ] Both teammates read the write-up and can explain every part.
- [ ] Look up A1823T on the University of Washington BRCA1 site or ClinVar. Only state a cause if a source confirms it.
- [x] Look up ProteinGym's published ESM-2 650M result for this table and compare it with 0.506. Done about 7 PM: published 0.515, 16th of 97 models, best 0.560 (`outputs/results/reference_proteingym.json`). An AlphaMissense reference row was added at the same time (`outputs/results/reference_alphamissense.json`): AUROC 0.9155 against 0.8995 for ESM-2 on the same held-out rows. No model, cutoff or stage 5 result changed.

---

## 12. What is in this pack, and what is not

**Included**

- `CONTEXT.md` (this file) and `RUBRIC_AND_RULES.md`.
- `docs/`: the organisers' three documents (project briefs, Hacker Guide, sample submission).
- `data/`: the three data files.
- `spec/BUILD_SPEC.md`, `spec/AGENTS.md`, `spec/START_PROMPT.md`: the audited version of the build specification.
- `spec/review_agents/`: the six review-agent definitions.
- `history/`: two earlier planning files, superseded.

**Not included, because they exist only on the build machine**

- The code (`src/`, `tests/`, `app/`), `notebook.ipynb`, `outputs/` and `writeup/`. Get them from the GitHub repository.

**One mismatch to know about**

The build ran on an earlier copy of `BUILD_SPEC.md`, with the important fixes passed in as human decisions. The copy in `spec/` is the audited version. Where they differ, the decisions in section 4 above and the repository's `writeup/decision_log.md` describe what was built: for example the heat map uses a −12 to 12 red-blue scale, and the reverse split was cut.

---

## 13. Sources

- Findlay et al. 2018, Nature: https://doi.org/10.1038/s41586-018-0461-z
- University of Washington BRCA1 database: https://sge.gs.washington.edu/BRCA1/
- ProteinGym: https://github.com/OATML-Markslab/ProteinGym
- ESM-2 model: https://huggingface.co/facebook/esm2_t33_650M_UR50D
- Protein Data Bank entry 1JM7 (zinc-holding letters): https://www.rcsb.org/structure/1JM7
