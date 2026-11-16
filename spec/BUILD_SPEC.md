# BRCA Spellcheck: build spec

This file is the complete specification for building the project. It stands alone: everything you need is written here. If another file in the repository disagrees with this one, this file wins.

Written Sunday Nov 15 2026, about 1:00 AM US Central. Revised about 1:45 AM after an outside audit. Hard deadline: Sunday Nov 15 2026, 10:59 PM US Central (11:59 PM US Eastern).

## Table of contents

- 0. How to work (instructions to the coding agent)
- 1. The project in plain words
- 2. Glossary
- 3. Locked rules
- 4. Repository layout
- 5. Environment
- 6. Data contracts
- 7. Constants
- 8. Module specs
- 9. Tests
- 10. Stages
- 11. Reviews
- 12. Write-up, form, video and Q&A
- 13. Known traps
- 14. Cut list and never-cut list
- 15. Troubleshooting
- 16. Definition of done

Words used with a fixed meaning in this file:

- **"the humans"** = the two team members. **"the operator"** is the one running you. **"the reviewer"** is the teammate who reviews at checkpoints. Never write either person's name in any file or commit.
- **"the test run"** = the single execution of `src/evaluate.py` on the end region (section 3, rule 8).
- **"estimate"** next to a number = my guess, not a fact. Do not assert on it in a test.
- **"verified"** next to a number = it was checked against the real files on Nov 15 before this spec was written.

---

## 0. How to work (instructions to the coding agent)

You are the coding agent for a two-person high school hackathon team. You write and run all the code. The humans review.

1. **Build stage by stage** (section 10). Do not start a stage until the previous one passes every acceptance criterion.
2. **At the start of each stage**, print the stage goal and the list of files you will create or change. Check the time (`Get-Date` in PowerShell) and print it.
3. **At the end of each stage**, run every command in the stage's "Verify" list yourself, read the output, and report:
   - files changed
   - commands run
   - the actual outputs (paste them, do not paraphrase numbers)
   - pass or fail for each acceptance criterion
   - three short bullet points a human could say to a judge who asks "why did you do it this way?"
4. **Never move past a failing test.** Never weaken, skip or delete a test to make it pass. If you believe a test is wrong, stop and explain why.
5. **Never change a locked rule** (section 3). If you think one is wrong, stop and say why.
6. **Never invent a number.** Results come only from files in `outputs/results/`. Until a result exists, write a placeholder such as `[X of 267]`.
7. **At every 🛑 checkpoint, stop completely.** Tell the operator exactly what to show the reviewer. Wait for the exact words given in the checkpoint before continuing. Checkpoint A is the one exception: it is automatic when its tests pass (stage 3).
8. **Never run `src/evaluate.py` without the exact go-ahead in checkpoint B** ("run the test"). The one exception is `python -m src.evaluate --selfcheck-start`, which never touches the end region (section 8.9).
9. **Never push to a git remote, log in to any service, or upload anything** without asking first.
10. **If you are stuck on one problem for 20 minutes**, stop and propose a simpler route, with what it costs.
11. **Do not add features, files or dependencies** that are not in this spec without asking.
12. **Code style:** simple and readable. Lowercase `snake_case` names. Short comments only where a line is not obvious. No clever one-liners. Type hints on every function. Every module has a one-paragraph docstring at the top saying what it does.
13. **Commit after each stage** with the commit message given in that stage. Use `git add -A` then `git commit -m "<message>"`. Never commit anything that looks like a token or password. Before your first commit, run `git config user.name` and `git config user.email`. They must print `BRCA Spellcheck Team` and `team@brca-spellcheck.invalid` (round 1 judging is blind, so commits must not carry a personal name or address). If they print anything else, stop and tell the operator.
14. **Windows and PowerShell are the main environment.** Run every command from the repository root. Colab notes appear only where Colab differs.
15. **If you are unsure about a library detail** (a function name, an argument), check the installed package (`python -c "help(...)"`) or its documentation. Do not guess. Several pinned libraries are recent major versions (transformers 5, Gradio 6, Plotly 7).
16. **If a download fails** (packages, model weights), report the exact error. Do not switch to a different model or source without asking.
17. **Do not read or print end-region lab labels before the test run.** Section 3, rule 9 says exactly what this means. This is the rule a coding agent is most likely to break by accident, for example with a quick `df.head()` or a debugging plot.

---

## 1. The project in plain words

**BRCA Spellcheck** is an entry in the Veritas AI International High School Hackathon 2026, AI + Health track, under "your own idea".

**The question.** When a person has a genetic test, the lab often reports a "variant of uncertain significance": a mutation nobody can yet call harmful or harmless. Most genes have no lab map of which mutations are harmful, so a genetic counselor often has only computer predictions to go on. Can a small, free model that only read protein sequences be trusted for that, and where does it go blind?

**The exam.** We test the model on BRCA1, a gene that helps repair DNA. Mutations that break it raise the risk of breast and ovarian cancer. BRCA1 is one of the few genes with a lab answer key: a 2018 study (Findlay et al., Nature) used CRISPR to make BRCA1 mutations in living human cells and measured whether each one still worked.

**BRCA1 is the exam, not the use case.** At the positions we test, the lab already measured 96.5% of all possible single-letter DNA changes. A counselor would use the lab result there, not our tool. The value of the project is what the exam tells us about genes that have no lab map. This framing must appear in the write-up, the notebook, the demo and the video.

**The AI.** ESM-2 is a protein language model. It was trained only on natural protein sequences, by hiding letters and learning to predict them. It never saw a patient, a diagnosis or a lab label. Using it this way, with no training on our labels, is called zero-shot. To score a mutation we hide the letter at that position, ask the model how likely each of the 20 letters is there, and compare the new letter with the normal one. A mutation the model finds very unlikely is a probable "typo".

**The four arms** (methods compared). Every arm gives each mutation a harm score, where higher means more likely harmful:

- Arm 0: always say "harmless". The simplest possible baseline.
- Arm 1: a 1992 lookup table of how unusual each letter swap is (BLOSUM62).
- Arm 2: ESM-2 zero-shot.
- Arm 3: a small logistic regression on top of the ESM score, the BLOSUM score and a few chemistry features.

**The split.** The lab data cover two separate regions of the protein. We tune every cutoff on the start region (positions 1–101) and test once on the end region (positions 1631–1855). Rows are never split at random, because mutations at the same position behave alike and a random split would leak the answer.

**What we add.** (1) A whole region of the protein held out, with the cutoff frozen in advance. (2) A test of whether labels from one end of the protein help at the other end. (3) Missed cases and false alarms counted separately. (4) Where it fails and why.

**What is already known. We must not overclaim.**

- ProteinGym, a public benchmark, already scores ESM-2 and many other models zero-shot on this exact table.
- A public Evo 2 notebook (Arc Institute, also published by NVIDIA) runs a zero-shot DNA model on the same study's 3,893 DNA changes.
- Stronger tools exist, such as AlphaMissense. As far as we know its predictions are free to download. We do not claim to be better or cheaper.

**The demo.** The user types a mutation such as `C61G`. The tool shows the model's verdict next to the lab's answer.

**The disclaimer**, word for word, in the write-up, the notebook and the demo:

> This project is a prototype and is not medical advice.

---

## 2. Glossary

| Term | Meaning in this project |
| --- | --- |
| Protein | A chain of building blocks, written as a string of letters. BRCA1 is 1,863 letters long. |
| Amino acid letter | One of 20 capital letters, each standing for one building block: `ACDEFGHIKLMNPQRSTVWY`. |
| Position | A 1-based index into the protein string. Position 1 is the first letter. |
| Wild type | The normal, unmutated protein, or the normal letter at a position. |
| Mutation notation (`C61G`) | Normal letter, position, new letter. `C61G` means position 61 normally holds `C` and the mutation swaps it for `G`. |
| Swap | One single-letter change. Each position has 19 possible swaps. |
| Region | One of two blocks of positions with lab data: `start` (1–101) or `end` (1631–1855). |
| RING | The biological name of the start region. |
| BRCT | The biological name of the end region. |
| Function score | The lab's measurement for a mutation. Higher means the protein still works. Range in our table: −3.72 to 1.31. |
| LOF | Loss of function: function score below −1.328. We call this **harmful**. |
| FUNC | Functional: function score above −0.748. We call this **harmless**. |
| INT | Intermediate: function score from −1.328 to −0.748. Left out of tuning and of the headline. |
| Harm score | A number an arm gives a mutation. Higher means more likely harmful. |
| Zero-shot | Using a model on a task without training it on any labels for that task. |
| Masked language model | A model trained to predict hidden pieces of a sequence. |
| Token | One unit of model input. For ESM-2, each protein letter is one token, plus one start token and one end token. |
| Logits | The raw numbers a model outputs for each possible token at a position. |
| Log-probability | The natural logarithm of a probability. Always zero or negative. |
| ESM score | log P(new letter) − log P(normal letter) at the hidden position. More negative means the model finds the swap more surprising. |
| Entropy | A measure of how unsure the model is at a position: −Σ p·ln(p) over the 20 letters. 0 means certain; ln(20) ≈ 2.996 means no idea. |
| BLOSUM62 | A 20 × 20 table of whole numbers from 1992. A high number means the two letters are often swapped in nature. |
| Cutoff | A harm-score threshold. A mutation is **flagged** when its harm score is at or above the cutoff. |
| Catch rate | Share of harmful mutations that are flagged. Also called recall or sensitivity. |
| Missed case | A harmful mutation that is not flagged. A false negative. |
| False alarm | A harmless mutation that is flagged. A false positive. |
| AUROC | A cutoff-free score from 0.5 (coin flip) to 1.0 (perfect): the chance that a random harmful mutation has a higher harm score than a random harmless one. |
| Bootstrap | Resample the test data with replacement many times and see how much a result moves. |
| Out-of-fold | A prediction for a row made by a model that was not trained on that row. |
| Grouped folds | Cross-validation where all rows with the same position stay in the same fold. |
| Leakage | Any way information from the test data reaches tuning. It makes results look better than they are. |
| Frozen file | `outputs/frozen.json`. Holds every tuned number, written before the test run and never changed after. |
| Band | The demo's three verdicts: "Looks harmful", "Can't tell", "Looks harmless". |
| Splicing | How a cell cuts and joins a gene's message before making protein. Some DNA changes are harmful because they break splicing, which a protein model cannot see. |
| VUS | Variant of uncertain significance: a mutation nobody can yet call harmful or harmless. |
| Window | A 500-letter slice of the protein given to the model. |
| Arm | One of the four methods being compared. |
| CRISPR | A lab technique for editing DNA. The 2018 study used it to make each mutation in living cells. |
| Hash (SHA-256) | A 64-character fingerprint of a file's bytes. If one byte changes, the hash changes. |
| Spearman correlation | A number from −1 to 1 saying how well two lists agree in rank order. Positive means they rise together. |
| In-sample | Measured on the same rows that were used for tuning. Such numbers flatter the method and are never reported as results. |
| Tuning | Any choice made by looking at labels: a cutoff, a model fit, a scaler, a feature choice. |
| Pipeline | A scikit-learn object that chains a scaler and a model so both are fit on the same rows. |
| Nats | The unit of entropy when the natural logarithm is used. |
| Gradio | A Python library that turns a function into a small web page with a text box and a button. |
| Hugging Face | A website that hosts models (we download ESM-2 from it) and small free web apps, which it calls Spaces. |
| Colab | Google Colab: free notebooks in the browser, with an optional GPU. Judges may open our notebook there. |
| Notebook | `notebook.ipynb`: a file of text cells and code cells with their outputs saved. |

---

## 3. Locked rules

Do not change any of these. If you think one is wrong, stop and ask.

### 3.1 The ten locked decisions

**Rule 1. Label.** Harmful = `LOF`. Harmless = `FUNC`. The 155 `INT` rows (59 in the start region, 96 in the end region) are left out of the headline and of **all** tuning: every cutoff, the arm 3 model fit, the scaler and the cross-validation. They are reported separately.

**Rule 2. Split.** Tune on the start region (positions 1–101). Test on the end region (positions 1631–1855), once. Never split rows at random.

**Rule 3. Four arms.** Every arm produces a harm score where higher means more likely harmful.

- Arm 0: always "harmless". It flags nothing. It has no cutoff.
- Arm 1: BLOSUM62 from Biopython. Harm score = minus the BLOSUM62 value of the swap.
- Arm 2: ESM-2 zero-shot. Harm score = minus the ESM score.
- Arm 3: logistic regression (scikit-learn defaults, features standardised) on eight features: the ESM score, the BLOSUM score, the model's entropy at that position, size change, hydrophobicity change, and flags for "to proline", "from cysteine", "from glycine". Harm score = predicted probability of harmful.

**Rule 4. Cutoff rule.** A mutation is flagged when its harm score is **at or above** the cutoff (`>=`).

- For each of arms 1, 2 and 3, on the start region, take the **highest** cutoff that still flags at least 90% of the 131 harmful mutations. That is the cutoff with the fewest flags. (90% of 131 is 117.9, so at least 118 must be flagged.)
- Mutations with equal harm scores always get the same call. So the cutoff is always one of the observed harm scores.
- For arm 3, choose the cutoff from out-of-fold predictions (5 folds, grouped by position), then refit the model on the whole start region. Otherwise the cutoff would be chosen on rows the model has already memorised.
- Freeze every cutoff and apply it unchanged to the end region.

**Rule 5. Headline.** One table: for each arm on the end region, harmful caught, harmful missed, false alarms, and the catch rate actually reached.

**Rule 6. Secondary results.**

- AUROC with harmful = 1 and the harm score as the prediction.
- A bootstrap interval for AUROC of ESM (arm 2) minus BLOSUM (arm 1): 1,000 resamples of end-region **positions** (not rows), 2.5th to 97.5th percentile.
- The reverse split (tune on the end region, test on the start region), clearly labelled as secondary.
- Results by kind of swap: to proline, loss of cysteine, charge change, other.
- What each arm says about the `INT` rows.
- Only if the study's full table is in `data/` and a human has confirmed its columns: the headline again without the rows flagged as splice-disrupting.

**Rule 7. No hyperparameter search.** No trying several settings and keeping the best. Every setting is fixed in this spec.

**Rule 8. One test run.** Before it, `outputs/frozen.json` is written with every cutoff, the feature list, the model settings, the fitted model numbers, the model name and a hash of `outputs/esm_scores.csv`. The test run is one script, `src/evaluate.py`, which reads `frozen.json` and refuses to run if the hash does not match. A rerun is allowed only to fix a code bug. Every run is logged in `outputs/results/run_log.txt` with the time and the reason.

**Rule 9. No peeking.** Before the test run, no code may compare scores with end-region lab labels. That means:

- no correlations, plots, tables, prints or `head()` of data frames that hold end-region `function_score`, `findlay_class` or `proteingym_bin` next to any harm score, ESM score or feature;
- no checks on M1775R or A1708E (two famous end-region mutations);
- no full-table correlation (it includes the end region).

Allowed before the test run: the end-region label **counts** already written in this spec (267 LOF, 96 INT, 899 FUNC), and the label-free letter test on end-region positions. Computing ESM scores for end-region positions is allowed, because it uses no labels.

**Rule 10. "Can't tell" band.** Demo only, arm 2 only, set on the start region.

- Cutoff A is arm 2's 90%-catch cutoff from rule 4.
- Cutoff B is the lowest cutoff that flags at most 10% of the start region's 385 harmless mutations (so at most 38).
- Harm score at or above B: "Looks harmful". Below A: "Looks harmless". In between: "Can't tell".
- If B is not above A there is no band. Say so, and use A alone: at or above A is "Looks harmful", below A is "Looks harmless".
- After the test run, report how many end-region mutations fall in each band.

### 3.2 Five further locked rules

**Rule 11. Start region only before the test run.** Before the test run, the heat map draws the start region only, and the demo is tried on start-region mutations only.

**Rule 12. The reverse split lives inside `src/evaluate.py`** and runs only after the headline table has been written to disk. It is the only code allowed to tune on end-region labels.

**Rule 13. `frozen.json` holds every fitted number:** arm 3's coefficients, intercept, scaler means and scaler scales, and cutoff B. `evaluate.py` never refits a model for the headline.

**Rule 14. BLOSUM62 ties are expected.** BLOSUM62 values are whole numbers, so many mutations share a harm score and arm 1's cutoff flags a large share. Report the share. Do not "fix" it by adding tie-breakers or noise.

**Rule 15. End-region lab labels are switched on by the run log.** The heat map's end-region lab markers and the demo's end-region lab results appear only after `outputs/results/run_log.txt` records a successful test run. The exact check, used everywhere, is the function `test_run_done()` in `src/load.py`:

1. If `outputs/results/run_log.txt` does not exist, return `False`. (The path is read from `config.RUN_LOG_PATH` at call time.)
2. Read it as UTF-8 text. For each non-empty line, split on the character `|` and strip spaces from each part.
3. If any line has a part exactly equal to `status=success`, return `True`.
4. Otherwise return `False`.

Only `src/evaluate.py` may write to `run_log.txt`. Never create or edit that file by hand.

A second function in `src/load.py`, `test_attempted()`, answers a different question: has the test been **started**, even if it did not finish? It returns `True` if any line of the run log has a part exactly equal to `status=started`, `status=success` or `status=failed`. Once it is `True`, end-region labels may already have been compared with scores, so nothing that was frozen may change: `src/freeze.py` refuses to run, `src/esm_score.py` refuses to recompute, and `src/evaluate.py` refuses to run again without a bug-fix reason. `test_run_done()` (success only) is what switches end-region labels **on**; `test_attempted()` is what locks tuning **off**.

### 3.3 Other locked settings

- **Seed 42 everywhere.** Every place a seed must be set:
  - Python: `random.seed(42)`
  - NumPy legacy: `numpy.random.seed(42)`
  - PyTorch: `torch.manual_seed(42)` and `torch.cuda.manual_seed_all(42)`
  - scikit-learn: `random_state=42` in `LogisticRegression` and in `StratifiedGroupKFold`
  - the bootstrap: `numpy.random.default_rng(42)`
  - `src/config.py` provides `set_seeds()` for the first three. Call it at the start of every command-line entry point.
- **INT rows are excluded from all tuning** (repeated here because it is easy to forget in a new function).
- **No model is ever fit, and no scaler is ever fit, on rows from more than one region.**
- **ESM scores are computed once and cached.** After `frozen.json` exists, `outputs/esm_scores.csv` must not change by a single byte.

### 3.4 Result sentences fixed in advance

Use the sentence that matches the result, word for word, and claim nothing stronger. The rules that choose the sentence are in section 8.9, step 10.

- **S1.** "A model that only read protein sequences ranked harmful BRCA1 mutations better than a classic substitution table."
- **S2.** "No detectable difference at this sample size."
- **S3.** "We could not detect a benefit from adding labels from the other end of the protein."
- **S4.** No fixed wording. Say plainly that ESM barely beat the rule; the finding is that a 1992 table is hard to beat.
- **S5.** "A cutoff tuned on one end of the protein did not carry over to the other: it caught X% instead of 90%."
- **S6.** "To catch X% of harmful mutations, the model also flagged Y of 899 harmless ones."

### 3.5 Never claim

Never write, in code comments, the notebook, the demo, the write-up or the video script, that:

- the method is new, or that grading a protein language model against CRISPR data is new;
- we are better than, or cheaper than, AlphaMissense or any other tool;
- the tool diagnoses anything;
- the tool classifies variants for clinical use.

The tool is a prototype that shows how far a sequence-only model can be trusted.

### 3.6 The organisers' rules

Breaking these costs the team its submission.

- All code is written this weekend. Do not copy a finished project from anywhere.
- Results must come from data the model did not see while being tuned.
- Report missed cases and false alarms separately, not only accuracy.
- Compare with a simple baseline and show both numbers.
- Show at least one case the model gets wrong and explain why.
- Public data only. No patient records.
- Honest results. Weak results are reported as they are.
- The write-up must contain the exact line: "This project is a prototype and is not medical advice."
- Credit every dataset, model, library and AI tool.
- Free tools only.
- Round 1 judging is blind: reviewers must not see names or schools. See section 12.2.

---

## 4. Repository layout

`C` = committed to git. `I` = ignored by git.

```
.
├── AGENTS.md                     C  standing rules for any coding agent. Exists already.
├── BUILD_SPEC.md                 C  this file. Exists already.
├── START_PROMPT.md               C  first message for the coding agent. Exists already.
├── .gitignore                    C  exists already. Contents below.
├── .gitattributes                C  exists already. Stops git changing line endings.
├── README.md                     C  stage 1 writes a stub; stage 5b turns it into version zero; stage 10 replaces it with the write-up.
├── requirements.txt              C  stage 1. Pinned versions (section 5.3).
├── pytest.ini                    C  stage 1. Makes `src` importable in tests (section 5.6).
├── notebook.ipynb                C  stage 9. The version judges open. Saved with outputs.
├── data/                            input files. Never modified by any code.
│   ├── BRCA1_HUMAN_Findlay_2018.csv              C  exists. The raw lab table from ProteinGym.
│   ├── brca1_missense_scores_slim.csv            C  exists. Compact table built from the raw one.
│   ├── BRCA1_HUMAN_Findlay_2018_wildtype.fasta   C  exists. The normal protein, 1,863 letters.
│   └── DMS_ProteinGym_substitutions/             I  the full ProteinGym download, about 1 GB. Not used.
├── src/
│   ├── __init__.py               C  stage 1. Empty file.
│   ├── config.py                 C  stage 1. Every constant and path.
│   ├── load.py                   C  stage 1. The only module that reads files in data/ and the score cache.
│   ├── baselines.py              C  stage 2. Arm 0 and arm 1.
│   ├── cutoffs.py                C  stage 2. The shared cutoff rule, cutoff B and the band.
│   ├── esm_score.py              C  stage 3. Runs ESM-2 and writes the score cache.
│   ├── features.py               C  stage 4. Arm 3 features, swap kinds, arm 3 probability.
│   ├── freeze.py                 C  stage 4. Tunes on the start region and writes frozen.json.
│   ├── evaluate.py               C  stage 4 (written), stage 5 (run once). The test run.
│   ├── failures.py               C  stage 6. Worst misses and false alarms.
│   ├── heatmap.py                C  stage 7. The heat map.
│   └── export_app.py             C  stage 8. Builds the two small files the demo reads.
├── tests/
│   ├── test_units.py             C  stage 1 onwards. Unit tests on made-up data.
│   └── test_sanity.py            C  stage 1 onwards. Sanity tests 1–7 on the real data.
├── tools/
│   └── build_notebook.py         C  stage 9. Generates notebook.ipynb.
├── outputs/
│   ├── esm_scores.csv            C  stage 3. 6,194 rows.
│   ├── frozen.json               C  stage 4.
│   ├── heatmap.html              C  stage 7.
│   ├── heatmap.png               C  stage 7.
│   ├── selfcheck_start/          I  stage 4. Scratch output of the start-region self-check.
│   └── results/                  C  stage 5 and 6. Every result file (section 6.6).
├── app/                             the Hugging Face Space. Self-contained: reads nothing outside app/.
│   ├── app.py                    C  stage 8.
│   ├── app_data.csv              C  stage 8. Built by src/export_app.py.
│   ├── app_config.json           C  stage 8. Built by src/export_app.py.
│   ├── requirements.txt          C  stage 8.
│   └── README.md                 C  stage 8. Has the Space header.
├── writeup/
│   ├── decision_log.md           C  stage 1 onwards. One line per choice and why.
│   ├── credits.md                C  stage 1 onwards. Everything we use, with links.
│   ├── writeup.md                C  stage 10.
│   ├── form_answers.md           C  stage 10.
│   ├── video_script.md           C  stage 11.
│   └── qa_sheet.md               C  stage 12.
├── docs/                         I  the organisers' three PDFs. Theirs, not ours to republish.
├── .venv/                        I  the virtual environment.
├── CLAUDE.md, .claude/           I  Claude Code configuration and the six review agents (section 11). Never edit or commit them.
├── MASTER_PROMPT.md, PROJECT_BRIEF.md   I  superseded planning files. Do not read, edit or commit them.
```

**`.gitignore`, in full (already in place; do not change it without asking):**

```
# environments and caches
.venv/
__pycache__/
*.pyc
.pytest_cache/
.ipynb_checkpoints/

# the full ProteinGym download (about 1 GB); only the BRCA1 file is used
data/DMS_ProteinGym_substitutions/

# organisers' PDFs: theirs, not ours to republish
docs/

# scratch output of the start-region self-check
outputs/selfcheck_start/

# planning files for a different tool; they also name the team (round 1 is blind)
CLAUDE.md
MASTER_PROMPT.md
PROJECT_BRIEF.md
.claude/
```

**`.gitattributes`, in full (already in place):**

```
*.csv   -text
*.json  -text
*.fasta -text
*.ipynb -text
```

Why: `-text` tells git never to convert line endings in these files. Without it, git on Windows can turn `LF` into `CRLF` on checkout, which changes the bytes of `outputs/esm_scores.csv` and makes its hash stop matching `frozen.json`.

**Why outputs are committed:** judges and the notebook load `outputs/esm_scores.csv`, `outputs/frozen.json` and `outputs/results/` from GitHub, so nobody needs a GPU to check the numbers.

---

## 5. Environment

### 5.1 Python and the virtual environment

Python **3.12** (3.12.10 is on the machine this spec was checked on). Python 3.11 is acceptable. Do not use 3.13 or later without asking.

PowerShell, from the repository root:

```powershell
python --version
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
```

Always call the environment's Python by its path, `.\.venv\Scripts\python.exe`, instead of activating the environment. Activation can fail under PowerShell's script policy, and an agent's shell does not keep state between commands. Everywhere this spec writes `python`, it means `.\.venv\Scripts\python.exe`.

### 5.2 PyTorch

Install PyTorch **first**, before `requirements.txt`, so the GPU build is used:

```powershell
.\.venv\Scripts\python.exe -m pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu128
```

The download is about 3 GB. This build (`2.11.0+cu128`) is verified to see the RTX 5060 on the GPU laptop.

Verify:

```powershell
.\.venv\Scripts\python.exe -c "import torch; print(torch.__version__); print(torch.cuda.is_available()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'no gpu'); x = torch.ones(1000, 1000, device='cuda' if torch.cuda.is_available() else 'cpu'); print((x @ x).sum().item())"
```

Expected on the GPU laptop, four lines:

```
2.11.0+cu128
True
NVIDIA GeForce RTX 5060 Laptop GPU
1000000000.0
```

On a machine with no NVIDIA GPU the second line is `False` and the third is `no gpu`. That is fine for every stage except stage 3 (section 5.5). If the machine has no NVIDIA GPU you may install the plain build instead: `pip install torch==2.11.0`.

### 5.3 `requirements.txt` (main project), in full

```
# Install torch first with the GPU index URL (see README), then: pip install -r requirements.txt
torch==2.11.0
transformers==5.23.0
numpy==2.4.4
pandas==2.3.3
scipy==1.17.1
scikit-learn==1.8.0
biopython==1.88
plotly==7.1.0
matplotlib==3.10.9
gradio==6.32.0
pytest==9.2.0
nbformat==5.11.1
nbconvert==7.17.2
ipykernel==7.4.0
```

Install:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

All fourteen pins exist on PyPI and resolve together without conflict (checked Nov 15 2026). Several are deliberately not the newest release: do not upgrade them. torch, transformers, numpy, pandas, scipy, scikit-learn, biopython and matplotlib were installed and imported together; the others were resolved by pip but not installed. If pip reports a conflict, stop and report the exact message. Do not unpin anything on your own.

Nothing else may be imported. In particular: no `kaleido`, no `seaborn`, no `statsmodels`, no `fair-esm`.

### 5.4 `app/requirements.txt` (the Hugging Face Space), in full

```
pandas==2.3.3
```

No torch and no transformers: the demo only reads a small CSV. Gradio is not listed because a Space installs the Gradio version named in `app/README.md` (`sdk_version: 6.32.0`).

### 5.5 The GPU note and the handoff

Stage 3 (running ESM-2) needs an NVIDIA GPU: the team's RTX 5060 laptop (8 GB) or a free Colab T4. Every other stage runs on any machine, because everything after stage 3 reads `outputs/esm_scores.csv`.

If the machine you are running on prints `False` for `torch.cuda.is_available()`:

1. Finish writing `src/esm_score.py` and its unit tests. Run `python -m src.esm_score --letter-test-only`, which needs no GPU and no model weights.
2. 🛑 Stop. Tell the operator: "Stage 3 needs a GPU. Choose route A or route B."
3. **Route A, the GPU laptop.** The humans get the current commit onto the GPU laptop (a git push and pull, which needs their approval, or copying the repository folder). There, from the repository root with the environment set up as in 5.1–5.3:
   `.\.venv\Scripts\python.exe -m src.esm_score --only C61G`, then `.\.venv\Scripts\python.exe -m src.esm_score`.
   Then commit `outputs/esm_scores.csv` and bring it back the same way.
4. **Route B, Colab.** A human opens a new Colab notebook, sets the runtime to T4 GPU, uploads a zip of the repository (or clones it once it is on GitHub), and runs:
   `!pip install -q biopython==1.88` then `!python -m src.esm_score --only C61G` then `!python -m src.esm_score`.
   Colab already has torch, transformers, pandas and numpy; do not install `requirements.txt` there. The human downloads `outputs/esm_scores.csv` and puts it in `outputs/` on the build machine.
5. However the file travels, it must arrive **byte for byte**. Do not open and re-save it in Excel or an editor.
6. On the build machine, run `python -m src.esm_score --validate-cache`, then continue stage 3 from its "Verify" list.

### 5.6 Making `src` importable in tests

Use `pytest.ini` in the repository root, with exactly this content:

```
[pytest]
pythonpath = .
testpaths = tests
addopts = -q
```

`src/__init__.py` is an empty file. Run modules as `python -m src.<name>` from the repository root, never as `python src/<name>.py`.

Run the tests with:

```powershell
.\.venv\Scripts\python.exe -m pytest tests/ -q
```

### 5.7 Colab differences (for the notebook, stage 9)

- Colab has its own Python and its own torch, transformers, pandas, numpy, scikit-learn, scipy, plotly and matplotlib. Do not install `requirements.txt` in Colab. The notebook installs only `biopython==1.88`.
- Paths use `/`. All code must build paths with `pathlib.Path`, never by joining strings with `\`.
- The notebook's first code cell clones the GitHub repository if it is running in Colab (section 10, stage 9).

---

## 6. Data contracts

All CSV files: UTF-8, comma-separated, a header row, `LF` line endings, no index column. Write them with `df.to_csv(path, index=False, lineterminator="\n")` on a file opened with `open(path, "w", newline="", encoding="utf-8")`. Read them with explicit dtypes.

### 6.1 `data/BRCA1_HUMAN_Findlay_2018.csv` (input, read-only)

The raw lab table as published by ProteinGym. 1,837 data rows. Size 3,463,590 bytes. SHA-256 `abde5d47e84b1feeac6c32b60c1f3b5d2881accf3ceccfbadccceeea73e81c43` (verified).

| Column | Type | Meaning | Allowed values | Nulls |
| --- | --- | --- | --- | --- |
| `mutant` | str | The mutation, such as `C61G` | one letter, 1–4 digits, one letter; unique | none |
| `mutated_sequence` | str | The full protein with this one change applied | 1,863 letters | none |
| `DMS_score` | float | The lab's function score | −3.718221153 to 1.306892088 | none |
| `DMS_score_bin` | int | ProteinGym's own yes/no label | 0 (474 rows) or 1 (1,363 rows) | none |

Example row: `M1I,IDLSALRVEEVQNVINAMQK…(1,863 letters),-2.502128,0`.

Only `src/load.py` reads this file, and only for the rebuild check (section 8.2).

### 6.2 `data/brca1_missense_scores_slim.csv` (input, read-only, already built)

The table every stage uses. 1,837 data rows, in the same order as the raw file. Size 71,925 bytes. SHA-256 `4eb3d5597a4a2fcad88fa371b55cee9d8028c6df1ea18b95ad23802484e1a0ac` (verified).

| Column | Type | Meaning | Allowed values | Nulls |
| --- | --- | --- | --- | --- |
| `mutant` | str | The mutation | unique | none |
| `wild_type` | str | Normal letter: the first character of `mutant` | one of the 20 letters | none |
| `position` | int | The digits in the middle of `mutant` | 1–101 or 1631–1855 | none |
| `new_letter` | str | The last character of `mutant` | one of the 20 letters, never equal to `wild_type` | none |
| `region` | str | Derived from `position` only | `start` if 1 ≤ position ≤ 101; `end` if 1631 ≤ position ≤ 1855 | none |
| `function_score` | float | Copied unchanged from `DMS_score` | as above | none |
| `findlay_class` | str | Derived from `function_score` | `LOF` if score < −1.328; `FUNC` if score > −0.748; otherwise `INT` | none |
| `proteingym_bin` | int | Copied from `DMS_score_bin`. Not used by any arm. | 0 or 1 | none |

Example rows (verified):

```
mutant,wild_type,position,new_letter,region,function_score,findlay_class,proteingym_bin
M1I,M,1,I,start,-2.502128,LOF,0
C61G,C,61,G,start,-1.738348,LOF,0
D2E,D,2,E,start,0.088696,FUNC,1
```

(The file stores `function_score` with more digits than shown here; do not round it.)

Verified facts. Use them as tests:

| Fact | Value |
| --- | --- |
| Rows | 1,837, all `mutant` values unique |
| Distinct positions | 326 = 101 (every position 1–101) + 225 (every position 1631–1855) |
| Start region classes | 131 LOF, 59 INT, 385 FUNC, total 575 |
| End region classes | 267 LOF, 96 INT, 899 FUNC, total 1,262 |
| Both regions | 398 LOF, 155 INT, 1,284 FUNC |
| Lab-tested swaps per position | 1 to 7, average 5.63. Not every swap was tested. |
| Scores exactly equal to a class cutoff | none |
| C61G, C64G (start) | both `LOF`. (M1775R and A1708E in the end region are also `LOF`, but **never write a test or a check on those two before the test run**; rule 9.) |
| Lab-tested swaps at position 61 | `C61F`, `C61G`, `C61R`, `C61S`, `C61W`, `C61Y`. `C61A` is **not** lab-tested. |

**Reading rule:** always read with `dtype={"mutant": str, "wild_type": str, "position": "int64", "new_letter": str, "region": str, "function_score": "float64", "findlay_class": str, "proteingym_bin": "int64"}`.

### 6.3 `data/BRCA1_HUMAN_Findlay_2018_wildtype.fasta` (input, read-only, already built)

The normal BRCA1 protein. Size 1,981 bytes. SHA-256 `47529e88620d4d071f3ebec5212a399d74e102bc48f1ef5a1e1215186a196c36` (verified).

- Line 1 is a header starting with `>`: `>BRCA1_HUMAN wild type, 1863 aa, rebuilt from ProteinGym BRCA1_HUMAN_Findlay_2018.csv`
- The following lines hold the sequence, 60 letters per line (the last line is shorter).
- To read it: drop every line that starts with `>`, strip whitespace from the rest, join them. The result is one string of exactly 1,863 letters.

Verified facts about the sequence (call it `seq`; Python indexing is 0-based, positions are 1-based, so position `p` is `seq[p - 1]`):

| Fact | Value |
| --- | --- |
| Length | 1,863 |
| First 20 letters | `MDLSALRVEEVQNVINAMQK` |
| Letters outside the 20 standard ones | none |
| `seq[60]` (position 61) | `C` |
| `seq[63]` (position 64) | `C` |
| `seq[1707]` (position 1708) | `A` |
| `seq[1774]` (position 1775) | `M` |
| Positions 1, 101, 1631, 1855 | `M`, `Y`, `S`, `I` |
| `seq[0:500]` first 5 and last 5 letters | `MDLSA` … `RPLTN` |
| `seq[1363:1863]` first 5 and last 5 letters | `NLGEA` … `PHSHY` |

### 6.4 `outputs/esm_scores.csv` (written by stage 3, then read-only)

One row per position × new letter: 326 positions × 19 swaps = **6,194 rows**. Sorted by `position` ascending, then by `new_letter` in the fixed letter order `ACDEFGHIKLMNPQRSTVWY`. There is no row where `new_letter` equals `wild_type`.

| Column | Type | Meaning | Allowed values | Nulls |
| --- | --- | --- | --- | --- |
| `position` | int | Protein position, 1-based | 1–101, 1631–1855 | none |
| `region` | str | From position | `start`, `end` | none |
| `wild_type` | str | Normal letter at this position | one of 20 | none |
| `new_letter` | str | The swapped-in letter | one of 20, not `wild_type` | none |
| `mutant` | str | `wild_type` + `position` + `new_letter` | such as `C61G`; unique | none |
| `window` | str | Which window the model saw | `A` for start, `B` for end | none |
| `token_index` | int | Index of this position in the model input | 1–101 for A, 268–492 for B | none |
| `token_at_index` | str | The token the tokenizer put at `token_index` | must equal `wild_type` | none |
| `esm_score` | float | log P(new) − log P(normal), natural log, rounded to 6 decimals | any real; mostly negative | none |
| `entropy` | float | Entropy over the 20 letters at this position, in nats, rounded to 6 decimals. Same for all 19 rows of a position. | 0 to 2.995733 | none |
| `is_lab_tested` | int | 1 if `mutant` is in the slim table, else 0 | 0 or 1; exactly 1,837 ones | none |
| `model_id` | str | Hugging Face model ID used | such as `facebook/esm2_t33_650M_UR50D`; same in every row | none |
| `transformers_version` | str | `transformers.__version__` at scoring time | such as `5.23.0`; same in every row | none |
| `precision` | str | Floating-point precision the model ran in | `float32` (or `float16` only if a human approved it) | none |

Example row. Index columns are verified; the two score values are made up for illustration:

```
position,region,wild_type,new_letter,mutant,window,token_index,token_at_index,esm_score,entropy,is_lab_tested,model_id,transformers_version,precision
61,start,C,G,C61G,A,61,C,-9.123456,0.412345,1,facebook/esm2_t33_650M_UR50D,5.23.0,float32
```

Write floats with `float_format="%.6f"`.

### 6.5 `outputs/frozen.json` (written by stage 4, then read-only)

UTF-8 JSON, written with `json.dump(obj, f, indent=2, sort_keys=False)` and a final newline. No intermediate features file is used anywhere; features are cheap and are rebuilt in memory from `esm_scores.csv`.

Every key is required. Values shown are examples of the format, not predictions, except where marked verified.

```json
{
  "schema_version": 1,
  "created_at": "2026-11-15T03:10:42-06:00",
  "code_commit": "0123abc…(40 hex characters)",
  "code_dirty": false,
  "seed": 42,
  "model_id": "facebook/esm2_t33_650M_UR50D",
  "inputs": {
    "esm_scores_path": "outputs/esm_scores.csv",
    "esm_scores_sha256": "…(64 hex characters)",
    "esm_scores_rows": 6194,
    "slim_csv_sha256": "4eb3d5597a4a2fcad88fa371b55cee9d8028c6df1ea18b95ad23802484e1a0ac"
  },
  "tuning_region": "start",
  "test_region": "end",
  "labels": {
    "harmful_class": "LOF",
    "harmless_class": "FUNC",
    "excluded_class": "INT",
    "lof_below": -1.328,
    "func_above": -0.748
  },
  "tuning_counts": {"harmful": 131, "harmless": 385, "int_excluded": 59},
  "catch_target_percent": 90,
  "false_alarm_target_percent": 10,
  "arms": {
    "arm0": {
      "name": "Always harmless",
      "cutoff": null
    },
    "arm1": {
      "name": "BLOSUM62 rule",
      "harm_score": "0 - blosum62(wild_type, new_letter)",
      "cutoff": 0.0,
      "tuning_caught": 120,
      "tuning_false_alarms": 217
    },
    "arm2": {
      "name": "ESM-2 zero-shot",
      "harm_score": "0 - esm_score",
      "cutoff": 6.543210,
      "tuning_caught": 118,
      "tuning_false_alarms": 99
    },
    "arm3": {
      "name": "Logistic regression on ESM-2 + BLOSUM62 + chemistry",
      "harm_score": "predicted probability of harmful",
      "features": ["esm_score", "blosum_score", "entropy", "size_change",
                   "hydrophobicity_change", "to_proline", "from_cysteine", "from_glycine"],
      "scaler_mean": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      "scaler_scale": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
      "coefficients": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      "intercept": 0.0,
      "cutoff": 0.123456,
      "cutoff_source": "out-of-fold probabilities, StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42), groups = position",
      "model_settings": {"C": 1.0, "solver": "lbfgs", "max_iter": 100, "random_state": 42},
      "converged": true,
      "out_of_fold_caught": 118,
      "out_of_fold_false_alarms": 99,
      "tuning_caught": 120,
      "tuning_false_alarms": 90
    }
  },
  "band": {
    "arm": "arm2",
    "cutoff_a": 6.543210,
    "cutoff_b": 9.876543,
    "has_band": true
  },
  "versions": {
    "python": "3.12.10",
    "numpy": "2.4.4",
    "pandas": "2.3.3",
    "scikit_learn": "1.8.0",
    "biopython": "1.88"
  }
}
```

Key by key:

| Key | Type | Meaning |
| --- | --- | --- |
| `schema_version` | int | Always 1. |
| `created_at` | str | Local time with UTC offset: `datetime.now().astimezone().isoformat(timespec="seconds")`. |
| `code_commit` | str | Output of `git rev-parse HEAD` at freeze time. |
| `code_dirty` | bool | `True` if `git status --porcelain` printed anything other than lines for `outputs/frozen.json`. |
| `seed` | int | 42. |
| `model_id` | str | Copied from the `model_id` column of `esm_scores.csv`. |
| `inputs.esm_scores_sha256` | str | SHA-256 of the bytes of `outputs/esm_scores.csv`, lowercase hex. |
| `inputs.esm_scores_rows` | int | 6194. |
| `inputs.slim_csv_sha256` | str | SHA-256 of the slim CSV. |
| `labels.*` | | The class rule, recorded so the file explains itself. |
| `tuning_counts` | object | Start-region counts used for tuning. Must equal 131 / 385 / 59. |
| `catch_target_percent` | int | 90. |
| `false_alarm_target_percent` | int | 10. |
| `arms.arm0.cutoff` | null | Arm 0 flags nothing. |
| `arms.arm1.cutoff` | float | Arm 1's frozen cutoff. Verified value: 0.0 (written as positive zero). |
| `arms.arm1.tuning_caught` | int | Start-region harmful flagged at that cutoff. Verified: 120. |
| `arms.arm1.tuning_false_alarms` | int | Start-region harmless flagged. Verified: 217. |
| `arms.arm2.cutoff` | float | Arm 2's frozen cutoff, on the harm scale (minus the ESM score). Unknown until run. |
| `arms.arm3.features` | list of str | Feature names in the fixed order. Must equal the list in section 7.8. |
| `arms.arm3.scaler_mean` | list of 8 floats | `StandardScaler.mean_` from the final fit, in feature order. |
| `arms.arm3.scaler_scale` | list of 8 floats | `StandardScaler.scale_` from the final fit. |
| `arms.arm3.coefficients` | list of 8 floats | `LogisticRegression.coef_[0]` from the final fit. |
| `arms.arm3.intercept` | float | `LogisticRegression.intercept_[0]`. |
| `arms.arm3.cutoff` | float | A probability between 0 and 1, chosen on out-of-fold probabilities. |
| `arms.arm3.converged` | bool | `True` if the final fit used fewer iterations than `max_iter`. |
| `arms.arm3.out_of_fold_*` | int | Counts at the cutoff on the out-of-fold probabilities. |
| `arms.arm3.tuning_*` | int | Counts at the same cutoff on the refit model's in-sample probabilities. Informational. |
| `band.cutoff_a` | float | Equal to `arms.arm2.cutoff`. |
| `band.cutoff_b` | float or null | Cutoff B, or `null` if no cutoff qualifies. |
| `band.has_band` | bool | `True` only if `cutoff_b` is not null and `cutoff_b > cutoff_a`. |
| `versions` | object | Library versions at freeze time. |

Write every float in full precision (`json.dump` default). Do not round.

### 6.6 Files in `outputs/results/` (written by stage 5, except where noted)

All percentages are stored as numbers from 0 to 100, rounded to 1 decimal.

**`headline.csv`** — four rows, one per arm, in order arm0 to arm3. End region, LOF and FUNC rows only (1,166 rows).

| Column | Type | Meaning |
| --- | --- | --- |
| `arm` | str | `arm0`, `arm1`, `arm2`, `arm3` |
| `arm_name` | str | The `name` from `frozen.json` |
| `cutoff` | float or empty | Frozen cutoff; empty for arm0 |
| `n_harmful` | int | 267 |
| `n_harmless` | int | 899 |
| `harmful_caught` | int | Harmful rows flagged |
| `harmful_missed` | int | `n_harmful − harmful_caught` |
| `false_alarms` | int | Harmless rows flagged |
| `harmless_cleared` | int | `n_harmless − false_alarms` |
| `catch_rate_percent` | float | `100 × harmful_caught / n_harmful` |
| `false_alarm_rate_percent` | float | `100 × false_alarms / n_harmless` |
| `accuracy_percent` | float | `100 × (harmful_caught + harmless_cleared) / 1166` |

The arm0 row is known in advance: `arm0,Always harmless,,267,899,0,267,0,899,0.0,0.0,77.1`.

**`headline.md`** — the same table as a Markdown table, for pasting into the write-up, with one extra last column `AUROC` (each arm's value from `auroc.csv`, shown with 3 decimals), so the table carries its own ranking.

**`auroc.csv`** — columns `region_tested` (`end`), `arm`, `n_harmful`, `n_harmless`, `auroc` (float, 4 decimals). Four rows. Arm 0's `auroc` is 0.5 by definition (a constant score).

**`bootstrap.json`**

```json
{
  "unit": "end-region positions",
  "n_positions": 225,
  "n_resamples_requested": 1000,
  "n_resamples_used": 1000,
  "n_resamples_skipped_one_class": 0,
  "seed": 42,
  "esm_minus_blosum": {"point": 0.0, "low": 0.0, "high": 0.0},
  "arm3_minus_esm": {"point": 0.0, "low": 0.0, "high": 0.0}
}
```

`point` is the difference in AUROC on the real end-region rows. `low` and `high` are the 2.5th and 97.5th percentiles of the resampled differences. `n_positions` is the number of distinct end-region positions that have at least one LOF or FUNC row (at most 225; the exact value is unknown until the run).

**`band_counts.csv`** — arm 2 only, end region. Columns `findlay_class` (`LOF`, `FUNC`, `INT`), `looks_harmful`, `cant_tell`, `looks_harmless`, `total` (all int). Three rows. If `has_band` is false, `cant_tell` is 0 in every row.

**`by_swap.csv`** — end region, LOF and FUNC rows. One row per arm (arm1, arm2, arm3) × swap kind (`to_proline`, `loss_of_cysteine`, `charge_change`, `other`): 12 rows. Columns `arm`, `swap_kind`, `n_harmful`, `n_harmless`, `harmful_caught`, `harmful_missed`, `false_alarms`, `catch_rate_percent`, `false_alarm_rate_percent`. A rate whose denominator is 0 is written as an empty cell.

**`int_rows.csv`** — what each arm says about the INT rows. Columns `region` (`start`, `end`), `arm` (arm1, arm2, arm3), `n_int` (59 or 96), `flagged`, `not_flagged`, `flagged_percent`. Six rows.

**`reverse_split.csv`** — secondary. Tuned on the end region, tested on the start region. Same columns as `headline.csv` plus a first column `direction` whose value is always `SECONDARY: tuned on end, tested on start`, and a last column `auroc`. Four rows.

**`predictions_end.csv`** — one row per lab-tested end-region mutation (1,262 rows, INT included).

| Column | Type | Meaning |
| --- | --- | --- |
| `mutant`, `position`, `wild_type`, `new_letter` | | as in the slim table |
| `function_score`, `findlay_class` | | the lab's answer |
| `is_harmful` | int or empty | 1 for LOF, 0 for FUNC, empty for INT |
| `swap_kind` | str | section 7.9 |
| `blosum_score` | int | raw BLOSUM62 value |
| `esm_score`, `entropy` | float | from the cache |
| `harm_arm1`, `harm_arm2`, `harm_arm3` | float | harm scores |
| `flag_arm1`, `flag_arm2`, `flag_arm3` | int | 1 if harm ≥ that arm's frozen cutoff |
| `band_arm2` | str | `Looks harmful`, `Can't tell` or `Looks harmless` |

**`result_sentences.txt`** — the fixed sentences that apply, one per line, each prefixed with its code (`S1:` and so on), with X and Y filled in from the tables.

**`sanity7.json`** — sanity test 7's report: `{"M1775R_esm_score": float, "A1708E_esm_score": float, "spearman_all_1837": float, "spearman_start_575": float, "spearman_end_1262": float}`.

**`run_log.txt`** — appended, never edited. One `refused` line for an attempt that failed the checks in step 2 of the test run (a step 1 refusal, a repeat without `--bugfix`, writes nothing); otherwise a `started` line followed by a `success` or `failed` line. Format, with fields separated by ` | ` (space, pipe, space):

```
2026-11-15T09:05:10-06:00 | run=1 | status=started | kind=first | reason=first run | commit=0123abc | esm_sha256=…64 hex… | frozen_sha256=…64 hex…
2026-11-15T09:05:12-06:00 | run=1 | status=success | kind=first | reason=first run | commit=0123abc | esm_sha256=…64 hex… | frozen_sha256=…64 hex…
```

A normal first run therefore leaves exactly two lines: `started`, then `success`.

| Field | Values |
| --- | --- |
| timestamp | local time with offset, to the second |
| `run=` | The run number: the count of earlier `status=started` lines, plus 1. A `started` line and the `success` or `failed` line that closes it carry the same number. A `refused` line carries the number the run would have had. |
| `status=` | `refused` (stopped before loading any end-region label); `started` (written immediately before end-region labels are loaded); `success` (all result files written); `failed` (an error or interruption after `started`) |
| `kind=` | `first` or `bugfix` |
| `reason=` | `first run`, or the text given to `--bugfix`, or for `refused` and `failed` the error message. Replace any `|` or newline in it with a space. |
| `commit=` | first 7 characters of `git rev-parse HEAD` |
| `esm_sha256=`, `frozen_sha256=` | hashes of the two input files as found on disk |

**`failures.csv` and `failures.md`** (stage 6) — ten rows: the five worst misses and the five worst false alarms of arm 2. Columns `kind` (`miss` or `false_alarm`), `rank` (1–5), `mutant`, `position`, `wild_type`, `new_letter`, `function_score`, `findlay_class`, `esm_score`, `harm_arm2`, `band_arm2`, `blosum_score`, `entropy`, `swap_kind`, `flag_arm1`, `flag_arm3`, `guess_why`. `guess_why` is free text that always begins with the word `Guess:`.

### 6.7 `app/app_data.csv` and `app/app_config.json` (written by stage 8)

The demo reads only these two files. Both are built by `src/export_app.py`.

**`app/app_data.csv`** — 6,194 rows, same order as `esm_scores.csv`.

| Column | Type | Meaning | Nulls |
| --- | --- | --- | --- |
| `position` | int | 1–101, 1631–1855 | none |
| `wild_type` | str | normal letter | none |
| `new_letter` | str | swapped-in letter | none |
| `esm_score` | float | 6 decimals | none |
| `lab_class` | str | `LOF`, `INT`, `FUNC`; empty if this swap was not lab-tested; **also empty for every end-region row until the test run is done** (and until then the app refuses every end-region question, section 8.13) | allowed |
| `lab_score` | float | `function_score` rounded to 3 decimals; empty under the same conditions | allowed |

Example row: `61,C,G,-9.123456,LOF,-1.738` (the `esm_score` is made up; the rest is verified).

**`app/app_config.json`**

```json
{
  "model_id": "facebook/esm2_t33_650M_UR50D",
  "cutoff_a": 6.543210,
  "cutoff_b": 9.876543,
  "has_band": true,
  "end_labels_included": false,
  "protein_length": 1863,
  "n_tested_positions": 326,
  "exam_counts": null
}
```

`end_labels_included` is the value `test_run_done()` returned when the export ran.

`exam_counts` is `null` until the test run is done. After it, it is an object with one key per verdict, each holding how many lab-tested end-region mutations got that verdict from arm 2, split by the lab's class:

```json
"exam_counts": {
  "Looks harmful":  {"LOF": 0, "FUNC": 0, "INT": 0},
  "Can't tell":     {"LOF": 0, "FUNC": 0, "INT": 0},
  "Looks harmless": {"LOF": 0, "FUNC": 0, "INT": 0}
}
```

(The zeros are placeholders.) These nine numbers must equal `outputs/results/band_counts.csv`. They let the demo say how often each verdict matched the lab on the region the cutoffs never saw.

---

## 7. Constants

All constants live in `src/config.py`. No other module may hard-code any of these values. `app/app.py` is the one exception: it cannot import `src`, so it keeps its own copies of the few it needs, each with a comment naming the constant it copies.

### 7.1 Regions, classes, targets

| Name | Value | Notes |
| --- | --- | --- |
| `SEED` | `42` | |
| `PROTEIN_LENGTH` | `1863` | |
| `START_REGION` | `(1, 101)` | inclusive, 1-based |
| `END_REGION` | `(1631, 1855)` | inclusive, 1-based |
| `LOF_BELOW` | `-1.328` | `LOF` if `score < -1.328` (strict) |
| `FUNC_ABOVE` | `-0.748` | `FUNC` if `score > -0.748` (strict) |
| | | `INT` if `-1.328 <= score <= -0.748` |
| `HARMFUL_CLASS` | `"LOF"` | |
| `HARMLESS_CLASS` | `"FUNC"` | |
| `EXCLUDED_CLASS` | `"INT"` | |
| `CATCH_TARGET_PERCENT` | `90` | an integer, so the comparison is exact |
| `FALSE_ALARM_TARGET_PERCENT` | `10` | an integer |
| `N_FOLDS` | `5` | |
| `N_BOOTSTRAP` | `1000` | |
| `BOOTSTRAP_PERCENTILES` | `(2.5, 97.5)` | |
| `N_ROWS` | `1837` | |
| `N_POSITIONS` | `326` | |
| `N_SCORE_ROWS` | `6194` | 326 × 19 |
| `EXPECTED_COUNTS` | `{"start": {"LOF": 131, "INT": 59, "FUNC": 385}, "end": {"LOF": 267, "INT": 96, "FUNC": 899}}` | |
| `SLIM_CSV_SHA256` | `"4eb3d5597a4a2fcad88fa371b55cee9d8028c6df1ea18b95ad23802484e1a0ac"` | |
| `RAW_CSV_SHA256` | `"abde5d47e84b1feeac6c32b60c1f3b5d2881accf3ceccfbadccceeea73e81c43"` | |
| `FASTA_SHA256` | `"47529e88620d4d071f3ebec5212a399d74e102bc48f1ef5a1e1215186a196c36"` | |
| `SEQUENCE_START` | `"MDLSALRVEEVQNVINAMQK"` | |
| `SPEARMAN_MIN_START` | `0.15` | sanity test 5 fails below this |
| `BARELY_BEATS_AUROC` | `0.05` | used only to choose sentence S4 |
| `CATCH_DRIFT_POINTS` | `5.0` | used only to choose sentence S5 |

### 7.2 Models

| Name | Value |
| --- | --- |
| `MODEL_ID` | `"facebook/esm2_t33_650M_UR50D"` |
| `FALLBACK_MODEL_IDS` | `("facebook/esm2_t30_150M_UR50D", "facebook/esm2_t12_35M_UR50D")` |
| `ESM_BATCH_SIZE` | `8` |

The fallback order is 650M, then 150M, then 35M. A fallback is used only after a human says so. The default is always `MODEL_ID`.

### 7.3 Windows and the token-index formula

The protein is longer than the model's comfortable input, so the model sees one of two 500-letter windows.

| Window | Covers region | Residues (1-based, inclusive) | Python slice | `WINDOW_START` (1-based) |
| --- | --- | --- | --- | --- |
| `A` | start | 1–500 | `seq[0:500]` | 1 |
| `B` | end | 1364–1863 | `seq[1363:1863]` | 1364 |

`WINDOWS = {"A": (1, 500), "B": (1364, 1863)}`.

The tokenizer puts one start token at index 0, then one token per letter, then one end token. So a window of 500 letters becomes 502 tokens, and:

```
window          = "A" if 1 <= position <= 101, "B" if 1631 <= position <= 1855
token_index     = position - WINDOW_START[window] + 1
```

- Window A: `token_index = position`. Positions 1–101 map to token indices 1–101.
- Window B: `token_index = position − 1363`. Positions 1631–1855 map to token indices 268–492.

### 7.4 The 20 letters, in one fixed order

`AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"`

Use this order everywhere: heat map rows, the 20 token IDs, the entropy sum, sorting. Never use `sorted(set(...))` or dictionary order instead.

Three-letter codes, for the demo's input parser:

`THREE_TO_ONE = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"}`

### 7.5 Kyte-Doolittle hydrophobicity

Source: Kyte J, Doolittle RF (1982), "A simple method for displaying the hydropathic character of a protein", Journal of Molecular Biology 157:105–132. No units. Higher means more water-avoiding. Name this source in a code comment and in `writeup/decision_log.md`.

`HYDROPHOBICITY`: `{"A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3}`

### 7.6 Amino-acid size

Source: Zamyatnin AA (1972), "Protein volume in solution", Progress in Biophysics and Molecular Biology 24:107–123. Units: cubic ångströms (Å³), the volume of the residue. Name this source in a code comment and in `writeup/decision_log.md`.

`VOLUME`: `{"A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9, "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7, "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4, "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6}`

### 7.7 Charge

`POSITIVE_LETTERS = "KR"`, `NEGATIVE_LETTERS = "DE"`. Every other letter is neutral. **Histidine (`H`) counts as neutral**, because it is mostly uncharged at the pH inside a cell. `charge(letter)` is `+1`, `−1` or `0`.

### 7.8 Arm 3 feature order

`FEATURE_ORDER = ["esm_score", "blosum_score", "entropy", "size_change", "hydrophobicity_change", "to_proline", "from_cysteine", "from_glycine"]`

Formulas are in section 8.6.

### 7.9 Swap kinds

`SWAP_KINDS = ["to_proline", "loss_of_cysteine", "charge_change", "other"]`

Each swap gets exactly one kind. Check in this order and take the first that matches:

1. `to_proline`: `new_letter == "P"`.
2. `loss_of_cysteine`: `wild_type == "C"`.
3. `charge_change`: `charge(wild_type) != charge(new_letter)`. This covers gaining a charge, losing one, and reversing one.
4. `other`: everything else.

Examples, written as letter pairs (normal letter → new letter): L→P is `to_proline`. C→G is `loss_of_cysteine`. C→R is `loss_of_cysteine` (rule 2 comes before rule 3). C→P is `to_proline` (rule 1 comes first). S→D is `charge_change`. L→K is `charge_change`. K→R is `other` (both positive). H→R is `charge_change` (H is neutral, R is positive). L→S is `other`.

### 7.10 Zinc letters

`ZINC_POSITIONS = [24, 27, 39, 41, 44, 47, 61, 64]`

The letters at these positions are C, C, C, H, C, C, C, C (verified in our sequence). They are the eight zinc-holding letters of the RING region, **confirmed** on Nov 15 against the solved structure (Protein Data Bank entry 1JM7, the BRCA1/BARD1 RING-domain heterodimer): one zinc is held by C24, C27, C44 and C47, the other by C39, H41, C61 and C64. The heat map marks them with the label "zinc-holding letters (PDB 1JM7)". The write-up, the video and the `guess_why` text may say that a position is one of these eight letters, which is a fact. Whether the model's scores stand out at these positions is a result: say it only if `outputs/heatmap.png` shows it, and describe what the picture shows, nothing more.

### 7.11 Heat map

| Name | Value |
| --- | --- |
| `HEATMAP_ZMIN` | `-15.0` |
| `HEATMAP_ZMAX` | `0.0` |
| `HEATMAP_COLORSCALE` | `"Reds_r"` (dark red = strongly negative = the model finds the swap very surprising; near white = zero) |

The colour scale is one-sided, because ESM scores are mostly negative. Scores below −15 are drawn in the darkest colour and scores above 0 in the lightest. This range is fixed here and must not be adjusted after looking at the data.

### 7.12 Paths

Build every path from `ROOT = Path(__file__).resolve().parent.parent` (the repository root).

| Name | Path |
| --- | --- |
| `DATA_DIR` | `ROOT / "data"` |
| `RAW_CSV_PATH` | `DATA_DIR / "BRCA1_HUMAN_Findlay_2018.csv"` |
| `SLIM_CSV_PATH` | `DATA_DIR / "brca1_missense_scores_slim.csv"` |
| `FASTA_PATH` | `DATA_DIR / "BRCA1_HUMAN_Findlay_2018_wildtype.fasta"` |
| `OUTPUTS_DIR` | `ROOT / "outputs"` |
| `ESM_SCORES_PATH` | `OUTPUTS_DIR / "esm_scores.csv"` |
| `FROZEN_PATH` | `OUTPUTS_DIR / "frozen.json"` |
| `RESULTS_DIR` | `OUTPUTS_DIR / "results"` |
| `RUN_LOG_PATH` | `RESULTS_DIR / "run_log.txt"` |
| `SELFCHECK_DIR` | `OUTPUTS_DIR / "selfcheck_start"` |
| `HEATMAP_HTML_PATH` | `OUTPUTS_DIR / "heatmap.html"` |
| `HEATMAP_PNG_PATH` | `OUTPUTS_DIR / "heatmap.png"` |
| `APP_DIR` | `ROOT / "app"` |
| `APP_DATA_PATH` | `APP_DIR / "app_data.csv"` |
| `APP_CONFIG_PATH` | `APP_DIR / "app_config.json"` |

### 7.13 Fixed text

| Name | Value |
| --- | --- |
| `DISCLAIMER` | `"This project is a prototype and is not medical advice."` |
| `DEMO_NOTE` | `"At lab-tested positions a counselor would use the lab result. This demo shows how the model does where we can check it."` |
| `VERDICT_HARMFUL` | `"Looks harmful"` |
| `VERDICT_CANT_TELL` | `"Can't tell"` |
| `VERDICT_HARMLESS` | `"Looks harmless"` |
| `ARM_NAMES` | `{"arm0": "Always harmless", "arm1": "BLOSUM62 rule", "arm2": "ESM-2 zero-shot", "arm3": "Logistic regression on ESM-2 + BLOSUM62 + chemistry"}` |

`"Can't tell"` uses a plain ASCII apostrophe.

---

## 8. Module specs

### 8.0 Call graph

An arrow means "imports". There are no cycles. Nothing imports `evaluate`, `freeze`, `failures`, `heatmap`, `export_app` or `esm_score`; they are entry points.

```
config        (imports nothing from src)
load          -> config
cutoffs       -> config
baselines     -> config
features      -> config, baselines, cutoffs
esm_score     -> config, load
freeze        -> config, load, baselines, cutoffs, features
evaluate      -> config, load, baselines, cutoffs, features
failures      -> config, load
heatmap       -> config, load
export_app    -> config, load, cutoffs
app/app.py    -> nothing from src (it must run with only the app/ folder present)
tools/build_notebook.py -> nothing from src (it only writes notebook cells as text)
```

Rules that keep imports cheap and safe:

- `torch` and `transformers` are imported **inside functions** in `src/esm_score.py`, never at the top of any module. Every other module must import without torch installed.
- `Bio` is imported only in `src/baselines.py`.
- `plotly` and `matplotlib` are imported only in `src/heatmap.py`.
- `gradio` is imported only in `app/app.py`.
- **Only `src/load.py` may read** `data/*.csv`, `data/*.fasta` and `outputs/esm_scores.csv`. Every other module gets its data by calling functions in `load`. This makes the no-peeking rule checkable: a reviewer can search for `read_csv` and expect to find it only in `load.py`, in `failures.py` (which reads `predictions_end.csv` after the test) and in `app/app.py` (which reads `app_data.csv`). `heatmap.py` and `export_app.py` never read a file in `outputs/results/`.

Every entry-point module has a function `main(argv: list[str] | None = None) -> int` that parses `argv` with `argparse` (`None` means the real command line), calls `config.set_seeds()` first, and returns 0 on success and 1 on any refusal or failure. The module ends with `if __name__ == "__main__": raise SystemExit(main())`. Refusals print the exact message given in this spec to standard error. Tests call `main([...])` directly and check the returned code.

**Use paths through the module, not by name.** Write `from src import config` and then `config.FROZEN_PATH`, `config.RUN_LOG_PATH` and so on at the point of use. Do not write `from src.config import FROZEN_PATH`, and do not use a path constant as a function's default argument. Otherwise the unit tests cannot point the code at a temporary folder, and a test could overwrite a real file.

### 8.1 `src/config.py`

**Purpose.** Hold every constant and path from section 7, plus two tiny helpers.

**Public names.** Everything listed in section 7, with exactly those names.

**Functions.**

`set_seeds(seed: int = SEED) -> None`
1. `random.seed(seed)`.
2. `numpy.random.seed(seed)`.
3. Try to import torch. If it imports: `torch.manual_seed(seed)` and `torch.cuda.manual_seed_all(seed)`. If the import fails, do nothing (most stages do not need torch).

`sha256_of_file(path: Path) -> str`
1. Open the file in binary mode (`"rb"`).
2. Read it in chunks of 1 MB and feed each to `hashlib.sha256()`.
3. Return the lowercase hex digest (64 characters).
4. If the file does not exist, raise `FileNotFoundError(f"Cannot hash {path}: file not found.")`.

**Must not.** Import anything from `src`. Read any data file. Import torch at the top of the module.

### 8.2 `src/load.py`

**Purpose.** The single door to the data. It reads the input files, checks them, and hands out labels only to callers that are allowed to have them.

**Design.** Labels are separated from everything else. A function that needs no labels calls `load_mutations()` and cannot leak anything. A function that needs labels calls `load_labels(region, purpose)`, which refuses to give end-region labels before the test run unless the caller is the test run itself.

#### Public functions

`parse_mutant(text: str) -> tuple[str, int, str]`

Parses the strict form used in the data files, such as `C61G`. (The demo has its own, more forgiving parser in `app/app.py`.)
1. If `text` is not a `str`, raise `TypeError("mutant must be a string")`.
2. Match the whole string against the regular expression `^([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]{0,3})([ACDEFGHIKLMNPQRSTVWY])$`. No stripping, no case change.
3. If it does not match, raise `ValueError(f"Cannot parse mutant '{text}': expected one letter, a position, one letter, for example C61G.")`.
4. If the two letters are equal, raise `ValueError(f"Mutant '{text}' does not change the letter.")`.
5. Return `(wild_type, int(position), new_letter)`.

`region_of(position: int) -> str`
1. If `1 <= position <= 101`, return `"start"`.
2. If `1631 <= position <= 1855`, return `"end"`.
3. Otherwise raise `ValueError(f"Position {position} is outside the two lab-tested regions (1-101 and 1631-1855).")`.

`class_of(function_score: float) -> str`
1. If `function_score < -1.328`, return `"LOF"`.
2. If `function_score > -0.748`, return `"FUNC"`.
3. Otherwise return `"INT"`.

`load_wild_type() -> str`
1. Read `FASTA_PATH` as UTF-8 text.
2. Drop lines starting with `>`. Strip the others. Join them.
3. Check: length is 1,863; it starts with `MDLSALRVEEVQNVINAMQK`; every character is in `AA_ORDER`. On failure raise `ValueError("Wild-type sequence check failed: <which check>.")`.
4. Return the string.

`load_mutations() -> pandas.DataFrame`

Label-free. Safe to call anywhere, any time.
1. Read `SLIM_CSV_PATH` with the dtypes in section 6.2.
2. Return only the columns `mutant`, `wild_type`, `position`, `new_letter`, `region` (1,837 rows, original order, a fresh copy).

`class_counts() -> dict[str, dict[str, int]]`

Counts only. Safe to call anywhere, any time.
1. Read the slim CSV.
2. Return `{"start": {"LOF": n, "INT": n, "FUNC": n}, "end": {...}}`.
3. Do not return, print or log any row.

`test_run_done(run_log_path: Path | None = None) -> bool`

Exactly the four steps in section 3, rule 15. If `run_log_path` is `None`, use `config.RUN_LOG_PATH`, looked up at call time.

`test_attempted(run_log_path: Path | None = None) -> bool`

Same reading steps, but returns `True` if any line has a part exactly equal to `status=started`, `status=success` or `status=failed`. `status=refused` lines do not count. Used by `freeze`, `esm_score` and `evaluate` to lock everything once the test has begun.

`load_labels(region: str, purpose: str) -> pandas.DataFrame`

Returns `mutant`, `wild_type`, `position`, `new_letter`, `region`, `function_score`, `findlay_class` for one region (`proteingym_bin` is never returned; no arm uses it).
1. If `region` is not `"start"` or `"end"`, raise `ValueError("region must be 'start' or 'end'.")`.
2. If `purpose` is not one of `"tuning"`, `"evaluate"`, `"post_test"`, raise `ValueError("purpose must be 'tuning', 'evaluate' or 'post_test'.")`.
3. Permission table:

| region | purpose | Allowed? |
| --- | --- | --- |
| `start` | any | always |
| `end` | `tuning` | never. Raise `PermissionError("End-region lab labels may not be used for tuning.")` |
| `end` | `evaluate` | only if `test_attempted()` is `True` (the `started` line has been written); otherwise raise `PermissionError("End-region lab labels may be loaded only inside the test run, after its started line is written.")`. **Only `src/evaluate.py` may pass this value, only inside the real test run, and at exactly one place in the code, written with the literal strings: `load.load_labels("end", "evaluate")`.** Never pass the region as a variable to a helper. |
| `end` | `post_test` | only if `test_run_done()` is `True`. Otherwise raise `PermissionError("End-region lab labels are locked until the one test run has succeeded (outputs/results/run_log.txt has no status=success line).")` |

4. Read the slim CSV, keep the rows of that region, return a fresh copy.

This is a tripwire against accidents, not a security system. Do not work around it. If it stops you, you are about to break rule 9.

`load_scores() -> pandas.DataFrame`

Reads the ESM score cache. Label-free.
1. If `ESM_SCORES_PATH` does not exist, raise `FileNotFoundError("outputs/esm_scores.csv not found. Run: python -m src.esm_score")`.
2. Read with dtypes: `position` int64, `token_index` int64, `is_lab_tested` int64, `esm_score` float64, `entropy` float64, the rest `str`. Pass `keep_default_na=False` for the string columns so that no letter is misread as missing.
3. Call `validate_scores(df)` and return `df`.

`validate_scores(df: pandas.DataFrame) -> None`

Raises `ValueError("esm_scores.csv check failed: <which check>.")` unless all of these hold:
1. Columns are exactly those in section 6.4, in that order.
2. 6,194 rows; `mutant` unique; 326 distinct positions; every position has exactly 19 rows.
3. No nulls. `esm_score` and `entropy` are finite.
4. `token_at_index == wild_type` in every row.
5. `wild_type` equals the wild-type sequence letter at `position` in every row.
6. `window` and `token_index` equal what section 7.3's formula gives.
7. `is_lab_tested` sums to 1,837, and the set of `mutant` values with `is_lab_tested == 1` equals the set of `mutant` values in the slim table.
8. `entropy` is constant within each position and lies in `[0, 2.9958]`.
9. `model_id`, `transformers_version` and `precision` each have exactly one distinct value.

`lab_tested_scores(region: str | None = None) -> pandas.DataFrame`

Label-free. `load_scores()` filtered to `is_lab_tested == 1`, and to one region if `region` is given. 575 rows for `start`, 1,262 for `end`, 1,837 for `None`.

`check_inputs() -> None`

Runs every check below and raises `ValueError("Input check failed: <which check>.")` on the first failure. Prints one line per passed check.
1. `sha256_of_file` of the raw CSV, the slim CSV and the FASTA equal the three constants in section 7.1.
2. The slim CSV has 1,837 rows, the eight columns of section 6.2 in order, no nulls, unique `mutant`.
3. For every row (on failure, report only which check failed and in how many rows, never the rows themselves): `parse_mutant(mutant)` equals `(wild_type, position, new_letter)`; `region == region_of(position)`; `findlay_class == class_of(function_score)`; `wild_type` equals the sequence letter at `position`.
4. `class_counts()` equals `EXPECTED_COUNTS`.
5. 326 distinct positions.
6. No `function_score` equals −1.328 or −0.748.

`rebuild_and_compare() -> None`

Proves the slim CSV and FASTA really come from the raw file. It never writes into `data/`.
1. Read `RAW_CSV_PATH` with dtypes `mutant` str, `mutated_sequence` str, `DMS_score` float64, `DMS_score_bin` int64.
2. Check: columns are exactly `mutant`, `mutated_sequence`, `DMS_score`, `DMS_score_bin`; 1,837 rows; no nulls; `mutant` unique.
3. Parse every `mutant` with `parse_mutant`.
4. Rebuild the wild type: take the first row's `mutated_sequence`, and put that row's wild-type letter back at its position (`index = position - 1`).
5. Check that **every** row agrees with it: for each row, `mutated_sequence` must equal the wild type with exactly `new_letter` at `position`, and the wild type's letter at `position` must equal the row's wild-type letter. If any row disagrees, raise `ValueError(f"Row {mutant} disagrees with the rebuilt wild-type sequence.")`.
6. Build the slim frame in memory: the eight columns of section 6.2, same row order.
7. Compare with the files on disk: the rebuilt sequence must equal `load_wild_type()`; the rebuilt frame must equal the slim CSV on every column (compare column by column with exact equality). If they differ, do **not** let pandas print the differing values, because they are lab labels: raise `ValueError(f"Rebuild check failed: column {name} differs in {n} rows.")`, naming only the column and the number of rows.
8. Print `REBUILD MATCHES: slim CSV and FASTA agree with the raw ProteinGym file.`

#### Command line

`python -m src.load --check` runs `check_inputs()` then prints the class-count table (counts only).
`python -m src.load --rebuild-check` runs `rebuild_and_compare()`.

#### Must not

- Write anything into `data/`.
- Print or log any row that contains `function_score` or `findlay_class` for the end region.
- Offer any function that returns end-region labels without the permission check.

### 8.3 `src/cutoffs.py`

**Purpose.** The one place cutoffs and bands are decided. Pure functions on numbers. No file reading.

All functions take NumPy arrays or anything `numpy.asarray` accepts. `harm_scores` is a 1-D float array. `is_harmful` is a 1-D array of 0 and 1 (or `False` and `True`) of the same length. Rows that are neither harmful nor harmless (INT) must be removed by the caller before calling. Each function begins by checking: same length, at least one harmful row, at least one harmless row, no NaN; otherwise it raises `ValueError` with a message naming the failed check.

`flag(harm_scores, cutoff: float) -> numpy.ndarray`

Returns a boolean array: `harm_scores >= cutoff`. This is the only place in the project where the flagging comparison is written. Everything else calls it.

`choose_catch_cutoff(harm_scores, is_harmful, target_percent: int = 90) -> float`

Returns the **highest** cutoff that flags at least `target_percent` % of the harmful rows.

1. Let `h` = the harm scores of the harmful rows, sorted from highest to lowest. Let `n` = `len(h)`.
2. Let `needed` = the smallest whole number `k` with `k * 100 >= target_percent * n`. Compute it with integers: `needed = -(-(target_percent * n) // 100)`. For `n = 131` and 90: `needed = 118`. For `n = 267` and 90: `needed = 241`.
3. Return `float(h[needed - 1])`: the `needed`-th highest harmful score.

Why this is right: any cutoff higher than `h[needed - 1]` flags fewer than `needed` harmful rows. A cutoff equal to it flags every harmful row with a score at least that high, which is `needed` or more. So it is the highest cutoff that reaches the target, and it is always an observed score. Because flagging is `>=`, rows with equal scores always get the same call.

Because of ties, the catch rate actually reached can be above the target. Report the rate reached, not the target.

**Worked example 1** (test `test_catch_cutoff_worked_example`). Six rows:

| Row | Harm score | Harmful? |
| --- | --- | --- |
| a | 0.9 | yes |
| b | 0.7 | no |
| c | 0.5 | yes |
| d | 0.5 | no |
| e | 0.5 | yes |
| f | 0.1 | no |

- Harmful scores sorted high to low: 0.9, 0.5, 0.5. `n = 3`.
- Target 90: `needed` = smallest `k` with `100k >= 270` = 3. Cutoff = third highest = **0.5**.
- Flagged at `>= 0.5`: a, b, c, d, e. Caught 3 of 3. False alarms 2 of 3 (b and d). Row d ties with c and e, so it gets the same call.
- Target 60 on the same rows: `needed` = smallest `k` with `100k >= 180` = 2. Cutoff = second highest = **0.5** again. The catch rate reached is 100%, not 60%, because the tie cannot be split.
- Target 30: `needed` = 1. Cutoff = **0.9**. Caught 1 of 3, false alarms 0.

`choose_false_alarm_cutoff(harm_scores, is_harmful, target_percent: int = 10) -> float | None`

Returns cutoff B: the **lowest** observed harm score `c` such that flagging at `>= c` flags at most `target_percent` % of the harmless rows. Returns `None` if no observed score qualifies.

1. Let `n_harmless` = number of harmless rows. Let `allowed = (target_percent * n_harmless) // 100` (integer division). For 385 and 10: `allowed = 38`.
2. Let `candidates` = the distinct values of `harm_scores` over **all** rows passed in (harmful and harmless), sorted from lowest to highest.
3. For each candidate `c` from lowest to highest: count harmless rows with score `>= c`. Return the first `c` whose count is `<= allowed`.
4. If none qualifies, return `None`.

**Worked example 2** (test `test_false_alarm_cutoff_worked_example`). Fourteen rows:

- Ten harmless rows with harm scores 8, 6, 5, 5, 4, 3, 3, 2, 1, 0.
- Four harmful rows with harm scores 9, 7, 6.5, 2.
- Target 10: `allowed = (10 × 10) // 100 = 1`.
- Candidates low to high: 0, 1, 2, 3, 4, 5, 6, 6.5, 7, 8, 9.
- Harmless rows flagged at each: 0→10, 1→9, 2→8, 3→7, 4→5, 5→4, 6→2, **6.5→1**. Stop.
- Cutoff B = **6.5**. It is a harmful row's score; that is allowed, because candidates come from all rows.
- Second case (test `test_false_alarm_cutoff_none`): harmless scores 5, 5, 5, 1 and harmful scores 4, 2, target 10. `allowed = 0`. The highest candidate, 5, still flags three harmless rows. Return `None`.

`has_band(cutoff_a: float, cutoff_b: float | None) -> bool`

Returns `True` only if `cutoff_b is not None and cutoff_b > cutoff_a`.

`band_verdict(harm_score: float, cutoff_a: float, cutoff_b: float | None) -> str`

1. If `has_band(cutoff_a, cutoff_b)`:
   - `harm_score >= cutoff_b` → `"Looks harmful"`
   - `harm_score < cutoff_a` → `"Looks harmless"`
   - otherwise → `"Can't tell"`
2. If there is no band:
   - `harm_score >= cutoff_a` → `"Looks harmful"`
   - otherwise → `"Looks harmless"`

Worked example (test `test_band_verdict`): A = 2.0, B = 5.0. Scores 5.0 → harmful; 7.1 → harmful; 4.99 → can't tell; 2.0 → can't tell; 1.99 → harmless. With B = `None` or B = 2.0 or B = 1.0 (no band): 2.0 → harmful; 1.99 → harmless.

`count_outcomes(harm_scores, is_harmful, cutoff: float | None) -> dict[str, int]`

Returns `{"n_harmful", "n_harmless", "harmful_caught", "harmful_missed", "false_alarms", "harmless_cleared"}`. If `cutoff` is `None` (arm 0), nothing is flagged.

`rates(counts: dict[str, int]) -> dict[str, float]`

Returns `catch_rate_percent`, `false_alarm_rate_percent`, `accuracy_percent`, each `round(100 * a / b, 1)`.

**Must not.** Read files. Know anything about regions or arms. Add noise or tie-breakers.

### 8.4 `src/baselines.py`

**Purpose.** Arm 0 and arm 1.

`load_blosum62() -> dict[tuple[str, str], int]`
1. `from Bio.Align import substitution_matrices`; `matrix = substitution_matrices.load("BLOSUM62")`. (Verified to work in Biopython 1.88. The matrix's alphabet is `ARNDCQEGHILKMFPSTWYVBZX*`, which is **not** `AA_ORDER`, so always index it by letter, never by number.)
2. Build a dictionary for all 20 × 20 pairs of letters in `AA_ORDER`: `table[(a, b)] = int(matrix[a, b])`. `matrix[a, b]` returns a float such as `-3.0`.
3. Checks, each raising `ValueError("BLOSUM62 check failed: <which check>.")`:
   - all 400 pairs are present;
   - every value is a whole number before `int()` (`float(v).is_integer()`);
   - the table is symmetric: `table[(a, b)] == table[(b, a)]`;
   - spot values: `("C", "G")` is −3, `("D", "E")` is 2, `("M", "I")` is 1, `("W", "W")` is 11;
   - over the 380 off-diagonal pairs, the minimum is −4 and the maximum is 3.
   (All verified.)
4. Cache the result in a module-level variable so the matrix loads once.

`blosum_score(wild_type: str, new_letter: str) -> int`

Returns `load_blosum62()[(wild_type, new_letter)]`. Raises `KeyError` with the pair if either letter is not one of the 20.

`harm_arm1(wild_type: str, new_letter: str) -> float`

Returns `0.0 - blosum_score(wild_type, new_letter)`. Write it as a subtraction from `0.0`, not as a unary minus: `-0.0` is a real floating-point value that prints as `-0.0` and looks wrong in `frozen.json`.

`add_arm1(df: pandas.DataFrame) -> pandas.DataFrame`

Returns a copy of `df` with two new columns: `blosum_score` (int) and `harm_arm1` (float), computed from the `wild_type` and `new_letter` columns.

`arm0_flags(n_rows: int) -> numpy.ndarray`

Returns `numpy.zeros(n_rows, dtype=bool)`. Arm 0 flags nothing and has no cutoff. Its harm score, where one is needed (AUROC), is the constant `0.0`.

**Command line.** `python -m src.baselines` loads start-region labels with `load.load_labels("start", "tuning")`, drops INT rows, computes arm 1's cutoff with `cutoffs.choose_catch_cutoff`, and prints the cutoff, the counts and the rates. It also prints arm 0's known end-region numbers from `load.class_counts()` (0 of 267 caught; right on 899 of 1,166; 77.1%). It writes no file.

**Must not.** Load end-region labels. Use any ESM score.

### 8.5 `src/esm_score.py`

**Purpose.** Run ESM-2 once over the 326 positions and write `outputs/esm_scores.csv`. This module uses no labels at all.

#### The index maths, with two worked examples

**Example 1: C61 (position 61, start region).**

| Step | Value |
| --- | --- |
| Protein position | 61 |
| Normal letter | `seq[60]` = `C` |
| Region → window | start → `A` |
| Window slice | `window_seq = seq[0:500]` (residues 1–500) |
| Residue number inside the window | `61 − 1 + 1 = 61` |
| Letter in the window string | `window_seq[60]` = `C` |
| Token index | 61 (index 0 is the start token) |
| Token at that index | `C` |

**Example 2: M1775 (position 1775, end region).**

| Step | Value |
| --- | --- |
| Protein position | 1775 |
| Normal letter | `seq[1774]` = `M` |
| Region → window | end → `B` |
| Window slice | `window_seq = seq[1363:1863]` (residues 1364–1863) |
| Residue number inside the window | `1775 − 1364 + 1 = 412` |
| Letter in the window string | `window_seq[411]` = `M` |
| Token index | 412 |
| Token at that index | `M` |

More verified pairs, for tests: position 1 → A, 1, `M`. Position 101 → A, 101, `Y`. Position 1631 → B, 268, `S`. Position 1855 → B, 492, `I`.

#### Facts about the tokenizer (verified with transformers 5.23.0 and `facebook/esm2_t33_650M_UR50D`)

- `AutoTokenizer.from_pretrained(MODEL_ID)` returns an `EsmTokenizer`.
- Vocabulary size 33. `<cls>` = 0 (start token), `<pad>` = 1, `<eos>` = 2 (end token), `<unk>` = 3, `<mask>` = 32.
- Letter IDs: L 4, A 5, G 6, V 7, S 8, E 9, R 10, T 11, I 12, D 13, P 14, K 15, Q 16, N 17, F 18, Y 19, M 20, H 21, W 22, C 23.
- `tokenizer(window_seq, return_tensors="pt")["input_ids"]` for a 500-letter window has shape `(1, 502)`: `<cls>`, 500 letters, `<eos>`.
- Window A's first six tokens: `<cls>`, `M`, `D`, `L`, `S`, `A`. Last three: `T`, `N`, `<eos>`.
- Window B's first six tokens: `<cls>`, `N`, `L`, `G`, `E`, `A`. Last three: `H`, `Y`, `<eos>`.
- Pass the window as a plain string with no spaces.

Do not hard-code any of these IDs. Look them up from the tokenizer and check them (the letter test below). They are listed so you know what to expect.

#### Public functions

`window_of(position: int) -> str`

`"A"` if `1 <= position <= 101`; `"B"` if `1631 <= position <= 1855`; otherwise raise `ValueError(f"Position {position} is not in a scored region.")`.

`token_index_of(position: int) -> int`

`position - WINDOW_START[window_of(position)] + 1`, where `WINDOW_START` is 1 for A and 1364 for B.

`window_sequence(seq: str, window: str) -> str`

`seq[0:500]` for `"A"`, `seq[1363:1863]` for `"B"`. Check the result has length 500.

`scored_positions() -> list[int]`

The 326 distinct positions in `load.load_mutations()`, sorted ascending. Check there are 326.

`letter_token_ids(tokenizer) -> list[int]`
1. For each letter in `AA_ORDER`, `tokenizer.convert_tokens_to_ids(letter)`.
2. Check: 20 distinct IDs; none equals `tokenizer.unk_token_id`; `tokenizer.convert_ids_to_tokens(id)` gives back the same letter for each.
3. On failure raise `ValueError("Letter token lookup failed: <detail>.")`.
4. Return the 20 IDs in `AA_ORDER`.

`run_letter_test(tokenizer, seq: str, positions: list[int]) -> pandas.DataFrame`

The label-free guard against off-by-one bugs. It needs the tokenizer only, not the model weights.
1. For each window `A` and `B`: tokenize `window_sequence(seq, window)`; get the list of tokens with `tokenizer.convert_ids_to_tokens`.
2. Check each window: exactly 502 tokens; token 0 is `tokenizer.cls_token`; token 501 is `tokenizer.eos_token`; for every `i` from 1 to 500, token `i` equals `window_seq[i - 1]`.
3. For each position in `positions`: compute `window_of`, `token_index_of`; read the token at that index; compare with `seq[position - 1]`.
4. Return a frame with columns `position`, `window`, `token_index`, `token_at_index`, `wild_type`, `ok` (bool), one row per position.
5. If any check in step 2 fails, or any `ok` is `False`, raise `RuntimeError(f"LETTER TEST FAILED at {n_bad} of {n} positions. First failure: position {p}, window {w}, token index {i}, token '{t}', expected '{x}'. Do not score anything. Stop and report.")`.
6. On success print `LETTER TEST PASSED: 326 of 326 positions, both windows 502 tokens.`

**If the letter test fails: stop.** Do not adjust indices by trial and error until it passes. Report the failure. The `index-debugger` review (section 11) comes next. Section 15 lists the eight checks that review makes.

`choose_device(allow_cpu: bool) -> str`
1. If `torch.cuda.is_available()`, return `"cuda"`.
2. Else if `allow_cpu`, print a warning that scoring will be slow and return `"cpu"`.
3. Else raise `RuntimeError("CUDA is not available. Stage 3 needs an NVIDIA GPU. See BUILD_SPEC.md section 5.5. To run on CPU anyway (slow), pass --allow-cpu after asking a human.")`.

`load_model(model_id: str, device: str)`
1. `from transformers import AutoTokenizer, EsmForMaskedLM`.
2. `tokenizer = AutoTokenizer.from_pretrained(model_id)`.
3. `model = EsmForMaskedLM.from_pretrained(model_id)`.
4. `model.eval()`, then `model.to(device)`. Full precision (float32). Do not call `.half()`.
5. Return `(tokenizer, model)`.

`score_positions(model, tokenizer, seq: str, positions: list[int], device: str, batch_size: int = ESM_BATCH_SIZE) -> pandas.DataFrame`

Returns the 6,194-row frame of section 6.4 (without `is_lab_tested`, `model_id`, `transformers_version`, `precision`, which `main` adds).

1. `ids20 = letter_token_ids(tokenizer)`. Look these up once.
2. `letter_test = run_letter_test(tokenizer, seq, positions)`. It raises on failure, so nothing is scored with a wrong index.
3. For each window in order A, B:
   1. `base_ids = tokenizer(window_seq, return_tensors="pt")["input_ids"]`, shape `(1, 502)`.
   2. Take this window's positions in ascending order. Process them in batches of `batch_size`.
   3. For a batch of `k` positions with token indices `i_1 … i_k`: build `batch = base_ids.repeat(k, 1)` (shape `(k, 502)`). For row `j`, **assert** `batch[j, i_j]` equals the token ID of the normal letter, then set `batch[j, i_j] = tokenizer.mask_token_id`. Each row has exactly one masked token.
   4. Move `batch` to `device`. Inside `with torch.no_grad():` run `logits = model(input_ids=batch).logits`, shape `(k, 502, 33)`. No attention mask is needed: every row has the same length and no padding.
   5. For row `j`, take `row_logits = logits[j, i_j]` — the logits **at the masked index**, the same index that was masked. Shape `(33,)`.
   6. `logp = torch.log_softmax(row_logits.float(), dim=-1)` over all 33 vocabulary entries.
   7. `logp20 = logp[ids20]`, the 20 letters' log-probabilities in `AA_ORDER`. Move to CPU as float64 NumPy.
   8. `wt_index = AA_ORDER.index(normal_letter)`.
   9. For each of the 19 other letters `b` (in `AA_ORDER`, skipping the normal letter): `esm_score = logp20[index of b] - logp20[wt_index]`.
   10. Entropy: `p = exp(logp20)`; `p = p / p.sum()` (rescale the 20 so they sum to 1); `entropy = -(p * log(p)).sum()` with natural log. Since softmax outputs are never exactly zero, no special case is needed; assert `0 <= entropy <= 2.9958`.
   11. Append 19 rows.
   12. Print progress every batch: `window A: 40/101 positions`.
4. Round `esm_score` and `entropy` to 6 decimals. Sort by `position`, then by `AA_ORDER` index of `new_letter`.
5. Return the frame.

**The sign.** `esm_score = log P(new letter) − log P(normal letter)`. The normal letter is usually what the model expects, so most scores are negative. A harmful swap should be strongly negative. The harm score for arm 2 is `0.0 − esm_score`, so it is usually positive. If most of your scores are positive, the subtraction is the wrong way round.

**Batching: yes.** Batch 8 masked copies of the same window per forward pass. All copies have the same length, so there is no padding and batching does not change the maths (results differ from unbatched only by floating-point noise, around 1e-5).

**Memory (estimate).** The 650M model in float32 holds about 2.6 GB of weights. A batch of 8 × 502 tokens with no gradients needs well under 2 GB more. The laptop GPU has 8 GB with about 1.4 GB used by the desktop, so batch 8 should fit.

**If the GPU runs out of memory**, in this order:
1. Retry with `--batch-size 1`. You may do this without asking.
2. Still out of memory: 🛑 stop and ask. The options a human may choose are half precision (`model.half()`, recorded as `precision = float16`), then the 150M model, then the 35M model. Say what each loses: half precision adds small rounding noise; smaller models are generally weaker predictors, and the write-up must name the model actually used.

**Runtime (estimate).** First run: downloading the weights (about 2.6 GB) takes most of the time. Scoring is 326 forward passes: roughly 1–3 minutes on the RTX 5060 or a T4, roughly 10–30 minutes on a CPU.

#### Command line

`python -m src.esm_score [--only MUTANT] [--letter-test-only] [--validate-cache] [--recompute] [--batch-size N] [--allow-cpu] [--model-id ID]`

| Flag | Behaviour |
| --- | --- |
| (none) | If the cache exists: validate it, print `Loaded cached scores: 6194 rows, model <id>.` and exit 0 without loading the model. If not: first check the lock: if `outputs/frozen.json` exists or `load.test_attempted()` is `True`, refuse with the `--recompute` message below and exit 1 (a missing cache at that point is a problem for a human, not something to regenerate silently). Otherwise run the letter test, score all 326 positions, write the cache, validate it, print its SHA-256. |
| `--letter-test-only` | Load the tokenizer only (no weights, no GPU), run the letter test on all 326 positions, print the result, exit. Writes nothing. |
| `--only C61G` | Load the model, run the letter test, score that one mutation (batch of one), print `C61G: position 61, window A, token index 61, token C, esm_score = <value>, entropy = <value>`. Writes nothing. **Refused for any end-region position while `load.test_run_done()` is `False`**, with exit code 1, before the model is loaded, and with: `Refusing to score a single end-region mutation before the test run. Use --letter-test-only to check end-region indexing.` Looking at the score of one famous end-region mutation is a check rule 9 forbids. (The M1775 worked example above is about token indices only.) |
| `--validate-cache` | `load.load_scores()` and print the row count, the model ID and the SHA-256. |
| `--recompute` | Score again and overwrite the cache. **Refused if `outputs/frozen.json` exists or `load.test_attempted()` is `True`**, with: `Refusing to recompute: the current scores are locked (outputs/frozen.json exists or the test has started). Ask a human.` |
| `--model-id` | Use a fallback model. Only after a human approves. |

Writing the cache:
1. Add `is_lab_tested` (1 if `mutant` is in `load.load_mutations()["mutant"]`, else 0), `model_id`, `transformers_version` (`transformers.__version__`), `precision` (`"float32"`).
2. Create `outputs/` if needed. Write to `outputs/esm_scores.csv.tmp` with `float_format="%.6f"`, close the file, then rename it to `outputs/esm_scores.csv` with `os.replace`. This guarantees the file is complete before anything hashes it.
3. Read it back with `load.load_scores()` (which validates), then print `sha256_of_file`.

#### Must not

- Read or import any label. This module never calls `load.load_labels`.
- Compute anything from end-region labels. (Scoring end-region positions is fine and required.)
- Print a summary of scores by region that could be compared with labels by eye; print row counts and the C61G line only.
- Change the cache once `frozen.json` exists.

### 8.6 `src/features.py`

**Purpose.** Build arm 3's features and the swap kinds, and compute arm 3's probability from frozen numbers. Label-free.

`charge(letter: str) -> int`

`+1` for K or R, `−1` for D or E, `0` otherwise (including H).

`swap_kind(wild_type: str, new_letter: str) -> str`

The four rules in section 7.9, in that priority order.

`build_feature_table(scores: pandas.DataFrame) -> pandas.DataFrame`

Input: rows from `load.lab_tested_scores(...)`. Output: a copy with these added columns. It never looks at labels.

| Column | Formula | Type |
| --- | --- | --- |
| `blosum_score` | `baselines.blosum_score(wild_type, new_letter)` | int |
| `harm_arm1` | `0.0 - blosum_score` | float |
| `harm_arm2` | `0.0 - esm_score` | float |
| `size_change` | `VOLUME[new_letter] - VOLUME[wild_type]` (signed, Å³) | float |
| `hydrophobicity_change` | `HYDROPHOBICITY[new_letter] - HYDROPHOBICITY[wild_type]` (signed) | float |
| `to_proline` | 1 if `new_letter == "P"` else 0 | int |
| `from_cysteine` | 1 if `wild_type == "C"` else 0 | int |
| `from_glycine` | 1 if `wild_type == "G"` else 0 | int |
| `swap_kind` | `swap_kind(wild_type, new_letter)` | str |

The eight model features, in `FEATURE_ORDER`, are: `esm_score` (the raw score, not the harm score), `blosum_score` (the raw value, not the harm score), `entropy`, `size_change`, `hydrophobicity_change`, `to_proline`, `from_cysteine`, `from_glycine`. Size change and hydrophobicity change are **signed** differences, new minus normal, not absolute values.

Worked example for `C61G` (test `test_features_c61g`): `blosum_score = -3`, `harm_arm1 = 3.0`, `size_change = 60.1 − 108.5 = −48.4`, `hydrophobicity_change = −0.4 − 2.5 = −2.9`, `to_proline = 0`, `from_cysteine = 1`, `from_glycine = 0`, `swap_kind = "loss_of_cysteine"`.

`feature_matrix(table: pandas.DataFrame) -> numpy.ndarray`

Returns `table[FEATURE_ORDER].to_numpy(dtype="float64")`, shape `(n_rows, 8)`. Check there is no NaN.

`arm3_probability(x: numpy.ndarray, scaler_mean, scaler_scale, coefficients, intercept: float) -> numpy.ndarray`

Arm 3's harm score from frozen numbers, with no scikit-learn object:
1. `z = (x - scaler_mean) / scaler_scale` (column-wise).
2. `logit = z @ coefficients + intercept`.
3. Return `1.0 / (1.0 + numpy.exp(-logit))`.

Check the four arrays have length 8 and `x` has 8 columns; raise `ValueError` otherwise.

`fit_arm3(x, y, groups) -> dict`

The one function that fits arm 3. It is specified in section 8.7, next to the code that calls it.

**Must not.** Call `load.load_labels` (labels arrive as the `y` argument of `fit_arm3`; this module never fetches them). Fit anything outside `fit_arm3`.

### 8.7 `src/freeze.py`

**Purpose.** Do all tuning, on the start region only, and write `outputs/frozen.json`.

#### Private helpers worth naming

`_start_tuning_table() -> pandas.DataFrame`
1. `labels = load.load_labels("start", "tuning")` (575 rows).
2. `scores = load.lab_tested_scores("start")` (575 rows).
3. Merge on `mutant` with `validate="one_to_one"`. Check 575 rows.
4. `table = features.build_feature_table(merged)`.
5. **Drop INT rows:** keep only `findlay_class` in `{"LOF", "FUNC"}`. Check 516 rows: 131 LOF and 385 FUNC.
6. Add `is_harmful` = 1 for LOF, 0 for FUNC.
7. Return it.

`_fit_arm3(table) -> dict`

A thin wrapper: `x = features.feature_matrix(table)`; `y = table["is_harmful"].to_numpy()`; `groups = table["position"].to_numpy()`; return `features.fit_arm3(x, y, groups)`.

The fitting procedure itself is `features.fit_arm3(x: numpy.ndarray, y: numpy.ndarray, groups: numpy.ndarray) -> dict`, in `src/features.py`, so that `evaluate` can reuse it for the reverse split without importing `freeze`. It takes arrays, not regions: it knows nothing about where its rows came from. **The caller is responsible for passing rows from one region only, with INT rows already removed.** It imports scikit-learn inside the function body: `sklearn.pipeline.make_pipeline`, `sklearn.preprocessing.StandardScaler`, `sklearn.linear_model.LogisticRegression`, `sklearn.model_selection.StratifiedGroupKFold`. It returns `{"scaler_mean", "scaler_scale", "coefficients", "intercept", "cutoff", "converged", "n_iter", "out_of_fold_probabilities"}`. Its steps:

1. Check: `x` has 8 columns and no NaN; `y` holds only 0 and 1 and both are present; `x`, `y`, `groups` have the same length.
2. Splitter: `StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)`. This keeps all rows of a position in one fold (no leakage between near-copies), keeps the share of harmful rows similar across folds, and is fully determined by the seed. Do not use plain `KFold`, `StratifiedKFold` or `train_test_split`.
3. Out-of-fold probabilities: make an array `oof` of length `len(y)` (516 for the start region) filled with NaN. For each `(train_index, test_index)` from `splitter.split(x, y, groups)`:
   1. Check no position appears in both `groups[train_index]` and `groups[test_index]`.
   2. Build a **new** pipeline: `make_pipeline(StandardScaler(), LogisticRegression(random_state=42))`. All other settings are scikit-learn defaults: `C=1.0`, solver `lbfgs`, `max_iter=100`. Do not pass `penalty` (it is deprecated in scikit-learn 1.8). Do not set `class_weight`.
   3. Fit it on the training fold only. The scaler is inside the pipeline, so it is fit on the training fold only.
   4. `oof[test_index] = pipeline.predict_proba(x[test_index])[:, 1]`. Column 1 is the probability of class 1 (harmful); confirm with `pipeline.classes_`, which must be `[0, 1]`.
4. Check `oof` has no NaN left.
5. `cutoff = cutoffs.choose_catch_cutoff(oof, y, 90)`.
6. Refit: a new pipeline with the same settings, fit on all the rows passed in (516 for the start region).
7. Read out `scaler.mean_`, `scaler.scale_`, `model.coef_[0]`, `model.intercept_[0]`, `model.n_iter_[0]`.
8. Self-check: `features.arm3_probability(x, mean, scale, coef, intercept)` must equal `pipeline.predict_proba(x)[:, 1]` to within `1e-9` everywhere. If not, raise `RuntimeError("Arm 3 manual probability does not match scikit-learn.")`.
9. `converged = n_iter_ < 100`. If not converged, **do not** change `max_iter` (rule 7). Record `converged: false`, print a warning, and tell the humans at checkpoint B.
10. Return the dictionary described above. `freeze` copies its numbers into `frozen.json`.

#### `main()` steps

`python -m src.freeze [--force]`

1. `set_seeds()`.
2. **Refuse once the test has started:** if `load.test_attempted()` (any `started`, `success` or `failed` line in the run log), print `Refusing to freeze: the test run has started. Cutoffs may not change.` and exit 1. `--force` does not override this. It holds even if the run failed part-way, because a failed run may already have written end-region results.
3. **Refuse to overwrite:** if `outputs/frozen.json` exists and `--force` was not given, print `Refusing to overwrite outputs/frozen.json. Pass --force only if a human asked for a re-freeze before the test run.` and exit 1.
4. `load.check_inputs()`.
5. `table = _start_tuning_table()`.
6. Arm 1: `cutoff_1 = choose_catch_cutoff(table.harm_arm1, table.is_harmful, 90)`. Add `0.0` to the result so a negative zero becomes `0.0`. Verified expectation: cutoff 0.0, 120 of 131 caught, 217 of 385 false alarms.
7. Arm 2: `cutoff_a = choose_catch_cutoff(table.harm_arm2, table.is_harmful, 90)`.
8. Cutoff B: `cutoff_b = choose_false_alarm_cutoff(table.harm_arm2, table.is_harmful, 10)`.
9. `band_exists = has_band(cutoff_a, cutoff_b)`.
10. Arm 3: `_fit_arm3(table)`.
11. For each arm and cutoff, `count_outcomes` on the start table, for the `tuning_*` keys.
12. Hash `outputs/esm_scores.csv` and the slim CSV with `config.sha256_of_file`.
13. Get the git commit with `subprocess.run(["git", "rev-parse", "HEAD"], ...)` and the dirty flag with `git status --porcelain`.
14. Build the dictionary of section 6.5. Convert every NumPy number to a plain Python `float`, `int` or `bool` before dumping (NumPy types are not JSON-serialisable).
15. Write to `outputs/frozen.json.tmp`, close, then `os.replace` to `outputs/frozen.json`.
16. Read it back with `json.load` and check it round-trips.
17. Print a summary: each arm's cutoff, start-region caught and false alarms with percentages, the share of all 516 rows flagged, whether there is a band, the eight coefficients next to their feature names, and `converged`.

**Must not.** Call `load.load_labels("end", ...)`. Use INT rows. Try more than one model setting. Round any stored number.

### 8.8 How the arms are applied to a table (shared by `freeze` and `evaluate`)

For any table from `features.build_feature_table`:

| Arm | Harm score | Flag |
| --- | --- | --- |
| arm0 | constant `0.0` | never |
| arm1 | `harm_arm1` | `cutoffs.flag(harm_arm1, frozen["arms"]["arm1"]["cutoff"])` |
| arm2 | `harm_arm2` | `cutoffs.flag(harm_arm2, frozen["arms"]["arm2"]["cutoff"])` |
| arm3 | `features.arm3_probability(feature_matrix(table), …frozen arm3 numbers…)` | `cutoffs.flag(harm_arm3, frozen["arms"]["arm3"]["cutoff"])` |

Put this in one function in `src/evaluate.py`, `apply_frozen_arms(table: pandas.DataFrame, frozen: dict) -> pandas.DataFrame`, which returns a copy with `harm_arm3`, `flag_arm1`, `flag_arm2`, `flag_arm3` and `band_arm2` added. It is pure (no file access), so it can be unit-tested on made-up data.

### 8.9 `src/evaluate.py`

**Purpose.** The one test run. It is the only code that ever puts end-region lab labels next to scores before `test_run_done()` is true.

#### Pure functions (no file access; unit-tested on made-up data)

- `apply_frozen_arms(table, frozen) -> DataFrame` — section 8.8.
- `headline_table(pred: DataFrame, frozen: dict) -> DataFrame` — the four rows of `headline.csv` from rows whose `findlay_class` is LOF or FUNC.
- `auroc_table(pred: DataFrame) -> DataFrame` — uses `sklearn.metrics.roc_auc_score(y_true, y_score)` with `y_true = is_harmful` (1 = LOF, 0 = FUNC) and `y_score = ` the arm's harm score. Higher harm score must mean more harmful; never pass the raw ESM score. Arm 0 is written as 0.5 without calling the function.
- `bootstrap_auroc_differences(pred: DataFrame, n_resamples: int, seed: int) -> dict` — below.
- `band_count_table(pred: DataFrame) -> DataFrame`
- `by_swap_table(pred: DataFrame) -> DataFrame`
- `int_table(pred_region: DataFrame, region: str) -> DataFrame`
- `choose_sentences(headline: DataFrame, bootstrap: dict) -> list[str]` — step 10 below.
- `verify_hash(frozen: dict, scores_path: Path) -> None` — recompute `config.sha256_of_file(scores_path)` and compare with `frozen["inputs"]["esm_scores_sha256"]`. On mismatch raise `RuntimeError("HASH MISMATCH: outputs/esm_scores.csv has changed since frozen.json was written. Expected <expected>, found <found>. Refusing to run.")`. Takes the path as an argument so the unit test can point it at a temporary copy.
- `prior_runs(run_log_path: Path) -> list[dict]` — parse the run log into dictionaries; an absent file gives an empty list.
- `may_run(run_log_path: Path, bugfix_reason: str | None) -> tuple[bool, str]` — the refusal rule of step 1.
- `locked_hashes(run_log_path: Path) -> tuple[str, str] | None` — the `esm_sha256` and `frozen_sha256` on the first `status=started` line of the log, or `None` if there is none.

**The bootstrap, exactly:**

1. Keep rows with `is_harmful` 0 or 1 (drop INT).
2. `positions` = the sorted distinct `position` values among those rows. `n = len(positions)`.
3. Pre-compute, for each position, the arrays of `is_harmful`, `harm_arm1`, `harm_arm2`, `harm_arm3` for its rows.
4. `rng = numpy.random.default_rng(42)`.
5. Repeat 1,000 times:
   1. `chosen = rng.choice(positions, size=n, replace=True)`. Positions are resampled, not rows.
   2. Build the resampled row set by concatenating the rows of each chosen position, **keeping duplicates**: a position drawn three times contributes all its rows three times.
   3. If the resampled labels contain only one class, count it as skipped and go to the next repeat. Do not redraw.
   4. Compute `roc_auc_score` for arm 1, arm 2 and arm 3 on the resampled rows. Store `auroc_arm2 − auroc_arm1` and `auroc_arm3 − auroc_arm2`.
6. `low, high = numpy.percentile(differences, [2.5, 97.5])` for each of the two lists.
7. `point` = the same differences computed on the real (not resampled) rows.
8. Return the dictionary of `bootstrap.json`.

#### Command line

`python -m src.evaluate` — the test run. Only after the words "run the test" at checkpoint B.
`python -m src.evaluate --bugfix "<reason>"` — a rerun to fix a code bug. Only when a human gives the reason.
`python -m src.evaluate --selfcheck-start` — a rehearsal on the start region. Safe at any time after `frozen.json` exists.

There is no dry-run mode on the end region, and you must not add one.

#### `--selfcheck-start`

Exercises every table-building function on real data without touching the end region.
1. Load `frozen.json`, run `verify_hash`.
2. `labels = load.load_labels("start", "tuning")`; merge with `load.lab_tested_scores("start")`; build features; `apply_frozen_arms`.
3. Build the headline table, AUROC table, bootstrap (use 200 resamples here), band counts, by-swap table and INT table **for the start region**.
4. Write them under `outputs/selfcheck_start/` (git-ignored). Never write to `outputs/results/`. Never write the run log. It builds **only** the tables named in step 3. It never runs the reverse split (step 8 of the test run) or the sanity test 7 report (step 9), and it never calls `load.load_labels` with `"end"`. Do not implement it as "the test run with the region swapped": write it as its own short function that calls the pure table functions.
5. Check: arm 1's row shows 120 caught of 131 and 217 false alarms of 385; arm 2's catch rate is at least 90%; arm 3's counts equal `tuning_caught` and `tuning_false_alarms` in `frozen.json`.
6. Print `SELFCHECK PASSED (start region, in-sample; these are not results).`

These numbers are in-sample and must never appear in the write-up.

#### The test run: exact order of operations

1. **Refuse a repeat.** Call `may_run(RUN_LOG_PATH, bugfix_reason)`, which implements this rule: `runs = prior_runs(RUN_LOG_PATH)`; if any earlier line has `status=started`, `status=success` or `status=failed` (that is, `load.test_attempted()` is `True`), then end-region labels may already have been compared with scores. Unless `--bugfix "<reason>"` was given with a non-empty reason, print `Refusing to run: outputs/results/run_log.txt shows the test has already run. A rerun is allowed only to fix a code bug: python -m src.evaluate --bugfix "<what the bug was>"` and exit 1 **without** writing to the log. Lines with `status=refused` do not count as runs.
2. **Load and verify.** Load `frozen.json`. Check `schema_version == 1` and that `arms.arm3.features` equals `FEATURE_ORDER`. Run `verify_hash`. Check `sha256_of_file(SLIM_CSV_PATH)` equals `frozen["inputs"]["slim_csv_sha256"]`. Run `load.check_inputs()`. **On a rerun** (the log already has a `started` line): `sha256_of_file` of `esm_scores.csv` and of `frozen.json`, **in that order**, must equal the pair (`esm_sha256`, `frozen_sha256`) returned by `locked_hashes`; if either differs, refuse with `Refusing to run: frozen.json or esm_scores.csv is not the file the first run used.` If anything in this step fails, append a `status=refused` line to the run log with the error as the reason, print the error and exit 1. No end-region label has been read at this point.
3. **Compute every arm on the end region, using only frozen numbers.**
   1. **First append a `status=started` line to the run log** (create `outputs/results/` if needed) and close the file. Only then: `labels = load.load_labels("end", "evaluate")` — 1,262 rows. This is the moment the test begins, and it is the single place in the code where this call appears. Everything after the `started` line, up to and including step 11, runs inside `try: … except BaseException:` (not `except Exception`, so that Ctrl+C is caught too), which appends a `status=failed` line and re-raises. If the process is killed outright, the `started` line alone remains, and it still blocks a silent second run.
   2. `scores = load.lab_tested_scores("end")`; merge on `mutant`, `validate="one_to_one"`; check 1,262 rows.
   3. `table = features.build_feature_table(merged)`; add `is_harmful` (1 LOF, 0 FUNC, empty for INT).
   4. `pred = apply_frozen_arms(table, frozen)`.
   5. Check: 267 LOF, 96 INT, 899 FUNC; no NaN in any harm score.
4. **Write the headline table first.** Create `outputs/results/`. Write `headline.csv`. (`headline.md` is written in step 5, because it carries the AUROC column.) Check the arm0 row is 0 caught, 267 missed, 0 false alarms, 899 cleared, 77.1% accuracy.
5. **AUROC.** Write `auroc.csv`, then `headline.md`.
6. **Bootstrap.** Write `bootstrap.json`.
7. **Band counts, by-swap table, INT table.** Write `band_counts.csv`, `by_swap.csv`. For `int_rows.csv`: the end-region rows come from `pred`; the start-region rows come from `load.load_labels("start", "evaluate")` put through the same feature and `apply_frozen_arms` steps. Also write `predictions_end.csv`.
8. **The reverse split, last, clearly labelled secondary.** Only now, after every file above is on disk:
   1. Tuning table = end-region LOF and FUNC rows (267 + 899 = 1,166). Test table = start-region LOF and FUNC rows (131 + 385 = 516).
   2. Arm 1 and arm 2: `choose_catch_cutoff` on the end-region rows (`needed` = 241 of 267).
   3. Arm 3: call `features.fit_arm3(x, y, groups)` on the end-region rows (the same function `freeze` used on the start region). This is the only place `evaluate` fits a model.
   4. Apply the three reverse cutoffs and the reverse arm 3 model to the start-region rows. Build a table with the columns of `headline.csv` plus `direction` and `auroc`.
   5. Write `reverse_split.csv`.
   6. These reverse numbers are never written to `frozen.json` and never used by the demo, the heat map or the headline.
9. **Sanity test 7 report.** Write `sanity7.json`: the ESM scores of `M1775R` and `A1708E`; Spearman correlation (`scipy.stats.spearmanr`) between `esm_score` and `function_score` over all 1,837 lab-tested rows, over the 575 start rows and over the 1,262 end rows (INT included in all three).
10. **Sentences.** `choose_sentences` and write `result_sentences.txt`. Rules, using `d = bootstrap["esm_minus_blosum"]` and `t = bootstrap["arm3_minus_esm"]`, and arm 2's row of the headline table:
    - **S1** applies if `d.low > 0` (the whole interval is above zero).
    - **S2** applies if `d.low <= 0 <= d.high` (the interval contains zero).
    - If `d.high < 0` (BLOSUM ranked better than ESM), neither S1 nor S2 applies. Write the line `NO FIXED SENTENCE: the interval is entirely below zero, so BLOSUM62 ranked better than ESM-2. Report this plainly and ask the humans for wording.`
    - **S4** applies in addition to S1 if `d.low > 0` and `d.point < 0.05`.
    - **S3** applies if `t.low <= 0` (arm 3 did not clearly beat arm 2).
    - **S5** applies to each of arms 1, 2 and 3 whose `catch_rate_percent` differs from 90 by more than 5.0 points (below 85.0 or above 95.0). Write one S5 line per such arm, naming the arm, with X filled in.
    - **S6** is always written for arm 2, with X = arm 2's `catch_rate_percent` and Y = arm 2's `false_alarms`.
    - Percentages are written with one decimal, for example `83.5%`.
11. **Append the closing line to the run log**, `status=success`, `kind=first` or `kind=bugfix`, with the same `run=` number as the `started` line. Open the file in append mode with `newline="\n"`, write one line ending in `\n`, close it. This is the last file written. Only after this line exists does `test_run_done()` return `True`.
12. **Print** the headline table, the AUROC table, the bootstrap interval, the band counts and the contents of `result_sentences.txt`.

A rerun with `--bugfix` overwrites the result files. **Before any `--bugfix` rerun, commit `outputs/results/` exactly as it stands** (message `results as left by run N before bug-fix rerun`), even if the run failed and the files are partial, so nothing that was seen is lost. The run log keeps every attempt. The write-up must state the number of runs: the count of `status=started` lines.

#### Must not

- Read any cutoff, coefficient or setting from anywhere but `frozen.json` for the headline.
- Fit anything before step 8.
- Change `frozen.json` or `esm_scores.csv`.
- Offer a way to see end-region results without writing the run log.
- Be imported and called from the notebook or the tests in a way that performs the run. The notebook may import its **pure** functions only.

### 8.10 `src/failures.py`

**Purpose.** Stage 6: list the cases the model got most wrong.

`python -m src.failures [--full-table PATH]`

1. If `load.test_run_done()` is `False`, print `Refusing to run: the test run is not done.` and exit 1.
2. Read `outputs/results/predictions_end.csv`.
3. **Worst misses** (arm 2): rows with `findlay_class == "LOF"` and `flag_arm2 == 0`. Sort by `harm_arm2` ascending (the model was most sure these were harmless), ties broken by `function_score` ascending, then `mutant` ascending. Take the first five. If there are fewer than five, take them all.
4. **Worst false alarms** (arm 2): rows with `findlay_class == "FUNC"` and `flag_arm2 == 1`. Sort by `harm_arm2` descending, ties broken by `function_score` descending, then `mutant` ascending. Take the first five.
5. Write `outputs/results/failures.csv` with the columns of section 6.6 and an empty `guess_why`. Print the table.
6. After the table is printed, you (the agent) write one plain-language guess for each row into `guess_why` and into `failures.md`. Every guess starts with `Guess:`. Use only what the columns show (the letters, the swap kind, the BLOSUM value, the entropy, whether arms 1 and 3 agreed). Do not state biology as fact. The one exception is section 7.10: you may say that a position is one of the eight zinc-holding letters.

**The optional full study table (splicing).** The study's full table of 3,893 DNA changes may or may not be present. Its file name and column names are **unknown**.
0. This branch runs only after the `test_run_done()` check at the top of this module has passed.
1. If `--full-table PATH` is given, or if `data/` holds a file that is not one of the three known input files and is not inside `DMS_ProteinGym_substitutions/`: read only its header and first three rows, print the file name, the column names and those three rows, then print `STOP: human decision needed. Tell the humans the column names above. Do not guess which column is the RNA score.` and exit 0.
2. Implement nothing further for this file until a human names the columns to use and the cutoff for a "low" RNA score.
3. If no such file exists, print `No full study table found. Splicing is stated as a limit we could not test.` The write-up then uses the best-explained protein-level miss as its wrong case.

**Must not.** Run before the test. Present a guess as a fact.

### 8.11 `src/heatmap.py`

**Purpose.** Stage 7: draw the ESM scores as a heat map.

`python -m src.heatmap`

**What is drawn.**

- One panel per region. X axis: protein position (start panel 1–101; end panel 1631–1855). Y axis: the 20 letters in `AA_ORDER`, with `A` at the top.
- Cell colour: `esm_score` for that position and new letter. The cell where the letter equals the normal letter has no score; leave it blank (NaN).
- Colour scale: `Reds_r`, `zmin = -15`, `zmax = 0` (no `zmid`). Colour-bar title: `ESM score (darker = model finds the swap more surprising)`.
- **Lab-tested marker:** a small grey dot (size 3) on every cell whose swap was lab-tested and is not LOF.
- **LOF marker:** a black dot (size 6) on every cell whose swap the lab found to be loss of function. Give it a thin white outline so it stays visible on dark cells.
- Title: `BRCA Spellcheck: ESM-2 scores for every swap at the lab-tested positions`. Subtitle text under it: the disclaimer.

**The before/after rule (section 3, rules 11 and 15).**

1. `done = load.test_run_done()`.
2. If `done` is `False`: draw the **start panel only**. Get labels with `load.load_labels("start", "tuning")`. Do not load, draw or print anything about the end region, not even its scores. Add the figure note `End region hidden until the one test run is done.`
3. If `done` is `True`: draw both panels. Get end-region labels with `load.load_labels("end", "post_test")`.

**Building the grid.** For a region: pivot `load.load_scores()` rows of that region to a 20 × n_positions array (rows in `AA_ORDER`, columns in position order). Check each column has exactly one blank cell and it is the normal letter's row.

**HTML output.** Plotly: `plotly.subplots.make_subplots(rows=2, cols=1)` (or one row before the test), one `plotly.graph_objects.Heatmap` per panel with `coloraxis` shared, and one `plotly.graph_objects.Scatter` (mode `markers`) per marker type per panel. Hover text for a cell: `C61G · ESM score −9.12 · lab: LOF` (or `lab: not tested`). Write with `fig.write_html(HEATMAP_HTML_PATH, include_plotlyjs="cdn")`.

**PNG output.** Always drawn with matplotlib, not exported from Plotly. Plotly's image export needs an extra package and a Chrome install, which is a known failure point; matplotlib needs neither. Use `matplotlib.use("Agg")` before importing `pyplot`, `imshow` with `cmap="Reds_r"`, `vmin=-15`, `vmax=0`, `aspect="auto"`, the same markers with `scatter`, and `fig.savefig(HEATMAP_PNG_PATH, dpi=200, bbox_inches="tight")`. The two files must show the same data.

**Zinc positions.** Draw small triangles above the start panel at the positions in `ZINC_POSITIONS`, with the legend label `zinc-holding letters (PDB 1JM7)` (section 7.10).

**Must not.** Show any end-region cell before the test run. Rescale colours after looking at the data.

### 8.12 `src/export_app.py`

**Purpose.** Build `app/app_data.csv` and `app/app_config.json`, the only two data files the demo reads. This module is where the demo's before/after rule is enforced: before the test run, end-region lab results are simply not in the file, so the app cannot show them.

`python -m src.export_app`

1. `scores = load.load_scores()` (6,194 rows).
2. `frozen = json.load(FROZEN_PATH)`. Run the same hash check as `evaluate.verify_hash` (re-implement the three-line comparison here; do not import `evaluate`).
3. `done = load.test_run_done()`.
4. Start-region labels: `load.load_labels("start", "tuning")`.
5. End-region labels: if `done`, `load.load_labels("end", "post_test")`; otherwise an empty frame.
6. Left-merge the labels onto `scores` by `mutant`. `lab_class = findlay_class`; `lab_score = round(function_score, 3)`. Rows with no label keep empty cells.
7. Write `app/app_data.csv` with the six columns of section 6.7.
8. `exam_counts`: if `done`, take the 1,262 lab-tested end-region rows (scores merged with the end-region labels of step 5), compute `cutoffs.band_verdict(0.0 - esm_score, cutoff_a, cutoff_b)` for each, and count rows by verdict and `findlay_class` into the object of section 6.7 (all nine counts present, zeros included). If not `done`, `exam_counts` is `None`. Do not read `outputs/results/band_counts.csv` here; stage 8 checks the two agree.
9. Write `app/app_config.json` with `cutoff_a`, `cutoff_b`, `has_band` from `frozen["band"]`, `model_id`, `end_labels_included = done`, `protein_length = 1863`, `n_tested_positions = 326`, `exam_counts`.
10. Print the row count, how many rows have a `lab_class` (575 before the test; 1,837 after), `end_labels_included`, and `exam_counts`.

Run it in stage 8, and again after the test run if stage 8 was started early.

### 8.13 `app/app.py`

**Purpose.** The Gradio demo for a Hugging Face Space. It uses arm 2 only. It loads no model and needs no GPU.

**Hard rules.**

- It reads only `app_data.csv` and `app_config.json`, found next to itself: `Path(__file__).resolve().parent`.
- It imports only `gradio`, `pandas`, `json`, `re`, `pathlib`. Nothing from `src`. No torch, no transformers.
- It never crashes on any input. Every input gives either an answer or a clear refusal.
- **Both fixed notes appear in every response**, including refusals.
- `demo.launch()` is called with no arguments, inside `if __name__ == "__main__":`. Never pass `share=True`.

#### Input formats accepted

After removing leading and trailing whitespace, and one optional leading `p.` (any case):

| Format | Example | Reading |
| --- | --- | --- |
| One-letter, any case | `C61G`, `c61g`, `C61g` | C, 61, G |
| With prefix | `p.C61G`, `P.c61g` | C, 61, G |
| Three-letter, any case | `Cys61Gly`, `p.Cys61Gly`, `CYS61GLY` | C, 61, G |

Regular expressions, tried in this order, case-insensitive, matching the whole string:
1. `^([A-Za-z])\s*(\d{1,5})\s*([A-Za-z])$`
2. `^([A-Za-z]{3})\s*(\d{1,5})\s*([A-Za-z]{3})$`

For pattern 1, upper-case both letters; each must be in the 20 standard letters. For pattern 2, upper-case and look up in `THREE_TO_ONE`. Anything else is unreadable. Not accepted (all unreadable): DNA notation such as `c.181T>G`, stop codons (`*`, `Ter`, `X`), deletions, more than one mutation.

#### Functions

`parse_input(text: str) -> tuple[str, int, str]`

Returns `(wild_type, position, new_letter)` or raises `ValueError(message)` where `message` is one of the refusal messages below. Order of checks:
1. `text` is `None`, or empty after stripping → refusal R1.
2. No pattern matches, or a letter is not one of the 20 (for example `X61G`, `hello`) → refusal R2.
3. Both letters are the same → refusal R3.
4. Position is less than 1 or greater than 1,863 → refusal R4.

`lookup(wild_type: str, position: int, new_letter: str) -> dict`

Uses the data loaded at start-up.
1. If `position` is not one of the 326 positions in `app_data.csv` → raise `ValueError` with refusal R5.
   Then, if `position` is in the end region (1631–1855) and `end_labels_included` in the config is `false` → raise `ValueError` with refusal R8. Before the test run the app gives no verdict and no score for any end-region mutation.
2. If `wild_type` is not the normal letter at that position → raise `ValueError` with refusal R6.
3. Find the row for `(position, new_letter)`. It always exists, because all 19 swaps are scored.
4. `harm = 0.0 - esm_score`. Verdict from the band rule (copy of `band_verdict`, section 8.3) with `cutoff_a`, `cutoff_b` from the config.
5. Return a dictionary with `mutant`, `esm_score`, `verdict`, `lab_class`, `lab_score`, `region`.

`check_mutation(text: str) -> str`

The function the button calls, and the function a reviewer will call directly. Returns one Markdown string. Wrap everything in `try/except`: a `ValueError` becomes its refusal message; any other exception becomes refusal R7. The two fixed notes are appended in every case.

#### Refusal messages, word for word

`{…}` marks a value filled in.

- **R1:** `Type a mutation such as C61G: the normal letter, the position, then the new letter.`
- **R2:** `I could not read "{text}". Use the form C61G: the normal letter, the position, then the new letter. Three-letter names such as Cys61Gly also work.`
- **R3:** `"{mutant}" is not a change: both letters are {letter}.`
- **R4:** `Position {position} is outside BRCA1, which has 1,863 positions.`
- **R5:** `Position {position} is not one of the 326 lab-tested positions (1 to 101 and 1631 to 1855). This demo only answers where its answers can be checked against the lab, so it will not guess here.`
- **R6:** `BRCA1 has {actual} at position {position}, not {given}. Did you mean {actual}{position}{new}?` If `actual == new`, end the message after the first sentence.
- **R7:** `Something went wrong reading that input. Try a mutation such as C61G.`
- **R8:** `End-region answers (positions 1631 to 1855) are hidden until the one test run is done.`

In R2, `{text}` is the user's input cut to 40 characters, with any backtick, `<`, `>`, `*`, `_`, `[` and `]` removed so it cannot break the Markdown.

#### Output layout for a successful answer

A Markdown string with these fields, in this order:

1. `## {mutant}` — the mutation in one-letter form, such as `C61G`.
2. `**Model verdict:** {verdict}` — one of `Looks harmful`, `Can't tell`, `Looks harmless`.
3. `**Model score:** {esm_score:.2f}` followed by ` (more negative means the model finds this swap more surprising)`.
4. `**How the verdict is set:**` then one sentence. With a band: `Scores at or below {-cutoff_b:.2f} look harmful; scores above {-cutoff_a:.2f} look harmless; in between the model can't tell. Both lines were set on the start region only.` Without a band: `Scores at or below {-cutoff_a:.2f} look harmful; anything above looks harmless. The line was set on the start region only. There is no "can't tell" band, because the two cutoffs did not leave a gap.` (The cutoffs are stored on the harm scale, so they are negated to show them on the ESM-score scale.)
5. `**Lab result (answer key):**` then:
   - `lab_class` is `LOF` → `Loss of function (harmful). Lab function score {lab_score:.3f}.`
   - `FUNC` → `Functional (harmless). Lab function score {lab_score:.3f}.`
   - `INT` → `Intermediate (in between). Lab function score {lab_score:.3f}.`
   - empty → `not lab-tested`
6. `**On the held-out exam:**` then one sentence, built from `exam_counts` in the config for this answer's verdict (`h` = its `LOF` count, `f` = its `FUNC` count):
   - `exam_counts` is `null` → `Not available until the one test run is done.`
   - `h + f` is 0 → `On the end region, the model gave this verdict to no mutation the lab could classify.`
   - otherwise → `On the end region, which the cutoffs never saw, the model gave this verdict to {h + f} mutations the lab could classify: the lab found {h} harmful and {f} harmless.`
   This line is what makes the verdict usable: it tells the user how far to trust it. It appears in every successful answer, for start-region and end-region mutations alike.
7. A horizontal rule, then the two fixed notes, each on its own line:
   - `At lab-tested positions a counselor would use the lab result. This demo shows how the model does where we can check it.`
   - `This project is a prototype and is not medical advice.`

A refusal is: the refusal message, a horizontal rule, then the same two fixed notes.

#### Gradio components

Use `gradio.Blocks`:
- `gr.Markdown` — title `# BRCA Spellcheck` and two sentences: what the tool is, and "BRCA1 is our exam, not the use case: it is one of the few genes with a lab answer key, so we can check the model here."
- `gr.Textbox(label="Mutation", placeholder="C61G")`
- `gr.Button("Check")`
- `gr.Markdown` — the output.
- `gr.Examples` with `["C61G", "D2E", "C61A", "p.Cys64Gly"]`. If `end_labels_included` is true in the config, add `M1775R`, so that one example comes from the held-out region.
- `gr.Markdown` — the two fixed notes, shown permanently at the bottom as well.

Wire both the button's `click` and the textbox's `submit` to `check_mutation`.

#### `app/README.md`

Starts with this header (the Space reads it), then one paragraph of description and the disclaimer:

```
---
title: BRCA Spellcheck
emoji: 🧬
colorFrom: red
colorTo: blue
sdk: gradio
sdk_version: 6.32.0
app_file: app.py
pinned: false
---
```

### 8.14 `tools/build_notebook.py`

**Purpose.** Generate `notebook.ipynb` with `nbformat` so the notebook is never hand-edited JSON. Details in section 10, stage 9.

---

## 9. Tests

Two files. Unit tests use made-up data and temporary folders only; they never read `outputs/` or end-region labels. Sanity tests check the real data.

Run all: `.\.venv\Scripts\python.exe -m pytest tests/ -q`

A test that cannot run yet must **skip with a reason**, using `pytest.skip("<reason>")` or `pytest.mark.skipif`. It must never pass silently. The allowed skip conditions are listed below; no others.

### 9.1 Sanity tests (`tests/test_sanity.py`)

| # | Test name | Added in stage | Exact assertion | Skip condition |
| --- | --- | --- | --- | --- |
| 1 | `test_sanity_1_counts` | 1 | `load.check_inputs()` raises nothing. `load.class_counts()` equals start 131 / 59 / 385 and end 267 / 96 / 899. `load_mutations()` has 1,837 rows, unique `mutant`, 326 distinct positions. | never |
| 2 | `test_sanity_2_wild_type` | 1 | For every row of `load_mutations()`, `wild_type == load_wild_type()[position - 1]`. The sequence has 1,863 letters and starts `MDLSALRVEEVQNVINAMQK`. | never |
| 3a | `test_sanity_3a_letter_test_cached` | 3 | In `load.load_scores()`, `token_at_index == wild_type` in all 6,194 rows, and `window`, `token_index` match `esm_score.window_of` and `esm_score.token_index_of`. | skip if `outputs/esm_scores.csv` does not exist: "scores not computed yet (stage 3)" |
| 3b | `test_sanity_3b_letter_test_live` | 3 | `esm_score.run_letter_test(tokenizer, seq, scored_positions())` raises nothing and returns 326 rows, all `ok`. | skip if `transformers` cannot be imported or the tokenizer cannot be loaded (no network and no cache): give the exception text as the reason |
| 4 | `test_sanity_4_known_harmful_start` | 3 | Start region only. ESM scores of `C61G` and `C64G` are both `< 0` (the only assertion). The test also **prints**, without asserting, whether each is `<=` the 25th percentile (`numpy.percentile(x, 25)`) of the ESM scores of the 575 lab-tested start-region rows. That part is a report about the model, not a check on the index maths, so it never fails the build. | skip if the score cache does not exist |
| 5 | `test_sanity_5_start_correlation` | 3 | Start region only. `scipy.stats.spearmanr(esm_score, function_score)` over the 575 lab-tested start-region rows (INT included) has correlation `> 0.15`. Labels come from `load.load_labels("start", "tuning")`. | skip if the score cache does not exist |
| 6 | `test_sanity_6_arm0` | 1 | Using `load.class_counts()` only: arm 0 on the end region catches 0 of 267 and is right on 899 of 1,166; `round(100 * 899 / 1166, 1) == 77.1`. | never |
| 7 | `test_sanity_7_after_test` | 5 | `outputs/results/sanity7.json` exists; `M1775R_esm_score < 0`; `A1708E_esm_score < 0`; `spearman_all_1837 > 0`. | skip if `load.test_run_done()` is `False`: "runs only after the one test run" |

Notes:
- Tests 4 and 5 are sanity checks, not tuning: they choose nothing. They may use INT rows.
- Test 5's threshold of 0.15 is a tripwire for index bugs. A correct pipeline is expected to land well above it (estimate: 0.3 to 0.5). If the value is between 0 and 0.15, do not lower the threshold: stop and report.
- Test 7 never asserts on a specific Spearman value. A human compares it with the number ProteinGym publishes for ESM-2 650M on this table. Never guess that number.
- No sanity test other than 7 may load end-region labels. Tests 1 and 6 use counts only.

### 9.2 Unit tests (`tests/test_units.py`)

| Test name | Added in stage | What it checks |
| --- | --- | --- |
| `test_parse_mutant_good` | 1 | `C61G` → `("C", 61, "G")`; `S1631A` → `("S", 1631, "A")`; `M1I` → `("M", 1, "I")`. |
| `test_parse_mutant_bad` | 1 | Each of `""`, `"c61g"`, `"C61"`, `"61G"`, `"C61GG"`, `"X61G"`, `"C0G"`, `"C61G "`, `"C61C"`, `"p.C61G"` raises `ValueError`. `None` raises `TypeError`. |
| `test_region_of` | 1 | 1 and 101 → `start`; 1631 and 1855 → `end`; 0, 102, 1630, 1856 raise `ValueError`. |
| `test_class_of_boundaries` | 1 | −1.329 → `LOF`; −1.328 → `INT`; −1.0 → `INT`; −0.748 → `INT`; −0.747 → `FUNC`. |
| `test_sha256_of_file` | 1 | A temporary file holding the bytes `b"abc"` hashes to `ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad`. |
| `test_test_attempted` | 1 | Temporary log files: missing file → `False`; only `status=refused` lines → `False`; a `status=started` line → `True`; a `status=failed` line → `True`; a `status=success` line → `True`. |
| `test_test_run_done` | 1 | Temporary log files: missing file → `False`; a `status=refused` line → `False`; a `status=started` line alone → `False`; a `status=failed` line → `False`; a line containing `status=success` → `True`; a line whose reason text contains the word `success` but whose status is `failed` → `False`. |
| `test_load_labels_end_locked` | 1 | `load_labels("end", "tuning")` raises `PermissionError`. `load_labels("end", "post_test")` raises `PermissionError` when `test_run_done` is patched to return `False` (use `monkeypatch`; do not touch the real log). `load_labels("end", "evaluate")` raises `PermissionError` when `test_attempted` is patched to return `False`. `load_labels("start", "tuning")` returns 575 rows. |
| `test_catch_cutoff_worked_example` | 2 | Section 8.3 worked example 1: targets 90, 60, 30 give 0.5, 0.5, 0.9. |
| `test_catch_cutoff_matches_brute_force` | 2 | For 200 random small arrays (seeded, with scores rounded to one decimal so ties are common), the function's answer equals a brute-force search over every observed score for the highest one that reaches the target. |
| `test_catch_cutoff_needed_counts` | 2 | With 131 harmful rows of distinct scores and target 90, exactly 118 are flagged. With 267, exactly 241. |
| `test_false_alarm_cutoff_worked_example` | 2 | Section 8.3 worked example 2 gives 6.5. |
| `test_false_alarm_cutoff_none` | 2 | The second case in worked example 2 gives `None`. |
| `test_band_verdict` | 2 | The worked example in section 8.3, with and without a band. |
| `test_flag_is_inclusive` | 2 | `flag([1.0, 2.0, 3.0], 2.0)` is `[False, True, True]`. |
| `test_count_outcomes_arm0` | 2 | With `cutoff=None`: nothing flagged, `harmful_missed == n_harmful`. |
| `test_blosum_loads` | 2 | `load_blosum62()` has 400 entries; the four spot values of section 8.4; symmetric. |
| `test_harm_arm1_no_negative_zero` | 2 | For a pair whose BLOSUM62 value is 0 (for example `("A", "T")`), `str(harm_arm1("A", "T")) == "0.0"`. |
| `test_token_index_worked_examples` | 3 | `window_of` and `token_index_of` for 61 → `A`, 61; 1775 → `B`, 412; 1 → `A`, 1; 101 → `A`, 101; 1631 → `B`, 268; 1855 → `B`, 492. 102 and 1630 raise `ValueError`. |
| `test_window_sequence` | 3 | With the real wild type: window A has 500 letters, starts `MDLSA`, ends `RPLTN`; window B has 500 letters, starts `NLGEA`, ends `PHSHY`; `window_sequence(seq, "A")[60] == "C"`; `window_sequence(seq, "B")[411] == "M"`. |
| `test_swap_kind_each_kind` | 4 | `("L", "P")` → `to_proline`; `("C", "G")` → `loss_of_cysteine`; `("C", "R")` → `loss_of_cysteine`; `("C", "P")` → `to_proline`; `("S", "D")` → `charge_change`; `("L", "K")` → `charge_change`; `("H", "R")` → `charge_change`; `("K", "R")` → `other`; `("L", "S")` → `other`. |
| `test_features_c61g` | 4 | The worked example in section 8.6. |
| `test_arm3_probability_matches_sklearn` | 4 | On random made-up data, fit a `StandardScaler` + `LogisticRegression` pipeline and check `arm3_probability` with its fitted numbers equals `predict_proba[:, 1]` to `1e-9`. |
| `test_fit_arm3_groups_do_not_leak` | 4 | On made-up data with 40 groups, `features.fit_arm3` returns out-of-fold probabilities with no NaN, and a patched splitter check confirms no group is in both train and test of any fold. |
| `test_apply_frozen_arms_synthetic` | 4 | With a hand-written `frozen` dictionary and six made-up rows, flags and bands equal hand-computed values. |
| `test_headline_table_synthetic` | 4 | On made-up predictions, counts in each column equal hand-computed values and `harmful_caught + harmful_missed == n_harmful`. |
| `test_auroc_direction` | 4 | Made-up rows where harmful rows have higher harm scores give AUROC 1.0; reversed gives 0.0. |
| `test_bootstrap_resamples_positions` | 4 | On made-up data with 10 positions: the result is identical for two calls with the same seed; with a patched `roc_auc_score` that records its input sizes, resample sizes vary (they would be constant if rows were resampled); a data set where some resamples contain one class reports a non-zero skipped count. |
| `test_choose_sentences_rules` | 4 | Hand-written bootstrap dictionaries: interval (0.02, 0.10) with point 0.06 → S1 and not S4; (0.01, 0.04) with point 0.02 → S1 and S4; (−0.02, 0.05) → S2; (−0.08, −0.01) → the `NO FIXED SENTENCE` line; arm 2 catch rate 83.0 → an S5 line; 91.0 → no S5 line; S6 always present. |
| `test_verify_hash_on_temp_copy` | 4 | Copy a small made-up CSV to `tmp_path`, put its hash in a made-up `frozen` dictionary: `verify_hash` passes. Append one byte to the copy: `verify_hash` raises `RuntimeError` containing `HASH MISMATCH`. Never touches the real files. |
| `test_may_run_refusal` | 4 | Temporary logs: no file → may run; only `refused` lines → may run; a lone `started` line and no bugfix reason → refused; a `success` line and no bugfix reason → refused; a `failed` line and no bugfix reason → refused; a `success` line with reason `"fixed by-swap bug"` → may run; an empty-string reason → refused. |
| `test_freeze_refuses_overwrite` | 4 | With paths patched to `tmp_path` and an existing `frozen.json` there, `freeze.main([])` exits non-zero; and with a temporary run log holding a `status=success` line, or only a `status=started` line, or a `status=failed` line, it returns non-zero even with `--force`. |
| `test_locked_hashes_block_changed_inputs` | 4 | With a temporary log whose `started` line records two hashes, `locked_hashes` returns them in the order (`esm_sha256`, `frozen_sha256`); `None` for a log with no `started` line. With paths patched to `tmp_path`, a `started` line in the log and a `frozen.json` whose hash differs from the logged one, `evaluate.main(["--bugfix", "x"])` returns non-zero, appends a `status=refused` line, and never calls `load_labels` (patch it to raise if called). |
| `test_only_refuses_end_region` | 3 | `esm_score.main(["--only", "S1631A"])` returns 1 without loading the model when `test_run_done` is patched to `False` (patch `load_model` to raise if it is called). |
| `test_app_parser_good` | 8 | `C61G`, `c61g`, `C61g`, ` C61G `, `p.C61G`, `P.c61g`, `Cys61Gly`, `p.Cys61Gly`, `CYS61GLY` all parse to `("C", 61, "G")`. |
| `test_app_refusals` | 8 | `""` and `None` → R1; `"hello"`, `"X61G"`, `"c.181T>G"`, `"C61*"`, `"Cys61Ter"` → R2; `"C61C"` → R3; `"C0G"` and `"C5000G"` → R4; `"C999G"` → R5; `"A61G"` → R6 with the text `Did you mean C61G?`. Compare against the exact messages. |
| `test_app_notes_always_present` | 8 | For each of `C61G`, `D2E`, `C61A`, `""`, `"hello"`, `"X61G"`, `"C999G"`: the returned text contains both fixed notes word for word, and the call does not raise. |
| `test_app_lab_display` | 8 | `C61G` shows `Loss of function (harmful)`; `C61A` shows `not lab-tested`. |
| `test_app_exam_line` | 8 | With `exam_counts` patched to a made-up object, the answer for `C61G` contains the exam sentence of section 8.13 with the right three numbers for its verdict; with a verdict whose `LOF` and `FUNC` counts are both 0, the "no mutation" sentence; with `exam_counts` patched to `None`, `Not available until the one test run is done.` Refusals contain no exam line. |
| `test_app_end_region_locked` | 8 | With `end_labels_included` patched to `false`, `S1631A` returns refusal R8 word for word, with no verdict and no score in the text, and still contains both fixed notes. |
| `test_app_imports_no_model` | 8 | After importing `app/app.py`, `"torch"` and `"transformers"` are not in `sys.modules` (run this check in a fresh subprocess). |

The app tests import `app/app.py` by path with `importlib.util.spec_from_file_location`, because `app` is not a package. They skip with the reason "app data not exported yet (stage 8)" if `app/app_data.csv` does not exist.

---

## 10. Stages

Twelve stages, in order, plus the short stage 5b (version zero). Each block lists what to build and how to prove it works. "Time" is a target. If a stage runs more than 50% over its target, say so plainly at the end of the stage and propose a cut from section 14.

**Clock targets (US Central).**

| When | Target |
| --- | --- |
| Overnight Saturday to Sunday | Stages 1–4, ending at checkpoint B. If it gets too late, stop at checkpoint A and do stage 4 first thing Sunday. |
| Sunday 9:00–9:30 AM | Stage 5, the test run |
| As soon as stage 5b is committed | **Version zero submitted** by the humans: the code link and the write-up link both point at the GitHub repository |
| Sunday 9:30 AM–noon | Stages 6, 7, 8 |
| **Sunday noon** | **Feature freeze.** After this, only the notebook, demo polish, write-up and video. |
| Sunday noon–4 PM | Stages 9–12 |
| **Sunday 8:00 PM** | **Version one submitted** |
| Sunday 10:59 PM | Hard deadline |

In every stage, "python" means `.\.venv\Scripts\python.exe`, run from the repository root.

Every stage ends with:
1. `python -m pytest tests/ -q` — everything passes or skips for an allowed reason.
2. One new line per decision in `writeup/decision_log.md`, and any new source added to `writeup/credits.md`.
3. `git add -A` and `git commit -m "<the stage's message>"`.
4. The end-of-stage report (section 0, point 3).
5. The `test-runner` review, if you are Claude Code (section 11). In any other tool, your own pytest run in step 1 stands in for it: do not wait.

---

### Stage 1: Setup

- **Goal:** a working environment, the repository skeleton, the data loader, and proof the input files match the known facts.
- **Machine:** any.
- **Time:** 30 minutes (most of it is the PyTorch download).
- **Already in place (do not recreate):** the git repository, `.gitignore`, `.gitattributes`, the three files in `data/`, `AGENTS.md`, `BUILD_SPEC.md`, `START_PROMPT.md`.
- **Files created:** `.venv/`, `requirements.txt`, `pytest.ini`, `README.md` (stub: project name, one sentence, the disclaimer, "work in progress"), `src/__init__.py`, `src/config.py`, `src/load.py`, `tests/test_units.py`, `tests/test_sanity.py`, `writeup/decision_log.md`, `writeup/credits.md`.

**Steps**
1. `Get-Date`. `git status` and `git log --oneline -3` to confirm the repository exists and is clean.
2. Create the environment and install PyTorch (sections 5.1, 5.2). Run the verify command and paste its four lines.
3. Write `requirements.txt` exactly as in section 5.3 and install it.
4. Write `pytest.ini` exactly as in section 5.6.
5. Write `src/config.py` (sections 7 and 8.1).
6. Write `src/load.py` (section 8.2).
7. Write the stage-1 unit tests and sanity tests 1, 2 and 6 (section 9).
8. Start `writeup/decision_log.md` with one line each for: the label rule; the region split; seed 42; the Kyte-Doolittle source; the Zamyatnin source; histidine counted as neutral; signed size and hydrophobicity changes; `StratifiedGroupKFold` with shuffle and seed 42; regions derived from position because the table has no region column of its own.
9. Start `writeup/credits.md` with the full list in section 12.6.

**Verify**
```powershell
.\.venv\Scripts\python.exe -c "import torch, transformers, sklearn, Bio, pandas, numpy, scipy, plotly, matplotlib, gradio, nbformat; print('imports ok')"
.\.venv\Scripts\python.exe -m src.load --check
.\.venv\Scripts\python.exe -m src.load --rebuild-check
.\.venv\Scripts\python.exe -m pytest tests/ -q
git status --short
```

**Expected output**
- `imports ok`
- `--check` prints one line per passed check and this table (exact, verified):
  ```
  region  LOF  INT  FUNC  total
  start   131   59   385    575
  end     267   96   899   1262
  ```
- `--rebuild-check` ends with `REBUILD MATCHES: slim CSV and FASTA agree with the raw ProteinGym file.`
- pytest: all stage-1 tests pass; none fail. Sanity tests 3a, 3b, 4, 5, 7 do not exist yet.
- `git status --short` shows no change to any file in `data/`.

**Acceptance criteria**
- [ ] `torch.cuda.is_available()` printed, and the result reported (True or False).
- [ ] All packages import.
- [ ] The three data files' SHA-256 values equal the constants in section 7.1.
- [ ] Class counts match the table exactly.
- [ ] Sanity tests 1, 2, 6 and all stage-1 unit tests pass.
- [ ] No file in `data/` was changed.

**Commit message:** `stage 1: environment, config, data loader, sanity tests 1 2 6`

---

### Stage 2: Baselines

- **Goal:** arm 0, arm 1, and the shared cutoff logic, with arm 1's cutoff chosen on the start region.
- **Machine:** any.
- **Time:** 20 minutes.
- **Files created:** `src/cutoffs.py`, `src/baselines.py`. **Changed:** `tests/test_units.py`.

**Steps**
1. Write `src/cutoffs.py` (section 8.3) and its unit tests. Run them before writing anything else.
2. Write `src/baselines.py` (section 8.4) and its unit tests.
3. Run `python -m src.baselines`.

**Verify**
```powershell
.\.venv\Scripts\python.exe -m pytest tests/ -q
.\.venv\Scripts\python.exe -m src.baselines
```

**Expected output** of `python -m src.baselines` (exact, verified against the real start-region data):
- Arm 1 cutoff on the harm scale: `0.0` (that is, flag every swap whose BLOSUM62 value is 0 or lower).
- Start region: caught **120 of 131** (91.6%). False alarms **217 of 385** (56.4%). Flagged 337 of 516 rows (65.3%).
- Arm 0, end region, from counts alone: caught 0 of 267; right on 899 of 1,166 (77.1%).

For reference, the start-region counts by arm 1 harm score (verified). Harmful: −2 → 4, −1 → 7, 0 → 8, 1 → 34, 2 → 40, 3 → 38. Harmless: −3 → 16, −2 → 56, −1 → 96, 0 → 59, 1 → 68, 2 → 59, 3 → 31. The 118th highest harmful score is 0, which is why the cutoff is 0 and why eight tied rows push the catch from 112 to 120.

**Acceptance criteria**
- [ ] All unit tests for `cutoffs` pass, including the brute-force comparison.
- [ ] Arm 1's cutoff, caught count and false-alarm count equal the verified values above.
- [ ] The cutoff prints as `0.0`, not `-0.0`.
- [ ] No end-region label was loaded (only `class_counts()` was used for arm 0).

**Commit message:** `stage 2: cutoff rule, arm 0, arm 1 (BLOSUM62)`

---

### Stage 3: ESM-2 scoring

- **Goal:** `outputs/esm_scores.csv` with 6,194 rows, and proof the index maths is right.
- **Machine:** **needs an NVIDIA GPU** for step 4 onwards. See section 5.5 if this machine has none.
- **Time:** 45 minutes (the model download is about 2.6 GB).
- **Files created:** `src/esm_score.py`, `outputs/esm_scores.csv`. **Changed:** both test files.

**Steps**
1. Write `src/esm_score.py` (section 8.5) and the stage-3 unit tests.
2. Run the letter test alone. It downloads only the tokenizer (a few kilobytes).
3. If this machine has no GPU, follow section 5.5 now.
4. Score C61G alone and read the output.
5. Score everything and write the cache.
6. Add sanity tests 3a, 3b, 4, 5 and run them.

**Verify**
```powershell
.\.venv\Scripts\python.exe -m src.esm_score --letter-test-only
.\.venv\Scripts\python.exe -m src.esm_score --only C61G
.\.venv\Scripts\python.exe -m src.esm_score
.\.venv\Scripts\python.exe -m src.esm_score --validate-cache
.\.venv\Scripts\python.exe -m pytest tests/ -q
```

**Expected output**
- `LETTER TEST PASSED: 326 of 326 positions, both windows 502 tokens.` (exact, verified)
- `C61G: position 61, window A, token index 61, token C, esm_score = <negative number>, entropy = <number between 0 and 3>`. The index part is exact and verified. The score is **unknown until run**; it must be negative. Estimate: clearly negative, somewhere around −5 to −15, because position 61 is a cysteine the protein cannot do without.
- Full run: progress lines, then `Wrote outputs/esm_scores.csv: 6194 rows` and a 64-character SHA-256.
- The C61G score in the cache agrees with the `--only` score to within 0.001.
- Sanity test 5's Spearman value: **unknown until run**. It must be above 0.15. Estimate: 0.3 to 0.5.

**Acceptance criteria**
- [ ] Letter test: 326 of 326.
- [ ] Cache has 6,194 rows and passes `validate_scores`.
- [ ] `C61G` and `C64G` scores are below zero (sanity test 4). Report whether each is in the lowest quarter of start-region scores; this does not block the stage.
- [ ] Start-region Spearman is above 0.15 (sanity test 5). Report the value.
- [ ] Sanity tests 3a and 3b pass.
- [ ] The run printed no end-region label and no score-versus-label comparison for the end region.
- [ ] `precision` is `float32` and `model_id` is `facebook/esm2_t33_650M_UR50D` (or a human approved otherwise, and the decision log says why).

**If sanity test 3, 4 or 5 fails:** stop. Do not try index offsets until something passes. Commit what you have on a clearly named work-in-progress commit only if a human asks, and report: "Sanity test N failed. The index-debugger review is needed (section 11)." Section 15 lists what that review checks.

**Commit message:** `stage 3: ESM-2 scores for 326 positions, sanity tests 3 4 5`

#### Checkpoint A (automatic)

Put these three things in your end-of-stage report, so the humans can read them later:
1. The exact `C61G` line printed by `--only C61G`.
2. The letter-test line (`326 of 326`).
3. The start-region Spearman value from sanity test 5, and whether `C61G` and `C64G` are in the lowest quarter.

If the letter test and sanity tests 3a, 3b, 4 and 5 all pass (3b may skip only for its allowed reason), **continue to stage 4 without waiting**. If any of them fails, stop completely and report, as described above.

---

### Stage 4: Cutoffs, arm 3, the frozen file

- **Goal:** all tuning done on the start region; `outputs/frozen.json` written; `src/evaluate.py` written and rehearsed on the start region, **not run on the end region**.
- **Machine:** any.
- **Time:** 60 minutes.
- **Files created:** `src/features.py`, `src/freeze.py`, `src/evaluate.py`, `outputs/frozen.json`. **Changed:** `tests/test_units.py`.

**Steps**
1. Write `src/features.py` (section 8.6) and its unit tests.
2. Write `src/freeze.py` (section 8.7) and its unit tests.
3. Write `src/evaluate.py` (sections 8.8, 8.9) and its unit tests. All of them use made-up data.
4. Commit the code before freezing, so that `code_commit` in the frozen file points at the code that produced it: `git add -A; git commit -m "stage 4a: features, freeze and evaluate code (not yet frozen)"`. Then `git status --short` must print nothing.
5. Run `python -m src.freeze`.
6. Run `python -m src.evaluate --selfcheck-start`.
7. Do **not** run `python -m src.evaluate`.

**Verify**
```powershell
.\.venv\Scripts\python.exe -m pytest tests/ -q
.\.venv\Scripts\python.exe -m src.freeze
.\.venv\Scripts\python.exe -m src.freeze
.\.venv\Scripts\python.exe -m src.evaluate --selfcheck-start
.\.venv\Scripts\python.exe -m src.esm_score --recompute
Get-Content outputs\frozen.json
Test-Path outputs\results\run_log.txt
```

**Expected output**
- First `freeze`: the summary. Arm 1: cutoff 0.0, 120 of 131, 217 of 385 (exact, verified). Arm 2 and arm 3 cutoffs and counts: **unknown until run**. Arm 2 must catch at least 118 of 131 in-sample. Arm 3 must catch at least 118 of 131 on the out-of-fold probabilities.
- Second `freeze`: `Refusing to overwrite outputs/frozen.json. …` and a non-zero exit code. This proves the guard works.
- `--selfcheck-start`: `SELFCHECK PASSED (start region, in-sample; these are not results).`
- `--recompute`: `Refusing to recompute: the current scores are locked (outputs/frozen.json exists or the test has started). Ask a human.` and a non-zero exit code.
- `Test-Path outputs\results\run_log.txt` prints `False`. The test has not run.

**Acceptance criteria**
- [ ] All unit tests pass, including the hash check on a temporary copy and the run-log refusal.
- [ ] `frozen.json` has every key in section 6.5; `tuning_counts` is 131 / 385 / 59; the feature list equals `FEATURE_ORDER`; four lists of exactly eight numbers.
- [ ] `band.cutoff_a` equals `arms.arm2.cutoff`.
- [ ] `inputs.esm_scores_sha256` equals the hash printed by `python -m src.esm_score --validate-cache`.
- [ ] The second freeze was refused, and the recompute was refused.
- [ ] The self-check passed.
- [ ] `outputs/results/run_log.txt` does not exist.
- [ ] Searching `src/` for `"end", "evaluate"` finds exactly one line, in `evaluate.py`, inside the test-run path and after the `started` line is written. The purpose string `"evaluate"` appears nowhere in `src/` outside `evaluate.py` and the permission check in `load.py`. No call to `load_labels(` passes the region as a variable.

**Commit message:** `stage 4: frozen.json written; evaluate.py written and NOT run`

#### 🛑 Checkpoint B

Stop completely. Show the humans:
1. The full contents of `outputs/frozen.json`.
2. One paragraph, in plain words, covering every choice in it: the label rule; that tuning used only the 516 start-region LOF and FUNC rows; each arm's cutoff and how many start-region mutations it flags; where arm 3's cutoff came from (out-of-fold, five folds grouped by position); whether arm 3 converged; whether there is a "Can't tell" band and where A and B sit; the model ID and precision; the hash.
3. The share of start-region rows arm 1 flags (65.3%), as rule 14 requires.
4. The `leak-auditor` report (section 11). If you are Claude Code, run that review now and paste its full report. In any other tool, say: "Stage 4 is committed. Please run leak-auditor in Claude Code."

Then wait. Two things must both happen before you continue:
- `leak-auditor` reported **CLEAN**. If it reports a leak, fix exactly what it names, re-run the stage-4 verify list, commit, and ask for it to be run again. If the fix changes anything `frozen.json` depends on, re-freeze with `python -m src.freeze --force` and say so.
- A human sends the exact words **"run the test"**.

---

### Stage 5: The one test run

- **Goal:** run `src/evaluate.py` once and save every result.
- **Machine:** any.
- **Time:** 30 minutes.
- **Files created:** everything in `outputs/results/` except `failures.*`. **Changed:** `tests/test_sanity.py` (add test 7).

**Steps**
1. `Get-Date`. Confirm you have the words "run the test" in the conversation.
2. `git status --short` must print nothing. If it prints anything, stop and ask.
3. `python -m src.evaluate` — once. This is the only time in the whole project that this bare command is typed.
4. If the log shows `status=started` without a `status=success` line after it (an error, an interruption or a crash): stop. Do not rerun. Commit `outputs/results/` exactly as it stands, report the error, and wait for a human to give a `--bugfix` reason.
5. Add sanity test 7 and run the tests.
6. Commit immediately, before looking for improvements. There are no improvements to make: nothing may be tuned after this point.

**Verify** (the run itself is step 3 above; it is **not** repeated here)
```powershell
Get-Content outputs\results\run_log.txt
Get-Content outputs\results\headline.md
Get-Content outputs\results\result_sentences.txt
Get-Content outputs\results\sanity7.json
.\.venv\Scripts\python.exe -m pytest tests/ -q
```

**Expected output**
- The arm0 row of the headline: 0 caught, 267 missed, 0 false alarms, 77.1% accuracy (exact, known in advance).
- Every other result: **unknown until run**. Do not predict them, and do not describe them as good or bad beyond the fixed sentences.
- `run_log.txt` has one `status=started` line and then one `status=success` line, both `run=1` and `kind=first`. Any other line is `status=refused` (a failed check before the run began).
- Do not run `python -m src.evaluate` a second time to "prove" the refusal. The refusal is proved by the unit test `test_may_run_refusal`.
- Sanity test 7 passes (both scores below zero, overall Spearman above zero). Report the three Spearman values.

**Acceptance criteria**
- [ ] The run log has one `status=started` line and one `status=success` line; any other line is `status=refused`.
- [ ] Every file of section 6.6 (except `failures.*`) exists and has the stated columns and row counts.
- [ ] In `headline.csv`, for every arm, `harmful_caught + harmful_missed == 267` and `false_alarms + harmless_cleared == 899`.
- [ ] `reverse_split.csv` has `SECONDARY` in its `direction` column.
- [ ] `frozen.json` and `esm_scores.csv` are byte-for-byte unchanged (`git status --short` shows them unmodified).
- [ ] Sanity test 7 passes.

**Commit message:** `stage 5: the one test run on the end region (run 1)`

#### Stage 5b: version zero (do this before you stop at checkpoint C)

- **Goal:** a first submittable version, so the team can submit hours before the deadline and resubmit later. A submission is complete once it has working links to the code and the write-up, and the organisers keep only the latest one.
- **Time:** 20 minutes. **File changed:** `README.md` only.

Replace the stub `README.md` with a short write-up that uses the headings of section 12.1, in order, with one to three sentences under each:

1. The title, `**Track:** Health`, `**Team:** [names — the humans fill this in]`, the disclaimer line, and the `**In one sentence:**` line.
2. **The problem** and **What we built**, from section 12.1, including the exam-not-use-case sentence.
3. **What we added:** `Our own idea.`, the two prior-work sentences, and the four things we added.
4. **How it works:** the data, the model, the two baselines and the split, in three sentences.
5. **Results:** the table from `outputs/results/headline.md`, pasted unchanged; the lines of `outputs/results/result_sentences.txt`, word for word; the number of test runs. Then this sentence: `Still to come in the next version: the failure analysis, the heat map, the demo and the notebook.`
6. **What we'd do next:** the two steps of section 12.1.
7. **Sources:** the contents of `writeup/credits.md`.

Rules: every number comes from a file in `outputs/results/`; nothing from section 3.5; no names except the placeholder on the `Team:` line; no adjectives about the results. Stage 10 replaces this file with the full write-up.

**Verify**
```powershell
Select-String -Path README.md -Pattern "This project is a prototype and is not medical advice." -SimpleMatch
git status --short
```

**Acceptance criteria**
- [ ] The disclaimer is present word for word.
- [ ] The table in `README.md` is identical to `outputs/results/headline.md`.
- [ ] No file other than `README.md` changed.

**Commit message:** `stage 5b: version zero README with the headline results`

The `claims-checker` review is due on `README.md` (section 11). If you are Claude Code, run it now and fix what it names before checkpoint C. In any other tool, ask for it as part of checkpoint C.

#### 🛑 Checkpoint C

Stop completely. Show the humans:
1. The headline table (`headline.md`).
2. The contents of `result_sentences.txt`, and which fixed sentence is the main one (S1 or S2, or the `NO FIXED SENTENCE` line).
3. The bootstrap interval and the three Spearman values. A human compares the overall Spearman value with the number ProteinGym publishes for ESM-2 650M on this table.
4. That version zero is committed (stage 5b), with this message: "Version zero is ready. To submit it now: fill in the Team line of README.md, push the repository to GitHub, and send the form with the repository link as both the code link and the write-up link. You can resubmit later; the organisers keep the latest version."

From now on no cutoff, feature, model setting or result file from stage 5 may be touched. Wait for the exact words **"go C"**.

---

### Stage 6: Failures

- **Goal:** the five worst misses and five worst false alarms, each with a labelled guess.
- **Machine:** any. **Time:** 45 minutes.
- **Files created:** `src/failures.py`, `outputs/results/failures.csv`, `outputs/results/failures.md`.

**Steps**
1. Write `src/failures.py` (section 8.10).
2. Run it. Read the table.
3. Write the ten guesses, each starting with `Guess:`.
4. If it printed `STOP: human decision needed`, relay the column names to the humans and continue with the rest of the stage while you wait. Splicing analysis is optional.

**Verify**
```powershell
.\.venv\Scripts\python.exe -m src.failures
Get-Content outputs\results\failures.md
.\.venv\Scripts\python.exe -m pytest tests/ -q
```

**Expected output:** ten rows (fewer only if arm 2 has fewer than five misses or five false alarms). Contents unknown until run.

**Acceptance criteria**
- [ ] Every miss has `findlay_class == LOF` and `flag_arm2 == 0`; every false alarm has `FUNC` and `flag_arm2 == 1`.
- [ ] Every `guess_why` begins with `Guess:`.
- [ ] A guess mentions zinc only for the eight positions of section 7.10.
- [ ] One miss is chosen as "the explained wrong case" for the write-up, with two sentences on why.

**Commit message:** `stage 6: worst misses and false alarms with labelled guesses`

---

### Stage 7: Heat map

- **Goal:** `outputs/heatmap.html` and `outputs/heatmap.png`.
- **Machine:** any. **Time:** 40 minutes.
- **Files created:** `src/heatmap.py`, the two outputs.

**Steps**
1. Write `src/heatmap.py` (section 8.11).
2. Run it. Because the test run is done, both panels are drawn.
3. Open the PNG and look at it: two panels, 20 rows each, markers visible, colour bar labelled, disclaimer present.

**Verify**
```powershell
.\.venv\Scripts\python.exe -m src.heatmap
Get-Item outputs\heatmap.html, outputs\heatmap.png | Select-Object Name, Length
.\.venv\Scripts\python.exe -m pytest tests/ -q
```

**Expected output:** both files exist; the HTML is a few hundred kilobytes; the PNG is non-empty. In the start panel, the column for position 61 has a blank cell in the `C` row and black LOF dots on `F`, `G`, `R`, `S`, `W`, `Y` only if the lab called those six LOF (check against `load_labels("start", "tuning")`; `C61G` is LOF, verified).

**Acceptance criteria**
- [ ] Each position column has exactly one blank cell, in the normal letter's row.
- [ ] Number of black dots equals the number of LOF rows drawn (398 for both panels); number of grey dots equals the lab-tested non-LOF rows (1,439).
- [ ] Colour range is exactly −15 to 0.
- [ ] Eight zinc triangles are drawn above the start panel, at positions 24, 27, 39, 41, 44, 47, 61 and 64.
- [ ] The HTML opens in a browser with no network errors other than loading Plotly from its CDN.

**Commit message:** `stage 7: heat map (interactive HTML and static PNG)`

*(If stage 7 is started before the test run because the humans asked for parallel work: only the start panel is drawn, 131 black dots and 444 grey dots, and the stage is re-run after stage 5.)*

---

### Stage 8: Gradio demo

- **Goal:** a working demo in `app/`, ready for a Hugging Face Space. No upload.
- **Machine:** any. **Time:** 60 minutes.
- **Files created:** `src/export_app.py`, `app/app.py`, `app/app_data.csv`, `app/app_config.json`, `app/requirements.txt`, `app/README.md`. **Changed:** `tests/test_units.py`.

**Steps**
1. Write `src/export_app.py` (section 8.12) and run it.
2. Write `app/app.py` (section 8.13), `app/requirements.txt` (section 5.4), `app/README.md`.
3. Write the app unit tests and run them.
4. Call `check_mutation` directly on the eight inputs below and paste every output.
5. Start the app locally once to confirm it launches, then stop it.

**Verify**
```powershell
.\.venv\Scripts\python.exe -m src.export_app
.\.venv\Scripts\python.exe -m pytest tests/ -q
.\.venv\Scripts\python.exe -c "import importlib.util as u; s=u.spec_from_file_location('app','app/app.py'); m=u.module_from_spec(s); s.loader.exec_module(m); [print(repr(t), '->', m.check_mutation(t), '\n') for t in ['C61G','c61g','C61g','D2E','X61G','C999G','','hello','C61A','A61G','p.Cys64Gly']]"
.\.venv\Scripts\python.exe app\app.py
```

**Expected output**
- `export_app` prints 6,194 rows; 1,837 rows with a lab class (the test run is done); `end_labels_included: True`.
- `C61G`, `c61g`, `C61g`: the same answer three times; lab result `Loss of function (harmful). Lab function score -1.738.` The model verdict is unknown until run.
- `D2E`: lab result `Functional (harmless). Lab function score 0.089.`
- `X61G`: refusal R2. `C999G`: refusal R5. `''`: refusal R1. `hello`: refusal R2. `A61G`: refusal R6 ending `Did you mean C61G?`
- `C61A`: an answer with `not lab-tested`.
- Every one of the eleven outputs contains both fixed notes.
- Every successful answer has the `On the held-out exam:` line with numbers, because the test run is done.
- `python app\app.py` prints a local address such as `http://127.0.0.1:7860` and **no** public `gradio.live` link.

**Acceptance criteria**
- [ ] All app unit tests pass.
- [ ] No input crashes `check_mutation`.
- [ ] Both fixed notes are in every output.
- [ ] `exam_counts` in `app/app_config.json` equals the nine counts in `outputs/results/band_counts.csv`. Check it with a one-line Python command and paste the output.
- [ ] `app/` holds `app.py`, `requirements.txt`, `README.md` with the header, `app_data.csv`, `app_config.json`.
- [ ] `app.py` imports neither torch nor transformers, and reads nothing outside `app/`.
- [ ] No share link was created.

**Commit message:** `stage 8: Gradio demo and Hugging Face Space files (not uploaded)`

#### 🛑 Checkpoint: before any upload or push

Do not log in to Hugging Face or GitHub, and do not push or upload. Give the operator these steps to do by hand:

1. On huggingface.co, create a new Space: SDK "Gradio", hardware "CPU basic" (free), visibility public.
2. In the Space's "Files" tab, choose "Add file → Upload files" and upload the five files inside `app/` (not the folder itself): `app.py`, `requirements.txt`, `README.md`, `app_data.csv`, `app_config.json`.
3. Wait for the build to finish, then open the Space in a private browser window and try `C61G`, `D2E`, `hello`.
4. Send the Space's address back so it can go into the write-up.

The Space name and the GitHub account name are visible to judges. Round 1 is blind, so the humans should avoid names that give away a school.

---

### Stage 9: Notebook

- **Goal:** `notebook.ipynb`, the version judges open. It runs top to bottom in Google Colab with no GPU, and is saved with its outputs.
- **Machine:** any. **Time:** 60 minutes.
- **Files created:** `tools/build_notebook.py`, `notebook.ipynb`.

**Before you start, ask the humans for the GitHub repository address** (for example `https://github.com/OWNER/REPO`). The notebook needs it. If they do not have one yet, use the placeholder `https://github.com/OWNER/REPO`, and list "replace the repository address in the notebook" in your end-of-stage report.

**How it is built.** `tools/build_notebook.py` uses `nbformat.v4.new_notebook`, `new_markdown_cell` and `new_code_cell` to write `notebook.ipynb`, with `REPO_URL` as a constant at the top of the script. Then the notebook is executed in place.

**Cells, in order** (M = Markdown, C = code):

1. M: title; one-paragraph summary; **the disclaimer**; "BRCA1 is the exam, not the use case" in two sentences; what is already known (ProteinGym, the Evo 2 notebook, AlphaMissense) in three bullets.
2. C: setup. If `google.colab` can be imported: `!git clone --depth 1 REPO_URL repo`, `%cd repo`, `!pip install -q biopython==1.88`. Otherwise assume the notebook is running from the repository root. Then `import sys; sys.path.insert(0, ".")` and `from src import config, load`.
3. M: the data, in plain words: 1,837 lab-tested mutations, two regions, three classes, with the count table.
4. C: `load.check_inputs()` and print `load.class_counts()` as a table.
5. M: how ESM-2 scoring works, in four sentences, with the C61 worked example of the index maths.
6. C: `RECOMPUTE_SCORES = False`. If `True`, call `esm_score.score_positions` (needs a GPU) and keep the result **in memory only**: the notebook never writes `outputs/esm_scores.csv` or any other file. If `False`, `scores = load.load_scores()` and print the row count and model ID. Show the 19 scores at position 61.
7. M: the four arms and the cutoff rule, with the six-row worked example.
8. C: load `outputs/frozen.json` and print each arm's cutoff and the start-region counts.
9. M: the one test run: what was frozen, that it ran once, and the run-log lines (`started`, then `success`).
10. C: display `outputs/results/run_log.txt` and `headline.csv`. Then **recompute the headline table from scratch** with the pure functions (`load.load_labels("end", "post_test")`, `load.lab_tested_scores("end")`, `features.build_feature_table`, `evaluate.apply_frozen_arms`, `evaluate.headline_table`) and `assert` it equals the saved file. This never writes anything and never touches the run log.
11. M: the fixed result sentence(s) from `result_sentences.txt`, copied word for word.
12. C: display `auroc.csv`, `bootstrap.json`, `by_swap.csv`, `int_rows.csv`, `band_counts.csv`.
13. M: the reverse split, labelled secondary.
14. C: display `reverse_split.csv`.
15. M: where it fails: the explained wrong case, with every guess labelled as a guess.
16. C: display `failures.csv`.
17. M: the heat map and how to read it.
18. C: display `outputs/heatmap.png` with `IPython.display.Image`.
19. M: limits (all of section 12.4); the groups-check sentence (section 12.2); sources and credits; **the disclaimer again**.

The notebook must not: call `evaluate.main`, write to `outputs/`, contain a school name, or contain any number typed by hand that is not read from a file.

**Verify**
```powershell
.\.venv\Scripts\python.exe tools\build_notebook.py
$env:CUDA_VISIBLE_DEVICES = ""
.\.venv\Scripts\python.exe -m jupyter nbconvert --to notebook --execute --inplace notebook.ipynb
Remove-Item Env:\CUDA_VISIBLE_DEVICES
git status --short
```

**Expected output:** nbconvert finishes with no error; `git status --short` shows `notebook.ipynb` and `tools/build_notebook.py` only (no file in `outputs/` or `data/` changed).

**Acceptance criteria**
- [ ] The notebook executes top to bottom with the GPU hidden.
- [ ] The recomputed headline equals the saved one (the `assert` in cell 10 passes).
- [ ] The disclaimer appears word for word.
- [ ] Outputs are saved inside the notebook file.
- [ ] Nothing in `outputs/` changed.

**Commit message:** `stage 9: notebook built and executed`

After the commit, the `repro-checker` review is due (section 11). It must report WORKS. The Colab clone step can only be tested after the humans push to GitHub; say so.

---

### Stage 10: Write-up and form answers

- **Goal:** `writeup/writeup.md` (about one page), `writeup/form_answers.md`, and `README.md` replaced by the write-up.
- **Machine:** any. **Time:** 60 minutes.
- Follow section 12 exactly. Every number is copied from a file in `outputs/results/`; next to each number in your draft, keep an HTML comment naming the file and column it came from (for example `<!-- headline.csv arm2 harmful_caught -->`), so the claims check is fast.

**Verify**
```powershell
Select-String -Path writeup\writeup.md, README.md -Pattern "This project is a prototype and is not medical advice." -SimpleMatch
Select-String -Path writeup\*.md, README.md -Pattern "\[X", "\[Y", "TODO", "OWNER/REPO"
```

**Expected output:** the first command finds the line in both files. The second finds nothing, except placeholders the humans must fill (team names, emails, schools, links), which you list in your report.

**Acceptance criteria:** the checklist in section 12.7.

**Commit message:** `stage 10: write-up and form answers`

After the commit, the `claims-checker` review and then the `judge` review are due (section 11). `claims-checker` must report CLEAN. `judge` runs once; its ten questions feed stage 12.

---

### Stage 11: Video script

- **Goal:** `writeup/video_script.md`, a two-minute script. **Time:** 20 minutes.
- Structure and rules in section 12.5.

**Acceptance criteria**
- [ ] Four parts: 15, 60, 30, 15 seconds. About 280–320 spoken words in total.
- [ ] The demo part types `C61G`, shows the verdict next to the lab result, and shows the explained wrong case.
- [ ] Every number matches `outputs/results/`. Nothing from section 3.5.
- [ ] The disclaimer is spoken or shown on screen.

**Commit message:** `stage 11: two-minute video script`

Then the `claims-checker` review is due again (section 11).

---

### Stage 12: Q&A prep

- **Goal:** `writeup/qa_sheet.md`. **Time:** 60 minutes.

**Steps**
1. Part 1: a table with one row per locked rule (rules 1–15 of section 3): the rule in one line, the reason in one line, and the judge question it answers.
2. Part 2: the ten hardest questions from the `judge` review (the humans paste them in), each with a strong answer of two to four sentences built only from what the submission shows. If the judge's questions are not available yet, use the ten in section 12.8.
3. Part 3: quiz the human at the keyboard with the ten questions, **one at a time**: ask, wait for their answer, then say what was strong and what was missing. Do not show the model answer first.

**Acceptance criteria**
- [ ] Fifteen rule rows.
- [ ] Ten questions with answers; no answer claims anything in section 3.5.
- [ ] The quiz was run, or the human chose to skip it.

**Commit message:** `stage 12: Q&A sheet`

---

## 11. Reviews

Six review agents are defined in `.claude/agents/`. They are Claude Code subagents: read-only reviewers that start with a fresh context.

- **If you are Claude Code:** run each review yourself, as a subagent, at the moment listed below. Paste its full report into your end-of-stage report and act on it. Never edit a file in `.claude/`, and never do a reviewer's job in your own context instead: the point is a second pair of eyes that has not seen your reasoning.
- **If you are any other tool:** you cannot run them. Commit, then tell the operator which review is due, so the reviewer can run it in Claude Code on a copy of the repository (the humans handle pushing or copying; never push yourself without being asked). Do not wait for `test-runner`. Do wait for every other review that this spec says must report before you continue.

| Agent | When | What must be committed first | What it must report |
| --- | --- | --- | --- |
| `test-runner` | After every stage | The stage's commit | `ALL PASS` (skips are listed with reasons) |
| `index-debugger` | Only if sanity test 3, 4 or 5 fails | The failing output, pasted in full, and the commit containing `src/esm_score.py` | `BUG FOUND` with the exact line and smallest fix, or `NO BUG FOUND` |
| `leak-auditor` | Before checkpoint B, and again after any later change to anything in `src/` | The stage-4 commit (`frozen.json` present, `run_log.txt` absent) | `CLEAN` — required before "run the test" |
| `repro-checker` | After stage 9, and again right before submitting | The stage-9 commit; everything the notebook and app read must be committed | `WORKS` |
| `judge` | Once, after stage 10 | The stage-10 commit: `writeup/`, `app/`, `outputs/results/`, `notebook.ipynb` | Scores, weak points and ten hard questions (these feed stage 12) |
| `claims-checker` | After stage 5b (version zero), after stage 10, after stage 11, and right before submitting | The commit with the changed text | `CLEAN` — required before submitting |

If a review finds a problem: fix exactly what it names, re-run the affected stage's verify list, commit with a message starting `fix:`, and ask for the review to be run again. A fix to `src/` after stage 5 must not change any stage-5 result file. If it would, stop: that is a bug-fix rerun and needs a human's decision and reason.

`test-runner` runs `.venv/Scripts/python.exe -m pytest tests/ -q`, so if reviews run on a second machine, its `.venv` must have the packages in `requirements.txt` installed.

---

## 12. Write-up, form, video and Q&A

### 12.1 The organisers' template (use these headings, in this order)

```
# BRCA Spellcheck: an AI second opinion on uncertain BRCA1 mutations

**Track:** Health
**Team:** [names — the humans fill this in]

*This project is a prototype and is not medical advice.*

**In one sentence:** [the main fixed sentence, then sentence S6]

## The problem
## What we built
## What we added
## How it works
## Results
## What we'd do next
## Sources
```

About one page (roughly 450–600 words, plus the one results table). What each section must contain:

- **In one sentence** (directly under the disclaimer). The main fixed sentence (S1 or S2; if neither applies, the plain statement the humans approve) followed by S6, both word for word from `result_sentences.txt`. This is the one takeaway a judge should leave with. Nothing else goes on this line.
- **The problem** (2–3 sentences). Who: a genetic counselor. The moment: sitting with a patient whose test came back "variant of uncertain significance". Why it matters: most genes have no lab map, so the counselor leans on computer predictions and needs to know how far to trust them.
- **What we built.** What goes in (a mutation such as `C61G`), what comes out (a verdict, a score, and the lab's answer next to it), with one real example copied from the demo's output. Then the framing sentence: BRCA1 is the exam, not the use case; at these positions the lab has already measured 96.5% of possible single-letter DNA changes, so a counselor would use the lab result there. What the exam tells us is how far to trust this kind of model on genes with no lab map.
- **What we added.** Start with the words `Our own idea.` Then prior work in **exactly two sentences, and nowhere else in the write-up**: (1) grading a protein language model against this lab data has been done before (ProteinGym; the Evo 2 BRCA1 notebook); (2) stronger tools such as AlphaMissense exist, and we do not claim to match them. Then the four things we added: a whole region held out with the cutoff frozen in advance; a test of whether labels from one end of the protein help at the other; missed cases and false alarms counted separately; and a worked example of where it fails.
- **How it works.** Data: Findlay et al. 2018 via ProteinGym; 1,837 lab-tested single-letter protein changes in two regions; classes from the study's own cutoffs; 155 in-between rows left out of tuning and the headline. "Cleaning" is: rebuilding the normal sequence from the table and checking every row against it. Model: ESM-2 (650M), zero-shot, with the hide-a-letter explanation in one sentence; the model read the BRCA1 sequence during its training but never any lab label. Baselines: always-harmless, and the BLOSUM62 table. The split: tuned on positions 1–101, tested once on 1631–1855, never at random, and why.
- **Results.** The headline table from `headline.md` (four arms: caught, missed, false alarms, catch rate, AUROC). The fixed sentence(s) from `result_sentences.txt`, word for word. Sentence S6 with both numbers together. One sentence on the groups check (12.2). One explained wrong case from stage 6, with the guess labelled as a guess. One sentence on **what didn't work** (12.2). How many times the test was run.
- **What we'd do next.** Two concrete steps: add a splicing predictor for the DNA-level blind spot; repeat the exam on a second gene that has a lab map.
- **Sources.** Every dataset, model, library, method paper and AI tool, with links, and what each AI tool was used for (12.6).

### 12.2 Requirements from the organisers' documents that are easy to miss

1. **Round 1 is blind.** Reviewers see the project but not names or schools. So: team members' names appear only on the `Team:` line of the write-up and on the submission form. No school name anywhere. No names in the notebook, the demo, code comments, commit messages you write, the README title, or the video script. Commit author names and addresses are public on GitHub too, so this repository's local git configuration uses a neutral identity (section 0, point 13). The GitHub account name and the Space name are visible as well.
2. **"What didn't work"** must be in the write-up. Use a true one from the results: for example the arm 3 result if S3 applies, a cutoff that did not carry over if S5 applies, or arm 1's flood of false alarms. Do not invent one.
3. **The groups check.** Health-track judges look for whether the model works equally well across groups. Write: results are broken down by protein region (the main test and the reverse test) and by kind of swap (`by_swap.csv`); say plainly where they differ. Then: "Our data are mutations, not people, so we cannot check performance by age, sex or ancestry."
4. **The user's moment** must be named: the counselor with a VUS result in front of a patient.
5. **AI tools with what each was used for**, in the write-up and on the form.
6. **Missed cases and false alarms reported separately**, not only accuracy.
7. **The baseline's score next to the model's score.**
8. **At least one case the model gets wrong, with why.**
9. **Results from data the model did not see during tuning.**
10. **Every link must open in a private browser window** without signing in. The humans check this.
11. **The exact disclaimer line** in the write-up.
12. Submission needs: a code link (GitHub repository or Colab notebook, outputs visible), the write-up (the GitHub README is acceptable), credits. The video is optional but recommended.

### 12.3 The fixed result sentences

Section 3.4, chosen by the rules in section 8.9 step 10. Copy them from `outputs/results/result_sentences.txt`. Never soften or strengthen them. Never add adjectives such as "strong", "impressive" or "accurate" to a result.

### 12.4 Limits (two go in the one-page write-up; the rest go in the notebook)

- **Cells in a dish, not people.** The lab truth measures whether cells survive, which is strong but indirect evidence about cancer risk. A 2025 follow-up preprint from the same group reports that results can depend on the cell type used.
- **Blind to splicing.** The lab edited DNA. Some changes are harmful because they break how the gene's message is cut and joined, not because of the letter swap. A protein model cannot see that. (Say "we could not test this" unless stage 6's optional analysis was done.)
- **Small coverage.** 326 of 1,863 positions, single-letter protein changes only, one gene.
- **The model has read the BRCA1 sequence** (not the labels) during its training.
- **The cutoffs are ours.** The 90% rule is a choice; a different rule changes the trade-off.
- **Stronger tools exist** (for example AlphaMissense). As far as we know its predictions are also free to download, so we do not claim to be better or cheaper.
- **Not new as a method.** ProteinGym already scores ESM-2 on this table, and the Evo 2 notebook runs a DNA model on the same study.
- **Redundant where it is tested.** At these positions the lab result already exists.

Use two in the write-up, in the Results section: blind to splicing; cells in a dish. Prior work and stronger tools are already covered by the two sentences in "What we added", so do not repeat them as limits. Everything else on this list goes in the notebook only.

### 12.5 The video script (two minutes)

| Part | Length | Content |
| --- | --- | --- |
| The problem | about 15 s | The counselor, the patient, the "uncertain" result. |
| The demo | about 60 s | Type `C61G`; show the verdict next to the lab result. Type a harmless one (`D2E`). Show the explained wrong case. Say the exam-not-use-case sentence. |
| How it works | about 30 s | It read proteins, never patients. Graded against CRISPR lab results. Compared with two baselines on a held-out region: give the headline numbers for arm 0, arm 1 and arm 2, misses and false alarms both. |
| Next steps | about 15 s | One limit (blind to splicing) and one next step. Show or say the disclaimer. |

A screen recording with a voice-over. No names or school on screen.

### 12.6 Credits (start `writeup/credits.md` with all of these)

**Data**
- Findlay GM et al. (2018). Accurate classification of BRCA1 variants with saturation genome editing. *Nature* 562:217–222. https://doi.org/10.1038/s41586-018-0461-z — the lab measurements. Function scores are the authors' (University of Washington).
- University of Washington BRCA1 saturation genome editing site (classes and cutoffs): https://sge.gs.washington.edu/BRCA1/
- Protein Data Bank entry 1JM7, the solution structure of the BRCA1/BARD1 RING-domain heterodimer (the eight zinc-holding letters, section 7.10): https://www.rcsb.org/structure/1JM7
- ProteinGym (Notin P et al., 2023), the benchmark our copy of the table comes from, file `BRCA1_HUMAN_Findlay_2018.csv`: https://github.com/OATML-Markslab/ProteinGym — also **prior work**: it already scores ESM-2 on this table.

**Models and methods**
- ESM-2 (Lin Z et al., 2023, *Science* 379:1123–1130), model `facebook/esm2_t33_650M_UR50D`: https://huggingface.co/facebook/esm2_t33_650M_UR50D and https://github.com/facebookresearch/esm
- The hide-a-letter scoring method (masked-marginal scoring): Meier J et al. (2021), "Language models enable zero-shot prediction of the effects of mutations on protein function", NeurIPS. https://doi.org/10.1101/2021.07.09.450648
- BLOSUM62: Henikoff S, Henikoff JG (1992), *PNAS* 89:10915–10919, loaded through Biopython: https://biopython.org
- Kyte J, Doolittle RF (1982), *Journal of Molecular Biology* 157:105–132 — hydrophobicity scale.
- Zamyatnin AA (1972), *Progress in Biophysics and Molecular Biology* 24:107–123 — residue volumes.

**Prior work we do not use but must cite**
- Evo 2 BRCA1 zero-shot notebook (Arc Institute; also published by NVIDIA): https://github.com/ArcInstitute/evo2 and https://docs.nvidia.com/bionemo-recipes/2.7.1/main/examples/bionemo-evo2/zeroshot_brca1/
- AlphaMissense (Cheng J et al., 2023, *Science* 381:eadg7492) — a stronger tool; named as a limit.
- 2025 follow-up preprint on cell types: https://www.medrxiv.org/content/10.1101/2025.08.11.25333423.full.pdf

**Libraries:** PyTorch, Hugging Face Transformers, NumPy, pandas, SciPy, scikit-learn, Biopython, Plotly, Matplotlib, Gradio, pytest, Jupyter (nbformat, nbconvert).

**AI tools, with roles**
- **Claude (Anthropic), in chat and in Claude Code:** planning, writing and auditing this build specification, and review (running the tests, auditing for test-set leakage, checking the index maths, checking claims in the write-up, a reproducibility check, and a mock judge).
- **The coding agent that wrote and ran the project code from the specification.** Name the tool and model actually used, exactly as the operator gave it in the first message (for example Claude Code with a named Claude model, or a named Gemini model in Google Antigravity). Do not guess, and do not credit a tool that was not used.
- Any other AI tool either human used this weekend (ask them; for example a chatbot for brainstorming the idea). List each with its role.

**The organisers:** the Veritas AI hackathon project briefs, Hacker Guide and sample submission.

Links marked above are as recorded by the team; a human should open each once before submitting.

### 12.7 Write-up acceptance checklist

- [ ] The seven headings of 12.1, in order; about one page.
- [ ] The disclaimer, word for word.
- [ ] The exam-not-use-case framing, with 96.5%.
- [ ] The one-sentence takeaway directly under the disclaimer.
- [ ] Prior work in exactly two sentences (ProteinGym and the Evo 2 notebook; AlphaMissense), not repeated elsewhere.
- [ ] "Our own idea." and the four things we added.
- [ ] The headline table with all four arms; misses and false alarms in separate columns.
- [ ] The fixed sentence(s), word for word; S6 with both numbers.
- [ ] One explained wrong case; the guess is labelled as a guess.
- [ ] What didn't work.
- [ ] The groups-check sentences.
- [ ] Two limits (section 12.4).
- [ ] The number of test runs.
- [ ] Every AI tool with its role; every dataset, model and library with a link.
- [ ] No name outside the `Team:` line; no school anywhere.
- [ ] Nothing from section 3.5 (new, better, cheaper, diagnoses, clinical classification).
- [ ] No leftover placeholder such as `[X of 267]`.
- [ ] Every number traceable to a file in `outputs/results/`.

`writeup/form_answers.md` holds one answer per form field:

| Form field | Answer |
| --- | --- |
| Project title | BRCA Spellcheck: an AI second opinion on uncertain BRCA1 mutations |
| Building as | Team |
| Team name | [the humans fill in] |
| Team members | [names, emails, schools — the humans fill in] |
| Track | AI + Health |
| Starting point | Our own idea |
| One-sentence pitch | BRCA Spellcheck helps a genetic counselor judge an "uncertain" mutation by taking the mutation and giving back whether a protein-reading AI thinks it looks harmful, with that AI graded in the open against CRISPR lab results for BRCA1. |
| Code link | [GitHub repository address] |
| Write-up link | [the README on GitHub, or a shared document] |
| Demo video link | [the humans fill in, or left blank] |
| What we added | The "What we added" paragraph of the write-up. |
| Datasets, models and AI tools used | One line from 12.6, including each AI tool's role. |
| New to AI? | [the humans answer] |

### 12.8 Ten likely judge questions (fallback list for stage 12)

1. How do you know it works?
2. The lab already tested these mutations. Who needs your tool?
3. Did the model see BRCA1 during training?
4. Why not split the data at random?
5. Isn't this already on the ProteinGym leaderboard?
6. How is this different from AlphaMissense?
7. Why 90%? What happens with a different cutoff?
8. What would make it fail?
9. Would you use it on a patient?
10. You used AI tools to build this. Which parts can you explain yourselves?

---

## 13. Known traps

Each trap is followed by how it is caught.

**Index maths**
1. **Off-by-one token index.** Position `p` in window A is token `p`, not `p − 1`: index 0 is the start token. → `run_letter_test`; unit test `test_token_index_worked_examples`.
2. **Wrong window for a position.** End-region positions must use window B. A position such as 1775 does not exist in window A. → `window_of` raises outside the two regions; the letter test.
3. **Forgetting the start token**, or assuming an end token is absent. A 500-letter window is 502 tokens. → the letter test asserts the length and both special tokens.
4. **Using the protein position as a Python index.** `seq[61]` is position 62. Position `p` is `seq[p − 1]`. → sanity test 2.
5. **Window B's offset.** `token_index = position − 1363`, not `− 1364` and not `− 1362`. → the worked example (1775 → 412).
6. **Reading logits at the wrong index**, for example at index 0 or at `token_index − 1`. Read at the same index that was masked. → the assert before masking; sanity tests 4 and 5.
7. **Masking in a shared tensor** so that masks pile up across positions. Each row of a batch starts from a fresh copy of the unmasked window and has exactly one mask. → assert `(batch == mask_id).sum(dim=1)` is all ones.
8. **Tokenizing with spaces between letters**, or tokenizing a list of letters. Pass one plain string. → the letter test's length check.

**Signs and directions**
9. **ESM score sign flipped.** It is new minus normal. Most values are negative. → sanity test 4 (C61G below zero); sanity test 5 (positive correlation with the function score).
10. **Harm score direction.** Arm 1 and arm 2 harm scores are the **negated** raw values. Arm 3 uses the **raw** ESM and BLOSUM values as features. → `test_features_c61g`.
11. **AUROC label direction flipped.** `y_true` is 1 for LOF; `y_score` is the harm score. An AUROC below 0.5 for arm 2 almost certainly means a flipped sign, not a bad model. → `test_auroc_direction`.
12. **`predict_proba` column.** Column 1 is harmful only if `classes_` is `[0, 1]`. → the check in `fit_arm3`.
13. **Negative zero.** `-0.0` prints oddly and reads as a bug. → `0.0 - x`; `test_harm_arm1_no_negative_zero`.

**The cutoff rule**
14. **Cutoff chosen as "lowest" instead of "highest".** The lowest cutoff flags everything. The rule is the highest cutoff that still catches 90%. → worked example 1; the stage-2 verified numbers.
15. **Ties split across a cutoff**, for example by sorting and taking the first 118 rows. Flagging is by value (`>=`), never by rank. → `test_catch_cutoff_matches_brute_force`.
16. **`>` instead of `>=`.** → `test_flag_is_inclusive`.
17. **Floating-point target.** `0.9 * n` then `ceil` can be off by one. Use integer arithmetic. → `test_catch_cutoff_needed_counts`.
18. **Cutoff B direction.** It is the **lowest** cutoff that flags at most 10% of harmless rows. → worked example 2.

**Leakage**
19. **Random row splits** (`train_test_split`, `KFold`, `StratifiedKFold`). Rows at one position are near-copies. → only `StratifiedGroupKFold` with `groups = position`; `test_fit_arm3_groups_do_not_leak`.
20. **Scaler fit on all data** before cross-validation, or on both regions. → the scaler lives inside the pipeline and is fit per fold; the frozen scaler is fit on start-region rows only.
21. **Arm 3 cutoff chosen in-sample.** It comes from out-of-fold probabilities. → `cutoff_source` in `frozen.json`; the two separate count pairs.
22. **INT rows leaking into tuning.** → `_start_tuning_table` asserts 516 rows.
23. **Peeking at end-region labels through a plot, a print, a `describe()` or a notebook cell.** → `load_labels` permission check; `heatmap` and `export_app` gated on `test_run_done()`; the leak audit before checkpoint B.
24. **Checking M1775R or A1708E "just to be sure" before the test.** They are end-region labels. → not allowed; sanity test 7 skips until the run is done.
25. **A full-table correlation before the test.** It includes the end region. → sanity test 5 uses the start region only.
26. **Tuning anything after the test run** (a cutoff, a colour range chosen to flatter, a feature). → stage 5 commits at once; `freeze` refuses once the test has started (`test_attempted()`).
27. **Running `evaluate.py` "to see if it works".** → use `--selfcheck-start`; the run log makes any real run visible.

**The bootstrap**
28. **Resampling rows instead of positions.** → `test_bootstrap_resamples_positions`.
29. **Dropping duplicates** after resampling positions. A position drawn twice counts twice. → the same test.
30. **Resamples with one class** crash `roc_auc_score`. → skip and count them.
31. **Using the legacy global NumPy seed** for the bootstrap. Use `numpy.random.default_rng(42)`.

**Data handling**
32. **pandas reading `position` as float** (after a merge that introduces NaN, or from a CSV with blanks). → explicit dtypes; `validate="one_to_one"` on merges; row-count checks after every merge.
33. **pandas reading the letter-pair `NA` as missing.** A column holding amino-acid letters is safe, but any string column can be hit. → `keep_default_na=False` for string columns.
34. **Sorting letters alphabetically by accident and it happening to match.** Always sort by `AA_ORDER.index`.
35. **Indexing BLOSUM62 by number.** Its alphabet order is not `AA_ORDER`. → index by letter; spot values.
36. **Assuming all 19 swaps at a position were lab-tested.** Only 1 to 7 were. `C61A` has a score but no lab result. → `is_lab_tested`; the app's `not lab-tested`.

**Files and platform**
37. **Windows paths.** Use `pathlib.Path`; never build paths with backslashes or string joins.
38. **Encoding.** Open every text file with `encoding="utf-8"`. Windows defaults to another encoding, which breaks on `Å`, `−` and the emoji in `app/README.md`.
39. **CRLF line endings changing a file hash.** → `newline=""` and `lineterminator="\n"` when writing; `.gitattributes`.
40. **The hash computed before the file is fully written.** → write to a `.tmp` file, close, `os.replace`, then hash from disk.
41. **Excel or an editor re-saving a CSV.** It changes bytes and number formats. Never open `esm_scores.csv` in Excel.
42. **NumPy types in JSON.** `json.dump` fails on `numpy.float64`. Convert first.
43. **PowerShell is not bash.** No `&&` in Windows PowerShell 5; use `;`. Environment variables are `$env:NAME = "value"`. No `export`, no `source`.
44. **Running a module by path** (`python src/load.py`) breaks imports. Use `python -m src.load`.

**The demo**
45. **Gradio launching with a share link** (`share=True`) publishes the app through a third-party tunnel. Never pass it.
46. **The app loading the model** when it should only read the CSV. → `test_app_imports_no_model`.
47. **The app reading from `outputs/` or `data/`.** On a Space only `app/` exists. → paths relative to `app.py`.
48. **A refusal that forgets the fixed notes.** → `test_app_notes_always_present`.
49. **User text injected into Markdown.** → strip the listed characters in R2.

**Honesty**
50. **A number typed from memory into the write-up.** → every number comes from a results file; the claims check.
51. **Describing in-sample start-region numbers as results.** Only stage 5's files are results.
52. **Softening a weak result or dropping it.** Weak results are reported as they are.

---

## 14. Cut list and never-cut list

**If we are behind, cut in this order, and only after telling the humans:**

1. Extra features in arm 3 (keep only `esm_score` and `blosum_score`). This is allowed **only before `python -m src.freeze` has been run for the first time**, and only to save time. It must never be decided after seeing any arm 3 number: that would be a settings search (rule 7). It needs a human's approval, because `FEATURE_ORDER` and the checks on it would have to change together.
2. The reverse split.
3. The "Can't tell" band (the demo then uses cutoff A alone).
4. The interactive heat map (keep the static PNG).
5. The video.

**Never cut:**

- the two baselines (arm 0 and arm 1);
- the held-out region;
- the separate counts of missed cases and false alarms;
- one explained wrong case;
- the disclaimer.

Never add scope beyond this spec without asking.

---

## 15. Troubleshooting

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `torch.cuda.is_available()` is `False` on a machine with an NVIDIA card | The CPU build of torch was installed (for example by `pip install -r requirements.txt` before the GPU install) | `pip uninstall torch -y`, then the command in section 5.2. Verify with the four-line check. |
| `CUDA error: no kernel image is available` or a warning that the GPU is not supported | The PyTorch build is too old for this GPU. The RTX 5060 needs a CUDA 12.8 build (`cu128`) | Reinstall with the exact command in 5.2. Do not try other index URLs without asking. |
| `torch.cuda.is_available()` is `False` and there is no NVIDIA card | Expected | Section 5.5. Every stage except 3 still works. |
| `CUDA out of memory` | Other programs are using the GPU, or the batch is too large | Close GPU-heavy programs; rerun with `--batch-size 1`. If it still fails, stop and ask (section 8.5). |
| Model download fails or stalls | Network, or a Hugging Face outage | Report the exact error. Retrying the same command resumes the download. Do not switch model or mirror without asking. No login is needed for this model. |
| A warning about symlinks from `huggingface_hub` on Windows | Windows without Developer Mode | Harmless. Ignore it. |
| Tokenizer gives a length other than 502 | The window was passed with spaces, as a list, or is not 500 letters | Pass `seq[0:500]` or `seq[1363:1863]` as one plain string. Check `len(window_seq) == 500`. |
| The letter test fails | An index bug | Stop. See the eight checks below. |
| `substitution_matrices.load("BLOSUM62")` fails | Biopython missing or a different version | `pip show biopython` must say 1.88. Reinstall from `requirements.txt`. |
| `ConvergenceWarning` from `LogisticRegression` | lbfgs did not converge in 100 iterations | Do not raise `max_iter` (rule 7). Record `converged: false` and report at checkpoint B. With eight standardised features this is unlikely. |
| A `FutureWarning` about `penalty` | Passing `penalty=` to `LogisticRegression` in scikit-learn 1.8 | Do not pass it. Use `LogisticRegression(random_state=42)`. |
| Plotly image export asks for Kaleido or Chrome | Something called `fig.write_image` | Do not use it. The PNG is drawn with matplotlib (section 8.11). |
| matplotlib tries to open a window, or errors about a display | The default backend | `matplotlib.use("Agg")` before importing `pyplot`. |
| Gradio: `OSError: Cannot find empty port` or the port is in use | An earlier app is still running | Stop the earlier process, or launch with `demo.launch(server_port=7861)` for a local check only. |
| pytest: `ModuleNotFoundError: No module named 'src'` | Run from the wrong folder, or `pytest.ini` missing | Run from the repository root; check `pytest.ini` (section 5.6) and `src/__init__.py`. |
| `git` refuses a large file or the commit is huge | The 1 GB ProteinGym folder or `.venv/` was added | Check `.gitignore` is intact; `git rm -r --cached <path>`; never commit `data/DMS_ProteinGym_substitutions/`. |
| `HASH MISMATCH` from `evaluate` | `esm_scores.csv` changed after the freeze: recomputed, re-saved, or line endings converted | Do not re-freeze to make it pass. Find out what changed (`git diff --stat`, `git log -- outputs/esm_scores.csv`), restore the committed file with `git checkout -- outputs/esm_scores.csv`, and report. |
| `PermissionError: End-region lab labels are locked …` | Code asked for end-region labels before the test run | That code path must not run yet. Do not change the purpose string to get around it. |
| nbconvert: `No such kernel named python3` | The kernel is not registered | `python -m ipykernel install --user --name python3` |
| PowerShell says running scripts is disabled | Trying to activate the environment | Do not activate. Call `.\.venv\Scripts\python.exe` directly. |

**If the letter test (or sanity test 4 or 5) fails, these are the eight checks the debugging review makes, in order.** You may run them yourself to gather evidence, using start-region labels only, but do not edit code by trial and error:

1. The FASTA sequence is 1,863 letters and starts `MDLSALRVEEVQNVINAMQK`.
2. Window A is `seq[0:500]` and window B is `seq[1363:1863]`. Print the first and last five letters of each (`MDLSA…RPLTN`, `NLGEA…PHSHY`).
3. The tokenizer adds exactly one start token at index 0 and one end token at the end. Print the first six tokens of each window. Check no letter is merged or split.
4. For five positions in each region (first, last and three in between), print: protein position, window, token index, the token at that index, and the `wild_type` letter from the CSV. They must match.
5. The mask is placed at the same index the logits are read from.
6. The sign: score = log P(new) − log P(normal), so a harmful swap should be negative.
7. Letters are looked up with `convert_tokens_to_ids` for single upper-case letters, and none maps to the unknown token.
8. Start region only: the Spearman correlation of the ESM score with `function_score`, and where `C61G` and `C64G` rank.

---

## 16. Definition of done

Everything below must be true before the humans submit.

**Code and tests**
- [ ] `python -m pytest tests/ -q` passes with no failures. The only skip that may remain is 3b on a machine with no network.
- [ ] `python -m src.load --check` and `--rebuild-check` pass.
- [ ] No file in `data/` has changed since the first commit.

**The test**
- [ ] `outputs/frozen.json` was committed before `outputs/results/run_log.txt` existed (visible in `git log`).
- [ ] `leak-auditor` reported CLEAN before the run, and again after the last change to `src/`.
- [ ] `run_log.txt` shows one `status=started` line and one `status=success` line, or every extra run has a real bug-fix reason and the write-up states the number of runs (the count of `status=started` lines).
- [ ] `esm_scores.csv` still matches the hash in `frozen.json`.

**Results**
- [ ] `outputs/results/` holds: `headline.csv`, `headline.md`, `auroc.csv`, `bootstrap.json`, `band_counts.csv`, `by_swap.csv`, `int_rows.csv`, `reverse_split.csv`, `predictions_end.csv`, `result_sentences.txt`, `sanity7.json`, `run_log.txt`, `failures.csv`, `failures.md`.
- [ ] Both baselines are in the headline table next to arm 2 and arm 3.
- [ ] Misses and false alarms are separate columns.
- [ ] One wrong case is explained, with the guess labelled.

**Deliverables**
- [ ] `notebook.ipynb` runs top to bottom with no GPU and is saved with outputs. `repro-checker` reported WORKS.
- [ ] `outputs/heatmap.png` exists (and `heatmap.html`, unless cut).
- [ ] The demo answers `C61G`, a harmless mutation, an unreadable input and an untested position without crashing, with both fixed notes every time.
- [ ] `app/` holds the five files a Space needs.
- [ ] `writeup/writeup.md` passes the checklist in 12.7 and is also the `README.md`.
- [ ] `writeup/form_answers.md`, `writeup/video_script.md`, `writeup/qa_sheet.md`, `writeup/decision_log.md`, `writeup/credits.md` exist.
- [ ] `claims-checker` reported CLEAN on the final text.

**Words**
- [ ] "This project is a prototype and is not medical advice." appears word for word in the write-up, the README, the notebook and the demo.
- [ ] The exam-not-use-case framing appears in the write-up, the notebook, the demo and the video script.
- [ ] ProteinGym and the Evo 2 notebook are cited as prior work; AlphaMissense is named as stronger.
- [ ] Nothing claims the method is new, better or cheaper than another tool, diagnostic, or a clinical classifier.
- [ ] Every dataset, model, library and AI tool is credited with a link, and each AI tool has its role.
- [ ] No school name anywhere; names only on the `Team:` line and the form.
- [ ] No placeholder is left (`[X of 267]`, `OWNER/REPO`, `TODO`).

**Things only the humans do**
- [ ] Push to GitHub; create and upload the Space; record and upload the video.
- [ ] Replace the repository address in the notebook, rebuild and re-execute it, and open it once in Colab.
- [ ] Open every link in a private browser window.
- [ ] Fill in the team name, names, emails, schools and "New to AI?" on the form.
- [ ] Look up ProteinGym's published ESM-2 650M Spearman value for this table. (The zinc positions are already confirmed: section 7.10.)
- [ ] Submit version zero as soon as stage 5b is committed, then resubmit the final version.
- [ ] Submit version one by 8:00 PM Central on Sunday.
