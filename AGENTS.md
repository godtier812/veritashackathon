# Rules for any coding agent in this repository

Project: FoldFlux, a hackathon entry. Deadline Sun Nov 15 2026, 10:59 PM US Central.

**`BUILD_SPEC.md` is the spec.** Read it in full before changing anything. If anything disagrees with it, it wins. Never change a locked rule (its section 3) without asking a human.

## The five rules that matter most

1. **Tune on the start region (positions 1–101). Test on the end region (1631–1855) exactly once.** Never split rows at random.
2. **No peeking.** Before the test run, no code may print, plot or compare end-region lab labels with any score. Counts already written in the spec are the only exception.
3. **Freeze first.** Every cutoff and fitted number goes into `outputs/frozen.json` before the test run. After that, `frozen.json` and `outputs/esm_scores.csv` must not change by one byte.
4. **Never invent a number.** Every result comes from a file in `outputs/results/`. Use a placeholder such as `[X of 267]` until it exists.
5. **Never move past a failing test,** and never weaken, skip or delete a test to make it pass.

## Checkpoints: stop completely and wait for the exact words

- **A**, after stage 3 (ESM-2 scores): wait for "go A".
- **B**, after stage 4 (`frozen.json` written, leak audit CLEAN): wait for "run the test".
- **C**, after stage 5 (the one test run): wait for "go C".
- Before any Hugging Face upload or git push: ask first.

## Always

- Never run `python -m src.evaluate` without "run the test". The only exception is `--selfcheck-start`.
- Never push to a remote, log in to any service, or upload anything without asking.
- Never write to `data/`. Never create or edit `outputs/results/run_log.txt` by hand.
- Never claim the method is new, that it beats or is cheaper than another tool, or that it diagnoses or classifies variants for clinical use.
- Keep this line word for word in the write-up, notebook and demo: "This project is a prototype and is not medical advice."
- No names outside the write-up's `Team:` line, and no school names anywhere (round 1 judging is blind).
- Do not read, edit or commit `CLAUDE.md`, `MASTER_PROMPT.md`, `PROJECT_BRIEF.md` or `.claude/`.
- Work stage by stage, commit after each stage, and report actual command output.
- Windows and PowerShell. Call Python as `.\.venv\Scripts\python.exe` from the repository root.
