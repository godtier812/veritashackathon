# First message to Gemini in Antigravity

Paste everything below the line as your first message.

---

You are the coding agent for a two-person hackathon team. The project is called BRCA Spellcheck. The deadline is tonight, Sunday Nov 15 2026, 10:59 PM US Central.

The complete specification is `BUILD_SPEC.md` in the repository root. It is long (about 2,800 lines) and it stands alone. Standing rules are in `AGENTS.md`.

Do these four things, in this order, and nothing else yet:

1. **Read `BUILD_SPEC.md` in full, top to bottom.** Do not skim, and do not stop at the table of contents. Then read `AGENTS.md`.

2. **Restate the locked rules in your own words.** Cover all fifteen rules in section 3, plus:
   - what counts as "peeking" at end-region labels, and which three things you may do with the end region before the test run;
   - exactly how a cutoff is chosen, using the six-row worked example in section 8.3;
   - how a protein position becomes a token index, using position 1775 as the example;
   - the three checkpoints, and the exact words you must wait for at each;
   - the list of things you must never claim.

3. **List every ambiguity you see:** anything in the spec that could be read two ways, anything that looks contradictory, and anything you would otherwise have to guess. For each one, say which reading you would choose and why. If you see none, say so. Do not resolve an ambiguity about a locked rule by yourself: ask.

4. **Then start stage 1 only** (section 10). Before writing any file, print the stage goal, the files you will create, and the current time. When stage 1's acceptance criteria all pass, commit with the stage's commit message, give the end-of-stage report, and stop for me to read it.

Things to remember from the first minute:

- Do not read, edit or commit `CLAUDE.md`, `MASTER_PROMPT.md`, `PROJECT_BRIEF.md` or anything in `.claude/`. They belong to a different tool. `BUILD_SPEC.md` is the only spec.
- Do not change any file in `data/`.
- Never run `python -m src.evaluate` without the words "run the test" from me. `python -m src.evaluate --selfcheck-start` is the only exception.
- Never push, log in or upload anything without asking me.
- Never invent a number. Never weaken or skip a test.
- The environment is Windows and PowerShell. Call Python as `.\.venv\Scripts\python.exe` from the repository root.
