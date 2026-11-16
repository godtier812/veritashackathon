# First message to the coding agent

For the humans, before you paste:

1. Open this folder as the project in your coding tool. Use the laptop with the NVIDIA GPU if you can, so stage 3 does not need a handoff.
2. Fill in the one bracket in the message below with the tool and model you are using.
3. Paste everything below the line as your first message.

Words you will need later:

- At checkpoint B, after you have read `outputs/frozen.json` and the leak audit says CLEAN: `run the test`
- At checkpoint C, after you have read the results: `go C`
- If the agent asks about a GPU route, a fallback model, a push or an upload, answer in plain words. It will not do any of those without asking.

---

You are the coding agent for a two-person hackathon team. The project is called BRCA Spellcheck. The deadline is tonight, Sunday Nov 15 2026, 10:59 PM US Central.

You are running as: [TOOL AND MODEL, for example "Claude Code, Claude Opus 5.5"]. Use exactly this text when you credit the coding agent (`BUILD_SPEC.md` section 12.6).

The complete specification is `BUILD_SPEC.md` in the repository root. It is long (about 2,900 lines) and it stands alone. Standing rules are in `AGENTS.md`.

Do these four things, in this order, and nothing else yet:

1. **Read `BUILD_SPEC.md` in full, top to bottom.** Do not skim, and do not stop at the table of contents. Then read `AGENTS.md`.

2. **Restate the locked rules in your own words.** Cover all fifteen rules in section 3, plus:
   - what counts as "peeking" at end-region labels, and which three things you may do with the end region before the test run;
   - exactly how a cutoff is chosen, using the six-row worked example in section 8.3;
   - how a protein position becomes a token index, using position 1775 as the example;
   - the checkpoints: which one is automatic, which two are full stops, and the exact words you must wait for at each stop;
   - how the six reviews of section 11 will be run in the tool you are;
   - the list of things you must never claim.

3. **List every ambiguity you see:** anything in the spec that could be read two ways, anything that looks contradictory, and anything you would otherwise have to guess. For each one, say which reading you would choose and why. If you see none, say so. Do not resolve an ambiguity about a locked rule by yourself: ask.

4. **Then start stage 1 only** (section 10). Before writing any file, print the stage goal, the files you will create, and the current time. When stage 1's acceptance criteria all pass, commit with the stage's commit message, give the end-of-stage report, and stop for me to read it. After stage 1 you may run stages 2, 3 and 4 back to back, stopping only where the spec says to stop.

Things to remember from the first minute:

- `BUILD_SPEC.md` is the only spec. Do not read `MASTER_PROMPT.md` or `PROJECT_BRIEF.md`: they are superseded. Never edit or commit `CLAUDE.md` or anything in `.claude/`.
- If you are Claude Code, the review agents in `.claude/agents/` are yours to run as subagents (section 11).
- Do not change any file in `data/`.
- Never run `python -m src.evaluate` without the words "run the test" from me. `python -m src.evaluate --selfcheck-start` is the only exception.
- Never push, log in or upload anything without asking me.
- Never invent a number. Never weaken or skip a test.
- Before your first commit, check `git config user.name` and `git config user.email` print the neutral team identity (section 0, point 13). No personal names in any file or commit.
- The environment is Windows and PowerShell. Call Python as `.\.venv\Scripts\python.exe` from the repository root.
