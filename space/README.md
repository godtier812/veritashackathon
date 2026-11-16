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

# BRCA Spellcheck

This project is a prototype and is not medical advice.

One app, three pages:

- `/` is the landing page: what the project is, the results and where it fails.
- `/tool/` is the simple tool: type one BRCA1 letter change (C61G, p.Cys61Gly) to see the model's call, the lab's answer and whether they agree.
- `/explorer/` is the same data on the BRCA1 structure, with every lab-tested swap, a position scrubber and the scoreboard. Each tool answer has an "Open <mutation> in 3D" button that opens the explorer on that mutation (`/explorer/?m=C61G`), and the explorer links back to `/tool/?q=C61G`.

BRCA1 is our exam, not the use case. Both pages read cached scores and lab results; no model runs here. Structures (PDB 1JM7 and 1T15) are shown for orientation only; the model saw sequence alone.

Also deploys on Render's free tier from the repository root (`render.yaml`: `pip install -r space/requirements.txt`, then `python space/app.py`, which reads `PORT`).

Space layout: `app.py` (this launcher), `app/` (the tool, unchanged when run on its own with `python app/app.py`), `explorer/` (the static explorer), `landing/` (the static landing page), `requirements.txt`.
