---
title: BRCA Spellcheck Explorer
sdk: static
app_file: index.html
pinned: false
---

# BRCA Spellcheck Explorer

A static, read-only view of the BRCA Spellcheck results: the BRCA1 structure in 3D, every lab-tested swap, the model's call and the lab's answer.

This project is a prototype and is not medical advice.

- The page reads only `data.json`, built by `build_data.py` from the project's result files. Rebuild with `python explorer/build_data.py` from the repository root.
- `check_parity.py` checks that the page's verdict rule gives the same call as the Gradio tool in `app/` for every lab-tested swap.
- Structures: PDB 1JM7 and 1T15, shown for orientation only. The model saw sequence alone.
- 3Dmol.js 2.4.2 is vendored in `vendor/` (BSD-3-Clause). No network access is needed at runtime.
