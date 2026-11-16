---
title: FoldFlux
emoji: 🧬
colorFrom: red
colorTo: blue
sdk: gradio
sdk_version: 6.32.0
app_file: app.py
pinned: false
---

# FoldFlux

FoldFlux is a zero-shot protein language model prototype (`facebook/esm2_t33_650M_UR50D`) that scores single-letter BRCA1 amino-acid substitutions from sequence alone and compares those predictions against the Findlay et al. (2018) saturation genome editing lab assay. BRCA1 is the exam, not the use case: because BRCA1 is one of the few human genes with a dense experimental answer key across hundreds of positions, we use it to test how far a sequence-only model can be trusted—and where it misses harmful variants or raises false alarms—before anyone relies on such a model for an under-studied gene with no lab table.

*This project is a prototype and is not medical advice.*
