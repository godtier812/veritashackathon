# Credits and Citations

## Data

- Findlay GM et al. (2018). Accurate classification of BRCA1 variants with saturation genome editing. *Nature* 562:217–222. https://doi.org/10.1038/s41586-018-0461-z — lab measurements and function scores (University of Washington).
- University of Washington BRCA1 saturation genome editing resource (classes and cutoffs): https://sge.gs.washington.edu/BRCA1/
- ProteinGym (Notin P et al., 2023), benchmark source of `BRCA1_HUMAN_Findlay_2018.csv`: https://github.com/OATML-Markslab/ProteinGym — also prior work, as it already evaluates ESM-2 on this table. Its published zero-shot leaderboard row for this table (`outputs/reference/proteingym_brca1_leaderboard.csv`, retrieved 2026-11-15 from https://github.com/OATML-Markslab/ProteinGym/tree/main/benchmarks/DMS_zero_shot/substitutions) is used as reference context only.
- AlphaMissense published predictions for BRCA1 (UniProt `P38398`), from the AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk/files/AF-P38398-F1-aa-substitutions.csv (stored unedited as `outputs/reference/alphamissense_P38398.csv`, retrieved 2026-11-15) — scored as a reference row only, after the one test run; not one of our arms.

## Models and methods

- ESM-2 (Lin Z et al., 2023, *Science* 379:1123–1130), model `facebook/esm2_t33_650M_UR50D`: https://huggingface.co/facebook/esm2_t33_650M_UR50D and https://github.com/facebookresearch/esm
- ESM-C 600M (ESM Cambrian; EvolutionaryScale, now Biohub), model `biohub/ESMC-600M`, run through the `esm` Python package: https://huggingface.co/biohub/ESMC-600M and https://github.com/evolutionaryscale/esm — used only for the exploratory `arm4`, added after the one test run.
- Masked-marginal ("hide-a-letter") scoring method: Meier J et al. (2021), "Language models enable zero-shot prediction of the effects of mutations on protein function", NeurIPS. https://doi.org/10.1101/2021.07.09.450648
- BLOSUM62: Henikoff S, Henikoff JG (1992), *PNAS* 89:10915–10919, loaded via Biopython: https://biopython.org
- Kyte J, Doolittle RF (1982), *Journal of Molecular Biology* 157:105–132 — amino acid hydrophobicity scale.
- Zamyatnin AA (1972), *Progress in Biophysics and Molecular Biology* 24:107–123 — amino acid residue volumes in Å³.
- BRCA1-BARD1 RING-domain structure, RCSB Protein Data Bank entry `1JM7`: https://www.rcsb.org/structure/1JM7 — confirms RING zinc-coordinating positions 24, 27, 39, 41, 44, 47, 61, and 64.

## Prior work

- Evo 2 BRCA1 zero-shot notebook (Arc Institute; also published by NVIDIA): https://github.com/ArcInstitute/evo2 and https://docs.nvidia.com/bionemo-recipes/2.7.1/main/examples/bionemo-evo2/zeroshot_brca1/
- AlphaMissense (Cheng J et al., 2023, *Science* 381:eadg7492) — a stronger variant effect predictor, noted as a limitation; we do not run it, and use only its published predictions as the reference row described under Data.
- Follow-up preprint on cell-type context (2025): https://www.medrxiv.org/content/10.1101/2025.08.11.25333423.full.pdf

## Libraries

PyTorch, Hugging Face Transformers, NumPy, pandas, SciPy, scikit-learn, Biopython, Plotly, Matplotlib, Gradio, pytest, and Jupyter (`nbformat`, `nbconvert`, `ipykernel`).

## AI tools, with roles

- **Claude Code (Anthropic):** planning, writing the build specification, and review (running tests, auditing for test-set leakage, checking index math, checking claims in the write-up, reproducibility checks, and mock review); also wrote and ran the post-test reference comparison (`src/reference.py` and its unit tests), which reads ProteinGym's leaderboard row and AlphaMissense's published BRCA1 predictions.
- **Gemini, as a coding agent in Google Antigravity:** writing and executing the project code from the specification.

## The organisers

Veritas AI hackathon project briefs, Hacker Guide, and sample submission.
- AlphaFold Protein Structure Database model of BRCA1 (AF-P38398-F1, v6; Jumper J et al., 2021, *Nature* 596:583–589): https://alphafold.ebi.ac.uk/files/AF-P38398-F1-model_v6.pdb — used only for the burial part of the exploratory `arm6`.
