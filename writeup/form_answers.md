# Submission Form Answers

| Form field | Answer |
| --- | --- |
| Project title | FoldFlux: an AI second opinion on uncertain BRCA1 mutations |
| Building as | Team |
| Team name | [the humans fill in] |
| Team members | [names, emails, schools — the humans fill this in] |
| Track | AI + Health |
| Starting point | Our own idea |
| One-sentence pitch | FoldFlux helps a genetic counselor judge an "uncertain" mutation by taking the mutation and giving back whether a protein-reading AI thinks it looks harmful, with that AI graded in the open against CRISPR lab results for BRCA1. |
| Code link | https://github.com/godtier812/veritashackathon |
| Write-up link | https://github.com/godtier812/veritashackathon#readme |
| Demo video link | [the humans fill in, or left blank] |
| What we added | Our own idea. Grading a language model on this lab data has been done before (ProteinGym; the Evo 2 BRCA1 notebook), and stronger tools such as AlphaMissense exist. We added four things: holding out a whole protein region with cutoffs frozen in advance; testing whether start-region labels help at the other end; counting misses and false alarms separately; and explaining a failure case. |
| Datasets, models and AI tools used | Data: Findlay et al. (2018) UW BRCA1 SGE dataset via ProteinGym (`BRCA1_HUMAN_Findlay_2018.csv`) and PDB `1JM7`; Models & methods: ESM-2 (`facebook/esm2_t33_650M_UR50D`, Lin et al. 2023) with masked-marginal scoring (Meier et al. 2021), BLOSUM62 (Henikoff & Henikoff 1992) via Biopython, Kyte-Doolittle (1982) hydropathy, Zamyatnin (1972) residue volumes, and scikit-learn logistic regression; Libraries: PyTorch, Hugging Face Transformers, NumPy, pandas, SciPy, scikit-learn, Biopython, Plotly, Matplotlib, Gradio, pytest, and Jupyter (`nbformat`, `nbconvert`, `ipykernel`); AI tools: Claude Code (Anthropic) for planning, writing the build specification, and reviewing tests, leakage, index maths, claims, reproducibility, and mock judging, plus the post-test reference comparison (ProteinGym leaderboard row and AlphaMissense published BRCA1 predictions from the AlphaFold Database, reference only), and Gemini in Google Antigravity for writing and executing the project code from the specification. |
| New to AI? | [the humans answer] |

## Copyable Text Blocks

### Project title
```text
FoldFlux: an AI second opinion on uncertain BRCA1 mutations
```

### One-sentence pitch
```text
FoldFlux helps a genetic counselor judge an "uncertain" mutation by taking the mutation and giving back whether a protein-reading AI thinks it looks harmful, with that AI graded in the open against CRISPR lab results for BRCA1.
```

### Code link
```text
https://github.com/godtier812/veritashackathon
```

### Write-up link
```text
https://github.com/godtier812/veritashackathon#readme
```

### What we added
```text
Our own idea. Grading a language model on this lab data has been done before (ProteinGym; the Evo 2 BRCA1 notebook), and stronger tools such as AlphaMissense exist. We added four things: holding out a whole protein region with cutoffs frozen in advance; testing whether start-region labels help at the other end; counting misses and false alarms separately; and explaining a failure case.
```

### Datasets, models and AI tools used
```text
Data: Findlay et al. (2018) UW BRCA1 SGE dataset via ProteinGym (BRCA1_HUMAN_Findlay_2018.csv) and PDB 1JM7; Models & methods: ESM-2 (facebook/esm2_t33_650M_UR50D, Lin et al. 2023) with masked-marginal scoring (Meier et al. 2021), BLOSUM62 (Henikoff & Henikoff 1992) via Biopython, Kyte-Doolittle (1982) hydropathy, Zamyatnin (1972) residue volumes, and scikit-learn logistic regression; Libraries: PyTorch, Hugging Face Transformers, NumPy, pandas, SciPy, scikit-learn, Biopython, Plotly, Matplotlib, Gradio, pytest, and Jupyter (nbformat, nbconvert, ipykernel); AI tools: Claude Code (Anthropic) for planning, writing the build specification, and reviewing tests, leakage, index maths, claims, reproducibility, and mock judging, plus the post-test reference comparison (ProteinGym leaderboard row and AlphaMissense published BRCA1 predictions from the AlphaFold Database, reference only), and Gemini in Google Antigravity for writing and executing the project code from the specification.
```
