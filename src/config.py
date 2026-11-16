"""Central configuration holding all constants, filesystem paths, and seed/hashing helpers for BRCA Spellcheck."""

import hashlib
import random
from pathlib import Path

import numpy as np

# 7.1 Regions, classes, targets
SEED: int = 42
PROTEIN_LENGTH: int = 1863
START_REGION: tuple[int, int] = (1, 101)
END_REGION: tuple[int, int] = (1631, 1855)
LOF_BELOW: float = -1.328
FUNC_ABOVE: float = -0.748
HARMFUL_CLASS: str = "LOF"
HARMLESS_CLASS: str = "FUNC"
EXCLUDED_CLASS: str = "INT"
CATCH_TARGET_PERCENT: int = 90
FALSE_ALARM_TARGET_PERCENT: int = 10
N_FOLDS: int = 5
N_BOOTSTRAP: int = 1000
BOOTSTRAP_PERCENTILES: tuple[float, float] = (2.5, 97.5)
N_ROWS: int = 1837
N_POSITIONS: int = 326
N_SCORE_ROWS: int = 6194
EXPECTED_COUNTS: dict[str, dict[str, int]] = {
    "start": {"LOF": 131, "INT": 59, "FUNC": 385},
    "end": {"LOF": 267, "INT": 96, "FUNC": 899},
}
SLIM_CSV_SHA256: str = "4eb3d5597a4a2fcad88fa371b55cee9d8028c6df1ea18b95ad23802484e1a0ac"
RAW_CSV_SHA256: str = "abde5d47e84b1feeac6c32b60c1f3b5d2881accf3ceccfbadccceeea73e81c43"
FASTA_SHA256: str = "47529e88620d4d071f3ebec5212a399d74e102bc48f1ef5a1e1215186a196c36"
SEQUENCE_START: str = "MDLSALRVEEVQNVINAMQK"
SPEARMAN_MIN_START: float = 0.15
BARELY_BEATS_AUROC: float = 0.05
CATCH_DRIFT_POINTS: float = 5.0

# 7.2 Models
MODEL_ID: str = "facebook/esm2_t33_650M_UR50D"
FALLBACK_MODEL_IDS: tuple[str, ...] = (
    "facebook/esm2_t30_150M_UR50D",
    "facebook/esm2_t12_35M_UR50D",
)
ESM_BATCH_SIZE: int = 8

# 7.3 Windows and token-index offsets (1-based inclusive)
WINDOWS: dict[str, tuple[int, int]] = {
    "A": (1, 500),
    "B": (1364, 1863),
}
WINDOW_START: dict[str, int] = {
    "A": 1,
    "B": 1364,
}

# 7.4 Fixed amino-acid alphabet and three-letter mapping
AA_ORDER: str = "ACDEFGHIKLMNPQRSTVWY"
THREE_TO_ONE: dict[str, str] = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

# 7.5 Hydrophobicity (unitless; higher means more hydrophobic)
# Source: Kyte J, Doolittle RF (1982), J. Mol. Biol. 157:105–132.
HYDROPHOBICITY: dict[str, float] = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

# 7.6 Amino-acid residue volume (units: Å³)
# Source: Zamyatnin AA (1972), Prog. Biophys. Mol. Biol. 24:107–123.
VOLUME: dict[str, float] = {
    "A": 88.6,
    "C": 108.5,
    "D": 111.1,
    "E": 138.4,
    "F": 189.9,
    "G": 60.1,
    "H": 153.2,
    "I": 166.7,
    "K": 168.6,
    "L": 166.7,
    "M": 162.9,
    "N": 114.1,
    "P": 112.7,
    "Q": 143.8,
    "R": 173.4,
    "S": 89.0,
    "T": 116.1,
    "V": 140.0,
    "W": 227.8,
    "Y": 193.6,
}

# 7.7 Charge (H is neutral at cellular pH)
POSITIVE_LETTERS: str = "KR"
NEGATIVE_LETTERS: str = "DE"

# 7.8 Arm 3 feature order
FEATURE_ORDER: list[str] = [
    "esm_score",
    "blosum_score",
    "entropy",
    "size_change",
    "hydrophobicity_change",
    "to_proline",
    "from_cysteine",
    "from_glycine",
]

# 7.9 Swap kinds
SWAP_KINDS: list[str] = [
    "to_proline",
    "loss_of_cysteine",
    "charge_change",
    "other",
]

# 7.10 Zinc sites (to confirm on UniProt P38398 before stating as fact)
ZINC_POSITIONS_TO_CONFIRM: list[int] = [24, 27, 39, 41, 44, 47, 61, 64]

# 7.11 Heat map
HEATMAP_ZMIN: float = -12.0
HEATMAP_ZMAX: float = 12.0
HEATMAP_COLORSCALE: str = "RdBu"

# 7.12 Paths
ROOT: Path = Path(__file__).resolve().parent.parent
DATA_DIR: Path = ROOT / "data"
RAW_CSV_PATH: Path = DATA_DIR / "BRCA1_HUMAN_Findlay_2018.csv"
SLIM_CSV_PATH: Path = DATA_DIR / "brca1_missense_scores_slim.csv"
FASTA_PATH: Path = DATA_DIR / "BRCA1_HUMAN_Findlay_2018_wildtype.fasta"
OUTPUTS_DIR: Path = ROOT / "outputs"
ESM_SCORES_PATH: Path = OUTPUTS_DIR / "esm_scores.csv"
FROZEN_PATH: Path = OUTPUTS_DIR / "frozen.json"
RESULTS_DIR: Path = OUTPUTS_DIR / "results"
RUN_LOG_PATH: Path = RESULTS_DIR / "run_log.txt"
SELFCHECK_DIR: Path = OUTPUTS_DIR / "selfcheck_start"
HEATMAP_HTML_PATH: Path = OUTPUTS_DIR / "heatmap.html"
HEATMAP_PNG_PATH: Path = OUTPUTS_DIR / "heatmap.png"
APP_DIR: Path = ROOT / "app"
APP_DATA_PATH: Path = APP_DIR / "app_data.csv"
APP_CONFIG_PATH: Path = APP_DIR / "app_config.json"

# Post-test context (human decision, added after the one test run; not an arm).
# Public third-party files, stored unedited. Nothing here feeds frozen.json.
REFERENCE_DIR: Path = OUTPUTS_DIR / "reference"
ALPHAMISSENSE_PATH: Path = REFERENCE_DIR / "alphamissense_P38398.csv"
ALPHAMISSENSE_SHA256: str = (
    "e9519aa837f538536dd46ccbec755a9089bcc57ba2e33b98032596c36436fa51"
)
ALPHAMISSENSE_URL: str = (
    "https://alphafold.ebi.ac.uk/files/AF-P38398-F1-aa-substitutions.csv"
)
PROTEINGYM_ROW_PATH: Path = REFERENCE_DIR / "proteingym_brca1_leaderboard.csv"
PROTEINGYM_URL: str = (
    "https://github.com/OATML-Markslab/ProteinGym/tree/main/benchmarks/"
    "DMS_zero_shot/substitutions"
)
PROTEINGYM_OUR_MODEL: str = "ESM2 (650M)"
REFERENCE_RETRIEVED: str = "2026-11-15"
REFERENCE_ALPHAMISSENSE_RESULT_PATH: Path = RESULTS_DIR / "reference_alphamissense.json"
REFERENCE_PROTEINGYM_RESULT_PATH: Path = RESULTS_DIR / "reference_proteingym.json"

# Exploratory arm 4 (human decision, added after the one test run). It has its own
# score file, frozen file and run log, so nothing from stages 3 to 5 is touched.
ARM4_MODEL_ID: str = "biohub/ESMC-600M"
ARM4_CANDIDATES: tuple[str, str] = ("esmc", "esm2_esmc_mean")
ARM4_CANDIDATE_NAMES: dict[str, str] = {
    "esmc": "ESM-C 600M zero-shot",
    "esm2_esmc_mean": "the mean of the ESM-2 650M and ESM-C 600M scores",
}
ARM4_SCORES_PATH: Path = OUTPUTS_DIR / "arm4_esmc_scores.csv"
ARM4_FROZEN_PATH: Path = OUTPUTS_DIR / "frozen_arm4.json"
ARM4_RUN_LOG_PATH: Path = RESULTS_DIR / "arm4_run_log.txt"
ARM4_RESULT_PATH: Path = RESULTS_DIR / "arm4.json"
ARM4_PARAGRAPH_PATH: Path = RESULTS_DIR / "arm4_paragraph.md"

# Exploratory arm 5 (human decision, added after the one test run). Same ESM-2 scores as
# arm 2; only the cutoff rule changes, chosen on the start region and frozen before one run.
ARM5_CUTOFF_RULE: str = (
    "highest start-region catch rate minus false-alarm rate; a tie goes to the higher cutoff"
)
ARM5_FROZEN_PATH: Path = OUTPUTS_DIR / "frozen_arm5.json"
ARM5_RUN_LOG_PATH: Path = RESULTS_DIR / "arm5_run_log.txt"
ARM5_RESULT_PATH: Path = RESULTS_DIR / "arm5.json"
ARM5_PARAGRAPH_PATH: Path = RESULTS_DIR / "arm5_paragraph.md"

# 7.13 Fixed text
DISCLAIMER: str = "This project is a prototype and is not medical advice."
DEMO_NOTE: str = (
    "At lab-tested positions a counselor would use the lab result. "
    "This demo shows how the model does where we can check it."
)
VERDICT_HARMFUL: str = "Looks harmful"
VERDICT_CANT_TELL: str = "Can't tell"
VERDICT_HARMLESS: str = "Looks harmless"
ARM_NAMES: dict[str, str] = {
    "arm0": "Always harmless",
    "arm1": "BLOSUM62 rule",
    "arm2": "ESM-2 zero-shot",
    "arm3": "Logistic regression on ESM-2 + BLOSUM62 + chemistry",
}


def set_seeds(seed: int = SEED) -> None:
    """Set random seeds for Python, NumPy, and PyTorch when available."""
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass


def sha256_of_file(path: Path) -> str:
    """Compute and return the 64-character lowercase SHA-256 hex digest of a file."""
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"Cannot hash {path}: file not found.")
    digest = hashlib.sha256()
    with open(file_path, "rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()
