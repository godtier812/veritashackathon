"""Self-contained Gradio demo for BRCA Spellcheck (arm 2 zero-shot ESM-2 lookup)."""

import json
from pathlib import Path
import re

import gradio as gr
import pandas as pd

# Copies of constants from src/config.py (app/app.py cannot import src)
PROTEIN_LENGTH: int = 1863  # Copy of config.PROTEIN_LENGTH
START_REGION: tuple[int, int] = (1, 101)  # Copy of config.START_REGION
END_REGION: tuple[int, int] = (1631, 1855)  # Copy of config.END_REGION
AA_ORDER: str = "ACDEFGHIKLMNPQRSTVWY"  # Copy of config.AA_ORDER
THREE_TO_ONE: dict[str, str] = {  # Copy of config.THREE_TO_ONE
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
HARMFUL_CLASS: str = "LOF"  # Copy of config.HARMFUL_CLASS
HARMLESS_CLASS: str = "FUNC"  # Copy of config.HARMLESS_CLASS
EXCLUDED_CLASS: str = "INT"  # Copy of config.EXCLUDED_CLASS
DISCLAIMER: str = (
    "This project is a prototype and is not medical advice."
)  # Copy of config.DISCLAIMER
DEMO_NOTE: str = (  # Copy of config.DEMO_NOTE
    "At lab-tested positions a counselor would use the lab result. "
    "This demo shows how the model does where we can check it."
)
VERDICT_HARMFUL: str = "Looks harmful"  # Copy of config.VERDICT_HARMFUL
VERDICT_CANT_TELL: str = "Can't tell"  # Copy of config.VERDICT_CANT_TELL
VERDICT_HARMLESS: str = "Looks harmless"  # Copy of config.VERDICT_HARMLESS

_PATTERN_ONE = re.compile(r"^([A-Za-z])\s*(\d{1,5})\s*([A-Za-z])$", re.IGNORECASE)
_PATTERN_THREE = re.compile(
    r"^([A-Za-z]{3})\s*(\d{1,5})\s*([A-Za-z]{3})$", re.IGNORECASE
)
_STRIP_MARKDOWN_RE = re.compile(r"[`<>*_\[\]]")

_R1_MESSAGE: str = (
    "Type a mutation such as C61G: the normal letter, the position, then the new letter."
)
_R7_MESSAGE: str = "Something went wrong reading that input. Try a mutation such as C61G."
_R8_MESSAGE: str = (
    "End-region answers (positions 1631 to 1855) are hidden until the one test run is done."
)

APP_DIR: Path = Path(__file__).resolve().parent
APP_DATA_PATH: Path = APP_DIR / "app_data.csv"
APP_CONFIG_PATH: Path = APP_DIR / "app_config.json"

with open(APP_CONFIG_PATH, "r", encoding="utf-8") as _cfg_handle:
    APP_CONFIG: dict = json.load(_cfg_handle)

APP_DATA: pd.DataFrame = pd.read_csv(
    APP_DATA_PATH,
    dtype={
        "position": "int64",
        "wild_type": str,
        "new_letter": str,
        "esm_score": "float64",
        "lab_class": str,
        "lab_score": str,
    },
    keep_default_na=False,
)

_WT_BY_POS: dict[int, str] = {}
_ROW_BY_POS_NEW: dict[tuple[int, str], dict] = {}
_END_LAB_ROWS: list[tuple[float, str]] = []
for _pos, _wt, _new, _esm, _lclass, _lscore in zip(
    APP_DATA["position"],
    APP_DATA["wild_type"],
    APP_DATA["new_letter"],
    APP_DATA["esm_score"],
    APP_DATA["lab_class"],
    APP_DATA["lab_score"],
):
    _pos_int = int(_pos)
    _esm_float = float(_esm)
    _lclass_str = str(_lclass)
    _WT_BY_POS[_pos_int] = str(_wt)
    _ROW_BY_POS_NEW[(_pos_int, str(_new))] = {
        "esm_score": _esm_float,
        "lab_class": _lclass_str,
        "lab_score": None if _lscore == "" else float(_lscore),
    }
    if END_REGION[0] <= _pos_int <= END_REGION[1] and _lclass_str != "":
        _END_LAB_ROWS.append((_esm_float, _lclass_str))


def _r2_message(raw_text: str) -> str:
    """Format Refusal R2 with 40-char truncation and Markdown characters stripped."""
    cleaned = _STRIP_MARKDOWN_RE.sub("", raw_text[:40])
    return (
        f'I could not read "{cleaned}". Use the form C61G: the normal letter, '
        "the position, then the new letter. Three-letter names such as Cys61Gly also work."
    )


def band_verdict(
    harm_score: float,
    cutoff_a: float,
    cutoff_b: float | None,
) -> str:
    """Map a single harm score and optional band cutoffs to the fixed verdict text."""
    score = float(harm_score)
    cut_a = float(cutoff_a)
    cut_b = None if cutoff_b is None else float(cutoff_b)
    if cut_b is not None and cut_b > cut_a:
        if score >= cut_b:
            return VERDICT_HARMFUL
        if score < cut_a:
            return VERDICT_HARMLESS
        return VERDICT_CANT_TELL
    if score >= cut_a:
        return VERDICT_HARMFUL
    return VERDICT_HARMLESS


def parse_input(text: str) -> tuple[str, int, str]:
    """Parse one-letter or three-letter mutation input or raise ValueError with R1-R4."""
    if text is None:
        raise ValueError(_R1_MESSAGE)
    if not isinstance(text, str):
        raise TypeError("Input must be a string.")
    stripped = text.strip()
    if stripped == "":
        raise ValueError(_R1_MESSAGE)

    if stripped[:2].lower() == "p.":
        body = stripped[2:].strip()
    else:
        body = stripped

    m1 = _PATTERN_ONE.fullmatch(body)
    if m1 is not None:
        wt_raw, pos_str, new_raw = m1.groups()
        wt = wt_raw.upper()
        new_letter = new_raw.upper()
        if wt not in AA_ORDER or new_letter not in AA_ORDER:
            raise ValueError(_r2_message(text))
    else:
        m2 = _PATTERN_THREE.fullmatch(body)
        if m2 is None:
            raise ValueError(_r2_message(text))
        wt_raw, pos_str, new_raw = m2.groups()
        wt_mapped = THREE_TO_ONE.get(wt_raw.upper())
        new_mapped = THREE_TO_ONE.get(new_raw.upper())
        if wt_mapped is None or new_mapped is None:
            raise ValueError(_r2_message(text))
        wt = wt_mapped
        new_letter = new_mapped

    position = int(pos_str)
    if wt == new_letter:
        raise ValueError(
            f'"{wt}{position}{new_letter}" is not a change: both letters are {wt}.'
        )
    protein_len = int(APP_CONFIG.get("protein_length", PROTEIN_LENGTH))
    if position < 1 or position > protein_len:
        raise ValueError(
            f"Position {position} is outside BRCA1, which has 1,863 positions."
        )
    return wt, position, new_letter


def lookup(wild_type: str, position: int, new_letter: str) -> dict:
    """Look up a parsed mutation in app_data.csv or raise ValueError with R5, R6, or R8."""
    if position not in _WT_BY_POS:
        raise ValueError(
            f"Position {position} is not one of the 326 lab-tested positions "
            "(1 to 101 and 1631 to 1855). This demo only answers where its answers "
            "can be checked against the lab, so it will not guess here."
        )

    in_end_region = END_REGION[0] <= position <= END_REGION[1]
    if in_end_region and not bool(APP_CONFIG.get("end_labels_included", False)):
        raise ValueError(_R8_MESSAGE)

    actual = _WT_BY_POS[position]
    if wild_type != actual:
        if actual == new_letter:
            raise ValueError(
                f"BRCA1 has {actual} at position {position}, not {wild_type}."
            )
        raise ValueError(
            f"BRCA1 has {actual} at position {position}, not {wild_type}. "
            f"Did you mean {actual}{position}{new_letter}?"
        )

    row = _ROW_BY_POS_NEW[(position, new_letter)]
    esm_score = float(row["esm_score"])
    harm = 0.0 - esm_score

    cutoff_a = float(APP_CONFIG["cutoff_a"])
    raw_cutoff_b = APP_CONFIG.get("cutoff_b")
    cutoff_b = None if raw_cutoff_b is None else float(raw_cutoff_b)
    if not bool(APP_CONFIG.get("has_band", cutoff_b is not None and cutoff_b > cutoff_a)):
        cutoff_b = None

    verdict = band_verdict(harm, cutoff_a, cutoff_b)
    region = "start" if START_REGION[0] <= position <= START_REGION[1] else "end"

    return {
        "mutant": f"{wild_type}{position}{new_letter}",
        "esm_score": esm_score,
        "verdict": verdict,
        "lab_class": str(row["lab_class"]),
        "lab_score": row["lab_score"],
        "region": region,
    }


def _append_notes(content: str) -> str:
    """Append a horizontal rule and both fixed notes to any output or refusal."""
    return f"{content}\n\n---\n{DEMO_NOTE}\n{DISCLAIMER}"


def _end_region_verdict_line(
    verdict: str,
    cutoff_a: float,
    cutoff_b: float | None,
) -> str:
    """Format the held-out end-region track-record line for the given verdict."""
    n_lof = 0
    n_func = 0
    n_int = 0
    for esm_val, lclass in _END_LAB_ROWS:
        if band_verdict(0.0 - esm_val, cutoff_a, cutoff_b) == verdict:
            if lclass == HARMFUL_CLASS:
                n_lof += 1
            elif lclass == HARMLESS_CLASS:
                n_func += 1
            elif lclass == EXCLUDED_CLASS:
                n_int += 1

    if verdict == VERDICT_HARMFUL:
        clear_total = n_lof + n_func
        pct = 100.0 * n_lof / clear_total if clear_total > 0 else 0.0
        return (
            '**End-region track record:** On the held-out end region, "Looks harmful" '
            f"matched the lab (LOF, harmful) on {n_lof} of {clear_total} clear-class "
            f"mutations ({pct:.1f}%; {n_func} were FUNC harmless, plus {n_int} INT in between)."
        )
    if verdict == VERDICT_HARMLESS:
        clear_total = n_lof + n_func
        pct = 100.0 * n_func / clear_total if clear_total > 0 else 0.0
        return (
            '**End-region track record:** On the held-out end region, "Looks harmless" '
            f"matched the lab (FUNC, harmless) on {n_func} of {clear_total} clear-class "
            f"mutations ({pct:.1f}%; {n_lof} were LOF harmful, plus {n_int} INT in between)."
        )
    return (
        '**End-region track record:** On the held-out end region, mutations given '
        f'"Can\'t tell" were {n_lof} LOF harmful, {n_func} FUNC harmless, and {n_int} INT '
        f"in between in the lab ({n_lof + n_func + n_int} total)."
    )


def check_mutation(text: str) -> str:
    """Validate, look up, and format the Markdown response for a single mutation input."""
    try:
        wild_type, position, new_letter = parse_input(text)
        info = lookup(wild_type, position, new_letter)

        cutoff_a = float(APP_CONFIG["cutoff_a"])
        raw_cutoff_b = APP_CONFIG.get("cutoff_b")
        cutoff_b = None if raw_cutoff_b is None else float(raw_cutoff_b)
        has_band_flag = (
            bool(APP_CONFIG.get("has_band", False))
            and cutoff_b is not None
            and cutoff_b > cutoff_a
        )

        if has_band_flag and cutoff_b is not None:
            how_text = (
                f"Scores at or below {-cutoff_b:.2f} look harmful; "
                f"scores above {-cutoff_a:.2f} look harmless; "
                "in between the model can't tell. Both lines were set on the start region only."
            )
        else:
            how_text = (
                f"Scores at or below {-cutoff_a:.2f} look harmful; "
                "anything above looks harmless. The line was set on the start region only. "
                'There is no "can\'t tell" band, because the two cutoffs did not leave a gap.'
            )

        lab_class = info["lab_class"]
        lab_score = info["lab_score"]
        if lab_class == HARMFUL_CLASS and lab_score is not None:
            lab_text = f"Loss of function (harmful). Lab function score {lab_score:.3f}."
        elif lab_class == HARMLESS_CLASS and lab_score is not None:
            lab_text = f"Functional (harmless). Lab function score {lab_score:.3f}."
        elif lab_class == EXCLUDED_CLASS and lab_score is not None:
            lab_text = f"Intermediate (in between). Lab function score {lab_score:.3f}."
        else:
            lab_text = "not lab-tested"

        lines = [
            f"## {info['mutant']}",
            f"**Model verdict:** {info['verdict']}",
            (
                f"**Model score:** {info['esm_score']:.2f} "
                "(more negative means the model finds this swap more surprising)"
            ),
            f"**How the verdict is set:** {how_text}",
            f"**Lab result (answer key):** {lab_text}",
        ]
        if bool(APP_CONFIG.get("end_labels_included", False)):
            lines.append(
                _end_region_verdict_line(
                    info["verdict"],
                    cutoff_a,
                    cutoff_b if has_band_flag else None,
                )
            )
        return _append_notes("\n\n".join(lines))
    except ValueError as exc:
        return _append_notes(str(exc))
    except Exception:
        return _append_notes(_R7_MESSAGE)


with gr.Blocks(title="BRCA Spellcheck") as demo:
    gr.Markdown(
        "# BRCA Spellcheck\n\n"
        "Look up how a protein language model (`ESM-2`) scores a single-letter change "
        "in BRCA1 from sequence alone, and compare its call against the lab answer key.\n\n"
        "BRCA1 is our exam, not the use case: it is one of the few genes with a lab "
        "answer key, so we can check the model here."
    )
    mutation_input = gr.Textbox(label="Mutation", placeholder="C61G")
    check_button = gr.Button("Check")
    result_output = gr.Markdown()
    gr.Examples(
        examples=["C61G", "D2E", "C61A", "p.Cys64Gly"],
        inputs=mutation_input,
    )
    gr.Markdown(f"{DEMO_NOTE}\n\n{DISCLAIMER}")

    check_button.click(fn=check_mutation, inputs=mutation_input, outputs=result_output)
    mutation_input.submit(fn=check_mutation, inputs=mutation_input, outputs=result_output)


if __name__ == "__main__":
    demo.launch()
