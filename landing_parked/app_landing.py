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


def track_record_text(verdict: str, n_lof: int, n_func: int, n_int: int) -> str:
    total = n_lof + n_func + n_int
    if verdict == VERDICT_HARMFUL:
        return (
            f'On the held-out end region the model said "Looks harmful" {total} times. '
            f"The lab agreed {n_lof} times, {n_func} were false alarms (the protein still worked), "
            f"and {n_int} landed in the lab's grey zone, where it could not call them either way."
        )
    if verdict == VERDICT_HARMLESS:
        return (
            f'On the held-out end region the model said "Looks harmless" {total} times. '
            f"The lab agreed {n_func} times, {n_lof} were misses (the protein was broken), "
            f"and {n_int} landed in the lab's grey zone, where it could not call them either way."
        )
    return (
        f'On the held-out end region the model said "Can\'t tell" {total} times. '
        f"The lab found {n_lof} of those broken, {n_func} working, and {n_int} in its grey zone."
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


lof_below: float = -1.328
func_above: float = -0.748
na: str = "not available"
results_dir: Path = APP_DIR / "results"
repo_results_dir: Path = APP_DIR.parent / "outputs" / "results"
bundle_names: list[str] = [
    "headline.csv",
    "band_counts.csv",
    "failures.csv",
    "auroc.csv",
    "bootstrap.json",
    "result_sentences.txt",
]
aa_names: dict[str, str] = {
    "A": "alanine",
    "C": "cysteine",
    "D": "aspartate",
    "E": "glutamate",
    "F": "phenylalanine",
    "G": "glycine",
    "H": "histidine",
    "I": "isoleucine",
    "K": "lysine",
    "L": "leucine",
    "M": "methionine",
    "N": "asparagine",
    "P": "proline",
    "Q": "glutamine",
    "R": "arginine",
    "S": "serine",
    "T": "threonine",
    "V": "valine",
    "W": "tryptophan",
    "Y": "tyrosine",
}


def find_file(name: str) -> Path | None:
    for folder in (repo_results_dir, results_dir):
        path = folder / name
        if path.exists():
            return path
    return None


def heatmap_path() -> Path | None:
    for path in (APP_DIR.parent / "outputs" / "heatmap.png", results_dir / "heatmap.png"):
        if path.exists():
            return path
    return None


def check_bundle() -> list[str]:
    differ = []
    pairs = [(repo_results_dir / n, results_dir / n) for n in bundle_names]
    pairs.append((APP_DIR.parent / "outputs" / "heatmap.png", results_dir / "heatmap.png"))
    for repo_path, app_path in pairs:
        if repo_path.exists() and app_path.exists():
            if repo_path.read_bytes() != app_path.read_bytes():
                differ.append(app_path.name)
    if differ:
        print("warning: app/results copies differ from outputs:", ", ".join(differ))
    return differ


def read_table(name: str) -> pd.DataFrame | None:
    path = find_file(name)
    if path is None:
        return None
    try:
        return pd.read_csv(path, keep_default_na=False)
    except Exception:
        return None


def read_json(name: str) -> dict | None:
    path = find_file(name)
    if path is None:
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def read_sentences() -> dict[str, str]:
    path = find_file("result_sentences.txt")
    found = {}
    if path is None:
        return found
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" in line:
            key, text = line.split(":", 1)
            found[key.strip()] = text.strip()
    return found


def read_blosum() -> dict[tuple[str, str], int]:
    table = read_table("blosum62.csv")
    found = {}
    if table is None:
        return found
    for _, row in table.iterrows():
        for letter in AA_ORDER:
            found[(str(row["letter"]), letter)] = int(row[letter])
    return found


bundle_differs = check_bundle()
headline = read_table("headline.csv")
band_counts = read_table("band_counts.csv")
failures = read_table("failures.csv")
auroc = read_table("auroc.csv")
bootstrap = read_json("bootstrap.json")
sentences = read_sentences()
blosum = read_blosum()

esm_values = APP_DATA["esm_score"].astype(float)
esm_lo = float(int(esm_values.min()) - 1)
esm_hi = float(int(esm_values.max()) + 1)
lab_rows = APP_DATA[APP_DATA["lab_class"] != ""]
lab_values = lab_rows["lab_score"].astype(float)
lab_lo = float(int(lab_values.min()) - 1) if len(lab_values) else -4.0
lab_hi = float(int(lab_values.max()) + 1) if len(lab_values) else 2.0
lab_tested_list = [
    f"{wt}{int(pos)}{new}"
    for pos, wt, new in zip(lab_rows["position"], lab_rows["wild_type"], lab_rows["new_letter"])
]
start_lab = lab_rows[lab_rows["position"] <= START_REGION[1]]
start_n_harmful = int((start_lab["lab_class"] == HARMFUL_CLASS).sum())
start_n_harmless = int((start_lab["lab_class"] == HARMLESS_CLASS).sum())
scored_positions = sorted(_WT_BY_POS)


def cutoffs() -> tuple[float, float | None]:
    cut_a = float(APP_CONFIG["cutoff_a"])
    raw_b = APP_CONFIG.get("cutoff_b")
    cut_b = None if raw_b is None else float(raw_b)
    if not bool(APP_CONFIG.get("has_band", False)) or cut_b is None or cut_b <= cut_a:
        cut_b = None
    return cut_a, cut_b


def esc(text) -> str:
    out = str(text)
    for old, new in (("&", "&amp;"), ("<", "&lt;"), (">", "&gt;"), ('"', "&quot;"), ("'", "&#x27;")):
        out = out.replace(old, new)
    return out


def signed(value: float, places: int) -> str:
    text = f"{value:.{places}f}"
    if text.startswith("-"):
        return "−" + text[1:]
    if float(text) > 0:
        return "+" + text
    return text


def plain_signed(value: float, places: int) -> str:
    return f"{value:.{places}f}"


def numspan(value: float, places: int, show_plus: bool = True) -> str:
    text = signed(value, places) if show_plus else signed(value, places).lstrip("+")
    return f'<span class="num" data-count="{value}" data-dec="{places}">{text}</span>'


icon_paths: dict[str, str] = {
    "harmful": '<path d="M8 2.5 14.2 13.5H1.8Z"/>',
    "harmless": '<circle cx="8" cy="8" r="5.6"/>',
    "cant": '<path d="M8 1.8 14.2 8 8 14.2 1.8 8Z"/>',
    "untested": '<circle cx="8" cy="8" r="5.6" stroke-dasharray="2.2 2.2"/>',
    "agree": '<path d="M2.8 8.6 6.4 12.2 13.2 4.4"/>',
    "disagree": '<path d="M3.5 3.5 12.5 12.5M12.5 3.5 3.5 12.5"/>',
    "grey": '<path d="M2.5 8h11"/>',
    "copy": '<rect x="5" y="5" width="8.5" height="8.5" rx="1.5"/><path d="M11 5V3.5A1.5 1.5 0 0 0 9.5 2h-6A1.5 1.5 0 0 0 2 3.5v6A1.5 1.5 0 0 0 3.5 11H5"/>',
    "chev": '<path d="M6 3.5 10.5 8 6 12.5"/>',
}


def icon(kind: str, extra: str = "") -> str:
    return (
        f'<svg class="ico {extra}" viewBox="0 0 16 16" aria-hidden="true">'
        f"{icon_paths[kind]}</svg>"
    )


verdict_style: dict[str, tuple[str, str]] = {
    VERDICT_HARMFUL: ("harmful", "st-harmful"),
    VERDICT_HARMLESS: ("harmless", "st-harmless"),
    VERDICT_CANT_TELL: ("cant", "st-cant"),
}

verdict_plain: dict[str, str] = {
    VERDICT_HARMFUL: "The model thinks this change probably breaks the protein.",
    VERDICT_HARMLESS: "The model thinks the protein probably copes with this change.",
    VERDICT_CANT_TELL: "The model is not sure about this one, so it does not make a call.",
}


def agreement(verdict: str, lab_class: str) -> dict:
    if lab_class == "":
        return {
            "key": "untested",
            "label": "Not lab-tested",
            "icon": "untested",
            "css": "st-untested",
            "mean": "Nobody has measured this exact swap in the lab, so there is no answer key to check the model against.",
        }
    if verdict == VERDICT_CANT_TELL:
        return {
            "key": "cant",
            "label": "Model can't tell",
            "icon": "cant",
            "css": "st-cant",
            "mean": "The score landed in the middle zone, where the model is wrong too often to make a call. Here the lab result is the answer.",
        }
    if lab_class == EXCLUDED_CLASS:
        return {
            "key": "grey",
            "label": "Lab result in between",
            "icon": "grey",
            "css": "st-untested",
            "mean": "The lab measured something between broken and working, so the model's call can't be marked right or wrong.",
        }
    if (verdict == VERDICT_HARMFUL and lab_class == HARMFUL_CLASS) or (
        verdict == VERDICT_HARMLESS and lab_class == HARMLESS_CLASS
    ):
        return {
            "key": "agree",
            "label": "Model and lab agree",
            "icon": "agree",
            "css": "st-agree",
            "mean": "The model's call and the lab's answer point the same way. That is one match, not a promise about the next one.",
        }
    if verdict == VERDICT_HARMLESS and lab_class == HARMFUL_CLASS:
        return {
            "key": "miss",
            "label": "Model missed this one",
            "icon": "disagree",
            "css": "st-disagree",
            "mean": "The lab found this change breaks the protein, but the model said it looked fine. This is the kind of mistake that matters most.",
        }
    return {
        "key": "alarm",
        "label": "Model raised a false alarm",
        "icon": "disagree",
        "css": "st-disagree",
        "mean": "The model flagged this change, but in the lab the protein kept working normally.",
    }


def region_words(region: str) -> str:
    if region == "start":
        return "the RING region at the start of BRCA1"
    return "the BRCT region near the end of BRCA1"


def model_sentence(esm: float, verdict: str) -> str:
    cut_a, cut_b = cutoffs()
    line_a = -cut_a
    line_b = -cut_b if cut_b is not None else line_a
    score = numspan(esm, 2)
    if verdict == VERDICT_HARMFUL:
        gap = abs(esm - line_b)
        return (
            f"The model scored it {score}, which is <span class=\"num\">{gap:.2f}</span> past the harmful line "
            f"at <span class=\"num\">{signed(line_b, 2)}</span>, so its call is harmful."
        )
    if verdict == VERDICT_HARMLESS:
        gap = abs(esm - line_a)
        return (
            f"The model scored it {score}, <span class=\"num\">{gap:.2f}</span> above the harmless line "
            f"at <span class=\"num\">{signed(line_a, 2)}</span>, so its call is harmless."
        )
    return (
        f"The model scored it {score}, between the harmless line "
        f"(<span class=\"num\">{signed(line_a, 2)}</span>) and the harmful line "
        f"(<span class=\"num\">{signed(line_b, 2)}</span>), so it makes no call."
    )


def lab_sentence(lab_class: str, lab_score: float | None, key: str) -> str:
    if lab_class == "" or lab_score is None:
        return "Nobody tested this exact swap in the lab, so we cannot check the call."
    score = f'<span class="num">{signed(lab_score, 3)}</span>'
    tails = {
        "agree": ", and the model agreed.",
        "miss": ", so the model missed it.",
        "alarm": ", so the model raised a false alarm.",
        "cant": "; the model made no call.",
        "grey": ", so there is no clear answer to grade the model against.",
    }
    if lab_class == HARMFUL_CLASS:
        start = f"In the lab this change knocked out the protein's function (score {score})"
    elif lab_class == HARMLESS_CLASS:
        start = f"In the lab the protein kept working (score {score})"
    else:
        start = f"The lab's score ({score}) fell in its grey zone between broken and working"
    return start + tails.get(key, ".")


def pct(value: float, low: float, high: float) -> float:
    return max(0.0, min(100.0, (value - low) / (high - low) * 100.0))


def track(title: str, value_text: str, low: float, high: float, left_cut: float, right_cut: float,
          zone_fills: list[str], zone_names: list[tuple[str, str, str]], cut_labels: list[str],
          value: float, marker_class: str) -> str:
    p_left = pct(left_cut, low, high)
    p_right = pct(right_cut, low, high)
    p_val = pct(value, low, high)
    rects = (
        f'<rect x="0" y="0" width="{p_left:.3f}" height="12" style="fill:{zone_fills[0]}"/>'
        f'<rect x="{p_left:.3f}" y="0" width="{p_right - p_left:.3f}" height="12" style="fill:{zone_fills[1]}"/>'
        f'<rect x="{p_right:.3f}" y="0" width="{100 - p_right:.3f}" height="12" style="fill:{zone_fills[2]}"/>'
    )
    names = (
        f'<span class="{zone_names[0][1]}">{icon(zone_names[0][2])}{zone_names[0][0]}</span>'
        f'<span class="{zone_names[1][1]}">{icon(zone_names[1][2])}{zone_names[1][0]}</span>'
        f'<span class="{zone_names[2][1]}">{icon(zone_names[2][2])}{zone_names[2][0]}</span>'
    )
    return (
        f'<div class="track"><div class="track-label"><span>{title}</span><b class="mono">{value_text}</b></div>'
        f'<div class="bar"><svg class="zones" viewBox="0 0 100 12" preserveAspectRatio="none" aria-hidden="true">{rects}</svg>'
        f'<div class="cut cut-l" style="left:{p_left:.2f}%"><span>{cut_labels[0]}</span></div>'
        f'<div class="cut cut-r" style="left:{p_right:.2f}%"><span>{cut_labels[1]}</span></div>'
        f'<div class="mk-layer" style="--p:{p_val:.2f}%"><div class="mk {marker_class}"></div></div></div>'
        f'<div class="zone-names">{names}</div></div>'
    )


def gauge(esm: float, lab_class: str, lab_score: float | None, mutant: str) -> str:
    cut_a, cut_b = cutoffs()
    line_a = -cut_a
    line_b = -cut_b if cut_b is not None else line_a
    model_track = track(
        "Model's call (ESM-2 score)",
        signed(esm, 2),
        esm_lo,
        esm_hi,
        line_b,
        line_a,
        ["var(--harmful-soft)", "var(--cant-soft)", "var(--harmless-soft)"],
        [("harmful", "st-harmful", "harmful"), ("can't tell", "st-cant", "cant"), ("harmless", "st-harmless", "harmless")],
        [signed(line_b, 2), signed(line_a, 2)],
        esm,
        "",
    )
    if lab_class == "" or lab_score is None:
        lab_track = (
            '<div class="track"><div class="track-label"><span>Lab\'s answer (function score)</span>'
            '<b class="mono">not tested</b></div>'
            f'<div class="untested-track">{icon("untested", "st-untested")}&nbsp; No lab measurement for this swap.</div></div>'
        )
    else:
        lab_track = track(
            "Lab's answer (function score)",
            signed(lab_score, 3),
            lab_lo,
            lab_hi,
            lof_below,
            func_above,
            ["var(--harmful-soft)", "var(--untested-soft)", "var(--harmless-soft)"],
            [("broken (LOF)", "st-harmful", "harmful"), ("in between", "st-untested", "grey"), ("working (FUNC)", "st-harmless", "harmless")],
            [signed(lof_below, 3), signed(func_above, 3)],
            lab_score,
            "lab",
        )
    return (
        f'<div class="gauge" role="img" aria-label="Score gauge for {esc(mutant)}">'
        f'<p class="fig caps">Figure. Where {esc(mutant)} lands, left is worse</p>'
        f"{model_track}{lab_track}</div>"
    )


def band_column(verdict: str) -> str:
    return {
        VERDICT_HARMFUL: "looks_harmful",
        VERDICT_HARMLESS: "looks_harmless",
        VERDICT_CANT_TELL: "cant_tell",
    }[verdict]


def track_record_from_file(verdict: str) -> str:
    if band_counts is None:
        return na
    try:
        column = band_column(verdict)
        counts = {str(r["findlay_class"]): int(r[column]) for _, r in band_counts.iterrows()}
        return track_record_text(
            verdict,
            counts[HARMFUL_CLASS],
            counts[HARMLESS_CLASS],
            counts[EXCLUDED_CLASS],
        )
    except Exception:
        return na


def arm_row(arm_id: str) -> dict | None:
    if headline is None:
        return None
    rows = headline[headline["arm"] == arm_id]
    if len(rows) == 0:
        return None
    return rows.iloc[0].to_dict()


def blosum_text(wild_type: str, new_letter: str) -> tuple[str, str]:
    value = blosum.get((wild_type, new_letter))
    if value is None:
        return na, ""
    arm1 = arm_row("arm1")
    if arm1 is None or str(arm1.get("cutoff", "")) == "":
        return f'<span class="num">{value}</span>', "BLOSUM62 rule cutoff not available"
    cut = float(arm1["cutoff"])
    harm = 0.0 - value
    call = "harmful" if harm >= cut else "harmless"
    note = (
        f"Harm = −BLOSUM62 = {signed(harm, 1).lstrip('+')}. "
        f"The BLOSUM62 rule flags harm at or above {cut:.1f} (headline.csv), so it calls this {call}."
    )
    return f'<span class="num">{value}</span>', note


def copy_summary(info: dict, agree: dict) -> str:
    cut_a, cut_b = cutoffs()
    lab_class = info["lab_class"]
    lab_score = info["lab_score"]
    if lab_class and lab_score is not None:
        lab = f"{lab_class}, function score {lab_score:.3f} (LOF < -1.328, FUNC > -0.748)"
    else:
        lab = "not lab-tested"
    value = blosum.get((info["mutant"][0], info["mutant"][-1]))
    lines = [
        f"BRCA Spellcheck: {info['mutant']}",
        f"Model's call: {info['verdict']}",
        f"ESM-2 score (log P(new) - log P(old)): {info['esm_score']:.6f}",
        f"Cutoffs on harm scale, set on start region 1-101: A = {cut_a:.6f}, B = {'none' if cut_b is None else f'{cut_b:.6f}'}",
        f"Lab's answer: {lab}",
        f"Agreement: {agree['label']}",
        f"Region: {'start (tuning)' if info['region'] == 'start' else 'end (held-out test)'}",
        f"BLOSUM62: {na if value is None else value}",
        f"Model: {APP_CONFIG.get('model_id', na)}",
        DISCLAIMER,
    ]
    return "\n".join(lines)


def numbers_panel(info: dict, agree: dict) -> str:
    cut_a, cut_b = cutoffs()
    esm = info["esm_score"]
    lab_class = info["lab_class"]
    lab_score = info["lab_score"]
    wt = info["mutant"][0]
    new = info["mutant"][-1]
    if lab_class and lab_score is not None:
        lab_dd = (
            f'<span class="num">{plain_signed(lab_score, 3)}</span> &middot; {esc(lab_class)}'
            f'<span class="sub">LOF below −1.328 &middot; INT from −1.328 to −0.748 &middot; FUNC above −0.748 (Findlay et al. 2018)</span>'
        )
    else:
        lab_dd = "not lab-tested"
    if info["region"] == "start":
        region_dd = "Start region, positions 1–101. This is the tuning region: the cutoffs were set here."
    else:
        region_dd = "End region, positions 1631–1855. This is the held-out test region, scored once with frozen cutoffs."
    cut_b_text = "none (no band)" if cut_b is None else f"{cut_b:.6f}"
    b_value, b_note = blosum_text(wt, new)
    rows = [
        (
            "Model score",
            f'<span class="num">{plain_signed(esm, 6)}</span><span class="sub">log P(new letter) − log P(old letter), {esc(APP_CONFIG.get("model_id", na))}</span>',
        ),
        (
            "Harm score",
            f'<span class="num">{plain_signed(0.0 - esm, 6)}</span><span class="sub">harm = −score; higher means more likely harmful</span>',
        ),
        (
            "Cutoffs",
            f'A = <span class="num">{cut_a:.6f}</span> &middot; B = <span class="num">{cut_b_text}</span>'
            '<span class="sub">From app_config.json, set on the start region (1–101) only. Harm ≥ B: looks harmful. Harm &lt; A: looks harmless. In between: can\'t tell.</span>',
        ),
        ("Lab's answer", lab_dd),
        ("Region", region_dd),
        ("BLOSUM62 baseline", f'{b_value}<span class="sub">{b_note}</span>'),
        ("Track record of this call", esc(track_record_from_file(info["verdict"])) + '<span class="sub">From band_counts.csv, held-out end region only.</span>'),
    ]
    body = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in rows)
    summary = esc(copy_summary(info, agree))
    return (
        '<details class="numbers"><summary>'
        f'{icon("chev", "chev")}Show the numbers</summary>'
        f'<dl class="kv">{body}</dl>'
        f'<button type="button" class="copy-btn">{icon("copy")}<span>Copy result as text</span></button>'
        f'<pre class="copy-text" aria-hidden="true">{summary}</pre>'
        "</details>"
    )


def result_html(info: dict) -> str:
    verdict = info["verdict"]
    lab_class = info["lab_class"]
    lab_score = info["lab_score"]
    position = int(info["mutant"][1:-1])
    wt = info["mutant"][0]
    new = info["mutant"][-1]
    vkind, vcss = verdict_style[verdict]
    agree = agreement(verdict, lab_class)
    where = "start region &middot; tuning" if info["region"] == "start" else "end region &middot; held-out test"
    story = " ".join(
        [
            f"{esc(info['mutant'])} swaps {aa_names[wt]} ({wt}) for {aa_names[new]} ({new}) at position "
            f'<span class="num">{position}</span>, in {region_words(info["region"])}.',
            model_sentence(info["esm_score"], verdict),
            lab_sentence(lab_class, lab_score, agree["key"]),
        ]
    )
    return (
        '<article class="card result-card" aria-live="polite">'
        '<div class="result-top">'
        f'<h3 class="mutant">{esc(info["mutant"])}</h3>'
        f'<p class="where caps">position {position} &middot; {where}</p></div>'
        f'<span class="verdict {vcss}">{icon(vkind)}{esc(verdict)}</span>'
        f'<p class="plain">{verdict_plain[verdict]}</p>'
        '<div class="agree">'
        f'<span class="agree-badge {agree["css"]}">{icon(agree["icon"])}{agree["label"]}</span>'
        f'<p class="agree-mean"><span class="q">What does this mean?</span> {agree["mean"]}</p></div>'
        f'<p class="story">{story}</p>'
        f'{gauge(info["esm_score"], lab_class, lab_score, info["mutant"])}'
        f"{numbers_panel(info, agree)}"
        "</article>"
    )


refusal_titles: list[tuple[str, str]] = [
    ("Type a mutation", "Nothing to check yet."),
    ("I could not read", "That doesn't look like a protein change."),
    ("is not a change", "That one doesn't change anything."),
    ("is outside BRCA1", "That position is off the end of the protein."),
    ("is not one of the", "We didn't score that position."),
    ("BRCA1 has", "Wrong starting letter."),
    ("End-region answers", "Not yet."),
]


def nearest_suggestion(message: str) -> list[str]:
    found = re.search(r"Position (\d+) is not one of", message)
    if found is None or not scored_positions:
        return []
    asked = int(found.group(1))
    best = min(scored_positions, key=lambda p: abs(p - asked))
    for mutant in lab_tested_list:
        if int(mutant[1:-1]) == best:
            return [mutant]
    return []


def refusal_html(message: str) -> str:
    title = "Something broke on our side."
    for prefix, text in refusal_titles:
        if prefix in message:
            title = text
            break
    picks = []
    meant = re.search(r"Did you mean (\w+)\?", message)
    if meant:
        picks.append(meant.group(1))
    picks += nearest_suggestion(message)
    for default in ("C61G", "D2E", "p.Cys64Gly"):
        if default not in picks:
            picks.append(default)
    chips = "".join(
        f'<button type="button" class="try-chip" data-q="{esc(q)}">{esc(q)}</button>' for q in picks[:4]
    )
    return (
        '<article class="card refusal" role="status">'
        f"<h3>{esc(title)}</h3><p>{esc(message)}</p>"
        f'<div class="try"><span class="caps">Try one of these</span>{chips}</div></article>'
    )


def empty_html() -> str:
    return (
        '<div class="empty"><h3>Nothing checked yet.</h3>'
        "<p>Type a change on the left, or pick one of the examples. C61G is the usual first try.</p>"
        "<p>Two of the examples are cases where the model got it wrong. We left them in on purpose.</p></div>"
    )


def render_result(text) -> str:
    try:
        wild_type, position, new_letter = parse_input(text)
        info = lookup(wild_type, position, new_letter)
        return result_html(info)
    except ValueError as exc:
        return refusal_html(str(exc))
    except Exception:
        return refusal_html(_R7_MESSAGE)


def failure_row(kind: str) -> dict | None:
    if failures is None:
        return None
    rows = failures[(failures["kind"] == kind) & (failures["rank"].astype(str) == "1")]
    if len(rows) == 0:
        return None
    return rows.iloc[0].to_dict()


def stat_html(value: str, label: str, raw=None, k: int = 0) -> str:
    shown = value
    if raw is not None:
        shown = f'<span class="num" data-count="{raw}" data-dec="0" data-fmt="int">{value}</span>'
    return f'<div class="stat rv" style="--d:{300 + k * 90}ms"><dt>{shown}</dt><dd>{label}</dd></div>'


def ticker_html() -> str:
    start = "".join(_WT_BY_POS[p] for p in scored_positions if p <= START_REGION[1])
    end = "".join(_WT_BY_POS[p] for p in scored_positions if p >= END_REGION[0])
    rows = []
    for k, seq in enumerate((start, end)):
        text = esc(seq + " " + seq)
        rows.append(f'<div class="tick-row r{k}"><span>{text}</span><span>{text}</span></div>')
    return '<div class="ticker" aria-hidden="true">' + "".join(rows) + "</div>"


def title_html(text: str) -> str:
    letters = "".join(
        f'<span class="ch" style="--i:{i}">{"&nbsp;" if c == " " else esc(c)}</span>' for i, c in enumerate(text)
    )
    return f'<h1 aria-label="{esc(text)}"><span class="title-in" aria-hidden="true">{letters}</span></h1>'


def hero_html() -> str:
    n_positions = APP_CONFIG.get("n_tested_positions")
    end_n = bootstrap.get("n_positions") if bootstrap else None
    stats = (
        stat_html(esc(n_positions) if n_positions is not None else na, "positions scored", n_positions, 0)
        + stat_html(f"{len(lab_tested_list):,}", "swaps with a lab answer", len(lab_tested_list), 1)
        + stat_html(
            esc(end_n) if end_n is not None else na,
            "held-out test positions<br>1631–1855, a spec constant",
            end_n,
            2,
        )
    )
    legend_items = [
        ("harmful", "st-harmful", "Looks harmful", "the model's call"),
        ("harmless", "st-harmless", "Looks harmless", "the model's call"),
        ("cant", "st-cant", "Can't tell", "score in the middle zone"),
        ("untested", "st-untested", "Not lab-tested", "no answer key"),
        ("agree", "st-agree", "Agree", "model and lab match"),
        ("disagree", "st-disagree", "Disagree", "a miss or a false alarm"),
    ]
    legend = "".join(
        f'<li class="rv" style="--d:{350 + j * 70}ms"><span class="{css}">{icon(kind)}</span><span><b>{label}</b> &middot; {note}</span></li>'
        for j, (kind, css, label, note) in enumerate(legend_items)
    )
    return (
        '<div class="scroll-progress" aria-hidden="true"><i></i></div>'
        '<header class="hero">'
        f"{ticker_html()}"
        "<div>"
        '<p class="eyebrow caps rv" style="--d:0ms">A prototype &middot; BRCA1 &middot; ESM-2</p>'
        f"{title_html('BRCA Spellcheck')}"
        '<p class="lede rv" style="--d:180ms">Type one letter change in the BRCA1 protein. We show the model\'s call, the lab\'s answer, and whether the two agree.</p>'
        '<p class="exam rv" style="--d:240ms">BRCA1 is our exam, not the use case: it is one of the few genes with a lab answer key, so we can check the model here.</p>'
        f'<dl class="stats">{stats}</dl></div>'
        f'<aside class="legend"><h2 class="caps rv" style="--d:300ms">How to read the marks</h2><ul>{legend}</ul></aside>'
        "</header>"
    )


def input_intro_html() -> str:
    return (
        '<div class="input-intro rv" style="--d:420ms"><h2>Check a mutation</h2>'
        "<p>Old letter, position, new letter. C61G, c61g, Cys61Gly and p.Cys61Gly all work. "
        "We scored positions 1–101 and 1631–1855.</p></div>"
    )


def chip_note_html() -> str:
    miss = failure_row("miss")
    alarm = failure_row("false_alarm")
    parts = []
    if miss:
        parts.append(f"<b>{esc(miss['mutant'])}</b> is our worst miss")
    if alarm:
        parts.append(f"<b>{esc(alarm['mutant'])}</b> our worst false alarm")
    if not parts:
        return ""
    return '<p class="chip-note">' + ", ".join(parts) + ".</p>"


def meter(label: str, count: int, total: int, css: str, k: int = 0) -> str:
    share = count / total if total else 0.0
    return (
        f'<div class="metric"><div class="m-top"><span>{label}</span>'
        f'<span><span class="num" data-count="{count}" data-dec="0" data-fmt="int">{count}</span> <span class="of">of {total}</span></span></div>'
        f'<div class="meter" role="img" aria-label="{label}: {count} of {total}"><i class="{css}" style="--v:{share:.4f};--d:{k * 120}ms"></i></div></div>'
    )


arm_titles: dict[str, str] = {
    "arm0": "Always says harmless",
    "arm1": "BLOSUM62 rule",
    "arm2": "ESM-2",
    "arm3": "ESM-2 plus extras",
}


def scoreboard_html() -> str:
    if headline is None:
        return '<div class="empty"><h3>Scoreboard not available.</h3><p>headline.csv is missing, so we show nothing rather than guess.</p></div>'
    first = headline.iloc[0]
    n_harmful = int(first["n_harmful"])
    n_harmless = int(first["n_harmless"])
    rows = []
    for row_k, (_, arm) in enumerate(headline.iterrows()):
        arm_id = str(arm["arm"])
        sub = esc(arm["arm_name"])
        if str(arm["cutoff"]) != "":
            sub += f' &middot; cutoff <span class="num">{float(arm["cutoff"]):.6f}</span>'
        else:
            sub += " &middot; no cutoff"
        if arm_id == "arm2":
            sub += " &middot; the model in this demo"
        bars = (
            meter("Harmful caught", int(arm["harmful_caught"]), int(arm["n_harmful"]), "caught", 0)
            + meter("Harmful missed", int(arm["harmful_missed"]), int(arm["n_harmful"]), "missed", 1)
            + meter("False alarms", int(arm["false_alarms"]), int(arm["n_harmless"]), "alarm", 2)
        )
        featured = " featured" if arm_id == "arm2" else ""
        rows.append(
            f'<div class="arm{featured} rv" style="--d:{row_k * 110}ms"><div><h3>{esc(arm_titles.get(arm_id, arm["arm_name"]))}</h3><p class="arm-sub">{sub}</p></div>'
            f'<div class="arm-bars">{bars}</div>'
            f'<div class="arm-rate"><span class="big"><span class="num" data-count="{float(arm["catch_rate_percent"])}" data-dec="1" data-fmt="pct">{float(arm["catch_rate_percent"]):.1f}%</span></span><span class="lbl">catch rate</span></div></div>'
        )
    sentence = sentences.get("S6")
    sentence_html = (
        f'<p class="sentence rv">{esc(sentence)}</p><p class="small rv" style="--d:120ms">Result sentence S6, worded before the test run and filled in from it.</p>'
        if sentence
        else f'<p class="sentence">{na}</p>'
    )
    rank_rows = ""
    if auroc is not None:
        for _, r in auroc.iterrows():
            rank_rows += f'<tr><td>{esc(arm_titles.get(str(r["arm"]), r["arm"]))}</td><td class="num">{float(r["auroc"]):.4f}</td></tr>'
    else:
        rank_rows = f"<tr><td>{na}</td><td></td></tr>"
    boot = ""
    if bootstrap and "esm_minus_blosum" in bootstrap:
        b = bootstrap["esm_minus_blosum"]
        boot = (
            f'<p class="small">ESM-2 minus BLOSUM62: <span class="num">{b["point"]:.3f}</span> '
            f'(2.5th to 97.5th percentile <span class="num">{b["low"]:.3f}</span> to <span class="num">{b["high"]:.3f}</span>, '
            f'{int(bootstrap.get("n_resamples_used", 0)):,} resamples of {esc(bootstrap.get("unit", "positions"))}).</p>'
        )
    ranking = (
        '<div class="rv" style="--d:160ms"><table class="rank-table"><caption class="caps">Ranking quality, end region (AUROC, 0.5 is a coin flip)</caption>'
        f"<tbody>{rank_rows}</tbody></table>{boot}</div>"
    )
    head = (
        '<div class="section-head"><h2 class="rv">How the model did on the end region it never saw.</h2>'
        f'<p class="rv" style="--d:120ms">We set every cutoff on positions'
        f" 1–101, froze them, then scored positions 1631–1855 once. "
        f'There the lab gives a clear answer for <span class="num">{n_harmful}</span> harmful and <span class="num">{n_harmless}</span> harmless mutations. '
        "Misses and false alarms are shown side by side, in grey, not hidden.</p></div>"
    )
    return (
        head
        + f'<div class="board">{"".join(rows)}</div>'
        + f'<div class="board-note"><div>{sentence_html}</div>{ranking}</div>'
        + failure_html()
    )


def failure_html() -> str:
    miss = failure_row("miss")
    if miss is None:
        return ""
    position = int(miss["position"])
    near = 0
    if failures is not None:
        misses = failures[failures["kind"] == "miss"]
        near = int(((misses["position"].astype(int) - position).abs() <= 1).sum())
        total_misses = len(misses)
    others = []
    if str(miss.get("flag_arm1")) == "0":
        others.append("the BLOSUM62 rule")
    if str(miss.get("flag_arm3")) == "0":
        others.append("ESM-2 plus extras")
    joined = " and ".join(others)
    also = f" {joined[:1].upper() + joined[1:]} missed it too." if others else ""
    cluster = (
        f" {near} of our {total_misses} worst misses sit at positions {position} and {position + 1}."
        if near > 1
        else ""
    )
    return (
        '<section class="failure rv"><div>'
        f'<p class="tag">{icon("disagree")}Failure: worst miss on the end region</p>'
        f'<h3 class="mutant">{esc(miss["mutant"])}</h3>'
        f'<p class="small">From failures.csv, rank 1 of the misses.</p></div>'
        "<div>"
        f'<p>The lab measured a function score of <span class="num">{signed(float(miss["function_score"]), 3)}</span>, '
        f'well into broken. The model scored it <span class="num">{signed(float(miss["esm_score"]), 2)}</span>, '
        f'an ordinary-looking swap, and said "{esc(miss["band_arm2"])}".{also}{cluster}</p>'
        '<p class="guess"><b>Our guess, not tested:</b> the lab edited DNA, not protein. This DNA change may upset '
        "how the gene's message is spliced, and a model that only reads protein letters cannot see that.</p>"
        "</div></section>"
    )


def heat_head_html() -> str:
    return (
        '<div class="section-head"><h2 class="rv">Every swap we scored, in one picture.</h2>'
        '<p class="rv" style="--d:120ms">'
        "Columns are positions, rows are the 20 possible new letters. Red means the model found the swap "
        "surprising (a negative score), blue means it did not. Black dots mark swaps the lab found broken; "
        "grey dots mark the other lab-tested swaps. The colour scale stops at ±12.</p></div>"
    )


def heat_caption_html() -> str:
    return (
        '<p class="small heat rv"><span class="caps">Figure 1.</span> ESM-2 scores at the lab-tested positions, '
        "start region above, end region below. Saved by the pipeline as outputs/heatmap.png.</p>"
    )


dash: str = ' stroke-dasharray="3 3"'
qfill: str = "var(--accent-soft)"
step_svgs: list[str] = [
    '<svg class="dia" viewBox="0 0 220 64" aria-hidden="true">'
    + "".join(
        f'<g class="pc" style="--j:{i}"><rect x="{4 + i * 30}" y="16" width="26" height="30" rx="4" fill="{qfill if c == "?" else "none"}" stroke="currentColor" stroke-width="1.25"{dash if c == "?" else ""}/>'
        f'<text x="{17 + i * 30}" y="36" text-anchor="middle">{c}</text></g>'
        for i, c in enumerate("LKCE?QR")
    )
    + "</svg>",
    '<svg class="dia" viewBox="0 0 220 64" aria-hidden="true">'
    '<line class="pc" style="--j:0" x1="10" y1="54" x2="210" y2="54" stroke="currentColor" stroke-width="1"/>'
    '<rect class="pc grow" style="--j:1" x="40" y="18" width="40" height="36" fill="var(--harmless-soft)" stroke="var(--harmless)" stroke-width="1.25"/>'
    '<rect class="pc grow" style="--j:3" x="130" y="44" width="40" height="10" fill="var(--harmful-soft)" stroke="var(--harmful)" stroke-width="1.25"/>'
    '<text class="pc" style="--j:2" x="60" y="13" text-anchor="middle">old C</text>'
    '<text class="pc" style="--j:4" x="150" y="38" text-anchor="middle">new G</text>'
    "</svg>",
    '<svg class="dia" viewBox="0 0 220 64" aria-hidden="true">'
    '<g class="pc wipe" style="--j:0"><rect x="6" y="26" width="80" height="10" fill="var(--harmful-soft)"/>'
    '<rect x="86" y="26" width="40" height="10" fill="var(--cant-soft)"/>'
    '<rect x="126" y="26" width="88" height="10" fill="var(--harmless-soft)"/></g>'
    '<g class="pc drop" style="--j:4"><line x1="86" y1="16" x2="86" y2="46" stroke="currentColor" stroke-dasharray="3 3"/>'
    '<text x="86" y="60" text-anchor="middle">B</text></g>'
    '<g class="pc drop" style="--j:6"><line x1="126" y1="16" x2="126" y2="46" stroke="currentColor" stroke-dasharray="3 3"/>'
    '<text x="126" y="60" text-anchor="middle">A</text></g>'
    "</svg>",
    '<svg class="dia" viewBox="0 0 220 64" aria-hidden="true">'
    '<g class="pc" style="--j:0"><rect x="6" y="20" width="60" height="24" rx="4" fill="none" stroke="currentColor" stroke-width="1.25"/>'
    '<text x="36" y="36" text-anchor="middle">1–101</text></g>'
    '<g class="pc wipe" style="--j:2"><path d="M72 32h50" stroke="currentColor" stroke-width="1.25" fill="none"/><path d="M116 27l7 5-7 5" stroke="currentColor" stroke-width="1.25" fill="none"/></g>'
    '<g class="pc" style="--j:5"><rect x="128" y="20" width="86" height="24" rx="4" fill="var(--accent-soft)" stroke="var(--accent)" stroke-width="1.25"/>'
    '<text x="171" y="36" text-anchor="middle">1631–1855</text></g>'
    "</svg>",
]


def how_html() -> str:
    c61g = _ROW_BY_POS_NEW.get((61, "G"))
    example = (
        f' C61G scores <span class="num">{signed(c61g["esm_score"], 2)}</span>.' if c61g else ""
    )
    steps = [
        (
            "Hide one letter",
            "We give ESM-2, a protein language model, a 500-letter stretch of BRCA1 with one position blanked out, "
            "and ask how likely each of the 20 letters is to belong there.",
        ),
        (
            "Compare the old and new letter",
            "The score is log P(new letter) minus log P(old letter). Below zero means the model finds the new letter "
            "more surprising than the one nature put there." + example,
        ),
        (
            "Draw two lines on the start region",
            f"On positions 1–101 we set cutoff A to still catch at least 90% of the {start_n_harmful} harmful mutations there, "
            f"and cutoff B to flag at most 10% of the {start_n_harmless} harmless ones. Between the lines the model says \"can't tell\".",
        ),
        (
            "Grade once on the end region",
            "Then we froze both lines and scored positions 1631–1855 one time, against the lab. "
            "The scoreboard tab is that one run, mistakes included.",
        ),
    ]
    items = "".join(
        f'<li class="rv" style="--d:{i * 90}ms"><span class="n">0{i + 1}</span><div><h3>{t}</h3><p>{p}</p></div>{step_svgs[i]}</li>'
        for i, (t, p) in enumerate(steps)
    )
    head = (
        '<div class="section-head"><h2 class="rv">Four steps, no training.</h2>'
        '<p class="rv" style="--d:120ms">'
        "The model never saw a lab label. It only learned which letters tend to sit where in natural proteins. "
        "We borrowed the lab's answers to grade it.</p></div>"
    )
    return head + f'<ol class="steps">{items}</ol>'


def credits_html() -> str:
    limits = [
        (
            "Blind to splicing",
            "The lab edited DNA. Some DNA changes break how the gene's message is spliced, which a model that "
            "reads only protein letters cannot see. We could not test this.",
        ),
        (
            "Cells in a dish, not people",
            "The answer key comes from cells grown in a lab. A 2025 preprint notes that results can depend on cell type.",
        ),
        (
            "Stronger tools exist",
            "AlphaMissense is free to download and is a stronger variant effect predictor. Zero-shot scoring like ours "
            "is not a new method; ProteinGym already evaluates ESM-2 on this same table.",
        ),
    ]
    limit_items = "".join(
        f'<li class="rv" style="--d:{k * 100}ms"><h3>{t}</h3><p>{p}</p></li>' for k, (t, p) in enumerate(limits)
    )

    def link(text: str, url: str) -> str:
        return f'<a href="{url}" target="_blank" rel="noopener">{text}</a>'

    groups = [
        (
            "Data",
            [
                "Findlay GM et al. (2018), Nature 562:217–222, saturation genome editing of BRCA1: the lab answer key. "
                + link("doi", "https://doi.org/10.1038/s41586-018-0461-z"),
                "BRCA1 saturation genome editing resource (classes and cutoffs). "
                + link("sge.gs.washington.edu", "https://sge.gs.washington.edu/BRCA1/"),
                "ProteinGym (Notin P et al., 2023), source of BRCA1_HUMAN_Findlay_2018.csv. "
                + link("github", "https://github.com/OATML-Markslab/ProteinGym"),
            ],
        ),
        (
            "Models and methods",
            [
                "ESM-2 (Lin Z et al., 2023, Science 379:1123–1130), "
                + esc(APP_CONFIG.get("model_id", na))
                + ". "
                + link("model card", "https://huggingface.co/facebook/esm2_t33_650M_UR50D"),
                "Masked-marginal scoring: Meier J et al. (2021), NeurIPS. "
                + link("doi", "https://doi.org/10.1101/2021.07.09.450648"),
                "BLOSUM62: Henikoff S, Henikoff JG (1992), PNAS 89:10915–10919.",
                "Hydrophobicity: Kyte J, Doolittle RF (1982). Residue volumes: Zamyatnin AA (1972).",
            ],
        ),
        (
            "Prior work we build on",
            [
                "ProteinGym, which already scores ESM-2 on this table.",
                "Evo 2 BRCA1 zero-shot notebook (Arc Institute, NVIDIA). "
                + link("github", "https://github.com/ArcInstitute/evo2"),
                "AlphaMissense (Cheng J et al., 2023, Science 381:eadg7492).",
                "Cell-type context preprint (2025). "
                + link("medRxiv", "https://www.medrxiv.org/content/10.1101/2025.08.11.25333423.full.pdf"),
            ],
        ),
        (
            "AI tools",
            [
                "Claude Code (Anthropic): planning, the build spec, review, and this demo's redesign.",
                "Gemini in Google Antigravity: writing and running the project code from the spec.",
            ],
        ),
        (
            "Software and type",
            [
                "PyTorch, Transformers, NumPy, pandas, SciPy, scikit-learn, Biopython, Plotly, Matplotlib, Gradio, pytest, Jupyter.",
                "Instrument Serif, IBM Plex Sans and JetBrains Mono, via Google Fonts.",
            ],
        ),
    ]
    credit_blocks = "".join(
        f'<div class="rv" style="--d:{k * 80}ms"><h3 class="caps">{g}</h3><ul>' + "".join(f"<li>{i}</li>" for i in items) + "</ul></div>"
        for k, (g, items) in enumerate(groups)
    )
    head = (
        '<div class="section-head"><h2 class="rv">What it can\'t do, and whose work it stands on.</h2>'
        '<p class="rv" style="--d:120ms">'
        "This is a prototype that shows how far a sequence-only model can be trusted on one well-tested gene. "
        "It does not diagnose anything.</p></div>"
    )
    return head + (
        f'<div class="cols"><div><ol class="limits">{limit_items}</ol></div>'
        f'<div class="credits">{credit_blocks}</div></div>'
    )


def notes_html() -> str:
    return (
        '<div class="notes-bar" role="note"><div class="inner">'
        f'<span class="disc">{esc(DISCLAIMER)}</span>'
        f"<span>{esc(DEMO_NOTE)}</span></div></div>"
    )


examples: list[str] = ["C61G", "D2E", "p.Cys64Gly"]
for kind in ("miss", "false_alarm"):
    row = failure_row(kind)
    if row is not None:
        examples.append(str(row["mutant"]))


def pick_example(q: str) -> tuple[str, str]:
    return q, render_result(q)


def surprise() -> tuple[str, str]:
    if not lab_tested_list:
        return "", refusal_html("There are no lab-tested swaps loaded, so there is nothing to pick from.")
    q = str(pd.Series(lab_tested_list).sample(1).iloc[0])
    return q, render_result(q)


def on_load(request: gr.Request):
    q = ""
    tab = "check"
    try:
        q = (request.query_params.get("q") or "").strip()[:40]
        tab = request.query_params.get("tab") or "check"
    except Exception:
        pass
    if tab not in ("check", "score", "heat", "how", "limits"):
        tab = "check"
    shown = render_result(q) if q else empty_html()
    return q, shown, gr.Tabs(selected=tab)


font_link: str = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">'
    '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
    '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&'
    'family=Instrument+Serif&family=JetBrains+Mono:wght@400;500&display=swap">'
    '<meta name="description" content="A prototype that compares a protein language model\'s call on BRCA1 '
    'letter changes with lab measurements.">'
)

head_script: str = """
<script>
(function () {
  var root = document.documentElement;
  root.classList.add("js");
  try {
    if (new URLSearchParams(location.search).get("motion") === "off") root.classList.add("no-motion");
    if (new URLSearchParams(location.search).get("fx") === "max") root.classList.add("fx");
  } catch (e) {}
  function still() {
    return root.classList.contains("no-motion") || window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  }
  function fmt(v, d, kind) {
    if (kind === "int") return Math.round(v).toLocaleString("en-US");
    if (kind === "pct") return v.toFixed(d) + "%";
    var s = Math.abs(v).toFixed(d);
    if (v < 0 && parseFloat(s) !== 0) return "\\u2212" + s;
    if (v > 0 && parseFloat(s) !== 0) return "+" + s;
    return s;
  }
  function count(el) {
    if (el.dataset.done) return;
    el.dataset.done = "1";
    var to = parseFloat(el.dataset.count), d = parseInt(el.dataset.dec || "0", 10);
    if (isNaN(to) || still()) return;
    var t0 = performance.now() + (parseInt(el.dataset.delay || "0", 10)), dur = el.dataset.fmt ? 900 : 450, last = el.textContent;
    function step(t) {
      var p = Math.max(0, Math.min(1, (t - t0) / dur));
      p = 1 - Math.pow(1 - p, 3);
      el.textContent = p < 1 ? fmt(to * p, d, el.dataset.fmt) : last;
      if (p < 1) requestAnimationFrame(step);
    }
    requestAnimationFrame(step);
  }
  var io = "IntersectionObserver" in window ? new IntersectionObserver(function (entries) {
    entries.forEach(function (en) {
      if (en.isIntersecting) {
        en.target.classList.add("in");
        io.unobserve(en.target);
        if (en.target.dataset && en.target.dataset.count !== undefined) count(en.target);
        if (en.target.querySelectorAll) en.target.querySelectorAll("[data-count]").forEach(count);
      }
    });
  }, { threshold: 0.2 }) : null;
  function scan(node) {
    if (!node || !node.querySelectorAll) return;
    node.querySelectorAll(".result-card [data-count]").forEach(count);
    var reveal = [];
    if (node.matches && node.matches(".rv")) reveal.push(node);
    node.querySelectorAll(".rv:not(.seen)").forEach(function (el) { reveal.push(el); });
    reveal.forEach(function (el) {
      el.classList.add("seen");
      if (io && !still()) io.observe(el); else { el.classList.add("in"); }
    });
    var card = node.matches && node.matches(".result-card, .refusal") ? node : node.querySelector(".result-card, .refusal");
    if (card && window.innerWidth < 900 && document.readyState === "complete") {
      card.scrollIntoView({ behavior: still() ? "auto" : "smooth", block: "start" });
    }
    node.querySelectorAll(".meter i:not(.seen)").forEach(function (el) {
      el.classList.add("seen");
      if (io && !still()) io.observe(el); else el.classList.add("in");
    });
  }
  function light() {
    if (document.body && document.body.classList.contains("dark")) document.body.classList.remove("dark");
  }
  var mo = new MutationObserver(function (list) {
    light();
    list.forEach(function (m) { m.addedNodes.forEach(scan); });
  });
  var bar = null, ticking = false;
  function progress() {
    ticking = false;
    if (!bar) bar = document.querySelector(".scroll-progress i");
    if (!bar) return;
    var h = document.documentElement.scrollHeight - innerHeight;
    var v = h > 0 ? Math.min(1, Math.max(0, scrollY / h)) : 0;
    bar.style.transform = "scaleX(" + v.toFixed(4) + ")";
  }
  window.addEventListener("scroll", function () {
    if (!ticking) { ticking = true; requestAnimationFrame(progress); }
  }, { passive: true });
  function fx() {
    var hero = null, title = null, legend = null, tick = null, pending = false;
    function para() {
      pending = false;
      hero = hero || document.querySelector(".hero");
      if (!hero) return;
      title = title || hero.querySelector("h1");
      legend = legend || hero.querySelector(".legend");
      tick = tick || document.querySelector(".ticker");
      var y = Math.min(scrollY, 900);
      if (title) { title.style.transform = "translate3d(0," + (y * -0.25).toFixed(1) + "px,0)"; title.style.opacity = Math.max(0, 1 - y / 420).toFixed(3); }
      if (legend) legend.style.transform = "translate3d(0," + (y * 0.1).toFixed(1) + "px,0)";
      if (tick) tick.style.transform = "translate3d(0," + (y * 0.4).toFixed(1) + "px,0)";
    }
    window.addEventListener("scroll", function () {
      if (!pending) { pending = true; requestAnimationFrame(para); }
    }, { passive: true });
    document.addEventListener("pointermove", function (ev) {
      var card = ev.target.closest && ev.target.closest(".result-card, .failure, .arm");
      if (!card) return;
      var r = card.getBoundingClientRect();
      card.style.setProperty("--mx", (ev.clientX - r.left) + "px");
      card.style.setProperty("--my", (ev.clientY - r.top) + "px");
      if (card.classList.contains("result-card")) {
        var rx = ((ev.clientY - r.top) / r.height - 0.5) * -3;
        var ry = ((ev.clientX - r.left) / r.width - 0.5) * 3;
        card.style.transform = "perspective(1400px) rotateX(" + rx.toFixed(2) + "deg) rotateY(" + ry.toFixed(2) + "deg)";
      }
    });
    document.addEventListener("pointerout", function (ev) {
      var card = ev.target.closest && ev.target.closest(".result-card");
      if (card && !card.contains(ev.relatedTarget)) card.style.transform = "";
    });
  }
  function start() {
    light();
    root.classList.add("ready");
    if (root.classList.contains("fx") && !still()) fx();
    mo.observe(document.body, { childList: true, subtree: true, attributes: true, attributeFilter: ["class"] });
    scan(document.body);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start); else start();
  document.addEventListener("click", function (ev) {
    var chip = ev.target.closest && ev.target.closest(".try-chip");
    if (chip) {
      var box = document.querySelector("#mut-input textarea, #mut-input input");
      var btn = document.querySelector("#check-btn");
      if (box && btn) {
        var proto = Object.getPrototypeOf(box);
        var setter = Object.getOwnPropertyDescriptor(proto, "value").set;
        setter.call(box, chip.dataset.q);
        box.dispatchEvent(new Event("input", { bubbles: true }));
        setTimeout(function () { btn.click(); }, 30);
        box.focus();
      }
      return;
    }
    var copy = ev.target.closest && ev.target.closest(".copy-btn");
    if (copy) {
      var pre = copy.parentElement.querySelector(".copy-text");
      var label = copy.querySelector("span");
      var text = pre ? pre.textContent : "";
      var done = function () { label.textContent = "Copied"; setTimeout(function () { label.textContent = "Copy result as text"; }, 1600); };
      var fallback = function () {
        var area = document.createElement("textarea");
        area.value = text;
        area.setAttribute("readonly", "");
        area.style.position = "fixed";
        area.style.opacity = "0";
        document.body.appendChild(area);
        area.select();
        var ok = false;
        try { ok = document.execCommand("copy"); } catch (e) {}
        document.body.removeChild(area);
        if (ok) done(); else label.textContent = "Copy failed, select the numbers by hand";
      };
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(done, fallback);
      } else {
        fallback();
      }
    }
  });
})();
</script>
"""

theme = gr.themes.Base(
    primary_hue=gr.themes.colors.slate,
    neutral_hue=gr.themes.colors.stone,
    radius_size=gr.themes.sizes.radius_md,
    font=["IBM Plex Sans", "-apple-system", "Segoe UI", "Roboto", "sans-serif"],
    font_mono=["JetBrains Mono", "ui-monospace", "Consolas", "monospace"],
).set(
    body_background_fill="#f4f2ec",
    body_text_color="#1c2024",
    block_background_fill="#fbfaf6",
    block_border_color="#d6d2c8",
    input_background_fill="#fbfaf6",
    button_primary_background_fill="#1f4a73",
    button_primary_text_color="#f7f6f2",
)

css_text: str = (APP_DIR / "style.css").read_text(encoding="utf-8") if (APP_DIR / "style.css").exists() else ""


with gr.Blocks(title="BRCA Spellcheck") as demo:
    gr.HTML(hero_html())
    with gr.Tabs(selected="check", elem_classes=["tabs-main"]) as tabs:
        with gr.Tab("Check", id="check"):
            with gr.Row(elem_classes=["check-row"]):
                with gr.Column(scale=4, min_width=280, elem_classes=["input-col"]):
                    gr.HTML(input_intro_html())
                    mutation_input = gr.Textbox(
                        label="Mutation",
                        placeholder="C61G",
                        max_lines=1,
                        elem_id="mut-input",
                        elem_classes=["mut-input", "reveal", "d2"],
                    )
                    with gr.Row(elem_classes=["btn-row"]):
                        check_button = gr.Button("Check this swap", elem_id="check-btn", elem_classes=["btn-primary"], scale=2)
                        surprise_button = gr.Button("Surprise me", elem_classes=["btn-quiet"], scale=1)
                    gr.HTML('<p class="chip-head caps">Examples</p>')
                    with gr.Row(elem_classes=["chip-row"]):
                        chip_buttons = [
                            gr.Button(q, size="sm", elem_classes=["chip"], min_width=0, scale=0) for q in examples
                        ]
                    gr.HTML(chip_note_html())
                with gr.Column(scale=7, min_width=300):
                    result_output = gr.HTML(empty_html(), elem_classes=["result-out"], js_on_load=None)
        with gr.Tab("Scoreboard", id="score"):
            gr.HTML(scoreboard_html(), js_on_load=None)
        with gr.Tab("Heat map", id="heat"):
            gr.HTML(heat_head_html(), js_on_load=None)
            heat_file = heatmap_path()
            if heat_file is not None:
                gr.Image(
                    value=str(heat_file),
                    show_label=False,
                    interactive=False,
                    elem_classes=["heat-img"],
                    buttons=["fullscreen"],
                )
                gr.HTML(heat_caption_html(), js_on_load=None)
            else:
                gr.HTML('<div class="empty"><h3>Heat map not available.</h3><p>heatmap.png is missing from outputs/ and app/results/.</p></div>')
        with gr.Tab("How it works", id="how"):
            gr.HTML(how_html(), js_on_load=None)
        with gr.Tab("Limits", id="limits"):
            gr.HTML(credits_html(), js_on_load=None)
    gr.HTML(notes_html(), elem_classes=["notes-wrap"], js_on_load=None)

    check_button.click(fn=render_result, inputs=mutation_input, outputs=result_output)
    mutation_input.submit(fn=render_result, inputs=mutation_input, outputs=result_output)
    surprise_button.click(fn=surprise, inputs=None, outputs=[mutation_input, result_output])
    for chip, q in zip(chip_buttons, examples):
        chip.click(fn=lambda q=q: pick_example(q), inputs=None, outputs=[mutation_input, result_output])
    demo.load(fn=on_load, inputs=None, outputs=[mutation_input, result_output, tabs])


if __name__ == "__main__":
    demo.launch(
        theme=theme,
        css=css_text,
        head=font_link + head_script,
        footer_links=[],
    )
