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
explorer_url: str = ""
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


def check_bundle() -> list[str]:
    differ = []
    pairs = [(repo_results_dir / n, results_dir / n) for n in bundle_names]
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
    return f'<span class="num">{text}</span>'


icon_paths: dict[str, str] = {
    "harmful": '<path d="M8 2.5 14.2 13.5H1.8Z"/>',
    "harmless": '<circle cx="8" cy="8" r="5.6"/>',
    "cant": '<path d="M8 1.8 14.2 8 8 14.2 1.8 8Z"/>',
    "untested": '<circle cx="8" cy="8" r="5.6" stroke-dasharray="2.2 2.2"/>',
    "agree": '<path d="M2.8 8.6 6.4 12.2 13.2 4.4"/>',
    "disagree": '<path d="M3.5 3.5 12.5 12.5M12.5 3.5 3.5 12.5"/>',
    "grey": '<path d="M2.5 8h11"/>',
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


squig_paths: dict[str, str] = {
    "harmful": "M0 4 Q2 0 4 4 T8 4 T12 4 T16 4 T20 4 T24 4",
    "harmless": "M1 4H23",
    "cant": "M1 4H23",
    "probe": "M1 4H23",
}


def seq_html(position: int, new_letter: str | None = None, kind: str = "probe") -> str:
    span = [p for p in range(position - 6, position + 7) if p in _WT_BY_POS]
    cells = []
    for i, p in enumerate(span):
        wt = _WT_BY_POS[p]
        if p != position:
            cells.append(f'<span class="aa" style="--i:{i}">{esc(wt)}</span>')
            continue
        if new_letter:
            inner = f'<span class="old">{esc(wt)}</span><span class="new">{esc(new_letter)}</span>'
        else:
            inner = f'<span class="new">{esc(wt)}</span>'
        cells.append(
            f'<span class="aa hit k-{kind}" style="--i:{i}">{inner}'
            f'<svg class="squig" viewBox="0 0 24 8" preserveAspectRatio="none" aria-hidden="true"><path d="{squig_paths[kind]}"/></svg>'
            f'<span class="pos">{p}</span></span>'
        )
    letters = "".join(_WT_BY_POS[p] for p in span)
    if new_letter:
        label = f"BRCA1 letters {span[0]} to {span[-1]}: {letters}. The {_WT_BY_POS[position]} at {position} becomes {new_letter}."
    else:
        label = f"BRCA1 letters {span[0]} to {span[-1]}: {letters}."
    return (
        f'<div class="seq-wrap"><div class="seq" role="img" aria-label="{esc(label)}">{"".join(cells)}</div>'
        f'<p class="seq-cap">BRCA1, positions {span[0]}–{span[-1]}</p></div>'
    )


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
        f'<p class="fig">Where {esc(mutant)} lands. Further left means more likely harmful.</p>'
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
    return (
        '<details class="numbers"><summary>'
        f'{icon("chev", "chev")}Show the numbers</summary>'
        f'<dl class="kv">{body}</dl>'
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
    where = "start region, where we tuned the cutoffs" if info["region"] == "start" else "end region, the held-out test"
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
        f'<p class="where">Position {position}, {where}</p></div>'
        f'{seq_html(position, new, vkind)}'
        f'<span class="verdict {vcss}">{icon(vkind)}{esc(verdict)}</span>'
        f'<p class="plain">{verdict_plain[verdict]}</p>'
        '<div class="agree">'
        f'<span class="agree-badge {agree["css"]}">{icon(agree["icon"])}{agree["label"]}</span>'
        f'<p class="agree-mean"><span class="q">What does this mean?</span> {agree["mean"]}</p></div>'
        f'<p class="story">{story}</p>'
        f'{gauge(info["esm_score"], lab_class, lab_score, info["mutant"])}'
        f"{numbers_panel(info, agree)}"
        f"{explorer_link(info['mutant'])}"
        "</article>"
    )


def explorer_link(mutant: str) -> str:
    if not explorer_url:
        return ""
    return (
        f'<a class="to-3d" href="{esc(explorer_url)}?m={esc(mutant)}" target="_blank" rel="noopener">'
        f"Open {esc(mutant)} in 3D</a>"
        '<p class="to-3d-note">See where this swap sits on the BRCA1 structure (orientation only; some positions are not in the crystal structure), with every other letter the model scored at this position.</p>'
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
        f'<div class="try"><span class="try-lead">Try one of these</span>{chips}</div></article>'
    )


def empty_html() -> str:
    probe = 61 if 61 in _WT_BY_POS else (scored_positions[6] if len(scored_positions) > 6 else None)
    strip = seq_html(probe) if probe is not None else ""
    return (
        f'<div class="empty">{strip}'
        "<p>Each letter is one amino acid, a building block of the protein. A mutation swaps one letter for another: "
        "C61G turns the C at position 61 into a G.</p>"
        "<p>Type a change above, or pick an example. We show what the model calls it and what the lab measured.</p></div>"
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


def header_html() -> str:
    return (
        '<header class="top">'
        "<h1>BRCA Spellcheck</h1>"
        "<p>Type one letter change in the BRCA1 protein to see the model's call, the lab's answer, "
        "and whether the two agree. BRCA1 is our exam, not the use case.</p>"
        "</header>"
    )


def input_hint_html() -> str:
    return (
        '<p class="hint">Old letter, position, new letter: C61G, c61g, Cys61Gly or p.Cys61Gly. '
        f"We scored positions {START_REGION[0]}–{START_REGION[1]} and {END_REGION[0]}–{END_REGION[1]}.</p>"
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
    return '<p class="chip-note">' + ", ".join(parts) + ". We left them in on purpose.</p>"


arm_titles: dict[str, str] = {
    "arm0": "Always says harmless",
    "arm1": "BLOSUM62 rule",
    "arm2": "ESM-2 (this demo)",
    "arm3": "ESM-2 plus extras",
}


def cell(count: int, total: int, css: str, k: int) -> str:
    share = count / total if total else 0.0
    return (
        f'<td class="num"><span class="cnt">{count}</span>'
        f'<span class="mini" aria-hidden="true"><i class="{css}" style="--v:{share:.4f};--k:{k}"></i></span></td>'
    )


def scoreboard_html() -> str:
    if headline is None:
        return '<div class="plain-box"><h2>Scoreboard not available.</h2><p>headline.csv is missing, so we show nothing rather than guess.</p></div>'
    first = headline.iloc[0]
    n_harmful = int(first["n_harmful"])
    n_harmless = int(first["n_harmless"])
    body = ""
    for k, (_, arm) in enumerate(headline.iterrows()):
        arm_id = str(arm["arm"])
        cut = f'{float(arm["cutoff"]):.6f}' if str(arm["cutoff"]) != "" else "none"
        mark = ' class="ours"' if arm_id == "arm2" else ""
        body += (
            f'<tr{mark}><th scope="row">{esc(arm_titles.get(arm_id, arm["arm_name"]))}'
            f'<span class="sub">{esc(arm["arm_name"])}, cutoff <span class="num">{cut}</span></span></th>'
            f'{cell(int(arm["harmful_caught"]), n_harmful, "caught", k)}'
            f'{cell(int(arm["harmful_missed"]), n_harmful, "missed", k)}'
            f'{cell(int(arm["false_alarms"]), n_harmless, "alarm", k)}'
            f'<td class="num rate">{float(arm["catch_rate_percent"]):.1f}%</td></tr>'
        )
    table = (
        '<div class="table-wrap"><table class="board">'
        "<thead><tr><th scope=\"col\">Method</th>"
        f'<th scope="col">Harmful caught<span class="sub">of {n_harmful}</span></th>'
        f'<th scope="col">Harmful missed<span class="sub">of {n_harmful}</span></th>'
        f'<th scope="col">False alarms<span class="sub">of {n_harmless}</span></th>'
        '<th scope="col">Catch rate</th></tr></thead>'
        f"<tbody>{body}</tbody></table></div>"
    )
    sentence = sentences.get("S6")
    sentence_html = (
        f'<p class="sentence">{esc(sentence)}</p><p class="small">Our headline result. We wrote this sentence before the test run and filled in the numbers from it.</p>'
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
        '<table class="rank"><caption>Ranking quality on the end region (AUROC, 0.5 is a coin flip)</caption>'
        f"<tbody>{rank_rows}</tbody></table>{boot}"
    )
    head = (
        "<h2>The end region the model never saw</h2>"
        f"<p class=\"intro\">We set every cutoff on positions {START_REGION[0]}–{START_REGION[1]}, froze them, then scored "
        f"positions {END_REGION[0]}–{END_REGION[1]} once. There the lab gives a clear answer for "
        f'<span class="num">{n_harmful}</span> harmful and <span class="num">{n_harmless}</span> harmless mutations. '
        "Misses and false alarms sit side by side.</p>"
    )
    return f'<section class="score">{sentence_html}{head}{table}{ranking}{failure_html()}</section>'


def failure_html() -> str:
    miss = failure_row("miss")
    if miss is None:
        return ""
    position = int(miss["position"])
    near = 0
    total_misses = 0
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
        '<div class="failure">'
        f'<h3>{icon("disagree", "st-disagree")}Worst miss: <span class="num">{esc(miss["mutant"])}</span></h3>'
        f'<p>The lab measured a function score of <span class="num">{signed(float(miss["function_score"]), 3)}</span>, '
        f'well into broken. The model scored it <span class="num">{signed(float(miss["esm_score"]), 2)}</span>, '
        f'an ordinary-looking swap, and said "{esc(miss["band_arm2"])}".{also}{cluster}</p>'
        '<p class="guess"><b>Our guess, not tested:</b> the lab edited DNA, not protein. This DNA change may upset '
        "how the gene's message is spliced, and a model that only reads protein letters cannot see that.</p>"
        '<p class="small">From failures.csv, rank 1 of the misses.</p>'
        "</div>"
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


def on_load(request: gr.Request):
    q = ""
    tab = "check"
    try:
        q = (request.query_params.get("q") or "").strip()[:40]
        tab = request.query_params.get("tab") or "check"
    except Exception:
        pass
    if tab not in ("check", "score"):
        tab = "check"
    shown = render_result(q) if q else empty_html()
    return q, shown, gr.Tabs(selected=tab)


font_link: str = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">'
    '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
    '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&'
    'family=JetBrains+Mono:wght@400;500&display=swap">'
    '<meta name="description" content="A prototype that compares a protein language model\'s call on BRCA1 '
    'letter changes with lab measurements.">'
)

head_script: str = """
<script>
(function () {
  function light() {
    if (document.body && document.body.classList.contains("dark")) document.body.classList.remove("dark");
  }
  function start() {
    light();
    var watch = new MutationObserver(light);
    watch.observe(document.body, { attributes: true, attributeFilter: ["class"] });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start); else start();
  var asked = 0;
  document.addEventListener("click", function (ev) {
    if (ev.target.closest && ev.target.closest("#check-btn, .chip, .try-chip")) asked = Date.now();
  }, true);
  document.addEventListener("keydown", function (ev) {
    if (ev.key === "Enter" && ev.target.closest && ev.target.closest("#mut-input")) asked = Date.now();
  }, true);
  new MutationObserver(function () {
    if (Date.now() - asked > 4000 || window.innerWidth > 760) return;
    var card = document.querySelector(".result-out .result-card, .result-out .refusal");
    if (!card) return;
    asked = 0;
    var calm = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    card.scrollIntoView({ behavior: calm ? "auto" : "smooth", block: "start" });
  }).observe(document.documentElement, { childList: true, subtree: true });
  document.addEventListener("click", function (ev) {
    var chip = ev.target.closest && ev.target.closest(".try-chip");
    if (!chip) return;
    var box = document.querySelector("#mut-input textarea, #mut-input input");
    var btn = document.querySelector("#check-btn");
    if (!box || !btn) return;
    var setter = Object.getOwnPropertyDescriptor(Object.getPrototypeOf(box), "value").set;
    setter.call(box, chip.dataset.q);
    box.dispatchEvent(new Event("input", { bubbles: true }));
    setTimeout(function () { btn.click(); }, 30);
    box.focus();
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
    body_background_fill="#f3f5f8",
    body_text_color="#151a21",
    block_background_fill="#ffffff",
    block_border_color="#dbe0e6",
    input_background_fill="#ffffff",
    button_primary_background_fill="#1f4a73",
    button_primary_text_color="#f7f6f2",
)

css_text: str = (APP_DIR / "style.css").read_text(encoding="utf-8") if (APP_DIR / "style.css").exists() else ""


with gr.Blocks(title="BRCA Spellcheck") as demo:
    gr.HTML(header_html())
    with gr.Tabs(selected="check", elem_classes=["tabs-main"]) as tabs:
        with gr.Tab("Check", id="check"):
            with gr.Row(elem_classes=["ask-row"]):
                mutation_input = gr.Textbox(
                    label="Mutation",
                    placeholder="C61G",
                    max_lines=1,
                    elem_id="mut-input",
                    elem_classes=["mut-input"],
                    scale=3,
                    min_width=200,
                )
                check_button = gr.Button(
                    "Check this swap", elem_id="check-btn", elem_classes=["btn-primary"], scale=1, min_width=160
                )
            gr.HTML(input_hint_html())
            with gr.Row(elem_classes=["chip-row"]):
                chip_buttons = [
                    gr.Button(q, size="sm", elem_classes=["chip"], min_width=0, scale=0) for q in examples
                ]
            gr.HTML(chip_note_html())
            result_output = gr.HTML(empty_html(), elem_classes=["result-out"])
        with gr.Tab("Scoreboard", id="score"):
            gr.HTML(scoreboard_html())
    gr.HTML(notes_html(), elem_classes=["notes-wrap"])

    check_button.click(fn=render_result, inputs=mutation_input, outputs=result_output)
    mutation_input.submit(fn=render_result, inputs=mutation_input, outputs=result_output)
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
