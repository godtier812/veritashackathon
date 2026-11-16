import argparse
import json
import shutil
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent
root = here.parent
app_dir = root / "app"
na = "not available"
disclaimer = "This project is a prototype and is not medical advice."
demo_note = (
    "At lab-tested positions a counselor would use the lab result. "
    "This demo shows how the model does where we can check it."
)


def find(name):
    for folder in (root / "outputs" / "results", app_dir / "results"):
        if (folder / name).exists():
            return folder / name
    return None


def table(name):
    path = find(name)
    return pd.read_csv(path, keep_default_na=False) if path else None


def esc(text):
    out = str(text)
    for old, new in (("&", "&amp;"), ("<", "&lt;"), (">", "&gt;"), ('"', "&quot;")):
        out = out.replace(old, new)
    return out


def signed(value, places):
    text = f"{value:.{places}f}"
    if text.startswith("-"):
        return "−" + text[1:]
    return ("+" + text) if float(text) > 0 else text


config = json.loads((app_dir / "app_config.json").read_text(encoding="utf-8"))
data = pd.read_csv(app_dir / "app_data.csv", keep_default_na=False)
headline = table("headline.csv")
auroc = table("auroc.csv")
failures = table("failures.csv")
boot_path = find("bootstrap.json")
bootstrap = json.loads(boot_path.read_text(encoding="utf-8")) if boot_path else None
sentences = {}
sent_path = find("result_sentences.txt")
if sent_path:
    for line in sent_path.read_text(encoding="utf-8").splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            sentences[k.strip()] = v.strip()

start_region = (1, 101)
end_region = (1631, 1855)
protein_length = int(config.get("protein_length", 1863))
wt = dict(zip(data["position"].astype(int), data["wild_type"]))
start_positions = sorted(p for p in wt if p <= start_region[1])
harm_by_pos = (-data.groupby("position")["esm_score"].mean()).to_dict()
lab = data[data["lab_class"] != ""]
n_lab = len(lab)
start_lab = lab[lab["position"] <= start_region[1]]
start_lof = int((start_lab["lab_class"] == "LOF").sum())
start_func = int((start_lab["lab_class"] == "FUNC").sum())
row61 = data[(data["position"] == 61) & (data["new_letter"] == "G")]
c61g_esm = float(row61["esm_score"].iloc[0]) if len(row61) else None
c61g_lab = float(row61["lab_score"].iloc[0]) if len(row61) and row61["lab_score"].iloc[0] != "" else None


def arm(arm_id):
    if headline is None:
        return None
    rows = headline[headline["arm"] == arm_id]
    return rows.iloc[0].to_dict() if len(rows) else None


def fail(kind):
    if failures is None:
        return None
    rows = failures[(failures["kind"] == kind) & (failures["rank"].astype(str) == "1")]
    return rows.iloc[0].to_dict() if len(rows) else None


def hero_seq():
    order = sorted(start_positions, key=lambda p: harm_by_pos[p])
    rank = {p: k / (len(order) - 1) for k, p in enumerate(order)} if len(order) > 1 else {}
    cells = []
    for i, p in enumerate(start_positions):
        h = rank.get(p, 0) ** 2.2
        letter = esc(wt[p])
        beat = round(rank.get(p, 0) * 100)
        tip = f"Position {p} ({wt[p]}). Changes here look more harmful than at {beat}% of these positions."
        if p == 61:
            cells.append(
                f'<span class="lt hit" style="--i:{i};--h:{h:.3f}" data-tip="{esc(tip)}">'
                f'<span class="old">{letter}</span><span class="new">G</span><i></i></span>'
            )
        else:
            cells.append(f'<span class="lt" style="--i:{i};--h:{h:.3f}" data-tip="{esc(tip)}">{letter}<i></i></span>')
    letters = "".join(wt[p] for p in start_positions)
    return (
        f'<div class="seq" role="img" aria-label="BRCA1 positions {start_positions[0]} to {start_positions[-1]}: {letters}. '
        f'The cysteine at position 61 is replaced by glycine.">{"".join(cells)}</div>'
    )


def strip(position, new, kind):
    span = [p for p in range(position - 5, position + 6) if p in wt]
    cells = []
    for i, p in enumerate(span):
        if p == position:
            cells.append(
                f'<span class="aa hit k-{kind}" style="--i:{i}"><span class="old">{esc(wt[p])}</span>'
                f'<span class="new">{esc(new)}</span><i></i><span class="pos">{p}</span></span>'
            )
        else:
            cells.append(f'<span class="aa" style="--i:{i}">{esc(wt[p])}</span>')
    return f'<div class="strip" data-anim aria-label="BRCA1 positions {span[0]} to {span[-1]}, {esc(wt[position])}{position}{esc(new)}" role="img">{"".join(cells)}</div>'


def dots(total, marked, css_marked, css_rest, label):
    items = []
    for i in range(total):
        css = css_marked if i < marked else css_rest
        items.append(f'<i class="{css}" style="--i:{i}"></i>')
    return f'<div class="dots" data-anim role="img" aria-label="{esc(label)}">{"".join(items)}</div>'


def results_block():
    a2 = arm("arm2")
    if a2 is None:
        return f"<p>{na}</p>"
    n_harm = int(a2["n_harmful"])
    n_ok = int(a2["n_harmless"])
    caught = int(a2["harmful_caught"])
    missed = int(a2["harmful_missed"])
    alarms = int(a2["false_alarms"])
    cleared = int(a2["harmless_cleared"])
    left = (
        '<figure class="field"><figcaption>'
        f'<b>{n_harm}</b> mutations the lab found harmful</figcaption>'
        + dots(n_harm, caught, "d-caught", "d-missed", f"{caught} caught, {missed} missed")
        + f'<ul class="key"><li><span class="sw d-caught"></span><span class="tick" data-to="{caught}" data-of="{caught}">{caught}</span> caught</li>'
        f'<li><span class="sw d-missed"></span><span class="tick" data-to="{missed}" data-from="{caught}" data-of="{missed}">{missed}</span> missed</li></ul></figure>'
    )
    right = (
        '<figure class="field"><figcaption>'
        f'<b>{n_ok}</b> mutations the lab found harmless</figcaption>'
        + dots(n_ok, alarms, "d-alarm", "d-cleared", f"{alarms} false alarms, {cleared} cleared")
        + f'<ul class="key"><li><span class="sw d-alarm"></span><span class="tick" data-to="{alarms}" data-of="{alarms}">{alarms}</span> false alarms</li>'
        f'<li><span class="sw d-cleared"></span><span class="tick" data-to="{cleared}" data-from="{alarms}" data-of="{cleared}">{cleared}</span> correctly left alone</li></ul></figure>'
    )
    return f'<div class="fields">{left}{right}</div>'


arm_names = {
    "arm0": "Always says harmless",
    "arm1": "BLOSUM62 rule",
    "arm2": "ESM-2, our model",
    "arm3": "ESM-2 plus extras",
}


def arms_table():
    if headline is None:
        return f"<p>{na}</p>"
    auc = {}
    if auroc is not None:
        auc = {str(r["arm"]): float(r["auroc"]) for _, r in auroc.iterrows()}
    rows = ""
    for _, r in headline.iterrows():
        aid = str(r["arm"])
        mark = ' class="ours"' if aid == "arm2" else ""
        a = f"{auc[aid]:.3f}" if aid in auc else na
        rows += (
            f'<tr{mark}><th scope="row">{esc(arm_names.get(aid, r["arm_name"]))}</th>'
            f'<td>{int(r["harmful_caught"])}</td><td>{int(r["harmful_missed"])}</td>'
            f'<td>{int(r["false_alarms"])}</td><td>{float(r["catch_rate_percent"]):.1f}%</td><td>{a}</td></tr>'
        )
    first = headline.iloc[0]
    return (
        '<div class="tw"><table class="arms"><thead><tr><th scope="col">Method</th>'
        f'<th scope="col">Caught <small>of {int(first["n_harmful"])}</small></th>'
        f'<th scope="col">Missed <small>of {int(first["n_harmful"])}</small></th>'
        f'<th scope="col">False alarms <small>of {int(first["n_harmless"])}</small></th>'
        '<th scope="col">Catch rate</th><th scope="col">AUROC</th></tr></thead>'
        f"<tbody>{rows}</tbody></table></div>"
    )


def boot_line():
    if not bootstrap or "esm_minus_blosum" not in bootstrap:
        return ""
    b = bootstrap["esm_minus_blosum"]
    return (
        f'<p class="fine">AUROC ranks how well each method orders harmful above harmless; 0.5 is a coin flip. '
        f'ESM-2 minus BLOSUM62: {signed(b["point"], 3)}, with a 2.5th to 97.5th percentile range of '
        f'{signed(b["low"], 3)} to {signed(b["high"], 3)} over {int(bootstrap.get("n_resamples_used", 0)):,} '
        f'resamples of {esc(bootstrap.get("unit", "positions"))}.</p>'
    )


def failure_block():
    miss = fail("miss")
    alarm = fail("false_alarm")
    parts = []
    if miss:
        parts.append(
            '<article class="case">'
            + strip(int(miss["position"]), miss["new_letter"], "harmless")
            + f'<h3>{esc(miss["mutant"])}, our worst miss</h3>'
            f'<p>The lab measured {signed(float(miss["function_score"]), 3)}, well into broken. '
            f'The model scored it {signed(float(miss["esm_score"]), 2)} and called it "{esc(miss["band_arm2"])}". '
            + ("The BLOSUM62 rule and ESM-2 plus extras missed it too." if str(miss.get("flag_arm1")) == "0" and str(miss.get("flag_arm3")) == "0" else "")
            + "</p>"
            '<p class="guess"><b>Our guess, not tested:</b> the lab edited DNA, not protein. This change may upset '
            "how the gene's message is spliced, which a model that only reads protein letters cannot see.</p></article>"
        )
    if alarm:
        parts.append(
            '<article class="case">'
            + strip(int(alarm["position"]), alarm["new_letter"], "harmful")
            + f'<h3>{esc(alarm["mutant"])}, our worst false alarm</h3>'
            f'<p>The model scored it {signed(float(alarm["esm_score"]), 2)} and called it "{esc(alarm["band_arm2"])}". '
            f'In the lab the protein kept working (score {signed(float(alarm["function_score"]), 3)}).</p>'
            + ("<p>The BLOSUM62 rule and ESM-2 plus extras flagged it too.</p>" if str(alarm.get("flag_arm1")) == "1" and str(alarm.get("flag_arm3")) == "1" else "")
            + '<p class="guess"><b>Our guess, not tested:</b> a proline swap usually bends a protein out of shape, which is why every method '
            "flagged it. But this position varies a lot across related proteins, and here the lab found the protein tolerated the change.</p></article>"
        )
    return '<div class="cases">' + "".join(parts) + "</div>"


def jfile(name):
    path = find(name)
    return json.loads(path.read_text(encoding="utf-8")) if path else None


def card(big, title, text, tag=""):
    label = f'<span class="win-tag">{esc(tag)}</span>' if tag else ""
    return f'<article class="win">{label}<p class="win-big">{big}</p><h3>{esc(title)}</h3><p class="win-fine">{text}</p></article>'


def wins_block():
    cards = []
    a2 = arm("arm2")
    auc = {}
    if auroc is not None:
        auc = {str(r["arm"]): float(r["auroc"]) for _, r in auroc.iterrows()}
    bands = table("band_counts.csv")
    b = {}
    if bands is not None:
        b = {str(r["findlay_class"]): r for _, r in bands.iterrows()}
    if a2 is not None:
        split = ""
        if "LOF" in b and "FUNC" in b:
            split = (
                f" ({int(b['LOF']['looks_harmful'])} \u201clooks harmful\u201d, {int(b['LOF']['cant_tell'])} \u201ccan\u2019t tell\u201d)"
            )
        alarm_split = f" ({int(b['FUNC']['cant_tell'])} of them \u201ccan\u2019t tell\u201d)" if "FUNC" in b else ""
        cards.append(card(
            f"{float(a2['catch_rate_percent']):.1f}%",
            "of harmful mutations caught, on a part of the protein we kept hidden until one final test.",
            f"{int(a2['harmful_caught'])} of {int(a2['n_harmful'])} not cleared as harmless{split}, positions "
            f"{end_region[0]}\u2013{end_region[1]}, tested once. It also failed to clear {int(a2['false_alarms'])} of "
            f"{int(a2['n_harmless'])} harmless ones{alarm_split}.",
        ))
    if b:
        if "LOF" in b and "FUNC" in b:
            right = int(b["FUNC"]["looks_harmless"])
            wrong = int(b["LOF"]["looks_harmless"])
            n_func = int(b["FUNC"]["total"])
            n_clear = n_func + int(b["LOF"]["total"])
            cards.append(card(
                f"{100.0 * right / (right + wrong):.1f}%",
                "of the time, when the model said a change looks harmless, the lab agreed.",
                f"{right} of {right + wrong} \u201clooks harmless\u201d calls on the held-out test with a clear lab "
                f"answer. It made that call for {right + wrong} of {n_clear:,} changes; "
                f"{100.0 * n_func / n_clear:.0f}% of test changes are harmless. When unsure, it says "
                "\u201ccan\u2019t tell\u201d instead of guessing.",
            ))
    if "arm2" in auc:
        extra = f" A classic substitution table (BLOSUM62) scores {auc['arm1'] * 100:.0f}% on the same measure." if "arm1" in auc else ""
        cards.append(card(
            f"{auc['arm2'] * 100:.0f}%",
            "of the time, shown one harmful and one harmless change, the model ranks the harmful one as worse.",
            f"AUROC {auc['arm2']:.4f} on the held-out test; 50% would be a coin flip, and ties count as half.{extra}",
        ))
    cards.append(card(
        "0",
        "lab results the model learned from. It only read natural protein sequences.",
        "Zero-shot: lab answers were used only to set the cutoffs on the start region and to grade the model.",
    ))
    am = jfile("reference_alphamissense.json")
    if am and "end" in am:
        d = am["alphamissense_minus_esm"]
        cards.append(card(
            f"{am['end']['auroc_arm2_same_rows'] * 100:.0f}% <small>vs</small> {am['end']['auroc'] * 100:.0f}%",
            "close to AlphaMissense\u2019s ranking, just below it, using fewer kinds of data.",
            f"AUROC {am['end']['auroc_arm2_same_rows']:.4f} vs {am['end']['auroc']:.4f} on the same held-out rows; "
            f"the gap ({signed(d['point'], 3)}, 95% range {signed(d['low'], 3)} to {signed(d['high'], 3)}) includes zero. "
            "Our model reads only the protein sequence; AlphaMissense also uses 3D structure and population data. We do not claim to be better.",
            "checked after the test",
        ))
    if not cards:
        return f"<p>{na}</p>"
    return '<div class="wins">' + "".join(cards) + "</div>"


def protein_map():
    def pct(p):
        return (p - 1) / (protein_length - 1) * 100

    s0, s1 = pct(start_region[0]), pct(start_region[1])
    e0, e1 = pct(end_region[0]), pct(end_region[1])
    return (
        '<figure class="pmap" data-anim role="img" aria-label="BRCA1 is 1863 letters long. '
        f'We tuned on positions {start_region[0]} to {start_region[1]} and tested once on {end_region[0]} to {end_region[1]}.">'
        '<div class="rail"><i class="body"></i>'
        f'<i class="reg tune" style="left:{s0:.2f}%;width:{max(s1 - s0, 0.8):.2f}%"></i>'
        f'<i class="reg test" style="left:{e0:.2f}%;width:{e1 - e0:.2f}%"></i></div>'
        '<div class="ticks"><span>1</span>'
        f'<span>{protein_length:,}</span></div>'
        '<div class="tags">'
        f'<div class="tag tune"><b>Tuned here</b>Positions {start_region[0]}–{start_region[1]}, the RING region. '
        f'{start_lof} harmful and {start_func} harmless in the lab.</div>'
        f'<div class="tag test"><b>Tested once</b>Positions {end_region[0]}–{end_region[1]}, the BRCT region. '
        f'{int(arm("arm2")["n_harmful"]) if arm("arm2") else na} harmful and {int(arm("arm2")["n_harmless"]) if arm("arm2") else na} harmless in the lab.</div>'
        "</div></figure>"
    )


step_svgs = [
    '<svg viewBox="0 0 220 64" aria-hidden="true">'
    + "".join(
        f'<g class="pc" style="--j:{i}"><rect x="{4 + i * 30}" y="16" width="26" height="30" rx="5" fill="{"var(--accent-soft)" if c == "?" else "none"}" stroke="currentColor" stroke-width="1.25"{" stroke-dasharray=&quot;3 3&quot;" if c == "?" else ""}/>'
        f'<text x="{17 + i * 30}" y="36" text-anchor="middle">{c}</text></g>'
        for i, c in enumerate("".join(wt.get(p, "?") for p in range(58, 61)) + "?" + "".join(wt.get(p, "?") for p in range(62, 65)))
    )
    + "</svg>",
    '<svg viewBox="0 0 220 64" aria-hidden="true">'
    '<line class="pc" style="--j:0" x1="10" y1="54" x2="210" y2="54" stroke="currentColor"/>'
    '<rect class="pc grow" style="--j:1" x="40" y="14" width="40" height="40" rx="2" fill="var(--harmless-soft)" stroke="var(--harmless)" stroke-width="1.25"/>'
    '<rect class="pc grow" style="--j:3" x="130" y="46" width="40" height="8" rx="2" fill="var(--harmful-soft)" stroke="var(--harmful)" stroke-width="1.25"/>'
    '<text class="pc" style="--j:2" x="60" y="9" text-anchor="middle">C</text>'
    '<text class="pc" style="--j:4" x="150" y="40" text-anchor="middle">G</text>'
    "</svg>",
    '<svg viewBox="0 0 220 64" aria-hidden="true">'
    '<g class="pc wipe" style="--j:0"><rect x="6" y="26" width="80" height="10" fill="var(--harmful-soft)"/>'
    '<rect x="86" y="26" width="40" height="10" fill="var(--cant-soft)"/>'
    '<rect x="126" y="26" width="88" height="10" fill="var(--harmless-soft)"/></g>'
    '<g class="pc drop" style="--j:4"><line x1="86" y1="16" x2="86" y2="46" stroke="currentColor" stroke-dasharray="3 3"/>'
    '<text x="86" y="60" text-anchor="middle">B</text></g>'
    '<g class="pc drop" style="--j:6"><line x1="126" y1="16" x2="126" y2="46" stroke="currentColor" stroke-dasharray="3 3"/>'
    '<text x="126" y="60" text-anchor="middle">A</text></g>'
    "</svg>",
    '<svg viewBox="0 0 220 64" aria-hidden="true">'
    '<g class="pc" style="--j:0"><rect x="6" y="20" width="60" height="24" rx="5" fill="none" stroke="currentColor" stroke-width="1.25"/>'
    f'<text x="36" y="36" text-anchor="middle">{start_region[0]}–{start_region[1]}</text></g>'
    '<g class="pc wipe" style="--j:2"><path d="M72 32h50M116 27l7 5-7 5" stroke="currentColor" stroke-width="1.25" fill="none"/></g>'
    '<g class="pc" style="--j:5"><rect x="128" y="20" width="86" height="24" rx="5" fill="var(--accent-soft)" stroke="var(--accent)" stroke-width="1.25"/>'
    f'<text x="171" y="36" text-anchor="middle">{end_region[0]}–{end_region[1]}</text></g>'
    "</svg>",
]


def steps_block():
    example = f" C61G scores {signed(c61g_esm, 2)}." if c61g_esm is not None else ""
    steps = [
        ("Hide one letter",
         "We give ESM-2, a protein language model, a 500-letter stretch of BRCA1 with one position blanked out, "
         "and ask how likely each of the 20 amino acids is to belong there."),
        ("Compare the old and new letter",
         "The score is log P(new) minus log P(old). Below zero means the new letter surprises the model more than "
         "the one nature put there." + example),
        ("Draw two lines on the start region",
         f"On positions {start_region[0]}–{start_region[1]} we set one line to still catch at least 90% of the "
         f"{start_lof} harmful mutations, and another to flag at most 10% of the {start_func} harmless ones. "
         "Between the lines the model says it can't tell."),
        ("Grade once on the end region",
         f"Then we froze both lines and scored positions {end_region[0]}–{end_region[1]} one time against the lab. "
         "The results below are that single run, mistakes included."),
    ]
    items = "".join(
        f'<li data-anim><span class="n">{i + 1}</span><div><h3>{t}</h3><p>{p}</p></div><div class="dia">{step_svgs[i]}</div></li>'
        for i, (t, p) in enumerate(steps)
    )
    return f'<ol class="steps">{items}</ol>'


def credits_block():
    def link(text, url):
        return f'<a href="{url}" target="_blank" rel="noopener">{text}</a>'

    groups = [
        ("Data", [
            "Findlay GM et al. (2018), Nature 562:217–222, saturation genome editing of BRCA1. " + link("doi", "https://doi.org/10.1038/s41586-018-0461-z"),
            "ProteinGym (Notin P et al., 2023), source of the BRCA1 table we used. " + link("github", "https://github.com/OATML-Markslab/ProteinGym"),
        ]),
        ("Models and methods", [
            f"ESM-2 (Lin Z et al., 2023, Science 379:1123–1130), {esc(config.get('model_id', na))}.",
            "Masked-marginal scoring: Meier J et al. (2021), NeurIPS.",
            "BLOSUM62: Henikoff S, Henikoff JG (1992), PNAS 89:10915–10919.",
        ]),
        ("Prior work", [
            "ProteinGym already scores ESM-2 on this table.",
            "AlphaMissense (Cheng J et al., 2023, Science 381:eadg7492) is a stronger, freely available predictor.",
        ]),
        ("AI tools", [
            "Claude Code (Anthropic): planning, review, and the demo and page design.",
            "Gemini in Google Antigravity: writing and running the project code from our spec.",
        ]),
    ]
    return "".join(
        f'<div><h3>{g}</h3><ul>' + "".join(f"<li>{i}</li>" for i in items) + "</ul></div>"
        for g, items in groups
    )


def build(demo_url, explorer_url=""):
    a2 = arm("arm2")
    a1 = arm("arm1")
    s6 = sentences.get("S6", na)
    s1 = sentences.get("S1", "")
    auc = {}
    if auroc is not None:
        auc = {str(r["arm"]): float(r["auroc"]) for _, r in auroc.iterrows()}
    s1_line = ""
    if s1 and "arm2" in auc and "arm1" in auc:
        s1_line = (
            f'<p class="lead2">{esc(s1)} Its AUROC was {auc["arm2"]:.3f} against {auc["arm1"]:.3f} for BLOSUM62.</p>'
        )
    cta = (
        f'<a class="btn" href="{esc(demo_url)}">Try the demo</a>'
        + (f' <a class="btn btn-alt" href="{esc(explorer_url)}">Explore in 3D</a>' if explorer_url else "")
        if demo_url
        else ""
    )
    lab_line = (
        f"In the lab C61G scored {signed(c61g_lab, 3)}, well into broken. The model's score was {signed(c61g_esm, 2)}."
        if c61g_lab is not None and c61g_esm is not None
        else ""
    )
    html = TEMPLATE
    swaps = {
        "{{cta}}": cta,
        "{{hero_seq}}": hero_seq(),
        "{{protein_length}}": f"{protein_length:,}",
        "{{n_lab}}": f"{n_lab:,}",
        "{{c61g_line}}": lab_line,
        "{{steps}}": steps_block(),
        "{{wins}}": wins_block(),
        "{{pmap}}": protein_map(),
        "{{s6}}": esc(s6),
        "{{fields}}": results_block(),
        "{{s1_line}}": s1_line,
        "{{arms}}": arms_table(),
        "{{boot}}": boot_line(),
        "{{cases}}": failure_block(),
        "{{credits}}": credits_block(),
        "{{disclaimer}}": esc(disclaimer),
        "{{note}}": esc(demo_note),
        "{{css}}": (here / "style.css").read_text(encoding="utf-8"),
        "{{js}}": (here / "motion.js").read_text(encoding="utf-8"),
    }
    for k, v in swaps.items():
        html = html.replace(k, v)
    heat = root / "landing_parked" / "heatmap.png"
    if heat.exists():
        shutil.copyfile(heat, here / "heatmap.png")
        html = html.replace("{{heat}}", (
            '<figure class="heat" data-anim><img src="heatmap.png" alt="Heat map of ESM-2 scores for every swap at the lab-tested positions, '
            'start region above and end region below." loading="lazy">'
            "<figcaption>Columns are positions and rows are the 20 possible new letters. Red means the model found the swap surprising, "
            "blue means it did not. Black dots mark swaps the lab found broken.</figcaption></figure>"
        ))
    else:
        html = html.replace("{{heat}}", f"<p>{na}</p>")
    (here / "index.html").write_text(html, encoding="utf-8", newline="\n")
    return here / "index.html"


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BRCA Spellcheck</title>
<meta name="description" content="A prototype that checks a protein language model's calls on BRCA1 letter changes against lab measurements.">
<script>document.documentElement.classList.add("js")</script>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=JetBrains+Mono:wght@400;500&display=swap">
<style>{{css}}</style>
</head>
<body>
<a class="skip" href="#main">Skip to content</a>
<header class="bar">
  <div class="bar-in">
    <a class="name" href="#top">BRCA Spellcheck</a>
    <nav aria-label="Sections">
      <a href="#wins">Highlights</a>
      <a href="#how">How it works</a>
      <a href="#results">Results</a>
      <a href="#fails">Where it fails</a>
      <a href="#limits">Limits</a>
    </nav>
    {{cta}}
  </div>
</header>
<main id="main">
<section class="hero" id="top">
  <h1>A spellchecker for proteins, graded against the lab.</h1>
  <p class="lead">BRCA1 is written in {{protein_length}} letters, one per amino acid. We asked a protein language model which single-letter changes look like typos, then checked its calls against {{n_lab}} lab measurements.</p>
  {{hero_seq}}
  <p class="cap">The first 101 letters of BRCA1. The red underline marks how harmful the model finds changes at each position, on average across all 19 possible swaps. Darker means more harmful, compared with the other positions here. Then C61G: the C at position 61 becomes a G. {{c61g_line}}</p>
  <div class="hero-cta">{{cta}}<a class="link" href="#results">See the results</a></div>
  <p class="exam">BRCA1 is our exam, not the use case. It is one of the few genes with a lab answer key, so we can check the model here before anyone trusts this kind of model on a gene without one.</p>
</section>

<section class="sec" id="wins">
  <h2>What we can show.</h2>
  <p class="sub">Every number is read from our results files, from the one held-out test unless tagged.</p>
  {{wins}}
</section>

<section class="sec" id="how">
  <h2>Four steps, no training.</h2>
  <p class="sub">The model never saw a lab label. It only learned which letters tend to sit where in natural proteins. We borrowed the lab's answers to grade it.</p>
  {{steps}}
</section>

<section class="sec" id="exam">
  <h2>Tune on one end. Test on the other, once.</h2>
  <p class="sub">We never split mutations at random. Every cutoff was set at the start of the protein, frozen, and then used one time on the far end.</p>
  {{pmap}}
</section>

<section class="sec" id="results">
  <h2 class="big">{{s6}}</h2>
  <p class="sub">That is the held-out end region, scored once. Each dot is one mutation the lab gave a clear answer for.</p>
  {{fields}}
  {{s1_line}}
  {{arms}}
  {{boot}}
</section>

<section class="sec" id="fails">
  <h2>Where it fails.</h2>
  <p class="sub">We show our two worst mistakes on purpose. Both are guesses about why, and we did not test them.</p>
  {{cases}}
</section>

<section class="sec" id="map">
  <h2>Every swap we scored.</h2>
  {{heat}}
</section>

<section class="sec" id="limits">
  <h2>What it can't do.</h2>
  <ol class="limits">
    <li><h3>Blind to splicing</h3><p>The lab edited DNA. Some DNA changes break how the gene's message is spliced, which a model that reads only protein letters cannot see.</p></li>
    <li><h3>Cells in a dish, not people</h3><p>The answer key comes from cells grown in a lab. A 2025 preprint notes that results can depend on cell type.</p></li>
    <li><h3>Stronger tools exist</h3><p>Zero-shot scoring like ours is not a new method, and AlphaMissense is a stronger predictor. Our point is to show how far a sequence-only model can be trusted, and where it breaks.</p></li>
  </ol>
  <div class="credits">{{credits}}</div>
</section>
</main>
<footer class="notes" role="note"><div class="notes-in"><b>{{disclaimer}}</b><span>{{note}}</span></div></footer>
<script>{{js}}</script>
</body>
</html>
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo-url", default="http://127.0.0.1:7861/tool/")
    parser.add_argument("--explorer-url", default="http://127.0.0.1:7861/explorer/")
    args = parser.parse_args()
    print(build(args.demo_url, args.explorer_url))
