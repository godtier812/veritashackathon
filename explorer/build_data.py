import csv
import json
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parent
app_dir = root / "app"
results = root / "outputs" / "results"
app_results = app_dir / "results"


def pick(name):
    for folder in (results, app_results):
        path = folder / name
        if path.exists():
            return path
    return None


def read_rows(name):
    path = pick(name)
    if path is None:
        return None
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def num(text):
    if text is None or str(text).strip() == "":
        return None
    return float(text)


def read_config_constants():
    found = {}
    for line in (root / "src" / "config.py").read_text(encoding="utf-8").splitlines():
        for key in ("LOF_BELOW", "FUNC_ABOVE"):
            if line.startswith(key + ":"):
                found[key.lower()] = float(line.split("=", 1)[1].strip())
    return found


three = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I",
    "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S",
    "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def read_model_one(path):
    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            atoms.append({
                "name": line[12:16].strip(),
                "res": line[17:20].strip(),
                "chain": line[21],
                "resi": int(line[22:26]),
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "el": line[76:78].strip() or line[12:14].strip(),
            })
    return atoms


def check_structure(name, chain, positions):
    atoms = read_model_one(here / "structures" / f"{name}.pdb")
    letters = {}
    for a in atoms:
        if a["chain"] == chain and a["name"] == "CA" and a["res"] in three:
            letters[a["resi"]] = three[a["res"]]
    region = "start" if name == "1jm7" else "end"
    scored = [p for p in positions.values() if p["region"] == region]
    present = [p["p"] for p in scored if p["p"] in letters]
    match = [p for p in present if letters[p] == positions[p]["wt"]]
    missing = [p["p"] for p in scored if p["p"] not in letters]
    zinc = []
    if name == "1jm7":
        zn = [a for a in atoms if a["res"] == "ZN" and a["chain"] == chain]
        for a in atoms:
            if a["chain"] != chain:
                continue
            if (a["res"] == "CYS" and a["name"] == "SG") or (a["res"] == "HIS" and a["name"] in ("ND1", "NE2")):
                for z in zn:
                    d = sum((u - v) ** 2 for u, v in zip(a["xyz"], z["xyz"])) ** 0.5
                    if d <= 2.6:
                        zinc.append({"resi": a["resi"], "res": a["res"], "atom": a["name"],
                                     "zn": z["resi"], "dist": round(d, 2)})
    return {
        "id": name.upper(),
        "chain": chain,
        "n_scored": len(scored),
        "n_with_coords": len(present),
        "n_letter_match": len(match),
        "first": min(letters),
        "last": max(letters),
        "missing": missing,
        "zinc": zinc,
        "passed": len(present) > 0 and len(match) == len(present),
    }


def main():
    config = json.loads((app_dir / "app_config.json").read_text(encoding="utf-8"))
    consts = read_config_constants()

    with open(app_dir / "app_data.csv", newline="", encoding="utf-8") as handle:
        app_rows = list(csv.DictReader(handle))
    with open(root / "outputs" / "esm_scores.csv", newline="", encoding="utf-8") as handle:
        esm_rows = list(csv.DictReader(handle))

    entropy = {}
    esm_check = {}
    for row in esm_rows:
        pos = int(row["position"])
        entropy[pos] = float(row["entropy"])
        esm_check[(pos, row["new_letter"])] = float(row["esm_score"])

    blosum = {}
    blosum_rows = read_rows("blosum62.csv")
    if blosum_rows:
        for row in blosum_rows:
            for key, value in row.items():
                if key != "letter":
                    blosum[(row["letter"], key)] = int(value)

    positions = {}
    mismatch = 0
    for row in app_rows:
        pos = int(row["position"])
        new = row["new_letter"]
        esm = float(row["esm_score"])
        if abs(esm - esm_check.get((pos, new), 1e9)) > 1e-6:
            mismatch += 1
        entry = positions.setdefault(pos, {
            "p": pos,
            "wt": row["wild_type"],
            "region": "start" if pos <= 101 else "end",
            "entropy": entropy.get(pos),
            "swaps": [],
        })
        entry["swaps"].append([
            new,
            esm,
            row["lab_class"],
            num(row["lab_score"]),
            blosum.get((row["wild_type"], new)),
        ])

    headline = read_rows("headline.csv") or []
    arms = []
    for row in headline:
        arms.append({
            "arm": row["arm"],
            "name": row["arm_name"],
            "cutoff": num(row["cutoff"]),
            "n_harmful": int(row["n_harmful"]),
            "n_harmless": int(row["n_harmless"]),
            "caught": int(row["harmful_caught"]),
            "missed": int(row["harmful_missed"]),
            "false_alarms": int(row["false_alarms"]),
            "cleared": int(row["harmless_cleared"]),
            "catch_rate": float(row["catch_rate_percent"]),
            "false_alarm_rate": float(row["false_alarm_rate_percent"]),
        })

    auroc = {r["arm"]: float(r["auroc"]) for r in (read_rows("auroc.csv") or [])}
    bands = {r["findlay_class"]: {k: int(v) for k, v in r.items() if k != "findlay_class"}
             for r in (read_rows("band_counts.csv") or [])}

    sentences = {}
    path = pick("result_sentences.txt")
    if path:
        for line in path.read_text(encoding="utf-8").splitlines():
            if ":" in line:
                key, text = line.split(":", 1)
                sentences[key.strip()] = text.strip()

    failures = []
    for row in read_rows("failures.csv") or []:
        failures.append({
            "kind": row["kind"],
            "rank": int(row["rank"]),
            "mutant": row["mutant"],
            "position": int(row["position"]),
            "lab_score": float(row["function_score"]),
            "lab_class": row["findlay_class"],
            "esm_score": float(row["esm_score"]),
            "guess": row["guess_why"],
        })

    structures = {
        "start": check_structure("1jm7", "A", positions),
        "end": check_structure("1t15", "A", positions),
    }
    assumed = [24, 27, 39, 41, 44, 47, 61, 64]
    derived = sorted({z["resi"] for z in structures["start"]["zinc"]})
    structures["start"]["zinc_positions"] = derived
    structures["start"]["zinc_matches_assumed"] = derived == assumed

    data = {
        "structures": structures,
        "config": {
            "cutoff_a": config["cutoff_a"],
            "cutoff_b": config.get("cutoff_b"),
            "has_band": config.get("has_band", False),
            "end_labels_included": config.get("end_labels_included", False),
            "protein_length": config.get("protein_length"),
            "model_id": config.get("model_id"),
            "n_tested_positions": config.get("n_tested_positions"),
            "start": [1, 101],
            "end": [1631, 1855],
            "lof_below": consts.get("lof_below"),
            "func_above": consts.get("func_above"),
        },
        "positions": [positions[k] for k in sorted(positions)],
        "arms": arms,
        "auroc": auroc,
        "bands": bands,
        "sentences": sentences,
        "failures": failures,
    }
    (here / "data.json").write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")

    ok = True
    c61 = positions[61]
    c61g = [s for s in c61["swaps"] if s[0] == "G"][0]
    line = f"c61g: wt={c61['wt']} esm={c61g[1]:.2f} lab={c61g[3]:.3f} class={c61g[2]}"
    good = c61["wt"] == "C" and f"{c61g[1]:.2f}" == "-10.50" and f"{c61g[3]:.3f}" == "-1.738"
    print(line, "ok" if good else "FAIL")
    ok = ok and good

    with open(pick("headline.csv"), newline="", encoding="utf-8") as handle:
        raw = list(csv.DictReader(handle))
    for row, arm in zip(raw, arms):
        same = (int(row["harmful_caught"]) == arm["caught"]
                and int(row["harmful_missed"]) == arm["missed"]
                and int(row["false_alarms"]) == arm["false_alarms"])
        print(f"{arm['arm']}: caught={arm['caught']} missed={arm['missed']} "
              f"false_alarms={arm['false_alarms']} headline={'match' if same else 'FAIL'}")
        ok = ok and same

    n_swaps = sum(len(p["swaps"]) for p in data["positions"])
    n_lab = sum(1 for p in data["positions"] for s in p["swaps"] if s[2])
    print(f"positions={len(positions)} swaps={n_swaps} lab_tested={n_lab} "
          f"esm_vs_app_mismatch={mismatch} entropy_missing={sum(1 for p in positions.values() if p['entropy'] is None)}")
    ok = ok and mismatch == 0

    for key, s in structures.items():
        print(f"{s['id']} chain {s['chain']}: residues {s['first']}-{s['last']}, "
              f"scored {s['n_scored']}, with coords {s['n_with_coords']}, "
              f"letter match {s['n_letter_match']} of {s['n_with_coords']} -> {'PASS' if s['passed'] else 'FAIL'}")
        print(f"  scored positions without coords: {s['missing']}")
        ok = ok and s["passed"]
    for z in structures["start"]["zinc"]:
        print(f"  zinc ligand: {z['res']}{z['resi']} {z['atom']} to ZN {z['zn']} at {z['dist']} A")
    print(f"zinc positions from file: {derived}; assumed: {assumed}; "
          f"{'same' if derived == assumed else 'DIFFERENT'}")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
