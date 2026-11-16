import json
import re
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here.parent / "app"))

import app

data = json.loads((here / "data.json").read_text(encoding="utf-8"))
script = """
const rules = require(process.argv[1]);
const data = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
const out = [];
for (const p of data.positions) for (const s of p.swaps) if (s[2]) {
  const v = rules.verdict(s[1], data.config);
  out.push([p.wt + p.p + s[0], v, rules.agreement(v, s[2])]);
}
const parsed = ['C61G', 'c61g', 'p.Cys61Gly', 'xyz', 'C61C', 'A9999G', ''].map(t => [t, rules.parse(t, 1863)]);
console.log(JSON.stringify({ out, parsed }));
"""
result = subprocess.run(
    ["node", "-e", script, str(here / "verdict.js"), str(here / "data.json")],
    capture_output=True, text=True, check=True,
)
js = json.loads(result.stdout)

same = 0
same_agree = 0
diff = []
for mutant, verdict, agree in js["out"]:
    text = app.check_mutation(mutant)
    found = re.search(r"\*\*Model verdict:\*\* (.+)", text)
    app_verdict = found.group(1).strip() if found else None
    info = app.lookup(mutant[0], int(mutant[1:-1]), mutant[-1])
    app_agree = app.agreement(info["verdict"], info["lab_class"])["key"]
    if app_verdict == verdict:
        same += 1
    else:
        diff.append((mutant, verdict, app_verdict))
    if app_agree == agree:
        same_agree += 1

total = len(js["out"])
print(f"verdict parity: {same} of {total} lab-tested swaps agree ({100.0 * same / total:.1f}%)")
print(f"agreement-state parity: {same_agree} of {total}")
for text, parsed in js["parsed"]:
    try:
        py = app.parse_input(text)
        py_text = f"{py[0]}{py[1]}{py[2]}"
    except ValueError as exc:
        py_text = str(exc)
    js_text = parsed.get("error") or f"{parsed['wt']}{parsed['p']}{parsed['nw']}"
    print(f"parse {text!r}: {'match' if py_text == js_text else 'DIFF'} -> {js_text[:60]}")
if diff or same_agree != total:
    print("differences:", diff[:10])
    sys.exit(1)
