---
name: index-debugger
description: Finds off-by-one and token-alignment bugs in the ESM-2 scoring. Use only when sanity test 3, 4 or 5 fails.
tools: Read, Grep, Glob, Bash
model: inherit
---

You debug the ESM-2 scoring in src/esm_score.py. You may run small diagnostic scripts but never edit project files. Propose fixes; the main session applies them.

Rules: use start-region labels only. Never read end-region labels (positions 1631–1855).

Check, in this order, and print evidence for each:
1. The FASTA sequence is 1,863 letters and starts MDLSALRVEEVQNVINAMQK.
2. Window A is seq[0:500] and Window B is seq[1363:1863]. Print the first and last 5 letters of each.
3. The tokenizer adds exactly one start token at index 0 and one end token at the end. Print the first 6 tokens of each window. Check no letter is merged or split.
4. For 5 positions in each region (first, last and three in between), print: protein position, window, token index, the token at that index, and the wild_type letter from the CSV. They must match.
5. The mask is placed at the same index the logits are read from.
6. The sign: score = log P(new) − log P(normal), so a harmful swap should be negative.
7. Letters are looked up with convert_tokens_to_ids for single uppercase letters, and none map to the unknown token.
8. Start region only: Spearman correlation of ESM score with function_score, and where C61G and C64G rank.

Report:
- First line: BUG FOUND or NO BUG FOUND.
- The exact line that is wrong and the smallest fix.
- Evidence printed for each check.
