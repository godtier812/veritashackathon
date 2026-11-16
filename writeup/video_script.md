# BRCA Spellcheck — Two-Minute Video Script

**Total duration:** 2:00 (120 seconds)
**Format:** Screen recording with voice-over. No personal names or school names on screen or spoken.

---

## Part 1 — The problem (0:00–0:15 | 15 seconds)

[On screen: Title card reading "BRCA Spellcheck: an AI second opinion on uncertain BRCA1 mutations" with subtitle "This project is a prototype and is not medical advice." Transition to the Gradio demo landing page showing the mutation input box.]

**Voice-over:** When a patient's genetic test returns a "variant of uncertain significance," a genetic counselor faces a mutation nobody can call harmful or harmless. Most genes lack a lab map, so counselors need to know how far to trust computer predictions.

---

## Part 2 — The demo (0:15–1:15 | 60 seconds)

[On screen: In the Gradio app, type `C61G` into the Mutation box and click Check. Highlight the returned verdict (`Looks harmful`), model score (`-10.50`), start-region cutoffs (`-7.75` and `-3.03`), lab result (`Loss of function (harmful). Lab function score -1.738.`), and the held-out end-region track record line (`196 of 288`, `68.1%`). Then type `D2E` and click Check, highlighting `Looks harmless`, score `0.05`, lab result `Functional (harmless). Lab function score 0.089.`, and the track record line (`462 of 477`, `96.9%`). Finally type `A1823S` and click Check, highlighting `Looks harmless` (`+0.75`) next to `Loss of function (harmful). Lab function score -3.643.` and the `outputs/results/failures.md` table showing positions `1823–1824`.]

**Voice-over:** In BRCA Spellcheck, typing C61G returns "Looks harmful," score -10.50, beside the CRISPR lab result: loss of function, harmful, score -1.738. On the held-out end region, "Looks harmful" matched the lab on 196 of 288 clear-class mutations, 68.1%, with 92 harmless and 29 intermediate. Typing D2E returns "Looks harmless," score 0.05, matching the lab's functional, harmless score of 0.089. On the held-out end region, "Looks harmless" matched the lab on 462 of 477 clear-class mutations, 96.9%, with 15 harmful and 20 intermediate. Typing A1823S shows our explained wrong case—three of the five worst misses are at positions 1823 and 1824. All three arms missed A1823S, scoring +0.75 ("Looks harmless") against the lab's harmful -3.643; our guess is DNA splicing disruption invisible to a protein model. BRCA1 is the exam, not the use case: Findlay et al. measured 96.5% of single-letter DNA changes here, so counselors use the lab map on BRCA1 and this exam for unmapped genes.

---

## Part 3 — How it works & held-out results (1:15–1:45 | 30 seconds)

[On screen: Show `outputs/heatmap.png` (start region 1–101 and held-out end region 1631–1855), then display the headline table from `outputs/results/headline.md` comparing `arm0`, `arm1`, `arm2`, and `arm3` on the 1,166 held-out end-region mutations.]

**Voice-over:** ESM-2 read natural protein sequences, never patients or lab labels. We froze cutoffs on positions 1 to 101 and tested once on held-out positions 1631 to 1855: 267 harmful and 899 harmless mutations. Always-harmless arm 0 caught 0, missed 267, with 0 false alarms. BLOSUM62 arm 1 caught 232 (86.9%), missed 35, with 571 of 899 false alarms (63.5%). Zero-shot ESM-2 arm 2 caught 252 (94.4%), missed 15, with 437 of 899 false alarms (48.6%).

---

## Part 4 — Limits, next steps & disclaimer (1:45–2:00 | 15 seconds)

[On screen: Closing slide listing the splicing blind spot, prior work and AlphaMissense, the two next steps, and the full-screen disclaimer: "This project is a prototype and is not medical advice."]

**Voice-over:** Protein models are blind to DNA splicing—we could not test this—and stronger tools like AlphaMissense exist. Next, we would add a splicing predictor and repeat the exam on another gene. This project is a prototype and is not medical advice.
