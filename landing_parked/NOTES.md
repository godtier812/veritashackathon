# Parked landing-page parts

**The three parked files were removed from the tree on 2026-11-15 to tidy the repository.** They were never deployed and nothing used them. This note is kept as a map of what they held. To get one back, restore it from git history:

```
git show 4fe0ff5:landing_parked/app_landing.py > app_landing.py
git show 4fe0ff5:landing_parked/style_landing.css > style_landing.css
```

`heatmap.png` was a byte copy of `outputs/heatmap.png`; `landing/build.py` now reads `outputs/heatmap.png` directly.

The original description follows.

---

These parts came out of `app/` when the demo became a plain tool. Nothing here is deployed, and nothing in `app/` points here. They are kept so a later landing page can reuse them.

## Files

- `app_landing.py`: the full previous `app/app.py`, unchanged. The parked parts are its functions and strings, listed below.
- `style_landing.css`: the full previous `app/style.css`, unchanged.
- `heatmap.png`: byte copy of `outputs/heatmap.png`, which used to live in `app/results/`.

## Parts in `app_landing.py`

| Part | Where | What it does |
|---|---|---|
| Hero | `hero_html`, `title_html` | Eyebrow, big serif title split into per-letter spans, lede, "BRCA1 is our exam" line, the legend of six marks |
| Stat chips | `stat_html`, called from `hero_html` | 326 positions (app_config `n_tested_positions`), lab-tested swap count (counted from app_data.csv), 225 end-region positions (bootstrap.json `n_positions`) |
| Sequence ticker | `ticker_html` | Two rows of real BRCA1 amino-acid letters (start and end regions) for a marquee behind the hero. Only shown with `?fx=max` |
| How it works | `how_html`, `step_svgs` | Four steps with small inline SVG diagrams. Reads C61G's score and start-region class counts from app_data.csv |
| Limits and credits | `credits_html` | Three limits plus grouped credits with links |
| Heat map tab | `heat_head_html`, `heat_caption_html`, `heatmap_path`, the `gr.Image` block in the Blocks layout | Shows heatmap.png with a caption |
| Scoreboard bars | `meter`, old `scoreboard_html` | Arm rows with animated bars and count-up catch rates |
| Surprise me | `surprise` | Picks a random lab-tested swap |
| Copy button | `copy_summary`, copy part of `numbers_panel`, click handler in `head_script` | Copies a plain-text summary of a result |
| Animation JS | `head_script` | Count-up, IntersectionObserver reveals (`.rv`), MutationObserver re-scan, scroll-progress bar, mobile scroll-into-view, `?fx=max` parallax and cursor light, `?motion=off` |

## Parts in `style_landing.css`

- Tokens on `:root` at the top of the file. The landing parts need `--bg`, `--surface`, `--sunk`, `--ink`, `--ink-2`, `--ink-3`, `--rule`, `--rule-strong`, `--accent`, `--accent-soft`, the six state colours and their `-soft` variants, `--s1` to `--s7`, `--radius`, `--radius-sm`, `--shadow`, `--ease`, `--t1` to `--t3`, `--serif`, `--sans`, `--mono`. The animation section adds `--spring` and `--ease-out`.
- Hero, stats and legend: `.hero`, `.stat`, `.stats`, `.legend`.
- Section heads, steps, limits, credits: `.section-head`, `.steps`, `.limits`, `.credits`, `.cols`.
- Heat map: `.heat-img`, `.heat`.
- Scoreboard bars: `.arm`, `.arm-bars`, `.metric`, `.meter`, `.arm-rate`, `.board-note`.
- Animations: everything after the second `:root` block (`.title-in`, `.rv`, `.scroll-progress`, the result-card sequence, gauge wipe and marker slide, `.steps .pc`, keyframes `letter`, `panel`, `card-in`, `wipe`, `drop`, `draw`, `curtain`, `settle-img`, `rise`, `settle`, `slide`), plus the reduced-motion and `.no-motion` overrides.
- `?fx=max` extras at the end: `.ticker`, `.tick-row`, the cursor light, tilt, `beacon`, `ripple`, `sheen`, `marquee`. Gradio repeats every rule under a more specific container selector, so `display` and `color` on the ticker need `!important`.

## Fonts

Instrument Serif (hero and headings), IBM Plex Sans (body), JetBrains Mono (values), all loaded from Google Fonts in `font_link`.
