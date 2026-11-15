"""Stage 7 static Matplotlib heat map of ESM-2 scores and lab-tested mutation classes."""

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src import config, load

FIGURE_TITLE: str = (
    "BRCA Spellcheck: ESM-2 scores for every swap at the lab-tested positions"
)
COLORBAR_LABEL: str = "ESM score (negative = model finds the swap surprising)"
HIDDEN_END_NOTE: str = "End region hidden until the one test run is done."
ZINC_LEGEND_LABEL: str = "zinc-holding letters (PDB 1JM7)"
LOF_LEGEND_LABEL: str = "LOF (harmful)"
NON_LOF_LEGEND_LABEL: str = "Lab-tested (FUNC / INT)"


def build_region_panel(
    scores_region: pd.DataFrame,
    labels_region: pd.DataFrame,
) -> dict[str, object]:
    """Pivot one region's ESM scores into a 20 x n_positions matrix and compute marker coordinates."""
    if scores_region.empty:
        raise ValueError("scores_region must not be empty.")

    positions = sorted(int(pos) for pos in scores_region["position"].unique())
    n_positions = len(positions)
    pos_to_col = {pos: idx for idx, pos in enumerate(positions)}
    aa_to_row = {aa: idx for idx, aa in enumerate(config.AA_ORDER)}

    grid = np.full((len(config.AA_ORDER), n_positions), np.nan, dtype=np.float64)
    wt_by_pos: dict[int, str] = {}

    for pos, wt, new_aa, esm_val in zip(
        scores_region["position"],
        scores_region["wild_type"],
        scores_region["new_letter"],
        scores_region["esm_score"],
    ):
        pos_int = int(pos)
        wt_str = str(wt)
        new_str = str(new_aa)
        if pos_int in wt_by_pos and wt_by_pos[pos_int] != wt_str:
            raise ValueError(f"Inconsistent wild_type at position {pos_int}.")
        wt_by_pos[pos_int] = wt_str
        row_idx = aa_to_row[new_str]
        col_idx = pos_to_col[pos_int]
        grid[row_idx, col_idx] = float(esm_val)

    for pos in positions:
        col_idx = pos_to_col[pos]
        col_vals = grid[:, col_idx]
        nan_mask = np.isnan(col_vals)
        if int(nan_mask.sum()) != 1:
            raise ValueError(
                f"Position {pos} column must have exactly 1 NaN cell, found {int(nan_mask.sum())}."
            )
        wt_letter = wt_by_pos[pos]
        wt_row = aa_to_row[wt_letter]
        if not np.isnan(grid[wt_row, col_idx]):
            raise ValueError(
                f"Position {pos} wild-type cell ({wt_letter}) is not NaN."
            )

    merged = labels_region.merge(
        scores_region[["mutant", "position", "new_letter", "esm_score"]],
        on=["mutant", "position", "new_letter"],
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != len(labels_region):
        raise ValueError(
            f"Label merge mismatch: expected {len(labels_region)} rows, got {len(merged)}."
        )

    lof_xs: list[int] = []
    lof_ys: list[int] = []
    non_lof_xs: list[int] = []
    non_lof_ys: list[int] = []

    for pos, new_aa, fclass in zip(
        merged["position"],
        merged["new_letter"],
        merged["findlay_class"],
    ):
        col_idx = pos_to_col[int(pos)]
        row_idx = aa_to_row[str(new_aa)]
        if np.isnan(grid[row_idx, col_idx]):
            raise ValueError(
                f"Lab-tested swap at position {pos} ({new_aa}) landed on a NaN wild-type cell."
            )
        if str(fclass) == config.HARMFUL_CLASS:
            lof_xs.append(col_idx)
            lof_ys.append(row_idx)
        else:
            non_lof_xs.append(col_idx)
            non_lof_ys.append(row_idx)

    return {
        "positions": positions,
        "pos_to_col": pos_to_col,
        "wt_by_pos": wt_by_pos,
        "grid": grid,
        "lof_xs": lof_xs,
        "lof_ys": lof_ys,
        "non_lof_xs": non_lof_xs,
        "non_lof_ys": non_lof_ys,
    }


def _select_position_ticks(positions: list[int], step: int) -> tuple[list[int], list[str]]:
    """Select clean X-axis tick indices and position labels including endpoints."""
    tick_indices: list[int] = [0]
    for idx, pos in enumerate(positions):
        if pos % step == 0 and idx not in tick_indices:
            if idx > 2 and (len(positions) - 1 - idx) > 2:
                tick_indices.append(idx)
    if (len(positions) - 1) not in tick_indices:
        tick_indices.append(len(positions) - 1)
    tick_indices.sort()
    tick_labels = [str(positions[i]) for i in tick_indices]
    return tick_indices, tick_labels


def _draw_panel(
    ax: plt.Axes,
    panel: dict[str, object],
    cmap: matplotlib.colors.Colormap,
    panel_title: str,
    x_label: str,
    tick_step: int,
    draw_zinc: bool,
) -> matplotlib.image.AxesImage:
    """Render one region's heatmap matrix, lab-tested markers, and optional zinc markers."""
    grid = np.asarray(panel["grid"], dtype=np.float64)
    positions: list[int] = list(panel["positions"])  # type: ignore[arg-type]
    pos_to_col: dict[int, int] = dict(panel["pos_to_col"])  # type: ignore[arg-type]
    lof_xs: list[int] = list(panel["lof_xs"])  # type: ignore[arg-type]
    lof_ys: list[int] = list(panel["lof_ys"])  # type: ignore[arg-type]
    non_lof_xs: list[int] = list(panel["non_lof_xs"])  # type: ignore[arg-type]
    non_lof_ys: list[int] = list(panel["non_lof_ys"])  # type: ignore[arg-type]

    im = ax.imshow(
        grid,
        cmap=cmap,
        vmin=config.HEATMAP_ZMIN,
        vmax=config.HEATMAP_ZMAX,
        aspect="auto",
        interpolation="nearest",
    )

    ax.scatter(
        non_lof_xs,
        non_lof_ys,
        s=15,
        c="#666666",
        marker="o",
        edgecolors="white",
        linewidths=0.25,
        label=NON_LOF_LEGEND_LABEL,
        zorder=3,
    )
    ax.scatter(
        lof_xs,
        lof_ys,
        s=34,
        c="#000000",
        marker="o",
        edgecolors="white",
        linewidths=0.45,
        label=LOF_LEGEND_LABEL,
        zorder=4,
    )

    if draw_zinc:
        zinc_cols = [
            pos_to_col[pos]
            for pos in config.ZINC_POSITIONS_TO_CONFIRM
            if pos in pos_to_col
        ]
        zinc_ys = [-1.05] * len(zinc_cols)
        ax.scatter(
            zinc_cols,
            zinc_ys,
            s=52,
            c="#d95f02",
            marker="v",
            edgecolors="#000000",
            linewidths=0.6,
            clip_on=False,
            label=ZINC_LEGEND_LABEL,
            zorder=5,
        )

    ax.set_xlim(-0.5, len(positions) - 0.5)
    ax.set_ylim(len(config.AA_ORDER) - 0.5, -0.5)

    ax.set_yticks(range(len(config.AA_ORDER)))
    ax.set_yticklabels(list(config.AA_ORDER), fontsize=8.5)
    ax.set_ylabel("New amino acid", fontsize=10)

    tick_indices, tick_labels = _select_position_ticks(positions, step=tick_step)
    ax.set_xticks(tick_indices)
    ax.set_xticklabels(tick_labels, fontsize=8.5)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_title(panel_title, fontsize=11, fontweight="bold", pad=18 if draw_zinc else 8)

    return im


def render_heatmap(output_path: Path | None = None) -> dict[str, object]:
    """Generate the static PNG heat map obeying the before/after test-run visibility rule."""
    target_path = config.HEATMAP_PNG_PATH if output_path is None else Path(output_path)
    done = load.test_run_done()
    all_scores = load.load_scores()

    start_scores = all_scores[all_scores["region"] == "start"].reset_index(drop=True)
    start_labels = load.load_labels("start", "tuning")
    start_panel = build_region_panel(start_scores, start_labels)

    start_black = len(start_panel["lof_xs"])  # type: ignore[arg-type]
    start_grey = len(start_panel["non_lof_xs"])  # type: ignore[arg-type]
    expected_start_black = config.EXPECTED_COUNTS["start"][config.HARMFUL_CLASS]
    expected_start_grey = (
        config.EXPECTED_COUNTS["start"][config.HARMLESS_CLASS]
        + config.EXPECTED_COUNTS["start"][config.EXCLUDED_CLASS]
    )
    if start_black != expected_start_black or start_grey != expected_start_grey:
        raise ValueError(
            f"Start panel marker count mismatch: got {start_black} black and {start_grey} grey, "
            f"expected {expected_start_black} and {expected_start_grey}."
        )

    cmap = matplotlib.colormaps[config.HEATMAP_COLORSCALE].copy()
    cmap.set_bad("white")

    if not done:
        fig, ax_start = plt.subplots(
            nrows=1,
            ncols=1,
            figsize=(14, 5.8),
             constrained_layout=False,
        )
        fig.subplots_adjust(top=0.76, bottom=0.17, left=0.06, right=0.88)
        im = _draw_panel(
            ax=ax_start,
            panel=start_panel,
            cmap=cmap,
            panel_title="Start region (RING domain, positions 1–101)",
            x_label="BRCA1 protein position (start region: 1–101)",
            tick_step=10,
            draw_zinc=True,
        )
        handles, labels = ax_start.get_legend_handles_labels()
        order = [
            labels.index(LOF_LEGEND_LABEL),
            labels.index(NON_LOF_LEGEND_LABEL),
            labels.index(ZINC_LEGEND_LABEL),
        ]
        fig.legend(
            [handles[i] for i in order],
            [labels[i] for i in order],
            loc="upper center",
            bbox_to_anchor=(0.47, 0.875),
            ncol=3,
            frameon=True,
            fontsize=9.5,
        )
        cbar_ax = fig.add_axes([0.90, 0.17, 0.018, 0.59])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label(COLORBAR_LABEL, fontsize=9.5)

        fig.suptitle(FIGURE_TITLE, fontsize=13, fontweight="bold", y=0.965)
        fig.text(
            0.47,
            0.92,
            config.DISCLAIMER,
            ha="center",
            va="center",
            fontsize=10,
            style="italic",
        )
        fig.text(
            0.47,
            0.03,
            HIDDEN_END_NOTE,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )
        total_black = start_black
        total_grey = start_grey
        panels_drawn = ["start"]
    else:
        end_scores = all_scores[all_scores["region"] == "end"].reset_index(drop=True)
        end_labels = load.load_labels("end", "post_test")
        end_panel = build_region_panel(end_scores, end_labels)

        end_black = len(end_panel["lof_xs"])  # type: ignore[arg-type]
        end_grey = len(end_panel["non_lof_xs"])  # type: ignore[arg-type]
        expected_end_black = config.EXPECTED_COUNTS["end"][config.HARMFUL_CLASS]
        expected_end_grey = (
            config.EXPECTED_COUNTS["end"][config.HARMLESS_CLASS]
            + config.EXPECTED_COUNTS["end"][config.EXCLUDED_CLASS]
        )
        if end_black != expected_end_black or end_grey != expected_end_grey:
            raise ValueError(
                f"End panel marker count mismatch: got {end_black} black and {end_grey} grey, "
                f"expected {expected_end_black} and {expected_end_grey}."
            )

        total_black = start_black + end_black
        total_grey = start_grey + end_grey
        if total_black != 398 or total_grey != 1439:
            raise ValueError(
                f"Total marker count mismatch: got {total_black} black and {total_grey} grey, "
                "expected 398 and 1439."
            )

        fig, (ax_start, ax_end) = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(16, 10.2),
            constrained_layout=False,
        )
        fig.subplots_adjust(top=0.82, bottom=0.08, left=0.055, right=0.89, hspace=0.36)

        im = _draw_panel(
            ax=ax_start,
            panel=start_panel,
            cmap=cmap,
            panel_title="Start region (RING domain, positions 1–101)",
            x_label="BRCA1 protein position (start region: 1–101)",
            tick_step=10,
            draw_zinc=True,
        )
        _draw_panel(
            ax=ax_end,
            panel=end_panel,
            cmap=cmap,
            panel_title="End region (BRCT domain, positions 1631–1855)",
            x_label="BRCA1 protein position (end region: 1631–1855)",
            tick_step=20,
            draw_zinc=False,
        )

        handles, labels = ax_start.get_legend_handles_labels()
        order = [
            labels.index(LOF_LEGEND_LABEL),
            labels.index(NON_LOF_LEGEND_LABEL),
            labels.index(ZINC_LEGEND_LABEL),
        ]
        fig.legend(
            [handles[i] for i in order],
            [labels[i] for i in order],
            loc="upper center",
            bbox_to_anchor=(0.47, 0.905),
            ncol=3,
            frameon=True,
            fontsize=9.5,
        )

        cbar_ax = fig.add_axes([0.905, 0.12, 0.016, 0.66])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label(COLORBAR_LABEL, fontsize=10)

        fig.suptitle(FIGURE_TITLE, fontsize=13.5, fontweight="bold", y=0.97)
        fig.text(
            0.47,
            0.932,
            config.DISCLAIMER,
            ha="center",
            va="center",
            fontsize=10.5,
            style="italic",
        )
        panels_drawn = ["start", "end"]

    target_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target_path.with_name(f"{target_path.stem}.tmp{target_path.suffix}")
    try:
        fig.savefig(tmp_path, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        os.replace(tmp_path, target_path)
    finally:
        plt.close(fig)
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    pos61_col = int(start_panel["pos_to_col"][61])  # type: ignore[index]
    start_grid = np.asarray(start_panel["grid"], dtype=np.float64)
    c_row = config.AA_ORDER.index("C")
    pos61_c_is_nan = bool(np.isnan(start_grid[c_row, pos61_col]))
    pos61_lof_letters = sorted(
        config.AA_ORDER[y]
        for x, y in zip(start_panel["lof_xs"], start_panel["lof_ys"])  # type: ignore[arg-type]
        if x == pos61_col
    )

    return {
        "output_path": target_path,
        "panels_drawn": panels_drawn,
        "start_black": start_black,
        "start_grey": start_grey,
        "total_black": total_black,
        "total_grey": total_grey,
        "pos61_c_is_nan": pos61_c_is_nan,
        "pos61_lof_letters": pos61_lof_letters,
    }


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for generating outputs/heatmap.png."""
    config.set_seeds()
    parser = argparse.ArgumentParser(
        description="Render the BRCA Spellcheck static ESM-2 heat map PNG."
    )
    parser.parse_args(argv)
    try:
        summary = render_heatmap()
        out_path = Path(summary["output_path"])  # type: ignore[arg-type]
        size_bytes = out_path.stat().st_size
        panels_str = ", ".join(summary["panels_drawn"])  # type: ignore[arg-type]
        pos61_letters = ",".join(summary["pos61_lof_letters"])  # type: ignore[arg-type]
        print(
            f"Wrote {out_path} ({size_bytes} bytes) | panels: {panels_str} | "
            f"black LOF dots: {summary['total_black']} | "
            f"grey non-LOF dots: {summary['total_grey']} | "
            f"pos 61 C is NaN: {summary['pos61_c_is_nan']} | "
            f"pos 61 LOF swaps: {pos61_letters}"
        )
    except (
        FileNotFoundError,
        ValueError,
        TypeError,
        PermissionError,
        RuntimeError,
        OSError,
    ) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
