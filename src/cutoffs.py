"""Pure cutoff selection, flagging, band verdict, and outcome accounting rules for BRCA Spellcheck."""

import numpy as np
from numpy.typing import ArrayLike

from src import config


def _to_harmful_mask(raw_harmful: np.ndarray) -> np.ndarray:
    """Convert a validated 1-D 0/1 or boolean array into a boolean mask."""
    if raw_harmful.dtype == np.bool_:
        return raw_harmful.copy()
    if np.issubdtype(raw_harmful.dtype, np.integer):
        if not ((raw_harmful == 0) | (raw_harmful == 1)).all():
            raise ValueError("is_harmful must contain only 0, 1, or boolean values.")
        return raw_harmful == 1
    if np.issubdtype(raw_harmful.dtype, np.floating):
        if np.isnan(raw_harmful).any():
            raise ValueError("is_harmful must not contain NaN.")
        if not ((raw_harmful == 0.0) | (raw_harmful == 1.0)).all():
            raise ValueError("is_harmful must contain only 0, 1, or boolean values.")
        return raw_harmful == 1.0
    if raw_harmful.dtype == object:
        mask: list[bool] = []
        for val in raw_harmful.tolist():
            if (
                val is None
                or (isinstance(val, (float, np.floating)) and np.isnan(val))
                or type(val).__name__ == "NAType"
            ):
                raise ValueError("is_harmful must not contain NaN.")
            if isinstance(val, (bool, np.bool_)):
                mask.append(bool(val))
            elif isinstance(val, (int, np.integer, float, np.floating)) and val in (0, 1):
                mask.append(val == 1)
            else:
                raise ValueError("is_harmful must contain only 0, 1, or boolean values.")
        return np.asarray(mask, dtype=bool)
    raise ValueError("is_harmful must contain only 0, 1, or boolean values.")


def _validate_inputs(
    harm_scores: ArrayLike,
    is_harmful: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate 1-D harm scores and binary labels and return (float64 array, bool mask)."""
    harm_arr = np.asarray(harm_scores, dtype=np.float64)
    raw_harmful = np.asarray(is_harmful)

    if harm_arr.ndim != 1:
        raise ValueError("harm_scores must be a 1-D array.")
    if raw_harmful.ndim != 1:
        raise ValueError("is_harmful must be a 1-D array.")
    if len(harm_arr) != len(raw_harmful):
        raise ValueError("harm_scores and is_harmful must have the same length.")
    if len(harm_arr) == 0:
        raise ValueError("harm_scores and is_harmful must be non-empty.")
    if np.isnan(harm_arr).any():
        raise ValueError("harm_scores must not contain NaN.")

    harmful_mask = _to_harmful_mask(raw_harmful)
    if int(harmful_mask.sum()) < 1:
        raise ValueError("is_harmful must contain at least one harmful row.")
    if int((~harmful_mask).sum()) < 1:
        raise ValueError("is_harmful must contain at least one harmless row.")
    return harm_arr, harmful_mask


def flag(harm_scores: ArrayLike, cutoff: float) -> np.ndarray:
    """Return a boolean array marking rows with harm_scores >= cutoff."""
    harm_arr = np.asarray(harm_scores, dtype=np.float64)
    if harm_arr.ndim != 1:
        raise ValueError("harm_scores must be a 1-D array.")
    if np.isnan(harm_arr).any():
        raise ValueError("harm_scores must not contain NaN.")
    if cutoff is None or isinstance(cutoff, bool):
        raise ValueError("cutoff must be a numeric float value.")
    cutoff_val = float(cutoff)
    if np.isnan(cutoff_val):
        raise ValueError("cutoff must not be NaN.")
    return harm_arr >= cutoff_val


def choose_catch_cutoff(
    harm_scores: ArrayLike,
    is_harmful: ArrayLike,
    target_percent: int = 90,
) -> float:
    """Return the highest observed harm score flagging at least target_percent% of harmful rows."""
    harm_arr, harmful_mask = _validate_inputs(harm_scores, is_harmful)
    if (
        isinstance(target_percent, bool)
        or not isinstance(target_percent, (int, np.integer))
        or not (1 <= int(target_percent) <= 100)
    ):
        raise ValueError("target_percent must be an integer between 1 and 100.")

    target = int(target_percent)
    harmful_sorted = np.sort(harm_arr[harmful_mask])[::-1]
    n_harmful = len(harmful_sorted)
    needed = -(-(target * n_harmful) // 100)
    val = float(harmful_sorted[needed - 1])
    return 0.0 if val == 0.0 else val


def choose_false_alarm_cutoff(
    harm_scores: ArrayLike,
    is_harmful: ArrayLike,
    target_percent: int = 10,
) -> float | None:
    """Return the lowest observed harm score flagging at most target_percent% of harmless rows."""
    harm_arr, harmful_mask = _validate_inputs(harm_scores, is_harmful)
    if (
        isinstance(target_percent, bool)
        or not isinstance(target_percent, (int, np.integer))
        or not (0 <= int(target_percent) <= 100)
    ):
        raise ValueError("target_percent must be an integer between 0 and 100.")

    target = int(target_percent)
    harmless_scores = harm_arr[~harmful_mask]
    n_harmless = len(harmless_scores)
    allowed = (target * n_harmless) // 100
    candidates = np.unique(harm_arr)

    for candidate in candidates:
        cand_val = float(candidate)
        cand_val = 0.0 if cand_val == 0.0 else cand_val
        if int(flag(harmless_scores, cand_val).sum()) <= allowed:
            return cand_val
    return None


def has_band(cutoff_a: float, cutoff_b: float | None) -> bool:
    """Return True only when cutoff_b is present and strictly above cutoff_a."""
    return cutoff_b is not None and cutoff_b > cutoff_a


def band_verdict(
    harm_score: float,
    cutoff_a: float,
    cutoff_b: float | None,
) -> str:
    """Map a single harm score and optional band cutoffs to the fixed verdict text."""
    score = float(harm_score)
    cut_a = float(cutoff_a)
    cut_b = None if cutoff_b is None else float(cutoff_b)
    if np.isnan(score) or np.isnan(cut_a) or (cut_b is not None and np.isnan(cut_b)):
        raise ValueError("harm_score and cutoffs must not be NaN.")

    if has_band(cut_a, cut_b):
        assert cut_b is not None
        if bool(flag([score], cut_b)[0]):
            return config.VERDICT_HARMFUL
        if not bool(flag([score], cut_a)[0]):
            return config.VERDICT_HARMLESS
        return config.VERDICT_CANT_TELL

    if bool(flag([score], cut_a)[0]):
        return config.VERDICT_HARMFUL
    return config.VERDICT_HARMLESS


def count_outcomes(
    harm_scores: ArrayLike,
    is_harmful: ArrayLike,
    cutoff: float | None,
) -> dict[str, int]:
    """Count caught, missed, false alarm, and cleared rows at the given cutoff."""
    harm_arr, harmful_mask = _validate_inputs(harm_scores, is_harmful)
    if cutoff is None:
        flagged = np.zeros(len(harm_arr), dtype=bool)
    else:
        flagged = flag(harm_arr, cutoff)

    n_harmful = int(harmful_mask.sum())
    n_harmless = int((~harmful_mask).sum())
    harmful_caught = int((flagged & harmful_mask).sum())
    harmful_missed = n_harmful - harmful_caught
    false_alarms = int((flagged & (~harmful_mask)).sum())
    harmless_cleared = n_harmless - false_alarms

    return {
        "n_harmful": n_harmful,
        "n_harmless": n_harmless,
        "harmful_caught": harmful_caught,
        "harmful_missed": harmful_missed,
        "false_alarms": false_alarms,
        "harmless_cleared": harmless_cleared,
    }


def rates(counts: dict[str, int]) -> dict[str, float]:
    """Compute catch rate, false-alarm rate, and accuracy percentages rounded to 1 decimal."""
    n_harmful = counts["n_harmful"]
    n_harmless = counts["n_harmless"]
    if n_harmful <= 0 or n_harmless <= 0:
        raise ValueError("n_harmful and n_harmless must be positive.")

    total = n_harmful + n_harmless
    catch_rate = round(100.0 * counts["harmful_caught"] / n_harmful, 1)
    false_alarm_rate = round(100.0 * counts["false_alarms"] / n_harmless, 1)
    accuracy = round(
        100.0 * (counts["harmful_caught"] + counts["harmless_cleared"]) / total,
        1,
    )
    return {
        "catch_rate_percent": catch_rate,
        "false_alarm_rate_percent": false_alarm_rate,
        "accuracy_percent": accuracy,
    }
