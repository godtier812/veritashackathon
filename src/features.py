"""Arm 3 biochemical and sequence feature engineering, swap-kind classification, and out-of-fold logistic regression fitting."""

from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from src import baselines, config, cutoffs


def charge(letter: str) -> int:
    """Return formal charge at cellular pH (+1 for K/R, -1 for D/E, 0 otherwise including neutral H)."""
    if not isinstance(letter, str) or len(letter) != 1 or letter not in config.AA_ORDER:
        raise ValueError(f"Unknown amino-acid letter: {letter!r}.")
    if letter in config.POSITIVE_LETTERS:
        return 1
    if letter in config.NEGATIVE_LETTERS:
        return -1
    # Histidine (H) is counted as neutral (0) at intracellular pH.
    return 0


def swap_kind(wild_type: str, new_letter: str) -> str:
    """Classify a single-amino-acid substitution into one of the four priority-ordered swap kinds."""
    wt_charge = charge(wild_type)
    new_charge = charge(new_letter)
    if new_letter == "P":
        return "to_proline"
    if wild_type == "C":
        return "loss_of_cysteine"
    if wt_charge != new_charge:
        return "charge_change"
    return "other"


def build_feature_table(scores: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of scores with Arm 1/2 harm scores, biochemical features, and swap_kind added.

    Sources for biochemical scales:
    - Hydrophobicity: Kyte J, Doolittle RF (1982), J. Mol. Biol. 157:105-132.
    - Residue volume (size, Å³): Zamyatnin AA (1972), Prog. Biophys. Mol. Biol. 24:107-123.
    """
    missing = [
        col
        for col in ("wild_type", "new_letter", "esm_score")
        if col not in scores.columns
    ]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}.")

    result = scores.copy()
    wt_list = result["wild_type"].tolist()
    new_list = result["new_letter"].tolist()
    esm_list = result["esm_score"].tolist()

    blosum_vals = [
        baselines.blosum_score(wt, new) for wt, new in zip(wt_list, new_list)
    ]
    harm_arm1_vals = [0.0 - float(val) for val in blosum_vals]
    harm_arm2_vals = [0.0 - float(val) for val in esm_list]
    # Signed differences: new_letter minus wild_type (Zamyatnin 1972 volume; Kyte-Doolittle 1982 hydropathy).
    size_changes = [
        float(config.VOLUME[new] - config.VOLUME[wt])
        for wt, new in zip(wt_list, new_list)
    ]
    hydro_changes = [
        float(config.HYDROPHOBICITY[new] - config.HYDROPHOBICITY[wt])
        for wt, new in zip(wt_list, new_list)
    ]
    to_proline_vals = [1 if new == "P" else 0 for new in new_list]
    from_cys_vals = [1 if wt == "C" else 0 for wt in wt_list]
    from_gly_vals = [1 if wt == "G" else 0 for wt in wt_list]
    swap_kinds = [swap_kind(wt, new) for wt, new in zip(wt_list, new_list)]

    result["blosum_score"] = pd.Series(blosum_vals, index=result.index, dtype="int64")
    result["harm_arm1"] = pd.Series(harm_arm1_vals, index=result.index, dtype="float64")
    result["harm_arm2"] = pd.Series(harm_arm2_vals, index=result.index, dtype="float64")
    result["size_change"] = pd.Series(size_changes, index=result.index, dtype="float64")
    result["hydrophobicity_change"] = pd.Series(
        hydro_changes, index=result.index, dtype="float64"
    )
    result["to_proline"] = pd.Series(to_proline_vals, index=result.index, dtype="int64")
    result["from_cysteine"] = pd.Series(
        from_cys_vals, index=result.index, dtype="int64"
    )
    result["from_glycine"] = pd.Series(from_gly_vals, index=result.index, dtype="int64")
    result["swap_kind"] = pd.Series(swap_kinds, index=result.index, dtype=object)
    return result


def feature_matrix(table: pd.DataFrame) -> np.ndarray:
    """Extract the 8-column float64 feature matrix in config.FEATURE_ORDER and verify no NaNs."""
    missing = [col for col in config.FEATURE_ORDER if col not in table.columns]
    if missing:
        raise ValueError(f"Missing required feature column(s): {', '.join(missing)}.")
    mat = table[config.FEATURE_ORDER].to_numpy(dtype="float64")
    if mat.ndim != 2 or mat.shape[1] != len(config.FEATURE_ORDER):
        raise ValueError(
            f"Expected feature matrix with {len(config.FEATURE_ORDER)} columns, got shape {mat.shape}."
        )
    if np.isnan(mat).any():
        raise ValueError("Feature matrix must not contain NaN.")
    return mat


def arm3_probability(
    x: ArrayLike,
    scaler_mean: ArrayLike,
    scaler_scale: ArrayLike,
    coefficients: ArrayLike,
    intercept: float,
) -> np.ndarray:
    """Compute Arm 3 predicted probabilities directly from frozen standardization and logistic parameters."""
    n_features = len(config.FEATURE_ORDER)
    x_arr = np.asarray(x, dtype=np.float64)
    mean_arr = np.asarray(scaler_mean, dtype=np.float64)
    scale_arr = np.asarray(scaler_scale, dtype=np.float64)
    coef_arr = np.asarray(coefficients, dtype=np.float64)

    if x_arr.ndim != 2 or x_arr.shape[1] != n_features:
        raise ValueError(f"x must be a 2-D array with {n_features} columns.")
    if mean_arr.ndim != 1 or len(mean_arr) != n_features:
        raise ValueError(f"scaler_mean must have length {n_features}.")
    if scale_arr.ndim != 1 or len(scale_arr) != n_features:
        raise ValueError(f"scaler_scale must have length {n_features}.")
    if coef_arr.ndim != 1 or len(coef_arr) != n_features:
        raise ValueError(f"coefficients must have length {n_features}.")
    if isinstance(intercept, bool) or not isinstance(
        intercept, (int, float, np.number)
    ):
        raise ValueError("intercept must be a numeric float.")
    intercept_val = float(intercept)

    if (
        np.isnan(x_arr).any()
        or np.isnan(mean_arr).any()
        or np.isnan(scale_arr).any()
        or np.isnan(coef_arr).any()
        or np.isnan(intercept_val)
    ):
        raise ValueError("Inputs to arm3_probability must not contain NaN.")
    if (scale_arr == 0.0).any():
        raise ValueError("scaler_scale must not contain zero.")

    z = (x_arr - mean_arr) / scale_arr
    logit = z @ coef_arr + intercept_val
    return 1.0 / (1.0 + np.exp(-logit))


def fit_arm3(
    x: ArrayLike,
    y: ArrayLike,
    groups: ArrayLike,
) -> dict[str, Any]:
    """Fit Arm 3 out-of-fold probabilities for cutoff selection, then refit on all rows."""
    import sklearn.model_selection as model_selection
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    n_features = len(config.FEATURE_ORDER)
    x_arr = np.asarray(x, dtype=np.float64)
    y_raw = np.asarray(y)
    groups_arr = np.asarray(groups)

    if x_arr.ndim != 2 or x_arr.shape[1] != n_features:
        raise ValueError(f"x must be a 2-D array with {n_features} columns.")
    if np.isnan(x_arr).any():
        raise ValueError("x must not contain NaN.")
    if y_raw.ndim != 1 or groups_arr.ndim != 1:
        raise ValueError("y and groups must be 1-D arrays.")
    if len(x_arr) == 0 or len(x_arr) != len(y_raw) or len(x_arr) != len(groups_arr):
        raise ValueError("x, y, and groups must be non-empty and have the same length.")
    if pd.isna(y_raw).any() or pd.isna(groups_arr).any():
        raise ValueError("y and groups must not contain NaN.")

    if y_raw.dtype == np.bool_ or not (
        np.issubdtype(y_raw.dtype, np.integer)
        or np.issubdtype(y_raw.dtype, np.floating)
    ):
        raise ValueError("y must hold only 0 and 1.")
    if not ((y_raw == 0) | (y_raw == 1)).all():
        raise ValueError("y must hold only 0 and 1.")
    y_arr = y_raw.astype(np.int64)
    if set(np.unique(y_arr).tolist()) != {0, 1}:
        raise ValueError("y must contain both class 0 and class 1.")

    splitter = model_selection.StratifiedGroupKFold(
        n_splits=config.N_FOLDS,
        shuffle=True,
        random_state=config.SEED,
    )

    oof = np.full(len(y_arr), np.nan, dtype=np.float64)
    for train_index, test_index in splitter.split(x_arr, y_arr, groups_arr):
        train_groups = set(groups_arr[train_index].tolist())
        test_groups = set(groups_arr[test_index].tolist())
        if train_groups & test_groups:
            raise ValueError(
                "Group leakage detected: position appears in both train and test folds."
            )

        fold_pipeline = make_pipeline(
            StandardScaler(),
            LogisticRegression(random_state=config.SEED),
        )
        fold_pipeline.fit(x_arr[train_index], y_arr[train_index])
        if list(fold_pipeline.classes_) != [0, 1]:
            raise RuntimeError(
                f"Expected pipeline classes_ [0, 1], got {list(fold_pipeline.classes_)}."
            )
        oof[test_index] = fold_pipeline.predict_proba(x_arr[test_index])[:, 1]

    if np.isnan(oof).any():
        raise RuntimeError("Out-of-fold probabilities contain NaN.")

    cutoff = (
        float(cutoffs.choose_catch_cutoff(oof, y_arr, config.CATCH_TARGET_PERCENT))
        + 0.0
    )

    final_pipeline = make_pipeline(
        StandardScaler(),
        LogisticRegression(random_state=config.SEED),
    )
    final_pipeline.fit(x_arr, y_arr)

    scaler = final_pipeline.named_steps["standardscaler"]
    model = final_pipeline.named_steps["logisticregression"]
    if list(model.classes_) != [0, 1]:
        raise RuntimeError(
            f"Expected final model classes_ [0, 1], got {list(model.classes_)}."
        )

    scaler_mean = [float(v) for v in scaler.mean_.tolist()]
    scaler_scale = [float(v) for v in scaler.scale_.tolist()]
    coefficients = [float(v) for v in model.coef_[0].tolist()]
    intercept = float(model.intercept_[0])
    n_iter = int(model.n_iter_[0])

    manual_probs = arm3_probability(
        x_arr, scaler_mean, scaler_scale, coefficients, intercept
    )
    sklearn_probs = final_pipeline.predict_proba(x_arr)[:, 1]
    if not np.allclose(manual_probs, sklearn_probs, atol=1e-9, rtol=0.0):
        raise RuntimeError("Arm 3 manual probability does not match scikit-learn.")

    converged = bool(n_iter < 100)
    if not converged:
        print(
            "Warning: Arm 3 LogisticRegression reached max_iter=100 without converging."
        )

    return {
        "scaler_mean": scaler_mean,
        "scaler_scale": scaler_scale,
        "coefficients": coefficients,
        "intercept": intercept,
        "cutoff": cutoff,
        "converged": converged,
        "n_iter": n_iter,
        "out_of_fold_probabilities": oof,
    }
