"""Label-free ESM-2 zero-shot scorer, window-index verifier, and score cache CLI for BRCA1 missense mutations."""

import argparse
import os
import sys
from typing import Any

import numpy as np
import pandas as pd

from src import config, load

_BASE_COLUMNS: list[str] = [
    "position",
    "region",
    "wild_type",
    "new_letter",
    "mutant",
    "window",
    "token_index",
    "token_at_index",
    "esm_score",
    "entropy",
]

_FULL_COLUMNS: list[str] = _BASE_COLUMNS + [
    "is_lab_tested",
    "model_id",
    "transformers_version",
    "precision",
]

_LOCK_REFUSAL_MSG: str = (
    "Refusing to recompute: the current scores are locked "
    "(outputs/frozen.json exists or the test has started). Ask a human."
)

_ONLY_END_REFUSAL_MSG: str = (
    "Refusing to score a single end-region mutation before the test run. "
    "Use --letter-test-only to check end-region indexing."
)


def window_of(position: int) -> str:
    """Return 'A' for start-region positions or 'B' for end-region positions."""
    if isinstance(position, bool) or not isinstance(position, (int, np.integer)):
        raise ValueError(f"Position {position} is not in a scored region.")
    pos = int(position)
    start_lo, start_hi = config.START_REGION
    end_lo, end_hi = config.END_REGION
    if start_lo <= pos <= start_hi:
        return "A"
    if end_lo <= pos <= end_hi:
        return "B"
    raise ValueError(f"Position {position} is not in a scored region.")


def token_index_of(position: int) -> int:
    """Return the 1-based-offset ESM token index (after <cls>) for a scored position."""
    window = window_of(position)
    return int(position) - config.WINDOW_START[window] + 1


def window_sequence(seq: str, window: str) -> str:
    """Slice the 500-residue window ('A' or 'B') from the full BRCA1 wild-type sequence."""
    if window not in config.WINDOWS:
        raise ValueError(
            f"Unknown window '{window}': expected one of {sorted(config.WINDOWS)}."
        )
    start, end = config.WINDOWS[window]
    subseq = seq[start - 1 : end]
    expected_len = end - start + 1
    if len(subseq) != 500 or len(subseq) != expected_len:
        raise ValueError(
            f"Window {window} slice must have length 500, got {len(subseq)}."
        )
    return subseq


def scored_positions() -> list[int]:
    """Return the 326 distinct lab-tested positions sorted ascending."""
    mutations = load.load_mutations()
    positions = sorted({int(pos) for pos in mutations["position"]})
    if len(positions) != config.N_POSITIONS:
        raise ValueError(
            f"Expected {config.N_POSITIONS} distinct positions, got {len(positions)}."
        )
    return positions


def letter_token_ids(tokenizer: Any) -> list[int]:
    """Look up and validate the 20 amino-acid token IDs in fixed AA_ORDER."""
    unk_id = getattr(tokenizer, "unk_token_id", None)
    ids: list[int] = []
    for aa in config.AA_ORDER:
        raw_id = tokenizer.convert_tokens_to_ids(aa)
        if raw_id is None or isinstance(raw_id, bool) or not isinstance(raw_id, (int, np.integer)):
            raise ValueError(
                f"Letter token lookup failed: invalid token ID {raw_id!r} for '{aa}'."
            )
        tid = int(raw_id)
        if tid < 0:
            raise ValueError(
                f"Letter token lookup failed: negative token ID {tid} for '{aa}'."
            )
        if unk_id is not None and tid == int(unk_id):
            raise ValueError(
                f"Letter token lookup failed: '{aa}' mapped to unk_token_id ({unk_id})."
            )
        round_trip = tokenizer.convert_ids_to_tokens(tid)
        if round_trip != aa:
            raise ValueError(
                f"Letter token lookup failed: ID {tid} for '{aa}' mapped back to {round_trip!r}."
            )
        ids.append(tid)

    if len(set(ids)) != len(config.AA_ORDER):
        raise ValueError(
            f"Letter token lookup failed: expected {len(config.AA_ORDER)} distinct IDs, "
            f"got {len(set(ids))}."
        )
    return ids


def _tokenize_and_verify_window(tokenizer: Any, seq: str, window: str) -> list[str]:
    """Tokenize one 500-letter window and verify its 502 tokens match the sequence."""
    w_seq = window_sequence(seq, window)
    raw_ids = tokenizer(w_seq, return_tensors="pt")["input_ids"][0]
    input_ids = raw_ids.tolist() if hasattr(raw_ids, "tolist") else list(raw_ids)
    tokens = list(tokenizer.convert_ids_to_tokens(input_ids))

    def _fail(detail: str) -> RuntimeError:
        return RuntimeError(
            f"LETTER TEST FAILED: window {window} tokenization check failed: {detail}. "
            "Do not score anything. Stop and report."
        )

    if len(tokens) != 502:
        raise _fail(f"expected 502 tokens, got {len(tokens)}")
    if tokens[0] != tokenizer.cls_token:
        raise _fail(
            f"expected cls_token {tokenizer.cls_token!r} at index 0, got {tokens[0]!r}"
        )
    if tokens[501] != tokenizer.eos_token:
        raise _fail(
            f"expected eos_token {tokenizer.eos_token!r} at index 501, got {tokens[501]!r}"
        )
    for idx in range(1, 501):
        if tokens[idx] != w_seq[idx - 1]:
            raise _fail(
                f"token at index {idx} was {tokens[idx]!r}, expected {w_seq[idx - 1]!r}"
            )
    return tokens


def run_letter_test(tokenizer: Any, seq: str, positions: list[int]) -> pd.DataFrame:
    """Verify window tokenization and token-to-wild-type alignment across positions."""
    window_tokens = {
        window: _tokenize_and_verify_window(tokenizer, seq, window)
        for window in ("A", "B")
    }

    rows: list[dict[str, Any]] = []
    for pos in positions:
        w = window_of(pos)
        idx = token_index_of(pos)
        tok = window_tokens[w][idx]
        wt = seq[int(pos) - 1]
        rows.append(
            {
                "position": int(pos),
                "window": w,
                "token_index": idx,
                "token_at_index": tok,
                "wild_type": wt,
                "ok": bool(tok == wt),
            }
        )

    columns = ["position", "window", "token_index", "token_at_index", "wild_type", "ok"]
    df = pd.DataFrame(rows, columns=columns)
    if df.empty:
        df["ok"] = df["ok"].astype(bool)
    else:
        bad = df[~df["ok"]]
        if not bad.empty:
            first = bad.iloc[0]
            raise RuntimeError(
                f"LETTER TEST FAILED at {len(bad)} of {len(df)} positions. "
                f"First failure: position {int(first['position'])}, "
                f"window {first['window']}, token index {int(first['token_index'])}, "
                f"token '{first['token_at_index']}', expected '{first['wild_type']}'. "
                "Do not score anything. Stop and report."
            )

    print(f"LETTER TEST PASSED: {len(positions)} of {len(positions)} positions, both windows 502 tokens.")
    return df


def choose_device(allow_cpu: bool) -> str:
    """Return 'cuda' if available, else 'cpu' when allow_cpu is True, or raise."""
    import torch

    if torch.cuda.is_available():
        return "cuda"
    if allow_cpu:
        print(
            "Warning: CUDA is not available; scoring on CPU will be slow.",
            file=sys.stderr,
        )
        return "cpu"
    raise RuntimeError(
        "CUDA is not available. Stage 3 needs an NVIDIA GPU. "
        "See BUILD_SPEC.md section 5.5. To run on CPU anyway (slow), "
        "pass --allow-cpu after asking a human."
    )


def load_model(model_id: str, device: str) -> tuple[Any, Any]:
    """Load the ESM-2 tokenizer and masked language model in float32 eval mode."""
    from transformers import AutoTokenizer, EsmForMaskedLM

    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = EsmForMaskedLM.from_pretrained(model_id)
    model.eval()
    model.to(device)
    return tokenizer, model


def _compute_entropy(logp20: np.ndarray) -> float:
    """Compute Shannon entropy (in nats) over the 20 amino-acid log-probabilities."""
    probs = np.exp(logp20)
    probs = probs / probs.sum()
    nonzero = probs[probs > 0.0]
    entropy = float(-(nonzero * np.log(nonzero)).sum()) + 0.0
    if not (0.0 <= entropy <= 2.9958):
        raise AssertionError(f"Entropy {entropy} outside expected [0.0, 2.9958] bounds.")
    return entropy


def score_positions(
    model: Any,
    tokenizer: Any,
    seq: str,
    positions: list[int],
    device: str,
    batch_size: int = config.ESM_BATCH_SIZE,
) -> pd.DataFrame:
    """Run masked-LM scoring over the requested positions and return 19 rows per position."""
    import torch

    if isinstance(batch_size, bool) or not isinstance(batch_size, (int, np.integer)) or batch_size < 1:
        raise ValueError(f"batch_size must be a positive integer, got {batch_size!r}.")
    batch_size_int = int(batch_size)

    if len({int(p) for p in positions}) != len(positions):
        raise ValueError("positions must not contain duplicates.")

    ids20 = letter_token_ids(tokenizer)
    run_letter_test(tokenizer, seq, positions)

    rows: list[dict[str, Any]] = []
    for window in ("A", "B"):
        win_positions = sorted(int(pos) for pos in positions if window_of(pos) == window)
        if not win_positions:
            continue

        w_seq = window_sequence(seq, window)
        base_ids = tokenizer(w_seq, return_tensors="pt")["input_ids"]
        if not isinstance(base_ids, torch.Tensor):
            base_ids = torch.as_tensor(base_ids, dtype=torch.long)
        if base_ids.ndim == 1:
            base_ids = base_ids.unsqueeze(0)

        total = len(win_positions)
        done = 0
        for start_idx in range(0, total, batch_size_int):
            batch_positions = win_positions[start_idx : start_idx + batch_size_int]
            k = len(batch_positions)
            batch = base_ids.repeat(k, 1).clone()

            batch_meta: list[tuple[int, int, str, int]] = []
            for j, pos in enumerate(batch_positions):
                t_idx = token_index_of(pos)
                wt = seq[pos - 1]
                wt_index = config.AA_ORDER.index(wt)
                wt_tok_id = ids20[wt_index]
                observed_tok_id = int(batch[j, t_idx].item())
                if observed_tok_id != wt_tok_id:
                    raise AssertionError(
                        f"Token ID mismatch at position {pos} (index {t_idx}): "
                        f"got {observed_tok_id}, expected {wt_tok_id}."
                    )
                batch[j, t_idx] = tokenizer.mask_token_id
                batch_meta.append((pos, t_idx, wt, wt_index))

            mask_ok = bool((batch == tokenizer.mask_token_id).sum(dim=1).eq(1).all().item())
            if not mask_ok:
                raise AssertionError("Each batch row must contain exactly one masked token.")

            batch_device = batch.to(device)
            with torch.no_grad():
                logits = model(input_ids=batch_device).logits

            for j, (pos, t_idx, wt, wt_index) in enumerate(batch_meta):
                row_logits = logits[j, t_idx]
                logp = torch.log_softmax(row_logits.float(), dim=-1)
                logp20 = logp[ids20].detach().cpu().numpy().astype(np.float64)

                entropy = _compute_entropy(logp20)
                rounded_entropy = round(entropy, 6) + 0.0
                region = load.region_of(pos)
                wt_logp = float(logp20[wt_index])

                for b_idx, new_aa in enumerate(config.AA_ORDER):
                    if new_aa == wt:
                        continue
                    esm_score = round(float(logp20[b_idx]) - wt_logp, 6) + 0.0
                    rows.append(
                        {
                            "position": pos,
                            "region": region,
                            "wild_type": wt,
                            "new_letter": new_aa,
                            "mutant": f"{wt}{pos}{new_aa}",
                            "window": window,
                            "token_index": t_idx,
                            "token_at_index": wt,
                            "esm_score": esm_score,
                            "entropy": rounded_entropy,
                        }
                    )

            done += k
            print(f"window {window}: {done}/{total} positions")

    if not rows:
        return pd.DataFrame(columns=_BASE_COLUMNS)

    df = pd.DataFrame(rows, columns=_BASE_COLUMNS)
    aa_rank = {aa: idx for idx, aa in enumerate(config.AA_ORDER)}
    df["_aa_rank"] = df["new_letter"].map(aa_rank)
    df = df.sort_values(["position", "_aa_rank"], kind="mergesort").drop(
        columns=["_aa_rank"]
    )
    return df.reset_index(drop=True)


def _write_scores_cache(
    base_df: pd.DataFrame,
    model_id: str,
) -> pd.DataFrame:
    """Attach metadata columns, enforce Section 6.4 schema, write atomically, and validate."""
    import transformers

    tested_mutants = set(load.load_mutations()["mutant"])
    out_df = base_df.copy()
    out_df["is_lab_tested"] = out_df["mutant"].isin(tested_mutants).astype("int64")
    out_df["model_id"] = str(model_id)
    out_df["transformers_version"] = str(transformers.__version__)
    out_df["precision"] = "float32"
    out_df = out_df[_FULL_COLUMNS]

    scores_path = config.ESM_SCORES_PATH
    scores_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = scores_path.parent / f"{scores_path.name}.tmp"
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        out_df.to_csv(
            handle,
            index=False,
            lineterminator="\n",
            float_format="%.6f",
        )
    os.replace(tmp_path, scores_path)
    return load.load_scores()


def _is_cache_locked() -> bool:
    """Return True once frozen.json exists or the test run has started."""
    return config.FROZEN_PATH.exists() or load.test_attempted()


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for letter testing, cache validation, single-mutant checks, and full scoring."""
    config.set_seeds()

    parser = argparse.ArgumentParser(
        description="Run ESM-2 zero-shot scoring and manage outputs/esm_scores.csv."
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Score one mutation (e.g. C61G) without writing the cache.",
    )
    parser.add_argument(
        "--letter-test-only",
        action="store_true",
        help="Run the tokenizer letter test without loading model weights.",
    )
    parser.add_argument(
        "--validate-cache",
        action="store_true",
        help="Validate outputs/esm_scores.csv and print its row count, model ID, and SHA-256.",
    )
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="Recompute all scores and overwrite outputs/esm_scores.csv if unlocked.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=config.ESM_BATCH_SIZE,
        help="Batch size for ESM-2 forward passes.",
    )
    parser.add_argument(
        "--allow-cpu",
        action="store_true",
        help="Allow CPU execution when CUDA is unavailable (slow).",
    )
    parser.add_argument(
        "--model-id",
        type=str,
        default=config.MODEL_ID,
        help="Hugging Face ESM-2 model ID.",
    )
    args = parser.parse_args(argv)

    try:
        if args.letter_test_only:
            from transformers import AutoTokenizer

            tokenizer = AutoTokenizer.from_pretrained(args.model_id)
            run_letter_test(tokenizer, load.load_wild_type(), scored_positions())
            return 0

        if args.validate_cache:
            df = load.load_scores()
            digest = config.sha256_of_file(config.ESM_SCORES_PATH)
            print(
                f"Validated cache: {len(df)} rows, model {df['model_id'].iloc[0]}, "
                f"sha256 {digest}"
            )
            return 0

        if args.only is not None:
            wt, pos, _ = load.parse_mutant(args.only)
            region = load.region_of(pos)
            if region == "end" and not load.test_run_done():
                print(_ONLY_END_REFUSAL_MSG, file=sys.stderr)
                return 1

            seq = load.load_wild_type()
            if seq[pos - 1] != wt:
                raise ValueError(
                    f"Wild-type mismatch for {args.only}: FASTA has '{seq[pos - 1]}' "
                    f"at position {pos}, expected '{wt}'."
                )

            device = choose_device(args.allow_cpu)
            tokenizer, model = load_model(args.model_id, device)
            run_letter_test(tokenizer, seq, scored_positions())
            df = score_positions(model, tokenizer, seq, [pos], device, batch_size=1)
            matches = df.loc[df["mutant"] == args.only]
            if len(matches) != 1:
                raise ValueError(f"Expected 1 scored row for {args.only}, got {len(matches)}.")
            row = matches.iloc[0]
            print(
                f"{args.only}: position {pos}, window {row['window']}, "
                f"token index {row['token_index']}, token {row['token_at_index']}, "
                f"esm_score = {row['esm_score']:.6f}, entropy = {row['entropy']:.6f}"
            )
            return 0

        if args.recompute and _is_cache_locked():
            print(_LOCK_REFUSAL_MSG, file=sys.stderr)
            return 1

        if not args.recompute:
            if config.ESM_SCORES_PATH.exists():
                df = load.load_scores()
                print(
                    f"Loaded cached scores: {len(df)} rows, model {df['model_id'].iloc[0]}."
                )
                return 0
            if _is_cache_locked():
                print(_LOCK_REFUSAL_MSG, file=sys.stderr)
                return 1

        device = choose_device(args.allow_cpu)
        seq = load.load_wild_type()
        positions = scored_positions()
        tokenizer, model = load_model(args.model_id, device)
        base_df = score_positions(
            model, tokenizer, seq, positions, device, batch_size=args.batch_size
        )
        validated = _write_scores_cache(base_df, model_id=args.model_id)
        print(f"Wrote outputs/esm_scores.csv: {len(validated)} rows")
        print(config.sha256_of_file(config.ESM_SCORES_PATH))
        return 0
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
