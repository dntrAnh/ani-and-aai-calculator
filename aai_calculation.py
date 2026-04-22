#!/usr/bin/env python3
"""Phage-compatible Average Amino Acid Identity (AAI) using bidirectional best hits.

Approach:
  1. k-mer (k=4) pre-filter to find top candidate matches per protein (fast).
  2. BioPython local alignment on candidates to get true percent identity.
  3. Bidirectional best-hit (BBH) pairs → average identity = AAI.

This works for any organism including phages, with no external tools required.

Usage (standalone):
    python aai_calculation.py proteins/ --output outputs/aai_heatmap.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner

PROTEIN_EXTENSIONS = {".faa", ".fa", ".fas", ".fasta", ".fna", ".pep"}

# k-mer size for pre-filtering candidates
_KMER_K = 4
# Number of top k-mer candidates to align per query protein
_TOP_N = 5


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_aligner() -> PairwiseAligner:
    """Return a BioPython local aligner tuned for protein BBH."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def _protein_kmers(seq: str, k: int = _KMER_K) -> frozenset[str]:
    if len(seq) < k:
        return frozenset()
    return frozenset(seq[i : i + k] for i in range(len(seq) - k + 1))


def _kmer_jaccard(a: frozenset[str], b: frozenset[str]) -> float:
    if not a or not b:
        return 0.0
    union = len(a | b)
    return len(a & b) / union if union else 0.0


def _top_candidates(query_kmers: frozenset[str],
                    target_kmers: list[frozenset[str]],
                    top_n: int = _TOP_N) -> list[int]:
    """Return indices of the top-n targets by k-mer Jaccard similarity."""
    scored = sorted(
        range(len(target_kmers)),
        key=lambda j: _kmer_jaccard(query_kmers, target_kmers[j]),
        reverse=True,
    )
    return scored[:top_n]


def _pct_identity(seq_a: str, seq_b: str, aligner: PairwiseAligner) -> float:
    """Local-alignment percent identity: identities / alignment_length × 100.

    Returns 0.0 if the alignment covers less than 50% of the shorter sequence —
    this filters out spurious short local matches (e.g. a 2-residue hit that
    would otherwise report 100% identity).
    """
    if not seq_a or not seq_b:
        return 0.0
    la, lb = len(seq_a), len(seq_b)
    # Skip pairs whose lengths differ by more than 4× (very unlikely to be BBH)
    if min(la, lb) / max(la, lb) < 0.25:
        return 0.0
    try:
        aln = next(iter(aligner.align(seq_a, seq_b)))
    except StopIteration:
        return 0.0
    if aln.length == 0:
        return 0.0
    # Coverage filter: alignment must span ≥50% of the shorter sequence.
    # Without this, a 2-residue perfect match reports 100% identity.
    if aln.length / min(la, lb) < 0.50:
        return 0.0
    counts = aln.counts()
    return counts.identities / aln.length * 100


def _one_directional_best(seqs_query: list[str],
                          seqs_target: list[str],
                          target_kmers: list[frozenset[str]],
                          aligner: PairwiseAligner) -> dict[int, tuple[int, float]]:
    """Return {query_idx: (best_target_idx, pct_id)} using k-mer pre-filter."""
    best: dict[int, tuple[int, float]] = {}
    for i, seq_q in enumerate(seqs_query):
        q_kmers = _protein_kmers(seq_q)
        candidates = _top_candidates(q_kmers, target_kmers)
        best_id, best_j = 0.0, -1
        for j in candidates:
            pid = _pct_identity(seq_q, seqs_target[j], aligner)
            if pid > best_id:
                best_id, best_j = pid, j
        if best_j >= 0 and best_id > 0:
            best[i] = (best_j, best_id)
    return best


def _genome_aai(seqs_a: list[str], seqs_b: list[str],
                aligner: PairwiseAligner) -> float:
    """BBH-based AAI between two proteomes (returns 0–100)."""
    if not seqs_a or not seqs_b:
        return 0.0

    kmers_a = [_protein_kmers(s) for s in seqs_a]
    kmers_b = [_protein_kmers(s) for s in seqs_b]

    # A→B and B→A best hits
    a_to_b = _one_directional_best(seqs_a, seqs_b, kmers_b, aligner)
    b_to_a = _one_directional_best(seqs_b, seqs_a, kmers_a, aligner)

    # Collect bidirectional best hits
    bbh_ids: list[float] = []
    for i, (j, pid_ab) in a_to_b.items():
        if b_to_a.get(j, (-1,))[0] == i:
            pid_ba = b_to_a[j][1]
            bbh_ids.append((pid_ab + pid_ba) / 2.0)

    return float(np.mean(bbh_ids)) if bbh_ids else 0.0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_proteomes(files: list[Path]) -> dict[str, list[str]]:
    """Read protein FASTA files → {genome_stem: [seq1, seq2, ...]}."""
    proteomes: dict[str, list[str]] = {}
    for f in files:
        seqs = [str(r.seq).upper() for r in SeqIO.parse(f, "fasta") if r.seq]
        if not seqs:
            raise ValueError(f"No sequences found in {f}")
        proteomes[f.stem] = seqs
    return proteomes


def compute_aai_matrix(files: list[Path]) -> pd.DataFrame:
    """Compute all-vs-all AAI for a list of protein FASTA files.

    Returns a symmetric DataFrame (values 0–100) with genome stems as
    both index and columns.  Diagonal entries are 100.0.
    """
    if len(files) < 2:
        raise ValueError("Need at least 2 protein FASTA files for AAI.")

    proteomes = load_proteomes(files)
    names = list(proteomes.keys())
    n = len(names)
    matrix = np.full((n, n), np.nan)
    np.fill_diagonal(matrix, 100.0)

    aligner = _build_aligner()

    for i in range(n):
        for j in range(i + 1, n):
            aai = _genome_aai(proteomes[names[i]], proteomes[names[j]], aligner)
            matrix[i, j] = aai
            matrix[j, i] = aai  # symmetric

    return pd.DataFrame(matrix, index=names, columns=names)


# ---------------------------------------------------------------------------
# CLI entry-point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute all-vs-all AAI for protein FASTA files (phage-compatible)."
    )
    p.add_argument("inputs", nargs="+", type=Path,
                   help="Protein FASTA files or directories.")
    p.add_argument("-o", "--output", type=Path,
                   default=Path("outputs/aai_heatmap.png"),
                   help="Output heatmap PNG (default: outputs/aai_heatmap.png).")
    return p.parse_args()


def _collect(inputs: list[Path]) -> list[Path]:
    files: list[Path] = []
    for p in inputs:
        if not p.exists():
            raise FileNotFoundError(p)
        if p.is_dir():
            files.extend(
                c.resolve() for c in sorted(p.iterdir())
                if c.is_file() and c.suffix.lower() in PROTEIN_EXTENSIONS
            )
        elif p.suffix.lower() in PROTEIN_EXTENSIONS:
            files.append(p.resolve())
        else:
            raise ValueError(f"Not a supported file: {p}")
    return list(dict.fromkeys(files))


if __name__ == "__main__":
    import generate_heatmap

    args = _parse_args()
    files = _collect(args.inputs)
    if len(files) < 2:
        print("Error: need at least 2 protein FASTA files.", file=sys.stderr)
        sys.exit(2)

    df = compute_aai_matrix(files)
    saved = generate_heatmap.generate_from_df(
        df, args.output, title="AAI Identity (%)", vmin=0, vmax=100
    )
    print(f"Heatmap saved: {saved}")
