"""Genomic dual-track ASCII visualizer for Log R Ratio and B-Allele Frequency."""

from __future__ import annotations

import numpy as np

from geno_hmm.models.states import CNVState
from geno_hmm.algorithms.viterbi import CNVSegment


def render_ascii_profile(
    positions_bp: np.ndarray,
    lrr: np.ndarray,
    baf: np.ndarray,
    states: np.ndarray,
    width: int = 70,
    height_lrr: int = 10,
) -> str:
    """Render a dual-track ASCII visualization of LRR intensity and BAF allelic split.

    Parameters
    ----------
    positions_bp : np.ndarray
        Probe genomic positions.
    lrr : np.ndarray
        LRR values.
    baf : np.ndarray
        BAF values.
    states : np.ndarray
        Inferred hidden states.
    width : int
        Width of character canvas.
    height_lrr : int
        Height of LRR plot track.

    Returns
    -------
    chart_str : str
        Formatted multi-line text plot.
    """
    N = len(positions_bp)
    if N == 0:
        return "[Empty Profile]"

    # Bin probes into `width` bins
    bins = np.linspace(positions_bp[0], positions_bp[-1], width + 1)
    bin_idx = np.digitize(positions_bp, bins) - 1
    bin_idx = np.clip(bin_idx, 0, width - 1)

    bin_lrr = np.zeros(width)
    bin_state = np.full(width, 2)
    bin_counts = np.zeros(width)

    for i in range(N):
        b = bin_idx[i]
        bin_lrr[b] += lrr[i]
        bin_counts[b] += 1
        bin_state[b] = states[i]

    for b in range(width):
        if bin_counts[b] > 0:
            bin_lrr[b] /= bin_counts[b]

    lines: list[str] = []
    lines.append("=== Log R Ratio (LRR Intensity Profile) ===")
    lrr_min, lrr_max = -2.0, 1.0

    canvas = [[" " for _ in range(width)] for _ in range(height_lrr)]

    # Draw zero baseline
    zero_row = int(np.clip((lrr_max - 0.0) / (lrr_max - lrr_min) * (height_lrr - 1), 0, height_lrr - 1))
    for c in range(width):
        canvas[zero_row][c] = "."

    # Plot binned LRR points
    for c in range(width):
        val = bin_lrr[c]
        row = int(np.clip((lrr_max - val) / (lrr_max - lrr_min) * (height_lrr - 1), 0, height_lrr - 1))
        st = bin_state[c]
        if st in (0, 1):
            symbol = "v"  # Deletion
        elif st in (4, 5):
            symbol = "^"  # Duplication
        elif st == 3:
            symbol = "x"  # CN-LOH
        else:
            symbol = "*"  # Normal
        canvas[row][c] = symbol

    for r in range(height_lrr):
        y_val = lrr_max - r * (lrr_max - lrr_min) / (height_lrr - 1)
        row_str = "".join(canvas[r])
        lines.append(f"{y_val:+4.1f} |{row_str}|")

    lines.append("     +" + "-" * width + "+")
    lines.append(f"      {positions_bp[0] / 1e6:6.2f} Mb" + " " * (width - 25) + f"{positions_bp[-1] / 1e6:6.2f} Mb")

    # State call track
    lines.append("")
    lines.append("=== Genomic State Track (0=HOM_DEL, 1=HET_DEL, 2=NORMAL, 3=LOH, 4=DUP, 5=AMP) ===")
    state_chars = {
        0: "0",
        1: "-",
        2: "=",
        3: "L",
        4: "+",
        5: "#",
    }
    track_str = "".join(state_chars.get(int(bin_state[c]), "?") for c in range(width))
    lines.append("CALL |" + track_str + "|")

    return "\n".join(lines)


def format_clinical_table(segments: list[CNVSegment]) -> str:
    """Format detected CNVs into a clinical summary table."""
    if not segments:
        return "No high-confidence Copy Number Alterations detected."

    header = (
        f"{'State':<18} | {'Coordinates':<28} | {'Span (kb)':<10} | {'Probes':<8} | "
        f"{'Mean LRR':<10} | {'Score':<8}"
    )
    separator = "-" * len(header)
    rows = [header, separator]

    for seg in segments:
        coord = f"{seg.chromosome}:{seg.start_bp:,}-{seg.end_bp:,}"
        span_kb = f"{seg.length_bp / 1e3:.1f}"
        rows.append(
            f"{seg.state_name:<18} | {coord:<28} | {span_kb:<10} | {seg.num_probes:<8} | "
            f"{seg.mean_lrr:+8.3f}   | {seg.confidence_score:6.1f}"
        )

    return "\n".join(rows)
