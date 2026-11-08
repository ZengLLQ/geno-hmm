"""Copy Number Variation (CNV) 6-state definitions and theoretical signal parameters."""

from __future__ import annotations

from enum import IntEnum
from dataclasses import dataclass
import numpy as np


class CNVState(IntEnum):
    """6 discrete latent genomic copy number states."""

    HOM_DEL = 0       # State 0: Homozygous Deletion (CN = 0, nullisomy)
    HET_DEL = 1       # State 1: Hemizygous Deletion (CN = 1, monosomy)
    NORMAL = 2        # State 2: Normal Diploid Copy Neutral (CN = 2, euploid)
    CN_LOH = 3        # State 3: Copy-Neutral Loss of Heterozygosity (CN = 2, uniparental disomy)
    SINGLE_DUP = 4    # State 4: Single Copy Duplication (CN = 3, trisomy)
    DOUBLE_DUP = 5    # State 5: Double Copy Duplication / Amplification (CN = 4, tetrasomy)


NUM_STATES = 6


@dataclass(frozen=True)
class StateParameters:
    """Theoretical signal parameters for a copy number state."""

    state: CNVState
    copy_number: int
    name: str
    lrr_mean: float
    lrr_std: float
    baf_peaks: tuple[float, ...]  # Expected BAF peaks for possible genotypes


STATE_CONFIGS: dict[CNVState, StateParameters] = {
    CNVState.HOM_DEL: StateParameters(
        state=CNVState.HOM_DEL,
        copy_number=0,
        name="Homozygous Deletion",
        lrr_mean=-3.20,
        lrr_std=0.45,
        baf_peaks=(),  # Degenerate uniform noise
    ),
    CNVState.HET_DEL: StateParameters(
        state=CNVState.HET_DEL,
        copy_number=1,
        name="Hemizygous Deletion",
        lrr_mean=-0.55,
        lrr_std=0.22,
        baf_peaks=(0.0, 1.0),  # Genotypes {A, B} -> LOH
    ),
    CNVState.NORMAL: StateParameters(
        state=CNVState.NORMAL,
        copy_number=2,
        name="Normal Diploid",
        lrr_mean=0.00,
        lrr_std=0.18,
        baf_peaks=(0.0, 0.50, 1.0),  # Genotypes {AA, AB, BB}
    ),
    CNVState.CN_LOH: StateParameters(
        state=CNVState.CN_LOH,
        copy_number=2,
        name="Copy-Neutral LOH",
        lrr_mean=0.00,
        lrr_std=0.18,
        baf_peaks=(0.0, 1.0),  # Genotypes {AA, BB} only (no heterozygous AB band)
    ),
    CNVState.SINGLE_DUP: StateParameters(
        state=CNVState.SINGLE_DUP,
        copy_number=3,
        name="Single Duplication",
        lrr_mean=0.40,
        lrr_std=0.20,
        baf_peaks=(0.0, 0.33, 0.67, 1.0),  # Genotypes {AAA, AAB, ABB, BBB}
    ),
    CNVState.DOUBLE_DUP: StateParameters(
        state=CNVState.DOUBLE_DUP,
        copy_number=4,
        name="Double Duplication",
        lrr_mean=0.75,
        lrr_std=0.25,
        baf_peaks=(0.0, 0.25, 0.50, 0.75, 1.0),  # Genotypes {AAAA, AAAB, AABB, ABBB, BBBB}
    ),
}
