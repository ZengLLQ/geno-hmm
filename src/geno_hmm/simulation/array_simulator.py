"""Microarray signal simulator with synthetic chromosomal alterations."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from geno_hmm.models.states import CNVState, STATE_CONFIGS


@dataclass
class MicroarrayData:
    """Simulated or clinical SNP microarray chromosomal data."""

    chromosome: str
    positions_bp: np.ndarray
    lrr: np.ndarray
    baf: np.ndarray
    true_states: np.ndarray | None = None


class MicroarraySimulator:
    """Generates synthetic high-density SNP microarray profiles with ground truth alterations."""

    def __init__(self, random_seed: int | None = 42) -> None:
        self.rng = np.random.default_rng(random_seed)

    def generate_chromosome(
        self,
        chromosome: str = "chr22",
        start_bp: int = 15_000_000,
        end_bp: int = 25_000_000,
        avg_probe_spacing_bp: int = 2_000,
        noise_level: float = 1.0,
    ) -> MicroarrayData:
        """Create baseline diploid chromosome track."""
        num_probes = (end_bp - start_bp) // avg_probe_spacing_bp
        # Generate non-uniform probe spacing
        spacings = self.rng.integers(
            int(avg_probe_spacing_bp * 0.5),
            int(avg_probe_spacing_bp * 1.5),
            size=num_probes,
        )
        positions = start_bp + np.cumsum(spacings)
        positions = positions[positions <= end_bp]
        N = len(positions)

        # Baseline: NORMAL (diploid, State 2)
        true_states = np.full(N, int(CNVState.NORMAL), dtype=np.int32)
        lrr = np.zeros(N, dtype=np.float64)
        baf = np.zeros(N, dtype=np.float64)

        norm_cfg = STATE_CONFIGS[CNVState.NORMAL]
        lrr += self.rng.normal(norm_cfg.lrr_mean, norm_cfg.lrr_std * noise_level, size=N)

        # Simulate diploid BAF (peaks at 0.0, 0.5, 1.0)
        # Genotype probabilities: AA (0.25), AB (0.50), BB (0.25)
        genotypes = self.rng.choice([0.0, 0.5, 1.0], p=[0.25, 0.50, 0.25], size=N)
        baf_noise = self.rng.normal(0.0, 0.035 * noise_level, size=N)
        baf = np.clip(genotypes + baf_noise, 0.0, 1.0)

        return MicroarrayData(
            chromosome=chromosome,
            positions_bp=positions,
            lrr=lrr,
            baf=baf,
            true_states=true_states,
        )

    def inject_alteration(
        self,
        data: MicroarrayData,
        start_bp: int,
        end_bp: int,
        state: CNVState,
        noise_level: float = 1.0,
    ) -> MicroarrayData:
        """Inject a specific CNV event into the chromosomal dataset."""
        mask = (data.positions_bp >= start_bp) & (data.positions_bp <= end_bp)
        idx = np.where(mask)[0]
        if len(idx) == 0:
            return data

        cfg = STATE_CONFIGS[state]
        if data.true_states is not None:
            data.true_states[idx] = int(state)

        # Generate LRR under state parameters
        data.lrr[idx] = self.rng.normal(cfg.lrr_mean, cfg.lrr_std * noise_level, size=len(idx))

        # Generate BAF under state peaks
        if len(cfg.baf_peaks) == 0:
            # Nullisomy: uniform BAF
            data.baf[idx] = self.rng.uniform(0.0, 1.0, size=len(idx))
        else:
            chosen_peaks = self.rng.choice(cfg.baf_peaks, size=len(idx))
            baf_noise = self.rng.normal(0.0, 0.035 * noise_level, size=len(idx))
            data.baf[idx] = np.clip(chosen_peaks + baf_noise, 0.0, 1.0)

        return data
