"""Joint LRR and BAF emission probability distribution models."""

from __future__ import annotations

import math
import numpy as np

from geno_hmm.models.states import CNVState, NUM_STATES, STATE_CONFIGS


class JointEmissionModel:
    """Joint continuous emission probability density for Log R Ratio and B-Allele Frequency.

    Parameters
    ----------
    baf_std : float
        Standard deviation of BAF Gaussian genotype clusters (default 0.04).
    uniform_background_fraction : float
        Small uniform noise floor for artifactual probe readings (default 0.02).
    """

    def __init__(
        self,
        baf_std: float = 0.04,
        uniform_background_fraction: float = 0.02,
    ) -> None:
        self.baf_std = float(baf_std)
        self.bg_frac = float(uniform_background_fraction)
        self.configs = STATE_CONFIGS

    def evaluate_lrr_log_pdf(self, lrr: np.ndarray, state: CNVState) -> np.ndarray:
        """Evaluate Gaussian log-likelihood for LRR observations."""
        cfg = self.configs[state]
        mu = cfg.lrr_mean
        sig = cfg.lrr_std
        var = sig**2
        # -0.5 * log(2*pi*var) - (x - mu)^2 / (2*var)
        norm_const = -0.5 * (math.log(2.0 * math.pi) + 2.0 * math.log(sig))
        diff_sq = (lrr - mu) ** 2
        return norm_const - (diff_sq / (2.0 * var))

    def evaluate_baf_log_pdf(self, baf: np.ndarray, state: CNVState) -> np.ndarray:
        """Evaluate Gaussian Mixture Model log-likelihood for BAF observations."""
        cfg = self.configs[state]
        peaks = cfg.baf_peaks

        if len(peaks) == 0:
            # Nullisomy / Homozygous deletion: BAF is uniform noise in [0, 1]
            return np.zeros_like(baf)

        # GMM across genotype allelic peaks
        sig = self.baf_std
        var = sig**2
        norm_factor = 1.0 / (math.sqrt(2.0 * math.pi) * sig)
        weight_per_peak = (1.0 - self.bg_frac) / len(peaks)

        # Sum of Gaussian components
        # baf shape: (N,) -> broadcast with peaks: (N, num_peaks)
        baf_col = baf[:, np.newaxis]
        peaks_row = np.array(peaks)[np.newaxis, :]
        diffs = baf_col - peaks_row
        gaussians = norm_factor * np.exp(-0.5 * (diffs**2) / var)
        gmm_density = weight_per_peak * np.sum(gaussians, axis=1) + self.bg_frac * 1.0

        return np.log(np.maximum(1e-12, gmm_density))

    def compute_log_emission_matrix(self, lrr: np.ndarray, baf: np.ndarray) -> np.ndarray:
        """Compute complete log-emission probability matrix of shape (N, NUM_STATES).

        Parameters
        ----------
        lrr : np.ndarray
            (N,) array of Log R Ratio intensity values.
        baf : np.ndarray
            (N,) array of B-Allele Frequency values in [0, 1].

        Returns
        -------
        log_emissions : np.ndarray
            (N, NUM_STATES) matrix where element (i, s) = log P(LRR_i, BAF_i | state s).
        """
        N = len(lrr)
        log_emissions = np.zeros((N, NUM_STATES))

        for s in range(NUM_STATES):
            st = CNVState(s)
            log_p_lrr = self.evaluate_lrr_log_pdf(lrr, st)
            log_p_baf = self.evaluate_baf_log_pdf(baf, st)
            log_emissions[:, s] = log_p_lrr + log_p_baf

        return log_emissions
