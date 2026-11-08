"""Genomic distance-dependent transition matrix modeling probe spacing decay."""

from __future__ import annotations

import numpy as np

from geno_hmm.models.states import CNVState, NUM_STATES


DEFAULT_STATIONARY_PRIORS = np.array([
    1e-5,   # HOM_DEL (very rare)
    1e-3,   # HET_DEL (rare)
    0.995,  # NORMAL (common baseline)
    2e-4,   # CN_LOH (rare uniparental disomy)
    1e-3,   # SINGLE_DUP (rare)
    1e-5,   # DOUBLE_DUP (very rare)
], dtype=np.float64)
# Normalize priors
DEFAULT_STATIONARY_PRIORS /= np.sum(DEFAULT_STATIONARY_PRIORS)


class GenomicTransitionModel:
    """Distance-dependent transition matrix model.

    As physical genomic distance d between adjacent markers increases, the state transition
    relaxes exponentially from the identity matrix towards the stationary prior distribution:

        T_{s, s'}(d) = (1 - exp(-d / D)) * pi_{s'} + exp(-d / D) * delta_{s, s'}

    Parameters
    ----------
    characteristic_distance : float
        Genomic distance parameter D in base pairs (default 100,000 bp).
    stationary_priors : np.ndarray | None
        Stationary probability distribution across 6 states.
    """

    def __init__(
        self,
        characteristic_distance: float = 100_000.0,
        stationary_priors: np.ndarray | None = None,
    ) -> None:
        self.char_dist = float(characteristic_distance)
        if stationary_priors is None:
            self.priors = DEFAULT_STATIONARY_PRIORS.copy()
        else:
            p = np.asarray(stationary_priors, dtype=np.float64)
            self.priors = p / np.sum(p)

    def transition_matrix(self, distance_bp: float) -> np.ndarray:
        """Compute 6x6 transition probability matrix for a given physical distance in base pairs.

        Returns
        -------
        T : np.ndarray
            (6, 6) transition matrix where T[s, s'] = P(state_{i+1} = s' | state_i = s).
        """
        d = max(1.0, float(distance_bp))
        decay = np.exp(-d / self.char_dist)
        # Identity component + prior relaxation
        identity = np.eye(NUM_STATES, dtype=np.float64)
        priors_row = self.priors[np.newaxis, :]  # (1, 6)
        matrix = decay * identity + (1.0 - decay) * priors_row
        # Row-normalize to guarantee numerical sum to 1.0
        return matrix / np.sum(matrix, axis=1, keepdims=True)

    def log_transition_matrix(self, distance_bp: float) -> np.ndarray:
        """Compute log 6x6 transition matrix."""
        t_mat = self.transition_matrix(distance_bp)
        return np.log(np.maximum(1e-30, t_mat))
