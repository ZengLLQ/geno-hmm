"""Forward-backward algorithm for computing posterior marginal state probabilities."""

from __future__ import annotations

import numpy as np

from geno_hmm.models.states import CNVState, NUM_STATES
from geno_hmm.models.transitions import GenomicTransitionModel
from geno_hmm.models.emissions import JointEmissionModel


from scipy.special import logsumexp


class ForwardBackwardAlgorithm:
    """Computes posterior marginal probabilities for each genomic probe locus."""

    def __init__(
        self,
        emission_model: JointEmissionModel | None = None,
        transition_model: GenomicTransitionModel | None = None,
    ) -> None:
        self.emission_model = emission_model or JointEmissionModel()
        self.transition_model = transition_model or GenomicTransitionModel()

    def compute_posteriors(
        self,
        positions_bp: np.ndarray,
        lrr: np.ndarray,
        baf: np.ndarray,
    ) -> tuple[np.ndarray, float]:
        """Compute gamma_t(s) = P(state_t = s | Y) and overall sequence log-likelihood.

        Parameters
        ----------
        positions_bp : np.ndarray
            Sorted (N,) genomic probe positions.
        lrr : np.ndarray
            (N,) LRR observations.
        baf : np.ndarray
            (N,) BAF observations.

        Returns
        -------
        posteriors : np.ndarray
            (N, NUM_STATES) normalized posterior probability distribution.
        total_log_likelihood : float
            Total evidence log P(Y).
        """
        N = len(positions_bp)
        if N == 0:
            return np.zeros((0, NUM_STATES)), 0.0

        log_emissions = self.emission_model.compute_log_emission_matrix(lrr, baf)
        log_priors = np.log(self.transition_model.priors)
        distances = np.diff(positions_bp) if N > 1 else np.array([])

        # Forward pass (alpha) in log space
        log_alpha = np.zeros((N, NUM_STATES), dtype=np.float64)
        log_alpha[0] = log_priors + log_emissions[0]

        for t in range(1, N):
            d = distances[t - 1]
            log_trans = self.transition_model.log_transition_matrix(d)
            # log_alpha[t, s] = logsumexp_k(log_alpha[t-1, k] + log_trans[k, s]) + log_emissions[t, s]
            prev_plus_trans = log_alpha[t - 1, :, np.newaxis] + log_trans
            summed_prev = logsumexp(prev_plus_trans, axis=0)
            log_alpha[t] = summed_prev + log_emissions[t]

        total_log_likelihood = float(logsumexp(log_alpha[N - 1]))

        # Backward pass (beta) in log space
        log_beta = np.zeros((N, NUM_STATES), dtype=np.float64)
        log_beta[N - 1] = 0.0  # log(1.0) = 0.0

        for t in range(N - 2, -1, -1):
            d = distances[t]
            log_trans = self.transition_model.log_transition_matrix(d)
            # log_beta[t, s] = logsumexp_k(log_trans[s, k] + log_emissions[t+1, k] + log_beta[t+1, k])
            trans_emiss_beta = log_trans + log_emissions[t + 1, np.newaxis, :] + log_beta[t + 1, np.newaxis, :]
            log_beta[t] = logsumexp(trans_emiss_beta, axis=1)

        # Posterior marginals: log_gamma = log_alpha + log_beta - log P(Y)
        log_gamma = log_alpha + log_beta - total_log_likelihood
        posteriors = np.exp(np.clip(log_gamma, -700.0, 0.0))
        # Ensure row-wise normalization
        posteriors /= np.sum(posteriors, axis=1, keepdims=True)

        return posteriors, total_log_likelihood
