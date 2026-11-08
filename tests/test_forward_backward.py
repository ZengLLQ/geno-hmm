"""Unit tests for Forward-Backward marginal probability computation."""

import numpy as np
from geno_hmm.models.states import CNVState
from geno_hmm.algorithms.forward_backward import ForwardBackwardAlgorithm


def test_forward_backward_posteriors() -> None:
    fb = ForwardBackwardAlgorithm()
    pos = np.array([1000, 2000, 3000, 4000], dtype=np.int64)
    lrr = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    baf = np.array([0.5, 0.5, 0.5, 0.5], dtype=np.float64)

    posteriors, log_lik = fb.compute_posteriors(pos, lrr, baf)
    assert posteriors.shape == (4, 6)
    # Check probabilities sum to 1.0
    assert np.allclose(np.sum(posteriors, axis=1), 1.0)
    # Normal state should dominate
    assert np.all(posteriors[:, int(CNVState.NORMAL)] > 0.8)
    assert not np.isnan(log_lik)
