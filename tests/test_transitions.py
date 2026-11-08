"""Unit tests for genomic distance decay transition matrix."""

import numpy as np
from geno_hmm.models.transitions import GenomicTransitionModel


def test_transition_properties() -> None:
    model = GenomicTransitionModel(characteristic_distance=100_000.0)

    # For tiny distance, should be near identity
    t_near = model.transition_matrix(10.0)
    assert np.allclose(np.diag(t_near), 1.0, atol=0.01)
    assert np.allclose(np.sum(t_near, axis=1), 1.0)

    # For huge distance, rows should relax to priors
    t_far = model.transition_matrix(10_000_000.0)
    for row in t_far:
        assert np.allclose(row, model.priors, atol=0.01)
