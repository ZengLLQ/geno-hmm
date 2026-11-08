"""Unit tests for emission model evaluation."""

import numpy as np
from geno_hmm.models.states import CNVState
from geno_hmm.models.emissions import JointEmissionModel


def test_emission_log_pdf() -> None:
    model = JointEmissionModel()
    lrr = np.array([0.0, -0.55, 0.40])
    baf = np.array([0.5, 0.0, 0.33])

    mat = model.compute_log_emission_matrix(lrr, baf)
    assert mat.shape == (3, 6)
    assert not np.isnan(mat).any()
    assert not np.isinf(mat).any()

    # Normal state should have highest likelihood for LRR=0.0, BAF=0.5
    assert np.argmax(mat[0]) == int(CNVState.NORMAL)
