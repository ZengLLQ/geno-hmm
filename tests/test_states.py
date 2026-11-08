"""Unit tests for CNV states and configuration parameters."""

from geno_hmm.models.states import CNVState, NUM_STATES, STATE_CONFIGS


def test_num_states() -> None:
    assert NUM_STATES == 6
    assert len(CNVState) == 6
    assert len(STATE_CONFIGS) == 6


def test_state_configs() -> None:
    normal = STATE_CONFIGS[CNVState.NORMAL]
    assert normal.copy_number == 2
    assert normal.lrr_mean == 0.0
    assert 0.50 in normal.baf_peaks

    hom_del = STATE_CONFIGS[CNVState.HOM_DEL]
    assert hom_del.copy_number == 0
    assert hom_del.lrr_mean < -2.0
    assert len(hom_del.baf_peaks) == 0

    het_del = STATE_CONFIGS[CNVState.HET_DEL]
    assert het_del.copy_number == 1
    assert het_del.lrr_mean < 0.0
    assert het_del.baf_peaks == (0.0, 1.0)
