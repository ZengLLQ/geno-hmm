"""Unit tests for Viterbi decoding and segmentation."""

import numpy as np
from geno_hmm.models.states import CNVState
from geno_hmm.simulation.array_simulator import MicroarraySimulator
from geno_hmm.algorithms.viterbi import ViterbiDecoder


def test_viterbi_decoding_and_segmentation() -> None:
    sim = MicroarraySimulator(random_seed=42)
    data = sim.generate_chromosome(
        chromosome="chr22",
        start_bp=19_000_000,
        end_bp=21_000_000,
        avg_probe_spacing_bp=2_000,
    )
    # Inject clear hemizygous deletion
    data = sim.inject_alteration(
        data,
        start_bp=19_500_000,
        end_bp=20_500_000,
        state=CNVState.HET_DEL,
    )

    decoder = ViterbiDecoder()
    states = decoder.decode(data.positions_bp, data.lrr, data.baf)
    assert len(states) == len(data.positions_bp)

    segments = decoder.call_segments(
        chromosome="chr22",
        positions_bp=data.positions_bp,
        lrr=data.lrr,
        baf=data.baf,
        min_probes=10,
    )

    assert len(segments) >= 1
    main_seg = max(segments, key=lambda s: s.num_probes)
    assert main_seg.state == CNVState.HET_DEL
    assert main_seg.copy_number == 1
    assert main_seg.start_bp <= 19_600_000
    assert main_seg.end_bp >= 20_400_000
    assert main_seg.mean_lrr < -0.3
