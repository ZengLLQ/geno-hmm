"""Example: Detection and scoring of focal ERBB2/HER2 oncogene amplification."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from geno_hmm import (
    MicroarraySimulator,
    ViterbiDecoder,
    ACMGScorer,
    CNVState,
    render_ascii_profile,
    format_clinical_table,
)


def main() -> None:
    print("Simulating breast cancer biopsy chromosome 17 with focal ERBB2 amplification...")
    sim = MicroarraySimulator(random_seed=789)
    data = sim.generate_chromosome(
        chromosome="chr17",
        start_bp=38_000_000,
        end_bp=42_000_000,
        avg_probe_spacing_bp=2_000,
        noise_level=1.0,
    )
    # Inject focal tetrasomy (CN = 4) amplification
    data = sim.inject_alteration(
        data,
        start_bp=39_600_000,
        end_bp=40_000_000,
        state=CNVState.DOUBLE_DUP,
    )

    decoder = ViterbiDecoder()
    states = decoder.decode(data.positions_bp, data.lrr, data.baf)
    segments = decoder.call_segments(
        chromosome=data.chromosome,
        positions_bp=data.positions_bp,
        lrr=data.lrr,
        baf=data.baf,
    )

    print("\nASCII Profile:")
    print(render_ascii_profile(data.positions_bp, data.lrr, data.baf, states, width=65))

    print("\nSegment Calls:")
    print(format_clinical_table(segments))

    scorer = ACMGScorer()
    for seg in segments:
        report = scorer.classify_segment(seg)
        print(f"\nReport: {report.clinical_summary}")


if __name__ == "__main__":
    main()
