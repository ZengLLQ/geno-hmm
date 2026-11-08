"""Example: Detection and ACMG classification of de novo 22q11.2 DiGeorge microdeletion."""

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
    print("Simulating patient chromosome 22 with 22q11.2 microdeletion...")
    sim = MicroarraySimulator(random_seed=123)
    data = sim.generate_chromosome(
        chromosome="chr22",
        start_bp=18_000_000,
        end_bp=22_000_000,
        avg_probe_spacing_bp=2_000,
        noise_level=0.9,
    )
    # Inject 22q11.2 classic 1.5 Mb deletion encompassing TBX1
    data = sim.inject_alteration(
        data,
        start_bp=19_200_000,
        end_bp=20_700_000,
        state=CNVState.HET_DEL,
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

    print("\nClinical Interpretation:")
    scorer = ACMGScorer()
    for seg in segments:
        report = scorer.classify_segment(seg, is_de_novo=True, phenotype_match=True)
        print(f"Outcome: {report.acmg_class.value}")
        print(f"Summary: {report.clinical_summary}")


if __name__ == "__main__":
    main()
