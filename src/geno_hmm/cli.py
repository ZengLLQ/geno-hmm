"""Command line interface for GenoHMM clinical CNV detection and interpretation."""

from __future__ import annotations

import argparse
import sys
import numpy as np

from geno_hmm.models.states import CNVState
from geno_hmm.simulation.array_simulator import MicroarraySimulator
from geno_hmm.algorithms.viterbi import ViterbiDecoder
from geno_hmm.clinical.acmg_scorer import ACMGScorer
from geno_hmm.visualizer import render_ascii_profile, format_clinical_table


def handle_demo(args: argparse.Namespace) -> int:
    """Run demonstration with synthetic microdeletion or duplication."""
    print("=" * 78)
    print("  GenoHMM: Clinical CNV Detection & ACMG Classification Demonstration")
    print("=" * 78)

    sim = MicroarraySimulator(random_seed=args.seed)
    scenario = args.scenario.lower()

    if scenario == "digeorge":
        print("\n[Scenario] Simulating Chromosome 22 with 22q11.2 Microdeletion (TBX1 locus)...")
        data = sim.generate_chromosome(
            chromosome="chr22",
            start_bp=18_000_000,
            end_bp=22_000_000,
            avg_probe_spacing_bp=2_500,
            noise_level=args.noise,
        )
        # Inject hemizygous deletion (22q11.2 classic 3Mb or 1.5Mb region)
        data = sim.inject_alteration(
            data,
            start_bp=19_000_000,
            end_bp=20_500_000,
            state=CNVState.HET_DEL,
            noise_level=args.noise,
        )
    elif scenario == "her2":
        print("\n[Scenario] Simulating Chromosome 17 with Focal HER2 / ERBB2 Amplification...")
        data = sim.generate_chromosome(
            chromosome="chr17",
            start_bp=37_000_000,
            end_bp=42_000_000,
            avg_probe_spacing_bp=2_500,
            noise_level=args.noise,
        )
        data = sim.inject_alteration(
            data,
            start_bp=39_500_000,
            end_bp=40_200_000,
            state=CNVState.DOUBLE_DUP,
            noise_level=args.noise,
        )
    else:
        print(f"Unknown scenario: {args.scenario}")
        return 1

    print(f"Synthesized {len(data.positions_bp):,} microarray SNP probes.")
    print("Running 6-State Genomic Viterbi HMM Decoder...")

    decoder = ViterbiDecoder()
    states = decoder.decode(data.positions_bp, data.lrr, data.baf)
    segments = decoder.call_segments(
        chromosome=data.chromosome,
        positions_bp=data.positions_bp,
        lrr=data.lrr,
        baf=data.baf,
        min_probes=5,
    )

    print("\n" + render_ascii_profile(data.positions_bp, data.lrr, data.baf, states))
    print("\n" + format_clinical_table(segments))

    # ACMG Clinical Classification
    scorer = ACMGScorer()
    print("\n" + "=" * 78)
    print("  ACMG / ClinGen Technical Standard Classification Dossier")
    print("=" * 78)

    for i, seg in enumerate(segments, 1):
        report = scorer.classify_segment(seg, is_de_novo=True, phenotype_match=True)
        print(f"\n[CNV #{i}] {report.clinical_summary}")
        print(f"  Classification: {report.acmg_class.value}")
        print(f"  Total Score:    {report.total_score:+.2f} points")
        print("  Evidence Items:")
        for ev in report.evidence_breakdown:
            print(f"    - [{ev.section} / {ev.code}] {ev.description} ({ev.points:+.2f} pts)")

    return 0


def handle_benchmark(args: argparse.Namespace) -> int:
    """Run sensitivity and specificity benchmark over simulated CNVs."""
    print("=" * 78)
    print(f"  GenoHMM Benchmark: Evaluating on {args.num_runs} Injected Chromosomes")
    print("=" * 78)

    sim = MicroarraySimulator(random_seed=42)
    decoder = ViterbiDecoder()

    tp_probes = 0
    fp_probes = 0
    fn_probes = 0
    tn_probes = 0

    for run in range(args.num_runs):
        data = sim.generate_chromosome(
            chromosome="chr22",
            start_bp=18_000_000,
            end_bp=22_000_000,
            avg_probe_spacing_bp=3_000,
            noise_level=1.0,
        )
        # Random CNV location
        c_start = 19_000_000 + (run * 200_000) % 1_500_000
        c_end = c_start + 500_000
        st = CNVState.HET_DEL if run % 2 == 0 else CNVState.SINGLE_DUP
        data = sim.inject_alteration(data, c_start, c_end, st)

        pred_states = decoder.decode(data.positions_bp, data.lrr, data.baf)
        true_alt = data.true_states != int(CNVState.NORMAL)
        pred_alt = pred_states != int(CNVState.NORMAL)

        tp_probes += int(np.sum(true_alt & pred_alt))
        fp_probes += int(np.sum((~true_alt) & pred_alt))
        fn_probes += int(np.sum(true_alt & (~pred_alt)))
        tn_probes += int(np.sum((~true_alt) & (~pred_alt)))

    sens = tp_probes / max(1, (tp_probes + fn_probes))
    spec = tn_probes / max(1, (tn_probes + fp_probes))
    f1 = 2 * tp_probes / max(1, (2 * tp_probes + fp_probes + fn_probes))

    print(f"\nProbe-Level Results across {args.num_runs} runs:")
    print(f"  True Positives:  {tp_probes:,}")
    print(f"  False Positives: {fp_probes:,}")
    print(f"  False Negatives: {fn_probes:,}")
    print(f"  True Negatives:  {tn_probes:,}")
    print(f"  Sensitivity (Recall): {sens * 100:.2f}%")
    print(f"  Specificity:          {spec * 100:.2f}%")
    print(f"  F1 Score:             {f1:.4f}")

    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI main entry point."""
    parser = argparse.ArgumentParser(
        prog="geno-hmm",
        description="Clinical Copy Number Variation (CNV) Detection and ACMG Pathogenicity Scoring",
    )
    subparsers = parser.add_subparsers(dest="command", help="Subcommand to run")

    # Demo
    demo_parser = subparsers.add_parser("demo", help="Run interactive CNV demonstration")
    demo_parser.add_argument(
        "--scenario",
        choices=["digeorge", "her2"],
        default="digeorge",
        help="Simulated clinical syndrome (default: digeorge)",
    )
    demo_parser.add_argument("--noise", type=float, default=1.0, help="Array noise level multiplier")
    demo_parser.add_argument("--seed", type=int, default=42, help="Random number generator seed")

    # Benchmark
    bench_parser = subparsers.add_parser("benchmark", help="Run accuracy benchmarks")
    bench_parser.add_argument("--num-runs", type=int, default=10, help="Number of simulated chromosomes")

    args = parser.parse_args(argv)
    if args.command == "demo":
        return handle_demo(args)
    elif args.command == "benchmark":
        return handle_benchmark(args)
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
