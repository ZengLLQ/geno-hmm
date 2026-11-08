"""ClinGen curated dosage sensitivity map and reference gene database."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DosageGene:
    """ClinGen curated gene with dosage sensitivity scores."""

    symbol: str
    chromosome: str
    start_bp: int
    end_bp: int
    haploinsufficiency_score: int  # 0: no evidence, 1: emerging, 2: some, 3: sufficient evidence
    triplosensitivity_score: int   # 0: no evidence, 1: emerging, 2: some, 3: sufficient evidence
    associated_disorder: str
    omim_id: str


# High-impact curated dosage sensitive genes across human genome
CLINGEN_DOSAGE_GENES: list[DosageGene] = [
    DosageGene(
        symbol="TBX1",
        chromosome="chr22",
        start_bp=19_750_000,
        end_bp=19_780_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=0,
        associated_disorder="22q11.2 Deletion / DiGeorge Syndrome",
        omim_id="188400",
    ),
    DosageGene(
        symbol="ELN",
        chromosome="chr7",
        start_bp=73_440_000,
        end_bp=73_485_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=0,
        associated_disorder="Williams-Beuren Syndrome",
        omim_id="194050",
    ),
    DosageGene(
        symbol="SHANK3",
        chromosome="chr22",
        start_bp=51_110_000,
        end_bp=51_170_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=0,
        associated_disorder="Phelan-McDermid Syndrome / 22q13.3",
        omim_id="606230",
    ),
    DosageGene(
        symbol="PMP22",
        chromosome="chr17",
        start_bp=15_229_000,
        end_bp=15_265_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=3,
        associated_disorder="HNPP (del) / CMT1A (dup)",
        omim_id="118220",
    ),
    DosageGene(
        symbol="MECP2",
        chromosome="chrX",
        start_bp=153_287_000,
        end_bp=153_363_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=3,
        associated_disorder="Rett Syndrome (del) / MECP2 Duplication",
        omim_id="312750",
    ),
    DosageGene(
        symbol="PTEN",
        chromosome="chr10",
        start_bp=89_623_000,
        end_bp=89_731_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=0,
        associated_disorder="Cowden Syndrome / PTEN Hamartoma",
        omim_id="158350",
    ),
    DosageGene(
        symbol="TP53",
        chromosome="chr17",
        start_bp=7_668_000,
        end_bp=7_687_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=0,
        associated_disorder="Li-Fraumeni Syndrome",
        omim_id="191170",
    ),
    DosageGene(
        symbol="MYC",
        chromosome="chr8",
        start_bp=127_735_000,
        end_bp=127_742_000,
        haploinsufficiency_score=0,
        triplosensitivity_score=3,
        associated_disorder="Somatic Burkitt Lymphoma / Carcinoma Amplification",
        omim_id="190080",
    ),
    DosageGene(
        symbol="ERBB2",
        chromosome="chr17",
        start_bp=39_688_000,
        end_bp=39_728_000,
        haploinsufficiency_score=0,
        triplosensitivity_score=3,
        associated_disorder="HER2 Amplification in Breast / Gastric Carcinoma",
        omim_id="164870",
    ),
    DosageGene(
        symbol="RB1",
        chromosome="chr13",
        start_bp=48_303_000,
        end_bp=48_481_000,
        haploinsufficiency_score=3,
        triplosensitivity_score=0,
        associated_disorder="Retinoblastoma Susceptibility",
        omim_id="180200",
    ),
]


def find_overlapping_dosage_genes(
    chromosome: str,
    start_bp: int,
    end_bp: int,
) -> list[DosageGene]:
    """Find all curated dosage sensitive genes intersecting the interval [start_bp, end_bp]."""
    matches: list[DosageGene] = []
    chrom_norm = chromosome.lower()
    if not chrom_norm.startswith("chr"):
        chrom_norm = "chr" + chrom_norm

    for gene in CLINGEN_DOSAGE_GENES:
        gene_chr = gene.chromosome.lower()
        if gene_chr == chrom_norm:
            # Overlap test: max(start1, start2) <= min(end1, end2)
            if max(start_bp, gene.start_bp) <= min(end_bp, gene.end_bp):
                matches.append(gene)

    return matches
