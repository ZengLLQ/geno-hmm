"""Unit tests for ACMG/ClinGen classification logic."""

from geno_hmm.models.states import CNVState
from geno_hmm.algorithms.viterbi import CNVSegment
from geno_hmm.clinical.dosage_curated import find_overlapping_dosage_genes
from geno_hmm.clinical.acmg_scorer import ACMGScorer, ACMGClass


def test_dosage_overlap() -> None:
    # 22q11.2 deletion overlaps TBX1
    genes = find_overlapping_dosage_genes("chr22", 19_000_000, 20_000_000)
    symbols = [g.symbol for g in genes]
    assert "TBX1" in symbols


def test_acmg_pathogenic_classification() -> None:
    seg = CNVSegment(
        chromosome="chr22",
        start_bp=19_000_000,
        end_bp=21_500_000,
        length_bp=2_500_000,
        num_probes=1000,
        state=CNVState.HET_DEL,
        copy_number=1,
        state_name="Hemizygous Deletion",
        mean_lrr=-0.52,
        median_baf=0.02,
        confidence_score=45.0,
    )

    scorer = ACMGScorer()
    report = scorer.classify_segment(seg, is_de_novo=True, phenotype_match=True)
    assert report.acmg_class == ACMGClass.PATHOGENIC
    assert report.total_score >= 0.99
    assert any(g.symbol == "TBX1" for g in report.dosage_genes)
