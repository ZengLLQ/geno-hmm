"""ACMG/ClinGen Technical Standards quantitative classification rubric for CNVs."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from geno_hmm.models.states import CNVState
from geno_hmm.algorithms.viterbi import CNVSegment
from geno_hmm.clinical.dosage_curated import DosageGene, find_overlapping_dosage_genes


class ACMGClass(str, Enum):
    """ACMG/ClinGen 5-tier clinical classification."""

    PATHOGENIC = "Pathogenic (Class 5)"
    LIKELY_PATHOGENIC = "Likely Pathogenic (Class 4)"
    VUS = "Variant of Uncertain Significance (Class 3)"
    LIKELY_BENIGN = "Likely Benign (Class 2)"
    BENIGN = "Benign (Class 1)"


@dataclass
class EvidenceItem:
    """Individual evidence item contributing to ACMG score."""

    section: str
    code: str
    description: str
    points: float


@dataclass
class ACMGClassificationReport:
    """Detailed clinical interpretation and evidence dossier for a CNV segment."""

    segment: CNVSegment
    total_score: float
    acmg_class: ACMGClass
    dosage_genes: list[DosageGene]
    evidence_breakdown: list[EvidenceItem]
    clinical_summary: str


class ACMGScorer:
    """Evaluates CNVs under the ACMG/ClinGen quantitative point rubric (Riggs et al., 2020)."""

    def classify_segment(
        self,
        segment: CNVSegment,
        is_de_novo: bool = False,
        phenotype_match: bool = False,
    ) -> ACMGClassificationReport:
        """Score a genomic segment and determine pathogenicity classification.

        Parameters
        ----------
        segment : CNVSegment
            The identified CNV event.
        is_de_novo : bool
            Whether pedigree confirms de novo occurrence in patient (Section 4).
        phenotype_match : bool
            Whether patient's clinical indications match known syndrome (Section 4).

        Returns
        -------
        report : ACMGClassificationReport
            Complete clinical dossier.
        """
        evidence: list[EvidenceItem] = []
        is_loss = segment.state in (CNVState.HOM_DEL, CNVState.HET_DEL)
        is_gain = segment.state in (CNVState.SINGLE_DUP, CNVState.DOUBLE_DUP)
        genes = find_overlapping_dosage_genes(segment.chromosome, segment.start_bp, segment.end_bp)

        # 1. Section 1 / 2: Genomic content and Dosage Sensitivity Overlap
        if is_loss:
            hi_genes = [g for g in genes if g.haploinsufficiency_score == 3]
            if hi_genes:
                evidence.append(
                    EvidenceItem(
                        section="Section 2",
                        code="2A",
                        description=f"Complete/partial overlap with established haploinsufficient gene(s): {', '.join(g.symbol for g in hi_genes)} (HI=3)",
                        points=1.00,
                    )
                )
            elif len(genes) > 0:
                evidence.append(
                    EvidenceItem(
                        section="Section 2",
                        code="2C",
                        description=f"Overlaps gene(s) without established HI evidence: {', '.join(g.symbol for g in genes)}",
                        points=0.00,
                    )
                )
            else:
                if segment.length_bp >= 5_000_000:
                    evidence.append(
                        EvidenceItem(
                            section="Section 1",
                            code="1B",
                            description="Large gene-desert microdeletion >= 5 Mb",
                            points=0.45,
                        )
                    )
                else:
                    evidence.append(
                        EvidenceItem(
                            section="Section 1",
                            code="1A",
                            description="No dosage-sensitive genes involved in deletion",
                            points=0.00,
                        )
                    )

        elif is_gain:
            ts_genes = [g for g in genes if g.triplosensitivity_score == 3]
            if ts_genes:
                evidence.append(
                    EvidenceItem(
                        section="Section 2",
                        code="2A_dup",
                        description=f"Overlaps established triplosensitive gene(s): {', '.join(g.symbol for g in ts_genes)} (TS=3)",
                        points=1.00,
                    )
                )
            elif len(genes) > 0:
                evidence.append(
                    EvidenceItem(
                        section="Section 2",
                        code="2C_dup",
                        description=f"Overlaps gene(s) with benign duplication tolerance or unknown TS: {', '.join(g.symbol for g in genes)}",
                        points=0.00,
                    )
                )
            else:
                evidence.append(
                    EvidenceItem(
                        section="Section 1",
                        code="1A_dup",
                        description="Duplication does not encompass established triplosensitive loci",
                        points=0.00,
                    )
                )
        else:
            # CN-LOH
            evidence.append(
                EvidenceItem(
                    section="Section 1",
                    code="LOH",
                    description="Copy-neutral LOH / Uniparental Isodisomy candidate",
                    points=0.30 if segment.length_bp >= 10_000_000 else 0.00,
                )
            )

        # 2. Section 3: Genomic size and probe density
        if segment.length_bp >= 3_000_000:
            evidence.append(
                EvidenceItem(
                    section="Section 3",
                    code="3C",
                    description=f"Large chromosomal segment (span = {segment.length_bp / 1e6:.2f} Mb, {segment.num_probes} probes)",
                    points=0.45,
                )
            )
        elif segment.length_bp >= 1_000_000:
            evidence.append(
                EvidenceItem(
                    section="Section 3",
                    code="3B",
                    description=f"Moderate chromosomal segment ({segment.length_bp / 1e6:.2f} Mb)",
                    points=0.20,
                )
            )

        # 3. Section 4: Detailed Case Evidence & Inheritance
        if is_de_novo:
            evidence.append(
                EvidenceItem(
                    section="Section 4",
                    code="4B",
                    description="Confirmed de novo event (biparental origin verified)",
                    points=0.45,
                )
            )
        if phenotype_match:
            evidence.append(
                EvidenceItem(
                    section="Section 4",
                    code="4E",
                    description="Phenotype highly specific and consistent with known microdeletion/microduplication syndrome",
                    points=0.30,
                )
            )

        # Total points aggregation
        total_score = sum(e.points for e in evidence)

        # Determine 5-tier classification
        if total_score >= 0.99:
            acmg_class = ACMGClass.PATHOGENIC
        elif total_score >= 0.90:
            acmg_class = ACMGClass.LIKELY_PATHOGENIC
        elif total_score >= -0.89:
            acmg_class = ACMGClass.VUS
        elif total_score >= -0.98:
            acmg_class = ACMGClass.LIKELY_BENIGN
        else:
            acmg_class = ACMGClass.BENIGN

        # Summary text
        genes_str = ", ".join(g.symbol for g in genes) if genes else "none"
        summary = (
            f"Genomic {segment.state_name} on {segment.chromosome}:{segment.start_bp:,}-{segment.end_bp:,} "
            f"({segment.length_bp / 1e3:.1f} kb, {segment.num_probes} probes). "
            f"Overlapping genes: {genes_str}. Total ACMG score = {total_score:+.2f} points -> {acmg_class.value}."
        )

        return ACMGClassificationReport(
            segment=segment,
            total_score=total_score,
            acmg_class=acmg_class,
            dosage_genes=genes,
            evidence_breakdown=evidence,
            clinical_summary=summary,
        )
