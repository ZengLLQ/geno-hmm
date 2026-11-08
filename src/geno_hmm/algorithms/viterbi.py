"""Log-space Viterbi dynamic programming decoder and CNV segment extraction."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from geno_hmm.models.states import CNVState, NUM_STATES, STATE_CONFIGS
from geno_hmm.models.transitions import GenomicTransitionModel
from geno_hmm.models.emissions import JointEmissionModel


@dataclass
class CNVSegment:
    """A contiguous genomic copy number alteration call."""

    chromosome: str
    start_bp: int
    end_bp: int
    length_bp: int
    num_probes: int
    state: CNVState
    copy_number: int
    state_name: str
    mean_lrr: float
    median_baf: float
    confidence_score: float  # Log-odds improvement over diploid normal baseline


class ViterbiDecoder:
    """Dynamic programming Viterbi decoder operating strictly in log-space."""

    def __init__(
        self,
        emission_model: JointEmissionModel | None = None,
        transition_model: GenomicTransitionModel | None = None,
    ) -> None:
        self.emission_model = emission_model or JointEmissionModel()
        self.transition_model = transition_model or GenomicTransitionModel()

    def decode(
        self,
        positions_bp: np.ndarray,
        lrr: np.ndarray,
        baf: np.ndarray,
    ) -> np.ndarray:
        """Decode the maximum a posteriori (MAP) sequence of copy number states.

        Parameters
        ----------
        positions_bp : np.ndarray
            Sorted (N,) 1D array of probe chromosomal base pair coordinates.
        lrr : np.ndarray
            (N,) array of Log R Ratio values.
        baf : np.ndarray
            (N,) array of B-Allele Frequency values.

        Returns
        -------
        best_path : np.ndarray
            (N,) array of integer states (0..5).
        """
        N = len(positions_bp)
        if N == 0:
            return np.array([], dtype=int)
        if N == 1:
            log_emissions = self.emission_model.compute_log_emission_matrix(lrr, baf)
            log_priors = np.log(self.transition_model.priors)
            return np.array([int(np.argmax(log_emissions[0] + log_priors))])

        log_emissions = self.emission_model.compute_log_emission_matrix(lrr, baf)
        log_priors = np.log(self.transition_model.priors)

        # DP trellis matrices
        viterbi_trellis = np.zeros((N, NUM_STATES), dtype=np.float64)
        backpointers = np.zeros((N, NUM_STATES), dtype=np.int32)

        # Initialization at t = 0
        viterbi_trellis[0] = log_priors + log_emissions[0]

        # Forward DP recursion
        distances = np.diff(positions_bp)
        for t in range(1, N):
            d = distances[t - 1]
            log_trans = self.transition_model.log_transition_matrix(d)
            # For each state s at step t:
            # viterbi[t, s] = max_k (viterbi[t-1, k] + log_trans[k, s]) + log_emissions[t, s]
            # Vectorized across all previous states k (shape: 6, 6)
            # prev (6, 1) + trans (6, 6) -> max along axis 0
            prev_plus_trans = viterbi_trellis[t - 1, :, np.newaxis] + log_trans
            best_prev = np.argmax(prev_plus_trans, axis=0)
            best_vals = np.max(prev_plus_trans, axis=0)

            viterbi_trellis[t] = best_vals + log_emissions[t]
            backpointers[t] = best_prev

        # Backtracking
        best_path = np.zeros(N, dtype=np.int32)
        best_path[N - 1] = int(np.argmax(viterbi_trellis[N - 1]))

        for t in range(N - 2, -1, -1):
            best_path[t] = backpointers[t + 1, best_path[t + 1]]

        return best_path

    def call_segments(
        self,
        chromosome: str,
        positions_bp: np.ndarray,
        lrr: np.ndarray,
        baf: np.ndarray,
        min_probes: int = 3,
        min_length_bp: int = 1_000,
    ) -> list[CNVSegment]:
        """Decode path and aggregate contiguous CNV segments differing from NORMAL (State 2).

        Parameters
        ----------
        chromosome : str
            Chromosome identifier (e.g. "chr22").
        positions_bp : np.ndarray
            Marker coordinates.
        lrr : np.ndarray
            LRR measurements.
        baf : np.ndarray
            BAF measurements.
        min_probes : int
            Minimum number of contiguous probes to report a CNV call.
        min_length_bp : int
            Minimum genomic span in base pairs.

        Returns
        -------
        segments : list[CNVSegment]
            Filtered list of high-confidence genomic alteration segments.
        """
        states = self.decode(positions_bp, lrr, baf)
        N = len(states)
        if N == 0:
            return []

        log_emissions = self.emission_model.compute_log_emission_matrix(lrr, baf)
        normal_state_idx = int(CNVState.NORMAL)

        segments: list[CNVSegment] = []
        seg_start = 0

        while seg_start < N:
            cur_state = CNVState(states[seg_start])
            seg_end = seg_start
            while seg_end + 1 < N and states[seg_end + 1] == cur_state:
                seg_end += 1

            # Only report non-normal states (or CN_LOH)
            if cur_state != CNVState.NORMAL:
                num_probes = seg_end - seg_start + 1
                start_coord = int(positions_bp[seg_start])
                end_coord = int(positions_bp[seg_end])
                length = end_coord - start_coord

                if num_probes >= min_probes and length >= min_length_bp:
                    seg_lrr = lrr[seg_start : seg_end + 1]
                    seg_baf = baf[seg_start : seg_end + 1]
                    cfg = STATE_CONFIGS[cur_state]

                    # Confidence score: log-likelihood difference vs normal diploid
                    seg_emiss_call = log_emissions[seg_start : seg_end + 1, int(cur_state)]
                    seg_emiss_norm = log_emissions[seg_start : seg_end + 1, normal_state_idx]
                    log_odds_score = float(np.sum(seg_emiss_call - seg_emiss_norm))

                    segments.append(
                        CNVSegment(
                            chromosome=chromosome,
                            start_bp=start_coord,
                            end_bp=end_coord,
                            length_bp=length,
                            num_probes=num_probes,
                            state=cur_state,
                            copy_number=cfg.copy_number,
                            state_name=cfg.name,
                            mean_lrr=float(np.mean(seg_lrr)),
                            median_baf=float(np.median(seg_baf)),
                            confidence_score=max(0.1, log_odds_score),
                        )
                    )

            seg_start = seg_end + 1

        return segments
