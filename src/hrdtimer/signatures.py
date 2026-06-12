"""Signature fitting: SigProfilerAssignment baseline + MuSiCal refit.

Port of ``run_signature_analysis`` / ``calculate_exposures`` / ``save_probabilities``
from ``runSignatureAnalysis.py``, reorganised around a :class:`SignatureAnalyzer`
that holds the genome build / tumor type and lazily builds a single
:class:`~hrdtimer.catalog.SignatureCatalog`.
"""

from __future__ import annotations

import os
import sys

import pandas as pd

from .catalog import SignatureCatalog
from .constants import CLONALITY_CLASSES, DEFAULT_TUMOR_TYPE

__all__ = ["SignatureAnalyzer"]


class SignatureAnalyzer:
    """Run COSMIC signature assignment and MuSiCal refitting over a cohort tree."""

    def __init__(
        self, genome_build: str, tumor_type: str = DEFAULT_TUMOR_TYPE
    ) -> None:
        self.genome_build = genome_build
        self.tumor_type = tumor_type
        self._catalog: SignatureCatalog | None = None

    @property
    def catalog(self) -> SignatureCatalog:
        if self._catalog is None:
            self._catalog = SignatureCatalog(tumor_type=self.tumor_type)
        return self._catalog

    def run(self, parent_folder: str) -> None:
        """For each clonality subfolder: assign signatures, refit, save probs."""
        from SigProfilerAssignment import Analyzer as Analyze

        for subfolder in CLONALITY_CLASSES:
            subfolder_path = os.path.join(parent_folder, subfolder)
            if not os.path.exists(subfolder_path):
                continue

            alexandrov_output = os.path.join(subfolder_path, "Alexandrov")
            Analyze.cosmic_fit(
                subfolder_path,
                alexandrov_output,
                input_type="vcf",
                context_type="96",
                collapse_to_SBS96=False,
                cosmic_version=3.2,
                exome=False,
                genome_build=self.genome_build,
                make_plots=False,
                verbose=False,
            )

            activities_path = os.path.join(
                alexandrov_output,
                "Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
            )
            mutation_matrix_path = os.path.join(
                subfolder_path, "output/SBS/Input_vcffiles.SBS96.all"
            )
            musical_out_dir = os.path.join(parent_folder, "SignatureFitting", subfolder)

            if not os.path.exists(activities_path):
                sys.stdout.write(
                    f"Warning: SigProfiler output not found for {subfolder}. "
                    "Skipping MuSiCal.\n"
                )
                continue

            h_norm, w_reduced = self._calculate_exposures(
                activities_path,
                mutation_matrix_path,
                os.path.join(musical_out_dir, "exposures"),
            )
            self._save_probabilities(
                h_norm, w_reduced, os.path.join(musical_out_dir, "probabilities")
            )

    def _calculate_exposures(
        self, sigprofiler_activities_file: str, mutation_file: str, output_folder: str
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Refit with MuSiCal and persist SigProfiler/MuSiCal comparison matrices."""
        exposures = pd.read_csv(sigprofiler_activities_file, delimiter="\t", index_col=0)
        colnames = exposures.index.tolist()
        exposures = exposures.T
        exposures.columns = colnames

        mutation_matrix = pd.read_csv(mutation_file, delimiter="\t", index_col=0)
        numeric_cols = mutation_matrix.select_dtypes(include=["number"]).columns
        keep = [c for c in numeric_cols if not (mutation_matrix[c] == 0).all()]
        mutation_matrix = mutation_matrix[keep]

        h_naive, _ = self.catalog.refit(mutation_matrix, method="thresh_naive", thresh=0)
        h_naive.columns = colnames
        h_lik, model = self.catalog.refit(mutation_matrix)
        h_lik.columns = colnames

        os.makedirs(output_folder, exist_ok=True)
        exposures.to_csv(os.path.join(output_folder, "exposures_SigProfiler.csv"))
        h_naive.to_csv(os.path.join(output_folder, "exposures_MuSiCal_naive.csv"))
        h_lik.to_csv(os.path.join(output_folder, "exposures_MuSiCal_lik.csv"))
        model.H_reduced_normalized.to_csv(
            os.path.join(output_folder, "exposures_MuSiCal_lik_reduced_norm.csv")
        )
        model.W_reduced.to_csv(
            os.path.join(output_folder, "W_MuSiCal_lik_reduced_norm.csv")
        )
        return model.H_reduced_normalized, model.W_reduced

    @staticmethod
    def _save_probabilities(
        h_reduced_normalized: pd.DataFrame, w_reduced: pd.DataFrame, output_dir: str
    ) -> None:
        """Write per-sample mutation attribution probability matrices."""
        os.makedirs(output_dir, exist_ok=True)
        for col_name in h_reduced_normalized.columns:
            weighted = w_reduced.mul(h_reduced_normalized[col_name], axis=1)
            probabilities = weighted.div(weighted.sum(axis=1), axis=0)
            probabilities.to_csv(
                os.path.join(output_dir, f"prob_{col_name}.txt"), sep="\t", index=True
            )
