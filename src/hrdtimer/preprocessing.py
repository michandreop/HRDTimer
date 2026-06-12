"""Mutation preprocessing: split VCFs by clonality and assemble timing inputs.

Consolidates the three near-duplicate ``process_vcfs_early_late*`` functions
into the single subclonal-aware version (the one the production pipeline calls),
and keeps the vectorised ``prepare_samples_for_timing``.
"""

from __future__ import annotations

import os
import shutil
import sys

import pandas as pd
from tqdm import tqdm

from . import vcf as vcf_io
from .constants import CLONALITY_CLASSES

__all__ = ["MutationPreprocessor"]


class MutationPreprocessor:
    """Split MutationTimeR VCFs by organ and clonality, ready for signature fit."""

    def __init__(self, metadata_csv: str) -> None:
        metadata = pd.read_csv(metadata_csv)
        self._organ_lookup = metadata.set_index("sample")["organ"].to_dict()
        self._organs = metadata["organ"].unique()

    def split_by_clonality(
        self,
        input_folder: str,
        output_folder: str,
        time_analysis: bool = False,
        verbose: bool = False,
    ) -> None:
        """Categorise mutations by organ and clonality into a directory tree.

        Layout: ``output_folder/<organ>/<timing|all_mut>/<class>/<id>_<class>.vcf``.
        Empty organ trees created up front are removed at the end.
        """
        mode = "timing" if time_analysis else "all_mut"

        organ_folders = {
            organ: {
                cls: os.path.join(output_folder, organ, mode, cls)
                for cls in CLONALITY_CLASSES
            }
            for organ in self._organs
        }
        for folders in organ_folders.values():
            for folder in folders.values():
                os.makedirs(folder, exist_ok=True)

        written_folders: set[str] = set()
        vcf_files = [f for f in os.listdir(input_folder) if f.endswith(".vcf")]

        for filename in tqdm(vcf_files, desc="Processing VCFs", file=sys.stdout):
            aliquot_id = filename.split(".")[0]
            organ = self._organ_lookup.get(aliquot_id)
            if not organ:
                continue

            sample = vcf_io.vcf_to_dataframe(os.path.join(input_folder, filename))
            sample = vcf_io.process_vcf_dataframe(sample, time_analysis=time_analysis)

            for cls in CLONALITY_CLASSES:
                subset = sample[sample["CLS"] == cls.lower()]
                if subset.empty:
                    continue
                out_vcf = os.path.join(
                    organ_folders[organ][cls], f"{aliquot_id}_{cls.lower()}.vcf"
                )
                vcf_io.dataframe_to_vcf(subset, out_vcf)
                written_folders.add(os.path.dirname(out_vcf))

        self._cleanup_empty(output_folder, written_folders, verbose)

    def _cleanup_empty(
        self, output_folder: str, written_folders: set[str], verbose: bool
    ) -> None:
        for organ in self._organs:
            base_dir = os.path.join(output_folder, organ)
            keep = any(w.startswith(base_dir) for w in written_folders)
            if not keep and os.path.exists(base_dir):
                shutil.rmtree(base_dir)
                if verbose:
                    print(f"Removed empty folder: {base_dir}")

    @staticmethod
    def prepare_for_timing(vcf_folder_path: str) -> dict[str, pd.DataFrame] | None:
        """Merge per-mutation genomic info with signature posterior probabilities.

        Returns ``{sample_id: DataFrame}`` or ``None`` if nothing was processed.
        Expects a ``SignatureFitting`` directory produced by the signature step.
        """
        sample_dfs: dict[str, pd.DataFrame] = {}
        info_fields = ["context", "MajCN", "MinCN", "MutCN", "pSingle", "pGain", "VAF"]

        for subfolder in CLONALITY_CLASSES:
            subfolder_path = os.path.join(vcf_folder_path, subfolder)
            if not os.path.exists(subfolder_path):
                continue

            print(f"Processing {subfolder} samples:")
            vcf_files = [f for f in os.listdir(subfolder_path) if f.endswith(".vcf")]

            for filename in tqdm(vcf_files, desc=f"Merging {subfolder}", file=sys.stdout):
                sample_id = filename.split("_")[0]
                vcf_file_path = os.path.join(subfolder_path, filename)
                if os.path.getsize(vcf_file_path) == 0:
                    continue

                try:
                    df = pd.read_csv(vcf_file_path, sep="\t", header=None, comment="#",
                                     keep_default_na=False, na_values=[])
                except pd.errors.EmptyDataError:
                    continue

                df = df[[0, 1, 3, 4, 6]].copy()
                df.columns = ["CHROM", "POS", "REF", "ALT", "INFO"]
                df["INFO"] = df["INFO"].apply(vcf_io.info_field_to_dict)
                df["CHROM"] = df["CHROM"].astype(str)

                for field in info_fields:
                    col = "SBS96" if field == "context" else field
                    df[col] = df["INFO"].apply(lambda info: info.get(field))

                df["Classification"] = subfolder

                prob_path = os.path.join(
                    vcf_folder_path,
                    "SignatureFitting",
                    subfolder,
                    "probabilities",
                    f"prob_{sample_id}_{subfolder.lower()}.txt",
                )
                if os.path.exists(prob_path):
                    prob_df = pd.read_csv(prob_path, sep="\t", index_col=0)
                    prob_df.columns = [f"prob_{c}" for c in prob_df.columns]
                    df = df.join(prob_df, on="SBS96")

                if sample_id not in sample_dfs:
                    sample_dfs[sample_id] = df
                else:
                    sample_dfs[sample_id] = pd.concat(
                        [sample_dfs[sample_id], df], ignore_index=True
                    )

        return sample_dfs or None
