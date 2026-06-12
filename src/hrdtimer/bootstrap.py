"""Bootstrap resampling of mutation signatures.

Port of the vectorised ``generate_bootstraps`` (from ``timing.py``). Each
replicate resamples mutations per sample/class via a multinomial draw, refits
exposures, and writes per-sample annotated mutations plus cohort exposure
matrices. A ``seed`` argument is added for reproducibility (defaults to
non-deterministic, matching the original behaviour).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

from .catalog import SignatureCatalog
from .constants import DEFAULT_TUMOR_TYPE, TIMING_CLASSES

__all__ = ["Bootstrapper"]


class Bootstrapper:
    """Generate bootstrap replicates of signature posterior probabilities."""

    def __init__(
        self, tumor_type: str = DEFAULT_TUMOR_TYPE, seed: int | None = None
    ) -> None:
        self.catalog = SignatureCatalog(tumor_type=tumor_type)
        self._rng = np.random.default_rng(seed)

    def generate(
        self,
        samples: dict[str, pd.DataFrame],
        n_bootstraps: int,
        output_dir: str,
    ) -> None:
        """Write ``n_bootstraps`` replicates under ``output_dir/bootstrap_<i>/``."""
        catalog = self.catalog
        mutation_types = catalog.mutation_types
        os.makedirs(output_dir, exist_ok=True)

        for i in tqdm(range(1, n_bootstraps + 1), desc="Bootstrapping", file=sys.stdout):
            boot_dir = os.path.join(output_dir, f"bootstrap_{i}")
            os.makedirs(boot_dir, exist_ok=True)
            group_exposures: dict[str, list[pd.Series]] = {
                label: [] for label in TIMING_CLASSES
            }

            for sample_id, df in samples.items():
                merged = []
                for label in TIMING_CLASSES:
                    group = df[df["Classification"] == label].copy()
                    if group.empty:
                        continue

                    # Multinomial resample of mutation contexts.
                    counts = group["SBS96"].value_counts()
                    type_probs = counts.reindex(mutation_types, fill_value=0) / len(group)
                    sampled_counts = self._rng.multinomial(len(group), type_probs.values)
                    sampled_df = pd.DataFrame(
                        {sample_id: sampled_counts}, index=mutation_types
                    )
                    sampled_df.index.name = "Type"

                    exposures, _ = catalog.refit(sampled_df)
                    exp_series = exposures.iloc[:, 0]
                    exp_series.name = sample_id
                    group_exposures[label].append(exp_series)

                    total = exp_series.sum()
                    exposures_norm = exp_series / total if total > 0 else exp_series
                    prob_df = catalog.mutation_probabilities(exposures_norm).rename(
                        columns=lambda c: f"prob_{c}_boot"
                    )
                    merged.append(group.join(prob_df, on="SBS96"))

                if merged:
                    pd.concat(merged, ignore_index=True).to_csv(
                        os.path.join(boot_dir, f"{sample_id}.csv"), index=False
                    )

            for label, exposure_list in group_exposures.items():
                if exposure_list:
                    pd.concat(exposure_list, axis=1).T.to_csv(
                        os.path.join(boot_dir, f"exposures_{label}.csv")
                    )
