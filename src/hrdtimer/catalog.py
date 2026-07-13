"""COSMIC signature catalog wrapper.

In the original code the COSMIC catalog was loaded and restricted *inside* every
function that needed it (``calculate_exposures``, ``generate_bootstraps``,
``add_signature_probabilities``, ``get_signature_exposures`` ...). This class
loads/restricts it **once** and exposes the reused operations (refit, posterior
mutation probabilities, exposure aggregation) as methods, so the expensive
catalog setup happens a single time per analysis.
"""

from __future__ import annotations

import pandas as pd

from .constants import (
    COSMIC_CATALOG,
    DEFAULT_TUMOR_TYPE,
    REFIT_METHOD,
    REFIT_THRESH,
    TIMING_CLASSES,
)

__all__ = ["SignatureCatalog"]


class SignatureCatalog:
    """A restricted COSMIC SBS catalog plus the refit operations built on it."""

    def __init__(
        self,
        tumor_type: str = DEFAULT_TUMOR_TYPE,
        catalog_name: str = COSMIC_CATALOG,
        is_mmrd: bool = False,
        is_ppd: bool = False,
    ) -> None:
        from ._optional import require  # lazy: only needed when signatures are fit

        musical = require("musical")
        catalog = musical.load_catalog(catalog_name)
        catalog.restrict_catalog(tumor_type=tumor_type, is_MMRD=is_mmrd, is_PPD=is_ppd)
        self.tumor_type = tumor_type
        self.W: pd.DataFrame = catalog.W
        self.mutation_types: list[str] = self.W.index.tolist()

    # -- core refit ----------------------------------------------------------
    def refit(
        self,
        counts: pd.DataFrame,
        method: str = REFIT_METHOD,
        thresh: float = REFIT_THRESH,
    ):
        """Refit signature exposures for ``counts`` (mutation types x samples).

        Returns the ``(H, model)`` tuple from :func:`musical.refit.refit`.
        """
        from ._optional import require

        musical = require("musical")
        W = self.W.reindex(counts.index)
        return musical.refit.refit(counts, W, method=method, thresh=thresh)

    # -- derived quantities --------------------------------------------------
    def mutation_probabilities(self, exposures_norm: pd.Series) -> pd.DataFrame:
        """Posterior per-context signature probabilities (each row sums to 1)."""
        weighted = self.W.mul(exposures_norm, axis=1)
        return weighted.div(weighted.sum(axis=1), axis=0).fillna(0)

    def context_counts(self, group: pd.DataFrame) -> pd.Series:
        """Count SBS96 contexts in ``group``, reindexed to the full catalog."""
        return (
            group["SBS96"].value_counts().reindex(self.mutation_types, fill_value=0)
        )

    def aggregate_exposures(
        self, samples: dict[str, pd.DataFrame], labels=TIMING_CLASSES
    ) -> pd.DataFrame:
        """Cohort exposure matrix: per-sample exposures summed over ``labels``.

        Port of ``get_signature_exposures``. Returns a (signatures x samples)
        matrix with all-zero signature rows dropped.
        """
        summed: dict[str, pd.Series] = {}
        for sample_id, df in samples.items():
            exposures_sum = None
            for label in labels:
                group = df[df["Classification"] == label]
                if group.empty:
                    continue
                counts = self.context_counts(group)
                count_df = pd.DataFrame({sample_id: counts})
                count_df.index.name = "Type"
                exposures, _ = self.refit(count_df)
                raw = exposures.iloc[:, 0]
                exposures_sum = raw if exposures_sum is None else exposures_sum.add(
                    raw, fill_value=0
                )
            if exposures_sum is None:
                exposures_sum = pd.Series(
                    0, index=self.mutation_types, name=sample_id
                )
            summed[sample_id] = exposures_sum

        matrix = pd.DataFrame(summed)
        return matrix.loc[(matrix != 0).any(axis=1)]
