"""HRDTimer: timing of homologous recombination deficiency (HRD) and
whole-genome doubling (WGD) in tumours.

The pipeline is exposed as a small set of stage objects:

    >>> from hrdtimer import (
    ...     MutationPreprocessor, SignatureAnalyzer, Bootstrapper,
    ...     TimingEstimator, AgeModel,
    ... )

Stages
------
1. ``MutationPreprocessor``  - split MutationTimeR VCFs by organ/clonality and
   assemble per-sample timing inputs.
2. ``SignatureAnalyzer``     - COSMIC signature assignment + MuSiCal refit.
3. ``Bootstrapper``          - bootstrap resampling of signature probabilities.
4. ``TimingEstimator``       - WGD/HRD molecular timing from the bootstraps.
5. ``AgeModel``              - Monte Carlo crossing-age inference (gradual +
   two-step) from molecular timing and the SBS1 clock.
"""

from __future__ import annotations

from . import constants
from .age_model import AgeModel, SlopeDistribution, genome_scaling_factor, run_monte_carlo
from .bootstrap import Bootstrapper
from .catalog import SignatureCatalog
from .preprocessing import MutationPreprocessor
from .signatures import SignatureAnalyzer
from .timing import TimingEstimator
from .vcf import (
    dataframe_to_vcf,
    dict_to_vcf_info,
    extract_clonality,
    info_field_to_dict,
    process_vcf_dataframe,
    vcf_to_dataframe,
)

__version__ = "0.1.0"

__all__ = [
    "MutationPreprocessor",
    "SignatureAnalyzer",
    "SignatureCatalog",
    "Bootstrapper",
    "TimingEstimator",
    "AgeModel",
    "SlopeDistribution",
    "run_monte_carlo",
    "genome_scaling_factor",
    "process_vcf_dataframe",
    "vcf_to_dataframe",
    "dataframe_to_vcf",
    "extract_clonality",
    "info_field_to_dict",
    "dict_to_vcf_info",
    "constants",
]
