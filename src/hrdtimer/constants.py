"""Shared constants for the HRDTimer pipeline.

Centralising these values removes the magic strings/numbers that were copied
across the original ``utils.py`` (CpG context lists, the COSMIC catalog name,
the SBS1 reference slope, etc.) and makes them easy to audit in one place.
"""

from __future__ import annotations

# --- Mutational contexts -----------------------------------------------------
# C>T mutations in a CpG context: the molecular clock used for SBS1 timing.
CPG_CONTEXTS: tuple[str, ...] = ("A[C>T]G", "C[C>T]G", "T[C>T]G", "G[C>T]G")

# --- Clonality / timing classes ---------------------------------------------
# Classes produced by the preprocessing split. ``Subclonal`` is written to disk
# but only the first three carry interpretable molecular timing information.
CLONALITY_CLASSES: tuple[str, ...] = ("Early", "Late", "NA", "Subclonal")
TIMING_CLASSES: tuple[str, ...] = ("Early", "Late", "NA")

# Minor copy-number states evaluated during timing.
MIN_CN_STATES: tuple[int, ...] = (0, 1, 2)

# --- Signature catalog -------------------------------------------------------
COSMIC_CATALOG: str = "COSMIC_v3p2_SBS_WGS"
DEFAULT_TUMOR_TYPE: str = "Breast.AdenoCA"
REFIT_METHOD: str = "likelihood_bidirectional"
REFIT_THRESH: float = 0.001

# --- Age model (genome scaling) ---------------------------------------------
# Genome size scaling (3000 Mb) used to normalise SBS1 burden.
GENOME_SCALE_MB: int = 3000

# --- Monte Carlo age model ---------------------------------------------------
# Configuration for the sigmoid / two-step crossing-age inference (AgeModel).
# Subtype tumour-duration priors (years) with 95% CIs.
DEFAULT_AGE_CONFIG: dict = {
    # Monte Carlo
    "n_iter": 1000,
    "max_attempts_factor": 20,
    "rng_seed": 42,
    # Sigmoid transition constraint
    "p_transition": 0.90,
    "min_transition_years": 3.0,
    # Constant SBS1 burden added to the MRCA estimate (late accumulation,
    # estimated from scNanoSeq).
    "late_sbs1_burden": 0.0191667,
    # Output
    "output_dir": "hrdtime_plots",
    # Subtype tumour-duration priors (years)
    "subtype_durations": {
        "ER+": {"duration_mean": 16.33, "duration_CI": [13.55, 19.10]},
        "TN": {"duration_mean": 10.36, "duration_CI": [6.60, 14.32]},
    },
}
