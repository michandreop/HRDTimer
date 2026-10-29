"""Smoke tests: the package installs, imports, and its timing arithmetic is correct.

These run on a clean machine with only the PyPI dependencies declared in
pyproject.toml. They deliberately avoid MuSiCal, SigProfiler and the reference
genomes, so a failure here means the base install is broken rather than an
optional component being absent.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

import hrdtimer
from hrdtimer import TimingEstimator, constants, dict_to_vcf_info, info_field_to_dict

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_public_api_resolves():
    """Every name advertised in __all__ is importable."""
    missing = [n for n in hrdtimer.__all__ if not hasattr(hrdtimer, n)]
    assert not missing, f"declared in __all__ but absent: {missing}"


def test_version_matches_pyproject():
    text = (REPO_ROOT / "pyproject.toml").read_text()
    declared = re.search(r'^version = "([^"]+)"', text, re.M).group(1)
    assert hrdtimer.__version__ == declared


def test_constants_are_sane():
    assert constants.TIMING_CLASSES == ("Early", "Late", "NA")
    assert set(constants.MIN_CN_STATES) == {0, 1, 2}
    # The CpG clock channels used for WGD timing.
    assert all(re.fullmatch(r"[ACGT]\[C>T\]G", c) for c in constants.CPG_CONTEXTS)


def test_vcf_info_roundtrip():
    """INFO fields survive encode/decode, with numbers coerced back to numbers.

    Booleans are not round-tripped: they decode as the strings "True"/"False",
    which is why callers compare against the MutationTimeR spelling directly.
    """
    original = {"MajCN": 2, "MinCN": 0, "pGain": 0.97, "CLS": "clonal [early]"}
    assert info_field_to_dict(dict_to_vcf_info(original)) == original


def _cpg_frame(p_gain: float, n: int = 50) -> pd.DataFrame:
    """n CpG mutations in a (2,0) segment, each with the same multiplicity posterior."""
    return pd.DataFrame(
        {
            "SBS96": ["A[C>T]G"] * n,
            "MajCN": [2] * n,
            "MinCN": [0] * n,
            "pSingle": [1.0 - p_gain] * n,
            "pGain": [p_gain] * n,
            "prob_SBS1_boot": [1.0] * n,
        }
    )


@pytest.mark.parametrize(
    "p_gain, expected",
    [
        (1.0, 1.0),        # every mutation pre-duplication -> WGD at the end
        (0.0, 0.0),        # every mutation post-duplication -> WGD at the start
        (0.5, 2 / 3),      # t = 2*pi2 / (pi1 + 2*pi2)
    ],
)
def test_wgd_timing_matches_closed_form(p_gain, expected):
    """WGD molecular time follows Eq. 5 exactly for a uniform multiplicity posterior."""
    t, _ = TimingEstimator.wgd_time_from_replicate(_cpg_frame(p_gain))
    assert t == pytest.approx(expected, abs=1e-12)


def test_timing_is_an_expectation_not_a_hard_call():
    """Mutations of uncertain multiplicity are weighted, not committed to a phase.

    Half the mutations certainly pre-duplication and half certainly post gives
    the same answer as every mutation being 50/50, because molecular time is an
    expectation over the multiplicity posterior (Methods, Eq. 1).
    """
    mixed = pd.concat([_cpg_frame(1.0, 25), _cpg_frame(0.0, 25)], ignore_index=True)
    uncertain = _cpg_frame(0.5, 50)
    t_mixed, _ = TimingEstimator.wgd_time_from_replicate(mixed)
    t_uncertain, _ = TimingEstimator.wgd_time_from_replicate(uncertain)
    assert t_mixed == pytest.approx(t_uncertain, abs=1e-12)


def test_private_aliases_still_resolve():
    """Scripts written against the pre-rename private names keep working."""
    assert TimingEstimator._wgd_single_bootstrap is TimingEstimator.wgd_time_from_replicate
    assert TimingEstimator._hrd_single_bootstrap is TimingEstimator.hrd_time_from_replicate
