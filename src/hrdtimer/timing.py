"""WGD and HRD molecular timing from bootstrap replicates.

This consolidates the timing logic from ``timing.py``:

* WGD timing  -> ``calculate_wgd_time_bootstrapping``
* HRD per-replicate -> the **second** ``calculate_hrd_time_single_boot`` (the
  one using ``np.nan`` placeholders and NaN-safe ``c_avg`` weighting)
* HRD aggregation -> ``calculate_hrd_time_single_sample``
* cohort driver -> ``run_hrd_wgd_timing_analysis`` (the second, DataFrame-returning
  copy from ``utils.py``)

Efficiency change: the original read every ``bootstrap_<i>/<sample>.csv`` file
**twice** (once for WGD, once for HRD). :class:`TimingEstimator` reads each file
**once** and computes both, halving timing-stage I/O. The per-replicate maths
are byte-for-byte the originals (verified numerically against them).

This module depends only on numpy/pandas/tqdm -- no signature libraries.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

from .constants import CPG_CONTEXTS, MIN_CN_STATES

__all__ = ["TimingEstimator"]


class TimingEstimator:
    """Estimate WGD/HRD timing for a cohort from pre-computed bootstrap files."""

    def __init__(self, bootstraps_dir: str, num_bootstrap: int = 200) -> None:
        self.bootstraps_dir = bootstraps_dir
        self.num_bootstrap = num_bootstrap

    # -- public API ----------------------------------------------------------
    def run(self, sample_ids, output_csv_path: str | None = None) -> pd.DataFrame:
        """Time every sample in ``sample_ids`` and return a results DataFrame."""
        rows = [
            self._estimate_sample(sample_id)
            for sample_id in tqdm(sample_ids, desc="Analyzing Timing", file=sys.stdout)
        ]
        df = pd.DataFrame(rows)
        if output_csv_path:
            os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
            df.to_csv(output_csv_path, index=False)
        return df

    # -- per-sample driver (single file pass) --------------------------------
    def _estimate_sample(self, sample_id: str) -> dict:
        wgd_means: list[float] = []
        hrd_means: list[float] = []
        n_mut_cpg_total = np.zeros(3)
        n_mut_all = None
        c_avg_vals: list[float] = []
        stats = {
            key: {mc: [] for mc in MIN_CN_STATES}
            for key in (
                "pi_2_SBS1", "pi_2_SBS3", "pi_1_SBS1", "pi_1_SBS3",
                "Nt_SBS1", "Nt_SBS3", "c_val",
            )
        }

        for i in range(1, self.num_bootstrap + 1):
            path = os.path.join(self.bootstraps_dir, f"bootstrap_{i}", f"{sample_id}.csv")
            if not os.path.exists(path):
                continue
            df = pd.read_csv(path)

            # --- WGD (single read) ---
            wgd_mean, n_cpg = self._wgd_single_bootstrap(df)
            n_mut_cpg_total += n_cpg
            if not np.isnan(wgd_mean):
                wgd_means.append(wgd_mean)

            # --- HRD (same read) ---
            (t_values, n_mut, cavg, c_list, pi2_1, pi2_3,
             nt1, nt3, pi1_1, pi1_3) = self._hrd_single_bootstrap(df)
            t_arr = np.asarray(t_values, dtype=float)
            n_arr = np.asarray(n_mut, dtype=float)
            mask = ~np.isnan(t_arr)
            if mask.any():
                hrd_means.append(np.average(t_arr[mask], weights=n_arr[mask]))
            if n_mut_all is None:
                n_mut_all = n_arr
            c_avg_vals.append(cavg)
            for mc in MIN_CN_STATES:
                stats["pi_2_SBS1"][mc].append(pi2_1[mc])
                stats["pi_2_SBS3"][mc].append(pi2_3[mc])
                stats["pi_1_SBS1"][mc].append(pi1_1[mc])
                stats["pi_1_SBS3"][mc].append(pi1_3[mc])
                stats["Nt_SBS1"][mc].append(nt1[mc])
                stats["Nt_SBS3"][mc].append(nt3[mc])
                stats["c_val"][mc].append(c_list[mc])

        return self._aggregate(
            sample_id, wgd_means, hrd_means, n_mut_cpg_total, n_mut_all,
            c_avg_vals, stats,
        )

    # -- per-replicate WGD ---------------------------------------------------
    @staticmethod
    def _wgd_single_bootstrap(df: pd.DataFrame) -> tuple[float, np.ndarray]:
        """One replicate's CpG-clock WGD time; also returns per-state CpG counts."""
        df = df[df["SBS96"].isin(CPG_CONTEXTS)]
        t_values = []
        n_mut = []
        for min_cn in MIN_CN_STATES:
            sub = df[df["MinCN"] == min_cn]
            n_mut.append(len(sub))
            if len(sub) == 0:
                t_values.append(np.nan)
                continue
            probs = sub["prob_SBS1_boot"]
            sum_probs = probs.sum()
            if sum_probs == 0:
                t_values.append(np.nan)
                continue
            pi_1 = (probs * sub["pSingle"]).sum() / sum_probs
            pi_2 = (probs * sub["pGain"]).sum() / sum_probs
            denom = pi_1 + 2 * pi_2
            if denom == 0:
                t = np.nan
            elif pi_2 == 0:
                t = 0.0
            else:
                t = ((3 if min_cn == 1 else 2) * pi_2) / denom
            t_values.append(t)

        t_arr = np.asarray(t_values, dtype=float)
        n_arr = np.asarray(n_mut, dtype=float)
        mask = ~np.isnan(t_arr)
        weighted = (
            np.average(t_arr[mask], weights=n_arr[mask]) if mask.any() else np.nan
        )
        return weighted, n_arr

    # -- per-replicate HRD ---------------------------------------------------
    @staticmethod
    def _hrd_single_bootstrap(df: pd.DataFrame):
        """One replicate's HRD timing across MinCN states (SBS3-corrected SBS1).

        Faithful port of the second ``calculate_hrd_time_single_boot``.
        """
        if "prob_SBS3_boot" not in df.columns:
            df = df.assign(prob_SBS3_boot=0)

        res = {
            mc: {
                "t": np.nan, "n": 0, "c": np.nan,
                "pi_2_SBS1": 0, "pi_2_SBS3": 0, "Nt_SBS1": 0, "Nt_SBS3": 0,
                "pi_1_SBS1": 0, "pi_1_SBS3": 0,
            }
            for mc in MIN_CN_STATES
        }
        c_avg = 0

        # Process 0 and 2 before 1: state 1 needs the weighted c_avg.
        for min_cn in (0, 2, 1):
            sub = df[df["MinCN"] == min_cn]
            res[min_cn]["n"] = len(sub)
            if len(sub) == 0:
                continue

            def get_pi(prob_col, p_col):
                weights = sub[prob_col]
                total = weights.sum()
                return (weights * sub[p_col]).sum() / total if total > 0 else np.nan

            pi_1_SBS1 = float(get_pi("prob_SBS1_boot", "pSingle"))
            pi_2_SBS1 = float(get_pi("prob_SBS1_boot", "pGain"))
            pi_1_SBS3 = float(get_pi("prob_SBS3_boot", "pSingle"))
            pi_2_SBS3 = float(get_pi("prob_SBS3_boot", "pGain"))
            Nt_SBS1 = float(sub["prob_SBS1_boot"].sum())
            Nt_SBS3 = float(sub["prob_SBS3_boot"].sum())

            if min_cn in (0, 2):
                pi_2_SBS1_prime = (
                    pi_2_SBS1 - (pi_2_SBS3 / pi_1_SBS3 * pi_1_SBS1)
                    if pi_1_SBS3 > 0 else np.nan
                )
                if pi_1_SBS3 > 0 and Nt_SBS3 > 0:
                    res[min_cn]["c"] = (pi_1_SBS1 * Nt_SBS1) / (pi_1_SBS3 * Nt_SBS3)
                if min_cn == 2:
                    num = sum(
                        res[k]["c"] * res[k]["n"]
                        for k in (0, 2) if not np.isnan(res[k]["c"])
                    )
                    den = sum(
                        res[k]["n"] for k in (0, 2) if not np.isnan(res[k]["c"])
                    )
                    c_avg = num / den if den > 0 else 0
            else:  # min_cn == 1
                pi_2_SBS1_prime = (
                    pi_2_SBS1 - (pi_2_SBS3 * c_avg * (Nt_SBS3 / Nt_SBS1))
                    if (Nt_SBS1 > 0 and pi_1_SBS3 > 0) else np.nan
                )

            denom = pi_1_SBS1 + 2 * pi_2_SBS1
            if not np.isnan(pi_2_SBS1_prime) and denom != 0:
                mult = 3 if min_cn == 1 else 2
                res[min_cn]["t"] = (mult * pi_2_SBS1_prime) / denom

            res[min_cn].update(
                pi_2_SBS1=pi_2_SBS1, pi_2_SBS3=pi_2_SBS3,
                Nt_SBS1=Nt_SBS1, Nt_SBS3=Nt_SBS3,
                pi_1_SBS1=pi_1_SBS1, pi_1_SBS3=pi_1_SBS3,
            )

        order = list(MIN_CN_STATES)
        return (
            [res[k]["t"] for k in order],
            [res[k]["n"] for k in order],
            c_avg,
            [res[k]["c"] for k in order],
            [res[k]["pi_2_SBS1"] for k in order],
            [res[k]["pi_2_SBS3"] for k in order],
            [res[k]["Nt_SBS1"] for k in order],
            [res[k]["Nt_SBS3"] for k in order],
            [res[k]["pi_1_SBS1"] for k in order],
            [res[k]["pi_1_SBS3"] for k in order],
        )

    # -- aggregation ---------------------------------------------------------
    @staticmethod
    def _ci95_halfwidth(data) -> float:
        d = np.asarray(data, dtype=float)
        d = d[~np.isnan(d)]
        if len(d) == 0:
            return np.nan
        lo, hi = np.percentile(d, [2.5, 97.5])
        return (hi - lo) / 2

    def _aggregate(
        self, sample_id, wgd_means, hrd_means, n_mut_cpg_total, n_mut_all,
        c_avg_vals, stats,
    ) -> dict:
        # --- WGD aggregate ---
        wgd_means = np.asarray(wgd_means, dtype=float)
        if len(wgd_means) == 0:
            wgd_time = np.nan
            w_ci_hi = w_ci_lo = w_iqr_hi = w_iqr_lo = 0
        else:
            wgd_time = np.nanmean(wgd_means)
            p = np.nanpercentile(wgd_means, [2.5, 25, 75, 97.5])
            w_ci_hi, w_ci_lo = p[3] - wgd_time, wgd_time - p[0]
            w_iqr_hi, w_iqr_lo = p[2] - wgd_time, wgd_time - p[1]

        # --- HRD aggregate ---
        hrd_means = np.asarray(hrd_means, dtype=float)
        if len(hrd_means) == 0:
            nan3 = [np.nan, np.nan, np.nan]
            return {
                "ID": sample_id,
                "HRDTime": np.nan, "HRDTime_ci_hi": np.nan, "HRDTime_ci_lo": np.nan,
                "HRDTime_ci_IQR_hi": np.nan, "HRDTime_ci_IQR_lo": np.nan,
                "WGDTime": wgd_time, "WGDTime_ci_hi": w_ci_hi, "WGDTime_ci_lo": w_ci_lo,
                "WGDTime_ci_IQR_hi": w_iqr_hi, "WGDTime_ci_IQR_lo": w_iqr_lo,
                "pi2SBS1": nan3, "pi2SBS1_ci": nan3, "pi2SBS3": nan3, "pi2SBS3_ci": nan3,
                "pi1SBS1": nan3, "pi1SBS1_ci": nan3, "pi1SBS3": nan3, "pi1SBS3_ci": nan3,
                "c": nan3, "c_avg": np.nan, "Nt_SBS1": nan3, "Nt_SBS3": nan3,
                "N_mut_CpG": n_mut_cpg_total.tolist(),
                "N_mut_all": n_mut_all.tolist() if n_mut_all is not None else [],
            }

        hrd_time = np.nanmean(hrd_means)
        p_ci = np.nanpercentile(hrd_means, [2.5, 97.5])
        p_iqr = np.nanpercentile(hrd_means, [25, 75])

        def mean3(key):
            return [np.nanmean(stats[key][mc]) for mc in MIN_CN_STATES]

        def err3(key):
            return [self._ci95_halfwidth(stats[key][mc]) for mc in MIN_CN_STATES]

        c_val_mean = [
            0 if mc == 1 else np.nanmean(stats["c_val"][mc]) for mc in MIN_CN_STATES
        ]

        return {
            "ID": sample_id,
            "HRDTime": hrd_time,
            "HRDTime_ci_hi": p_ci[1] - hrd_time,
            "HRDTime_ci_lo": hrd_time - p_ci[0],
            "HRDTime_ci_IQR_hi": p_iqr[1] - hrd_time,
            "HRDTime_ci_IQR_lo": hrd_time - p_iqr[0],
            "WGDTime": wgd_time,
            "WGDTime_ci_hi": w_ci_hi,
            "WGDTime_ci_lo": w_ci_lo,
            "WGDTime_ci_IQR_hi": w_iqr_hi,
            "WGDTime_ci_IQR_lo": w_iqr_lo,
            "pi2SBS1": mean3("pi_2_SBS1"), "pi2SBS1_ci": err3("pi_2_SBS1"),
            "pi2SBS3": mean3("pi_2_SBS3"), "pi2SBS3_ci": err3("pi_2_SBS3"),
            "pi1SBS1": mean3("pi_1_SBS1"), "pi1SBS1_ci": err3("pi_1_SBS1"),
            "pi1SBS3": mean3("pi_1_SBS3"), "pi1SBS3_ci": err3("pi_1_SBS3"),
            "c": c_val_mean,
            "c_avg": np.nanmean(c_avg_vals),
            "Nt_SBS1": mean3("Nt_SBS1"),
            "Nt_SBS3": mean3("Nt_SBS3"),
            "N_mut_CpG": n_mut_cpg_total.tolist(),
            "N_mut_all": n_mut_all.tolist() if n_mut_all is not None else [],
        }
