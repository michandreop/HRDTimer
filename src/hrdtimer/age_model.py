"""Monte Carlo age model: HRD crossing-age inference from SBS1 burden.

This replaces the old deterministic line-intersection age model
(``MolecularTime_to_Age``). For each tumour it runs a Monte Carlo over the
normal-tissue slope, subtype-specific tumour duration, and the HRD molecular
time, fitting **two** crossing-age models per iteration:

* **gradual**  - a logistic (sigmoid) transition of the SBS1 rate, and
* **two-step** - a piecewise-linear (normal -> elevated) rate.

The inference core here depends only on numpy/pandas/scipy -- the figures from
the original notebook live in a separate (optional) plotting module so the
pipeline does not pull in matplotlib.

Ported from the ``HRDTimerAnalysis`` notebook class; the per-iteration maths is
verified numerically identical to the original.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from scipy.stats import gaussian_kde
from scipy.stats import norm as sp_norm

from .constants import DEFAULT_AGE_CONFIG

__all__ = [
    "SlopeDistribution",
    "AgeModel",
    "run_monte_carlo",
    "genome_scaling_factor",
    "sigmoid_curve",
    "sigmoid_derivative",
    "two_step_curve",
]


# --- genome scaling (G) ------------------------------------------------------
def genome_scaling_factor(df: pd.DataFrame) -> float:
    """SBS1 genome-scaling factor ``G`` (vectorised).

    ``G = sum(prob_SBS1) / sum(prob_SBS1 * MutCN / (MajCN + MinCN))``. Kept as a
    utility for producing ``scaled_SBS1`` upstream of the age model.
    """
    numerator = df["prob_SBS1"].sum()
    denominator = (df["prob_SBS1"] * df["MutCN"] / (df["MajCN"] + df["MinCN"])).sum()
    return numerator / denominator


# --- normal-tissue slope sampler --------------------------------------------
class SlopeDistribution:
    """KDE-based sampler for the normal-tissue SBS1 accumulation slope ``s_n``."""

    def __init__(self, csv_path: str, bw_method: float = 0.5) -> None:
        df = pd.read_csv(csv_path, index_col=0, na_values=["None"])
        df = df[df["driver"].isna()]
        x, y = df["age"].values, df["scaled_SBS1"].values
        x, y = x[x > 0], y[x > 0]

        self._s_array = y / x
        self.mean = float(np.mean(self._s_array))
        self.low = float(np.percentile(self._s_array, 5))
        self.high = float(np.percentile(self._s_array, 95))

        s_grid = np.linspace(self._s_array.min(), self._s_array.max(), 1000)
        pdf = gaussian_kde(self._s_array, bw_method=bw_method)(s_grid)
        pdf[[0, -1]] = 0
        pdf /= np.trapz(pdf, s_grid)
        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]
        self._inv_cdf = interp1d(
            cdf, s_grid, bounds_error=False,
            fill_value=(self._s_array.min(), self._s_array.max()),
        )

    def sample(self, n: int, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Draw *n* slope samples. Pass ``rng`` for reproducibility."""
        u = rng.uniform(0, 1, n) if rng is not None else np.random.uniform(0, 1, n)
        return self._inv_cdf(u)


# --- sigmoid (gradual) model -------------------------------------------------
def _sigmoid_integral(v: np.ndarray, mid: float, k: float) -> np.ndarray:
    """Anti-derivative of the logistic: (1/k) * log(1 + exp(k*(v-mid)))."""
    return (1.0 / k) * np.log1p(np.exp(k * (v - mid)))


def sigmoid_curve(x, mid: float, k: float, s_n: float, rd: float):
    """Cumulative SBS1 burden under a sigmoid rate model."""
    return s_n * x + rd * (_sigmoid_integral(x, mid, k) - _sigmoid_integral(0, mid, k))


def sigmoid_derivative(x, mid: float, k: float, s_n: float, rd: float):
    """Instantaneous mutation rate dy/dx at age *x*."""
    return s_n + rd / (1.0 + np.exp(-k * (x - mid)))


def fit_midpoint(
    k: float, s_n: float, s_late: float, age: float, y_obs: float
) -> Tuple[Optional[float], Optional[float]]:
    """Solve for the sigmoid midpoint through (age, y_obs) matching ``s_late``."""
    rd = s_late - s_n
    if rd <= 0:
        return None, None

    def residuals(m):
        return [
            sigmoid_curve(age, m[0], k, s_n, rd) - y_obs,
            sigmoid_derivative(age * 0.95, m[0], k, s_n, rd) - s_late,
        ]

    try:
        result = least_squares(
            residuals, [age * 0.5], bounds=(0, age), method="trf", max_nfev=200
        )
        return (float(result.x[0]), rd) if result.cost < 1e-3 else (None, None)
    except Exception:
        return None, None


def find_burden_crossing(
    mid: float, k: float, s_n: float, rd: float, y_target: float, age: float
) -> float:
    """Age at which the cumulative sigmoid curve first equals ``y_target``."""
    try:
        result = least_squares(
            lambda a: sigmoid_curve(a[0], mid, k, s_n, rd) - y_target,
            [age * 0.5], bounds=(0, age), method="trf", max_nfev=200,
        )
        return float(result.x[0]) if result.cost < 1e-4 else np.nan
    except Exception:
        return np.nan


# --- two-step (piecewise-linear) model --------------------------------------
def two_step_curve(x, x_div: float, s_n: float, s_late: float):
    """Cumulative SBS1 burden: normal slope until ``x_div``, then elevated."""
    return np.where(x <= x_div, s_n * x, s_n * x_div + s_late * (x - x_div))


def find_two_step_crossing(
    x_div: float, s_n: float, s_late: float, y_target: float, age: float
) -> float:
    """Analytic HRD crossing age under the two-step model (nan if out of range)."""
    if s_late <= 0:
        return np.nan
    y_at_div = s_n * x_div
    if y_target < y_at_div:
        a = y_target / s_n if s_n > 0 else np.nan
    else:
        a = x_div + (y_target - y_at_div) / s_late
    return float(a) if (not np.isnan(a) and 0 <= a <= age) else np.nan


def sample_steepness(
    dur: float, dur_mean: float, dur_sd: float, s_late: float, s_n: float,
    k_max: float, rng: np.random.Generator, p_transition: float,
) -> Tuple[float, float, float]:
    """Draw a sigmoid steepness ``k`` consistent with the duration prior."""
    logit_p = np.log(p_transition / (1 - p_transition))
    k_min = 2 * logit_p / dur
    k_min_eff = min(k_min, k_max)
    k_range = k_max - k_min_eff

    p_combined = (1 - sp_norm.cdf(dur, loc=dur_mean, scale=dur_sd)) * np.clip(
        (s_late - s_n) / s_n, 0, 40
    ) / 40.0
    k_centre = k_min_eff + p_combined * k_range
    k_draw = float(rng.normal(k_centre, max(k_range * 0.1, 1e-6)))
    return float(max(k_draw, k_min_eff)), float(k_centre), float(k_min)


# --- single-sample Monte Carlo ----------------------------------------------
def run_monte_carlo(
    row: pd.Series,
    slope_dist: SlopeDistribution,
    k_max: float,
    subtype_durations: Dict,
    n_iter: int = 1000,
    max_attempts_factor: int = 20,
    rng: Optional[np.random.Generator] = None,
    p_transition: float = 0.90,
    late_sbs1_burden: float = 0.0191667,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int], int]:
    """Monte Carlo crossing-age inference for one tumour (gradual + two-step).

    ``row`` must contain: age, scaled_SBS1, type, HRDTime, HRDTime_ci_hi,
    HRDTime_ci_lo. Returns ``(df_converged, df_all, fail_counts, n_attempts)``.
    """
    if rng is None:
        rng = np.random.default_rng()

    age = float(row["age"])
    y_obs = float(row["scaled_SBS1"]) + late_sbs1_burden
    stype = row["type"]
    dur_mean = subtype_durations[stype]["duration_mean"]
    dur_sd = (
        subtype_durations[stype]["duration_CI"][1]
        - subtype_durations[stype]["duration_CI"][0]
    ) / (2 * 1.96)
    hrd_mu = float(row["HRDTime"])
    hrd_sd = (float(row["HRDTime_ci_hi"]) + float(row["HRDTime_ci_lo"])) / (2 * 1.96)

    fail_counts = {
        "dur_invalid": 0, "s_late_leq_sn": 0, "mid_fit_failed": 0, "crossing_nan": 0
    }
    converged_records: List[Dict] = []
    all_records: List[Dict] = []
    max_attempts = n_iter * max_attempts_factor
    attempt_idx = iter_idx = 0

    while len(converged_records) < n_iter and attempt_idx < max_attempts:
        s_n = float(slope_dist.sample(1)[0])
        dur = rng.normal(dur_mean, dur_sd)
        x_div = age - dur
        hrd_f = float(np.clip(rng.normal(hrd_mu, hrd_sd), 0.01, 0.99))
        y_tgt = hrd_f * y_obs
        y_tgt_l = np.clip(hrd_f - abs(rng.normal(0, hrd_sd * 0.5)), 0, 1) * y_obs
        y_tgt_h = np.clip(hrd_f + abs(rng.normal(0, hrd_sd * 0.5)), 0, 1) * y_obs

        record: Dict = {
            "attempt_idx": attempt_idx, "iter_idx": np.nan,
            "s_n": s_n, "dur": dur, "x_div": x_div,
            "hrd_frac": hrd_f, "y_tgt": y_tgt, "y_tgt_l": y_tgt_l, "y_tgt_h": y_tgt_h,
            "s_late": np.nan, "rel_excess": np.nan,
            "k_centre": np.nan, "k": np.nan, "k_min_iter": np.nan,
            "mid": np.nan, "rd": np.nan,
            "HRD_crossing_age": np.nan, "HRD_crossing_age_low": np.nan,
            "HRD_crossing_age_high": np.nan, "HRD_Duration": np.nan,
            "TS_crossing_age": np.nan, "TS_crossing_age_low": np.nan,
            "TS_crossing_age_high": np.nan, "TS_Duration": np.nan,
            "converged": False, "fail_reason": None,
        }

        if dur <= 0 or x_div <= 0:
            fail_counts["dur_invalid"] += 1
            record["fail_reason"] = "dur_invalid"
            all_records.append(record)
            attempt_idx += 1
            continue

        s_late = (y_obs - s_n * x_div) / dur
        record["s_late"] = s_late
        record["rel_excess"] = float(np.clip((s_late - s_n) / s_n, 0, 40))

        if s_late <= s_n:
            fail_counts["s_late_leq_sn"] += 1
            record["fail_reason"] = "s_late_leq_sn"
            all_records.append(record)
            attempt_idx += 1
            continue

        # Two-step crossing (computed before the sigmoid fit).
        ts_ac = find_two_step_crossing(x_div, s_n, s_late, y_tgt, age)
        record.update({
            "TS_crossing_age": ts_ac,
            "TS_crossing_age_low": find_two_step_crossing(x_div, s_n, s_late, y_tgt_l, age),
            "TS_crossing_age_high": find_two_step_crossing(x_div, s_n, s_late, y_tgt_h, age),
            "TS_Duration": age - ts_ac if not np.isnan(ts_ac) else np.nan,
        })

        # Sigmoid fit.
        k, k_centre, k_min_iter = sample_steepness(
            dur, dur_mean, dur_sd, s_late, s_n, k_max, rng, p_transition
        )
        record.update({"k": k, "k_centre": k_centre, "k_min_iter": k_min_iter})

        mid, rd = fit_midpoint(k, s_n, s_late, age, y_obs)
        record["mid"] = mid if mid is not None else np.nan
        record["rd"] = rd if rd is not None else np.nan

        if mid is None:
            fail_counts["mid_fit_failed"] += 1
            record["fail_reason"] = "mid_fit_failed"
            all_records.append(record)
            attempt_idx += 1
            continue

        ac = find_burden_crossing(mid, k, s_n, rd, y_tgt, age)
        record.update({
            "HRD_crossing_age": ac,
            "HRD_crossing_age_low": find_burden_crossing(mid, k, s_n, rd, y_tgt_l, age),
            "HRD_crossing_age_high": find_burden_crossing(mid, k, s_n, rd, y_tgt_h, age),
            "HRD_Duration": age - ac if not np.isnan(ac) else np.nan,
        })

        if np.isnan(ac):
            fail_counts["crossing_nan"] += 1
            record["fail_reason"] = "crossing_nan"
            all_records.append(record)
            attempt_idx += 1
            continue

        record.update({"converged": True, "fail_reason": "none", "iter_idx": iter_idx})
        iter_idx += 1
        converged_records.append(record)
        all_records.append(record)
        attempt_idx += 1

    return (
        pd.DataFrame(converged_records),
        pd.DataFrame(all_records),
        fail_counts,
        attempt_idx,
    )


def _percentile_summary(values: np.ndarray, prefix: str) -> Dict:
    if len(values) == 0:
        return {f"{prefix}_median": np.nan, f"{prefix}_ci05": np.nan, f"{prefix}_ci95": np.nan}
    return {
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_ci05": float(np.percentile(values, 5)),
        f"{prefix}_ci95": float(np.percentile(values, 95)),
    }


def build_sample_summary(
    cohort: str, row: pd.Series, df_conv: pd.DataFrame,
    fail_counts: Dict, n_attempts: int, k_max: float,
) -> Dict:
    """Assemble the one-row per-sample summary."""
    s: Dict = {
        "cohort": cohort,
        "sample": row["sample.display"],
        "type": row["type"],
        "age_at_diagnosis": row["age"],
        "n_converged": len(df_conv),
        "n_attempts": n_attempts,
        "convergence_rate": 100 * len(df_conv) / max(n_attempts, 1),
        "k_max_used": k_max,
        **{f"fail_{k}": v for k, v in fail_counts.items()},
    }
    gradual_cols = [
        "HRD_crossing_age", "HRD_Duration", "k", "k_min_iter",
        "s_n", "dur", "s_late", "k_centre", "rel_excess",
    ]
    for col in gradual_cols + ["TS_crossing_age", "TS_Duration"]:
        vals = df_conv[col].dropna().values if col in df_conv.columns else np.array([])
        s.update(_percentile_summary(vals, col))
    return s


# --- orchestration -----------------------------------------------------------
class AgeModel:
    """End-to-end Monte Carlo crossing-age analysis over a merged cohort table.

    Parameters
    ----------
    data
        Merged DataFrame with one row per sample: must contain ``sample.display``,
        ``cohort``, ``type``, ``age``, ``scaled_SBS1``, ``HRDTime``,
        ``HRDTime_ci_hi``, ``HRDTime_ci_lo``.
    normal_tissue_csv
        Path to the normal-tissue SBS1 reference (for the slope distribution).
    config
        Overrides for :data:`hrdtimer.constants.DEFAULT_AGE_CONFIG`.
    """

    def __init__(
        self, data: pd.DataFrame, normal_tissue_csv: str, config: Optional[Dict] = None
    ) -> None:
        self.cfg = {**DEFAULT_AGE_CONFIG, **(config or {})}
        p = self.cfg["p_transition"]
        self._logit_p = np.log(p / (1 - p))
        self.k_max = 2 * self._logit_p / self.cfg["min_transition_years"]

        self.data = data[data["type"].isin(self.cfg["subtype_durations"])].copy()
        self.slope_dist = SlopeDistribution(normal_tissue_csv)
        self.rng = np.random.default_rng(self.cfg["rng_seed"])

        sd = self.cfg["subtype_durations"]
        self.dur_mean_map = {s: d["duration_mean"] for s, d in sd.items()}
        self.dur_sd_map = {
            s: (d["duration_CI"][1] - d["duration_CI"][0]) / (2 * 1.96)
            for s, d in sd.items()
        }

        self.summary_df: Optional[pd.DataFrame] = None
        self.iter_df: Optional[pd.DataFrame] = None
        self.diag_df: Optional[pd.DataFrame] = None
        self._conv_map: Dict[str, pd.DataFrame] = {}
        self._diag_conv: Optional[pd.DataFrame] = None

    def run(self, verbose: bool = True) -> "AgeModel":
        """Run the Monte Carlo across all samples; populate result DataFrames."""
        summary_rows: List[Dict] = []
        iter_rows: List[pd.DataFrame] = []
        diag_rows: List[pd.DataFrame] = []
        global_fail = {
            "dur_invalid": 0, "s_late_leq_sn": 0, "mid_fit_failed": 0, "crossing_nan": 0
        }
        total_conv = total_att = 0

        for cohort, cdata in self.data.groupby("cohort"):
            for _, row in cdata.iterrows():
                sid = row["sample.display"]
                df_conv, df_all, fail_counts, n_att = run_monte_carlo(
                    row,
                    slope_dist=self.slope_dist,
                    k_max=self.k_max,
                    subtype_durations=self.cfg["subtype_durations"],
                    n_iter=self.cfg["n_iter"],
                    max_attempts_factor=self.cfg["max_attempts_factor"],
                    rng=self.rng,
                    p_transition=self.cfg["p_transition"],
                    late_sbs1_burden=self.cfg["late_sbs1_burden"],
                )

                for df in (df_conv, df_all):
                    df["cohort"] = cohort
                    df["sample"] = sid
                    df["type"] = row["type"]
                    df["age_at_dx"] = row["age"]

                if verbose:
                    self._report_sample(sid, row, fail_counts, len(df_conv), n_att)

                for key in global_fail:
                    global_fail[key] += fail_counts.get(key, 0)
                total_conv += len(df_conv)
                total_att += n_att

                iter_rows.append(df_conv)
                diag_rows.append(df_all)
                self._conv_map[sid] = df_conv
                summary_rows.append(
                    build_sample_summary(cohort, row, df_conv, fail_counts, n_att, self.k_max)
                )

        if verbose:
            self._report_global(global_fail, total_conv, total_att)

        self.summary_df = pd.DataFrame(summary_rows)
        self.iter_df = pd.concat(iter_rows, ignore_index=True)
        self.diag_df = pd.concat(diag_rows, ignore_index=True)

        diag_conv = self.diag_df[self.diag_df["converged"]].copy()
        diag_conv["p_dur"] = diag_conv.apply(
            lambda r: sp_norm.cdf(
                r["dur"], loc=self.dur_mean_map[r["type"]], scale=self.dur_sd_map[r["type"]]
            ),
            axis=1,
        )
        diag_conv["p_excess"] = diag_conv["rel_excess"].clip(0, 40) / 40.0
        diag_conv["p_combined"] = (1 - diag_conv["p_dur"]) * diag_conv["p_excess"]
        self._diag_conv = diag_conv
        return self

    def save(
        self,
        summary_path: str,
        iter_path: Optional[str] = None,
        diag_path: Optional[str] = None,
    ) -> None:
        """Write the summary (and optionally iterations + diagnostics) to CSV."""
        self._ensure_run()
        os.makedirs(os.path.dirname(summary_path) or ".", exist_ok=True)
        self.summary_df.to_csv(summary_path, index=False)
        if iter_path:
            self.iter_df.to_csv(iter_path, index=False)
        if diag_path:
            self.diag_df.to_csv(diag_path, index=False)

    def _ensure_run(self) -> None:
        if self.summary_df is None:
            raise RuntimeError("Call .run() before accessing results.")

    @staticmethod
    def _report_sample(sid, row, fail_counts, n_converged, n_attempts, stream=sys.stderr):
        pct = 100 * n_converged / max(n_attempts, 1)
        print(f"\nSAMPLE {sid:<25s}  type={row['type']}  age={row['age']:.1f}", file=stream)
        print(f"  converged: {n_converged:>4d} / {n_attempts:<4d}  ({pct:5.1f}%)", file=stream)
        for reason, cnt in fail_counts.items():
            if cnt:
                print(f"    {reason:<22s}  {cnt}", file=stream)

    @staticmethod
    def _report_global(global_fail, total_conv, total_att, stream=sys.stderr):
        pct = 100 * total_conv / max(total_att, 1)
        sep = "=" * 52
        print(f"\n{sep}\n  GLOBAL CONVERGENCE: {total_conv}/{total_att} ({pct:.1f}%)\n{sep}",
              file=stream)
        for reason, cnt in global_fail.items():
            print(f"    {reason:<22s}  {cnt}", file=stream)

    # -- plotting (delegates to age_model_plots; matplotlib imported lazily) --
    def plot_cohort(self, cohort: str, show: bool = True):
        """Burden + sigmoid grid figures for *cohort*."""
        self._ensure_run()
        from . import age_model_plots as plots

        cdata = self.data[self.data["cohort"] == cohort]
        conv_map = {sid: self._conv_map[sid]
                    for sid in cdata["sample.display"] if sid in self._conv_map}
        plots.set_mpl_style()
        fig_burden = plots.plot_cohort_burden(cohort, cdata, conv_map, self.slope_dist)
        fig_sigmoid = plots.plot_cohort_sigmoid(cohort, cdata, conv_map)
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig_burden, fig_sigmoid

    def plot_two_step_cohort(self, cohort: str, show: bool = True):
        """Two-step burden grid figure for *cohort*."""
        self._ensure_run()
        from . import age_model_plots as plots

        cdata = self.data[self.data["cohort"] == cohort]
        conv_map = {sid: self._conv_map[sid]
                    for sid in cdata["sample.display"] if sid in self._conv_map}
        plots.set_mpl_style()
        fig = plots.plot_cohort_two_step(cohort, cdata, conv_map, self.slope_dist)
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig

    def plot_sample(self, sample_display: str, show: bool = True):
        """4-panel gradual diagnostic figure for one sample."""
        self._ensure_run()
        from . import age_model_plots as plots

        row = self.data[self.data["sample.display"] == sample_display]
        if row.empty:
            raise ValueError(f"Sample '{sample_display}' not found.")
        df_conv = self._conv_map.get(sample_display, pd.DataFrame())
        plots.set_mpl_style()
        fig = plots.make_sample_figure(
            row.iloc[0], df_conv, self.slope_dist, self.dur_mean_map, self.dur_sd_map
        )
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig

    def plot_two_step_sample(self, sample_display: str, show: bool = True):
        """3-panel two-step diagnostic figure for one sample."""
        self._ensure_run()
        from . import age_model_plots as plots

        row = self.data[self.data["sample.display"] == sample_display]
        if row.empty:
            raise ValueError(f"Sample '{sample_display}' not found.")
        df_conv = self._conv_map.get(sample_display, pd.DataFrame())
        plots.set_mpl_style()
        fig = plots.make_two_step_sample_figure(row.iloc[0], df_conv, self.slope_dist)
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig

    def plot_p_combined(self, output_path: Optional[str] = None, show: bool = True):
        """Global p_combined histogram grid across converged samples."""
        self._ensure_run()
        from . import age_model_plots as plots

        plots.set_mpl_style()
        fig = plots.plot_p_combined_grid(self._diag_conv)
        if output_path:
            fig.savefig(output_path, bbox_inches="tight")
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig

    def compare_models(
        self, metric: str = "HRD_Duration", groupby: str = "type",
        output_path: Optional[str] = None, show: bool = True,
    ):
        """Boxplot comparing gradual vs two-step per-sample medians."""
        self._ensure_run()
        from . import age_model_plots as plots

        plots.set_mpl_style()
        fig = plots.plot_model_comparison(self.iter_df, metric=metric, groupby=groupby)
        if output_path:
            fig.savefig(output_path, bbox_inches="tight")
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig

    def plot_hrd_onset(
        self,
        sample_order: Optional[List[str]] = None,
        meta_df=None,
        annotations: Optional[List[str]] = None,
        annotation_colors=None,
        sep_lines=None,
        sep_labels=None,
        title: str = "HRD onset vs age at diagnosis",
        output_path: Optional[str] = None,
        show: bool = True,
        **kwargs,
    ):
        """Lollipop plot of HRD crossing-age (median + 90% CI) vs age at diagnosis.

        Passes through to :func:`hrdtimer.age_model_plots.plot_hrd_onset`.
        Supply *meta_df* (indexed by sample name) and *annotations* to add a
        colour-coded annotation bar beneath the main panel.
        """
        self._ensure_run()
        from . import age_model_plots as plots
        import matplotlib.pyplot as plt

        plots.set_mpl_style()
        fig = plots.plot_hrd_onset(
            self.summary_df,
            sample_order=sample_order,
            meta_df=meta_df,
            annotations=annotations,
            annotation_colors=annotation_colors,
            sep_lines=sep_lines,
            sep_labels=sep_labels,
            title=title,
            **kwargs,
        )
        if output_path:
            fig.savefig(output_path, bbox_inches="tight")
        if show:
            plt.show()
        return fig

    def save_cohort_pdfs(self, output_dir: str) -> None:
        """Write burden + sigmoid + two-step cohort grids per cohort."""
        self._ensure_run()
        from . import age_model_plots as plots
        import matplotlib.pyplot as plt

        os.makedirs(output_dir, exist_ok=True)
        for cohort, cdata in self.data.groupby("cohort"):
            conv_map = {sid: self._conv_map[sid]
                        for sid in cdata["sample.display"] if sid in self._conv_map}
            plots.set_mpl_style()
            safe = cohort.replace(" ", "_")
            plots.plot_cohort_burden(cohort, cdata, conv_map, self.slope_dist).savefig(
                f"{output_dir}/HRD_crossing_ages_{safe}.pdf")
            plots.plot_cohort_sigmoid(cohort, cdata, conv_map).savefig(
                f"{output_dir}/HRD_sigmoid_transitions_{safe}.pdf", bbox_inches="tight")
            plots.plot_cohort_two_step(cohort, cdata, conv_map, self.slope_dist).savefig(
                f"{output_dir}/HRD_two_step_{safe}.pdf", bbox_inches="tight")
            plt.close("all")

    def save_sample_pdfs(self, output_dir: str) -> None:
        """Write gradual + two-step diagnostic PDFs per sample."""
        self._ensure_run()
        from . import age_model_plots as plots
        import matplotlib.pyplot as plt

        os.makedirs(output_dir, exist_ok=True)
        plots.set_mpl_style()
        for _, row in self.data.iterrows():
            sid = row["sample.display"]
            df_conv = self._conv_map.get(sid, pd.DataFrame())
            safe = sid.replace("/", "_").replace(" ", "_")
            fig_g = plots.make_sample_figure(
                row, df_conv, self.slope_dist, self.dur_mean_map, self.dur_sd_map)
            fig_g.savefig(f"{output_dir}/{safe}_gradual.pdf", bbox_inches="tight")
            plt.close(fig_g)
            fig_ts = plots.make_two_step_sample_figure(row, df_conv, self.slope_dist)
            fig_ts.savefig(f"{output_dir}/{safe}_two_step.pdf", bbox_inches="tight")
            plt.close(fig_ts)

    @property
    def cohorts(self) -> List[str]:
        return sorted(self.data["cohort"].unique().tolist())

    @property
    def samples(self) -> List[str]:
        return sorted(self.data["sample.display"].unique().tolist())
