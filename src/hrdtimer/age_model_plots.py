"""Figures for the Monte Carlo age model.

Optional plotting layer for :mod:`hrdtimer.age_model`. Kept in a separate module
so the inference/pipeline path never imports matplotlib. Ported verbatim (logic)
from the ``000_finalFig5_code`` notebook; only the hard-coded late-SBS1 burden
and the imports were tidied.

The :class:`hrdtimer.age_model.AgeModel` plotting methods delegate here.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import norm as sp_norm

from .age_model import SlopeDistribution, sigmoid_curve, sigmoid_derivative, two_step_curve
from .constants import DEFAULT_AGE_CONFIG

_LATE = DEFAULT_AGE_CONFIG["late_sbs1_burden"]

# ── Palette / style ─────────────────────────────────────────────────────────
PALETTE = {
    "initiation": "#27AE60",
    "midpoint": "#E67E22",
    "diagnosis": "#2C3E50",
    "median_curve": "firebrick",
    "tumour_window": "#2ECC71",
    "curve_sample": "#212C3E",
    "normal_band": "lightgrey",
    "hrd_band": "lightsalmon",
    "er_plus": "#185FA5",
    "tn": "#A32D2D",
    "two_step": "#8E44AD",
}
STYPE_COLORS = {"ER+": PALETTE["er_plus"], "TN": PALETTE["tn"]}


def set_mpl_style() -> None:
    """Apply the shared matplotlib style used across all age-model figures."""
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "font.family": "Sans",
        "font.size": 8,
        "text.color": "black",
        "axes.labelcolor": "black",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.5,
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })


# ── Shared panel drawing ────────────────────────────────────────────────────
def draw_sigmoid_panel(ax, df_conv: pd.DataFrame, age: float, n_show: int = 200):
    """Draw sampled + median sigmoid rate transitions on *ax*."""
    if df_conv.empty:
        return None

    x_sig = np.linspace(0, age * 1.08, 400)
    idx = np.random.choice(len(df_conv), min(n_show, len(df_conv)), replace=False)
    for ii in idx:
        r = df_conv.iloc[ii]
        ax.plot(x_sig, sigmoid_derivative(x_sig, r["mid"], r["k"], r["s_n"], r["rd"]),
                color=PALETTE["curve_sample"], lw=0.2, alpha=0.10)

    med_idx = (df_conv["k"] - df_conv["k"].median()).abs().idxmin()
    med = df_conv.loc[med_idx]
    x_d_med = age - med["dur"]
    m_med, k_med = med["mid"], med["k"]
    s_n_med = med["s_n"]
    s_l_med = s_n_med + med["rd"]

    ax.axvspan(x_d_med, age, color=PALETTE["tumour_window"], alpha=0.10, zorder=0)
    ax.plot(x_sig, sigmoid_derivative(x_sig, m_med, k_med, s_n_med, med["rd"]),
            color=PALETTE["median_curve"], lw=1.0, alpha=0.9, zorder=5)
    ax.axvline(x_d_med, color=PALETTE["initiation"], lw=0.9, ls="--", alpha=0.9)
    ax.axvline(m_med, color=PALETTE["midpoint"], lw=0.9, ls="-", alpha=0.9)
    ax.axvline(age, color=PALETTE["diagnosis"], lw=0.9, ls=":", alpha=0.9)
    ax.axhline(s_n_med, color="steelblue", lw=0.5, ls=":", alpha=0.7)
    ax.axhline(s_l_med, color=PALETTE["median_curve"], lw=0.5, ls=":", alpha=0.7)
    return med, x_d_med, m_med, k_med, s_n_med, s_l_med


def draw_two_step_burden_panel(
    ax, df_conv: pd.DataFrame, age: float, y_obs: float, hrd_f: float,
    hrd_sd: float, slope_dist: SlopeDistribution, n_show: int = 50,
) -> None:
    """Draw two-step (piecewise-linear) cumulative burden curves on *ax*."""
    y_tgt = hrd_f * y_obs
    y_tgt_l = np.clip(hrd_f - hrd_sd, 0, 1) * y_obs
    y_tgt_h = np.clip(hrd_f + hrd_sd, 0, 1) * y_obs
    x_full = np.linspace(0, 90, 400)
    x_samp = np.linspace(0, age, 300)

    ax.fill_between(x_full, slope_dist.low * x_full, slope_dist.high * x_full,
                    color=PALETTE["normal_band"], alpha=0.15)
    ax.plot(x_full, slope_dist.mean * x_full, ":", color="black", lw=0.6, alpha=0.6)
    ax.fill_between([0, 90], [y_tgt_l] * 2, [y_tgt_h] * 2, color="orange", alpha=0.15)
    ax.plot([0, 90], [y_tgt] * 2, color="darkorange", ls="--", lw=0.5, alpha=0.7)

    if not df_conv.empty:
        idx = np.random.choice(len(df_conv), min(n_show, len(df_conv)), replace=False)
        for ii in idx:
            r = df_conv.iloc[ii]
            ax.plot(x_samp, two_step_curve(x_samp, r["x_div"], r["s_n"], r["s_late"]),
                    color=PALETTE["two_step"], lw=0.3, alpha=0.25)
        ts_vals = df_conv["TS_crossing_age"].dropna().values
        if len(ts_vals):
            ax.axvspan(np.percentile(ts_vals, 5), np.percentile(ts_vals, 95),
                       color=PALETTE["two_step"], alpha=0.12)
            ax.axvline(np.median(ts_vals), color=PALETTE["two_step"], lw=0.9,
                       ls="--", label=f"TS median={np.median(ts_vals):.1f}yr")

    ax.scatter([age], [y_obs], color="k", s=20, zorder=10)
    ax.set(xlim=(0, 90), ylim=(0, 0.15), xlabel="Age [yr]",
           ylabel="SBS1 / Gx3000 [Mb]", title="Two-step cumulative burden")
    ax.tick_params(labelsize=7)


def _sigmoid_legend_handles() -> List[Line2D]:
    return [
        Line2D([0], [0], color=PALETTE["initiation"], lw=1, ls="--",
               label="Tumour initiation (A-d)"),
        Line2D([0], [0], color=PALETTE["midpoint"], lw=1, ls="-",
               label="Midpoint m (fitted)"),
        Line2D([0], [0], color=PALETTE["diagnosis"], lw=1, ls=":", label="Diagnosis A"),
        Line2D([0], [0], color=PALETTE["median_curve"], lw=1.5, label="Median rate curve"),
        plt.Rectangle((0, 0), 1, 1, color=PALETTE["tumour_window"], alpha=0.15,
                      label="Tumour window [A-d, A]"),
    ]


# ── Per-sample diagnostic panels ────────────────────────────────────────────
def _draw_burden_panel(ax, df_conv, age, y_obs, hrd_f, hrd_sd, slope_dist):
    y_tgt = hrd_f * y_obs
    y_tgt_l = np.clip(hrd_f - hrd_sd, 0, 1) * y_obs
    y_tgt_h = np.clip(hrd_f + hrd_sd, 0, 1) * y_obs
    x_full = np.linspace(0, 90, 400)
    x_samp = np.linspace(0, age, 200)

    ax.fill_between(x_full, slope_dist.low * x_full, slope_dist.high * x_full,
                    color=PALETTE["normal_band"], alpha=0.15)
    ax.plot(x_full, slope_dist.mean * x_full, ":", color="black", lw=0.6, alpha=0.6)
    ax.fill_between([0, 90], [y_tgt_l] * 2, [y_tgt_h] * 2, color="orange", alpha=0.15)
    ax.plot([0, 90], [y_tgt] * 2, color="darkorange", ls="--", lw=0.5, alpha=0.7)

    if not df_conv.empty:
        idx = np.random.choice(len(df_conv), min(50, len(df_conv)), replace=False)
        for ii in idx:
            r = df_conv.iloc[ii]
            ax.plot(x_samp, sigmoid_curve(x_samp, r["mid"], r["k"], r["s_n"], r["rd"]),
                    color=PALETTE["curve_sample"], lw=0.3, alpha=0.25)
        ac_vals = df_conv["HRD_crossing_age"].dropna().values
        if len(ac_vals):
            ax.axvspan(np.percentile(ac_vals, 5), np.percentile(ac_vals, 95),
                       color=PALETTE["hrd_band"], alpha=0.15)
            ax.axvline(np.median(ac_vals), color=PALETTE["median_curve"], lw=0.8)

    ax.scatter([age], [y_obs], color="k", s=20, zorder=10)
    ax.set(xlim=(0, 90), ylim=(0, 0.15), xlabel="Age [yr]",
           ylabel="SBS1 / Gx3000 [Mb]", title="Cumulative SBS1 burden")
    ax.tick_params(labelsize=7)


def _draw_sigmoid_panel_labeled(ax, df_conv, age, n_show):
    result = draw_sigmoid_panel(ax, df_conv, age, n_show=n_show)
    if result is not None:
        med, x_d_med, m_med, k_med, s_n_med, s_l_med = result
        fold = med["rd"] / s_n_med if s_n_med > 0 else np.nan
        txt = (f"m={m_med:.1f}yr  k={k_med:.2f}  dur={med['dur']:.1f}yr\n"
               f"x_d={x_d_med:.1f}yr  ctr={age - med['dur'] / 2:.1f}yr\n"
               f"s_n={s_n_med:.4f}  x{fold:.1f} -> s_late")
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, fontsize=6, va="top",
                ha="left", family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8, lw=0.3))
        leg = [
            Line2D([0], [0], color=PALETTE["initiation"], lw=1, ls="--", label="Initiation (A-d)"),
            Line2D([0], [0], color=PALETTE["midpoint"], lw=1, ls="-", label="Midpoint m"),
            Line2D([0], [0], color=PALETTE["diagnosis"], lw=1, ls=":", label="Diagnosis A"),
            Line2D([0], [0], color=PALETTE["median_curve"], lw=1.5, label="Median curve"),
        ]
        ax.legend(handles=leg, fontsize=5.5, frameon=False, loc="upper left",
                  bbox_to_anchor=(0, 0.62))
    ax.set_xlim(0, age * 1.08)
    ax.set(xlabel="Age [yr]", ylabel="dy/dx [SBS1/yr]", title="Sigmoid rate transitions")
    ax.tick_params(labelsize=7)


def _draw_k_distribution_panel(ax, df_conv):
    if not df_conv.empty:
        k_vals = df_conv["k"].dropna().values
        kmin_vals = df_conv["k_min_iter"].dropna().values
        ax.hist(k_vals, bins=30, color=PALETTE["diagnosis"], edgecolor="white",
                lw=0.3, alpha=0.8, label="sampled k")
        ax.axvline(np.median(k_vals), color=PALETTE["median_curve"], lw=1.0, ls="--",
                   label=f"median={np.median(k_vals):.3f}")
        ax.axvline(np.median(kmin_vals), color=PALETTE["initiation"], lw=1.0, ls=":",
                   label=f"k_min={np.median(kmin_vals):.3f}")
        ax.legend(fontsize=6, frameon=False)
    ax.set(xlabel="k", ylabel="Count", title="Sampled k distribution")
    ax.tick_params(labelsize=7)


def _draw_p_combined_panel(ax, df_conv, stype, dur_mean_map, dur_sd_map):
    if not df_conv.empty:
        p_dur = df_conv["dur"].apply(
            lambda d: sp_norm.cdf(d, loc=dur_mean_map[stype], scale=dur_sd_map[stype]))
        p_excess = df_conv["rel_excess"].clip(0, 40) / 40.0
        p_combined = (1 - p_dur) * p_excess
        ax.hist(p_combined.dropna(), bins=30, color=STYPE_COLORS.get(stype, "steelblue"),
                edgecolor="white", lw=0.3, alpha=0.8)
        ax.axvline(p_combined.median(), color="black", lw=1.0, ls="--",
                   label=f"median={p_combined.median():.2f}")
        ax.legend(fontsize=6, frameon=False)
    ax.set(xlim=(0, 1), xlabel="p_combined", ylabel="Count", title="p_combined distribution")
    ax.tick_params(labelsize=7)


# ── Per-sample figures ──────────────────────────────────────────────────────
def make_sample_figure(
    row: pd.Series, df_conv: pd.DataFrame, slope_dist: SlopeDistribution,
    dur_mean_map: Dict, dur_sd_map: Dict, n_show: int = 40,
) -> plt.Figure:
    """4-panel diagnostic figure for the gradual (sigmoid) model."""
    age = float(row["age"])
    y_obs = float(row["scaled_SBS1"]) + _LATE
    stype = row["type"]
    sample_id = row["sample.display"]
    hrd_f = float(row["HRDTime"])
    hrd_sd = (float(row["HRDTime_ci_hi"]) + float(row["HRDTime_ci_lo"])) / (2 * 1.96)

    fig, axes = plt.subplots(1, 4, figsize=(12, 3.2))
    fig.suptitle(f"{sample_id}  ({stype})  -  age={age:.1f} yr  [Gradual model]",
                 fontsize=10, y=1.01)
    _draw_burden_panel(axes[0], df_conv, age, y_obs, hrd_f, hrd_sd, slope_dist)
    _draw_sigmoid_panel_labeled(axes[1], df_conv, age, n_show)
    _draw_k_distribution_panel(axes[2], df_conv)
    _draw_p_combined_panel(axes[3], df_conv, stype, dur_mean_map, dur_sd_map)
    plt.tight_layout()
    return fig


def make_two_step_sample_figure(
    row: pd.Series, df_conv: pd.DataFrame, slope_dist: SlopeDistribution,
) -> plt.Figure:
    """3-panel diagnostic figure for the two-step model."""
    age = float(row["age"])
    y_obs = float(row["scaled_SBS1"]) + _LATE
    stype = row["type"]
    sample_id = row["sample.display"]
    hrd_f = float(row["HRDTime"])
    hrd_sd = (float(row["HRDTime_ci_hi"]) + float(row["HRDTime_ci_lo"])) / (2 * 1.96)

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2))
    fig.suptitle(f"{sample_id}  ({stype})  -  age={age:.1f} yr  [Two-step model]",
                 fontsize=10, y=1.01)
    draw_two_step_burden_panel(axes[0], df_conv, age, y_obs, hrd_f, hrd_sd, slope_dist)

    ax = axes[1]
    if not df_conv.empty:
        ts_vals = df_conv["TS_crossing_age"].dropna().values
        if len(ts_vals):
            ax.hist(ts_vals, bins=30, color=PALETTE["two_step"], edgecolor="white",
                    lw=0.3, alpha=0.8)
            ax.axvline(np.median(ts_vals), color="black", lw=1.0, ls="--",
                       label=f"median={np.median(ts_vals):.1f}yr")
            ax.axvline(np.percentile(ts_vals, 5), color="black", lw=0.6, ls=":",
                       alpha=0.6, label="5-95th pct")
            ax.axvline(np.percentile(ts_vals, 95), color="black", lw=0.6, ls=":", alpha=0.6)
            ax.legend(fontsize=6, frameon=False)
    ax.set(xlabel="HRD crossing age [yr]", ylabel="Count", title="TS crossing age distribution")
    ax.tick_params(labelsize=7)

    ax = axes[2]
    if not df_conv.empty:
        dur_vals = df_conv["TS_Duration"].dropna().values
        if len(dur_vals):
            ax.hist(dur_vals, bins=30, color=PALETTE["two_step"], edgecolor="white",
                    lw=0.3, alpha=0.8)
            ax.axvline(np.median(dur_vals), color="black", lw=1.0, ls="--",
                       label=f"median={np.median(dur_vals):.1f}yr")
            ax.legend(fontsize=6, frameon=False)
    ax.set(xlabel="HRD duration [yr]", ylabel="Count", title="TS HRD duration distribution")
    ax.tick_params(labelsize=7)

    plt.tight_layout()
    return fig


# ── Cohort grid figures ─────────────────────────────────────────────────────
def plot_cohort_burden(
    cohort: str, cdata: pd.DataFrame, conv_map: Dict[str, pd.DataFrame],
    slope_dist: SlopeDistribution, n_cols: int = 5,
) -> plt.Figure:
    """Grid of cumulative SBS1 burden (gradual model) for a cohort."""
    n_rows = int(np.ceil(len(cdata) / n_cols))
    set_mpl_style()
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 1.8, n_rows * 2.2), sharey=True)
    fig.suptitle(cohort, fontsize=13, y=0.975)
    axes_flat = np.atleast_1d(axes).flatten()

    i = -1
    for i, (_, row) in enumerate(cdata.iterrows()):
        ax = axes_flat[i]
        df_it = conv_map[row["sample.display"]]
        age = float(row["age"])
        y_obs = float(row["scaled_SBS1"]) + _LATE
        hrd_f = float(row["HRDTime"])
        hrd_sd = (float(row["HRDTime_ci_hi"]) + float(row["HRDTime_ci_lo"])) / (2 * 1.96)
        y_tgt = hrd_f * y_obs
        y_tgt_l = np.clip(hrd_f - hrd_sd, 0, 1) * y_obs
        y_tgt_h = np.clip(hrd_f + hrd_sd, 0, 1) * y_obs
        x_full = np.linspace(0, 90, 400)
        x_samp = np.linspace(0, age, 200)

        ax.fill_between(x_full, slope_dist.low * x_full, slope_dist.high * x_full,
                        color=PALETTE["normal_band"], alpha=0.15)
        ax.plot(x_full, slope_dist.mean * x_full, ":", color="black", lw=0.5, alpha=0.6)
        ax.fill_between([0, 90], [y_tgt_l] * 2, [y_tgt_h] * 2, color="orange", alpha=0.15)
        ax.plot([0, 90], [y_tgt] * 2, color="darkorange", ls="--", lw=0.4, alpha=0.6)

        if not df_it.empty:
            for ii in np.random.choice(len(df_it), min(50, len(df_it)), replace=False):
                r = df_it.iloc[ii]
                ax.plot(x_samp, sigmoid_curve(x_samp, r["mid"], r["k"], r["s_n"], r["rd"]),
                        color=PALETTE["curve_sample"], lw=0.3, alpha=0.25)
            ac_vals = df_it["HRD_crossing_age"].dropna().values
            if len(ac_vals):
                ax.axvspan(np.percentile(ac_vals, 5), np.percentile(ac_vals, 95),
                           color=PALETTE["hrd_band"], alpha=0.15)
                ax.axvline(np.median(ac_vals), color=PALETTE["median_curve"], lw=0.5)

        ax.scatter([age], [y_obs], color="k", s=15, zorder=10)
        ax.set(xlim=(0, 90), ylim=(0, 0.15), xlabel="Age [yr]",
               ylabel="SBS1 / Gx3000 [Mb]", title=f"{row['sample.display']} ({row['type']})")
        ax.tick_params(labelsize=7)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def plot_cohort_sigmoid(
    cohort: str, cdata: pd.DataFrame, conv_map: Dict[str, pd.DataFrame], n_cols: int = 5,
) -> plt.Figure:
    """Grid of sigmoid rate transitions for a cohort."""
    n_rows = int(np.ceil(len(cdata) / n_cols))
    set_mpl_style()
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 2.8, n_rows * 3.0), sharey=False)
    fig.suptitle(f"{cohort}  -  sigmoid rate transitions", fontsize=11, y=0.975)
    axes_flat = np.atleast_1d(axes).flatten()

    i = -1
    for i, (_, row) in enumerate(cdata.iterrows()):
        ax = axes_flat[i]
        df_it = conv_map[row["sample.display"]]
        age = float(row["age"])
        stype = row["type"]

        if df_it.empty:
            ax.set_title(f"{row['sample.display']} ({stype})", fontsize=7)
            ax.axis("off")
            continue

        result = draw_sigmoid_panel(ax, df_it, age)
        if result is not None:
            med, x_d_med, m_med, k_med, s_n_med, _ = result
            fold = med["rd"] / s_n_med if s_n_med > 0 else np.nan
            txt = (f"m={m_med:.1f}yr  k={k_med:.2f}\n"
                   f"dur={med['dur']:.1f}yr  x_d={x_d_med:.1f}yr\n"
                   f"s_n={s_n_med:.4f}  x{fold:.1f}->s_late")
            ax.text(0.03, 0.97, txt, transform=ax.transAxes, fontsize=5.2, va="top",
                    ha="left", family="monospace",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", alpha=0.8, lw=0.3))

        ax.set_xlim(0, age * 1.08)
        ax.set(xlabel="Age [yr]", ylabel="dy/dx [SBS1/yr]",
               title=f"{row['sample.display']} ({stype})")
        ax.tick_params(labelsize=6)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.legend(handles=_sigmoid_legend_handles(), loc="lower center", ncol=5,
               fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.01))
    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    return fig


def plot_cohort_two_step(
    cohort: str, cdata: pd.DataFrame, conv_map: Dict[str, pd.DataFrame],
    slope_dist: SlopeDistribution, n_cols: int = 5,
) -> plt.Figure:
    """Grid of two-step cumulative burden curves for a cohort."""
    n_rows = int(np.ceil(len(cdata) / n_cols))
    set_mpl_style()
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 1.8, n_rows * 2.2), sharey=True)
    fig.suptitle(f"{cohort}  -  Two-step model", fontsize=13, y=0.975)
    axes_flat = np.atleast_1d(axes).flatten()

    i = -1
    for i, (_, row) in enumerate(cdata.iterrows()):
        ax = axes_flat[i]
        df_it = conv_map[row["sample.display"]]
        age = float(row["age"])
        y_obs = float(row["scaled_SBS1"]) + _LATE
        hrd_f = float(row["HRDTime"])
        hrd_sd = (float(row["HRDTime_ci_hi"]) + float(row["HRDTime_ci_lo"])) / (2 * 1.96)
        draw_two_step_burden_panel(ax, df_it, age, y_obs, hrd_f, hrd_sd, slope_dist, n_show=50)
        ax.set_title(f"{row['sample.display']} ({row['type']})", fontsize=7, pad=2)
        ax.tick_params(labelsize=6)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def plot_p_combined_grid(diag_df_conv: pd.DataFrame, n_cols: int = 6) -> plt.Figure:
    """Grid of p_combined histograms across all converged samples."""
    samples = diag_df_conv["sample"].unique()
    n_rows = int(np.ceil(len(samples) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 2.2, n_rows * 2.0))
    axes_flat = np.atleast_1d(axes).flatten()

    i = -1
    for i, sid in enumerate(samples):
        ax = axes_flat[i]
        d = diag_df_conv[diag_df_conv["sample"] == sid]
        vals = d["p_combined"].dropna().values
        ax.hist(vals, bins=30, color=STYPE_COLORS.get(d["type"].iloc[0], "steelblue"),
                edgecolor="white", lw=0.3, alpha=0.8)
        ax.axvline(np.median(vals), color="black", lw=0.8, ls="--")
        ax.set_title(sid, fontsize=7, pad=2)
        ax.set_xlim(0, 1)
        ax.tick_params(labelsize=6)
        if i % n_cols == 0:
            ax.set_ylabel("Count", fontsize=7)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.legend(
        handles=[plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.8, label=s)
                 for s, c in STYPE_COLORS.items()],
        loc="lower center", ncol=2, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02),
    )
    plt.suptitle("p_combined distribution per sample (converged iterations)", fontsize=11, y=1.01)
    plt.tight_layout()
    return fig


def plot_model_comparison(
    iter_df: pd.DataFrame, metric: str = "HRD_Duration", groupby: str = "type",
) -> plt.Figure:
    """Side-by-side boxplot comparing gradual vs two-step per-sample medians."""
    ts_col = "TS_Duration" if metric == "HRD_Duration" else "TS_crossing_age"
    groups = sorted(iter_df[groupby].dropna().unique())
    n_cols = len(groups)

    set_mpl_style()
    fig, axes = plt.subplots(1, n_cols, figsize=(n_cols * 3.5, 4.5), sharey=True)
    axes = np.atleast_1d(axes)
    ylabel = "HRD duration [yr]" if metric == "HRD_Duration" else "HRD crossing age [yr]"

    for ax, grp in zip(axes, groups):
        sub = iter_df[iter_df[groupby] == grp]
        gradual_medians = sub.groupby("sample")[metric].median().dropna().values
        ts_medians = sub.groupby("sample")[ts_col].median().dropna().values

        bp = ax.boxplot(
            [gradual_medians, ts_medians], positions=[1, 2], widths=0.5,
            patch_artist=True, medianprops=dict(color="white", lw=1.5),
            whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8),
            flierprops=dict(marker="o", markersize=3, alpha=0.4),
        )
        bp["boxes"][0].set_facecolor(PALETTE["median_curve"])
        bp["boxes"][1].set_facecolor(PALETTE["two_step"])

        for pos, vals, col in zip([1, 2], [gradual_medians, ts_medians],
                                  [PALETTE["median_curve"], PALETTE["two_step"]]):
            jitter = np.random.uniform(-0.12, 0.12, len(vals))
            ax.scatter(pos + jitter, vals, color=col, s=18, alpha=0.5, zorder=5,
                       edgecolors="white", lw=0.3)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Gradual", "Two-step"], fontsize=9)
        ax.set_title(grp, fontsize=10)
        ax.tick_params(labelsize=8)
        if ax is axes[0]:
            ax.set_ylabel(ylabel, fontsize=9)

    title_metric = "HRD Duration" if metric == "HRD_Duration" else "HRD Crossing Age"
    fig.suptitle(f"Gradual vs Two-step model - {title_metric} (per-sample medians)",
                 fontsize=11, y=1.02)
    fig.legend(
        handles=[mpatches.Patch(color=PALETTE["median_curve"], label="Gradual (sigmoid)"),
                 mpatches.Patch(color=PALETTE["two_step"], label="Two-step (piecewise)")],
        loc="lower center", ncol=2, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.04),
    )
    plt.tight_layout()
    return fig


def plot_hrd_onset(
    summary_df: pd.DataFrame,
    sample_order: Optional[List[str]] = None,
    meta_df: Optional[pd.DataFrame] = None,
    annotations: Optional[List[str]] = None,
    annotation_colors: Optional[Dict[str, List[str]]] = None,
    sep_lines: Optional[List[float]] = None,
    sep_labels: Optional[List[Tuple]] = None,
    title: str = "HRD onset vs age at diagnosis",
    ci_color: str = "#A7C7E7",
    lw: float = 0.25,
    fig_width_per_sample: float = 0.18,
    main_height: float = 3.5,
) -> plt.Figure:
    """Lollipop plot of HRD crossing-age (median + 90% CI) vs age at diagnosis.

    Parameters
    ----------
    summary_df:
        Output of ``AgeModel.summary_df``.  Must contain ``sample``,
        ``age_at_diagnosis``, ``HRD_crossing_age_median``,
        ``HRD_crossing_age_ci05``, ``HRD_crossing_age_ci95``.
    sample_order:
        Order samples appear on the x-axis.  Defaults to the order in
        *summary_df*.
    meta_df:
        Optional metadata DataFrame indexed by sample name (matching
        ``summary_df["sample"]``).  When supplied, an annotation colour-bar
        panel is added below the main panel.
    annotations:
        Columns from *meta_df* to draw as annotation rows (top -> bottom).
    annotation_colors:
        ``{column: [color0, color1, ...]}`` mapping.  Values for each column
        are sorted alphabetically and assigned colours in order.
    sep_lines:
        x-positions of vertical separator lines (e.g. cohort boundaries).
    sep_labels:
        ``[(x_centre, label, colour), ...]`` text labels above separator lines.
    """
    col_med  = "HRD_crossing_age_median"
    col_low  = "HRD_crossing_age_ci05"
    col_high = "HRD_crossing_age_ci95"

    if sample_order is None:
        sample_order = summary_df["sample"].tolist()

    df = summary_df.set_index("sample").loc[sample_order].reset_index()
    n  = len(df)

    # Bar width shrinks with more samples, capped between 0.08 and 0.2
    bar_width = float(np.clip(0.6 / max(n, 1) ** 0.4, 0.08, 0.2))

    has_ann    = meta_df is not None and bool(annotations)
    n_ann_rows = len(annotations) if has_ann and annotations else 0
    ann_height = main_height / 3
    fig_h      = main_height + (ann_height if has_ann else 0)
    fig_w      = max(2, n * fig_width_per_sample) + 1.0  # +1 for right-side legends

    if has_ann:
        fig, axes = plt.subplots(
            2, 1, figsize=(fig_w, fig_h),
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        )
        ax_main, ax_ann = axes
    else:
        fig, ax_main = plt.subplots(figsize=(fig_w, fig_h))
        ax_ann = None

    # ── Main lollipop panel ──────────────────────────────────────────────────
    for i, row in df.iterrows():
        ax_main.plot([i, i], [row["age_at_diagnosis"], row[col_med]],
                     color="k", lw=0.5, zorder=1)
        ax_main.add_patch(plt.Rectangle(
            (i - bar_width / 2, row[col_low]),
            bar_width, row[col_high] - row[col_low],
            color=ci_color, zorder=0, linewidth=0,
        ))
        ax_main.plot(i, row["age_at_diagnosis"], "o",
                     color="grey", ms=4, zorder=2,
                     markeredgewidth=lw, markeredgecolor="grey",
                     label="Age at diagnosis" if i == 0 else "")
        ax_main.plot(i, row[col_med], "o",
                     color="#191970", ms=4, zorder=3,
                     markeredgewidth=lw, markeredgecolor="#191970",
                     label="HRD onset" if i == 0 else "")

    ax_main.set_xlim(-0.5, n - 0.5)
    ax_main.set_ylim(bottom=0)
    y_top = ax_main.get_ylim()[1]

    for sx in (sep_lines or []):
        ax_main.axvline(sx, color="black", lw=lw, zorder=5)
    for sx, sl, sc in (sep_labels or []):
        ax_main.text(sx, y_top, sl, ha="center", va="bottom", fontsize=7, color=sc)

    ax_main.set_xticks(range(n))
    ax_main.set_xticklabels([] if has_ann else sample_order, rotation=90, fontsize=5)
    ax_main.set_ylabel("Age (years)")
    ax_main.set_title(title, y=1.02, fontsize=9)
    for sp in ax_main.spines.values():
        sp.set_linewidth(lw)
    ax_main.tick_params(width=lw, labelsize=6)

    # Main legend — collected at the end so annotation patches can be appended
    _main_handles, _main_labels = ax_main.get_legend_handles_labels()

    # ── Annotation panel ─────────────────────────────────────────────────────
    if has_ann and ax_ann is not None and meta_df is not None and annotations:
        ann_colors = annotation_colors or {}
        ordered_vals: Dict[str, List] = {}
        for col in annotations:
            if col not in meta_df.columns:
                ordered_vals[col] = []
                continue
            uniq = sorted(v for v in meta_df[col].dropna().unique() if v != "otherHRD")
            if "otherHRD" in meta_df[col].values:
                uniq.append("otherHRD")
            ordered_vals[col] = uniq

        for i, sid in enumerate(sample_order):
            if sid not in meta_df.index:
                continue
            for j, col in enumerate(annotations):
                val = meta_df.at[sid, col] if col in meta_df.columns else np.nan
                uv  = ordered_vals[col]
                rc  = ann_colors.get(col, [])
                idx = uv.index(val) if pd.notna(val) and val in uv else -1
                color = rc[idx] if 0 <= idx < len(rc) else "white"
                ax_ann.add_patch(mpatches.Rectangle(
                    (i - 0.5, j), 1, 1,
                    facecolor=color, edgecolor="lightgrey", linewidth=0.0,
                ))

        for sx in (sep_lines or []):
            ax_ann.axvline(sx, color="black", lw=lw, zorder=5)

        ax_ann.set_xlim(-0.5, n - 0.5)
        ax_ann.set_ylim(0, len(annotations))
        ax_ann.invert_yaxis()
        ax_ann.set_xticks(range(n))
        ax_ann.set_xticklabels(sample_order, rotation=90, fontsize=5)
        ax_ann.set_yticks([j + 0.5 for j in range(len(annotations))])
        ax_ann.set_yticklabels(annotations, fontsize=7)
        for sp in ax_ann.spines.values():
            sp.set_linewidth(lw)
        ax_ann.tick_params(width=lw, labelsize=6)

        # Annotation legend — outside to the right, one section per annotation row
        leg_patches = []
        for col in annotations:
            if not ordered_vals[col]:
                continue
            leg_patches.append(mpatches.Patch(color="none", label=f"— {col} —"))
            for idx, val in enumerate(ordered_vals[col]):
                rc = ann_colors.get(col, [])
                color = rc[idx] if idx < len(rc) else "white"
                leg_patches.append(mpatches.Patch(color=color, label=str(val)))

        # Combined legend on ax_main: main items first, then annotation sections
        all_handles = list(_main_handles) + [
            mpatches.Patch(color="none", label="")  # spacer
        ] + leg_patches
        all_labels = list(_main_labels) + [""] + [
            p.get_label() for p in leg_patches
        ]

        leg_combined = ax_main.legend(
            handles=all_handles,
            labels=all_labels,
            loc="upper left", bbox_to_anchor=(1.01, 1.0),
            fontsize=5, frameon=True, borderaxespad=0, handlelength=1.0,
        )
        leg_combined.get_frame().set_linewidth(lw)
        for text, patch in zip(leg_combined.get_texts(), leg_combined.legend_handles):
            label = text.get_text()
            if label.startswith("— "):
                text.set_fontweight("bold")
                text.set_fontsize(6)
                patch.set_visible(False)
            elif label in _main_labels:
                text.set_fontsize(7)

    else:
        # No annotation panel — just draw the main legend
        h, lbl = ax_main.get_legend_handles_labels()
        leg = ax_main.legend(
            dict(zip(lbl, h)).values(), dict(zip(lbl, h)).keys(),
            loc="upper left", bbox_to_anchor=(1.01, 1.0),
            fontsize=7, frameon=True, borderaxespad=0,
        )
        leg.get_frame().set_linewidth(lw)

    #plt.tight_layout()
    return fig
