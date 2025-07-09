import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_filled_lines(ax, y, WGD, HRD, W_err_hi, W_err_lo, H_err_hi, H_err_lo):
    y_val = y.values[0]
    ax.fill_betweenx(
        y=[(WGD - W_err_lo) * y_val, (WGD + W_err_hi) * y_val],
        x1=ax.get_xlim()[0], x2=85, color='#e68f32', alpha=0.5
    )
    ax.axhline(y=WGD * y_val, color='#e68f32', linestyle='--')

    ax.fill_betweenx(
        y=[(HRD - H_err_lo) * y_val, (HRD + H_err_hi) * y_val],
        x1=ax.get_xlim()[0], x2=85, color='#21255e', alpha=0.5
    )
    ax.axhline(y=HRD * y_val, color='#21255e', linestyle='--')
    
def MolecularTime_to_Age(SBS1_exposures, moleculartimedf, metadata, late_SBS1_burden=0, output_fig=None):
    """
    Estimate chronological timing of WGD and HRD events using SBS1 molecular clock data.

    This function merges molecular timing data with SBS1 exposures and metadata, fits multiple 
    slope models to approximate SBS1 accumulation over age, and computes estimated ages at 
    WGD and HRD events under different assumptions (linear, intermediate, and maximum slopes). 
    Optionally, it generates a multi-panel plot showing the fits for all samples.

    Parameters
    ----------
    SBS1_exposures : pd.DataFrame
        DataFrame containing SBS1 exposure information with at least the following columns:
        - 'ID': sample identifier
        - 'scaled_SBS1': scaled SBS1 mutation burden

    moleculartimedf : pd.DataFrame
        DataFrame containing molecular timing data for each sample, including columns:
        - 'ID': sample identifier
        - 'WGDTime', 'WGDTime_ci_lo', 'WGDTime_ci_hi': timing and confidence intervals for WGD event
        - 'HRDTime', 'HRDTime_ci_lo', 'HRDTime_ci_hi': timing and confidence intervals for HRD event
        - Other molecular timing related columns

    metadata : pd.DataFrame
        DataFrame with sample metadata, providing at least:
        - 'sample': sample identifier (used to merge with SBS1_exposures 'ID')
        - 'age': chronological age of the sample
        - 'type': subtype classification (e.g., 'TN', 'ER+')

    late_SBS1_burden : float, optional (default=0)
        Additional SBS1 burden to add to molecular clock estimates to account for late SBS1 accumulation 
        (e.g., post-WGD accumulation). Estimated from scNanoSeq

    output_fig : str or None, optional (default=None)
        Path to save the output plot (PDF). If None, the figure is not saved.

    Returns
    -------
    pd.DataFrame
        DataFrame containing estimated ages for WGD and HRD events under multiple timing models:
        linear, intermediate, and maximum slopes, including their confidence intervals.
    """

    results = []
    IDs = moleculartimedf['ID'].unique()
    nrows, ncols = -(-len(IDs) // 5), 5
    fig, axs = plt.subplots(nrows, ncols, figsize=(20, nrows * 4))
    axs = axs.flatten()

    def intersect_lines(m1, b1, m2, b2):
        return (b2 - b1) / (m1 - m2) if abs(m1 - m2) >= 1e-14 else np.nan

    def min_positive(*args):
        pos = [x for x in args if x is not None and not pd.isna(x)]
        return max(0, min(pos)) if pos else np.nan

    def zero_floor(val):
        return max(0, val) if pd.notna(val) else val

    x_vals = np.linspace(0, 90, 1000)
    ref_slope = 0.0002569
    slope_max = 30 * ref_slope

    SBS1_exposures = SBS1_exposures.drop(columns=['age', 'type'], errors='ignore')
    SBS1_exposures = SBS1_exposures.merge(
        metadata[['sample', 'age', 'type']],
        left_on='ID', right_on='sample',
        how='left'
    ).drop(columns=['sample'])

    merged_df = moleculartimedf.merge(SBS1_exposures, on='ID', how='inner')
    merged_df = merged_df[merged_df['age'].notna() & (merged_df['age'] > 0)]


    for idx, row in merged_df.iterrows():
        sid = row['ID']
        x_breast = row['age']
        y_breast = row['scaled_SBS1']
        y_breast_corrected = y_breast + late_SBS1_burden

        WGD = row['WGDTime']
        HRD = row['HRDTime']
        W_err_hi = row['WGDTime_ci_hi']
        W_err_lo = row['WGDTime_ci_lo']
        H_err_hi = row['HRDTime_ci_hi']
        H_err_lo = row['HRDTime_ci_lo']

        ax = axs[idx]

        plot_filled_lines(ax, pd.Series(y_breast), WGD, HRD, W_err_hi, W_err_lo, H_err_hi, H_err_lo)
        ax.plot(x_vals, ref_slope * x_vals, color='grey', linestyle='-')

        timing_class = row['type'] if 'type' in row else 'Unknown'
        duration = 10.5 if timing_class == 'TN' else 16.3 if timing_class == 'ER+' else 0
        x_div = max(x_breast - duration, 1e-3)
        y_div = ref_slope * x_div

        slope_intermediate = (y_breast_corrected - y_div) / (x_breast - x_div)
        intercept_intermediate = y_breast_corrected - slope_intermediate * x_breast
        intercept_max = y_breast_corrected - slope_max * x_breast
        straight_slope = y_breast_corrected / x_breast

        ax.plot(x_vals, slope_intermediate * x_vals + intercept_intermediate, '--', color='#658E9C')
        ax.plot(x_vals, slope_max * x_vals + intercept_max, '--', color='k')
        ax.plot(x_vals, straight_slope * x_vals, '--', color='grey')

        ax.scatter(x_breast, y_breast_corrected, color='#C9A690', edgecolors='black',
                   linewidths=0.5, s=70, zorder=20, label='SBS1 incl. late')
        ax.scatter(x_breast, y_breast, color='#D36582', edgecolors='black',
                   linewidths=0.5, s=30, zorder=20, label='SBS1 at MRCA')
        ax.scatter(x_div, y_div, color='#658E9C', s=70, edgecolors='black',
                   linewidths=0.5, zorder=20, label='Divergence')

        y_scenarios = {
            'HRD': HRD * y_breast,
            'HRD_low': (HRD - H_err_lo) * y_breast,
            'HRD_high': (HRD + H_err_hi) * y_breast,
            'WGD': WGD * y_breast,
            'WGD_low': (WGD - W_err_lo) * y_breast,
            'WGD_high': (WGD + W_err_hi) * y_breast
        }

        def compute_intersections(y_val):
            return {
                'linear': y_val / straight_slope if straight_slope != 0 else np.nan,
                'norm_int': min_positive(
                    intersect_lines(ref_slope, 0, 0, y_val),
                    intersect_lines(slope_intermediate, intercept_intermediate, 0, y_val)
                ),
                'norm_max': min_positive(
                    intersect_lines(ref_slope, 0, 0, y_val),
                    intersect_lines(slope_max, intercept_max, 0, y_val)
                )
            }

        HRD_results = {k: compute_intersections(v) for k, v in y_scenarios.items() if 'HRD' in k}
        WGD_results = {k: compute_intersections(v) for k, v in y_scenarios.items() if 'WGD' in k}

        results.append({
            'ID': sid,
            'Age': x_breast,
            'HRD_linear': zero_floor(HRD_results['HRD']['linear']),
            'HRD_linear_low': zero_floor(HRD_results['HRD_low']['linear']),
            'HRD_linear_high': zero_floor(HRD_results['HRD_high']['linear']),
            'WGD_linear': zero_floor(WGD_results['WGD']['linear']),
            'WGD_linear_low': zero_floor(WGD_results['WGD_low']['linear']),
            'WGD_linear_high': zero_floor(WGD_results['WGD_high']['linear']),
            'HRD_norm_int': zero_floor(HRD_results['HRD']['norm_int']),
            'HRD_norm_int_low': zero_floor(HRD_results['HRD_low']['norm_int']),
            'HRD_norm_int_high': zero_floor(HRD_results['HRD_high']['norm_int']),
            'WGD_norm_int': zero_floor(WGD_results['WGD']['norm_int']),
            'WGD_norm_int_low': zero_floor(WGD_results['WGD_low']['norm_int']),
            'WGD_norm_int_high': zero_floor(WGD_results['WGD_high']['norm_int']),
            'HRD_norm_max': zero_floor(HRD_results['HRD']['norm_max']),
            'HRD_norm_max_low': zero_floor(HRD_results['HRD_low']['norm_max']),
            'HRD_norm_max_high': zero_floor(HRD_results['HRD_high']['norm_max']),
            'WGD_norm_max': zero_floor(WGD_results['WGD']['norm_max']),
            'WGD_norm_max_low': zero_floor(WGD_results['WGD_low']['norm_max']),
            'WGD_norm_max_high': zero_floor(WGD_results['WGD_high']['norm_max']),
        })

        ax.set_title(f"{sid} ({timing_class})", fontsize=9)
        ax.set_xlim(0, 85)
        ax.set_ylim(0, 0.13)
        ax.set_xlabel('Age')
        ax.set_ylabel('SBS1 / G x 3000 Mb')

    plt.tight_layout()

    if output_fig:
        plt.savefig(output_fig, format='pdf', bbox_inches='tight', dpi=600)

    plt.show()

    return pd.DataFrame(results)

