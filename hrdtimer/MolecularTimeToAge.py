import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import linregress

def calculate_line(x, slope, intercept):
    return slope * x + intercept

def plot_filled_lines(ax, y, WGD, HRD, W_err, H_err):
    y_val = y.values[0]
    for t, e, c in [(WGD, W_err, '#e68f32'), (HRD, H_err, '#21255e')]:
        for s in [1, -1]:
            ax.fill_betweenx(y=[t * y_val, (t + s * e) * y_val], x1=ax.get_xlim()[0], x2=85, color=c, alpha=0.5)
        ax.axhline(y=t * y_val, color=c, linestyle='--')

def polynomial_model(x, *params):
    powers = [1,2,3,4,5,6,7,8,9,10,11,12,13,20,22,23,24,40]
    return sum(p * (x ** pow) for p, pow in zip(params, powers))

def normal_line(x, x_breast, y_breast):
    return 0.0002569 * x * x_breast / y_breast

def max_acc_line(x, slope, intercept, x_breast, y_breast):
    return slope * x * (x_breast / y_breast) + intercept / y_breast

def objective(params, x_breast, y_breast, slope, intercept):
    x = np.linspace(0, 1, 10000)
    y_model = polynomial_model(x, *params)
    nl = normal_line(x, x_breast, y_breast)
    ml = max_acc_line(x, slope, intercept, x_breast, y_breast)
    penalties = (
        np.sum(np.clip(-y_model, 0, None)) +
        np.sum(np.clip(nl - y_model, 0, None)) * 10 +
        np.sum(np.clip(ml - y_model, 0, None)) * 10 +
        np.sum(np.clip(y_model[x > 0.9] - ml[x > 0.9], 0, None)) +
        np.sum((y_model[x <= 0.2] - nl[x <= 0.2])**2) * 100
    )
    return penalties

def analyze_cohort_timing(cohort_name, bulk_tumor_df, hrd_timer_df, normal_tissue_df, output_csv, output_fig=None):
    samples = bulk_tumor_df[bulk_tumor_df['Cohort'] == cohort_name]
    merged = pd.merge(samples, hrd_timer_df, left_on='aliquot_id', right_on='ID')
    merged = merged.drop(columns='ID')

    x_norm, y_norm = normal_tissue_df['age'], normal_tissue_df['scaled_SBS1']
    slope_norm, intercept_norm, *_ = linregress(x_norm, y_norm)
    x_values_norm = np.linspace(0, 90, 1000)
    y_values_norm = slope_norm * x_values_norm

    IDs = merged['aliquot_id'].unique()
    nrows, ncols = -(-len(IDs) // 5), 5
    fig, axs = plt.subplots(nrows, ncols, figsize=(20, nrows * 4))
    axs = axs.flatten()
    results = pd.DataFrame(columns=[
        'ID','Age',
        'WGD','WGD_low','WGD_high',
        'HRD','HRD_low','HRD_high',
        'WGD_linear','WGD_linear_low','WGD_linear_high',
        'HRD_linear','HRD_linear_low','HRD_linear_high'
    ])

    for idx, sid in enumerate(IDs):
        d = merged[merged['aliquot_id'] == sid]
        if d['scaled_SBS1'].empty: continue

        WGD, HRD = d[['WGDTime', 'HRDTime']].values[0]
        W_err, H_err = d[['WGDTime_ci', 'HRDTime_ci']].values[0]
        x_breast, y_breast = d['age'].values[0], d['scaled_SBS1']
        ax = axs[idx]

        plot_filled_lines(ax, y_breast, WGD, HRD, W_err, H_err)

        slope = 30 * 0.0002569
        intercept = -slope * x_breast + y_breast
        ax.plot(x_values_norm, calculate_line(x_values_norm, slope, intercept.mean()), '--', color='k')
        ax.plot(x_values_norm, y_values_norm, color='k')
        ax.scatter(x_breast, y_breast, color='k')

        init_guess = [0.1] * 18
        bnds = [(0, None)] * 18
        res = minimize(objective, init_guess, args=(x_breast, y_breast.values[0], slope, intercept.mean()), method='L-BFGS-B', bounds=bnds)
        params = res.x

        x = np.linspace(0, 1, 10000)
        x_scaled = x * x_breast
        y_scaling = y_breast.values[0]
        y_fit = polynomial_model(x, *params) * y_scaling
        ax.plot(x_scaled, y_fit, '-', linewidth=0.6, color='k')

        def get_mean_x(target):
            ix = np.where(np.abs(y_fit - target) < 1e-3)[0]
            return np.mean(x_scaled[ix]) if ix.size else "NA"
        
        # Polynomial
        HRD_pow = get_mean_x(float(HRD * y_breast))
        HRD_high = get_mean_x(float((HRD + H_err) * y_breast))
        HRD_low = get_mean_x(float((HRD - H_err) * y_breast)) or 0
        WGD_pow = get_mean_x(float(WGD * y_breast))
        WGD_high = get_mean_x(float((WGD + W_err) * y_breast))
        WGD_low = get_mean_x(float((WGD - W_err) * y_breast))

        # Linear model-based
        HRD_linear = "NA" if (HRD * x_breast) < 0 else float(HRD * x_breast)
        HRD_linear_high = "NA" if ((HRD + H_err) * x_breast) < 0 else float((HRD + H_err) * x_breast)
        HRD_linear_low  = "NA" if ((HRD - H_err) * x_breast) < 0 else float((HRD - H_err) * x_breast)
        WGD_linear = "NA" if (WGD * x_breast) < 0 else float(WGD * x_breast)
        WGD_linear_high = "NA" if ((WGD + W_err) * x_breast) < 0 else float((WGD + W_err) * x_breast)
        WGD_linear_low = "NA" if ((WGD - W_err) * x_breast) < 0 else float((WGD - W_err) * x_breast)

        results.loc[len(results)] = [
            sid, x_breast,
            WGD_pow, WGD_low, WGD_high,
            HRD_pow, HRD_low, HRD_high,
            WGD_linear, WGD_linear_low, WGD_linear_high,
            HRD_linear, HRD_linear_low, HRD_linear_high
        ]

        ax.set(title=sid, xlim=(0, 85), ylim=(0, 0.1), xlabel='Age', ylabel='SBS1 / G x 3000 Mb')

    results.to_csv(output_csv, index=False)
    plt.tight_layout()
    #if output_fig:
        #plt.savefig(output_fig, dpi=400)
    plt.show()