# Standard Library
import os
import sys

# Third-Party Libraries
import numpy as np
import pandas as pd
from tqdm import tqdm

# Domain-Specific Libraries
import musical

# -------------------------- HRDTimer pipeline (bootstrapping) ---------------------------------
# ------------ Changing only signature posterior probabilities in each bootstrap ---------------

def generate_bootstraps(samples_dict, n_bootstraps, output_dir, tumor_type='Breast.AdenoCA'):
    """
    Performs bootstrap resampling to estimate signature exposure variance and mutation attribution.

    For each bootstrap replicate, this function:
    1. Resamples mutations per sample and classification group (Early, Late, NA) using multinomial sampling.
    2. Infers signature exposures for the resampled data using MuSiCal's likelihood refitting.
    3. Calculates posterior mutation probabilities (bootstrap-specific) for each signature.
    4. Saves annotated mutation lists and exposure matrices for each replicate.

    Args:
        samples_dict (dict[str, pd.DataFrame]): Dictionary mapping sample IDs to DataFrames 
            containing 'SBS96' and 'Classification' columns.
        n_bootstraps (int): Number of bootstrap iterations to perform.
        output_dir (str): Base directory for saving bootstrap results.
        tumor_type (str): Tumor type for MuSiCal catalog restriction. Defaults to 'Breast.AdenoCA'.

    Returns:
        None: Results are written to disk under 'bootstrap_i' subdirectories.

    Output Structure:
        output_dir/
        └── bootstrap_i/
            ├── [sample_id].csv        # Annotated mutations with 'prob_[Sig]_boot' columns
            ├── exposures_Early.csv    # Combined exposure matrix for Early mutations
            ├── exposures_Late.csv     # Combined exposure matrix for Late mutations
            └── exposures_NA.csv       # Combined exposure matrix for NA mutations
    """
    # Initialize MuSiCal Catalog
    catalog = musical.load_catalog('COSMIC_v3p2_SBS_WGS')
    catalog.restrict_catalog(tumor_type=tumor_type, is_MMRD=False, is_PPD=False)
    mutation_types = catalog.W.index.tolist()
    W = catalog.W.reindex(mutation_types)
    
    os.makedirs(output_dir, exist_ok=True)
    labels = ['Early', 'Late', 'NA']
    
    for i in tqdm(range(1, n_bootstraps + 1), desc="Bootstrapping", file=sys.stdout):
        boot_dir = os.path.join(output_dir, f'bootstrap_{i}')
        os.makedirs(boot_dir, exist_ok=True)

        # Dictionary to accumulate exposures for the entire cohort in this bootstrap
        group_exposure_dict = {label: [] for label in labels}
        
        for sample_id, df in samples_dict.items():
            sample_merged_data = []

            for label in labels:
                group = df[df['Classification'] == label].copy()
                if group.empty:
                    continue
                
                # 1. Multinomial Sampling
                counts = group['SBS96'].value_counts()
                type_probs = counts.reindex(mutation_types, fill_value=0) / len(group)
                
                sampled_counts = np.random.multinomial(len(group), type_probs.values)
                
                # Format for MuSiCal refit (expects DataFrame with samples as columns)
                sampled_df = pd.DataFrame({sample_id: sampled_counts}, index=mutation_types)
                sampled_df.index.name = 'Type'

                # 2. MuSiCal Refit
                exposures, _ = musical.refit.refit(sampled_df, W, method='likelihood_bidirectional', thresh=0.001)
                
                # Capture exposures for cohort-level file
                exp_series = exposures.iloc[:, 0]
                exp_series.name = sample_id
                group_exposure_dict[label].append(exp_series)

                # 3. Calculate Mutation Attribution Probabilities
                exposures_norm = exp_series / exp_series.sum() if exp_series.sum() > 0 else exp_series
                W_prob = W.mul(exposures_norm, axis=1)
                
                # Normalize across signatures (rows sum to 1)
                prob_df = W_prob.div(W_prob.sum(axis=1), axis=0).fillna(0)
                prob_df = prob_df.rename(columns=lambda c: f'prob_{c}_boot')

                # 4. Map probabilities back to the original group mutations
                group = group.join(prob_df, on='SBS96')
                sample_merged_data.append(group)

            # Save individual sample file for this bootstrap
            if sample_merged_data:
                pd.concat(sample_merged_data, ignore_index=True).to_csv(
                    os.path.join(boot_dir, f"{sample_id}.csv"), index=False
                )

        # Save cohort-level exposure files for this bootstrap
        for label, exposure_list in group_exposure_dict.items():
            if exposure_list:
                df_exposures = pd.concat(exposure_list, axis=1).T
                df_exposures.to_csv(os.path.join(boot_dir, f"exposures_{label}.csv"))

def calculate_wgd_time_bootstrapping(sample_id, base_dir, num_bootstrap=200):
    """
    Estimates the timing of Whole Genome Doubling (WGD) using bootstrap replicates.

    This function calculates WGD timing based on the molecular clock of CpG > TpG 
    mutations (SBS1). It computes a weighted mean of time estimates across different 
    minor copy number states (MinCN 0, 1, 2) for each bootstrap, then aggregates 
    them to provide a final estimate with confidence intervals and IQR.

    Args:
        sample_id (str): Unique identifier for the sample.
        base_dir (str): Path to the directory containing 'bootstrap_i' subfolders.
        num_bootstrap (int): Number of bootstrap replicates to process. Defaults to 200.

    Returns:
        tuple: A collection of timing metrics:
            - N_mut_CpG_total (np.array): Total CpG mutations found for MinCN 0, 1, 2.
            - weighted_means (np.array): The WGD time estimate for each bootstrap.
            - WGD_time (float): The mean WGD time across all bootstraps.
            - WGD_time_CI_hi/lo (float): 95% Confidence Interval bounds.
            - WGD_time_IQR_hi/lo (float): Interquartile Range bounds.
    """
    weighted_means = []
    N_mut_CpG_total = np.zeros(3)
    cpg_contexts = ['A[C>T]G', 'C[C>T]G', 'T[C>T]G', 'G[C>T]G']

    for i in range(1, num_bootstrap + 1):
        file_path = os.path.join(base_dir, f"bootstrap_{i}", f"{sample_id}.csv")
        if not os.path.exists(file_path):
            continue

        sample_df = pd.read_csv(file_path)
        # Filter for CpG mutations immediately
        sample_df = sample_df[sample_df['SBS96'].isin(cpg_contexts)]

        t_values_replicate = []
        n_mut_replicate = []

        for min_cn in range(3):
            sub = sample_df[sample_df['MinCN'] == min_cn]
            n_count = len(sub)
            n_mut_replicate.append(n_count)
            N_mut_CpG_total[min_cn] += n_count

            if n_count == 0:
                t_values_replicate.append(np.nan)
                continue

            # Vectorized calculation of pi_1 and pi_2
            # sum(prob * p) / sum(prob)
            probs = sub['prob_SBS1_boot']
            sum_probs = probs.sum()
            
            if sum_probs == 0:
                t_values_replicate.append(np.nan)
                continue

            pi_1 = (probs * sub['pSingle']).sum() / sum_probs
            pi_2 = (probs * sub['pGain']).sum() / sum_probs

            # Time calculation logic
            denom = pi_1 + 2 * pi_2
            if denom == 0:
                t_value = np.nan
            elif pi_2 == 0:
                t_value = 0
            else:
                multiplier = 3 if min_cn == 1 else 2
                t_value = (multiplier * pi_2) / denom
            
            t_values_replicate.append(t_value)

        # Calculate weighted mean for this bootstrap replicate
        t_values_replicate = np.array(t_values_replicate)
        n_mut_replicate = np.array(n_mut_replicate)
        
        mask = ~np.isnan(t_values_replicate)
        if np.any(mask):
            w_mean = np.average(t_values_replicate[mask], weights=n_mut_replicate[mask])
            weighted_means.append(w_mean)

    # Aggregate results across all bootstraps
    weighted_means = np.array(weighted_means)
    if len(weighted_means) == 0:
        return N_mut_CpG_total, np.array([]), np.nan, 0, 0, 0, 0

    wgd_time = np.nanmean(weighted_means)
    
    # Percentiles for CI and IQR
    p = np.nanpercentile(weighted_means, [2.5, 25, 75, 97.5])
    
    return (
        N_mut_CpG_total,
        weighted_means,
        wgd_time,
        p[3] - wgd_time,  # CI hi
        wgd_time - p[0],  # CI lo
        p[2] - wgd_time,  # IQR hi
        wgd_time - p[1]   # IQR lo
    )

def calculate_hrd_time_single_boot(sample_df):
    """
    Calculates HRD timing values for a sample across Minor Copy Number states.

    Args:
        sample_df (pd.DataFrame): DataFrame containing 'prob_SBS1_boot', 
            'prob_SBS3_boot', 'pSingle', 'pGain', and 'MinCN'.

    Returns:
        tuple: A collection of calculated metrics ordered by [MinCN 0, 1, 2]:
            t_values, N_mut, c_avg, c_list, pi_2_SBS1_vals, pi_2_SBS3_vals, 
            Nt_SBS1_vals, Nt_SBS3_vals, pi_1_SBS1_vals, pi_1_SBS3_vals.
    """
    if 'prob_SBS3_boot' not in sample_df.columns:
        sample_df['prob_SBS3_boot'] = 0

    # Initialize storage for indices [0, 1, 2]
    res = {k: {"t": np.nan, "n": 0, "c": None, "pi_2_SBS1": 0, "pi_2_SBS3": 0, 
               "Nt_SBS1": 0, "Nt_SBS3": 0, "pi_1_SBS1": 0, "pi_1_SBS3": 0} 
           for k in [0, 1, 2]}
    
    c_avg = 0

    # We process 0 and 2 first to calculate the weighted c_avg needed for state 1
    for min_cn in [0, 2, 1]:
        sub = sample_df[sample_df['MinCN'] == min_cn]
        n_count = len(sub)
        res[min_cn]["n"] = n_count
        
        if n_count == 0:
            continue

        # Vectorized weighted means for pi calculations
        def get_pi(prob_col, p_col):
            weights = sub[prob_col]
            total_w = weights.sum()
            return (weights * sub[p_col]).sum() / total_w if total_w > 0 else np.nan

        pi_1_SBS1 = get_pi('prob_SBS1_boot', 'pSingle')
        pi_2_SBS1 = get_pi('prob_SBS1_boot', 'pGain')
        pi_1_SBS3 = get_pi('prob_SBS3_boot', 'pSingle')
        pi_2_SBS3 = get_pi('prob_SBS3_boot', 'pGain')

        Nt_SBS1 = sub['prob_SBS1_boot'].sum()
        Nt_SBS3 = sub['prob_SBS3_boot'].sum()

        # <-- PASTE HERE
        pi_1_SBS1 = float(pi_1_SBS1) if pi_1_SBS1 is not None else np.nan
        pi_2_SBS1 = float(pi_2_SBS1) if pi_2_SBS1 is not None else np.nan
        pi_1_SBS3 = float(pi_1_SBS3) if pi_1_SBS3 is not None else np.nan
        pi_2_SBS3 = float(pi_2_SBS3) if pi_2_SBS3 is not None else np.nan
        Nt_SBS1 = float(Nt_SBS1) if Nt_SBS1 is not None else np.nan
        Nt_SBS3 = float(Nt_SBS3) if Nt_SBS3 is not None else np.nan


        # Calculation for MinCN 0 and 2
        if min_cn in [0, 2]:
            pi_2_SBS1_prime = pi_2_SBS1 - (pi_2_SBS3 / pi_1_SBS3 * pi_1_SBS1) if pi_1_SBS3 > 0 else np.nan
            if pi_1_SBS3 > 0 and Nt_SBS3 > 0:
                res[min_cn]["c"] = (pi_1_SBS1 * Nt_SBS1) / (pi_1_SBS3 * Nt_SBS3)
            
            # Update c_avg once both 0 and 2 are processed
            if min_cn == 2:
                num = sum((res[k]["c"] * res[k]["n"]) for k in [0, 2] if res[k]["c"] is not None)
                den = sum(res[k]["n"] for k in [0, 2] if res[k]["c"] is not None)
                c_avg = num / den if den > 0 else 0

        # Special Calculation for MinCN 1 using the weighted c_avg
        else: 
            pi_2_SBS1_prime = pi_2_SBS1 - (pi_2_SBS3 * c_avg * (Nt_SBS3 / Nt_SBS1)) if (Nt_SBS1 > 0 and pi_1_SBS3 > 0) else np.nan

        # Final t_value assignment
        denom = pi_1_SBS1 + 2 * pi_2_SBS1
        if not np.isnan(pi_2_SBS1_prime) and denom != 0:
            mult = 3 if min_cn == 1 else 2
            res[min_cn]["t"] = (mult * pi_2_SBS1_prime) / denom

        # Save intermediates back to dictionary
        res[min_cn].update({
            "pi_2_SBS1": pi_2_SBS1, "pi_2_SBS3": pi_2_SBS3, 
            "Nt_SBS1": Nt_SBS1, "Nt_SBS3": Nt_SBS3, 
            "pi_1_SBS1": pi_1_SBS1, "pi_1_SBS3": pi_1_SBS3
        })

    # Prepare flat lists for return to maintain compatibility
    order = [0, 1, 2]
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
        [res[k]["pi_1_SBS3"] for k in order]
    )

def calculate_hrd_time_single_boot(sample_df):
    """
    Calculates HRD timing values for a sample across Minor Copy Number states.
    """

    if 'prob_SBS3_boot' not in sample_df.columns:
        sample_df['prob_SBS3_boot'] = 0

    # Initialize storage for indices [0, 1, 2]
    res = {k: {"t": np.nan, "n": 0, "c": np.nan,  # <-- changed None -> np.nan
               "pi_2_SBS1": 0, "pi_2_SBS3": 0, 
               "Nt_SBS1": 0, "Nt_SBS3": 0, 
               "pi_1_SBS1": 0, "pi_1_SBS3": 0} 
           for k in [0, 1, 2]}
    
    c_avg = 0

    # We process 0 and 2 first to calculate the weighted c_avg needed for state 1
    for min_cn in [0, 2, 1]:
        sub = sample_df[sample_df['MinCN'] == min_cn]
        n_count = len(sub)
        res[min_cn]["n"] = n_count
        
        if n_count == 0:
            continue

        # Vectorized weighted means for pi calculations
        def get_pi(prob_col, p_col):
            weights = sub[prob_col]
            total_w = weights.sum()
            return (weights * sub[p_col]).sum() / total_w if total_w > 0 else np.nan

        pi_1_SBS1 = get_pi('prob_SBS1_boot', 'pSingle')
        pi_2_SBS1 = get_pi('prob_SBS1_boot', 'pGain')
        pi_1_SBS3 = get_pi('prob_SBS3_boot', 'pSingle')
        pi_2_SBS3 = get_pi('prob_SBS3_boot', 'pGain')

        Nt_SBS1 = sub['prob_SBS1_boot'].sum()
        Nt_SBS3 = sub['prob_SBS3_boot'].sum()

        pi_1_SBS1 = float(pi_1_SBS1) if pi_1_SBS1 is not None else np.nan
        pi_2_SBS1 = float(pi_2_SBS1) if pi_2_SBS1 is not None else np.nan
        pi_1_SBS3 = float(pi_1_SBS3) if pi_1_SBS3 is not None else np.nan
        pi_2_SBS3 = float(pi_2_SBS3) if pi_2_SBS3 is not None else np.nan
        Nt_SBS1 = float(Nt_SBS1) if Nt_SBS1 is not None else np.nan
        Nt_SBS3 = float(Nt_SBS3) if Nt_SBS3 is not None else np.nan

        # Calculation for MinCN 0 and 2
        if min_cn in [0, 2]:
            pi_2_SBS1_prime = pi_2_SBS1 - (pi_2_SBS3 / pi_1_SBS3 * pi_1_SBS1) if pi_1_SBS3 > 0 else np.nan
            if pi_1_SBS3 > 0 and Nt_SBS3 > 0:
                res[min_cn]["c"] = (pi_1_SBS1 * Nt_SBS1) / (pi_1_SBS3 * Nt_SBS3)
            
            # Update c_avg once both 0 and 2 are processed
            if min_cn == 2:
                # <-- changed to skip NaNs
                num = sum((res[k]["c"] * res[k]["n"]) for k in [0, 2] if not np.isnan(res[k]["c"]))
                den = sum(res[k]["n"] for k in [0, 2] if not np.isnan(res[k]["c"]))
                c_avg = num / den if den > 0 else 0

        # Special Calculation for MinCN 1 using the weighted c_avg
        else: 
            pi_2_SBS1_prime = pi_2_SBS1 - (pi_2_SBS3 * c_avg * (Nt_SBS3 / Nt_SBS1)) if (Nt_SBS1 > 0 and pi_1_SBS3 > 0) else np.nan

        # Final t_value assignment
        denom = pi_1_SBS1 + 2 * pi_2_SBS1
        if not np.isnan(pi_2_SBS1_prime) and denom != 0:
            mult = 3 if min_cn == 1 else 2
            res[min_cn]["t"] = (mult * pi_2_SBS1_prime) / denom

        # Save intermediates back to dictionary
        res[min_cn].update({
            "pi_2_SBS1": pi_2_SBS1, "pi_2_SBS3": pi_2_SBS3, 
            "Nt_SBS1": Nt_SBS1, "Nt_SBS3": Nt_SBS3, 
            "pi_1_SBS1": pi_1_SBS1, "pi_1_SBS3": pi_1_SBS3
        })

    # Prepare flat lists for return to maintain compatibility
    order = [0, 1, 2]
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
        [res[k]["pi_1_SBS3"] for k in order]
    )


def calculate_hrd_time_single_sample(sample_id, base_dir, num_bootstrap=200):
    """
    Aggregates HRD timing across multiple bootstrap replicates.

    This function iterates through bootstrap CSV files, calls the HRD timing 
    calculation for each, and aggregates the results to provide robust mean 
    estimates, confidence intervals (95%), and interquartile ranges (IQR). 
    It also summarizes intermediate parameters like mutation counts (Nt) 
    and proportions (pi) for both SBS1 and SBS3.

    Args:
        sample_id (str): Unique identifier for the sample.
        base_dir (str): Path to directory containing 'bootstrap_i' subfolders.
        num_bootstrap (int): Number of bootstrap replicates to process. Defaults to 200.

    Returns:
        tuple: A large collection of aggregated metrics including:
            - N_mut_all (list): Representative mutation counts per MinCN.
            - HRD_means (np.array): Weighted HRD time for every bootstrap replicate.
            - HRD_time (float): Mean HRD time across all replicates.
            - CI/IQR bounds (float): High/Low bounds for 95% CI and IQR.
            - Parameter means/errors (lists): Averaged pi and Nt values for SBS1 and SBS3.
    """
    # Initialization of storage structures
    HRD_means = []
    c_avg_val = []
    
    # Nested dict to store lists of intermediate values for each MinCN [0, 1, 2]
    stats = {
        'pi_2_SBS1': {0: [], 1: [], 2: []}, 'pi_2_SBS3': {0: [], 1: [], 2: []},
        'pi_1_SBS1': {0: [], 1: [], 2: []}, 'pi_1_SBS3': {0: [], 1: [], 2: []},
        'Nt_SBS1': {0: [], 1: [], 2: []}, 'Nt_SBS3': {0: [], 1: [], 2: []},
        'c_val': {0: [], 1: [], 2: []}
    }
    
    N_mut_all = None

    for i in range(1, num_bootstrap + 1):
        file_path = os.path.join(base_dir, f"bootstrap_{i}", f"{sample_id}.csv")
        if not os.path.exists(file_path):
            continue

        bootstrap_sample = pd.read_csv(file_path)

        # Call the core calculation function
        (t_values, N_mut, cavg, c_list, pi2_1, pi2_3, 
         nt1, nt3, pi1_1, pi1_3) = calculate_hrd_time_single_boot(bootstrap_sample)

        # 1. Calculate Weighted Mean Time for this bootstrap replicate
        t_values = np.array(t_values)
        N_mut = np.array(N_mut)
        mask = ~np.isnan(t_values)
        
        if np.any(mask):
            w_mean = np.average(t_values[mask], weights=N_mut[mask])
            HRD_means.append(w_mean)
        
        if N_mut_all is None:
            N_mut_all = N_mut

        # 2. Store intermediate values
        c_avg_val.append(cavg)
        for mc in range(3):
            stats['pi_2_SBS1'][mc].append(pi2_1[mc])
            stats['pi_2_SBS3'][mc].append(pi2_3[mc])
            stats['pi_1_SBS1'][mc].append(pi1_1[mc])
            stats['pi_1_SBS3'][mc].append(pi1_3[mc])
            stats['Nt_SBS1'][mc].append(nt1[mc])
            stats['Nt_SBS3'][mc].append(nt3[mc])
            stats['c_val'][mc].append(c_list[mc])

    # Convert to arrays for statistical aggregation
    HRD_means = np.array(HRD_means)
    if len(HRD_means) == 0:
        return [None]*19 # Return empty placeholders if no data found

    # Final Timing Aggregates
    hrd_time = np.nanmean(HRD_means)
    p_ci = np.nanpercentile(HRD_means, [2.5, 97.5])
    p_iqr = np.nanpercentile(HRD_means, [25, 75])

    # Helper for 95% CI semi-amplitude (standard for error bars)
    def get_err(data):
        d = np.array(data)
        d = d[~np.isnan(d)]
        if len(d) == 0: return np.nan
        bounds = np.percentile(d, [2.5, 97.5])
        return (bounds[1] - bounds[0]) / 2

    # Map intermediates to means and error bars
    results = {}
    for key in ['pi_2_SBS1', 'pi_2_SBS3', 'pi_1_SBS1', 'pi_1_SBS3']:
        results[f"{key}_mean"] = [np.nanmean(stats[key][mc]) for mc in range(3)]
        results[f"{key}_err"] = [get_err(stats[key][mc]) for mc in range(3)]

    Nt_SBS1_mean = [np.nanmean(stats['Nt_SBS1'][mc]) for mc in range(3)]
    Nt_SBS3_mean = [np.nanmean(stats['Nt_SBS3'][mc]) for mc in range(3)]
    c_val_mean = [0 if mc == 1 else np.nanmean(stats['c_val'][mc]) for mc in range(3)]
    c_avg = np.nanmean(c_avg_val)

    return (
        N_mut_all, HRD_means, hrd_time, 
        p_ci[1] - hrd_time, hrd_time - p_ci[0],  # CI hi, lo
        p_iqr[1] - hrd_time, hrd_time - p_iqr[0], # IQR hi, lo
        c_val_mean, c_avg, Nt_SBS1_mean, Nt_SBS3_mean,
        results['pi_2_SBS1_mean'], results['pi_2_SBS1_err'],
        results['pi_2_SBS3_mean'], results['pi_2_SBS3_err'],
        results['pi_1_SBS1_mean'], results['pi_1_SBS1_err'],
        results['pi_1_SBS3_mean'], results['pi_1_SBS3_err']
    )

def run_hrd_wgd_timing_analysis(hrd_wgd_timing_samples, bootstraps_dir, output_csv_path=None):
    """
    Executes the complete HRD and WGD timing analysis pipeline for a cohort.

    This function processes each sample to estimate WGD and HRD-adjusted timing 
    from bootstrapped mutation data. It aggregates timing means, 95% confidence 
    intervals, and IQRs along with mutational signature proportions (SBS1, SBS3).

    Args:
        hrd_wgd_timing_samples (dict): Dictionary where keys are Sample IDs.
        bootstraps_dir (str): Path to directory containing bootstrap results.
        output_csv_path (str, optional): Path to save the final results CSV.

    Returns:
        pd.DataFrame: A table of results for all processed samples.
    """
    all_results = []
    
    # Process each sample in a single pass to optimize file I/O
    for sample_id in tqdm(hrd_wgd_timing_samples.keys(), desc="Analyzing Timing", file=sys.stdout):
        
        # 1. Standard WGD Timing Analysis
        try:
            # Calls the CpG-based WGD timing function refactored earlier
            wgd_res = calculate_wgd_time_bootstrapping(sample_id, bootstraps_dir)
            n_mut_cpg, _, w_time, w_ci_hi, w_ci_lo, w_iqr_hi, w_iqr_lo = wgd_res
        except Exception as e:
            print(f"Warning: WGD timing failed for {sample_id}: {e}")
            w_time = w_ci_hi = w_ci_lo = w_iqr_hi = w_iqr_lo = "Error"
            n_mut_cpg = []

        # 2. HRD Timing Analysis
        try:
            # Calls the SBS1/SBS3-adjusted timing function refactored earlier
            hrd_res = calculate_hrd_time_single_sample(sample_id, bootstraps_dir)
            (n_mut_all, _, h_time, h_ci_hi, h_ci_lo, h_iqr_hi, h_iqr_lo, 
             c_list, c_avg, nt1, nt3, p2_1, p2_1_e, p2_3, p2_3_e, 
             p1_1, p1_1_e, p1_3, p1_3_e) = hrd_res
        except Exception as e:
            print(f"Warning: HRD timing failed for {sample_id}: {e}")
            h_time = h_ci_hi = h_ci_lo = h_iqr_hi = h_iqr_lo = "Error"
            n_mut_all = []
            c_list = c_avg = nt1 = nt3 = p2_1 = p2_1_e = p2_3 = p2_3_e = p1_1 = p1_1_e = p1_3 = p1_3_e = "Error"

        # 3. Compile Data Row
        all_results.append({
            "ID": sample_id,
            "HRDTime": h_time,
            "HRDTime_ci_hi": h_ci_hi,
            "HRDTime_ci_lo": h_ci_lo,
            "HRDTime_ci_IQR_hi": h_iqr_hi,
            "HRDTime_ci_IQR_lo": h_iqr_lo,
            "WGDTime": w_time,
            "WGDTime_ci_hi": w_ci_hi,
            "WGDTime_ci_lo": w_ci_lo,
            "WGDTime_ci_IQR_hi": w_iqr_hi,
            "WGDTime_ci_IQR_lo": w_iqr_lo,
            "pi2SBS1": p2_1,
            "pi2SBS1_ci": p2_1_e,
            "pi2SBS3": p2_3,
            "pi2SBS3_ci": p2_3_e,
            "pi1SBS1": p1_1,
            "pi1SBS1_ci": p1_1_e,
            "pi1SBS3": p1_3,
            "pi1SBS3_ci": p1_3_e,
            "c": c_list,
            "c_avg": c_avg,
            "Nt_SBS1": nt1,
            "Nt_SBS3": nt3,
            "N_mut_CpG": n_mut_cpg.tolist() if hasattr(n_mut_cpg, 'tolist') else n_mut_cpg,
            "N_mut_all": n_mut_all.tolist() if hasattr(n_mut_all, 'tolist') else n_mut_all
        })

    # Convert to DataFrame
    df_results = pd.DataFrame(all_results)

    # Save to CSV
    if output_csv_path:
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        df_results.to_csv(output_csv_path, index=False)

    return df_results
