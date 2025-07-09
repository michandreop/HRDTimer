# Standard Library
import os
import sys
import re
import ast
import csv
import shutil
import logging
from contextlib import redirect_stdout

# Third-Party Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from tqdm import tqdm
from tabulate import tabulate
from scipy import stats, odr
from scipy.stats import linregress
import vcfpy

# Domain-Specific Libraries
from SigProfilerAssignment import Analyzer as Analyze
import musical

logging.getLogger().setLevel(logging.WARNING)
pd.options.mode.chained_assignment = None  # Disable the warning
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.size'] = 13  

pd.set_option("display.max_columns", None)

# Create a wrapper around tqdm that sets file=sys.stdout if not provided
class TQDMStdout(tqdm):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('file', sys.stdout)
        super().__init__(*args, **kwargs)

# Monkey patch tqdm globally
import builtins
import tqdm as tqdm_module
tqdm_module.tqdm = TQDMStdout
builtins.tqdm = TQDMStdout  # Optional: if some packages use `tqdm` from builtins


# -------------------- Functions Definitions -------------------------

# --------- VCF parsing/handling -----------
def read_vcf(record):
    return {
        "Chromosome": record.CHROM,
        "Position": record.POS,
        "ID": record.ID,
        "Reference Allele": record.REF,
        "Alternate Alleles": record.ALT,
        "Quality Score": record.QUAL,
        "Info": record.INFO
    }

def dict_to_vcf_INFO(d):
    parts = []
    for k, v in d.items():
        if isinstance(v, list):
            v_str = ','.join(str(i) for i in v)
        else:
            v_str = str(v)
        parts.append(f"{k}:{v_str}")
    return ';'.join(parts)

def INFO_field_to_dict(s):
    d = {}
    for item in s.split(';'):
        if not item:
            continue
        key, val = item.split(':', 1)
        if ',' in val:
            val_list = val.split(',')
            # Try to infer types (int, float)
            try:
                val_list = [x for x in val_list]
            except:
                pass
            d[key] = val_list
        else:
            # Try to infer scalar types
            if val == 'TRUE':
                d[key] = True
            elif val == 'FALSE':
                d[key] = False
            else:
                try:
                    d[key] = int(val)
                except ValueError:
                    try:
                        d[key] = float(val)
                    except ValueError:
                        d[key] = val
    return d

def vcf_to_dataframe(vcf_file):
    data = []
    with open(vcf_file, 'r') as f:
        reader = vcfpy.Reader.from_stream(f)
        for record in reader:
            data.append(read_vcf(record))
    return pd.DataFrame(data)

def dataframe_to_vcf(df, output_vcf):
    cols = df.columns[:df.columns.get_loc('INFO') + 1]
    df = df[cols].copy()
    df['INFO'] = df['INFO'].apply(dict_to_vcf_INFO)
    with open(output_vcf, 'w') as vcf:
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tINFO\n")
    df.to_csv(output_vcf, sep="\t", mode='a', index=False, header=False)

def process_vcf_dataframe(df, time_analysis=False):
    df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'INFO']
    df['ID'] = df['ID'].apply(lambda x: x[0] if x else '.')
    df['ALT'] = df['ALT'].apply(lambda x: x[0].value)

    get_info = lambda info, key: info.get(key, '.')
    fields = ['MajCN', 'MinCN', 'MutCN', 'context', 'CLS', 'inWGDregion', 'powr']
    for f in fields:
        df[f] = df['INFO'].apply(lambda info: get_info(info, f))

    df['CLS'] = df['CLS'].str.extract(r"\[(.+?)\]")
    df = df[(df['context'] != '.') & df['CLS'].notna()]
    
    if time_analysis:
        df = df[(df['MajCN'] == 2) & (df['inWGDregion'] == 'TRUE')]

    df.rename(columns={'powr': 'pow'}, inplace=True)
    return df

# --------- HRDTimer related functions -----

def G(df):
    n = len(df)
    sum_pi = 0
    sum_den = 0 
    for index, row in df.iterrows(): 
        sum_den += row['prob_SBS1'] * row['MutCN'] / (row['MajCN'] + row['MinCN'])  
        sum_pi += row['prob_SBS1']
    G = sum_pi / sum_den
    return G

def process_vcfs_early_late(input_folder, output_folder, organ_csv_path, time_analysis=False, verbose=False):
    """
    Process all VCF files in the input folder and save the processed files in the output folder.

    Args:
    - input_folder (str): Path to the folder containing VCF files.
    - output_folder (str): Path to the folder where processed VCF files will be saved.
    - organ_csv_path (str): Path to the CSV file containing organ information.
    - time_analysis (bool): If True, create a 'timing' folder; otherwise, create 'all_mut' folder.
    """
    
    organ_df = pd.read_csv(organ_csv_path)
    organ_lookup = organ_df.set_index('sample')['organ'].to_dict()

    # Predefine output folders for each organ and CLS type
    organ_folders = {
        organ: {
            "timing" if time_analysis else "all_mut": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut"),
            "Early": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "Early"),
            "Late": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "Late"),
            "NA": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "NA")
        }
        for organ in organ_df['organ'].unique()
    }
    
    # Create all required folders
    for folders in organ_folders.values():
        for folder in folders.values():
            os.makedirs(folder, exist_ok=True)

    # Process VCF files in the input folder
    for filename in tqdm([f for f in os.listdir(input_folder) if f.endswith(".vcf")], desc="Processing VCFs"):
        vcf_path = os.path.join(input_folder, filename)
        aliquot_id = filename.split(".")[0]
        organ = organ_lookup.get(aliquot_id)

        if not organ:
            #print(f"Organ information not found for {filename}, skipping.")
            continue

        # Process VCF
        sample = vcf_to_dataframe(vcf_path)
        if time_analysis:
            sample = process_vcf_dataframe(sample, time_analysis=True)
        else:
            sample = process_vcf_dataframe(sample)

        cls_groups = {
            "Early": sample[sample['CLS'] == 'early'],
            "Late": sample[sample['CLS'] == 'late'],
            "NA": sample[sample['CLS'] == 'NA']
        }

        for cls, df in cls_groups.items():
            if not df.empty:
                output_vcf = os.path.join(organ_folders[organ]["timing" if time_analysis else "all_mut"], cls, f"{aliquot_id}_{cls.lower()}.vcf")
                dataframe_to_vcf(df, output_vcf)

def process_vcfs_early_late(input_folder, output_folder, organ_csv_path, time_analysis=False, verbose=False):
    """
    Process all VCF files in the input folder and save the processed files in the output folder.

    Args:
    - input_folder (str): Path to the folder containing VCF files.
    - output_folder (str): Path to the folder where processed VCF files will be saved.
    - organ_csv_path (str): Path to the CSV file containing organ information.
    - time_analysis (bool): If True, create a 'timing' folder; otherwise, create 'all_mut' folder.
    """

    organ_df = pd.read_csv(organ_csv_path)
    organ_lookup = organ_df.set_index('sample')['organ'].to_dict()
    created_folders = set()

    # Predefine output folders for each organ and CLS type
    organ_folders = {
        organ: {
            "timing" if time_analysis else "all_mut": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut"),
            "Early": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "Early"),
            "Late": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "Late"),
            "NA": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "NA")
        }
        for organ in organ_df['organ'].unique()
    }

    # Create all required folders
    for folders in organ_folders.values():
        for folder in folders.values():
            os.makedirs(folder, exist_ok=True)

    written_folders = set()

    # Process VCF files in the input folder
    for filename in tqdm([f for f in os.listdir(input_folder) if f.endswith(".vcf")], desc="Processing VCFs"):
        vcf_path = os.path.join(input_folder, filename)
        aliquot_id = filename.split(".")[0]
        organ = organ_lookup.get(aliquot_id)

        if not organ:
            continue

        # Process VCF
        sample = vcf_to_dataframe(vcf_path)
        if time_analysis:
            sample = process_vcf_dataframe(sample, time_analysis=True)
        else:
            sample = process_vcf_dataframe(sample)

        cls_groups = {
            "Early": sample[sample['CLS'] == 'early'],
            "Late": sample[sample['CLS'] == 'late'],
            "NA": sample[sample['CLS'] == 'NA']
        }

        for cls, df in cls_groups.items():
            if not df.empty:
                output_vcf = os.path.join(
                    organ_folders[organ]["timing" if time_analysis else "all_mut"], cls, f"{aliquot_id}_{cls.lower()}.vcf"
                )
                dataframe_to_vcf(df, output_vcf)
                written_folders.add(os.path.dirname(output_vcf))

    # Cleanup: remove any empty organ folders (and subfolders)
    for organ in organ_folders:
        base_dir = os.path.join(output_folder, organ)
        if not any(
            written_folder.startswith(os.path.join(output_folder, organ)) for written_folder in written_folders
        ):
            if os.path.exists(base_dir):
                shutil.rmtree(base_dir)
                if verbose:
                    print(f"Removed empty folder: {base_dir}")

def save_probabilities(H_reduced_normalized, W_reduced, output_dir):
    # Create a directory to store probability files if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    for col_name in H_reduced_normalized.columns:
        # Create an empty DataFrame with the dimensions of W
        result = pd.DataFrame(np.zeros((1, len(W_reduced.columns))), columns=W_reduced.columns)
        col_to_multiply = H_reduced_normalized[col_name]

        # Perform element-wise multiplication for each row in W with the selected column from H
        result = W_reduced.mul(col_to_multiply, axis=1)
        probabilities = result.div(result.sum(axis=1), axis=0)

        output_file = os.path.join(output_dir, f'prob_{col_name}.txt')
        probabilities.to_csv(output_file, sep='\t', index=True)

def calculate_exposures(activities_file, mutation_file, tumor_type, output_folder):
    """
    Process data and run the MuSiCal algorithm, then save the results as CSV files.
    
    Args:
    - activities_file (str): Path to the activities file (CSV format).
    - mutation_file (str): Path to the mutation file (CSV format).
    - tumor_type (str): Tumor type for restricting the catalog.
    - parent_folder (str): Path to the parent folder where CSV files will be saved.
    
    Returns:
    - exposures (DataFrame): The exposures dataframe from SigProfiler.
    - H1 (DataFrame): The exposures dataframe from MuSiCal naive method.
    - H2 (DataFrame): The exposures dataframe from MuSiCal likelihood method.
    - model2.H_reduced_normalized (DataFrame): The normalized exposures from MuSiCal likelihood reduced model.
    """
    
    # Read Activities file
    exposures = pd.read_csv(activities_file, delimiter='\t', index_col=0)
    colnames = exposures.index.tolist()
    exposures = exposures.T
    exposures.columns = colnames
 
    mutation_matrix = pd.read_csv(mutation_file, delimiter='\t', index_col=0)
    numeric_cols = mutation_matrix.select_dtypes(include=['number']).columns
    numeric_cols_filtered = [col for col in numeric_cols if not all(mutation_matrix[col] == 0)]

    mutation_matrix = mutation_matrix[numeric_cols_filtered]
    
    # Musical
    catalog = musical.load_catalog('COSMIC_v3p2_SBS_WGS')
    catalog.restrict_catalog(tumor_type=tumor_type, is_MMRD=False, is_PPD=False)
    W = catalog.W

    # Refit with different methods
    W = W.reindex(mutation_matrix.index)

    H1, model1 = musical.refit.refit(mutation_matrix, W, method='thresh_naive', thresh=0)
    H1.columns = colnames

    H2, model2 = musical.refit.refit(mutation_matrix, W, method='likelihood_bidirectional', thresh=0.001)
    H2.columns = colnames

    # Create the parent folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Save the DataFrames to CSV files within the parent folder
    exposures.to_csv(os.path.join(output_folder, "exposures_SigProfiler.csv"))
    H1.to_csv(os.path.join(output_folder, "exposures_MuSiCal_naive.csv"))
    H2.to_csv(os.path.join(output_folder, "exposures_MuSiCal_lik.csv"))
    model2.H_reduced_normalized.to_csv(os.path.join(output_folder, "exposures_MuSiCal_lik_reduced_norm.csv"))
    model2.W_reduced.to_csv(os.path.join(output_folder, "W_MuSiCal_lik_reduced_norm.csv"))

    return model2.H_reduced_normalized, model2.W_reduced

def run_Signature_Analysis(parent_folder, genome_build):
    """
    Perform COSMIC FIT analysis and MuSiCal analysis for each subfolder in the given parent folder.
    
    Args:
    - parent_folder (str): Path to the parent folder that contains subfolders 'Early', 'Late', 'NA'.
    """
    subfolders = ['Early', 'Late', 'NA']
    
    # Iterate through subfolders
    for subfolder in subfolders:
        subfolder_path = os.path.join(parent_folder, subfolder)
        output_path = os.path.join(parent_folder, subfolder + '/Alexandrov/')
        
        Analyze.cosmic_fit(subfolder_path, output_path,
                                   input_type="vcf", context_type="96",
                                   collapse_to_SBS96=False, cosmic_version=3.2, exome=False,
                                   genome_build=genome_build,
                                   exclude_signature_subgroups=None, export_probabilities=False,
                                   export_probabilities_per_mutation=False, make_plots=False,
                                   sample_reconstruction_plots=False, verbose=False)

        # Process data and run musical for the subfolder
        H_reduced_normalized, W_reduced = calculate_exposures(
            os.path.join(output_path, 'Assignment_Solution/Activities/Assignment_Solution_Activities.txt'),
            os.path.join(parent_folder, subfolder + '/output/SBS/Input_vcffiles.SBS96.all'),
            tumor_type='Breast.AdenoCA',
            output_folder=os.path.join(parent_folder, 'SignatureFitting', subfolder, 'exposures')
        )

        # Apply save_probabilities function
        save_probabilities(H_reduced_normalized, W_reduced, 
                           output_dir=os.path.join(parent_folder, 'SignatureFitting', subfolder, 'probabilities'))
        
def prepare_samples_for_timing(vcf_folder_path):
    """
    Processes SNV files, extracts SBS96 context, and appends relevant mutation information for each sample.

    Args:
        vcf_folder_path (str): Path to folder containing SNV files.

    Returns:
        dict: Dictionary of DataFrames with sample IDs as keys, containing SBS96 info and probabilities.
    """
        
    def extract_info(info_str, key):
        return info_str.get(key, None)

    # Initialize result dictionary
    sample_dfs = {}

    # Iterate over subfolders
    for subfolder in ['Early', 'Late', 'NA']:
        subfolder_path = os.path.join(vcf_folder_path, subfolder)
        if not os.path.exists(subfolder_path): continue

        print(f'Processing {subfolder} samples:')

        # Process VCF files
        for filename in tqdm([f for f in os.listdir(subfolder_path) if f.endswith(".vcf")], desc="Processing Files"):
            sample_id = filename.split('_')[0]
            vcf_file_path = os.path.join(subfolder_path, filename)

            # Skip empty VCF files
            if os.path.getsize(vcf_file_path) == 0: continue

            # Read VCF file and extract necessary info
            try:
                vcf_df = pd.read_csv(vcf_file_path, sep='\t', header=None, comment='#')
            except pd.errors.EmptyDataError:
                continue
            
            extracted_df = vcf_df[[0, 1, 3, 4, 6]].copy()  # Extract relevant columns
            extracted_df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'INFO']
            extracted_df['INFO'] = extracted_df['INFO'].apply(INFO_field_to_dict)
            extracted_df['CHROM'] = extracted_df['CHROM'].astype(str)

            # Apply extraction functions to extract mutation data
            extracted_df['SBS96'] = extracted_df['INFO'].apply(lambda x: extract_info(x, 'context'))
            extracted_df['MajCN'] = extracted_df['INFO'].apply(lambda x: extract_info(x, 'MajCN'))
            extracted_df['MinCN'] = extracted_df['INFO'].apply(lambda x: extract_info(x, 'MinCN'))
            extracted_df['MutCN'] = extracted_df['INFO'].apply(lambda x: extract_info(x, 'MutCN'))
            extracted_df['pSingle'] = extracted_df['INFO'].apply(lambda x: extract_info(x, 'pSingle'))
            extracted_df['pGain'] = extracted_df['INFO'].apply(lambda x: extract_info(x, 'pGain'))
            extracted_df['VAF'] = extracted_df['INFO'].apply(lambda x: extract_info(x, 'VAF'))

            # Add classification column based on subfolder
            extracted_df['Classification'] = subfolder

            # Load probabilities and merge with the extracted data
            probabilities_file = os.path.join(vcf_folder_path, "SignatureFitting", subfolder, "probabilities", f'prob_{sample_id}_{subfolder.lower()}.txt')
            if os.path.exists(probabilities_file):
                probabilities_df = pd.read_csv(probabilities_file, sep='\t', index_col=0)
                for index, row in extracted_df.iterrows():
                    sbs96_value = row['SBS96']
                    if sbs96_value in probabilities_df.index:
                        for col in probabilities_df.columns:
                            extracted_df.at[index, f'prob_{col}'] = probabilities_df.at[sbs96_value, col]

            # Append to the sample's DataFrame or initialize it if not present
            if sample_id not in sample_dfs:
                sample_dfs[sample_id] = extracted_df
            else:
                sample_dfs[sample_id] = pd.concat([sample_dfs[sample_id], extracted_df], ignore_index=True)

    # Return the result dictionary if data is found, otherwise None
    return sample_dfs if sample_dfs else None

def calculate_WGDtime_prob_bootstrapping(sample_df, num_bootstrap=200, sample_fraction=1):
    # Keep [C>T]pG mutations
    sample_df = sample_df[sample_df['SBS96'].isin(['A[C>T]G', 'C[C>T]G', 'T[C>T]G', 'G[C>T]G'])]

    N_mut_CpG_all = np.array([])
    for min_cn in range(3):
        filtered_df = sample_df[sample_df['MinCN'] == min_cn]
        N_mut_CpG_all = np.append(N_mut_CpG_all, filtered_df.shape[0])

    weighted_means = np.array([])
    # Perform bootstrap sampling
    for _ in range(num_bootstrap):
        # Sample data with replacement
        bootstrap_sample = sample_df.sample(frac=sample_fraction, replace=True)
        t_values_bootstrap = np.array([])
        N_mut_CpG = np.array([])
        
        for min_cn in range(3):
            filtered_df = bootstrap_sample[bootstrap_sample['MinCN'] == min_cn]

            N_mut_CpG = np.append(N_mut_CpG, filtered_df.shape[0])

            sum_num = 0
            sum_pi = 0
            for index, row in filtered_df.iterrows():
                sum_num += row['prob_SBS1'] * row['pSingle'] / (row['pSingle'] + row['pGain'])
                sum_pi += row['prob_SBS1']
            pi_1 = sum_num / sum_pi if sum_pi != 0 else 0

            # Calculate pi_2
            sum_num = 0
            sum_pi = 0
            for index, row in filtered_df.iterrows():
                sum_num += row['prob_SBS1'] * row['pGain'] / (row['pSingle'] + row['pGain'])
                sum_pi += row['prob_SBS1']
            pi_2 = sum_num / sum_pi if sum_pi != 0 else 0

            # Calculate t_value based on MinCN
            if pi_2 == 0 and pi_1 != 0:
                t_value = 0
            elif pi_1 + 2 * pi_2 == 0:
                t_value = np.nan
            else:
                t_value = (3 * pi_2) / (pi_1 + 2 * pi_2) if min_cn == 1 else (2 * pi_2) / (pi_1 + 2 * pi_2)

            t_values_bootstrap = np.append(t_values_bootstrap,t_value)
            
        nan_indices = np.isnan(t_values_bootstrap)
        
        if np.sum(nan_indices) == 1:
            non_nan_means = t_values_bootstrap[~nan_indices]
            non_nan_weights = N_mut_CpG[~nan_indices]
            weighted_mean = np.sum(non_nan_means * non_nan_weights) / np.sum(non_nan_weights)
        elif np.sum(nan_indices) == 2:
            non_nan_means = t_values_bootstrap[~nan_indices]
            non_nan_weights = N_mut_CpG[~nan_indices]
            weighted_mean = non_nan_means[0]
        else:
            weighted_mean = np.sum(t_values_bootstrap * N_mut_CpG) / np.sum(N_mut_CpG)

        weighted_means = np.append(weighted_means,weighted_mean)

    WGD_time = np.mean(weighted_means)

    # Calculate the 95th percentile CI
    lower_bound = np.percentile(weighted_means, 2.5)
    upper_bound = np.percentile(weighted_means, 97.5)
    WDG_time_CI = (upper_bound - lower_bound)/2

    return N_mut_CpG_all, weighted_means, WGD_time, WDG_time_CI

def calculate_HRD_time(sample_df):
    # Initialize result variables
    t_values_dict = {0: None, 1: None, 2: None}
    N_mut_dict = {0: None, 1: None, 2: None}
    pi_2_SBS1_val = {0: None, 1:None, 2:None}
    pi_2_SBS3_val = {0: None, 1:None, 2:None}
    pi_1_SBS1_val = {0: None, 1:None, 2:None}
    pi_1_SBS3_val = {0: None, 1:None, 2:None}
    c_dict = {0: None, 1: None, 2: None}
    c_avg_values = {0: None, 2: None}
    Nt_SBS1_val = {0: None, 1:None, 2:None}
    Nt_SBS3_val = {0: None, 1:None, 2:None}
    c_avg = 0

    # Define the order of MinCN values
    min_cn_order = [0, 2, 1]

    # Loop through MinCN values in the specified order
    for min_cn in min_cn_order:
        filtered_df = sample_df[sample_df['MinCN'] == min_cn]
        N_mut_dict[min_cn] = filtered_df.shape[0]

        # Calculate pi_1 and pi_2 for SBS1
        sum_num_SBS1 = 0
        sum_pi_SBS1 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS1 += row['prob_SBS1'] * row['pSingle'] / (row['pSingle'] + row['pGain'])
            sum_pi_SBS1 += row['prob_SBS1']
        pi_1_SBS1 = sum_num_SBS1 / sum_pi_SBS1 if sum_pi_SBS1 != 0 else np.nan

        sum_num_SBS1 = 0
        sum_pi_SBS1 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS1 += row['prob_SBS1'] * row['pGain'] / (row['pSingle'] + row['pGain'])
            sum_pi_SBS1 += row['prob_SBS1']
        pi_2_SBS1 = sum_num_SBS1 / sum_pi_SBS1 if sum_pi_SBS1 != 0 else np.nan

        # Calculate pi_1 and pi_2 for SBS3
        sum_num_SBS3 = 0
        sum_pi_SBS3 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS3 += row['prob_SBS3'] * row['pSingle'] / (row['pSingle'] + row['pGain'])
            sum_pi_SBS3 += row['prob_SBS3']
        pi_1_SBS3 = sum_num_SBS3 / sum_pi_SBS3 if sum_pi_SBS3 != 0 else np.nan

        sum_num_SBS3 = 0
        sum_pi_SBS3 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS3 += row['prob_SBS3'] * row['pGain'] / (row['pSingle'] + row['pGain'])
            sum_pi_SBS3 += row['prob_SBS3']
        pi_2_SBS3 = sum_num_SBS3 / sum_pi_SBS3 if sum_pi_SBS3 != 0 else np.nan

        Nt_SBS1 = np.sum(filtered_df['prob_SBS1'].tolist())
        Nt_SBS3 = np.sum(filtered_df['prob_SBS3'].tolist())

        # Adjust pi_2 for SBS1 using the new formula
        pi_2_SBS1_prime = pi_2_SBS1 - (pi_2_SBS3 / pi_1_SBS3) * pi_1_SBS1 if pi_1_SBS3 != 0 else np.nan

        # Calculate t_value based on MinCN using the adjusted pi_2_SBS1_prime
        if (pi_1_SBS1 + 2 * pi_2_SBS1_prime == 0) or np.isnan(pi_2_SBS1_prime):
            t_value = np.nan
        else:
            t_value = (3 * pi_2_SBS1_prime) / (pi_1_SBS1 + 2 * pi_2_SBS1) if min_cn == 1 else (2 * pi_2_SBS1_prime) / (pi_1_SBS1 + 2 * pi_2_SBS1)

        t_values_dict[min_cn] = t_value

        if min_cn == 0 or min_cn == 2:
            c_dict[min_cn] = (pi_1_SBS1 * Nt_SBS1) / (pi_1_SBS3 * Nt_SBS3)
            c_avg_values[min_cn] = c_dict[min_cn]

        # Calculate weighted average of c for min_cn 0 and 2
        if min_cn == 2:
            c_avg_numerator = 0
            c_avg_denominator = 0
            if c_avg_values[0] is not None:
                c_avg_numerator += c_avg_values[0] * N_mut_dict[0]
                c_avg_denominator += N_mut_dict[0]
            if c_avg_values[2] is not None:
                c_avg_numerator += c_avg_values[2] * N_mut_dict[2]
                c_avg_denominator += N_mut_dict[2]
            if c_avg_denominator != 0:
                c_avg = c_avg_numerator / c_avg_denominator

        if min_cn == 1:
            # Recalculate pi_2_SBS1_prime using the weighted average c0
            pi_2_SBS1_prime = pi_2_SBS1 - pi_2_SBS3 * c_avg * (Nt_SBS3 / Nt_SBS1) if pi_1_SBS3 != 0 else np.nan
            # Recalculate t_value for min_cn 1 using the adjusted pi_2_SBS1_prime
            if (pi_1_SBS1 + 2 * pi_2_SBS1_prime == 0) or np.isnan(pi_2_SBS1_prime):
                t_value = np.nan
            else:
                t_value = (3 * pi_2_SBS1_prime) / (pi_1_SBS1 + 2 * pi_2_SBS1)
            t_values_dict[min_cn] = t_value

        pi_2_SBS1_val[min_cn] = pi_2_SBS1
        pi_2_SBS3_val[min_cn] = pi_2_SBS3
        pi_1_SBS1_val[min_cn] = pi_1_SBS1
        pi_1_SBS3_val[min_cn] = pi_1_SBS3
        Nt_SBS1_val[min_cn] = Nt_SBS1
        Nt_SBS3_val[min_cn] = Nt_SBS3


    # Arrange the results in the order 0, 1, 2
    t_values = [t_values_dict[0], t_values_dict[1], t_values_dict[2]]
    N_mut = [N_mut_dict[0], N_mut_dict[1], N_mut_dict[2]]
    c = [c_dict[0], c_dict[1], c_dict[2]]
    pi_2_SBS1_values = [pi_2_SBS1_val[0], pi_2_SBS1_val[1], pi_2_SBS1_val[2]]
    pi_2_SBS3_values = [pi_2_SBS3_val[0], pi_2_SBS3_val[1], pi_2_SBS3_val[2]]
    pi_1_SBS1_values = [pi_1_SBS1_val[0], pi_1_SBS1_val[1], pi_1_SBS1_val[2]]
    pi_1_SBS3_values = [pi_1_SBS3_val[0], pi_1_SBS3_val[1], pi_1_SBS3_val[2]]
    Nt_SBS1_values = [Nt_SBS1_val[0], Nt_SBS1_val[1], Nt_SBS1_val[2]]
    Nt_SBS3_values = [Nt_SBS3_val[0], Nt_SBS3_val[1], Nt_SBS3_val[2]]

    return t_values, N_mut, c_avg, c, pi_2_SBS1_values, pi_2_SBS3_values, Nt_SBS1_values, Nt_SBS3_values, pi_1_SBS1_values, pi_1_SBS3_values

def calculate_HRDtime_prob_bootstrapping(sample_df, num_bootstrap=200, sample_fraction=1):
    # Filter the DataFrame based on the 'SBS96' column values
    # sample_df = sample_df[sample_df['SBS96'].isin(['A[C>T]G', 'C[C>T]G', 'T[C>T]G', 'G[C>T]G'])]

    HRD_means = np.array([])
    pi_2_SBS1 = {min_cn: [] for min_cn in range(3)}
    pi_2_SBS3 = {min_cn: [] for min_cn in range(3)}
    pi_1_SBS1 = {min_cn: [] for min_cn in range(3)}
    pi_1_SBS3 = {min_cn: [] for min_cn in range(3)}
    c_val = {min_cn: [] for min_cn in range(3)}
    c_avg_val = np.array([])
    Nt_SBS1 = {min_cn: [] for min_cn in range(3)}
    Nt_SBS3 = {min_cn: [] for min_cn in range(3)}

    N_mut_all = np.array([])
    for min_cn in range(3):
        filtered_df = sample_df[sample_df['MinCN'] == min_cn]
        N_mut_all = np.append(N_mut_all, filtered_df.shape[0])

    # Perform bootstrap sampling
    for _ in range(num_bootstrap):
        # Sample data with replacement
        bootstrap_sample = sample_df.sample(frac=sample_fraction, replace=True)

        t_values, N_mut, cavg, c, pi2SBS1, pi2SBS3, NtSBS1, NtSBS3, pi1SBS1, pi1SBS3= calculate_HRD_time(bootstrap_sample)

        t_values = np.array(t_values)
        N_mut = np.array(N_mut)
        nan_indices = np.isnan(t_values)

        if np.sum(nan_indices) == 1:
            non_nan_means = t_values[~nan_indices]
            non_nan_weights = N_mut[~nan_indices]
            weighted_mean = np.sum(non_nan_means * non_nan_weights) / np.sum(non_nan_weights)
        elif np.sum(nan_indices) == 2:
            non_nan_means = t_values[~nan_indices]
            non_nan_weights = N_mut[~nan_indices]
            weighted_mean = non_nan_means[0]
        else:
            weighted_mean = np.sum(t_values * N_mut) / np.sum(N_mut)

        HRD_means = np.append(HRD_means,weighted_mean)

        for min_cn in range(3):
            pi_2_SBS1[min_cn].append(pi2SBS1[min_cn])
            pi_2_SBS3[min_cn].append(pi2SBS3[min_cn])
            pi_1_SBS1[min_cn].append(pi1SBS1[min_cn])
            pi_1_SBS3[min_cn].append(pi1SBS3[min_cn])
            Nt_SBS1[min_cn].append(NtSBS1[min_cn])
            Nt_SBS3[min_cn].append(NtSBS3[min_cn])
            c_val[min_cn].append(c[min_cn])

        c_avg_val = np.append(c_avg_val, cavg)

    HRD_time = np.mean(HRD_means)

    # Calculate the 95th percentile CI
    lower_bound = np.percentile(HRD_means, 2.5)
    upper_bound = np.percentile(HRD_means, 97.5)
    HRD_time_CI = (upper_bound - lower_bound)/2
            
    # Assuming gauss_means and gauss_4stds are lists of lists and initialized elsewhere
    pi_2_SBS1_mean = [[] for _ in range(3)]
    pi_2_SBS3_mean = [[] for _ in range(3)]
    pi_1_SBS1_mean = [[] for _ in range(3)]
    pi_1_SBS3_mean = [[] for _ in range(3)]
    pi_2_SBS1_err = [[] for _ in range(3)]
    pi_2_SBS3_err = [[] for _ in range(3)]
    pi_1_SBS1_err = [[] for _ in range(3)]
    pi_1_SBS3_err = [[] for _ in range(3)]
    Nt_SBS1_mean = [[] for _ in range(3)]
    Nt_SBS3_mean = [[] for _ in range(3)]
    c_val_mean = [[] for _ in range(3)]

    for i in range(3):
        if i == 1:
            c_val_mean[i] = 0
        else:
            c_val_mean[i] = np.mean(c_val[i])
        pi_2_SBS1_mean[i] = np.mean(pi_2_SBS1[i])
        pi_2_SBS1_err[i] = np.std(pi_2_SBS1[i])
        pi_2_SBS3_mean[i] = np.mean(pi_2_SBS3[i])
        pi_2_SBS3_err[i] = np.std(pi_2_SBS3[i])
        pi_1_SBS1_mean[i] = np.mean(pi_1_SBS1[i])
        pi_1_SBS1_err[i] = np.std(pi_1_SBS1[i])
        pi_1_SBS3_mean[i] = np.mean(pi_1_SBS3[i])
        pi_1_SBS3_err[i] = np.std(pi_1_SBS3[i])
        Nt_SBS1_mean[i] = np.mean(Nt_SBS1[i])
        Nt_SBS3_mean[i] = np.mean(Nt_SBS3[i])

    c_avg = np.mean(c_avg_val)

    return N_mut_all, HRD_means, HRD_time, HRD_time_CI, c_val_mean, c_avg, Nt_SBS1_mean, Nt_SBS3_mean, pi_2_SBS1_mean, pi_2_SBS1_err, pi_2_SBS3_mean, pi_2_SBS3_err, pi_1_SBS1_mean, pi_1_SBS1_err, pi_1_SBS3_mean, pi_1_SBS3_err

def calculate_WGDtime_CTpGs(sample_df, N_mut_CpGs):
        sample_df = sample_df[sample_df['SBS96'].isin(['A[C>T]G', 'C[C>T]G', 'T[C>T]G', 'G[C>T]G'])]
        t_values = []
        
        for min_cn in range(3):
            filtered_df = sample_df[sample_df['MinCN'] == min_cn]

            # Calculate pi_1
            sum_num = 0
            sum_pi = 0
            for index, row in filtered_df.iterrows():
                sum_num += row['pSingle'] / (row['pSingle'] + row['pGain'])
                sum_pi += 1
                #row['prob_SBS1']
            pi_1 = sum_num / sum_pi if sum_pi != 0 else 0

            # Calculate pi_2
            sum_num = 0
            sum_pi = 0
            for index, row in filtered_df.iterrows():
                sum_num += row['pGain'] / (row['pSingle'] + row['pGain'])
                sum_pi += 1
            pi_2 = sum_num / sum_pi if sum_pi != 0 else 0

            # Calculate t_value based on MinCN
            if pi_2 == 0 and pi_1 != 0:
                t_value = 0
            elif pi_1 + 2 * pi_2 == 0:
                t_value = np.nan
            else:
                t_value = (3 * pi_2) / (pi_1 + 2 * pi_2) if min_cn == 1 else (2 * pi_2) / (pi_1 + 2 * pi_2)

            t_values.append(t_value)
        

        # Convert t_values to a Pandas Series to use fillna()
        df_values_filled = pd.Series(t_values).fillna(0)
        
        # Ensure N_mut_CpGs is also a Pandas Series
        df_weights_filled = pd.Series(N_mut_CpGs).fillna(0)

        # Normalize the weights for each row
        normalized_weights = df_weights_filled.div(df_weights_filled.sum(), axis=0)
        
        # Calculate the weighted mean
        weighted_means = (df_values_filled * normalized_weights)

        return weighted_means.to_dict()

def calculate_WGDtime_prob_bootstrapping_CTpG(sample_df, num_bootstrap=200, sample_fraction=1):
    # Keep [C>T]pG mutations
    sample_df = sample_df[sample_df['SBS96'].isin(['A[C>T]G', 'C[C>T]G', 'T[C>T]G', 'G[C>T]G'])]

    N_mut_CpG_all = np.array([])
    for min_cn in range(3):
        filtered_df = sample_df[sample_df['MinCN'] == min_cn]
        N_mut_CpG_all = np.append(N_mut_CpG_all, filtered_df.shape[0])

    weighted_means = np.array([])
    # Perform bootstrap sampling
    for _ in range(num_bootstrap):
        # Sample data with replacement
        bootstrap_sample = sample_df.sample(frac=sample_fraction, replace=True)
        t_values_bootstrap = np.array([])
        N_mut_CpG = np.array([])
        
        for min_cn in range(3):
            filtered_df = bootstrap_sample[bootstrap_sample['MinCN'] == min_cn]

            N_mut_CpG = np.append(N_mut_CpG, filtered_df.shape[0])

            sum_num = 0
            sum_pi = 0
            for index, row in filtered_df.iterrows():
                sum_num += row['pSingle'] / (row['pSingle'] + row['pGain'])
                sum_pi += 1
            pi_1 = sum_num / sum_pi if sum_pi != 0 else 0

            # Calculate pi_2
            sum_num = 0
            sum_pi = 0
            for index, row in filtered_df.iterrows():
                sum_num += row['pGain'] / (row['pSingle'] + row['pGain'])
                sum_pi += 1
            pi_2 = sum_num / sum_pi if sum_pi != 0 else 0

            # Calculate t_value based on MinCN
            if pi_2 == 0 and pi_1 != 0:
                t_value = 0
            elif pi_1 + 2 * pi_2 == 0:
                t_value = np.nan
            else:
                t_value = (3 * pi_2) / (pi_1 + 2 * pi_2) if min_cn == 1 else (2 * pi_2) / (pi_1 + 2 * pi_2)

            t_values_bootstrap = np.append(t_values_bootstrap,t_value)
            
        nan_indices = np.isnan(t_values_bootstrap)
        
        if np.sum(nan_indices) == 1:
            non_nan_means = t_values_bootstrap[~nan_indices]
            non_nan_weights = N_mut_CpG[~nan_indices]
            weighted_mean = np.sum(non_nan_means * non_nan_weights) / np.sum(non_nan_weights)
        elif np.sum(nan_indices) == 2:
            non_nan_means = t_values_bootstrap[~nan_indices]
            non_nan_weights = N_mut_CpG[~nan_indices]
            weighted_mean = non_nan_means[0]
        else:
            weighted_mean = np.sum(t_values_bootstrap * N_mut_CpG) / np.sum(N_mut_CpG)

        weighted_means = np.append(weighted_means,weighted_mean)

    WGD_time = np.mean(weighted_means)

    # Calculate the 95th percentile CI
    lower_bound = np.percentile(weighted_means, 2.5)
    upper_bound = np.percentile(weighted_means, 97.5)
    WDG_time_CI = (upper_bound - lower_bound)/2

    return N_mut_CpG_all, weighted_means, WGD_time, WDG_time_CI

def add_signature_probabilities(samples_dict):
    # Load and restrict the COSMIC signature catalog
    catalog = musical.load_catalog('COSMIC_v3p2_SBS_WGS')
    catalog.restrict_catalog(tumor_type='Breast.AdenoCA', is_MMRD=False, is_PPD=False)
    
    mutation_types = catalog.W.index.tolist()
    W = catalog.W.reindex(mutation_types)

    updated_samples_dict = {}

    for sample_id, df in samples_dict.items():
        merged = []

        for label in ['Early', 'Late', 'NA']:
            group = df[df['Classification'] == label].copy()
            if group.empty:
                continue

            # Count actual SBS96 mutation types directly
            counts = group['SBS96'].value_counts().reindex(mutation_types, fill_value=0)
            count_df = pd.DataFrame({sample_id: counts})
            count_df.index.name = 'Type'

            # Refit exposures to the actual mutation profile
            exposures, _ = musical.refit.refit(count_df, W, method='likelihood_bidirectional', thresh=0.001)
            exposures_norm = (exposures / exposures.sum()).iloc[:, 0]  # Normalize exposures

            # Compute signature probabilities per mutation type
            W_weighted = W.mul(exposures_norm, axis=1)
            W_probs = W_weighted.div(W_weighted.sum(axis=1), axis=0).fillna(0).reset_index()

            # Rename signature columns to have 'prob_' prefix
            W_probs = W_probs.rename(columns={col: f"prob_{col}" for col in W.columns})

            # Merge back into the group-level mutation data
            group = group.merge(W_probs, left_on='SBS96', right_on='Type', how='left')
            group.drop(columns=['Type'], inplace=True)

            merged.append(group)

        if merged:
            updated_samples_dict[sample_id] = pd.concat(merged, ignore_index=True)
        else:
            updated_samples_dict[sample_id] = df.copy()  # fallback if nothing was merged

    return updated_samples_dict

def get_signature_exposures(samples_dict):
    # Load and restrict the COSMIC signature catalog
    catalog = musical.load_catalog('COSMIC_v3p2_SBS_WGS')
    catalog.restrict_catalog(tumor_type='Breast.AdenoCA', is_MMRD=False, is_PPD=False)
    
    mutation_types = catalog.W.index.tolist()
    W = catalog.W.reindex(mutation_types)

    # Dictionary to accumulate summed exposures per sample
    summed_exposures = {}

    for sample_id, df in samples_dict.items():
        exposures_sum = None

        for classification in ['Early', 'Late', 'NA']:
            group = df[df['Classification'] == classification]
            if group.empty:
                continue

            counts = group['SBS96'].value_counts().reindex(mutation_types, fill_value=0)
            count_df = pd.DataFrame({sample_id: counts})
            count_df.index.name = 'Type'

            exposures, _ = musical.refit.refit(count_df, W, method='likelihood_bidirectional', thresh=0.001)
            exposures_raw = exposures.iloc[:, 0]

            if exposures_sum is None:
                exposures_sum = exposures_raw
            else:
                exposures_sum = exposures_sum.add(exposures_raw, fill_value=0)

        if exposures_sum is None:
            # No data for any classification, fill with zeros
            exposures_sum = pd.Series(0, index=mutation_types, name=sample_id)

        summed_exposures[sample_id] = exposures_sum

    exposures_matrix = pd.DataFrame(summed_exposures)
    # Keep only signatures with at least one non-zero exposure across samples
    exposures_matrix = exposures_matrix.loc[(exposures_matrix != 0).any(axis=1)]

    return exposures_matrix

# -------------------------- HRD Time New bootstrapping method ---------------------------------
# ------------ Changing only signature posterior probabilities in each bootstrap ---------------

def generate_bootstraps(samples_dict, n_bootstraps, output_dir):
    """
    Perform bootstrap resampling of SBS96 mutations per sample and time classification
    ('Early', 'Late', 'NA') to estimate signature exposures using the COSMIC v3.2 catalog.

    For each bootstrap iteration:
    - For each sample and time group, mutations are resampled using multinomial sampling
      based on observed SBS96 mutation type frequencies.
    - Signature exposures are inferred via likelihood-based refitting using `musical`.
    - Posterior mutation probabilities per signature are calculated and merged with
      the original group data.
    - Annotated sample data are saved to per-sample CSV files.
    - For each group label, combined exposure vectors across all samples are saved as CSV.

    Parameters:
    ----------
    samples_dict : dict
        Dictionary mapping sample IDs to DataFrames. Each DataFrame must contain:
        - 'SBS96': mutation types (str)
        - 'Classification': one of {'Early', 'Late', 'NA'}

    n_bootstraps : int
        Number of bootstrap replicates to perform.

    output_dir : str
        Path to directory where bootstrap results will be saved. Each bootstrap replicate
        will be saved in a subdirectory named 'bootstrap_<i>'.

    Output:
    -------
    - For each bootstrap replicate:
        - <output_dir>/bootstrap_<i>/<sample_id>.csv:
            Sample-level annotated mutations with posterior signature probabilities.
        - <output_dir>/bootstrap_<i>/exposures_<label>.csv:
            Combined exposure matrix for each group label.
    """
    
    catalog = musical.load_catalog('COSMIC_v3p2_SBS_WGS')
    catalog.restrict_catalog(tumor_type='Breast.AdenoCA', is_MMRD=False, is_PPD=False)
    mutation_types = catalog.W.index.tolist()
    W = catalog.W.reindex(mutation_types)
    
    os.makedirs(output_dir, exist_ok=True)
    
    for i in tqdm(range(1, n_bootstraps + 1), desc="Bootstrapping"):
        boot_dir = os.path.join(output_dir, f'bootstrap_{i}')
        os.makedirs(boot_dir, exist_ok=True)

        group_exposure_dict = {'Early': [], 'Late': [], 'NA': []}
        
        for sample_id, df in samples_dict.items():
            merged = []

            for label in ['Early', 'Late', 'NA']:
                group = df[df['Classification'] == label].copy()
                if group.empty:
                    continue
                
                type_probs = group['SBS96'].value_counts(normalize=True).reindex(mutation_types, fill_value=0.0)
                sampled = np.random.multinomial(len(group), type_probs.values)
                sampled_types = np.repeat(type_probs.index.values, sampled)
                sampled_df = pd.DataFrame({sample_id: pd.Series(sampled_types).value_counts()
                                           .reindex(mutation_types, fill_value=0)}, index=mutation_types)
                sampled_df.index.name = 'Type'

                exposures, _ = musical.refit.refit(sampled_df, W, method='likelihood_bidirectional', thresh=0.001)
                exposures_norm = (exposures / exposures.sum()).iloc[:, 0]
                exposures_norm.name = sample_id
                exposures.name = sample_id
                group_exposure_dict[label].append(exposures)

                W_prob = W.mul(exposures_norm, axis=1)
                prob_df = W_prob.div(W_prob.sum(axis=1), axis=0).fillna(0).rename(
                    columns=lambda c: f'prob_{c}_boot').reset_index()

                group['Type'] = group['SBS96']
                group = group.merge(prob_df, on='Type', how='left')
                merged.append(group)

            if merged:
                pd.concat(merged, ignore_index=True).to_csv(os.path.join(boot_dir, f"{sample_id}.csv"), index=False)

        for label, exposure_list in group_exposure_dict.items():
            if exposure_list:
                df_exposures = pd.concat(exposure_list, axis=1).T
                df_exposures.to_csv(os.path.join(boot_dir, f"exposures_{label}.csv"))

def calculate_WGDtime_prob_bootstrapping_p(sample_id, base_dir, num_bootstrap=200):
    weighted_means = np.array([])

    for i in range(1, num_bootstrap + 1):
        N_mut_CpG_all = np.array([0, 0, 0])  # MinCN = 0, 1, 2
        # Read bootstrap file
        file_path = os.path.join(base_dir, f"bootstrap_{i}", f"{sample_id}.csv")
        if not os.path.exists(file_path):
            continue  # or raise warning/log

        sample_df = pd.read_csv(file_path)

        # Filter for [C>T]pG mutations
        sample_df = sample_df[sample_df['SBS96'].isin(['A[C>T]G', 'C[C>T]G', 'T[C>T]G', 'G[C>T]G'])]

        t_values_bootstrap = np.array([])
        N_mut_CpG = np.array([])

        for min_cn in range(3):
            filtered_df = sample_df[sample_df['MinCN'] == min_cn]
            N_mut_CpG = np.append(N_mut_CpG, filtered_df.shape[0])
            N_mut_CpG_all[min_cn] += filtered_df.shape[0]

            sum_num = sum_pi = 0
            for _, row in filtered_df.iterrows():
                sum_num += row['prob_SBS1_boot'] * row['pSingle']
                sum_pi += row['prob_SBS1_boot']
            pi_1 = sum_num / sum_pi if sum_pi else 0

            sum_num = sum_pi = 0
            for _, row in filtered_df.iterrows():
                sum_num += row['prob_SBS1_boot'] * row['pGain']
                sum_pi += row['prob_SBS1_boot']
            pi_2 = sum_num / sum_pi if sum_pi else 0

            if pi_2 == 0 and pi_1 != 0:
                t_value = 0
            elif pi_1 + 2 * pi_2 == 0:
                t_value = np.nan
            else:
                t_value = (3 * pi_2) / (pi_1 + 2 * pi_2) if min_cn == 1 else (2 * pi_2) / (pi_1 + 2 * pi_2)

            t_values_bootstrap = np.append(t_values_bootstrap, t_value)

        nan_indices = np.isnan(t_values_bootstrap)

        if np.sum(nan_indices) == 1:
            non_nan_means = t_values_bootstrap[~nan_indices]
            non_nan_weights = N_mut_CpG[~nan_indices]
            weighted_mean = np.sum(non_nan_means * non_nan_weights) / np.sum(non_nan_weights)
        elif np.sum(nan_indices) == 2:
            weighted_mean = t_values_bootstrap[~nan_indices][0]
        else:
            weighted_mean = np.sum(t_values_bootstrap * N_mut_CpG) / np.sum(N_mut_CpG)

        weighted_means = np.append(weighted_means, weighted_mean)

    WGD_time = np.nanmean(weighted_means)
    lower_bound = np.nanpercentile(weighted_means, 2.5)
    upper_bound = np.nanpercentile(weighted_means, 97.5)

    WGD_time_CI_hi = upper_bound - WGD_time
    WGD_time_CI_lo = WGD_time - lower_bound

        # Interquartile range (IQR) bounds
    q25 = np.nanpercentile(weighted_means, 25)
    q75 = np.nanpercentile(weighted_means, 75)

    WGD_time_IQR_hi = q75 - WGD_time
    WGD_time_IQR_lo = WGD_time - q25

    return N_mut_CpG_all, weighted_means, WGD_time, WGD_time_CI_hi, WGD_time_CI_lo, WGD_time_IQR_hi, WGD_time_IQR_lo

def calculate_WGDtime_prob_bootstrapping_CTpG_p(sample_id, base_dir, num_bootstrap=200):
    weighted_means = np.array([])
    N_mut_CpG_all = np.array([0, 0, 0])

    for i in range(1, num_bootstrap + 1):
        file_path = os.path.join(base_dir, f"bootstrap_{i}", f"{sample_id}.csv")
        if not os.path.exists(file_path):
            continue

        sample_df = pd.read_csv(file_path)
        sample_df = sample_df[sample_df['SBS96'].isin(['A[C>T]G', 'C[C>T]G', 'T[C>T]G', 'G[C>T]G'])]

        t_values_bootstrap = np.array([])
        N_mut_CpG = np.array([])

        for min_cn in range(3):
            filtered_df = sample_df[sample_df['MinCN'] == min_cn]
            N_mut_CpG = np.append(N_mut_CpG, filtered_df.shape[0])
            N_mut_CpG_all[min_cn] += filtered_df.shape[0]

            pi_1 = np.nanmean(filtered_df['pSingle'] / (filtered_df['pSingle'] + filtered_df['pGain'])) if not filtered_df.empty else 0
            pi_2 = np.nanmean(filtered_df['pGain'] / (filtered_df['pSingle'] + filtered_df['pGain'])) if not filtered_df.empty else 0

            if pi_2 == 0 and pi_1 != 0:
                t_value = 0
            elif pi_1 + 2 * pi_2 == 0:
                t_value = np.nan
            else:
                t_value = (3 * pi_2) / (pi_1 + 2 * pi_2) if min_cn == 1 else (2 * pi_2) / (pi_1 + 2 * pi_2)

            t_values_bootstrap = np.append(t_values_bootstrap, t_value)

        nan_indices = np.isnan(t_values_bootstrap)

        if np.sum(nan_indices) == 1:
            weighted_mean = np.sum(t_values_bootstrap[~nan_indices] * N_mut_CpG[~nan_indices]) / np.sum(N_mut_CpG[~nan_indices])
        elif np.sum(nan_indices) == 2:
            weighted_mean = t_values_bootstrap[~nan_indices][0]
        else:
            weighted_mean = np.sum(t_values_bootstrap * N_mut_CpG) / np.sum(N_mut_CpG)

        weighted_means = np.append(weighted_means, weighted_mean)

    WGD_time = np.nanmean(weighted_means)
    lower_bound = np.nanpercentile(weighted_means, 2.5)
    upper_bound = np.nanpercentile(weighted_means, 97.5)

    WGD_time_CI_hi = upper_bound - WGD_time
    WGD_time_CI_lo = WGD_time - lower_bound

    return N_mut_CpG_all, weighted_means, WGD_time, WGD_time_CI_hi, WGD_time_CI_lo

def calculate_HRD_time_p(sample_df):

    # Ensure the prob_SBS3_boot_restricted column exists
    if 'prob_SBS3_boot' not in sample_df.columns:
        sample_df['prob_SBS3_boot'] = 0
        
    # Initialize result variables
    t_values_dict = {0: None, 1: None, 2: None}
    N_mut_dict = {0: None, 1: None, 2: None}
    pi_2_SBS1_val = {0: None, 1:None, 2:None}
    pi_2_SBS3_val = {0: None, 1:None, 2:None}
    pi_1_SBS1_val = {0: None, 1:None, 2:None}
    pi_1_SBS3_val = {0: None, 1:None, 2:None}
    c_dict = {0: None, 1: None, 2: None}
    c_avg_values = {0: None, 2: None}
    Nt_SBS1_val = {0: None, 1:None, 2:None}
    Nt_SBS3_val = {0: None, 1:None, 2:None}
    c_avg = 0

    # Define the order of MinCN values
    min_cn_order = [0, 2, 1]

    # Loop through MinCN values in the specified order
    for min_cn in min_cn_order:
        filtered_df = sample_df[sample_df['MinCN'] == min_cn]
        N_mut_dict[min_cn] = filtered_df.shape[0]


        # Calculate pi_1 and pi_2 for SBS1
        sum_num_SBS1 = 0
        sum_pi_SBS1 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS1 += row['prob_SBS1_boot'] * row['pSingle'] 
            sum_pi_SBS1 += row['prob_SBS1_boot']
        pi_1_SBS1 = sum_num_SBS1 / sum_pi_SBS1 if sum_pi_SBS1 != 0 else np.nan

        sum_num_SBS1 = 0
        sum_pi_SBS1 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS1 += row['prob_SBS1_boot'] * row['pGain'] 
            sum_pi_SBS1 += row['prob_SBS1_boot']
        pi_2_SBS1 = sum_num_SBS1 / sum_pi_SBS1 if sum_pi_SBS1 != 0 else np.nan

        # Calculate pi_1 and pi_2 for SBS3
        sum_num_SBS3 = 0
        sum_pi_SBS3 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS3 += row['prob_SBS3_boot'] * row['pSingle'] 
            sum_pi_SBS3 += row['prob_SBS3_boot']
        pi_1_SBS3 = sum_num_SBS3 / sum_pi_SBS3 if sum_pi_SBS3 != 0 else np.nan

        sum_num_SBS3 = 0
        sum_pi_SBS3 = 0
        for index, row in filtered_df.iterrows():
            sum_num_SBS3 += row['prob_SBS3_boot'] * row['pGain'] 
            sum_pi_SBS3 += row['prob_SBS3_boot']
        pi_2_SBS3 = sum_num_SBS3 / sum_pi_SBS3 if sum_pi_SBS3 != 0 else np.nan

        Nt_SBS1 = np.sum(filtered_df['prob_SBS1_boot'].tolist())
        Nt_SBS3 = np.sum(filtered_df['prob_SBS3_boot'].tolist())

        # Adjust pi_2 for SBS1 using the new formula
        pi_2_SBS1_prime = pi_2_SBS1 - (pi_2_SBS3 / pi_1_SBS3) * pi_1_SBS1 if pi_1_SBS3 != 0 else np.nan

        # Calculate t_value based on MinCN using the adjusted pi_2_SBS1_prime
        if (pi_1_SBS1 + 2 * pi_2_SBS1_prime == 0) or np.isnan(pi_2_SBS1_prime):
            t_value = np.nan
        else:
            t_value = (3 * pi_2_SBS1_prime) / (pi_1_SBS1 + 2 * pi_2_SBS1) if min_cn == 1 else (2 * pi_2_SBS1_prime) / (pi_1_SBS1 + 2 * pi_2_SBS1)

        t_values_dict[min_cn] = t_value

        if min_cn == 0 or min_cn == 2:
            c_dict[min_cn] = (pi_1_SBS1 * Nt_SBS1) / (pi_1_SBS3 * Nt_SBS3)
            c_avg_values[min_cn] = c_dict[min_cn]

        # Calculate weighted average of c for min_cn 0 and 2
        if min_cn == 2:
            c_avg_numerator = 0
            c_avg_denominator = 0
            if c_avg_values[0] is not None:
                c_avg_numerator += c_avg_values[0] * N_mut_dict[0]
                c_avg_denominator += N_mut_dict[0]
            if c_avg_values[2] is not None:
                c_avg_numerator += c_avg_values[2] * N_mut_dict[2]
                c_avg_denominator += N_mut_dict[2]
            if c_avg_denominator != 0:
                c_avg = c_avg_numerator / c_avg_denominator

        if min_cn == 1:
            # Recalculate pi_2_SBS1_prime using the weighted average c0
            pi_2_SBS1_prime = pi_2_SBS1 - pi_2_SBS3 * c_avg * (Nt_SBS3 / Nt_SBS1) if pi_1_SBS3 != 0 else np.nan
            # Recalculate t_value for min_cn 1 using the adjusted pi_2_SBS1_prime
            if (pi_1_SBS1 + 2 * pi_2_SBS1_prime == 0) or np.isnan(pi_2_SBS1_prime):
                t_value = np.nan
            else:
                t_value = (3 * pi_2_SBS1_prime) / (pi_1_SBS1 + 2 * pi_2_SBS1)
            t_values_dict[min_cn] = t_value

        pi_2_SBS1_val[min_cn] = pi_2_SBS1
        pi_2_SBS3_val[min_cn] = pi_2_SBS3
        pi_1_SBS1_val[min_cn] = pi_1_SBS1
        pi_1_SBS3_val[min_cn] = pi_1_SBS3
        Nt_SBS1_val[min_cn] = Nt_SBS1
        Nt_SBS3_val[min_cn] = Nt_SBS3


    # Arrange the bootstrapped_matrices in the order 0, 1, 2
    t_values = [t_values_dict[0], t_values_dict[1], t_values_dict[2]]
    N_mut = [N_mut_dict[0], N_mut_dict[1], N_mut_dict[2]]
    c = [c_dict[0], c_dict[1], c_dict[2]]
    pi_2_SBS1_values = [pi_2_SBS1_val[0], pi_2_SBS1_val[1], pi_2_SBS1_val[2]]
    pi_2_SBS3_values = [pi_2_SBS3_val[0], pi_2_SBS3_val[1], pi_2_SBS3_val[2]]
    pi_1_SBS1_values = [pi_1_SBS1_val[0], pi_1_SBS1_val[1], pi_1_SBS1_val[2]]
    pi_1_SBS3_values = [pi_1_SBS3_val[0], pi_1_SBS3_val[1], pi_1_SBS3_val[2]]
    Nt_SBS1_values = [Nt_SBS1_val[0], Nt_SBS1_val[1], Nt_SBS1_val[2]]
    Nt_SBS3_values = [Nt_SBS3_val[0], Nt_SBS3_val[1], Nt_SBS3_val[2]]

    return t_values, N_mut, c_avg, c, pi_2_SBS1_values, pi_2_SBS3_values, Nt_SBS1_values, Nt_SBS3_values, pi_1_SBS1_values, pi_1_SBS3_values

def calculate_HRDtime_prob_bootstrapping_from_dir(sample_id, base_dir, num_bootstrap=200):
    HRD_means = np.array([])
    pi_2_SBS1 = {min_cn: [] for min_cn in range(3)}
    pi_2_SBS3 = {min_cn: [] for min_cn in range(3)}
    pi_1_SBS1 = {min_cn: [] for min_cn in range(3)}
    pi_1_SBS3 = {min_cn: [] for min_cn in range(3)}
    c_val = {min_cn: [] for min_cn in range(3)}
    c_avg_val = np.array([])
    Nt_SBS1 = {min_cn: [] for min_cn in range(3)}
    Nt_SBS3 = {min_cn: [] for min_cn in range(3)}
    N_mut_all = np.array([])

    for i in range(num_bootstrap):
        file_path = os.path.join(base_dir, f"bootstrap_{i+1}", f"{sample_id}.csv")
        if not os.path.exists(file_path):
            continue  # Skip missing files

        bootstrap_sample = pd.read_csv(file_path)

        t_values, N_mut, cavg, c, pi2SBS1, pi2SBS3, NtSBS1, NtSBS3, pi1SBS1, pi1SBS3 = calculate_HRD_time_p(bootstrap_sample)

        t_values = np.array(t_values)
        N_mut = np.array(N_mut)
        nan_indices = np.isnan(t_values)

        if np.sum(nan_indices) == 1:
            non_nan_means = t_values[~nan_indices]
            non_nan_weights = N_mut[~nan_indices]
            weighted_mean = np.sum(non_nan_means * non_nan_weights) / np.sum(non_nan_weights)
        elif np.sum(nan_indices) == 2:
            non_nan_means = t_values[~nan_indices]
            weighted_mean = non_nan_means[0]
        else:
            weighted_mean = np.sum(t_values * N_mut) / np.sum(N_mut)

        HRD_means = np.append(HRD_means, weighted_mean)
        N_mut_all = N_mut if i == 0 else N_mut_all  # only need to store once

        for min_cn in range(3):
            pi_2_SBS1[min_cn].append(pi2SBS1[min_cn])
            pi_2_SBS3[min_cn].append(pi2SBS3[min_cn])
            pi_1_SBS1[min_cn].append(pi1SBS1[min_cn])
            pi_1_SBS3[min_cn].append(pi1SBS3[min_cn])
            Nt_SBS1[min_cn].append(NtSBS1[min_cn])
            Nt_SBS3[min_cn].append(NtSBS3[min_cn])
            c_val[min_cn].append(c[min_cn])

        c_avg_val = np.append(c_avg_val, cavg)

    #HRD_time = np.mean(HRD_means)
    #HRD_time_CI = (np.percentile(HRD_means, 97.5) - np.percentile(HRD_means, 2.5)) / 2
    #HRD_time_CI_hi = np.percentile(HRD_means, 97.5) - HRD_time
    #HRD_time_CI_lo = HRD_time - np.percentile(HRD_means, 2.5)

    
    HRD_time_CI = (np.nanpercentile(HRD_means, 97.5) - np.nanpercentile(HRD_means, 2.5)) / 2

    HRD_time = np.nanmean(HRD_means)
    HRD_time_CI_hi = np.nanpercentile(HRD_means, 97.5) - HRD_time
    HRD_time_CI_lo = HRD_time - np.nanpercentile(HRD_means, 2.5)

    # Interquartile range (IQR)
    q25 = np.nanpercentile(HRD_means, 25)
    q75 = np.nanpercentile(HRD_means, 75)

    HRD_time_IQR_hi = q75 - HRD_time
    HRD_time_IQR_lo = HRD_time - q25

    # Mean and std calculations
    pi_2_SBS1_mean, pi_2_SBS1_err = [], []
    pi_2_SBS3_mean, pi_2_SBS3_err = [], []
    pi_1_SBS1_mean, pi_1_SBS1_err = [], []
    pi_1_SBS3_mean, pi_1_SBS3_err = [], []
    Nt_SBS1_mean, Nt_SBS3_mean = [], []
    c_val_mean = []

    '''
    for i in range(3):
        c_val_mean.append(0 if i == 1 else np.nanmean(c_val[i]))
        pi_2_SBS1_mean.append(np.nanmean(pi_2_SBS1[i]))
        pi_2_SBS1_err.append(np.nanstd(pi_2_SBS1[i]))
        pi_2_SBS3_mean.append(np.nanmean(pi_2_SBS3[i]))
        pi_2_SBS3_err.append(np.nanstd(pi_2_SBS3[i]))
        pi_1_SBS1_mean.append(np.nanmean(pi_1_SBS1[i]))
        pi_1_SBS1_err.append(np.nanstd(pi_1_SBS1[i]))
        pi_1_SBS3_mean.append(np.nanmean(pi_1_SBS3[i]))
        pi_1_SBS3_err.append(np.nanstd(pi_1_SBS3[i]))
        Nt_SBS1_mean.append(np.nanmean(Nt_SBS1[i]))
        Nt_SBS3_mean.append(np.nanmean(Nt_SBS3[i]))

    c_avg = np.nanmean(c_avg_val)
    '''

    for i in range(3):
        c_val_mean.append(0 if i == 1 else np.nanmean(c_val[i]))

        def ci95_percentile(data):
            data = np.array(data)
            data = data[~np.isnan(data)]
            if len(data) == 0:
                return np.nan
            lower = np.percentile(data, 2.5)
            upper = np.percentile(data, 97.5)
            return (upper - lower) / 2

        pi_2_SBS1_mean.append(np.nanmean(pi_2_SBS1[i]))
        pi_2_SBS1_err.append(ci95_percentile(pi_2_SBS1[i]))

        pi_2_SBS3_mean.append(np.nanmean(pi_2_SBS3[i]))
        pi_2_SBS3_err.append(ci95_percentile(pi_2_SBS3[i]))

        pi_1_SBS1_mean.append(np.nanmean(pi_1_SBS1[i]))
        pi_1_SBS1_err.append(ci95_percentile(pi_1_SBS1[i]))

        pi_1_SBS3_mean.append(np.nanmean(pi_1_SBS3[i]))
        pi_1_SBS3_err.append(ci95_percentile(pi_1_SBS3[i]))

        Nt_SBS1_mean.append(np.nanmean(Nt_SBS1[i]))
        Nt_SBS3_mean.append(np.nanmean(Nt_SBS3[i]))

    c_avg = np.nanmean(c_avg_val)


    return N_mut_all, HRD_means, HRD_time, HRD_time_CI_hi, HRD_time_CI_lo,  HRD_time_IQR_hi, HRD_time_IQR_lo, c_val_mean, c_avg, Nt_SBS1_mean, Nt_SBS3_mean, pi_2_SBS1_mean, pi_2_SBS1_err, pi_2_SBS3_mean, pi_2_SBS3_err, pi_1_SBS1_mean, pi_1_SBS1_err, pi_1_SBS3_mean, pi_1_SBS3_err

def run_HRD_WGD_timing_analysis(hrd_wgd_timing_samples, base_dir, output_csv_path):
    """
    Run HRD and WGD timing analysis for a set of samples and save the results to a CSV file.

    This function processes each sample in the input dictionary `hrd_wgd_timing_samples` to estimate
    the timing of whole-genome duplication (WGD) and homologous recombination deficiency (HRD) events
    based on bootstrapped mutation data. It calculates means and confidence intervals for timing estimates,
    mutational signature proportions (SBS1, SBS3), and other related metrics. The results are printed
    in tabular format and saved to a CSV file.

    Parameters:
    ----------
    hrd_wgd_timing_samples : dict
        A dictionary where keys are sample IDs to be processed.
        
    base_dir : str
        The base directory path where input data for each sample is stored.
        
    output_csv_path : str
        Path to the output CSV file where the aggregated timing results will be saved.

    Outputs:
    -------
    - Prints a summary table of HRD and WGD timing estimates with associated statistics.
    - Saves a CSV file containing the timing results and signature metrics for each sample.

    For each sample, the following values are computed and recorded:
    - HRDTime, HRDTime_ci_hi, HRDTime_ci_lo
    - WGDTime, WGDTime_ci_hi, WGDTime_ci_lo
    - WGDTime_CpG
    - pi2SBS1, pi2SBS1_ci, pi2SBS3, pi2SBS3_ci
    - pi1SBS1, pi1SBS1_ci, pi1SBS3, pi1SBS3_ci
    - c (mean coefficient), c21 (adjusted coefficient)
    - Nt_SBS1, Nt_SBS3 (mutation counts)
    - N_mut(C>TpG), N_mut_all
    """

    # Initialize containers
    WGDTime_means, WGDTime_CpGs = {}, {}
    WGDTime_error_hi, WGDTime_error_lo = {}, {}
    WGDTime_CpGs_error_hi, WGDTime_CpGs_error_lo = {}, {}
    WGDTime_error_IQR_hi, WGDTime_error_IQR_lo = {}, {}

    NmutCpG = {}

    HRDTime_means, HRDTime_error_hi, HRDTime_error_lo = {}, {}, {}
    HRDTime_error_IQR_hi, HRDTime_error_IQR_lo = {}, {}
    pi2SBS1_means, pi2SBS1_err = {}, {}
    pi2SBS3_means, pi2SBS3_err = {}, {}
    pi1SBS1_means, pi1SBS1_err = {}, {}
    pi1SBS3_means, pi1SBS3_err = {}, {}
    c, c_avg, NtSBS1, NtSBS3, Nmutall = {}, {}, {}, {}, {}

    # WGD Time estimation loop
    for sample_id in tqdm(hrd_wgd_timing_samples.keys(), desc="Processing WGD Samples"):
        N_mut_CpG, _, WGDTime, WGDTime_CI_hi, WGDTime_CI_lo,  WGD_time_IQR_hi, WGD_time_IQR_lo = calculate_WGDtime_prob_bootstrapping_p(sample_id, base_dir)
        _, _, WGDTime_CpG, WGDTime_CpG_CI_hi, WGDTime_CpG_CI_lo = calculate_WGDtime_prob_bootstrapping_CTpG_p(sample_id, base_dir)

        WGDTime_means[sample_id] = WGDTime
        WGDTime_error_hi[sample_id] = WGDTime_CI_hi
        WGDTime_error_lo[sample_id] = WGDTime_CI_lo

        WGDTime_error_IQR_hi[sample_id] = WGD_time_IQR_hi
        WGDTime_error_IQR_lo[sample_id] = WGD_time_IQR_lo

        NmutCpG[sample_id] = N_mut_CpG.tolist()
        WGDTime_CpGs[sample_id] = WGDTime_CpG
        WGDTime_CpGs_error_hi[sample_id] = WGDTime_CpG_CI_hi
        WGDTime_CpGs_error_lo[sample_id] = WGDTime_CpG_CI_lo

    # HRD Time estimation loop
    for sample_id in tqdm(hrd_wgd_timing_samples.keys(), desc="Processing HRD Samples"):
        results = calculate_HRDtime_prob_bootstrapping_from_dir(sample_id, base_dir)
        N_mut_all, _, HRD_time, HRD_time_CI_hi, HRD_time_CI_lo, HRD_time_IQR_hi, HRD_time_IQR_lo, c_val_mean, cavg, Nt_SBS1, Nt_SBS3, \
        pi_2_SBS1_mean, pi_2_SBS1_err, pi_2_SBS3_mean, pi_2_SBS3_err, \
        pi_1_SBS1_mean, pi_1_SBS1_err, pi_1_SBS3_mean, pi_1_SBS3_err = results

        HRDTime_means[sample_id] = HRD_time
        HRDTime_error_hi[sample_id] = HRD_time_CI_hi
        HRDTime_error_lo[sample_id] = HRD_time_CI_lo

        HRDTime_error_IQR_hi[sample_id] = HRD_time_IQR_hi
        HRDTime_error_IQR_lo[sample_id] = HRD_time_IQR_lo


        pi2SBS1_means[sample_id] = pi_2_SBS1_mean
        pi2SBS1_err[sample_id] = pi_2_SBS1_err
        pi2SBS3_means[sample_id] = pi_2_SBS3_mean
        pi2SBS3_err[sample_id] = pi_2_SBS3_err
        pi1SBS1_means[sample_id] = pi_1_SBS1_mean
        pi1SBS1_err[sample_id] = pi_1_SBS1_err
        pi1SBS3_means[sample_id] = pi_1_SBS3_mean
        pi1SBS3_err[sample_id] = pi_1_SBS3_err

        c[sample_id] = c_val_mean 
        c_avg[sample_id] = cavg
        NtSBS1[sample_id] = Nt_SBS1
        NtSBS3[sample_id] = Nt_SBS3
        Nmutall[sample_id] = N_mut_all.tolist()

    # Aggregate results
    results = []
    for aliquot_id in hrd_wgd_timing_samples:
        results.append([
            aliquot_id,
            HRDTime_means.get(aliquot_id, "Not available"),
            HRDTime_error_hi.get(aliquot_id, "Not available"),
            HRDTime_error_lo.get(aliquot_id, "Not available"),
            HRDTime_error_IQR_hi.get(aliquot_id, "Not available"),
            HRDTime_error_IQR_lo.get(aliquot_id, "Not available"),

            WGDTime_means.get(aliquot_id, "Not available"),
            WGDTime_error_hi.get(aliquot_id, "Not available"),
            WGDTime_error_lo.get(aliquot_id, "Not available"),
            WGDTime_error_IQR_hi.get(aliquot_id, "Not available"),
            WGDTime_error_IQR_lo.get(aliquot_id, "Not available"),

            WGDTime_CpGs.get(aliquot_id, "Not available"),
            pi2SBS1_means.get(aliquot_id, "Not available"),
            pi2SBS1_err.get(aliquot_id, "Not available"),
            pi2SBS3_means.get(aliquot_id, "Not available"),
            pi2SBS3_err.get(aliquot_id, "Not available"),
            pi1SBS1_means.get(aliquot_id, "Not available"),
            pi1SBS1_err.get(aliquot_id, "Not available"),
            pi1SBS3_means.get(aliquot_id, "Not available"),
            pi1SBS3_err.get(aliquot_id, "Not available"),
            c.get(aliquot_id, "Not available"),
            c_avg.get(aliquot_id, "Not available"),
            NtSBS1.get(aliquot_id, "Not available"),
            NtSBS3.get(aliquot_id, "Not available"),
            NmutCpG.get(aliquot_id, "Not available"),
            Nmutall.get(aliquot_id, "Not available")
        ])

    # Print result table
    #print(tabulate(results, headers=[
    #    "ID", "HRDTime", "HRDTime_ci_hi", "HRDTime_ci_lo", "HRDTime_ci_IQR_hi", "HRDTime_ci_IQR_lo",
    #    "WGDTime", "WGDTime_ci_hi", "WGDTime_ci_lo", "WGDTime_ci_IQR_hi", "WGDTime_ci_IQR_lo",
    #    "WGDTime_CpG", 'pi2SBS1', 'pi2SBS1_ci', 'pi2SBS3', 'pi2SBS3_ci',
    #    'pi1SBS1', 'pi1SBS1_ci', 'pi1SBS3', 'pi1SBS3_ci',
    #    "c", "c21", "Nt_SBS1", "Nt_SBS3", "N_mut(C>TpG)", "N_mut_all"
    #], tablefmt="grid"))

    # Save results to CSV
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    with open(output_csv_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([
            "ID", "HRDTime", "HRDTime_ci_hi", "HRDTime_ci_lo", "HRDTime_ci_IQR_hi", "HRDTime_ci_IQR_lo",
            "WGDTime", "WGDTime_ci_hi", "WGDTime_ci_lo", "WGDTime_ci_IQR_hi", "WGDTime_ci_IQR_lo",
            "WGDTime_CpG", 'pi2SBS1', 'pi2SBS1_ci', 'pi2SBS3', 'pi2SBS3_ci',
            'pi1SBS1', 'pi1SBS1_ci', 'pi1SBS3', 'pi1SBS3_ci',
            "c", "c21", "Nt_SBS1", "Nt_SBS3", "N_mut(C>TpG)", "N_mut_all"
        ])
        writer.writerows(results)

def run_HRD_WGD_timing_analysis(hrd_wgd_timing_samples, bootstraps_dir, output_csv_path=None):
    """
    Run HRD and WGD timing analysis for a set of samples and save the results to a CSV file.

    This function processes each sample in the input dictionary `hrd_wgd_timing_samples` to estimate
    the timing of whole-genome duplication (WGD) and homologous recombination deficiency (HRD) events
    based on bootstrapped mutation data. It calculates means and confidence intervals for timing estimates,
    mutational signature proportions (SBS1, SBS3), and other related metrics. The results are printed
    in tabular format and saved to a CSV file.

    Parameters:
    ----------
    hrd_wgd_timing_samples : dict
        A dictionary where keys are sample IDs to be processed.
        
    bootstraps_dir : str
        The directory path where input bootstrapped data for each sample is stored.
        
    output_csv_path : str
        Path to the output CSV file where the aggregated timing results will be saved.

    Returns:
    -------
    pandas.DataFrame
        A DataFrame containing the full set of results per sample.
    """

    # Initialize containers
    WGDTime_means, WGDTime_CpGs = {}, {}
    WGDTime_error_hi, WGDTime_error_lo = {}, {}
    WGDTime_CpGs_error_hi, WGDTime_CpGs_error_lo = {}, {}
    WGDTime_error_IQR_hi, WGDTime_error_IQR_lo = {}, {}
    NmutCpG = {}

    HRDTime_means, HRDTime_error_hi, HRDTime_error_lo = {}, {}, {}
    HRDTime_error_IQR_hi, HRDTime_error_IQR_lo = {}, {}
    pi2SBS1_means, pi2SBS1_err = {}, {}
    pi2SBS3_means, pi2SBS3_err = {}, {}
    pi1SBS1_means, pi1SBS1_err = {}, {}
    pi1SBS3_means, pi1SBS3_err = {}, {}
    c, c_avg, NtSBS1, NtSBS3, Nmutall = {}, {}, {}, {}, {}

    # WGD Time estimation loop
    for sample_id in tqdm(hrd_wgd_timing_samples.keys(), desc="Timing WGD"):
        N_mut_CpG, _, WGDTime, WGDTime_CI_hi, WGDTime_CI_lo,  WGD_time_IQR_hi, WGD_time_IQR_lo = calculate_WGDtime_prob_bootstrapping_p(sample_id, bootstraps_dir)
        _, _, WGDTime_CpG, WGDTime_CpG_CI_hi, WGDTime_CpG_CI_lo = calculate_WGDtime_prob_bootstrapping_CTpG_p(sample_id, bootstraps_dir)

        WGDTime_means[sample_id] = WGDTime
        WGDTime_error_hi[sample_id] = WGDTime_CI_hi
        WGDTime_error_lo[sample_id] = WGDTime_CI_lo
        WGDTime_error_IQR_hi[sample_id] = WGD_time_IQR_hi
        WGDTime_error_IQR_lo[sample_id] = WGD_time_IQR_lo

        NmutCpG[sample_id] = N_mut_CpG.tolist()
        WGDTime_CpGs[sample_id] = WGDTime_CpG
        WGDTime_CpGs_error_hi[sample_id] = WGDTime_CpG_CI_hi
        WGDTime_CpGs_error_lo[sample_id] = WGDTime_CpG_CI_lo

    # HRD Time estimation loop
    for sample_id in tqdm(hrd_wgd_timing_samples.keys(), desc="Timing HRD"):
        results = calculate_HRDtime_prob_bootstrapping_from_dir(sample_id, bootstraps_dir)
        N_mut_all, _, HRD_time, HRD_time_CI_hi, HRD_time_CI_lo, HRD_time_IQR_hi, HRD_time_IQR_lo, c_val_mean, cavg, Nt_SBS1, Nt_SBS3, \
        pi_2_SBS1_mean, pi_2_SBS1_err, pi_2_SBS3_mean, pi_2_SBS3_err, \
        pi_1_SBS1_mean, pi_1_SBS1_err, pi_1_SBS3_mean, pi_1_SBS3_err = results

        HRDTime_means[sample_id] = HRD_time
        HRDTime_error_hi[sample_id] = HRD_time_CI_hi
        HRDTime_error_lo[sample_id] = HRD_time_CI_lo
        HRDTime_error_IQR_hi[sample_id] = HRD_time_IQR_hi
        HRDTime_error_IQR_lo[sample_id] = HRD_time_IQR_lo

        pi2SBS1_means[sample_id] = pi_2_SBS1_mean
        pi2SBS1_err[sample_id] = pi_2_SBS1_err
        pi2SBS3_means[sample_id] = pi_2_SBS3_mean
        pi2SBS3_err[sample_id] = pi_2_SBS3_err
        pi1SBS1_means[sample_id] = pi_1_SBS1_mean
        pi1SBS1_err[sample_id] = pi_1_SBS1_err
        pi1SBS3_means[sample_id] = pi_1_SBS3_mean
        pi1SBS3_err[sample_id] = pi_1_SBS3_err

        c[sample_id] = c_val_mean 
        c_avg[sample_id] = cavg
        NtSBS1[sample_id] = Nt_SBS1
        NtSBS3[sample_id] = Nt_SBS3
        Nmutall[sample_id] = N_mut_all.tolist()

    # Aggregate results
    results = []
    for aliquot_id in hrd_wgd_timing_samples:
        results.append([
            aliquot_id,
            HRDTime_means.get(aliquot_id, "Not available"),
            HRDTime_error_hi.get(aliquot_id, "Not available"),
            HRDTime_error_lo.get(aliquot_id, "Not available"),
            HRDTime_error_IQR_hi.get(aliquot_id, "Not available"),
            HRDTime_error_IQR_lo.get(aliquot_id, "Not available"),
            WGDTime_means.get(aliquot_id, "Not available"),
            WGDTime_error_hi.get(aliquot_id, "Not available"),
            WGDTime_error_lo.get(aliquot_id, "Not available"),
            WGDTime_error_IQR_hi.get(aliquot_id, "Not available"),
            WGDTime_error_IQR_lo.get(aliquot_id, "Not available"),
            WGDTime_CpGs.get(aliquot_id, "Not available"),
            pi2SBS1_means.get(aliquot_id, "Not available"),
            pi2SBS1_err.get(aliquot_id, "Not available"),
            pi2SBS3_means.get(aliquot_id, "Not available"),
            pi2SBS3_err.get(aliquot_id, "Not available"),
            pi1SBS1_means.get(aliquot_id, "Not available"),
            pi1SBS1_err.get(aliquot_id, "Not available"),
            pi1SBS3_means.get(aliquot_id, "Not available"),
            pi1SBS3_err.get(aliquot_id, "Not available"),
            c.get(aliquot_id, "Not available"),
            c_avg.get(aliquot_id, "Not available"),
            NtSBS1.get(aliquot_id, "Not available"),
            NtSBS3.get(aliquot_id, "Not available"),
            NmutCpG.get(aliquot_id, "Not available"),
            Nmutall.get(aliquot_id, "Not available")
        ])

    # Write to CSV
    headers = [
        "ID", "HRDTime", "HRDTime_ci_hi", "HRDTime_ci_lo", "HRDTime_ci_IQR_hi", "HRDTime_ci_IQR_lo",
        "WGDTime", "WGDTime_ci_hi", "WGDTime_ci_lo", "WGDTime_ci_IQR_hi", "WGDTime_ci_IQR_lo",
        "WGDTime_CpG", 'pi2SBS1', 'pi2SBS1_ci', 'pi2SBS3', 'pi2SBS3_ci',
        'pi1SBS1', 'pi1SBS1_ci', 'pi1SBS3', 'pi1SBS3_ci',
        "c", "c21", "Nt_SBS1", "Nt_SBS3", "N_mut(C>TpG)", "N_mut_all"
    ]
    if output_csv_path:
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        with open(output_csv_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            writer.writerows(results)

    # Return a DataFrame for further use
    return pd.DataFrame(results, columns=headers)
