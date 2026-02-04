import os
import sys
import shutil
import pandas as pd
from tqdm import tqdm

# Domain-Specific Libraries
import musical
from SigProfilerAssignment import Analyzer as Analyze


def save_probabilities(H_reduced_normalized, W_reduced, output_dir):
    """
    Calculates and saves mutation attribution probabilities for each sample.

    For each sample (column) in H, this function multiplies the signature 
    definitions (W) by the sample's exposures (H), normalizes the results 
    so each row sums to 1, and saves the resulting probability matrix to a file.

    Args:
        H_reduced_normalized (pd.DataFrame): Normalized exposure matrix where 
            rows are signatures and columns are samples.
        W_reduced (pd.DataFrame): Signature matrix where rows are mutation 
            types and columns are signatures.
        output_dir (str): Path to the directory where the probability 
            text files will be saved.

    Returns:
        None: Files are written directly to the disk in tab-separated format.

    Notes:
        The output files are named using the format: 'prob_{sample_name}.txt'.
    """
    # Create a directory to store probability files if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    for col_name in H_reduced_normalized.columns:
        # col_to_multiply represents the exposure of signatures for one sample
        col_to_multiply = H_reduced_normalized[col_name]

        # Element-wise multiplication: Weighting signatures by their exposure
        result = W_reduced.mul(col_to_multiply, axis=1)
        
        # Normalize rows to get probabilities (sum of each row = 1)
        probabilities = result.div(result.sum(axis=1), axis=0)

        output_file = os.path.join(output_dir, f'prob_{col_name}.txt')
        probabilities.to_csv(output_file, sep='\t', index=True)

def calculate_exposures(SigProfiler_activities_file, mutation_file, tumor_type, output_folder):
    """
    Refits mutation data using MuSiCal and compares results against SigProfiler.

    This function performs signature refitting using MuSiCal's likelihood-based 
    approach. It also processes SigProfiler activities and MuSiCal's naive refit 
    method to serve as baseline comparisons. All models (primary and comparison) 
    are saved to the output directory.

    Args:
        activities_file (str): Path to the baseline SigProfiler activities file 
            (used for comparison).
        mutation_file (str): Path to the mutation count matrix file (SBS format).
        tumor_type (str): Tumor type used to restrict the signature catalog.
        output_folder (str): Directory where the primary results and 
            comparison baseline CSVs will be saved.

    Returns:
        tuple: A tuple containing the primary results:
            - H_reduced_normalized (pd.DataFrame): Normalized exposure matrix (H) 
              from the MuSiCal likelihood reduced model.
            - W_reduced (pd.DataFrame): The reduced signature matrix (W).
    """
    
    # Read Activities file
    exposures = pd.read_csv(SigProfiler_activities_file, delimiter='\t', index_col=0)
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

def run_signature_analysis(parent_folder, genome_build, tumor_type='Breast.AdenoCA'):
    """
    Orchestrates a comparative signature analysis pipeline (SigProfiler + MuSiCal).

    For each clonality subfolder (Early, Late, NA, Subclonal), this function:
    1. Runs SigProfilerAssignment (cosmic_fit) to obtain baseline exposures.
    2. Performs MuSiCal refitting using the SigProfiler output as a comparison.
    3. Calculates and saves mutation attribution probabilities.

    Args:
        parent_folder (str): Path to the directory containing clonality subfolders.
        genome_build (str): Genome assembly version (e.g., 'GRCh37', 'GRCh38').
        tumor_type (str): COSMIC tumor type for MuSiCal catalog restriction. 
            Defaults to 'Breast.AdenoCA'.

    Returns:
        None: All results are written to 'SignatureFitting' subdirectories.
    """
    subfolders = ['Early', 'Late', 'NA', 'Subclonal']

    for subfolder in subfolders:
        subfolder_path = os.path.join(parent_folder, subfolder)
        if not os.path.exists(subfolder_path):
            continue

        # 1. SigProfilerAssignment
        alexandrov_output = os.path.join(subfolder_path, 'Alexandrov')

        Analyze.cosmic_fit(
            subfolder_path, alexandrov_output,
            input_type="vcf", context_type="96",
            collapse_to_SBS96=False, cosmic_version=3.2, exome=False,
            genome_build=genome_build,
            make_plots=False,
            verbose=False
        )

        # 2. MuSiCal refitting
        activities_path = os.path.join(
            alexandrov_output,
            'Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
        )
        mutation_matrix_path = os.path.join(
            subfolder_path,
            'output/SBS/Input_vcffiles.SBS96.all'
        )
        musical_out_dir = os.path.join(parent_folder, 'SignatureFitting', subfolder)

        if os.path.exists(activities_path):
            H_reduced_normalized, W_reduced = calculate_exposures(
                SigProfiler_activities_file=activities_path,
                mutation_file=mutation_matrix_path,
                tumor_type=tumor_type,
                output_folder=os.path.join(musical_out_dir, 'exposures')
            )

            save_probabilities(
                H_reduced_normalized,
                W_reduced,
                output_dir=os.path.join(musical_out_dir, 'probabilities')
            )
        else:
            sys.stdout.write(
                f"Warning: SigProfiler output not found for {subfolder}. "
                "Skipping MuSiCal.\n"
            )