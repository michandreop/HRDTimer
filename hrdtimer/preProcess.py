# -------------------- Imports -------------------------
import os
import sys
import shutil
import pandas as pd
from tqdm import tqdm

# Domain-Specific Libraries
import musical
from SigProfilerAssignment import Analyzer as Analyze

# Local Utilities
from . import utils_v2 as HRDTimerUtils

def G(df):
    n = len(df)
    sum_pi = 0
    sum_den = 0 
    for index, row in df.iterrows(): 
        sum_den += row['prob_SBS1'] * row['MutCN'] / (row['MajCN'] + row['MinCN'])  
        sum_pi += row['prob_SBS1']
    G = sum_pi / sum_den
    return G

def process_vcfs_early_late_subclonal(input_folder, output_folder, metadata_csv_path, time_analysis=False, verbose=False):
    """
    Batch processes VCF files to categorize mutations by organ and clonality status.

    This function iterates through VCFs in an input directory, looks up the 
    associated organ via a metadata CSV, and splits mutations into distinct 
    files based on their 'CLS' (Clonality) type: Early, Late, NA, or Subclonal. 
    It creates a structured directory tree for the output.

    Folder Structure:
        output_folder/
        └── [Organ]/
            └── [timing | all_mut]/
                ├── Early/
                ├── Late/
                ├── NA/
                └── Subclonal/

    Args:
        input_folder (str): Path to the directory containing raw VCF files.
        output_folder (str): Base path where the organ-stratified directory 
            structure will be created.
        metadata_csv_path (str): Path to a CSV file containing at least 
            'sample' and 'organ' columns for sample-to-organ mapping.
        time_analysis (bool): If True, filters for high-confidence timing 
            regions (e.g., MajCN=2) and saves to a 'timing' subfolder. 
            If False, saves to 'all_mut'. Defaults to False.
        verbose (bool): If True, prints messages regarding directory cleanup 
            and processing status. Defaults to False.

    Returns:
        None: Processed VCFs are written to disk. Empty organ directories 
            created during initialization are removed at the end of execution.

    Note:
        The function expects `HRDTimerUtils` to be available in the namespace 
        for VCF-to-DataFrame conversions and processing.
    """

    metadata_df = pd.read_csv(metadata_csv_path)
    organ_lookup = metadata_df.set_index('sample')['organ'].to_dict()

    # Predefine output folders for each organ and CLS type
    organ_folders = {
        organ: {
            "timing" if time_analysis else "all_mut": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut"),
            "Early": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "Early"),
            "Late": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "Late"),
            "NA": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "NA"),
            "Subclonal": os.path.join(output_folder, organ, "timing" if time_analysis else "all_mut", "Subclonal")
        }
        for organ in metadata_df['organ'].unique()
    }

    # Create all required folders
    for folders in organ_folders.values():
        for folder in folders.values():
            os.makedirs(folder, exist_ok=True)

    written_folders = set()

    # Process VCF files in the input folder
    for filename in tqdm([f for f in os.listdir(input_folder) if f.endswith(".vcf")], desc="Processing VCFs", file=sys.stdout):
        vcf_path = os.path.join(input_folder, filename)
        aliquot_id = filename.split(".")[0]
        organ = organ_lookup.get(aliquot_id)

        if not organ:
            continue

        sample = HRDTimerUtils.vcf_to_dataframe(vcf_path)
        if time_analysis:
            sample = HRDTimerUtils.process_vcf_dataframe(sample, time_analysis=True)
        else:
            sample = HRDTimerUtils.process_vcf_dataframe(sample)

        cls_groups = {
            "Early": sample[sample['CLS'] == 'early'],
            "Late": sample[sample['CLS'] == 'late'],
            "NA": sample[sample['CLS'] == 'NA'],
            "Subclonal": sample[sample['CLS'] == 'subclonal']
        }

        for cls, df in cls_groups.items():
            if not df.empty:
                output_vcf = os.path.join(
                    organ_folders[organ]["timing" if time_analysis else "all_mut"], cls, f"{aliquot_id}_{cls.lower()}.vcf"
                )
                HRDTimerUtils.dataframe_to_vcf(df, output_vcf)
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
 
def prepare_samples_for_timing(vcf_folder_path):
    """
    Integrates VCF mutation data with signature attribution probabilities.

    This function aggregates mutations from 'Early', 'Late', 'NA', and 'Subclonal' 
    VCFs. For each mutation, it extracts genomic info (CNV states, VAF, context) 
    and maps it to the specific signature probabilities calculated during the 
    Signature Analysis phase.

    Args:
        vcf_folder_path (str): Path to the directory containing clonality subfolders 
            and the 'SignatureFitting' results directory.

    Returns:
        dict[str, pd.DataFrame]: A dictionary where keys are Sample IDs and 
            values are DataFrames containing merged mutation info and 
            signature probabilities. Returns None if no data is processed.

    Note:
        Assumes that 'SignatureFitting' exists within vcf_folder_path and 
        contains probability files named 'prob_{sample_id}_{clonality}.txt'.
    """
    sample_dfs = {}
    clonality_types = ['Early', 'Late', 'NA', 'Subclonal']

    for subfolder in clonality_types:
        subfolder_path = os.path.join(vcf_folder_path, subfolder)
        if not os.path.exists(subfolder_path):
            continue

        print(f'Processing {subfolder} samples:')
        vcf_files = [f for f in os.listdir(subfolder_path) if f.endswith(".vcf")]

        for filename in tqdm(vcf_files, desc=f"Merging {subfolder}", file=sys.stdout):
            sample_id = filename.split('_')[0]
            vcf_file_path = os.path.join(subfolder_path, filename)

            if os.path.getsize(vcf_file_path) == 0:
                continue

            try:
                # Load VCF and parse INFO field
                df = pd.read_csv(vcf_file_path, sep='\t', header=None, comment='#')
                df = df[[0, 1, 3, 4, 6]].copy()
                df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'INFO']
                df['INFO'] = df['INFO'].apply(HRDTimerUtils.INFO_field_to_dict)
                
                # Vectorized extraction of INFO fields
                fields = ['context', 'MajCN', 'MinCN', 'MutCN', 'pSingle', 'pGain', 'VAF']
                for field in fields:
                    col_name = 'SBS96' if field == 'context' else field
                    df[col_name] = df['INFO'].apply(lambda x: x.get(field))

                df['Classification'] = subfolder
                df['CHROM'] = df['CHROM'].astype(str)

                # Vectorized Merge of Probabilities (Much faster than row-by-row)
                prob_path = os.path.join(vcf_folder_path, "SignatureFitting", subfolder, 
                                         "probabilities", f'prob_{sample_id}_{subfolder.lower()}.txt')
                
                if os.path.exists(prob_path):
                    prob_df = pd.read_csv(prob_path, sep='\t', index_col=0)
                    prob_df.columns = [f'prob_{c}' for c in prob_df.columns]
                    # Join on SBS96 context
                    df = df.join(prob_df, on='SBS96')

                # Group by Sample ID
                if sample_id not in sample_dfs:
                    sample_dfs[sample_id] = df
                else:
                    sample_dfs[sample_id] = pd.concat([sample_dfs[sample_id], df], ignore_index=True)

            except pd.errors.EmptyDataError:
                continue

    return sample_dfs if sample_dfs else None