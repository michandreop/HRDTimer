import os
import sys
import logging
import argparse
from contextlib import contextmanager, redirect_stdout

# 1. Add the parent directory to sys.path to find the hrdtimer package
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

try:
    from hrdtimer import preProcess
    from hrdtimer import runSignatureAnalysis as runSigAnalysis
except ImportError as e:
    print(f"Error: Could not find 'hrdtimer' modules in the parent directory: {e}")
    sys.exit(1)

def safe_print(message, outfile):
    """Prints to console and appends to the .out file safely"""
    print(message)
    try:
        with open(outfile, "a", encoding="utf-8") as f:
            f.write(message + "\n")
    except Exception as e:
        # Fallback if the file system is temporarily locked
        logging.error(f"Failed to write to .out file: {e}")

@contextmanager
def suppress_stdout():
    """Temporarily diverts stdout to devnull to hide tool-specific noise."""
    with open(os.devnull, 'w') as devnull:
        with redirect_stdout(devnull):
            yield

def find_all_mut_and_timing_dirs(base_path):
    """Recursively finds all 'all_mut' and 'timing' directories."""
    target_dirs = []
    for root, dirs, _ in os.walk(base_path):
        for d in dirs:
            if d in ['all_mut', 'timing']:
                target_dirs.append(os.path.join(root, d))
    return target_dirs

def main():
    parser = argparse.ArgumentParser(description="HRDTimer Preprocessing Pipeline")
    
    parser.add_argument("--input", type=str, 
                        default="../raw_input_data/PCAWG_MutationTimeR_output")
    parser.add_argument("--output", type=str, 
                        default="../HRDTimer_PreProcess/PCAWG_v3")
    parser.add_argument("--metadata", type=str, 
                        default="../data/metadata/pan_metadata_v5.csv")
    parser.add_argument("--log", type=str, 
                        default="output/preprocess.log")

    args = parser.parse_args()

    # Create directories before initializing loggers
    log_dir = os.path.dirname(args.log)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)

    # Path for the summary file
    summary_out = os.path.splitext(args.log)[0] + ".out"
    
    # Initialize/Clear the .out file
    with open(summary_out, "w", encoding="utf-8") as f:
        f.write("--- Pipeline Execution Started ---\n")

    # Initialize Technical Logging (.log file)
    logging.basicConfig(
        filename=args.log,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        force=True 
    )

    safe_print(f"Summary (.out) file: {os.path.abspath(summary_out)}", summary_out)
    safe_print(f"Detailed (.log) file: {os.path.abspath(args.log)}", summary_out)
    safe_print("-" * 50, summary_out)

    # Step 1: Data Output Directory
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)
        safe_print(f"Created data output directory: {args.output}", summary_out)

    # Step 2: VCF Processing
    safe_print("Step 1: Processing VCFs ...", summary_out)
    try:
        preProcess.process_vcfs_early_late_subclonal(args.input, args.output, args.metadata, time_analysis=False)
        preProcess.process_vcfs_early_late_subclonal(args.input, args.output, args.metadata, time_analysis=True)
        logging.info("VCF processing completed successfully.")
        safe_print("VCF processing complete.", summary_out)
    except Exception as e:
        logging.error(f"VCF processing failed: {e}", exc_info=True)
        safe_print(f"CRITICAL ERROR during VCF processing: {e}", summary_out)
        sys.exit(1)

    # Step 3: Signature Analysis
    safe_print("Step 2: Identifying directories for Signature Analysis...", summary_out)
    dirs_to_process = find_all_mut_and_timing_dirs(args.output)
    safe_print(f"Found {len(dirs_to_process)} target directories.", summary_out)

    for d in dirs_to_process:
        # Build logic
        genome_build = "GRCh38" if "INFORM" in d else "GRCh37"
        # Tumor type logic
        tumor_type = 'Ovary.AdenoCA' if "Ovary" in d else 'Breast.AdenoCA'
        
        msg = f"--> Running: {os.path.basename(os.path.dirname(d))}/{os.path.basename(d)} | Type: {tumor_type}"
        safe_print(msg, summary_out)
        
        try:
            runSigAnalysis.run_signature_analysis(
                d,
                genome_build=genome_build,
                tumor_type=tumor_type
            )
            logging.info(f"Successfully analyzed: {d}")
        except Exception as e:
            logging.error(f"Signature Analysis failed for {d}: {e}", exc_info=True)
            safe_print(f"FAILED: {d}. Check .log for details.", summary_out)
            # Continue to next sample despite failure
            continue

    safe_print("\n" + "="*30, summary_out)
    safe_print("FINISH: Pipeline completed successfully.", summary_out)
    safe_print("="*30, summary_out)

if __name__ == "__main__":
    main()

# -----------------------------------------------------------------------------------
# EXAMPLE RUN CASES
# -----------------------------------------------------------------------------------
#
# CASE 1: Default run (uses 'output/' folder within pipeline directory)
#   $ python preprocess.py
#   -> Creates: output/preprocess.log & output/preprocess.out
#
# CASE 2: Custom run with project-specific names
#   $ python preprocess.py --log logs/breast_cancer_run.log
#   -> Creates: logs/breast_cancer_run.log & logs/breast_cancer_run.out
#
# CASE 3: Full Custom Paths (Input, Output, and Logs in different drives/folders)
#   $ python preprocess.py \
#       --input /mnt/data/raw_vcfs \
#       --output /mnt/analysis/processed_data \
#       --metadata /mnt/data/samples.csv \
#       --log /home/user/logs/full_cohort_v1.log
#
#   -> Note: This will automatically create /home/user/logs/ and /mnt/analysis/processed_data/
#
# CASE 4: Background Execution (Recommended for Case 3)
#   $ nohup python preprocess.py --log logs/long_run.log &
#   -> You can safely close the terminal. Check progress later with: tail -f logs/long_run.out
# -----------------------------------------------------------------------------------