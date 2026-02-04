#!/usr/bin/env python3
"""
generate_slurm_jobs.py

Generate SLURM batch scripts for HRDTimer analysis.

This script:
- Finds samples in the input directory
- Filters samples based on metadata (HRDetect.isHRD & isWGD)
- Chunks them into batches
- Creates SLURM scripts ready for submission
- Handles logging directories

Example Usage:
$ python generate_slurm_jobs.py \
    --python-script /path/to/run_analysis.py \
    --input-dir /path/to/vcfs \
    --metadata /path/to/metadata.csv \
    --chunk-size 10
"""

import os
import sys
import argparse
import logging
import pandas as pd
from pathlib import Path

# ------------------
# Add parent directory to sys.path to find hrdtimer
# ------------------
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
try:
    from hrdtimer import preProcess
except ImportError as e:
    print(f"ERROR: Could not find 'hrdtimer' package: {e}")
    sys.exit(1)

# ------------------
# Helper functions
# ------------------
def chunk_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def safe_print(msg, out_file=None):
    print(msg)
    if out_file:
        try:
            with open(out_file, "a", encoding="utf-8") as f:
                f.write(msg + "\n")
        except Exception as e:
            logging.warning(f"Failed to write to {out_file}: {e}")

# ------------------
# Main
# ------------------
def main():
    parser = argparse.ArgumentParser(description="Generate SLURM jobs for HRDTimer analysis")
    parser.add_argument("--python-script", required=True, help="Path to the Python analysis script")
    parser.add_argument("--input-dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--metadata", type=str, default="../data/metadata/pan_metadata_v5.csv",
                        help="Path to metadata CSV file")
    parser.add_argument("--chunk-size", type=int, default=10, help="Number of samples per SLURM job")
    parser.add_argument("--slurm-dir", default="./slurm_jobs", help="Directory to save SLURM scripts")
    parser.add_argument("--time", default="0-06:00", help="SLURM job time")
    parser.add_argument("--mem", default="10G", help="SLURM job memory")
    parser.add_argument("--cpus", type=int, default=1, help="SLURM job CPUs")
    args = parser.parse_args()

    # ------------------
    # Prepare directories
    # ------------------
    slurm_dir = Path(args.slurm_dir)
    slurm_dir.mkdir(parents=True, exist_ok=True)

    log_dir = slurm_dir / "logs"
    log_dir.mkdir(exist_ok=True, parents=True)
    summary_out = log_dir / "slurm_job_summary.out"

    # Clear previous summary
    with open(summary_out, "w") as f:
        f.write("--- SLURM Job Generation Started ---\n")

    safe_print(f"Using input directory: {args.input_dir}", summary_out)
    safe_print(f"SLURM scripts will be created in: {slurm_dir}", summary_out)
    safe_print(f"Summary file: {summary_out}", summary_out)
    safe_print("-" * 50, summary_out)

    # ------------------
    # Load metadata and filter samples
    # ------------------
    metadata = pd.read_csv(args.metadata)

    filtered_metadata = metadata[
        (metadata['HRDetect.isHRD'] == True) &
        (metadata['isWGD'] == True)
    ]

    valid_sample_ids = filtered_metadata['sample'].tolist()
    safe_print(f"{len(valid_sample_ids)} samples passed HRD & WGD filters.", summary_out)

    # ------------------
    # Load samples from input directory
    # ------------------
    samples_to_time = preProcess.prepare_samples_for_timing(args.input_dir)

    # Only keep samples that pass the metadata filter
    sample_ids = [s for s in samples_to_time.keys() if s in valid_sample_ids]
    safe_print(f"{len(sample_ids)} samples will be processed after filtering.", summary_out)

    if len(sample_ids) == 0:
        safe_print("No samples to process after filtering. Exiting.", summary_out)
        return

    # ------------------
    # Generate SLURM scripts
    # ------------------
    for job_idx, chunk in enumerate(chunk_list(sample_ids, args.chunk_size)):
        batch_id = f"batch_{job_idx:03d}"
        job_name = f"hrdtimer_run_{batch_id}"
        slurm_file = slurm_dir / f"{job_name}.sbatch"
        sample_args = " ".join(chunk)

        slurm_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/{job_name}.out
#SBATCH --error={log_dir}/{job_name}.err
#SBATCH --time={args.time}
#SBATCH --mem={args.mem}
#SBATCH -c {args.cpus}
#SBATCH -A park_contrib
#SBATCH -p park

set -euo pipefail

source /home/mia218/miniforge3/etc/profile.d/conda.sh
conda activate python311_signatures

# Combine deterministic batch ID with SLURM job ID
JOB_ID={batch_id}_${{SLURM_JOB_ID}}

python {args.python_script} \\
  --job-id $JOB_ID \\
  --sample-ids {sample_args} \\
  --input-dir {args.input_dir}
"""
        slurm_file.write_text(slurm_text)
        safe_print(f"Created SLURM job script: {slurm_file}", summary_out)

    safe_print(f"\nGenerated {job_idx + 1} SLURM job files in {slurm_dir}", summary_out)
    safe_print("FINISHED: SLURM job generation.", summary_out)

# ------------------
# Entry point
# ------------------
if __name__ == "__main__":
    main()

# -----------------------------------------------------------------------------------
# EXAMPLE COMMANDS
# -----------------------------------------------------------------------------------
#
# Default run (creates ./slurm_jobs and logs/ subfolder)
# $ python generate_slurm_jobs.py \
#    --python-script /home/mia218/park_home_dir/HRDTimer/pipeline/run_timing_analysis.py \
#    --input-dir /home/mia218/park_home_dir/HRDTimer/HRDTimer_PreProcess/PCAWG_v3/Breast/timing
#
# Custom chunk size and SLURM folder
# $ python generate_slurm_jobs.py \
#     --python-script /home/mia218/park_home_dir/HRDTimer/pipeline/run_timing_analysis.py \
#     --input-dir /home/mia218/park_home_dir/HRDTimer/HRDTimer_PreProcess/PCAWG_v3/Breast/timing \
#     --chunk-size 15 \
#     --slurm-dir /home/mia218/park_home_dir/HRDTimer/slurm_jobs_v2
#
# Background execution (useful for long job generation)
# $ nohup python generate_slurm_jobs.py \
#     --python-script /path/to/run_analysis.py \
#     --input-dir /path/to/vcfs \
#     --chunk-size 10 &
# -----------------------------------------------------------------------------------
