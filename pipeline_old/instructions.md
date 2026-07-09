# HRDTimer Pipeline Instructions

Assuming you are in the `pipeline` directory, the following steps describe how to run the HRDTimer pipeline. This example uses the **PCAWG Breast and Ovarian cancer samples** after processing through the full MutationTimeR pipeline.

---

## Step 1: Preprocess VCF Files

**Input directories and files:**

- **VCF directory:**  
  MutationTimeR output folder stored within the `hrdtimer` directory:  
  ```
  ../raw_input_data/PCAWG_MutationTimeR_output
  ```
- **Metadata file:**  
  ```
  ../data/metadata/pan_metadata_v5.csv
  ```
  The metadata file must include the following **mandatory columns**:  
  `sample`, `organ`, `HRDetect.isHRD`, `isWGD`  

  > Note: The VCF filenames in the MutationTimeR output folder should follow the convention: `sample.vcf`.

- **Output folder for preprocessed VCFs and signature analysis results:**  
  ```
  ../HRDTimer_PreProcess/PCAWG_v3
  ```

**Run preprocessing with:**

```bash
python preprocess.py \
    --input ../raw_input_data/PCAWG_MutationTimeR_output \
    --output ../HRDTimer_PreProcess/PCAWG_v3 \
    --metadata ../data/metadata/pan_metadata_v5.csv
```

> This will create a **logs folder** within the pipeline directory to track preprocessing progress.

After preprocessing, results will be saved in:

```
/n/data1/hms/dbmi/park/michail_a/HRDTimer_PreProcess/PCAWG_v3
```

---

## Step 2: Generate Slurm Jobs for HRDTimer Timing Analysis

The next step focuses on the **timing set of mutations**. To reduce runtime, jobs are split into chunks of 10–15 samples for parallel execution.

**Generate Slurm jobs with:**

```bash
python generate_slurm_jobs.py \
    --python-script run_timing_analysis.py \
    --input-dir ../HRDTimer_PreProcess/PCAWG_v3/Breast/timing \
    --chunk-size 15 \
    --slurm-dir slurm_jobs
```

> This will create a **Slurm jobs folder** at the specified location.

---

## Step 3: Submit Slurm Jobs

To submit all `.sbatch` files in the Slurm jobs folder:

```bash
cd slurm_jobs

for file in *.sbatch; do
    sbatch "$file"
done
```

---

## Step 4: Outputs

Running the pipeline will generate two main outputs:

1. **Bootstraps folder:**  
   ```
   bootstraps
   ```  
   Contains bootstrap results and signature exposures per bootstrap, used for exploration and timing analysis.

2. **Results folder:**  
   ```
   results
   ```  
   Contains batch-wise results. Example:  
   ```
   results/batch_000_28263966.csv
   ```

> All CSVs can then be merged to obtain the complete results.
