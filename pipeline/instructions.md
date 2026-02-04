# HRDTimer Pipeline Instructions

Assuming you are in the **pipeline directory**, the following steps describe how to run the HRDTimer pipeline. These instructions are for analyzing **PCAWG Breast and Ovarian samples** after running the full MutationTimeR pipeline.

---

## Step 1: Preprocess VCF Files

**Input directories and files:**

- **VCF directory:**  
  ```
  /n/data1/hms/dbmi/park/michail_a/HRDTimer/raw_input_data/PCAWG_MutationTimeR_output
  ```
- **Metadata file:**  
  ```
  /n/data1/hms/dbmi/park/michail_a/HRDTimer/data/metadata/pan_metadata_v5.csv
  ```
  The metadata file must include the following **mandatory columns**:  
  `sample`, `organ`, `HRDetect.isHRD`, `isWGD`

- **Output folder for preprocessed VCFs and signature analysis results:**  
  ```
  /n/data1/hms/dbmi/park/michail_a/HRDTimer_PreProcess/PCAWG_v3
  ```

**Command to run preprocessing:**

```bash
python preprocess.py \
    --input /home/mia218/park_home_dir/HRDTimer/raw_input_data/PCAWG_MutationTimeR_output \
    --output /n/data1/hms/dbmi/park/michail_a/HRDTimer_PreProcess/PCAWG_v3 \
    --metadata /home/mia218/park_home_dir/HRDTimer/data/metadata/pan_metadata_v5.csv
```

> This will create a **logs folder** within the pipeline directory to track the progress of preprocessing jobs.  

After preprocessing, the following directory will be created:  
```
/n/data1/hms/dbmi/park/michail_a/HRDTimer_PreProcess/PCAWG_v3
```

---

## Step 2: Generate Slurm Jobs for HRDTimer Timing Analysis

For the HRDTimer pipeline, we focus on the **timing set of mutations**. To parallelize the analysis and reduce runtime, jobs are split per 10–15 samples.  

**Command to generate Slurm jobs:**

```bash
python generate_slurm_jobs.py \
    --python-script /home/mia218/park_home_dir/HRDTimer/pipeline/run_timing_analysis.py \
    --input-dir /n/data1/hms/dbmi/park/michail_a/HRDTimer_PreProcess/PCAWG_v3/Breast/timing \
    --chunk-size 15 \
    --slurm-dir /home/mia218/park_home_dir/HRDTimer/pipeline/slurm_jobs
```

> This will create a **Slurm jobs folder** at the specified location.

---

## Step 3: Submit Slurm Jobs

To submit all the generated `.sbatch` files:

```bash
cd /home/mia218/park_home_dir/HRDTimer/pipeline/slurm_jobs

for file in *.sbatch; do
    sbatch "$file"
done
```

---

## Step 4: Outputs

Running the pipeline will generate two main outputs:

1. **Bootstraps folder:**  
   ```
   /home/mia218/park_home_dir/HRDTimer/pipeline/bootstraps
   ```  
   Stores bootstrap results and signature exposures per bootstrap for exploration and timing.

2. **Results folder:**  
   ```
   /home/mia218/park_home_dir/HRDTimer/pipeline/results
   ```  
   Contains results split by batch. Example:  
   ```
   /home/mia218/park_home_dir/HRDTimer/pipeline/results/batch_000_28263966.csv
   ```

All CSVs can then be merged to obtain the full results.
