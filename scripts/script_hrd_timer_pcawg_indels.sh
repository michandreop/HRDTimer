#!/bin/bash
#SBATCH -c 8                               # Request 16 cores
#SBATCH -t 0-4:05                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8G                         # Memory total in MB (for all cores)
#SBATCH -o log/indel_timeR_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e log/indel_timeR_%j.err                 # File to which STDERR will be written, including job ID (%j)
#SBATCH --array=6-321                                      # 2-2779 2-321 for breast/ovarian

module load gcc/9.2.0 python/3.10.11 R/4.3.1
module load samtools
export R_LIBS_USER="~/R-4.3.1/library"

SAMPLE_TABLE_FN=~/projects/brca_timing/data/processed/breast_ovarian_sel_cols.csv
#SAMPLE_TABLE_FN=~/projects/rsignatures/data/processed/all_pcawg_sel_cols.csv
LNO=${SLURM_ARRAY_TASK_ID}
NTH_LINE=$(cat ${SAMPLE_TABLE_FN} | head -${LNO} | tail -1)

echo $NTH_LINE

# the line should be: tumor ID, tumor bam, normal ID, normal bam, (37|38)
ALIQUOT_ID=$(echo $NTH_LINE | awk -F"," '{ print $1}') 
ICGC_SPECIMEN_ID=$(echo $NTH_LINE | awk -F"," '{ print $2}')  
PURITY=$(echo $NTH_LINE | awk -F"," '{ print $3}')  

Rscript run_hrd_timer_pcawg.R \
"$ALIQUOT_ID" \
"$ICGC_SPECIMEN_ID" \
"$PURITY" \
"indels"
# should write the updated vcf files here: ~/park_dglodzik/TimeR_indels/