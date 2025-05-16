#!/bin/bash
#SBATCH -c 2                               # Request 16 cores
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -A park_contrib
#SBATCH -p park                            # Partition to run in
#SBATCH --mem=32G                         # Memory total in MB (for all cores)
#SBATCH -o log/timeR_icgc_%A_%a.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e log/timeR_icgc_%A_%a.err                 # File to which STDERR will be written, including job ID (%j)
#SBATCH --array=2-9                                     # 2-9%2   
#SBATCH -x compute-p-17-34,compute-p-17-43,compute-p-17-174,compute-p-17-01,compute-p-17-02,compute-p-17-03,compute-p-17-04,compute-p-17-05,compute-p-17-06,compute-p-17-07,compute-p-17-08,compute-p-17-26,compute-p-17-27,compute-p-17-28,compute-p-17-29,compute-p-17-30,compute-p-17-38,compute-p-17-39,compute-p-17-40,compute-p-17-41,compute-p-17-42,compute-p-17-44,compute-p-17-45,compute-p-17-46,compute-p-17-180,compute-p-17-181,compute-p-17-182,compute-p-17-183,compute-p-17-184,compute-p-17-185,compute-p-17-31,compute-p-17-32,compute-p-17-33,compute-p-17-35,compute-p-17-36,compute-p-17-37

module load gcc/9.2.0 python/3.10.11 R/4.3.1
module load samtools
export R_LIBS_USER="~/R-4.3.1/library"

SAMPLE_TABLE_FN=/home/dg204/park_dglodzik/ICGC/CGAP/sample_table.tsv
LNO=${SLURM_ARRAY_TASK_ID}
NTH_LINE=$(cat ${SAMPLE_TABLE_FN} | head -${LNO} | tail -1)
echo $NTH_LINE
# the line should be: tumor ID, tumor bam, normal ID, normal bam, (37|38)
ALIQUOT_ID=$(echo $NTH_LINE | awk -F" " '{ print $1}')  
PURITY=$(echo $NTH_LINE | awk -F" " '{ print $6}')  
TUMOR_ID=$(echo $NTH_LINE | awk -F" " '{ print $15}')  
RESULT_PATH="/home/dg204/park_dglodzik/TimeR_ICGC/"
VCF_FN="/home/dg204/park_dglodzik/ICGC/CGAP/SNVs_HC/${ALIQUOT_ID}_snv_high_conf.vcf.gz"
CN_FN="/home/dg204/park_dglodzik/ICGC/CGAP/purple_tumor_cnv_somatic_tsv/${ALIQUOT_ID}.tsv.gz"

echo "ALIQUOT_ID $ALIQUOT_ID"
echo "PURITY $PURITY"
echo "TUMOR_ID $TUMOR_ID"

Rscript run_timer.R \
    "$ALIQUOT_ID" \
    "$PURITY" \
    "$RESULT_PATH" \
    "$VCF_FN" \
    "$CN_FN" \
    "hg38" \
    "$TUMOR_ID"

RDATA="${RESULT_PATH}/data/${ALIQUOT_ID}.RData"
PLOT_FOLDER=${RESULT_PATH}/hrdtimer_plots/

if [ ! -d "$PLOT_FOLDER" ]; then
  mkdir -p "$PLOT_FOLDER"
fi
RESULT_HRDTIMER_PATH="$RESULT_PATH/hrdtimer_results/"
if [ ! -d "$RESULT_HRDTIMER_PATH" ]; then
  mkdir -p "$RESULT_HRDTIMER_PATH"
fi

Rscript run_hrdtimer.R \
"$ALIQUOT_ID" \
"$RDATA" \
"$PLOT_FOLDER" \
"$RESULT_HRDTIMER_PATH" \
"$PURITY"