#!/bin/bash
#SBATCH -c 8                               # Request 16 cores
#SBATCH -t 0-10:05                         # Runtime in D-HH:MM format
#SBATCH -A park_contrib
#SBATCH -p park                              # Partition to run in
#SBATCH --mem=8G                         # Memory total in MB (for all cores)
#SBATCH -o log/timeR_example%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e log/timeR_example%j.err                 # File to which STDERR will be written, including job ID (%j)

module load gcc/14.2.0 python/3.13.1 R/4.4.2
module load samtools
export R_LIBS_USER="~/park_dglodzik/Renvs/R-4.4.2-IRkernel/library";


# the line should be: tumor ID, tumor bam, normal ID, normal bam, (37|38)
ALIQUOT_ID="987528ac-437a-4eb8-a335-4f2076d5c006" 
ICGC_SPECIMEN_ID="SA20573" 
BASE_PURITY="0.818"
RESULT_PATH="/home/dg204/park_dglodzik/TimeR/example_take2/"
VCF_FN="/home/dg204/park_data/ICGC/SNV_indel_calls/final_consensus_12oct_passonly/snv_mnv/987528ac-437a-4eb8-a335-4f2076d5c006.consensus.20160830.somatic.snv_mnv.vcf.gz"
CN_FN="/home/dg204/park_data/ICGC/CNV_all_callers_included/987528ac-437a-4eb8-a335-4f2076d5c006.consensus.20170119.somatic.cna.annotated.txt"
TUMOR_ID=""

mkdir -p "$RESULT_PATH"

# runs vanilla MutationTimeR
Rscript run_timer.R \
    "${ALIQUOT_ID}" \
    "${BASE_PURITY}" \
    "${RESULT_PATH}" \
    "${VCF_FN}" \
    "${CN_FN}" \
    "hg37" \
    "${TUMOR_ID}"
# writes per sample results to the following sufolders of $RESULT_PATH: /data/




# this is the output from TimeR part (above)
RDATA="${RESULT_PATH}/data/${ALIQUOT_ID}.RData"

# create appropriate folders for the WGD part
PLOT_FOLDER=${RESULT_PATH}/hrdtimer_plots/
if [ ! -d "$PLOT_FOLDER" ]; then
  mkdir -p "$PLOT_FOLDER"
fi
RESULT_HRDTIMER_PATH="$RESULT_PATH/wgd_results/"
if [ ! -d "$RESULT_HRDTIMER_PATH" ]; then
  mkdir -p "$RESULT_HRDTIMER_PATH"
fi

# runs MutationTimeR's logic for timing genome doubling
echo "hrdtimer ..."
Rscript run_timer_wgd.R \
"$ALIQUOT_ID" \
"$RDATA" \
"$PLOT_FOLDER" \
"$RESULT_HRDTIMER_PATH" \
"$BASE_PURITY"
echo "completed hrdtimer ..."
