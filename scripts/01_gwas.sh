#!/usr/bin/env bash
#
# Author: Frederik Lassen
#
#$ -N gwas
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/promoter_variants
#$ -o logs/gwas.log
#$ -e logs/gwas.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qe
#$ -t 17

set -o errexit
set -o nounset
module purge
#module load OpenBLAS/0.3.8-GCC-9.2.0 # for linreg
#module load OpenBLAS/0.2.15-GCC-4.9.3-2.25-LAPACK-3.6.0

# this step is required for running regression with HAIL
module load OpenBLAS/0.3.1-GCC-7.3.0-2.30
export LD_PRELOAD=/apps/eb/skylake/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas.so

source utils/bash_utils.sh

# directories
readonly in_dir="data/imputed"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/gwas"

# hail script
readonly hail_script="utils/hail_export.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in="${in_dir}/ukb_imp_chr${chr}_v3.bgen"

# output path
readonly out_prefix="${out_dir}/gwas_ukb_wes_200k_phased_chr${chr}"
readonly out="${out_prefix}.mt"

# setup hail
set_up_hail
mkdir -p ${out_dir}
python3 ${hail_script} \
    --chrom ${chr} \
    --input_path ${in} \
    --input_type "bgen" \
    --pheno_path "data/nicky_phenotypes.csv" \
    --pheno 'Hand_grip_strength_(left)_combined_white_ritish_InvNorm' 'Hand_grip_strength_(right)_combined_white_ritish_InvNorm' 'Forced_vital_capacity_(FVC)_Z-score_combined_white_ritish_InvNorm' 'FEV1-FVC_ratio_Z-score_combined_white_ritish_InvNorm'\
    --get_unrelated \
    --get_wb \
    --out_prefix "${out_dir}/test_ukb_wes_200k_phased_chr${chr}"



