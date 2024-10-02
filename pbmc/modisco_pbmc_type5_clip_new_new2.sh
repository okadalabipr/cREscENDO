#!/usr/bin/env bash
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=0-4:00:00
#SBATCH --output=logs/slurm-%x_%A.out
#SBATCH --error=logs/slurm-%x_%A.error

source /software/2020/software/miniconda3/23.9.0-0/etc/profile.d/conda.sh
conda activate modiscoenv

modisco motifs -a /users/ken.murakami/workspace/pbmcnew/try1/deeplift5_integratedGradient_norm_new2.npy -s /users/ken.murakami/workspace/pbmcnew/try1/inputmtx5_integratedGradient_norm_new2.npy -n 50000 -o /users/ken.murakami/workspace/pbmcnew/try1/pbmc_modisco_clus5_new_new2.h5

module load meme/5.1.1-foss-2018b-python-3.6.6

mkdir /users/ken.murakami/workspace/pbmcnew/try1/pbmc_modisco_clus5_new_new2_modiscoreport
modisco report -i /users/ken.murakami/workspace/pbmcnew/try1/pbmc_modisco_clus5_new_new2.h5 -o /users/ken.murakami/workspace/pbmcnew/try1/pbmc_modisco_clus5_new_new2_modiscoreport -m /users/ken.murakami/workspace/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
