#!/bin/bash
#PBS -q SQUID-H
#PBS --group=K2313
#PBS -l elapstim_req=1:00:00
#PBS -l gpunum_job=1
#PBS -m eb
#PBS -M k-mrkm@protein.osaka-u.ac.jp
cd $PBS_O_WORKDIR
source /system/apps/rhel8/cpu/Anaconda3/2020.11/etc/profile.d/conda.sh

samplename="/sqfs2/cmc/1/work/K2313/u6b872/pbmc_10fold1"
inputname="pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"

conda activate test-env7
python 0_preparation.py ${samplename} ${inputname}
python 1_peak_extended.py ${samplename}
bash 2_fastafrombed_extend.sh ${samplename}
python 3_fastatoread_extend.py ${samplename}
python 4_make_genelist.py ${samplename}
python 5_testtrain_split.py ${samplename}
python 5_gex_preprocess.py ${samplename}
python 6_PCA_rna_embedding_raw.py ${samplename}
conda deactivate