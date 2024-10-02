#!/bin/bash
#PBS -q SQUID-H
#PBS --group=K2313
#PBS -l elapstim_req=0:15:00
#PBS -l gpunum_job=8
#PBS -m eb
#PBS -M k-mrkm@protein.osaka-u.ac.jp
cd $PBS_O_WORKDIR
source /system/apps/rhel8/cpu/Anaconda3/2020.11/etc/profile.d/conda.sh

samplename="/sqfs2/cmc/1/work/K2313/u6b872/pbmc_10fold1"

conda activate test-env7
python 11fm_DeepLift_fullcell_merge_t2_raw.py ${samplename}
conda deactivate