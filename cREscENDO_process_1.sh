samplename=$1
inputname=$2

python 0_preparation.py ${samplename} ${inputname}
python 1_peak_extended.py ${samplename}
bash 2_fastafrombed_extend.sh ${samplename}
python 3_fastatoread_extend.py ${samplename}
python 4_make_genelist.py ${samplename}
python 5_testtrain_split.py ${samplename}
python 5_gex_preprocess.py ${samplename} ${inputname}
python 6_PCA_rna_embedding_raw.py ${samplename}
