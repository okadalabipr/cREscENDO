samplename=$1
inputname=$2
genomepath=$3
prefix=$4
dir=$(cd $(dirname $0); pwd)

cd ${dir}
python glioma/0_preparation.py ${samplename} ${inputname} ${prefix}
bash glioma/0a_liftover.sh ${samplename}
python glioma/0a_hg19tohg38.py ${samplename}
python script/1_peak_extended.py ${samplename}
bash script/2_fastafrombed_extend.sh ${samplename} ${genomepath}
python script/3_fastatoread_extend.py ${samplename}
python script/4_make_genelist.py ${samplename}
python script/5_testtrain_split.py ${samplename}
python glioma/5_gex_preprocess.py ${samplename} ${inputname} ${prefix}
python script/6_PCA_rna_embedding_raw.py ${samplename}
