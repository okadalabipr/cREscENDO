samplename=$1
inputname=$2
genomepath=$3
dir=$(cd $(dirname $0); pwd)

cd ${dir}
python script/0_preparation_ratio50.py ${samplename} ${inputname}
python script/1_peak_extended.py ${samplename}
bash script/2_fastafrombed_extend.sh ${samplename} ${genomepath}
python script/3_fastatoread_extend.py ${samplename}
python script/4_make_genelist.py ${samplename}
python script/5_testtrain_split.py ${samplename}
python script/5_gex_preprocess_ratio50.py ${samplename} ${inputname}
python script/6_PCA_rna_embedding_raw.py ${samplename}
