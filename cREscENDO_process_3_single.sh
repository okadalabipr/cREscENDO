samplename=$1
foldidx=$2
dir=$(cd $(dirname $0); pwd)

cd ${dir}

python script/11fm_DeepLift_fullcell_merge_t2_raw_single.py ${samplename} ${foldidx}
python script/13_promenhatag.py ${samplename}
