samplename=$1
foldidx=$2
dir=$(cd $(dirname $0); pwd)

cd ${dir}
python script/7_1split_single.py ${samplename} ${foldidx}
