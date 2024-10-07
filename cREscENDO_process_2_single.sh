samplename=$1
dir=$(cd $(dirname $0); pwd)

cd ${dir}
python script/7_1split_single.py ${samplename}
