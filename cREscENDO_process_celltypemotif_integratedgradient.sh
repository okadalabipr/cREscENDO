dir=$(cd $(dirname $0); pwd)

cd ${dir}

samplename=$1
enhaclusfile=$2
cellclusfile=$3
celltype=$4
prefix=$5
motiffile=$6

python for_motif/14_tfmodiscoinput_clip_new2.py ${samplename} ${enhaclusfile} ${cellcluspath} ${celltype} ${prefix}
bash for_motif/modisco.sh ${samplename} ${prefix} ${motiffile}