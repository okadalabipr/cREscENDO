samplename=$1
cellcluspath=$2
dir=$(cd $(dirname $0); pwd)

cd ${dir}

python script/13_promenhatag.py ${samplename}
python for_motif/13_DAR_select_celltype.py ${samplename} ${cellcluspath}
python enhaselect_celltype.py ${samplename} 0
python enhaselect_celltype.py ${samplename} 1
python enhaselect_celltype.py ${samplename} 2
python enhaselect_celltype.py ${samplename} 3
python enhaselect_celltype.py ${samplename} 4
python enhaselect_celltype.py ${samplename} 5