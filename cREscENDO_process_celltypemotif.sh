samplename=$1
cellcluspath=$2
dir=$(cd $(dirname $0); pwd)

cd ${dir}

python for_motif/13_DAR_select_celltype.py ${samplename} ${cellcluspath}
python for_motif/enhaselect_celltype.py ${samplename} 0
python for_motif/enhaselect_celltype.py ${samplename} 1
python for_motif/enhaselect_celltype.py ${samplename} 2
python for_motif/enhaselect_celltype.py ${samplename} 3
python for_motif/enhaselect_celltype.py ${samplename} 4
python for_motif/enhaselect_celltype.py ${samplename} 5