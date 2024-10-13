samplename=$1

mkdir ${samplename}/SOX2
fimo -oc ${samplename}/SOX2 --text -motif MA0143.5 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt ${samplename}/peaks_extend.fasta > ${samplename}/SOX2/fimo.tsv

mkdir ${samplename}/ZEB1
fimo -oc ${samplename}/ZEB1 --text -motif MA0103.4 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt ${samplename}/peaks_extend.fasta > ${samplename}/ZEB1/fimo.tsv

mkdir ${samplename}/CREM
fimo -oc ${samplename}/CREM --text -motif MA0609.3 --thresh 0.001 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt ${samplename}/peaks_extend.fasta > ${samplename}/CREM/fimo.tsv


