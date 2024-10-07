samplename=$1
genomepath=$2

fastaFromBed -fi ${genomepath} -bed ${samplename}/peaks_extend.bed -fo ${samplename}/peaks_extend.fasta