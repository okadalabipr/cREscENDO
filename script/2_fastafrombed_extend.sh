samplename=$1

fastaFromBed -fi /sqfs/work/K2313/u6b872/hg38/genome.fa -bed ${samplename}/peaks_extend.bed -fo ${samplename}/peaks_extend.fasta