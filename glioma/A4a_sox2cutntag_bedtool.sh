samplename=$1

bedtools intersect -wa -u -f 1 -b ${samplename}/sox2cre.bed -a ${samplename}/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > ${samplename}/tmp30b.bed
bedtools intersect -wa -u -f 1  -b ${samplename}/sox2noncre.bed -a ${samplename}/GSM6008250_SMNB19_SOX2_3_broadPeak.bed > ${samplename}/tmp31b.bed