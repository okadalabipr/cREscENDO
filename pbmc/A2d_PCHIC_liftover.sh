samplename=$1

liftOver ${samplename}/bait_hg19_tag.bed ${samplename}/hg19ToHg38.over.chain ${samplename}/bait_hg38_tag.bed bait_unlifted.bed
liftOver ${samplename}/oe_hg19_tag.bed ${samplename}/hg19ToHg38.over.chain ${samplename}/oe_hg38_tag.bed oe_unlifted.bed