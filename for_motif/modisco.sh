samplename=$1
prefix=$2
motiffile=$3

modisco motifs -a ${samplename}/${prefix}contribution.npy -s ${samplename}/${prefix}inputmtx.npy -n 50000 -o ${samplename}/${prefix}modiscoout.h5
mkdir ${samplename}/${prefix}_modiscoreport
modisco report -i ${samplename}/${prefix}modiscoout.h5 -o ${samplename}/${prefix}_modiscoreport -m ${motiffile}