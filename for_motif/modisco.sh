modisco motifs -a integratedGradient.npy -s inputseq.npy -n 50000 -o modiscoout.h5
mkdir modiscoout_modiscoreport
modisco report -i modiscoout.h5 -o modiscoout_modiscoreport -m JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
