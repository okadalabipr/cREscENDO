# cREscENDO

cREscENDO is a deep learning-based computational tool for determining cRE regions and their activities at the single cell level from single sample 10X Single Cell Multiome ATAC + Gene Expression data.

Installation should be performed in the following manner

```
## We recommend that you run it in a virtual environment. For example, use the following command.
## conda create --name crescendo_env python=3.8
## conda activate crescendo_env
git clone https://github.com/okadalabipr/cREscENDO
cd cREscENDO
pip install -r requirements.txt
```

to install the necessary files.
In addition, please install bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html) , Meme suite (https://meme-suite.org/meme/doc/install.html?man_type=web) , sync_batchnorm (https://github.com/vacancy/Synchronized-BatchNorm-PyTorch), please refer to their respective websites for installation.<br>
("sync_batchnorm" should be placed inside the “script” directory.)



```
conda install bioconda::bedtools
conda install bioconda::meme
cd ..
git clone https://github.com/vacancy/Synchronized-BatchNorm-PyTorch
cp -r Synchronized-BatchNorm-PyTorch/sync_batchnorm cREscENDO/script/
```

### How to run

In order to run this tool, you must provide a TSS list in bed file format. hg38's TSS list is included in this repository.<br>
(TSS list should be placed inside the working directory.)

To execute, run the following command
```
bash cREscENDO_process_1.sh [path/to/working/directory] [path/to/10X_multiome_h5.file] [path/to/ref/genome/fasta]
```

Perform the learning. As with other deep learning tools, cREscENDO should be trained using a GPU.<br>
If cross-validation is required, execute the following commands
```
bash cREscENDO_process_2.sh [path/to/working/directory]
bash cREscENDO_process_3.sh [path/to/working/directory]
```
If cross-validation is not required, execute the following commands
```
bash cREscENDO_process_2_single.sh [path/to/working/directory] 0
bash cREscENDO_process_3_single.sh [path/to/working/directory] 0
```
Cross-validation requires a lot of calculation time. If you want to divide the job or run the calculation in parallel due to restrictions on the execution environment, you can do it as follows.
```
##Please execute the following commands in separate jobs.
bash cREscENDO_process_2_single.sh [path/to/working/directory] 0
bash cREscENDO_process_2_single.sh [path/to/working/directory] 1
bash cREscENDO_process_2_single.sh [path/to/working/directory] 2
bash cREscENDO_process_2_single.sh [path/to/working/directory] 3
bash cREscENDO_process_2_single.sh [path/to/working/directory] 4
bash cREscENDO_process_2_single.sh [path/to/working/directory] 5
bash cREscENDO_process_2_single.sh [path/to/working/directory] 6
bash cREscENDO_process_2_single.sh [path/to/working/directory] 7
bash cREscENDO_process_2_single.sh [path/to/working/directory] 8
bash cREscENDO_process_2_single.sh [path/to/working/directory] 9

##Finally, run the following command to aggregate the results of cross-validation.
bash cREscENDO_process_3.sh [path/to/working/directory]
```


### Output

The single-cell level cRE activity is provided in a file named "Deeplift_full_ver2_all.npy"<br>
This is a 3-dimensional matrix consisting of [genes, cells, peaks]<br>
<br>
The cRE matrix for the top 60% of cREs, starting from the one with the highest activity maximum, is provided as the file named “allgrad_ssep_max.npy”.<br>
This matrix is a two-dimensional matrix of [genes,peaks]<br>
"1" represents cRE, "0" represents non-cRE, and "-1" represents missing values in this matrix. <br>
(The number of peaks validated for each gene is different, resulting in missing values.)<br>
<br>
As for the peak dimension, the target promoter is always entered at index 0.<br>
Then, starting from index 1, the regions with the smallest genome coordinates are stored in order.<br>
<br>
For the correspondence between the gene index and the gene name, please refer to the file named "pair_300000.csv" (Gene indexes start with 0)<br>
The file “pair_300000.csv” also contains an index of which ATAC-seq peak is associated with each gene.<br>
The file “pair_promoter.csv” contains the index of the target promoter of each gene.<br>
The index of the ATAC-seq peak corresponds to the index (starting from 0) of the bed file ("peaks.bed").<br>
<br>

-------------------------------
In addition, the code used in the paper is stored in this github repository.

### for_motif
Code for motif analysis

### pbmc
Code for Fig.2 and Fig.3A

### glioma
Code for Fig.3B-J
