# cREscENDO

cREscENDO is a deep learning-based computational tool for determining cRE regions and their activities at the single cell level from single sample 10X Single Cell Multiome ATAC + Gene Expression data.

Installation should be performed in the following manner

```
pip install -r requirements.txt
```

to install the necessary files.
In addition, please install bedtools(https://bedtools.readthedocs.io/en/latest/content/installation.html),sync_batchnorm(https://github.com/vacancy/ Synchronized-BatchNorm-PyTorch), please refer to their respective websites for installation.



### How to run
To execute, run the following command
```
pip install -r requirements.txt
```

Perform the learning. As with other deep learning tools, cREscENDO should be trained using a GPU.
If cross-validation is required, execute the following command
```
pip install -r requirements.txt
```
If cross-validation is not required, execute the following commands
```
pip install -r requirements.txt
```
### Output
The single-cell level cRE activity is provided in a file named ""
This is a 3-dimensional matrix consisting of [genes, cells, peaks]
The cRE matrix for the top 60% of cREs, starting from the one with the highest activity maximum, is provided as the file named “”.
This matrix is a two-dimensional matrix of [genes,peaks]
1 represents cRE, 0 represents non-cRE, and -1 represents missing values. 
(The number of peaks validated for each gene is different, resulting in missing values.)


As for the peak dimension, the target promoter is always entered at index 0.
Then, starting from index 1, the regions with the smallest genome coordinates are stored in order.

For the correspondence between the gene index and the gene name, please refer to the file named "" (Gene indexes start with 0)
The file “” contains an index of which ATAC-seq peak is associated with each gene.
The file “” contains the index of the target promoter of each gene.
The index of the ATAC-seq peak corresponds to the index (starting from 0) of the bed file.

In addition, the code used in the paper is stored in this github repository.

### for_motif
Code for motif analysis

### pbmc
Code for Fig.2 and Fig.3A

### glioma
Code for Fig.3B-J
