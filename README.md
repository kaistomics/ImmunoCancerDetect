# ImmunoCancerDetect

Cancer detection by immunological features

---
### Description

The utility of peripheral blood features other than cell-free DNA in noninvasive multi-cancer screening remains unexplored. Here, we developed a deep learning model that extracts cancer-specific immunological features from 10,929 tumor and 10,845 normal tissue transcriptomes. This model achieves an ROC-AUC of 0.94 in discriminating between our own 1,055 blood samples from various types of cancer and 745 normal blood samples. Its capability of cancer detection is not undermined by immunological perturbation resulting from checkpoint blockade therapy. Model interpretation unveils the significance of features pertaining to the B cell receptor repertoire. Coupled with our proof-of-concept model that distinguishes cancer from infectious or autoimmune conditions, this method achieves an ultimate ROC-AUC of 0.91 when applied to our separate test dataset comprising 183 cancer or normal samples. Our results highlight the clinical implications of local and peripheral immune characteristics in the application of blood-based multi-cancer detection. 

---
### Prerequisites

* Required R packages:
  * ImmuCellAI
  * alakazam
  * dplyr
* Required python packages:
  * numpy
  * pandas
  * gseapy
  * matplotlib
* Required github code
  * <https://github.com/liulab-dfci/RIMA_pipeline/blob/0a5595070e8844a3eef65931b74aae6ee600bf6b/src/immune_repertoire/trust4_metric_functions.R>

---
### Input file preparation

* To generate gene expression matrices
  1. The count matrix file can be created from the raw bam files by using the featureCounts command with gencode v34. It is also necessary to convert the Ensemble gene ID to the HGNC gene name and remove any duplicate genes.

e.g.
````
featureCounts -T 30 -p -s 0 -t gene -g gene_id -a gencode.v34.annotation.gtf -o count.txt *.bam

````
````
import csv
import pandas as pd
import numpy as np
import sys
import os

raw = pd.read_csv("count.txt", sep='\t')
raw.rename(columns ={'Geneid':'ENSG'}, inplace=True)
raw.set_index('ENSG', inplace=True)
print("RAW data "+" shape: ", raw.shape)

## ref hg19 : 20961 genes
ref19 = pd.read_csv('gencode.v19.annotation.pcTRIG.gene_name.final.tsv', header=None, sep='\t')
ref19.rename(columns={0:"ENSG", 1:"Symbol"}, inplace=True)
ref19.set_index("ENSG", inplace=True)

mergeref19 = raw.join(ref19, how = "inner")
mergeref19.set_index('Symbol', inplace=True)
print("HGNC pcgene filtered data with hg19 shape: ", mergeref19.shape)

## filtered protein-coding gene reference(hg38) : 20677 genes
ref38 = pd.read_csv('gencode.v42.annotation.pcTRIG.gene_name.final.tsv', header=None, sep='\t')
ref38.rename(columns={0:"ENSG", 1:"Symbol"}, inplace=True)
ref38.set_index("ENSG", inplace=True)

mergeref38 = raw.join(ref38, how = "inner")
mergeref38.set_index('Symbol', inplace=True)
print("HGNC pcgene filtered data with hg38 shape: ", mergeref38.shape)

##save HGNC data
mergeref19.to_csv("input_pcTRIG_filtered_hg19.txt", sep='\t')
mergeref38.to_csv("input_pcTRIG_filtered_hg38.txt", sep='\t')
````
````
cut -f 1 input.txt | sort -n | uniq -c | sort | awk '$1 > 1 {print($2)}' > gene_dup.txt
grep -vf gene_dup.txt input.txt > input_rmdup.txt
````

  2. Next, normalize expression levels with TPM values.

e.g. R code
````
library(dplyr)
library(tidyr)

ftr.cnt <- read.table(f_in, sep="\t", stringsAsFactors=FALSE, header=TRUE)

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

ftr.tpm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(tpm=tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)
write.table(ftr.tpm, file=f_out, sep="\t", row.names=FALSE, quote=FALSE)
````


  3. Then, remove unnecessary columns from the TPM matrix file to create an input file with only numeric values for feature calculation.

e.g.

````
cut -f 2-6 --complement TPM.txt > TPM_valueonly.txt
````

* To generate input matrices for model training, combine all previously calculated features.

-> Input file format: 124 features * n samples input matrix csv file
