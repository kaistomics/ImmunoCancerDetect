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

  2. Next, normalize expression levels with TPM values.

  3. Then, remove unnecessary columns from the TPM matrix file to create an input file with only numeric values for feature calculation.

* To generate input matrices for model training, combine all previously calculated features.

-> Input file format: 124 features * n samples input matrix csv file
