### Steps before calculation


1. The count matrix file can be created from the raw bam files by using the featureCounts command with gencode v34. It is also necessary to convert the Ensemble gene ID to HGNC gene name and remove duplicate genes.

2. Next, normalize expression levels with TPM values.

3. Then, remove unnecessary columns from the TPM matrix file to create an input file with only numeric values for feature calculation.
