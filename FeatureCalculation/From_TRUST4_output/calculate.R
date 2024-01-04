# repertoire features
library(dplyr)

setwd('trust4_out')
source("src/immune_repertoire/trust4_metric_functions.R")
args = commandArgs(trailingOnly=TRUE)
filename <- args[1]
cdr3 <- read.table(file = filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# filter out count of cdr < 3 and add the phenotype info for downstream comparison
phenotype <- "phenotype"
sample <- str_remove(filename, "_report.tsv")
cdr3 <- subset(cdr3, count > 3) %>%
  mutate(V = as.character(V), J = as.character(J), C = as.character(C), CDR3aa = as.character(CDR3aa)) %>%
  mutate(clinic = as.character(phenotype), sample = as.character(sample))

# determine whether the cdr animo acid is complete or partial
cdr3$is_complete <- sapply(cdr3$CDR3aa, function(x) ifelse(x != "partial" && x != "out_of_frame" && !grepl("^_",x) && !grepl("^\\?", x),"Y","N"))

# exact the TCR and BCR
cdr3.bcr <- subset(cdr3, grepl("^IG",V) | grepl("^IG",J) | grepl("^IG",C))
cdr3.tcr <- subset(cdr3, grepl("^TR",V) | grepl("^TR",J) | grepl("^TR",C))

# add lib size and clinic traits
cdr3.bcr <- cdr3.bcr %>% mutate(lib.size = sum(count))
cdr3.tcr <- cdr3.tcr %>% mutate(lib.size = sum(count))

# split BCR into heavy chain and light chain
cdr3.bcr.heavy <- subset(cdr3.bcr, grepl("^IGH",V) | grepl("^IGH",J) | grepl("^IGH",C))
cdr3.bcr.light <- subset(cdr3.bcr, grepl("^IG[K|L]",V) | grepl("^IG[K|L]",J) | grepl("^IG[K|L]",C))

cdr3.tcr.trb <- subset(cdr3.tcr, grepl("^TRB",V) | grepl("^TRB",J) | grepl("^TRB",C))

sample_bcr_cluster <- BuildBCRlineage(sampleID = sample, Bdata = cdr3.bcr.heavy, start=3, end=10)

# TCR CPK
tcrcpk <- aggregate(CDR3aa ~ sample+clinic+lib.size, cdr3.tcr, function(x) length(unique(x))) %>%
  mutate(TCR_CPK = signif(CDR3aa/(lib.size/1000),4))
bcrcpk <- aggregate(CDR3aa ~ sample+clinic+lib.size, cdr3.bcr, function(x) length(unique(x))) %>%
  mutate(BCR_CPK = signif(CDR3aa/(lib.size/1000),4))

# BCR clonality and entropy
single_sample_bcr_clonality <- getClonality(sample, cdr3.bcr.heavy, start=3, end=10)

# TCR clonality and entropy
single_sample_tcr_clonality <- getClonalityTCR(sample,cdr3.tcr)

SHM.ratio <- getSHMratio(sample_bcr_cluster)

st.Ig <- cdr3.bcr.heavy %>%
  group_by(clinic,sample) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%   #as.numeric(sample.clones[filename,2])
  dplyr::filter(C != "*" & C !=".") %>%
  group_by(sample, C) %>%
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))
tmp <- t(st.Ig)
colnames(tmp) <- tmp[2,]

# chemical properties of amino acids
library(alakazam)
                    
bcrAA <- aminoAcidProperties(cdr3.bcr, property = c("length", "gravy", "bulk", "aliphatic", "polarity", "charge", "basic","acidic", "aromatic"),seq='CDR3aa', nt=FALSE)
tcrAA <- aminoAcidProperties(cdr3.tcr, property = c("length", "gravy", "bulk", "aliphatic", "polarity", "charge", "basic","acidic", "aromatic"),seq='CDR3aa', nt=FALSE)
igh <- subset(cdr3.bcr, grepl("^IGH",V))
igk <- subset(cdr3.bcr, grepl("^IGK",V))
igl <- subset(cdr3.bcr, grepl("^IGL",V))
tra <- subset(cdr3.tcr, grepl("^TRA",V))
trb <- subset(cdr3.tcr, grepl("^TRB",V))
trd <- subset(cdr3.tcr, grepl("^TRD",V))
trg <- subset(cdr3.tcr, grepl("^TRG",V))

# combine all features of receptor sequences
final <- cbind(data.frame(sample) %>% mutate(clinic =  as.character(phenotype), TCR.lib.size = tcrcpk$lib.size, TCR.CDR3.aa = tcrcpk$CDR3aa, TCR.CPK = tcrcpk$TCR_CPK) %>%
                 mutate(TCR.clonality = single_sample_tcr_clonality[2], TCR.entropy = single_sample_tcr_clonality[3]) %>%
                 mutate(BCR.lib.size = bcrcpk$lib.size, BCR.CDR3.aa = bcrcpk$CDR3aa, BCR.CPK = bcrcpk$BCR_CPK) %>%
                 mutate(BCR.clonality = single_sample_bcr_clonality[2], BCR.entropy = single_sample_bcr_clonality[3]) %>%
                 mutate(BCR.cluster.number = length(sample_bcr_cluster), BCR.SHM.rate = SHM.ratio) %>%
                 mutate(TCR.length.mean = mean(tcrAA$CDR3aa_aa_length, na.rm=TRUE), TCR.gravy.mean = mean(tcrAA$CDR3aa_aa_gravy, na.rm=TRUE),
                        TCR.bulk.mean = mean(tcrAA$CDR3aa_aa_bulk, na.rm=TRUE), TCR.polarity.mean = mean(tcrAA$CDR3aa_aa_polarity, na.rm=TRUE),
                        TCR.charge.mean = mean(tcrAA$CDR3aa_aa_charge, na.rm=TRUE), TCR.basic.mean = mean(tcrAA$CDR3aa_aa_basic, na.rm=TRUE),
                        TCR.aliphatic.mean = mean(tcrAA$CDR3aa_aa_aliphatic, na.rm=TRUE), TCR.acidic.mean = mean(tcrAA$CDR3aa_aa_acidic, na.rm=TRUE),
                        TCR.aromatic.mean = mean(tcrAA$CDR3aa_aa_aromatic, na.rm=TRUE), BCR.length.mean = mean(bcrAA$CDR3aa_aa_length, na.rm=TRUE),
                        BCR.gravy.mean = mean(bcrAA$CDR3aa_aa_gravy, na.rm=TRUE), BCR.bulk.mean = mean(bcrAA$CDR3aa_aa_bulk, na.rm=TRUE),
                        BCR.polarity.mean = mean(bcrAA$CDR3aa_aa_polarity, na.rm=TRUE), BCR.charge.mean = mean(bcrAA$CDR3aa_aa_charge, na.rm=TRUE),
                        BCR.basic.mean = mean(bcrAA$CDR3aa_aa_basic, na.rm=TRUE), BCR.aliphatic.mean = mean(bcrAA$CDR3aa_aa_aliphatic, na.rm=TRUE),
                        BCR.acidic.mean = mean(bcrAA$CDR3aa_aa_acidic, na.rm=TRUE), BCR.aromatic.mean = mean(bcrAA$CDR3aa_aa_aromatic, na.rm=TRUE))  %>%
                 mutate(IGH = sum(igh$count)/sum(cdr3.bcr$count), IGK = sum(igk$count)/sum(cdr3.bcr$count), IGL = sum(igl$count)/sum(cdr3.bcr$count),
                        TRA = sum(tra$count)/sum(cdr3.tcr$count), TRB = sum(trb$count)/sum(cdr3.tcr$count),
                        TRD = sum(trd$count)/sum(cdr3.tcr$count), TRG = sum(trg$count)/sum(cdr3.tcr$count)), data.frame(tmp)['Num.Ig',, drop=FALSE])

write.table(final, noquote(paste(sample, "_features.txt", sep = "")), sep="\t", row.names = F)
