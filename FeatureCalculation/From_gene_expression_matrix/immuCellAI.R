args=(commandArgs(TRUE))
f_name <- as.character(args[1])

f_in <- paste("./", f_name, "_input.txt", sep='')
f_out <- paste("./", f_name, "_output.txt", sep='')

library(ImmuCellAI)

# Open input files (Symbol must shown without index 'Symbol')
input <- read.csv(f_in, sep='\t', row.names=1, stringsAsFactors = FALSE, header=T)

# Run ImmuCellAI_new 
res <- ImmuCellAI_new(sample = input, data_type = "rnaseq", customer = 0, group_tag = 0, response_tag = 0)

# Save output
output <- res$Sample_abundance
write.table(output, file=f_out, sep='\t', row.names=T, col.names=NA, quote=F)
