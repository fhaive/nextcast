# INPUT:
# This script takes in input the phenodata table of the samples preprocessed in
# eUTOPIA, the preprocessed expression matrix and the results of the differential
# analysis unfiltered (meaning that the file contains all the genes profiled in
# the experiment and their corresponding log-fold-change and the p-values).

# WHAT THE SCRIPT DOES:
# It will convert these files in a format ready to use for the INfORM analysis.
# This script will generate the files required to create one co-expression network
# for each one of the comparison performed in the eUTOPIA analysis. 
# In this script, only a certain amount of genes will be included and they will
# be the same in all the datasets.
# Since in eUTOPIA multiple comparisons can be made, a gene can be or not 
# DEG in all of them. The strategy implemented in this script
# is to count how many times the a gene is a DEG one across the different
# comparisons, rank them based on how many times they are DEG and select only the 
# genes that are in the top N position of the list.

# DETAILS ON THE PARAMETERS:
# The threshold to define a gene as DEG (logFC and p.value) and the number of top positions (top_position)
# are set-up by the user. A gene is considered DEG if its absolute log-fold-change 
# is higher than logFC and its adjusted pvalue is less than p.value.
# In the example the logFC threshold is set to 0.58 that correspond to 
# a fold-change of 1.5

#OUTPUT:
# The script will create two excel files. The first will contain 
# one sheet for each expression matrix (corresponding to each experimental condition 
# (e.g. treatment) present in the data. The second will contain 
# one sheet for each phenodata table describing the samples in the first excel files.
# Expression matrices and phenodata files are matched by position in the two excel files

library(readr)
library(rio)

### READ EUTOPIA OUTPUTS
pheno <- read_delim("data/GSE146708_metadata.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
pheno = as.data.frame(pheno)
rownames(pheno) = pheno$sample_id

expMat <- read.delim("data/MQ_Expression_Matrix_Aggregated.txt")
expMat= as.matrix(expMat)
expMat = expMat[,rownames(pheno)]

deg <- import_list("data/MQ_Differential_Expression_Tables_UNFILTERED.xlsx")
deg = deg[2:length(deg)]
names(deg) = gsub(pattern = "[0-9]\\.", replacement = "",names(deg))

path_output_files = "data/"
contrast = names(deg)

# this is a boolean, if TRUE all the genes are used;
# if FALSE, only the DEG genes are used. DEG are selected by using logFC and p.value thresholds
all_genes  = FALSE

# INPUTS: 
logFC = 0.58 #threshold on the log-fold-change
p.value = 0.05# threshold on the adjusted p-values
# number of genes that will be included in the final results
top_positions = 1000

# column of the pheno data containing the grouping used when performing the comparisons
# and the sample ids.
pheno_group_col = 4
pheno_sample_col = 1

# index of the columns that contains the genes ID in the differential expression analysis results
gene_symbol_col = 8

# The nTimes vector is used to count how many times the genes are differentially 
# expressed across the different comparisons
nTimes = rep(0, nrow(deg[[1]]))
names(nTimes) =deg[[1]]$GeneSymbol

for(i in 1:length(deg)){
  di = deg[[i]]
  rownames(di) = di[, gene_symbol_col]
  di = di[abs(di$logFC)>logFC & di$adj.P.Val<p.value,]
  nTimes[rownames(di)] = nTimes[rownames(di)]+1
}

nTimes = sort(nTimes, decreasing = T)
# Only the genes that are in the first top position of the rank are selected
common_genes = names(nTimes)[1:top_positions]

for(cr in contrast){
  elem = unlist(strsplit(x = cr,split = "-"))
  # identifying the index of the samples used in this specific contrast
  idx = which(pheno[,pheno_group_col] %in% elem)
  ph_cr = pheno[idx,]
  ph_cr = ph_cr[order(ph_cr$dose),]
  exp_cr = expMat[,ph_cr[,pheno_sample_col]]
  di = deg[[cr]]
  rownames(di) = di[, gene_symbol_col]
  
  if(all_genes  == FALSE){
    di = di[common_genes,]
    exp_cr = exp_cr[common_genes,]
  }
  
  di = di[,c("logFC","P.Value")]
  write.table(exp_cr, paste(path_output_files, cr,"_exp.txt",sep = ""), quote = FALSE, sep = "\t")
  write.table(di, paste(path_output_files,cr,"_deg.txt",sep = ""), quote = FALSE, sep = "\t")
}


