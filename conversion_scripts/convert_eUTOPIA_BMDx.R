# INPUT:
# This script takes in input the phenodata table of the samples preprocessed in
# eUTOPIA, the preprocessed expression matrix and the results of the differential
# analysis unfiltered (meaning that the file contains all the genes profiled in
# the experiment and their corresponding log-fold-change and the p-values).

# WHAT THE SCRIPT DOES:
# It will convert these files in a format ready to use for the BMDx analysis.
# In this script, only the genes that are differentially expressed (named DEG genes) will be included.
# Since in eUTOPIA multiple comparisons can be made, a gene can be or not 
# DEG in all of them. The strategy implemented in this script
# is to count how many times the a gene is a DEG one across the different
# comparisons and keep in the data only those that are DEG in at least N comparisons.

# DETAILS ON THE PARAMETERS:
# The threshold to define a gene as DEG (log_fc_th and pval_th) and the number of comparison (N)
# are set-up by the user. A gene is considered DEG if its absolute log-fold-change 
# is higher than log_fc_th and its adjusted pvalue is less than pval_th.
# In the example the lof_fc_th threshold is set to 0.58 that correspond to 
# a fold-change of 1.5

#OUTPUT:
# The script will create two excel files. The first will contain 
# one sheet for each expression matrix (corresponding to each experimental condition 
# (e.g. treatment) present in the data. The second will contain 
# one sheet for each phenodata table describing the samples in the first excel files.
# Expression matrices and phenodata files are matched by position in the two excel files

library(readr)
library(rio)

# READ INPUTS: update paths and file names with your own

### Read phenodata table
pheno <- read_delim("data/GSE146708_metadata.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
pheno = as.data.frame(pheno)
rownames(pheno) = pheno$sample_id

### Read expression matrix as stored from eUTOPIA
expMat <- read.delim("data/MQ_Expression_Matrix_Aggregated.txt")
expMat= as.matrix(expMat)
expMat = expMat[,rownames(pheno)]

### Read differential analysis results are stored from eUTOPIA
deg <- import_list("data/MQ_Differential_Expression_Tables_UNFILTERED.xlsx")
deg = deg[2:length(deg)]
names(deg) = gsub(pattern = "[0-9]\\.", replacement = "",names(deg))
names(deg) # this shows the names of the comparisons performed in eUTOPIA

# PATH WHERE OUTPUT ARE GOING TO BE STORED: can be updated with your own path
path_output_files = "data/"

# INPUTS: 
log_fc_th = 0.58 #threshold on the log-fold-change
pval_th = 0.05 # threshold on the adjusted p-values

# number of comparison in which a gene has to be differentially expressed
# in order to be included in the files
number_of_comparisons = 4

# name of the experiment that will be used in the excel file sheet
exp_name = "MWCNT"

# names of the phenodata columns specifying dose, time points, the sample id, 
# and the type of exposure (e.g. MWCNT) 
pheno_dose_col = "dose"
pheno_time_col = "time_point"
pheno_sample_col = "sample_id"
pheno_exposure_col = "NA" 


# The nTimes vector is used to count how many times the genes are differentially 
# expressed across the different comparisons
nTimes = rep(0, nrow(deg[[1]]))
names(nTimes) =deg[[1]]$GeneSymbol

# for each comparison performed in the eUTOPIA analysis
for(i in 1:length(deg)){
  deg_i = deg[[i]]
  rownames(deg_i) = deg_i$GeneSymbol
  
  # identify differentially expressed genes for the i-th comparison
  deg_f = deg_i[abs(deg_i$logFC)> log_fc_th & deg_i$adj.P.Val< pval_th,]
  
  # Update gene counts
  nTimes[rownames(deg_f)] = nTimes[rownames(deg_f)]+1
}
nTimes = sort(nTimes, decreasing = T)

# Only the genes that are differentially expressed in at least number_of_comparisons
# will be included in the file
good_genes = names(nTimes)[nTimes>=number_of_comparisons]

# if in the phenodata table the column specifying the exposure is is not present, 
# an extra column is added. This is the case of the sample data we are using 
if(pheno_exposure_col=="NA"){
  pheno = cbind(pheno, rep(exp_name,nrow(pheno)))
  colnames(pheno)[ncol(pheno)] = "exposure_type"
  pheno_exposure_col = "exposure_type"
}

# if in the phenodata table the column specifying the time is is not present, 
# an extra column is added. This is not the case of the sample data we are using 
if(pheno_time_col=="NA"){
  pheno = cbind(pheno, rep(1:nrow(pheno)))
  colnames(pheno)[ncol(pheno)] = "time"
  pheno_time_col = "time"
}
colnames(expMat) = pheno[,pheno_sample_col]

# This for will create a list of expression matrices and a list of corresponding
# phenodata table. This will be stored in two excel files. The first will contain 
# one sheet for each expression matrix (corresponding to each treatment - in the 
# case of the sample data only one: MWCNT). The second will contain 
# one sheet for each phenodata table (corresponding to the samples of each treatment - 
# in the case of the sample data only one: MWCNT). 

exp_list = list()
pheno_list = list()
for(cr in unique(pheno[,pheno_exposure_col])){
  
  idx = which(pheno[,pheno_exposure_col] %in% cr)
  ph_cr = pheno[idx,c(pheno_sample_col,pheno_dose_col, pheno_time_col)]
  exp_cr = expMat[good_genes,ph_cr[,pheno_sample_col]]
  
  exp_list[[cr]] = as.data.frame(exp_cr)
  exp_list[[cr]] = cbind(" "=rownames(exp_cr), exp_list[[cr]])
  
  pheno_list[[cr]] = as.data.frame(ph_cr)
  
}

# NAMES OF THE STORED FILES: can be updated
exp_out_file = paste(path_output_files,"/MWCNT_expression_data_only_DEG_bmdx.xlsx",sep = "")
pheno_out_file = paste(path_output_files,"/MWCNT_pheno_data_bmdx.xlsx",sep = "")

# SAVING FILES
writexl::write_xlsx(x = exp_list,path = exp_out_file)
writexl::write_xlsx(x = pheno_list,path = pheno_out_file)
