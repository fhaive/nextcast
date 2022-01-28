# INPUT:
# file_names: vector of strings with path to the files containing the genes ranked by INfORM
# N: Number of top ranked genes to include in the output file
# experiment_names: vector of strings used to label the lists contained in the files
# group: vector of strings used to label the lists contained in the files

# WHAT THE SCRIPT DOES:
# It will convert these files in a format ready to use for the FunMappOne analysis.
# In this script, only the top N genes in the ranks will be selected


#OUTPUT:
# The script will create two excel files. The first will contain 
# one sheet for each expression matrix (corresponding to each experimental condition 
# (e.g. treatment) present in the data. The second will contain 
# one sheet for each phenodata table describing the samples in the first excel files.
# Expression matrices and phenodata files are matched by position in the two excel files

library(readxl)

gene_table_refactory = function(gene_table, n=N){
  gene_table = as.data.frame(gene_table)
  gene_table = gene_table[,c("SYMBOL","LFC")]
  gene_table = unique(gene_table)
  rownames(gene_table) = gene_table$SYMBOL
  gene_table = gene_table[1:n,]
  return(gene_table)
}

#### INPUT

# list of files with gene ranking to include in the FunMappOne analysis 

file_names = c("data/network_5_24/Gene_Tables_2021-12-15.xlsx",
              "data/network_5_48/Gene_Tables_2021-12-15.xlsx",
              "data/network_5_72/Gene_Tables_2021-12-16.xlsx",
              "data/network_10_24/Gene_Tables_2021-12-16.xlsx",
              "data/network_10_48/Gene_Tables_2021-12-16.xlsx",
              "data/network_10_72/Gene_Tables_2021-12-15.xlsx",
              "data/network_20_24/Gene_Tables_2021-12-15.xlsx",
              "data/network_20_48/Gene_Tables_2021-12-15.xlsx",
              "data/network_20_72/Gene_Tables_2021-12-15.xlsx")

# N is the number of genes that will be selected.
N = 200

# This is a vector of strings that will be used in the FunMappOne plot to name the
# different experiments. They need to follow the same order of the file names
experiment_names = c("MWCNT_5_24",
                     "MWCNT_5_48",
                     "MWCNT_5_72",
                     "MWCNT_10_24",
                     "MWCNT_10_48",
                     "MWCNT_10_72",
                     "MWCNT_20_24",
                     "MWCNT_20_48",
                     "MWCNT_20_72")


group =  c(1,1,1,2,2,2,3,3,3)

inputs = list()
for(i in 1:length(file_names)){
  inputs[[i]] = read_excel(file_names[[i]])
}

names(inputs) = experiment_names


for(i in 1:length(inputs)){
  inputs[[i]] = gene_table_refactory(inputs[[i]])
}

# change group assignment here.  
pheno = data.frame(samples = experiment_names, group = group)

data_frame_list = append(inputs, list("group" = pheno))

writexl::write_xlsx(data_frame_list, path =  "data/data_for_funmappone.xlsx")
