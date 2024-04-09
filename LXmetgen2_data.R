
library(openxlsx)

gene_list_file_example <- read.xlsx("gene_data.xlsx")
meta_pathways_MetabolAnalyst <- read.xlsx("meta_data.xlsx")
meta_pathways_Impala <- read.csv("meta_ORA_results.csv")

usethis::use_data(gene_list_file_example,overwrite = T)
usethis::use_data(meta_pathways_MetabolAnalyst,overwrite = T)
usethis::use_data(meta_pathways_Impala,overwrite = T)

rm(list=ls())

data(gene_list_file_example,overwrite = T)
data(meta_pathways_MetabolAnalyst,overwrite = T)
data(meta_pathways_Impala,overwrite = T)

