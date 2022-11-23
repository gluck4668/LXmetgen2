
setwd("D:/R-lin study/R packages/LXmetgen2")
library(openxlsx)

gene_list_file_example <- read.xlsx("gene_data.xlsx")
meta_data_file_example <- read.xlsx("meta_data.xlsx")

usethis::use_data(gene_list_file_example,overwrite = T)
usethis::use_data(meta_data_file_example,overwrite = T)

rm(list=ls())

data(gene_list_file_example)
data(meta_data_file_example)

