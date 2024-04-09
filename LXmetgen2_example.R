
if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXmetgen2")
library(LXmetgen2)

??LXmetgen2
#---------------------------------
data(gene_list_file_example) # a list of gene SYMBOL
data(meta_pathways_MetabolAnalyst) # obtained from the MetaboAnalyst online (https://www.metaboanalyst.ca/MetaboAnalyst/)
data(meta_pathways_Impala) # obtained from the Impala online (http://impala.molgen.mpg.de/)
#------------------------------

rm(list=ls())

gene_file="gene_data.txt"
meta_file="meta_data.xlsx"
meta_file="meta_ORA_results.csv"

group1="Model"
group2="Normal"

species= "rat"  # The species should be "human", "mouse", or "rat"

devtools::load_all()

LXmetgen2(gene_file,meta_file,group1,group2,species)




