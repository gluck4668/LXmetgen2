\name{LXmetgen2}
\alias{LXmetgen2}
\title{LXmetgen2::joint analysis for transcriptome and metabolome }
\usage{
 LXmetgen2 (gene_file,meta_file,group1,group2,species)
 }
\description{
LXmetgen2 can be used for the joint analysis of transcriptome and metabolome.
}

\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXmetgen2")
library(LXmetgen2)

??LXmetgen2
#---------------------------------
data(gene_list_file_example) # a list of gene SYMBOL
data(meta_data_file_example) # a metabolism pathwas data obtained from the MetaboAnalyst online (https://www.metaboanalyst.ca/MetaboAnalyst/)

#------------------------------

setwd("D:/Desktop/R_example/LXmetgen_example")

gene_file="gene_data.xlsx"
meta_file="meta_data.xlsx"

group1="Model"
group2="Normal"

species= "mouse"  # The species should be "human", "mouse", or "rat"


LXmetgen2 (gene_file,meta_file,group1,group2,species)

}
