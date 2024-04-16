
LXmetgen2 <- function(gene_file,meta_file,group1,group2,species){

#--------list all the packages that have been installed

  ist_packs <- function(){

  all_packages <- installed.packages()[,1]

# To judge whether a package was installed. If not, it will be installed.
  com_packs <- c("pak","devtools","BiocManager","plyr","openxlsx","dplyr","psych","ggplot2",
            "ggrepel","VennDiagram","ggvenn","RColorBrewer","ggthemes","rticles","httr","magrittr",
            "roxygen2","XML","RCurl","curl","stringr","patchwork","ggpubr","scales","readr")

  not_com <- com_packs[! com_packs %in% all_packages]


  if(length(not_com)>0){
    com_fun <- function(x){install.packages(x)}
    sapply(not_com,com_fun,simplify = T)
  }


#-------- Install the development version of "tidyverse" from GitHub
  if(!"tidyverse" %in% all_packages)
      pak::pak("tidyverse/tidyverse")
     library(tidyverse)


  bio_packs <- c("DOSE","clusterProfiler","do","enrichplot",
                   "pathview","BiocParallel","GO.db","KEGGREST",
                   "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db")
  # human: "org.Hs.eg.db"
  # mouse: "org.Mn.eg.db"
  # rat: "org.Rn.eg.db"

  not_bio <- bio_packs[! bio_packs %in% all_packages]

  if(length(not_bio)>0){
    bio_fun <- function(x){BiocManager::install(x)}
    sapply(not_bio,bio_fun,simplify = T)
  }


#----批量library
  lib_fun <- function(x){library(x,character.only = T)}
  sapply(c(com_packs,bio_packs),lib_fun)

}

ist_packs()

#------------------Gene enriched pathways analysis-----------------------------#
gene_path(gene_file,group1,group2,species)

#------------------metabolite pathways analysis-----------------------------#
meta_path(meta_file)

#---------Jiont analysis of the gene and metabolite enriched pathways---------#

jiont_analysis(kegg_gene_path,meta_pathway)

show_result

}



