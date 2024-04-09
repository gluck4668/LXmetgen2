
meta_path <- function(meta_file){

#---------metabolite pathways---------#
meta_type <- str_extract(meta_file,"(?<=[.]).*")
if(tolower(meta_type)=="txt")
  meta_path_0 <- read_table(meta_file) else
    meta_path_0 <- eval(str2expression(paste0("read.",meta_type,"(meta_file)")))

#----数据来源于IMPaLA------------------
is.impala <- grepl("overlap",colnames(meta_path_0),ignore.case = T) %>% any()
if(is.impala){
  meta_path <- dplyr::filter(meta_path_0,pathway_source=="KEGG")
  meta_path <- meta_path[,c("pathway_name","P_metabolites")]
  colnames(meta_path) <- c("pathway","log2p_meta")
  meta_path$pathway <- str_extract(meta_path$pathway,".*(?= -)") %>% trimws()
  meta_path$log2p_meta = -log2(meta_path$log2p_meta)
}


#----数据来源于MetabolAnalyst------------------
is.metabo <- grepl("impact",colnames(meta_path_0),ignore.case = T) %>% any()
if(is.metabo){
  meta_path <- meta_path_0[,c(1,5)]
  colnames(meta_path) <- c("pathway","log2p_meta")
  meta_path$log2p_meta = -log2(meta_path$log2p)
}


meta_pathway <<- data.frame(meta_path)

}



