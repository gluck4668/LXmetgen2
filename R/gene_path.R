
gene_path <- function(gene_file,group1,group2,species){

#------------------Gene enriched pathways analysis-----------------------------#
spe <- c("human","mouse","rat")
species <- tolower(species)

if(species %in% spe ==F)
   { spe_text <- paste("You choose the species is",species)
     print(spe_text)
     stop_text <- paste("The species of",species,"is not found; Note: the species should be human, mouse, or rat. Please check it!" )
     stop (stop_text)} else
       dir.file <<- dplyr::case_when ( tolower(species)== "human" ~ "analysis results (human)",
                                   tolower(species)== "mouse" ~ "analysis results (mouse)",
                                   tolower(species)== "rat" ~ "analysis results (rat)",
                                   TRUE ~ "analysis results") #设为全局变量

if (dir.exists(dir.file)==FALSE)
      dir.create(dir.file)

group <<- paste("(",group1,"VS",group2,")") # 全局变量

file_type <- str_extract(gene_file,"(?<=[.]).*")
if(tolower(file_type)=="txt")
  gene_df_0 <- read_table(gene_file) else
   gene_df_0 <- eval(str2expression(paste0("read.",file_type,"(gene_file)")))

colnames(gene_df_0)[1] <- "gene_id"
gene_df <- data.frame(distinct(gene_df_0, gene_id, .keep_all = TRUE))

gene_n <- dplyr::filter(gene_df,tolower(str_sub(gene_df$gene_id,1,1)) %in% 0:9) # 数字开头的基因
gene_L <- dplyr::filter(gene_df,tolower(str_sub(gene_df$gene_id,1,1)) %in% letters) # 字母开头的基因

# hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
gene_human <- c(gene_n$gene_id,toupper(gene_L$gene_id))

# mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
gene_mouse <- c( gene_n$gene_id, str_to_title(tolower(gene_L$gene_id)) )  # str_to_title()首字母大写
gene_rat <- c( gene_n$gene_id, str_to_title(tolower(gene_L$gene_id)) )

if(tolower(species)== "human")
      {gene_set <- gene_human
      orgdb = org.Hs.eg.db
      organ="hsa"}

if(tolower(species)== "mouse")
    {gene_set <- gene_mouse
    orgdb = org.Mm.eg.db
    organ="mmu"}


if(tolower(species)== "rat")
    {gene_set <- gene_rat
    orgdb = org.Rn.eg.db
    organ="rno"}

#---kegg pathways enrichment--------
gene_ENTREZID<-bitr(gene_set,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=orgdb)
gene_ENTREZID <- na.omit(gene_ENTREZID)
gene_ENT_h <- paste("The gene_ENTREZID_data","(",species,")",".xlsx")
gene_ENT_h <-paste0(dir.file,"/",gene_ENT_h)
write.xlsx(gene_ENTREZID,gene_ENT_h)
kegg_gene_df <- enrichKEGG(gene_ENTREZID$ENTREZID, organism = organ,
                               keyType = 'kegg',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH',
                               qvalueCutoff = 0.2,
                               minGSSize = 3,
                               maxGSSize = 3500,
                               use_internal_data = F)

pathways_geneID <- kegg_gene_df@result

symbol <- setReadable(kegg_gene_df, OrgDb = orgdb, keyType="ENTREZID")
pathways_geneSymbol <- symbol@result
kegg_all_pathways <- na.omit(pathways_geneSymbol) #全局变量

#--------------------------------------------
#colnames(kegg_all_pathways)[c(1,9,10)] <- c("pathwayID","geneSymbol","Count")

spe_type <- case_when(tolower(species)== "human" ~ "sapiens",
                      tolower(species)== "mouse" ~ "musculus",
                      tolower(species)== "rat" ~ "Rattus")

spe_exist <- grepl(spe_type,kegg_all_pathways$Description,ignore.case = T) %>% any()

if (spe_exist)
  kegg_all_pathways$Description <- str_extract(kegg_all_pathways$Description,".*(?= -)")

#把分数字符串变小数
str_fun <- function(i){eval(str2expression(i))}
kegg_all_pathways$GeneRatio <- purrr::map(kegg_all_pathways$GeneRatio,str_fun)
kegg_all_pathways$GeneRatio <- as.numeric(kegg_all_pathways$GeneRatio) %>% round(.,6)
kegg_all_pathways <- dplyr::arrange(kegg_all_pathways,desc(GeneRatio))

file_path_name <- paste("The KEGG pathways_data","(",species,")",".xlsx")
file_path_name <-paste0(dir.file,"/",file_path_name)
write.xlsx(kegg_all_pathways,file_path_name)

row_n <- nrow(kegg_all_pathways)

if(row_n<30)
  title_gene_text <- paste("The genes enriched pathways","(",group1,"VS",group2,")") else
  title_gene_text <- paste("TOP 30 genes enriched pathways","(",group1,"VS",group2,")")

title_size <- case_when(nrow(kegg_gene_df)>30 ~12,
                            nrow(kegg_gene_df)>20 ~12,
                            TRUE ~11)

xy_size <- case_when(nrow(kegg_gene_df)>30 ~9.5,
                         nrow(kegg_gene_df)>20 ~10.5,
                         TRUE ~11)

kegg_mytheme<-theme_bw()+
      theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size),
            panel.grid = element_blank(),
            panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
            axis.ticks.length = unit(1.5,units = "mm"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
      theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

kegg_xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size,angle =0,hjust=1))+
      theme(axis.text.y = element_text(face="bold",color="black",size=xy_size))+
      theme(legend.text=element_text(face="bold",color="black",size=xy_size))

if(row_n<30)
      path_n <- row_n else
        path_n <- 30

kegg_df <- kegg_all_pathways[1:path_n,]

#--------------kegg all pathways------------------------#
kegg_pathways <- ggplot(kegg_df)+
      geom_point(aes(x=GeneRatio,
                     y=fct_reorder(Description,GeneRatio),
                     color=-log10(pvalue),size=Count))+
      scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
      labs(x = 'GeneRatio', y = '',title=title_gene_text)+
      kegg_mytheme+kegg_xytheme

kegg_pathways

path_name <- paste("The KEGG pathways","(",species,")",".png")
path_name <-paste0(dir.file,"/",path_name)
ggsave(path_name, kegg_pathways,width=1200, height =1000, dpi=150,units = "px")

#------------kegg metabolism pathways----------------#
#筛选出含有"METABOL"字段的行
is.catetory <- grepl("category",names(kegg_all_pathways),ignore.case = T) %>% any()
if(is.catetory)
   meta_path_all <- dplyr::filter(kegg_all_pathways, category=="Metabolism") else
     meta_path_all <-kegg_all_pathways[grepl("metabol", kegg_all_pathways$Description),]

if(any(meta_path_all$Description=="Metabolic pathways"))
    meta_path_all <- meta_path_all[-grep("Metabolic pathways",meta_path_all$Description,ignore.case = T),]

kegg_meta_all <<- meta_path_all  # 设为全局变量

if(nrow(kegg_meta_all)>0){

    meta_name <- paste("The metabolism pathways_data","(",species,")",".xlsx")
    meta_name <-paste0(dir.file,"/",meta_name)
    write.xlsx(kegg_meta_all,meta_name)

    if(nrow(kegg_meta_all)>=30)
      {kegg_meta <- kegg_meta_all[1:30,]
       title_meta_text <- paste("Top 30 KEGG metabolism pathways","(",group1,"VS",group2,")")} else
       {kegg_meta <- kegg_meta_all
       title_meta_text <- paste("KEGG metabolism pathways","(",group1,"VS",group2,")")}

    kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)

    kegg_meta_pathways <- ggplot(kegg_meta)+
      geom_point(aes(x=GeneRatio,y=reorder(Description,GeneRatio),
                     color=-log10(pvalue),size=Count))+
      scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
      labs(x = 'GeneRatio', y = '',title=title_meta_text)+
      kegg_mytheme+kegg_xytheme

    kegg_meta_pathways

    meta_path_name <- paste("The metabolism pathways","(",species,")",".png")
    meta_path_name <-paste0(dir.file,"/",meta_path_name)

    ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")



       } else
print("There is no metabolism pathway!")

kegg_gene <- kegg_all_pathways[,c("Description","pvalue")]
colnames(kegg_gene) <- c("pathway","log2p_gene")
kegg_gene$log2p_gene = -log2(kegg_gene$log2p_gene)

kegg_gene_path <<- data.frame(kegg_gene)

}



