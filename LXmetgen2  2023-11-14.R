
LXmetgen2 <- function(gene_file,meta_file,group1,group2,species){

#--------list all the packages that have been installed

  all_packages <- installed.packages()[,1]

# To judge whether a package was installed. If not, it will be installed.
  com_packs <- c("pak","devtools","BiocManager","plyr","eoffice","openxlsx","dplyr","psych","ggplot2",
            "ggrepel","VennDiagram","ggvenn","RColorBrewer","ggthemes","rticles","httr","magrittr",
            "roxygen2","XML","RCurl","curl","stringr","patchwork","ggpubr","scales") 

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
  sapply(c(com_packs,bio_packs),lib_fun,simplify = T)


 rm(all_packages,com_packs,not_com,bio_packs,not_bio,lib_fun)
#------------------Gene enriched pathways analysis-----------------------------#

   spe <- c("human","mouse","rat")
   species <- tolower(species)

   if(species %in% spe ==F)
   { spe_text <- paste("You choose the species is",species)
     print(spe_text)
     stop_text <- paste("The species of",species,"is not found; Note: the species should be human, mouse, or rat. Please check it!" )
     stop (stop_text)} else
       dir.file <- dplyr::case_when ( tolower(species)== "human" ~ "analysis results (human)",
                                   tolower(species)== "mouse" ~ "analysis results (mouse)",
                                   tolower(species)== "rat" ~ "analysis results (rat)",
                                   TRUE ~ "analysis results")

    if (dir.exists(dir.file)==FALSE)
      dir.create(dir.file)

    group <- paste("(",group1,"VS",group2,")")

    gene_df_0 <- read.xlsx(gene_file)
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
    kegg_all_pathways <- cbind(pathways_geneID[,1:8],pathways_geneSymbol[,8],pathways_geneID[,9])
    
 #--------------------------------------------
    colnames(kegg_all_pathways)[c(1,9,10)] <- c("pathwayID","geneSymbol","Count")
    
    spe_type <- case_when(tolower(species)== "human" ~ "sapiens",
                          tolower(species)== "mouse" ~ "musculus",
                          tolower(species)== "rat" ~ "Rattus")
    
    spe_exist <- grepl(spe_type,kegg_all_pathways$Description,ignore.case = T) %>% table()
    
    if (grepl("TRUE",names(spe_exist)))
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
    kegg_meta <- kegg_gene_df

    #把Description转换成大写
    upper <- toupper(kegg_meta@result$Description)

    #筛选出含有"METABOL"字段的行
    kegg_meta <- dplyr::filter(kegg_meta@result,grepl("METABOL",upper))

    if(nrow(kegg_meta)>0)
    {kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)

    #把分数字符串变小数
    Gene_Ratio <- as.data.frame(apply(str_split(kegg_meta$GeneRatio,"/",simplify = T),2,as.numeric))
    Gene_Ratio <- Gene_Ratio[,1]/Gene_Ratio[,2]

    kegg_meta$Gene_Ratio <- Gene_Ratio
    kegg_meta <- kegg_meta[-1,]

    title_meta_text <- paste("KEGG metabolism pathways","(",group1,"VS",group2,")")

    kegg_meta_pathways <- ggplot(kegg_meta)+
      geom_point(aes(x=Gene_Ratio,y=reorder(Description,Gene_Ratio),
                     color=-log10(pvalue),size=Count))+
      scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
      labs(x = 'GeneRatio', y = '',title=title_meta_text)+
      kegg_mytheme+kegg_xytheme

    kegg_meta_pathways

    meta_path_name <- paste("The metabolism pathways","(",species,")",".png")
    meta_path_name <-paste0(dir.file,"/",meta_path_name)

    ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")

    #重新命名
    colnames(kegg_all_pathways)[1] <- "Pathways_ID"
    colnames(kegg_meta)[1] <- "Pathways_ID"

    kegg_meta01 <- kegg_meta[,c(1,10)]
    kegg_meta02 <- inner_join(kegg_all_pathways,kegg_meta01,by="Pathways_ID")

    #按Gene_Ratio降序排序
    kegg_meta02 <- kegg_meta02[order(-kegg_meta02$Gene_Ratio),]
    rownames(kegg_meta02) <- c(1:nrow(kegg_meta02))

    kegg_meta_path <- dplyr::arrange(kegg_meta02,desc(Gene_Ratio))

    meta_name <- paste("The metabolism pathways_data","(",species,")",".xlsx")
    meta_name <-paste0(dir.file,"/",meta_name)

    write.xlsx(kegg_meta_path,meta_name)} else
      print("There is no metabolism pathway!")

#---------Jiont analysis of the gene and metabolite enriched pathways---------#

if(file.exists(meta_file)==FALSE)
  {stop <- paste(meta_file,"is not found, or its format is csv (must be exlsx); please check it.")
  stop(stop)} else
    meta_path_0 <- read.xlsx(meta_file)

colnames(meta_path_0) <- c("pathways","Total","Expected","Hits",
                           "pvalue","-log10(p)","Holm.adjust","FDR","Impact")
meta_path <- meta_path_0[,c(1,5)]
colnames(meta_path) <- c("pathways","meta_pvalue")

kegg_meta <- kegg_all_pathways[,c(2,5)]
colnames(kegg_meta) <- c("pathways","kegg_pvalue")

kegg_meta_Venn <- inner_join(data.frame(kegg_meta),data.frame(meta_path),
                             by="pathways")

colnames(kegg_meta_Venn) <- c("Pathway","P_gene","P_meta")

gene_path <- kegg_meta_Venn[c("Pathway","P_gene")]
gene_P <- as.numeric(gene_path$P_gene)
class(gene_P)
gene_path$P_gene <- -log2(gene_P)

gene_path$types <- rep("genes",nrow(gene_path))
colnames(gene_path) <- c("Pathway","-log2(Pvalue)","types")

meta_path <- kegg_meta_Venn[c("Pathway","P_meta")]
meta_P <- as.numeric(meta_path$P_meta)
meta_path$P_meta <- -log2(meta_P)
meta_path$types <- rep("metabolites",nrow(meta_path))
colnames(meta_path) <- c("Pathway","-log2(Pvalue)","types")

gene_meta_path <- bind_rows(gene_path,meta_path)
colnames(gene_meta_path) <- c("Pathways","minus_log2_Pvalue","types")

g_m_name <- paste("The gene_meta_pathways_data","(",species,")",".xlsx")
g_m_name <-paste0(dir.file,"/", g_m_name)

write.xlsx(gene_meta_path, g_m_name)

y_p <- max(gene_meta_path$minus_log2_Pvalue)

height_y <- y_p*1.2

nrow_path <- nrow(gene_meta_path)/2

joint_title_size <- case_when(nrow_path>=30 ~12,
                              nrow_path>=20 ~12,
                              TRUE ~14)

joint_x_size <- case_when(nrow_path>=20 ~9,
                          nrow_path>=10 ~10,
                          TRUE ~12)

joint_y_size <- case_when(nrow_path>=20 ~12,
                          nrow_path>=10 ~12,
                          TRUE ~12)

joint_legend_size <- case_when(nrow_path>=30 ~12,
                               nrow_path>=20 ~12,
                               TRUE ~12)


bar_width <- case_when(nrow_path>=30 ~0.9,
                       nrow_path>=20 ~0.8,
                       TRUE ~0.7)

f1 <- ggplot(gene_meta_path, aes(x = Pathways, y = minus_log2_Pvalue,fill=types))+
  geom_bar(position = "dodge",stat = "identity",width = bar_width)+
  scale_fill_manual(values=c("#008b8b","#f08080"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(colour = "black", face="bold",size=12)) +
  labs(x="",y = ("-log2(Pvalue)"),title = paste('Joint-Pathways Analysis',group))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, height_y))

f1

log05 <- -log2(0.05)
log01 <- -log2(0.01)

line1 <- geom_hline(yintercept = c(log05),
                    linewidth = 0.6,
                    color = "blue",
                    lty = "dashed")
line2 <- geom_hline(yintercept = c(log01),
                    linewidth = 0.6,
                    color = "red",
                    lty = "dashed")

y1 <- geom_text(x=nrow(gene_meta_path)/2-2.5,y=log05+0.8,label = c("p<0.05"),
                size=4,color="blue",fontface="italic")
y2 <- geom_text(x=nrow(gene_meta_path)/2-2.5,y=log01+0.8,label = c("p<0.01"),
                size=4,color="blue",fontface="italic")

f2 <- f1+line1+line2+y1+y2
f2

mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =joint_title_size),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=2), "cm")
       )+
  theme(plot.title = element_text(hjust = 0.5))


xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=joint_x_size,angle =45,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=joint_y_size))

legend_theme <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = joint_legend_size, face = "bold"),
  legend.direction = "vertical",
  #legend.position = c(0.5,0.9),
  legend.background = element_blank()
)

f3 <- f2+mytheme+xytheme

f3

f3_name <- paste("The gene_metabolite_Joint_pathways 01","(",species,")",".png")
f3_name <-paste0(dir.file,"/", f3_name)

ggsave(f3_name,f3,width=1200, height =1000, dpi=150,units = "px")


f4 <- f3+
  labs(fill="")+
  theme(legend.direction = "horizontal",
        legend.position = c(0.5,0.92),
        legend.text = element_text(size=14,face = "bold") )

f4

f4_name <- paste("The gene_metabolite_Joint_pathways 02","(",species,")",".png")
f4_name <-paste0(dir.file,"/", f4_name)

ggsave(f4_name,f4,width=1200, height =1000, dpi=150,units = "px")

print("--------------------------------------------------------------")

print_text <- paste("The results can be found in the folder of",dir.file)

print(print_text)

f4

}



