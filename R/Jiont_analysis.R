
jiont_analysis <- function(kegg_gene_path,meta_pathway){

#---------Jiont analysis of the gene and metabolite enriched pathways---------#
kegg_meta_Venn <- inner_join(kegg_gene_path,meta_pathway,by="pathway")
names(kegg_meta_Venn)

gene_path <- kegg_meta_Venn[c("pathway","log2p_gene")]
gene_path$types <- rep("genes",nrow(gene_path))
colnames(gene_path) <- c("Pathway","-log2(Pvalue)","types")

meta_path <- kegg_meta_Venn[c("pathway","log2p_meta")]
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


#---------ggplot point----------------
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
  theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))


xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=10,angle =0,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=10))+
  theme(legend.text=element_text(face="bold",color="black",size=10))

gene_meta_df <- mutate(gene_meta_path,Ratio=c(kegg_meta_Venn$gene_Ratio,kegg_meta_Venn$meta_Ratio))
colnames(gene_meta_df) <- c("Pathways","log2pvalue","type","Ratio")

f5 <- ggplot(gene_meta_df)+
  geom_point(aes(x=Ratio,
                 y=Pathways,
                 shape=type,
                 color=log2pvalue,
                 size=Ratio))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'Ratio', y = 'Pathways',title=paste('Joint-Pathways Analysis',group),
       col="-log2(pvalue)",shape="Type",size="Ratio")+ # 修改图例名称
  mytheme+xytheme

f5

f5_name <- paste("The gene_metabolite_Joint_pathways 03","(",species,")",".png")
f5_name <-paste0(dir.file,"/", f5_name)
ggsave(f5_name,f5,width=1450, height =1200, dpi=150,units = "px")

#------------------
line05 <- geom_vline(xintercept = c(log05),
                    linewidth = 0.5,
                    color = "black",
                    lty = "dashed")

txt05 <- geom_text(x=log05+4,y=2.5,label = c("p<0.05"),
                size=5,color="blue",fontface="italic")

f6 <- ggplot(gene_meta_df)+
  geom_point(aes(x=log2pvalue,
                 y=Pathways,
                 shape=type,
                 color=log2pvalue,
                 size=Ratio))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = '-log2(pvalue)', y = 'Pathways',title=paste('Joint-Pathways Analysis',group),
       col="-log2(pvalue)",shape="Type",size="Ratio")+ # 修改图例名称
  mytheme+xytheme+
  line05+txt05

f6

f6_name <- paste("The gene_metabolite_Joint_pathways 04","(",species,")",".png")
f6_name <-paste0(dir.file,"/", f6_name)
ggsave(f6_name,f6,width=1450, height =1200, dpi=150,units = "px")

#-----------------------

print("--------------------------------------------------------------")

print_text <- paste("The results can be found in the folder of",dir.file)

print(print_text)

show_result <<- f4 #全局变量

}



