### Autor: Leandro de Mattos Pereira 12/02/2018
#### Fun??es usadas: Phyloseq, Melt and ggplot, biomformat,dplyr"
rm(list=ls())
############## Se necess?rio, instalar as bibliotecas a seguir com 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")

library("phyloseq")
library("biomformat")
library("dplyr")
library ("ggplot2")
library("ggsci")
library("gridExtra")
library("qiime2R")
library ("stringr")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")

metadata <- read_q2metadata("/home/mattoslmp/A_Erica_LEMM/qiime2/manifest")
SVs <-read_qza("/home/mattoslmp/A_Erica_LEMM/qiime2/OUTPUT_DADA_QIIME2/dada2_table_man2_250.qza")
taxonomy<-read_qza("/home/mattoslmp/A_Erica_LEMM/qiime2/OUTPUT_DADA_QIIME2/dada2_taxonomy_250pb.qza")
tree<-read_qza("/home/mattoslmp/A_Erica_LEMM/qiime2/OUTPUT_DADA_QIIME2/fasttree-tree_250.qza")
phylo <- qza_to_phyloseq(features="/home/mattoslmp/A_Erica_LEMM/qiime2/OUTPUT_DADA_QIIME2/dada2_table_man2_250.qza", tree="/home/mattoslmp/A_Erica_LEMM/qiime2/OUTPUT_DADA_QIIME2/fasttree-tree_250.qza", taxonomy = "/home/mattoslmp/A_Erica_LEMM/qiime2/OUTPUT_DADA_QIIME2/dada2_taxonomy_250pb.qza", "/home/mattoslmp/A_Erica_LEMM/qiime2/manifest")
source("/home/mattoslmp/A_Erica_LEMM/qiime2/Script_artigo_Raquel/scritps_v1/scripts/phyloseq_to_edgeR.R")


#######################################################
rank_names(phylo)
taxtable = tax_table(phylo)
phylo = prune_taxa(taxa_sums(phylo) > 0, phylo)

ps1a <- subset_taxa(phylo,Class!="Chloroplast")
ps2.1 <- subset_taxa(ps1a,Genus!="Mitochondria")


#ps3_SALINE <-subset_samples(ps2.1,Treatment =="SALINE")
#ps3_SALINE_CONTROL <-subset_samples(ps3_SALINE,Description =="Control_experiment")
#ps3_SALINE_CONTROL_T0 <-subset_samples(ps3_SALINE_CONTROL, Time =="T0")

ps3_BMC <-subset_samples(ps2.1,Treatment =="BMC")
ps3_BMC_temp <-subset_samples(ps3_BMC,Description =="Temp_experiment")

df <- ps3_BMC_temp %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 



#write.xlsx(df,"/home/mattoslmp/Desktop/OUTPUT_DADA_QIIME2/Resultado_tabular_qiime2.xlsx")
all <- df %>%
   dplyr:: select(Description,Treatment,Time, Genus, Abundance) %>%
  dplyr:: group_by(Description,Treatment,Time, Genus) %>%
   summarize(
     avg_abundance = mean(Abundance), SD = sd(Abundance)
   )  %>%
   arrange(-avg_abundance)
df2$Genus <- str_replace_all(df2$Genus, "Vibrio", "Pseudovibrio")
df2 <-data.frame(as.character(all$Description),as.character(all$Treatment),as.character(all$Time),
                   as.character(all$Genus),as.numeric(all$avg_abundance),as.numeric(all$SD), stringsAsFactors = FALSE)
colnames(df2) <- c("Description","Treatment","Time", "Genus", "avg_abundance","SD")
 as.factor(df2$Genus)
levels(df2$Genus) <- c(levels(df2$Genus), "Other") 

df2$Genus[df2$avg_abundance<=0.01]  <- "Others <1% Abundance" 
df2$Genus <- str_replace_all(df2$Genus, "Vibrio", "Pseudovibrio")
##45B344 #  #834399 ##6A71D4 #0000CC
Genus_colors <- c("#0F022F","#6B474C","#C2E66E", "#FF0FFA", "#5744B3", "#B3443E", "#94751E",
                  "#6864E6", "#44A7B3", "#993A6C", "#645999","#E67C63", "#569992", "#E6CE63", 
                  "#790000", "#FFFF00", "#E6B645", "#FF2E10", "#00EC00", "#FFA004", "#ED6F9B",
                  "#C27E6E", "#6C2333", "#005300", "#0000CC", "#D4CA6A", "#4F9E5D", "#4A6F87",
                  "#7D4F9E",  "black",    "#3A423A","#00005E", "#77B2D9", "#BE78F0", "#FF1240",
                  "#856E41", "#9C965D", "#77B24B")

gg_title <- expression(paste("\n\n Microbal communities at  Genus Level from " , italic ("Mussismilia hispida"),  " at Temperature Stress in different times"))



#  geom_segment(aes(x=0,xend=1,yend=Time))+
#geom_segment(aes(y=1,yend=0,xend=Time))
g.mid<-ggplot(df2,aes(x=0,y=`Time`))+geom_text(aes(label=`Time`)) +
  geom_segment(aes(x=0,xend=0,yend=`Time`)) +
  geom_segment(aes(x=0,xend=1,yend=`Time`)) +
  ggtitle("") +
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0,0))  +
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1, 0,1,-20), "mm"))

#top, right, bottom, left
g1  <- ggplot(df2)+  geom_col(mapping = aes(x = Time, y = avg_abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  theme(axis.title.x = element_blank())  +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
  ylab("Mean of Relative Abundance (% Phyla ) \n") +
   theme (plot.title =element_text(size=14, face="bold", hjust=0.45))   +
  theme(strip.text.x = element_text(size = 60, colour = "black", face="bold")) +
  theme(strip.text.y = element_text(size = 20, colour = "black", face = "bold")) +
  theme(legend.position = "bottom") +   theme(axis.title.y = element_blank(), axis.text.y = element_text(vjust=0.5,  size = 30, color= "black", face="bold")) +
  theme(legend.title = element_text(color = "black", size = 20),
  legend.text = element_text(color = "black", size = 40, face="bold"))  +
  #theme(plot.margin=unit(c(4,6,4,8),"cm"))  +
  # scale_fill_aaas() +
  scale_fill_manual(values = Genus_colors) +
  theme_minimal()   + coord_flip() + facet_wrap(~Treatment, strip.position="left") +
#  scale_y_reverse() +
    theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-2,1,0), "mm")) +
  theme(legend.position="bottom", legend.margin=margin(0,0,4,-26)) + theme(strip.background = element_rect(fill="grey"))
#top, right, bottom, left
#scale_y_reverse()   facet_wrap(~Time~Description, ncol=8)
#ggsave("phylo.pdf", phylo,height=20, width=40, units='in',limitsize = FALSE, dpi=300)
#scale_y_continuous(expand = c(0,0))
 ##########
ps3_BMC <-subset_samples(ps2.1,Treatment =="SALINE")
ps3_BMC_temp <-subset_samples(ps3_BMC,Description =="Temp_experiment")

df <- ps3_BMC_temp %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 

df$Treatment <- str_replace_all (df$Treatment, "SALINE", "PLACEBO")
#write.xlsx(df,"/home/mattoslmp/Desktop/OUTPUT_DADA_QIIME2/Resultado_tabular_qiime2.xlsx")
all <- df %>%
  dplyr:: select(Description,Treatment,Time, Genus, Abundance) %>%
  dplyr:: group_by(Description,Treatment,Time, Genus) %>%
  summarize(
    avg_abundance = mean(Abundance), SD = sd(Abundance)
  )  %>%
  arrange(-avg_abundance)

df2 <-data.frame(as.character(all$Description),as.character(all$Treatment),as.character(all$Time),
                 as.character(all$Genus),as.numeric(all$avg_abundance),as.numeric(all$SD), stringsAsFactors = FALSE)
colnames(df2) <- c("Description","Treatment","Time", "Genus", "avg_abundance","SD")
as.factor(df2$Genus)
levels(df2$Genus) <- c(levels(df2$Genus), "Other") 
df2$Genus[df2$avg_abundance<=0.01]  <- "Others <1% Abundance" 
df2$Genus <- str_replace_all(df2$Genus, "Vibrio", "Pseudovibrio")

##45B344""Black","#45B344" "#5744B3", "#B3443E"
#"#E3B7AB""#C27E6E", "#4A6F87", "#7D4F9E", "#DC9489" , "#78F08D",
#45B344"
##005300
Genus_colors <- c("#6B474C","#00A1E6","#C2E66E","#FF0FFA","#94751E",  "#6864E6",
                  "#993A6C", "#645999","#E67C63", "#E3B7AB", "#569992", "#E6CE63",
                  "#790000", "#FFFF00", "#C27E6E", "#4A6F87","#E6B645", "#00EC00",
                  "#FFA004", "#EE709C", "#6C2333", "#005300", "#0000CC", "#D4CA6A",
                  "#4F9E5D", "black","#3A423A", "#77B2D9", "#BE78F0",
                  "#FF1240","#856E41", "#9C965D", "#77B24B")

#28       
library(stringr)
df2$Treatment <- str_replace_all (df2$Treatment, "SALINE", "PLACEBO")
g2  <- ggplot(df2)+  geom_col(mapping = aes(x = Time, y = avg_abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  theme(axis.title.x = element_blank())  +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
  ylab("Mean of Relative Abundance (% Phyla ) \n") +
  theme (plot.title =element_text(size=14, face="bold", hjust=0.45))   +
  theme(strip.text.x = element_text(size = 60, colour = "black", face="bold")) +
  theme(strip.text.y = element_text(size = 20, colour = "black", face = "bold")) +
  theme(legend.position = "bottom") +   theme(axis.title.y = element_blank(), axis.text.y = element_text(vjust=0.5,  size = 30, color= "black", face="bold")) +
  theme(legend.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 40, face="bold"))  +
  #theme(plot.margin=unit(c(4,6,4,8),"cm"))  +
  # scale_fill_aaas() +
  scale_fill_manual(values = Genus_colors) +
  theme_minimal()   + coord_flip() + facet_wrap(~Treatment, strip.position="left") +
  #  scale_y_reverse() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-2,1,0), "mm")) +
  theme(legend.position="bottom", legend.margin=margin(0,0,4,-26)) + theme(strip.background = element_rect(fill="grey"))

#top, right, bottom, left
library(gridExtra)
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
#gg.mid <- ggplot_gtable(ggplot_build(g.mid))
fig1 <- grid.arrange(gg1,widths=c(4/9))
ggsave("Temp_Genus_T1_T2_T3_g1.tiff", fig1,height=8, width=15, units='in',limitsize = FALSE, dpi=300)
fig2 <- grid.arrange(gg2,widths=c(4/9))
ggsave("Temp_Genus_T1_T2_T3_g2.tiff", fig2,height=8, width=15, units='in',limitsize = FALSE, dpi=300)



