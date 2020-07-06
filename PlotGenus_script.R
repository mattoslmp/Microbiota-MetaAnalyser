### Autor: Leandro de Mattos Pereira 12/02/2018
#### Fun??es usadas: Phyloseq, Melt and ggplot, biomformat,dplyr"
rm(list=ls())
############## Se necess?rio, instalar as bibliotecas a seguir com 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("resh")

library("phyloseq")
library("biomformat")
library("dplyr")
library ("ggplot2")
library("ggsci")
library("gridExtra")
library("qiime2R")

setwd("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end")
metadata <- read_q2metadata("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/manifest2.1")
SVs <-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_man2.qza")
taxonomy<-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_taxonomy_man2.qza")
tree<-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/rooted-tree_cp2.qza")
#phylo <- qza_to_phyloseq(features="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_table_pubdata.qza", tree="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_rep_seqs_masked.rooted-tree.qza", taxonomy = "/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_pubdata_taxonomy.qza", metada="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/manifest3")
phylo <- qza_to_phyloseq(features="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_table_man2.qza", taxonomy = "/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_taxonomy_man2.qza", metada="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/manifest2.1", tree="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/rooted-tree_cp2.qza")

#######################################################
ps <- phylo 
rank_names(ps)
taxtable = tax_table(ps)
#ps = prune_taxa(taxa_sums(ps) > 0, ps)

######## Grafico a nivel de Filo
ps2.1 <- subset_taxa(ps,Class   != "Chloroplast" )
ps2.1 <- subset_taxa(ps,Kingdom   != "Archaea" )


df <- ps2.1 %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
 transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum) 



#View(df)
all <- df %>%
   dplyr:: select(Sample, Host, env_feature, isolation_source,geo_loc_name_country,Lat.Long,Temperature,Phylum, Abundance) %>%
  dplyr:: group_by(Sample,Host, env_feature, isolation_source,geo_loc_name_country,Temperature, Lat.Long,Phylum) %>%
   summarize(
     avg_abundance = mean(Abundance), SD = sd(Abundance)
   )  %>%
   arrange(-avg_abundance)

df2 <-data.frame(as.character(all$Sample),as.character(all$Host), as.character(all$env_feature), as.character(all$isolation_source),as.character(all$geo_loc_name_country), as.character(all$Lat.Long),
                 as.character (all$Temperature), as.character(all$Phylum), as.numeric(all$avg_abundance), stringsAsFactors = FALSE)
colnames(df2) <- c("Sample", "Host", "Env_feature", "Isolation Source","Geo Localization name country", "Lat-Long","Temperature","Phylum", "avg_abundance")
as.factor(df2$Phylum)

levels(df2$Phylum) <- c(levels(df2$Phylum), "Other") 

N = 50
GPN = prune_taxa(names(sort(taxa_sums(ps2.1), decreasing = TRUE)[1:N]), 
                 ps2.1)
p <- plot_tree(ps2.1, nodelabf = nodeplotblank, ladderize="left", color="Area", shape = "Host", label.tips = "Genus", text.size = 2.5)



df2$Phylum[df2$avg_abundance<=0.05]  <- "Others <5% Abundance" 
library(stringr)
df2$Env_feature<- str_replace_all(df2$Env_feature, "Shallow Sea", "Shallow Waters")


#factor <- c("S1", "S2", "S3", "S4","S5","S6","S7","S8","S9","S10",
 #           "S11","S12","S13","S14","S15","S16","S17","S18",
  #          "S19","S20","S21","S22", "S23", "S24", "S25", "S26",
   #         "S27", "S28", "S29","S30", "S31", "S32", "S33", "S34",
  #          "S35", "S36", "S37", "S38", "S39", "S40", "S41","S42",
 #           "S43","S44","S45","S46","S47","S48","S49","S50","S51", "S52")
#bar_order <- factor(df2$Sample, levels = factor) 



#phylum_colors <-c ("black", "#660000", "#994C00", "#FF3399", "#4C0099","#006633", "#CCFFCC", "#7F00FF", "#99004C","#808080", "#000080", "#7CFC00","orange","#89C5DA","#D1A33D", "#000000", "#3F4921", "yellow", "#228B22", "#0000CD", 
 #                  "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285", "red","#6DDE88","#808000", "#F08080", "#008B8B", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#FF7F50","#D1A33D", 
  #                 "#8B008B",	"#FFD700", "#FF1493", "#BDB76B","#696969","#FF4500", "#4169E1", "#8FBC8F","#E0FFFF", "#808000","#331900")
#phylum_colors <-c ("#000080", "#7CFC00","orange","#89C5DA","#D1A33D", "#000000", "#3F4921", "yellow", "#228B22", "#0000CD", 
 #                  "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285", "red","#6DDE88","#808000", "#F08080", "#008B8B", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#FF7F50","#D1A33D", 
  #                 "#8B008B",	"#FFD700", "#FF1493", "#BDB76B","#696969","#FF4500", "#4169E1", "#8FBC8F","#E0FFFF", "#808000")
phylum_colors <-c ("#45B344", 	"#708070", "#C2E66E", "#C6AAE9", "#5744B3", "#B3443E", "#B3A039", "#6864E6", "#44A7B3",
  "#993A6C", "#645999", "#E67C63", "#569992", "#E6CE63", "#834399", "#E6E63C", "#E6B645", "#E3B7AB",
  "#E39C6B", "#767CE6", "#628696", "#660000","#994C00","#FF3399","#4C0099","#006633","#CCFFCC","#7F00FF","#99004C","#808080", "#000080","#7CFC00",
  "#89C5DA","#D1A33D", "#000000", "#3F4921", "yellow", "#228B22", "#0000CD", "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285",  "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285", "red","#6DDE88","#808000", "#F08080", "#008B8B", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#FF7F50","#D1A33D", 
  "#8B008B")

gg_title <- expression(paste(""))
library(openxlsx)
#write.xlsx(df2, "phylo_avg_abundance_relative_deep_sea.xlsx")

phylo1  <- ggplot(df2)+  geom_col(mapping = aes(x = `Host`, y = avg_abundance, fill = Phylum), position = "fill", show.legend = TRUE) +
  theme(axis.title.x = element_blank()) +  scale_y_continuous(expand = c(0,0)) + facet_grid("~Env_feature", scales = "free" ) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
  ylab("Mean of Relative Abundance (% Phyla) \n") +
  ggtitle(gg_title) +
#  theme (plot.title =element_text(size=28, face="bold", hjust=0.45))   +
  theme(strip.text.x = element_text(size = 12, colour = "black", face="bold")) +
  theme(strip.text.y = element_text(size = 12, colour = "black", face = "bold")) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(vjust=0.5,  size = 8, color= "black", face="bold")) +
    theme(plot.margin=unit(c(4,6,4,8),"cm"))  +
 # scale_fill_aaas() +
  scale_fill_manual(values = phylum_colors)+
  theme_minimal() +  
  #theme(plot.margin=unit(c(1,11,2,16),"cm")) +
  theme(axis.title.y = element_text(face = "bold", size=12),
        axis.title.x = element_blank())  +
  theme(axis.text.x = element_text(color = "black", face = "bold.italic", size = 14, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(color = "black", face = "bold", size = 16)) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 12, face="bold"))
ggsave("/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/Figuras_paper/Figura_Deap_sea_phylo.png",height=12, width=15, units='in', dpi=300)

############################
rm(list=ls())
### Autor: Leandro de Mattos Pereira 12/02/2018
#### Fun??es usadas: Phyloseq, Melt and ggplot, biomformat,dplyr"
rm(list=ls())
############## Se necess?rio, instalar as bibliotecas a seguir com 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("resh")

library("phyloseq")
library("biomformat")
library("dplyr")
library ("ggplot2")
library("ggsci")
library("gridExtra")
library("qiime2R")

setwd("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end")
metadata <- read_q2metadata("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/manifest2.1")
SVs <-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_man2.qza")
taxonomy<-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_taxonomy_man2.qza")
tree<-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/rooted-tree_cp2.qza")
#phylo <- qza_to_phyloseq(features="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_table_pubdata.qza", tree="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_rep_seqs_masked.rooted-tree.qza", taxonomy = "/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_pubdata_taxonomy.qza", metada="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/manifest3")
phylo <- qza_to_phyloseq(features="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_table_man2.qza", taxonomy = "/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_taxonomy_man2.qza", metada="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/manifest2.1", tree="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/rooted-tree_cp2.qza")


#######################################################
ps <- phylo 
rank_names(ps)
taxtable = tax_table(ps)
#ps = prune_taxa(taxa_sums(ps) > 0, ps)

######## Grafico a nivel de Genus
######## Grafico a nivel de Genus
ps2.1 <- subset_taxa(ps,Class   != "Chloroplast" )
ps2.1 <- subset_taxa(ps,Kingdom   != "Archaea" )

df <-  ps2.1 %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt() %>%                                         # Melt to long format
  arrange(Genus) 


all <- df %>%
  dplyr:: select(Sample, Host, env_feature, isolation_source,geo_loc_name_country,Lat.Long,Temperature,Genus, Abundance) %>%
  dplyr:: group_by(Sample,Host, env_feature, isolation_source,geo_loc_name_country,Temperature, Lat.Long,Genus) %>%
  summarize(
    avg_abundance = mean(Abundance), SD = sd(Abundance)
  )  %>%
  arrange(-avg_abundance)

df2 <-data.frame(as.character(all$Sample),as.character(all$Host), as.character(all$env_feature), as.character(all$isolation_source),as.character(all$geo_loc_name_country), as.character(all$Lat.Long),
                 as.character (all$Temperature), as.character(all$Genus), as.numeric(all$avg_abundance), stringsAsFactors = FALSE)
colnames(df2) <- c("Sample", "Host", "Env_feature", "Isolation Source","Geo Localization name country", "Lat-Long","Temperature","Genus", "avg_abundance")
as.factor(df2$Genus)

levels(df2$Genus) <- c(levels(df2$Genus), "Other") 



df2$Genus[df2$avg_abundance<=0.05]  <- "Others <5% Abundance" 
library(stringr)
df2$Env_feature<- str_replace_all(df2$Env_feature, "Shallow Sea", "Shallow Waters")


#factor <- c("S1", "S2", "S3", "S4","S5","S6","S7","S8","S9","S10",
#           "S11","S12","S13","S14","S15","S16","S17","S18",
#          "S19","S20","S21","S22", "S23", "S24", "S25", "S26",
#         "S27", "S28", "S29","S30", "S31", "S32", "S33", "S34",
#          "S35", "S36", "S37", "S38", "S39", "S40", "S41","S42",
#           "S43","S44","S45","S46","S47","S48","S49","S50","S51", "S52")
#bar_order <- factor(df2$Sample, levels = factor) 



phylum_colors <-c ( "#45B344", 	"#708070", "#C2E66E", "#C6AAE9", "#5744B3", "#B3443E", "#B3A039", "#6864E6", "#44A7B3",
   "#993A6C", "#645999", "#E67C63", "#569992", "#E6CE63", "#834399", "#E6E63C", "#E6B645", "#E3B7AB",
   "#E39C6B", "#767CE6", "#628696", "#660000","#994C00","#FF3399","#4C0099","#006633","#CCFFCC","#7F00FF","#99004C","#808080", "#000080","#7CFC00",
   "#89C5DA","#D1A33D", "#000000", "#3F4921", "yellow", "#228B22", "#0000CD", "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285",  "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285", "red","#6DDE88","#808000", "#F08080", "#008B8B", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#FF7F50","#D1A33D", 
   "#8B008B",   "#903186", "#DFBCBB", "#335E40", "#9DE18B", "#E477EF", "#A97CC6", "#B4AAB7", "#747C4A" ,"#AADBDC",
   "#D3B36D", "#33BF8F", "#BFF2B3" ,"#D0DA98", "#792DBF", "#EA75BC", "#6F9ABC", "#A5EDC0", "#613772",
   "#C5D3C6","#6DA67E" ,"#73EDB7", "#AD29A9", "#9FAB77", "#80C5B0", "#ECAEE6", "#BBBFE5", "#EAB828",
   ,"#FBC177", "#EEC1E9", "#D1F3C6", "#D6D9EE", "#ABCEB0", "#4CE546", "#3380EE" ,"#D9F06B", "#E29133",
   "#789C42", "#3F6AE9" ,"#E688EC", "#EFC194", "#B59071", "#4E94ED", "#DD5252", "#8CF487", "#98C1D4",
   "#AE30F5", "#BACF6E", "#50E6F6", "#6FD8B4", "#C3A1F0", "#E3F03D", "#99F4AC", "#407194", "#CEF1F2",
   "#BE71B9", "#50D5C1", "#7349B0", "#F2E9A1" ,"#F3DE7E" ,"#5B46D5" ,"#AE8935", "#DA5DC0", "#BB84E8",
   "#D1F19C", "#DEE9ED", "#4F29E8" ,"#E9DA47" ,"#92E0EE", "#E15AE6", "#8E46F5", "#64D985", "#BD46E4",
   "#DDE3C6", "#8975E9", "#60B7F1", "#C49BCE", "#6CC230", "#95F1DC", "#B3E44D", "#EEE8BA", "#50CCF3",
   "#8B95E7", "#E9CEB1", "#B6BB38", "#E72AEA", "#EC9AC8", "#F1F7ED", "#BCF2E7", "#BCDCEE" ,"#BEBC96",
   "#7980E4", "#DB3B94", "#8067F4", "#CFF2D9", "#E9577E", "#E69BE9", "#BF8DB0", "#8FB6F0", "#EBB3C4",
   "#5BF4EC", "#78BECA", "#5FD1D0", "#B65882", "#E98E5F", "#EEE0D1", "#F4DFE7", "#A8AAE9", "#EE3AC2",
   "#7BE7EF" ,"#C6CACD" ,"#BA6BEE", "#E9CCE5", "#9AF828", "#AAB6A4", "#DA8576", "#F35927" ,"#DFBB4B",
   "#6DB46A", "#7C9795", "#737A88", "#51ADCD", "#6A7ABE" ,"#3F8F9B", "#4EF08F", "#608AC0", "#4EF6D7",
   "#ECF186", "#A3F569", "#B8F084", "#9BCC91", "#3EEEAF" ,"#8C7371", "#E07C9B", "#85CDEF", "#3FAF9D",
   "#459FE3", "#66B74B", "#98523B", "#D4B7F0" ,"#F2AFA2", "#9996C3" ,"#C98D93" ,"#98F7F6", "#32BCC9")

gg_title <- expression(paste(""))
library(openxlsx)
#write.xlsx(df2, "phylo_avg_abundance_relative_deep_sea.xlsx")

library(randomcoloR)
#n <- 144
#palette <- distinctColorPalette(n)

Genus  <- ggplot(df2)+  geom_col(mapping = aes(x = `Host`, y = avg_abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  theme(axis.title.x = element_blank()) +  scale_y_continuous(expand = c(0,0)) + facet_grid("~Env_feature", scales = "free" ) +
  theme(strip.text.x = element_text(size = 11, colour = "Black")) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
  ylab("Mean of Relative Abundance (% Phyla) \n") +
  ggtitle(gg_title) +
  #  theme (plot.title =element_text(size=28, face="bold", hjust=0.45))   +
  theme(strip.text.x = element_text(size = 12, colour = "black", face="bold")) +
  theme(strip.text.y = element_text(size = 12, colour = "black", face = "bold")) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(vjust=0.5,  size = 8, color= "black", face="bold")) +
  theme(plot.margin=unit(c(4,6,4,8),"cm"))  +
  # scale_fill_aaas() +
 # scale_fill_manual(values = palette) + 
  scale_fill_manual(values = phylum_colors)+
  theme_minimal() +  
  #theme(plot.margin=unit(c(1,11,2,16),"cm")) +
  theme(axis.title.y = element_text(face = "bold", size=12),
        axis.title.x = element_blank())  +
  theme(axis.text.x = element_text(color = "black", face = "bold.italic", size = 14, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(color = "black", face = "bold", size = 16)) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 12, face="bold"))



ggsave("Figura_Genus_Deap_sea_phylo.png",height=12, width=20, units='in', dpi=300)

##########
setwd("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end")
metadata <- read_q2metadata("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/manifest2.1")
SVs <-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_man2.qza")
taxonomy<-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_taxonomy_man2.qza")
tree<-read_qza("/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/rooted-tree_cp2.qza")
#phylo <- qza_to_phyloseq(features="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_table_pubdata.qza", tree="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_rep_seqs_masked.rooted-tree.qza", taxonomy = "/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/MERGED/merged_pubdata_taxonomy.qza", metada="/media/mattos/de7d744f-68db-41f9-9fed-15fd87338433/LEANDRO/Project1_Metabolical_dependence/mydatasets/Metabarcoding_Coral/qiime2/DeepSea/manifest3")
phylo <- qza_to_phyloseq(features="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_table_man2.qza", taxonomy = "/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/dada2_join_pair_rep-seqs_taxonomy_man2.qza", metada="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/manifest2.1", tree="/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end/rooted-tree_cp2.qza")


#######################################################
ps <- phylo 
rank_names(ps)
taxtable = tax_table(ps)
#ps = prune_taxa(taxa_sums(ps) > 0, ps)

######## Grafico a nivel de Genus
######## Grafico a nivel de Genus
ps2.1 <- subset_taxa(ps,Class   != "Chloroplast" )
ps2.1 <- subset_taxa(ps,Kingdom   != "Archaea" )
df <-  ps2.1 %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family) 

all <- df %>%
  dplyr:: select(Sample, Host, env_feature, isolation_source,geo_loc_name_country,Lat.Long,Temperature,Family, Abundance) %>%
  dplyr:: group_by(Sample,Host, env_feature, isolation_source,geo_loc_name_country,Temperature, Lat.Long,Family) %>%
  summarize(
    avg_abundance = mean(Abundance), SD = sd(Abundance)
  )  %>%
  arrange(-avg_abundance)

df2 <-data.frame(as.character(all$Sample),as.character(all$Host), as.character(all$env_feature), as.character(all$isolation_source),as.character(all$geo_loc_name_country), as.character(all$Lat.Long),
                 as.character (all$Temperature), as.character(all$Family), as.numeric(all$avg_abundance), stringsAsFactors = FALSE)
colnames(df2) <- c("Sample", "Host", "Env_feature", "Isolation Source","Geo Localization name country", "Lat-Long","Temperature","Family", "avg_abundance")
as.factor(df2$Family)

levels(df2$Family) <- c(levels(df2$Family), "Other") 



df2$Genus[df2$avg_abundance<=0.05]  <- "Others <5% Abundance" 
library(stringr)
df2$Env_feature<- str_replace_all(df2$Env_feature, "Shallow Sea", "Shallow Waters")


#factor <- c("S1", "S2", "S3", "S4","S5","S6","S7","S8","S9","S10",
#           "S11","S12","S13","S14","S15","S16","S17","S18",
#          "S19","S20","S21","S22", "S23", "S24", "S25", "S26",
#         "S27", "S28", "S29","S30", "S31", "S32", "S33", "S34",
#          "S35", "S36", "S37", "S38", "S39", "S40", "S41","S42",
#           "S43","S44","S45","S46","S47","S48","S49","S50","S51", "S52")
#bar_order <- factor(df2$Sample, levels = factor) 



phylum_colors <-c ( "#45B344", 	"#708070", "#C2E66E", "#C6AAE9", "#5744B3", "#B3443E", "#B3A039", "#6864E6", "#44A7B3",
                    "#993A6C", "#645999", "#E67C63", "#569992", "#E6CE63", "#834399", "#E6E63C", "#E6B645", "#E3B7AB",
                    "#E39C6B", "#767CE6", "#628696", "#660000","#994C00","#FF3399","#4C0099","#006633","#CCFFCC","#7F00FF","#99004C","#808080", "#000080","#7CFC00",
                    "#89C5DA","#D1A33D", "#000000", "#3F4921", "yellow", "#228B22", "#0000CD", "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285",  "#673770", "#D2691E", "#800000",  "#D7C1B1", "#689030", "#AD6F3B", "#FFE4B5","#CD9BCD","#D14285", "red","#6DDE88","#808000", "#F08080", "#008B8B", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#FF7F50","#D1A33D", 
                    "#8B008B",   "#903186", "#DFBCBB", "#335E40", "#9DE18B", "#E477EF", "#A97CC6", "#B4AAB7", "#747C4A" ,"#AADBDC",
                    "#D3B36D", "#33BF8F", "#BFF2B3" ,"#D0DA98", "#792DBF", "#EA75BC", "#6F9ABC", "#A5EDC0", "#613772",
                    "#C5D3C6","#6DA67E" ,"#73EDB7", "#AD29A9", "#9FAB77", "#80C5B0", "#ECAEE6", "#BBBFE5", "#EAB828",
                    ,"#FBC177", "#EEC1E9", "#D1F3C6", "#D6D9EE", "#ABCEB0", "#4CE546", "#3380EE" ,"#D9F06B", "#E29133",
                    "#789C42", "#3F6AE9" ,"#E688EC", "#EFC194", "#B59071", "#4E94ED", "#DD5252", "#8CF487", "#98C1D4",
                    "#AE30F5", "#BACF6E", "#50E6F6", "#6FD8B4", "#C3A1F0", "#E3F03D", "#99F4AC", "#407194", "#CEF1F2",
                    "#BE71B9", "#50D5C1", "#7349B0", "#F2E9A1" ,"#F3DE7E" ,"#5B46D5" ,"#AE8935", "#DA5DC0", "#BB84E8",
                    "#D1F19C", "#DEE9ED", "#4F29E8" ,"#E9DA47" ,"#92E0EE", "#E15AE6", "#8E46F5", "#64D985", "#BD46E4",
                    "#DDE3C6", "#8975E9", "#60B7F1", "#C49BCE", "#6CC230", "#95F1DC", "#B3E44D", "#EEE8BA", "#50CCF3",
                    "#8B95E7", "#E9CEB1", "#B6BB38", "#E72AEA", "#EC9AC8", "#F1F7ED", "#BCF2E7", "#BCDCEE" ,"#BEBC96",
                    "#7980E4", "#DB3B94", "#8067F4", "#CFF2D9", "#E9577E", "#E69BE9", "#BF8DB0", "#8FB6F0", "#EBB3C4",
                    "#5BF4EC", "#78BECA", "#5FD1D0", "#B65882", "#E98E5F", "#EEE0D1", "#F4DFE7", "#A8AAE9", "#EE3AC2",
                    "#7BE7EF" ,"#C6CACD" ,"#BA6BEE", "#E9CCE5", "#9AF828", "#AAB6A4", "#DA8576", "#F35927" ,"#DFBB4B",
                    "#6DB46A", "#7C9795", "#737A88", "#51ADCD", "#6A7ABE" ,"#3F8F9B", "#4EF08F", "#608AC0", "#4EF6D7",
                    "#ECF186", "#A3F569", "#B8F084", "#9BCC91", "#3EEEAF" ,"#8C7371", "#E07C9B", "#85CDEF", "#3FAF9D",
                    "#459FE3", "#66B74B", "#98523B", "#D4B7F0" ,"#F2AFA2", "#9996C3" ,"#C98D93" ,"#98F7F6", "#32BCC9")

gg_title <- expression(paste(""))
library(openxlsx)
#write.xlsx(df2, "phylo_avg_abundance_relative_deep_sea.xlsx")

library(randomcoloR)
#n <- 144
#palette <- distinctColorPalette(n)

Family  <- ggplot(df2)+  geom_col(mapping = aes(x = `Host`, y = avg_abundance, fill = Family), position = "fill", show.legend = TRUE) +
  theme(axis.title.x = element_blank()) +  scale_y_continuous(expand = c(0,0)) + facet_grid("~Env_feature", scales = "free" ) +
  theme(strip.text.x = element_text(size = 11, colour = "Black")) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
  ylab("Mean of Relative Abundance (% Phyla) \n") +
  ggtitle(gg_title) +
  #  theme (plot.title =element_text(size=28, face="bold", hjust=0.45))   +
  theme(strip.text.x = element_text(size = 12, colour = "black", face="bold")) +
  theme(strip.text.y = element_text(size = 12, colour = "black", face = "bold")) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(vjust=0.5,  size = 8, color= "black", face="bold")) +
  theme(plot.margin=unit(c(4,6,4,8),"cm"))  +
  # scale_fill_aaas() +
  # scale_fill_manual(values = palette) + 
  scale_fill_manual(values = phylum_colors)+
  theme_minimal() +  
  #theme(plot.margin=unit(c(1,11,2,16),"cm")) +
  theme(axis.title.y = element_text(face = "bold", size=12),
        axis.title.x = element_blank())  +
  theme(axis.text.x = element_text(color = "black", face = "bold.italic", size = 14, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(color = "black", face = "bold", size = 16)) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 12, face="bold"))


