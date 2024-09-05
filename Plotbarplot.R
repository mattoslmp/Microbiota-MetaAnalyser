# Data: 12/02/2018
# Descrição: Script para análise de dados de metagenômica 16S utilizando Phyloseq, Melt, ggplot e dplyr.
#            Agora com variáveis dinâmicas para facilitar a reutilização e melhorias na estrutura.

# Função para carregar pacotes automaticamente
carregar_pacotes <- function(pacotes) {
  for (pkg in pacotes) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Função para carregar dados com base em variáveis de nomes
carregar_dados_qiime2 <- function(diretorio, metadata_file, features_file, taxonomy_file, tree_file) {
  setwd(diretorio)  # Definir diretório de trabalho
  
  metadata <- read_q2metadata(metadata_file)  # Carregar metadados
  features <- read_qza(features_file)  # Carregar sequência de features
  taxonomy <- read_qza(taxonomy_file)  # Carregar dados de taxonomia
  tree <- read_qza(tree_file)  # Carregar a árvore filogenética
  
  # Retornar uma lista com todos os dados carregados
  return(list(metadata = metadata, features = features, taxonomy = taxonomy, tree = tree))
}

# Função para criar o objeto Phyloseq com dados dinâmicos
criar_phyloseq <- function(dados, features_file, taxonomy_file, metadata_file, tree_file) {
  phylo <- qza_to_phyloseq(
    features = features_file, 
    taxonomy = taxonomy_file, 
    metada = metadata_file, 
    tree = tree_file
  )
  return(phylo)
}

# Função para gerar gráficos com parâmetros dinâmicos
gerar_grafico_abundancia <- function(df, x_var, y_var, fill_var, facet_var, title, filename) {
  ggplot(df) +  
    geom_col(mapping = aes_string(x = x_var, y = y_var, fill = fill_var), position = "fill", show.legend = TRUE) +
    scale_y_continuous(expand = c(0,0)) + 
    facet_grid(paste0("~", facet_var), scales = "free") + 
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) + 
    ggtitle(title) +
    theme_minimal() +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 14)) +
    theme(legend.position = "bottom") +
    ggsave(filename, height = 12, width = 15, units = 'in', dpi = 300)
}

# Carregar pacotes necessários
pacotes_necessarios <- c("phyloseq", "biomformat", "dplyr", "ggplot2", "ggsci", "gridExtra", "qiime2R", "stringr", "openxlsx")
carregar_pacotes(pacotes_necessarios)

# Definir variáveis dinâmicas
diretorio <- "/home/mattoslmp/Desktop/METADADO_META_ANALISE/Resultados_16S/Resultados_pair_end"
metadata_file <- "manifest2.1"
features_file <- "dada2_join_pair_rep-seqs_man2.qza"
taxonomy_file <- "dada2_join_pair_rep-seqs_taxonomy_man2.qza"
tree_file <- "rooted-tree_cp2.qza"

# Carregar os dados
dados <- carregar_dados_qiime2(diretorio, metadata_file, features_file, taxonomy_file, tree_file)

# Criar o objeto Phyloseq
ps <- criar_phyloseq(dados, features_file, taxonomy_file, metadata_file, tree_file)

# Manipulação de dados e criação do dataframe para o gráfico
ps2.1 <- subset_taxa(ps, Class != "Chloroplast")
ps2.1 <- subset_taxa(ps, Kingdom != "Archaea")

df <- ps2.1 %>%
  tax_glom(taxrank = "Phylum") %>%  # Agregar por nível de filo
  transform_sample_counts(function(x) { x / sum(x) }) %>%  # Transformar em abundância relativa
  psmelt() %>%
  arrange(Phylum)

# Resumo dos dados e cálculo de média de abundância
all <- df %>%
  dplyr::select(Sample, Host, env_feature, isolation_source, geo_loc_name_country, Lat.Long, Temperature, Phylum, Abundance) %>%
  dplyr::group_by(Sample, Host, env_feature, isolation_source, geo_loc_name_country, Temperature, Lat.Long, Phylum) %>%
  summarize(avg_abundance = mean(Abundance), SD = sd(Abundance)) %>%
  arrange(-avg_abundance)

# Ajustar dataframe para o gráfico
df2 <- data.frame(as.character(all$Sample), as.character(all$Host), as.character(all$env_feature), as.character(all$isolation_source),
                  as.character(all$geo_loc_name_country), as.character(all$Lat.Long), as.character(all$Temperature),
                  as.character(all$Phylum), as.numeric(all$avg_abundance), stringsAsFactors = FALSE)

colnames(df2) <- c("Sample", "Host", "Env_feature", "Isolation Source", "Geo Localization name country", "Lat-Long", "Temperature", "Phylum", "avg_abundance")
df2$Phylum[df2$avg_abundance <= 0.05] <- "Others <5% Abundance"
df2$Env_feature <- str_replace_all(df2$Env_feature, "Shallow Sea", "Shallow Waters")

# Definir cores para o gráfico
phylum_colors <- c("#45B344", "#708070", "#C2E66E", "#C6AAE9", "#5744B3", "#B3443E", "#B3A039", "#6864E6", "#44A7B3")

# Gerar gráfico a nível de Phylum
gerar_grafico_abundancia(
  df = df2, 
  x_var = "Host", 
  y_var = "avg_abundance", 
  fill_var = "Phylum", 
  facet_var = "Env_feature", 
  title = "Mean of Relative Abundance (% Phyla)",
  filename = "phylo_avg_abundance_relative_deep_sea.png"
)

# Repetir para nível de Genus
df_genus <- ps2.1 %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) { x / sum(x) }) %>%
  psmelt() %>%
  arrange(Genus)

all_genus <- df_genus %>%
  dplyr::select(Sample, Host, env_feature, isolation_source, geo_loc_name_country, Lat.Long, Temperature, Genus, Abundance) %>%
  dplyr::group_by(Sample, Host, env_feature, isolation_source, geo_loc_name_country, Temperature, Lat.Long, Genus) %>%
  summarize(avg_abundance = mean(Abundance), SD = sd(Abundance)) %>%
  arrange(-avg_abundance)

df2_genus <- data.frame(as.character(all_genus$Sample), as.character(all_genus$Host), as.character(all_genus$env_feature),
                        as.character(all_genus$isolation_source), as.character(all_genus$geo_loc_name_country),
                        as.character(all_genus$Lat.Long), as.character(all_genus$Temperature),
                        as.character(all_genus$Genus), as.numeric(all_genus$avg_abundance), stringsAsFactors = FALSE)

colnames(df2_genus) <- c("Sample", "Host", "Env_feature", "Isolation Source", "Geo Localization name country", "Lat-Long", "Temperature", "Genus", "avg_abundance")
df2_genus$Genus[df2_genus$avg_abundance <= 0.05] <- "Others <5% Abundance"
df2_genus$Env_feature <- str_replace_all(df2_genus$Env_feature, "Shallow Sea", "Shallow Waters")

# Gerar gráfico a nível de Genus
gerar_grafico_abundancia(
  df = df2_genus, 
  x_var = "Host", 
  y_var = "avg_abundance", 
  fill_var = "Genus", 
  facet_var = "Env_feature", 
  title = "Mean of Relative Abundance (% Genus)",
  filename = "genus_avg_abundance_relative_deep_sea.png"
)

# Repetir para nível de Family
df_family <- ps2.1 %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) { x / sum(x) }) %>%
  psmelt() %>%
  arrange(Family)

all_family <- df_family %>%
  dplyr::select(Sample, Host, env_feature, isolation_source, geo_loc_name_country, Lat.Long, Temperature, Family, Abundance) %>%
  dplyr::group_by(Sample, Host, env_feature, isolation_source, geo_loc_name_country, Temperature, Lat.Long, Family) %>%
  summarize(avg_abundance = mean(Abundance), SD = sd(Abundance)) %>%
  arrange(-avg_abundance)

df2_family <- data.frame(as.character(all_family$Sample), as.character(all_family$Host), as.character(all_family$env_feature),
                         as.character(all_family$isolation_source), as.character(all_family$geo_loc_name_country),
                         as.character(all_family$Lat.Long), as.character(all_family$Temperature),
                         as.character(all_family$Family), as.numeric(all_family$avg_abundance), stringsAsFactors = FALSE)

colnames(df2_family) <- c("Sample", "Host", "Env_feature", "Isolation Source", "Geo Localization name country", "Lat-Long", "Temperature", "Family", "avg_abundance")
df2_family$Family[df2_family$avg_abundance <= 0.05] <- "Others <5% Abundance"
df2_family$Env_feature <- str_replace_all(df2_family$Env_feature, "Shallow Sea", "Shallow Waters")

# Gerar gráfico a nível de Family
gerar_grafico_abundancia(
  df = df2_family, 
  x_var = "Host", 
  y_var = "avg_abundance", 
  fill_var = "Family", 
  facet_var = "Env_feature", 
  title = "Mean of Relative Abundance (% Family)",
  filename = "family_avg_abundance_relative_deep_sea.png"
)

# Finalizar e limpar o ambiente
rm(list = ls())



