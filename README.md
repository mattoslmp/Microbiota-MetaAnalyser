# Readme
# Análise Metagenômica 16S

Este repositório contém um script em R para análise de dados de metagenômica 16S utilizando pacotes como Phyloseq, ggplot2 e dplyr. O script realiza a leitura de arquivos gerados pelo QIIME2 e faz a agregação e visualização de dados taxonômicos em diferentes níveis, como Filo, Gênero e Família.

## Funcionalidades

- **Carregamento de arquivos QIIME2**: O script aceita arquivos `.qza` gerados pelo QIIME2, como sequências, taxonomia e árvore filogenética.
- **Análise de Abundância Relativa**: Calcula a abundância relativa em níveis taxonômicos diferentes (Phylum, Genus, Family) e gera gráficos de barras empilhadas.
- **Flexibilidade com parâmetros**: O script é modular e pode ser reutilizado para diferentes datasets, bastando fornecer os nomes dos arquivos.
- **Visualizações customizadas**: Utiliza `ggplot2` para visualização dos resultados em gráficos com cores definidas para diferentes grupos taxonômicos.

## Estrutura do Projeto

- **`phylo_analysis.R`**: O script principal que realiza a análise de abundância relativa e gera os gráficos.
- **`data/`**: Diretório para armazenar os arquivos de dados utilizados (não incluído no repositório).
- **`figures/`**: Diretório para armazenar as figuras geradas pelos gráficos (não incluído no repositório).

## Requisitos

O script requer os seguintes pacotes R:

- `phyloseq`
- `biomformat`
- `dplyr`
- `ggplot2`
- `ggsci`
- `gridExtra`
- `qiime2R`
- `stringr`
- `openxlsx`

Para instalar todos os pacotes, execute o seguinte código no R:

```r
install.packages(c("phyloseq", "biomformat", "dplyr", "ggplot2", "ggsci", "gridExtra", "qiime2R", "stringr", "openxlsx"))
```
## Uso
## 1. Configuração do ambiente de trabalho
 Certifique-se de que todos os arquivos de dados estejam localizados no diretório correto. No início do script, você pode definir o diretório de trabalho e os arquivos de entrada:
diretorio <- "/caminho/para/seus/dados"
metadata_file <- "manifesto.qza"
features_file <- "sequencias.qza"
taxonomy_file <- "taxonomia.qza"
tree_file <- "arvore.qza"
## 2. Execução do Script

- O script é dividido em funções modulares que permitem carregar dados, processá-los e gerar gráficos de abundância relativa em diferentes níveis taxonômicos (Filo, Gênero e Família). Para rodar o script, basta:
source("phylo_analysis.R")

## 3. Geração de Gráficos

Os gráficos são gerados automaticamente e salvos no diretório especificado. Eles incluem:

- Abundância relativa por Filo (phylo_avg_abundance_relative_deep_sea.png)
- Abundância relativa por Gênero (genus_avg_abundance_relative_deep_sea.png)
- Abundância relativa por Família (family_avg_abundance_relative_deep_sea.png)


# Carregar o script
source("phylo_analysis.R")

# Definir os arquivos
diretorio <- "/home/user/metagenomics"
metadata_file <- "metadata.qza"
features_file <- "features.qza"
taxonomy_file <- "taxonomy.qza"
tree_file <- "tree.qza"

# Carregar os dados
dados <- carregar_dados_qiime2(diretorio, metadata_file, features_file, taxonomy_file, tree_file)

# Criar o objeto Phyloseq
ps <- criar_phyloseq(dados, features_file, taxonomy_file, metadata_file, tree_file)

# Gerar gráficos
gerar_grafico_abundancia(df2, "Host", "avg_abundance", "Phylum", "Env_feature", "Mean of Relative Abundance (% Phyla)", "phylo_abund.png")

## Script done for taxonomic profile analyes of coral microbiome:

- Coral Microbiome Manipulation Elicits metabolic and genetic restructuring to mitigate heat stress and evade Mortality. Santoro, E. P.; Borges, R. M.; Espinoza, J. L.; Freire., M.; Messias, C. S. M. A.; Villela, H. M. D.; Mattos, L. P.; Vilela, C. L. S.; Rosado, J. G.; Cardoso, P. M.; Rosado, P. M.; Assis, J. M.; Duarte, G. A. S.; Perna, G.; Rosado, A. S.; Macrae, A.; Dupont, C. L.; Nelson, K.E.; Sweet, M. J.; Voolstra, C. R.; Peixoto, R. S. Novembro de 2020. Science Advance. v. 7, p. eabg3088, 2021. Link to the Paper
