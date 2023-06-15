# Readme

This repository contains R code for analyzing 16S rRNA sequencing data using the Phyloseq, ggplot2, biomformat, and dplyr packages. The code processes the sequencing data, performs taxonomic classification, and generates visualizations.

## Installation

To run the code, you need to have R installed on your system. Additionally, the following R packages should be installed:

- Phyloseq
- ggplot2
- biomformat
- dplyr
- ggsci
- gridExtra
- qiime2R

You can install these packages using the following commands:

```R
install.packages("phyloseq")
install.packages("ggplot2")
install.packages("biomformat")
install.packages("dplyr")
install.packages("ggsci")
install.packages("gridExtra")
install.packages("qiime2R")

## Usage

- Clone or download the repository to your local machine.

- Set the working directory in the R code to the appropriate location where your input data files are located.

- Modify the file paths in the R code to match the locations of your input data files.

- Run the R code in an R environment or an R script editor.

- Note: Make sure to install the required R packages mentioned in the installation section.

## Description

The R code performs the following steps:

Reads the metadata file and input data files (SVs, taxonomy, and tree) using the read_qza and read_q2metadata functions from the qiime2R package.

Creates a Phyloseq object using the input data.

Filters the data, transforms it to relative abundance, and aggregates at the phylum level.

Calculates the mean abundance and standard deviation for each phylum in each sample.

Sorts the data by average abundance.

Creates a bar plot of the average abundance of phyla, grouped by host, environment feature, and isolation source.

Saves the plot as a PNG file.

Note: The code assumes specific file formats and data structure. Make sure to adapt it to your specific dataset.

## License
This code is released under the MIT License.

Please feel free to modify and use this code according to your needs.

## Contact
If you have any questions or suggestions, please feel free to contact the author:

Author: Leandro de Mattos Pereira
Email: mattoslmp@gmail.com
