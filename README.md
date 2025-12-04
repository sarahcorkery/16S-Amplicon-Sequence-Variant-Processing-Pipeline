# 16S-Amplicon-Variant-Sequence-Processing-Pipeline
The following a pipeline processing 16S amplicon variant sequence (ASV) data. The DADA2 portion of this script will allow us to process amplicon sequencing data to identify and quantify ASVs. The latter portion of this script uses PHYLOSEQ to visualize the relative abundance of our processed amplicon sequencing data.

# Installation Instructions 
## dada2 

This pipeline uses a plethora of dada2 functions to process and quantify amplicon sequence variants. Installation instructions can be found here: 

https://benjjneb.github.io/dada2/dada-installation.html

The following script from the installation link was pasted into the console to download dada2 for this pipeline. 

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2")
```
Be sure to load dada2's library prior to beginning work. 

```{r}
library(dada2); packageVersion("dada2")
```
## phyloseq

Phyloseq and a set of other packages are required for this processing pipeline. These will allow us to produce high-quality graphs from our processed amplicon sequencing data. 

Installation instructions for phyloseq can be found at the following link: 

https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

The following script from the installation link was pasted into the console to download phyloseq for this pipeline. 

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
```
The remaining packages required for this portion of the pipeline can be installed in the console using the script below. 

```{r}
install.packages("Biostrings")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("tidyverse")
```
Be sure to load these package's libraries prior to beginning the pipeline. 

```{r}
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
```

# Usage Instructions

The following code will outline the main steps of 16S Amplicon Sequence Variant processing pipeline. 

## Step One: Experimental Design 

Examples or demonstrations of how to use the project's features and functionalities.
Code snippets or commands illustrating common use cases.

## Step Two: Quality Control + processing raw reads 
## Step Three: Alignment
## Step Four: Identify ASVs 

ASVs = # of different types of microbes, you can set an OTU to be 100% and this would make it an ASV – they are different in the sense that an OTU does have 100% percent identity

## Step Five: Diversity Indices
## Step Six: Classification

# Manual Steps

## Features:
A list or description of the main functionalities and capabilities of the project.

## Technologies Used:
A list of programming languages, frameworks, libraries, and other tools utilized in the project.

## Known Issues or Limitations:
A section detailing any known bugs, limitations, or areas for improvement.

## License:
Information about the project's licensing, specifying how others can use and distribute the code.

## Contact Information and Acknowledgments:
Details on how to contact the project maintainers or authors.
Acknowledgments for any external contributions, resources, or inspirations.
