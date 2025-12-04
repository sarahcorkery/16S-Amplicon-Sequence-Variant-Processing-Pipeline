# 16S-Amplicon-Variant-Sequence-Processing-Pipeline
The following a pipeline processing 16S amplicon variant sequence (ASV) data. The DADA2 portion of this script will allow us to process amplicon sequencing data to identify and quantify ASVs. The latter portion of this script uses PHYLOSEQ to visualize the relative abundance of our processed amplicon sequencing data.

# Installation Instructions 
## dada2 

This pipeline uses a plethora of dada2 functions to process and quantify amplicon sequence variants. Installation instructions can be found at the link below. 

https://benjjneb.github.io/dada2/dada-installation.html

The following script from the installation link was pasted into the console to download dada2. 

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

Installation instructions for phyloseq can be found at the following link.

https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

The following script from the installation link was pasted into the console to download phyloseq. 

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
Be sure to load these package's libraries prior to beginning work. 

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

This step should be completed prior to 16S ASV processing. Please consider the following: 

### Sequencing Depth: 

How many sequencing runs (or depth) do you require per sample? The standard number of sequencing runs per sample is 10,000 or 1x, and tends to be sufficient (Lundin et al. 2012). Should you desire high sequencing diversity, you will likely require more runs. Check relevant existing literature to determine what a good depth for your work may be. This pipeline will use 1x sequencing depth. 

### Primer Bias: 

Different primers tend to work better for different organisms. Choose 16S primers based on the qualities of the organisms your are working with. 

### Amplicon Size:

We will be created paired-end reads, which are forward and reverse reads joined by a small overlap region. The longer these are, the better. Unfortunately, the length of your reads (or ASVs) is limited by your choice of sequencing method.

## Step Two: Quality Control

The input for the 16S ASV Processsing Pipeline and raw output from your sequencing center will be a .fastq file. These files contain your raw sequencing data and the quality scores associated with each base pair. It is not readable via your text editor. 

Most often, excluding those generated using Nanopore, your reads will be paired-end and need to be assembled. Prior to assembling our reads, we must trim our sequences to remove low quality base pairs. We will aim to trim away base pairs below a quality score of 30. This is imperative as poor quality data can yield inaccurate results. 

To check the quality of your forward and reverse reads, the following code is used: 

To visualize the quality profiles of our forward reads:

```{r}
plotQualityProfile(fnFs[1:2])
```
To visualize the quality profiles of our reverse reads:

```{r}
plotQualityProfile(fnRs[1:2])
```

We will trim and filter our reads using the filterAndTrim() function. We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns, the second it can't decipher a nucleotide it'll through the whole sequence out), truncQ=2, rm.phix=TRUE and maxEE=c(2,2) (maxEE allows two errors per read in the forward and reverse reads). The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. 

Be sure to change your truncLen parameter to reflect the lengths you want to truncate your forward and reverse reads to.

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows be sure to set multithread=FALSE!

head(out)
```

## Step Three: Evaluate Error Rates




## Step Three: Merge Paired-End Reads


Align thousands of 16S rRNA genes so that conserved regions line up 
Most often we will use a provided reference file for this purpose built off of the Silva reference database (more on Silva later) 
This is one of the most computational intensive steps ( but there are ways to make it easier*!) 

## Step Four: Evaluate and Identify ASVs 

ASVs = # of different types of microbes, you can set an OTU to be 100% and this would make it an ASV – they are different in the sense that an OTU does have 100% percent identity

## Step Six: Assign Taxonomy

# Manual Steps

# Features:
A list or description of the main functionalities and capabilities of the project.

# Technologies Used:
A list of programming languages, frameworks, libraries, and other tools utilized in the project.

# Known Issues or Limitations:
A section detailing any known bugs, limitations, or areas for improvement.

# License:
Information about the project's licensing, specifying how others can use and distribute the code.

# Contact Information and Acknowledgments:
Details on how to contact the project maintainers or authors.
Acknowledgments for any external contributions, resources, or inspirations.
