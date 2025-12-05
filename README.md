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

## SILVA

To assign taxonomy later in this pipeline, we will be using SILVA. 

Using the link below, we will download the most recent version of SILVA, version 138.2 to our computer.

https://zenodo.org/records/14169026

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

# Usage Instructions - dada2 Processing

The following code will outline the main steps of 16S Amplicon Sequence Variant processing pipeline. 

## Step One: Experimental Design 

This step should be completed prior to 16S ASV processing. Please consider the following: 

### Sequencing Depth: 

How many sequencing runs (or depth) do you require per sample? The standard number of sequencing runs per sample is 10,000 or 1x, and tends to be sufficient (Lundin et al. 2012). Should you desire high sequencing diversity, you will likely require more runs. Check relevant existing literature to determine what a good depth for your work may be. This pipeline will use 1x sequencing depth. 

### Primer Bias: 

Different primers tend to work better for different organisms. Choose 16S primers based on the qualities of the organisms your are working with. 

### Amplicon Size:

We will be created paired-end reads, which are forward and reverse reads joined by a small overlap region. The longer these are, the better. Unfortunately, the length of your reads (or ASVs) is limited by your choice of sequencing method.

## Step Two: Initial Preparation

First, you must saving your raw sequencing files to where you intend to set your working directory. Once they have been moved to your desired location, unzip your fastq.gz files. Ensure they are .fastq files before moving forward. 

Next, we will assign our working directory. Be sure to set your working directory to where your sequence files are located!

```{r}
setwd("/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq")
```
To make the location to these files more easily accessible for downstream commands, be sure to assign your working directory (i.e., where our files are located) to an object called path.

```{r}
path<-"/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq"
```

At this point, we must verify that our forward and reverse fastq files have the naming format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq. Below, we will tell R studio what formatting represents a forward and reverse read file. This is done by assigning forward and reverse reads to the objects fnFs for forward reads, and fnRs for reverse reads. 

```{r}
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
```

Finally, we will make a list of our sample names using our forward reads for downstream naming purposes. This is not done for the reverse read sample names as they will be identical to the forward read names.

```{r}    
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

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

Next, we will estimate error rates from our sequencing data as dada2's algorithm relies on a parametric error model to distinguish true biological sequences from sequencing errors. An error rates measures how likely any transition (i.e., A->C, A->G) occurred during sequencing of your data. dada2 does this through the # learnErrors() function. This will give us higher taxonomic resolution compared to what operational taxonomic units would provide.

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```
To visualize your error rates, use the script below. 

```{r}
plotErrors(errF, nominalQ=TRUE)
```

We can tell if our rates are good or bad based on whether the estimated error rates (red lines) are a good fit to the observed rates (black lines). If there is too much divergence, our data may not be useable. 

Finally, we will apply the core sample inference algorithm to our filtered and trimmed sequence data. With this command we take each sample's reads and use our error rates to infer which sequences are real biological variants and which are sequencing errors. dada() then infers the true error-free ASVs in each sample and how many reads belong to each.
 
```{r}
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
```

## Step Four: Merge Paired-End Reads

For our next step, we will merge our denoised forward and reverse reads to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, to construct our merged “contig” sequences. 

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

## Step Five: Evaluate and Identify ASVs 

An ASV is a unique, error-corrected DNA sequences from amplified marker genes that represent true biological variation. Amplicon sequence variants (ASVs) are grouped only when they share 100% sequence identity, offering greater precision than operational taxonomic units (OTUs), which require only a 97% similarity threshold.

Let's start step five with constructing an ASV table, a higher-resolution version of the OTU table produced by traditional methods.

``{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

The following code will allow us to inspect the distribution of sequence lengths.

```{r}
# From here we can inspect the distribution of sequence lengths. 
table(nchar(getSequences(seqtab)))

```

While the dada() function corrects substitutions, indel errors and chimeras remain. Fortunately, the accuracy of sequence variants after denoising makes identifying chimeric ASVs simpler than when dealing with OTUs. The following script will construct and remove chimeras from our collection of sequences. 

```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

Finally, we will divide the table containing no chimeras with the table containing chimeras, to get a percentage of chimeric sequences in our original list of ASVs. 

```{r}
sum(seqtab.nochim)/sum(seqtab)

As a final check of our progress, we will look at the number of reads that made it through each step in the pipeline using the commands below:

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

## Step Six: Assign Taxonomy

For step six, we will assign taxonomy to our sequences using the assignTaxonomy() function. This function requires appropriately formatted FASTA files containing taxonomically classified reference sequences to use as a training dataset. To facilitate this process, we will download the SILVA taxonomy database, which provides the necessary reference sequences for accurate taxonomic assignment.

We must state the path to this file on our computer in the assignTaxonomy() function below. Using the command below, we tell DADA2 to assign taxonomy to our seqtab.nochim matrix, which contains our non-chimeric ASVs and their abundances.

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/sarahcorkery/Desktop/6452 02 - Bioinformatics/R files/Silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
```

dada2 only makes species level assignments based on exact matching between ASVs and sequenced reference strains. Below, the addSpecies() function will assign species-level annotation to our taxonomic table.

```{r} 
taxa <- addSpecies(taxa, "/Users/sarahcorkery/Desktop/6452 02 - Bioinformatics/R files/Silva/silva_v138.2_assignSpecies.fa.gz")
```

Lastly, we will inspect the taxonomic assignmnet by removing sequence rownames for display purposes only

```{r}
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

## Step Seven: Save your Progress

Finally, we will save our taxa object, containing taxonomic assignments, our seqtab.nochim object, containing our ASVs and their abundances, and our track object, showing the remaining reads following each step of the DADA2 pipeline as .csv files. This step will allow us to save our progress, and later read in data, in the event our objects are lost or we terminate our R Studio session. 

write.csv(taxa, file = "/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5_taxa.csv")

write.csv(seqtab.nochim, file = "/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5_seqtabnochim.csv")

write.csv(track, file = "/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5_track.csv")

Using the commands below, we can read in the taxa, seqtab.nochim and track objects using the .csv files we made in the previous step. This is essential in the event our objects are lost or we terminate our R studio session:

taxa <- read.csv(file = "/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5_taxa.csv")
seqtab.nochim <- read.csv(file = "/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5_seqtabnochim.csv", header = FALSE)
track <- read.csv(file = "/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5_track.csv")

# Usage Instructions - phyloseq Processing

# In the next set of chunks, we will change the formatting of our seqtab.nochim table. This is necessary for downstream work with phyloseq, as well as an optional step to merge seqtab.nochim with our taxa table for easier viewing.

```{r}
# Step 1: To begin this process, we will first transpose our seqtab.nochim table. You must use the imported seqtab.nochim object (fromm the .csv) in the previous step for this step to run properly. 
flipped_seqtab.nochim <- as.data.frame(t(seqtab.nochim))
flipped_seqtab.nochim
```


# In this flipped csv, the header becomes a sequence (ASV) and moves the sample names to the first row (not the header).

```{r}
# Step 2: Next, we will copy the first row. 
colnames(flipped_seqtab.nochim) <- flipped_seqtab.nochim[1,]
# The command below allows us to verify that our formatting is correct thus far. It should look the same as it did following Step 1, apart from the header now stating the correct sample names. 
View(flipped_seqtab.nochim)
```

```{r}
# Step 3: Next, we will delete the first row of our transposed seqtab.nochim table. 
flipped_seqtab.nochim <- flipped_seqtab.nochim[-1,]
# The command below allows us to once again verify that our formatting is correct. The first row (not the header) should have been deleted. 
View(flipped_seqtab.nochim)
```

# Next we want to change the names of the sequences to "ASV(n)" so it is more digestable than the nucleotide sequence itself. 

```{r}
# Step 4: In the command below, we will take each column and its corresponding row name, and paste "ASV"1,2,3.. as a new column next to our nucleotide sequence column.  
rownames(flipped_seqtab.nochim) <- paste0("ASV", 1:nrow(flipped_seqtab.nochim))
```

```{r}
# Step 5: Next, we'll remove the nucleotide sequences column and save this new table as flipped_seqtab.nochim_forself for our ease of viewing.  
flipped_seqtab.nochim_forself <- flipped_seqtab.nochim[,-1]
```

```{r}
# Step 6: We will then save these two tables (the one with our nucleotide sequences and the other without) as .csv files. In the event our R Studio session is lost, we can always read these files back in. 
write.csv(flipped_seqtab.nochim, file = '/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5flipped_seqtab.nochim.csv')
write.csv(flipped_seqtab.nochim_forself, file ='/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5flipped_seqtab.nochim_forself.csv')
```

## The following steps are placed outside of a chunk as they are intended to be adjusted prior to running this script.

Manual Step Here: 

# Step 7: Next, as an optional step, we can use the cbind() function to save our flipped seqtab.nochim and taxa data as one csv. This approach offers the most comprehensive view of our ASV abundances and taxonomic classifications.

OTUabund<-cbind(flipped_seqtab.nochim,taxa)
write.csv(OTUabund,file='/Users/sarahcorkery/Desktop/ASV Processing/CorkeryV4V5/new_fastq/output_files/CorkeryV4V5OTUabund.csv')

```{r}
# Step 8: DADA2 generates the taxa object with sequences in the first column, a format incompatible with PHYLOSEQ. To ensure compatibility, we will remove this column using the command below. 
taxa <-taxa[-1]
# The command below allows us to verify that our formatting is correct. The first column containing nucleotide sequences should have been deleted.
View(taxa)

## Step One: Formatting for phyloseq


# Features:
A list or description of the main functionalities and capabilities of the project.

# Technologies Used:
A list of programming languages, frameworks, libraries, and other tools utilized in the project.

# Known Issues or Limitations:
A section detailing any known bugs, limitations, or areas for improvement.

- downloading packages

# License:
Information about the project's licensing, specifying how others can use and distribute the code.

# Contact Information and Acknowledgments:
Details on how to contact the project maintainers or authors.
Acknowledgments for any external contributions, resources, or inspirations.
