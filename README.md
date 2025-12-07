# 16S-Amplicon-Sequence-Variant-Processing-Pipeline
This README file provides an overview of the features and workflow steps included in the R Studio 16S amplicon sequence variant (ASV) processing pipeline script. The first section of the script utilizes dada2 to identify and quantify ASVs from 16S amplicon sequencing data. The subsequent section uses phyloseq to visualize the relative abundance of the processed amplicon sequencing data, facilitating effective interpretation of microbial community composition.

# Installation Instructions 
## dada2 

This pipeline uses a plethora of dada2 functions to process and quantify amplicon sequence variants. Installation instructions can be found at the link below. 

https://benjjneb.github.io/dada2/dada-installation.html

The following script from the installation link can be pasted into R Studio's console to download dada2.  

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2")
```
Be sure to load dada2's library into R Studio's environment prior to beginning work. 

```{r}
library(dada2); packageVersion("dada2")
```

## SILVA

To assign taxonomy later in the pipeline, we will be using SILVA. 

Using the link below, we will download the most recent version of SILVA, version 138.2 to our computer.

https://zenodo.org/records/14169026

Note its location on your computer as it will be necessary later in the pipeline. 

## phyloseq

Phyloseq and a set of other packages are required for this processing pipeline. These will allow us to produce high-quality graphs from our processed amplicon sequencing data. 

Installation instructions for phyloseq can be found at the following link.

https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

The following script from the installation link can be pasted into R Studio's console to download phyloseq. 

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
```
The remaining packages required for this portion of the pipeline can be installed by pasting the following script into R Studio's console.

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

The following code will outline the main steps of the R Studio 16S Amplicon Sequence Variant processing pipeline. 

## Step One: Experimental Design 

This step should be completed prior to 16S ASV processing. Please consider the following: 

### Sequencing Depth: 

How many sequencing runs (or what depth) do you require per sample? The standard number of sequencing runs per sample is 10,000 or 1x depth, which tends to be sufficient (Lundin et al. 2012). You will likely need more runs if you desire high sequencing diversity. Examine relevant literature to determine what a good depth for your work may be. 

### Primer Bias: 

Different primers work better with different organisms. Choose 16S primers based on the qualities of the bacteria (or archaea) you are working with. 

### Amplicon Size:

We will be creating paired-end reads, which are forward and reverse reads joined by a small overlap region. The longer these are, the better. Pay mind, as the length of your reads (or ASVs) is limited by your choice of sequencing method.

## Step Two: Initial Preparation

First, you must save your raw sequencing files to where you intend to set your working directory. Once they have been moved to your desired location, unzip your fastq.gz files. Ensure they have the extension .fastq, not .fastq.gz before continuing. 

Next, we will set our working directory. Be sure to do so where your sequence files are located!

```{r}
setwd("state_working_directory_here")
```

To make the location of these files more easily accessible for downstream commands, be sure to assign your working directory (i.e., where our files are located) to an object called path.

```{r}
path<-"state_working_directory_here"
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

The input for the 16S ASV Processsing Pipeline and raw output from your sequencing center will be a .fastq file. These files contain your raw sequencing data and the quality scores associated with each base pair. 

Most often, excluding those generated using Nanopore, your reads will be paired-end and need to be assembled. Prior to assembling our reads, we must trim our sequences to remove low quality base pairs. We will aim to trim away base pairs below a quality score of 30. This is imperative as poor quality data can yield inaccurate results. 

To check the quality of your forward and reverse reads, the following code can be used: 

To visualize the quality profiles of forward reads:

```{r}
plotQualityProfile(fnFs[1:2])
```
To visualize the quality profiles of reverse reads:

```{r}
plotQualityProfile(fnRs[1:2])
```

Reads can now be trimmed and filtered using the filterAndTrim() function. We will use standard filtering parameters: maxN=0, truncQ=2, rm.phix=TRUE and maxEE=c(2,2). The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores. Be sure to change your truncLen parameter to reflect the lengths you want to trim your forward and reverse reads to. 

When truncating reads, you need to make sure you're leaving at least a 50 base pair overlap between forward and reverse reads, as this is necessary for the merging of reads. If you have no sequences left after trimming, you were too stringent.  

```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows be sure to set multithread=FALSE!

head(out)
```

## Step Three: Evaluate Error Rates

The next step involves estimating error rates from our sequencing data using the learnErrors() function. dada2's algorithm uses a parametric error model to distinguish true biological sequences from sequencing errors. Error rates quantify the likelihood of specific nucleotide transitions (e.g., A→C, A→G) that may occur during sequencing. By accurately modeling these errors, dada2 achieves higher taxonomic resolution than traditional approaches based on operational taxonomic units (OTUs).

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

To visualize your error rates, use the script below. 

```{r}
plotErrors(errF, nominalQ=TRUE)
```

We can tell if our rates are good or bad based on whether the estimated error rates (red lines) follow observed error rate (black lines) trends. If there is too much divergence, our data may not be useable. 

Finally, we will apply dada2's core sample inference algorithm to our filtered and trimmed sequence data. With this command we take each sample's reads and use our error rates to infer which sequences are real biological variants and which are sequencing errors. dada() then infers the true error-free ASVs in each sample and how many reads belong to each.
 
```{r}
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
```

## Step Four: Merge Paired-End Reads

For our next step, we will merge our denoised forward and reverse reads to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads. This allows dada2 to construct merged “contig” sequences. 

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

## Step Five: Evaluate and Identify ASVs 

ASVs are unique, error-corrected DNA sequences from amplified marker genes that represent true biological variation. ASVs are grouped only when they share 100% sequence identity, offering greater precision than OTUs, which require only a 97% similarity threshold. Hence, dada2 only makes species level assignments based on exact matching between ASVs and sequenced reference strains.

Let's begin step five by constructing an ASV table, a higher-resolution version of the OTU table produced by traditional methods.

```{r}
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

The following code will allow us to inspect the distribution of sequence lengths.

```{r}
table(nchar(getSequences(seqtab)))
```

While the dada() function corrects substitutions, indel errors and chimeras remain. Fortunately, the accuracy of sequence variants after denoising makes identifying chimeric ASVs simpler than when dealing with OTUs. The following script will construct and remove chimeras from our collection of sequences. 

```{r}
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
```

To get a percentage of chimeric sequences in our original list of ASVs, we can divide the table containing no chimeras with the table containing chimeras.

```{r}
sum(seqtab.nochim)/sum(seqtab)
```

As a final check of our progress, we will look at the number of reads that made it through each step in the pipeline using the commands below:

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

## Step Six: Assign Taxonomy

For step six, we will assign taxonomy to our sequences using the assignTaxonomy() function. This function requires appropriately formatted FASTA files containing taxonomically classified reference sequences to use as a training dataset. To facilitate this process, we must download the SILVA taxonomy database, which provides the necessary reference sequences for accurate taxonomic assignment.

We must state the path to this file on our computer in the assignTaxonomy() function below. Using the script below, we tell dada2 to assign taxonomy to our seqtab.nochim matrix, which contains our non-chimeric ASVs and their abundances.

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/path_to_silva/Silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
```

Below, the addSpecies() function will assign species-level annotation to our taxonomic table.

```{r} 
taxa <- addSpecies(taxa, "/path_to_silva/Silva/silva_v138.2_assignSpecies.fa.gz")
```

Lastly, we will inspect the taxonomic assignment by removing sequence rownames for display purposes.

```{r}
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

## Step Seven: Save your Progress

Finally, we will save our taxa object, containing taxonomic assignments, our seqtab.nochim object, containing our ASVs and their abundances, and our track object, showing the remaining reads following each step of the dada2 pipeline as .csv files. This step will allow us to save our progress and later read in data, in the event our objects are lost or we terminate our R Studio session. 

```{r}
write.csv(taxa, file = "/state_working_directory_here/insert_your_file_name_here_TAXA.csv")

write.csv(seqtab.nochim, file = "/state_working_directory_here/insert_your_file_name_here_SEQTABNOCHIM.csv")

write.csv(track, file = "/state_working_directory_here/insert_your_file_name_here_TRACK.csv")
```

Using the commands below, we can read in the taxa, seqtab.nochim and track objects using the .csv files we made in the previous step.

```{r}
taxa <- read.csv(file = "/state_working_directory_here/insert_your_file_name_here_TAXA.csv")
seqtab.nochim <- read.csv(file = "/state_working_directory_here/insert_your_file_name_here_SEQTABNOCHIM.csv", header = FALSE)
track <- read.csv(file = "/state_working_directory_here/insert_your_file_name_here_TRACK.csv")
```

# Usage Instructions - phyloseq Processing

In the next set of chunks, we will change the formatting of our seqtab.nochim table. This is necessary for downstream work with phyloseq, as well as an optional step to merge seqtab.nochim with our taxa table for easier viewing.

## Step One: Formatting for phyloseq

### Step A: 

To begin this process, we will first transpose our seqtab.nochim table. You must use the seqtab.nochim object from our imported .csv file for this step to run properly.

```{r}
flipped_seqtab.nochim <- as.data.frame(t(seqtab.nochim))
```

### Step B: 

Next, we will copy the first row. 

```{r}
colnames(flipped_seqtab.nochim) <- flipped_seqtab.nochim[1,]
```

### Step C: 

Next, we will delete the first row of our transposed seqtab.nochim table.

```{r}
flipped_seqtab.nochim <- flipped_seqtab.nochim[-1,]
```

### Step D: 

In the command below, we will take each column and its corresponding row name, and paste "ASV"1,2,3.. as a new column next to our nucleotide sequence column.

```{r}
rownames(flipped_seqtab.nochim) <- paste0("ASV", 1:nrow(flipped_seqtab.nochim))
```

### Step E:

Next, we'll remove the nucleotide sequences column and save this new table as flipped_seqtab.nochim_forself for our ease of viewing.  

```{r}
flipped_seqtab.nochim_forself <- flipped_seqtab.nochim[,-1]
```

### Step F: 

We will then save these two tables (the one with our nucleotide sequences and the other without) as .csv files. In the event our R Studio session is lost, we can always read these files back in. 

```{r}
write.csv(flipped_seqtab.nochim, file = "/state_working_directory_here/insert_your_file_name_here_FLIPPEDSEQTABNOCHIM.csv")
write.csv(flipped_seqtab.nochim_forself, file ="/state_working_directory_here/insert_your_file_name_here_FLIPPEDSEQTABNOCHIM_forself.csv")
```

### Step G: 

Next, as an optional step, we can use the cbind() function to save our flipped seqtab.nochim and taxa data as one csv. This approach offers the most comprehensive view of our ASV abundances and taxonomic classifications.

```{r}
OTUabund<-cbind(flipped_seqtab.nochim,taxa)
write.csv(OTUabund,file="/state_working_directory_here/insert_your_file_name_here_OTUabund.csv")
```

### Step H: 

dada2 generates the taxa object with sequences in the first column, a format incompatible with phyloseq. To ensure compatibility, we will remove this column using the command below. 

```{r}
taxa <-taxa[-1]
```

### Step I: 

The script below allows us to verify that our formatting is correct. The first column containing nucleotide sequences should have been deleted.

```{r}
View(taxa)
```

## Step Two: Secondary phyloseq Formatting

Next, we will continue to format our data such that it can be used by phyloseq. 

First, we will create a phyloseq object called "taxmat", containing our taxa data. We will turn it into a matrix so it can be used by phyloseq.

```{r}
taxmat <- as.matrix(taxa)
```

Next, we will make the final formatting changes to our transposed seqtab.nochim table so it can be used by phyloseq. We will call this object otumat, and it will contain our ASV abundance data. Using the existing object "flipped_seqtab.nochim", we will delete the first column listing nucleotide sequences.

```{r}
otumat <- flipped_seqtab.nochim[,-1]
```

Before proceeding, we will ensure ensure both otumat and taxmat are matrices for downstream phyloseq processing using the script below. 

```{r}
class(otumat)
class(taxmat)
```

Next, we will ensure that the row names are ASVs for both matrices. This is also required for downstream phyloseq processing.  

```{r}
rownames(otumat) <- paste0("ASV", 1:nrow(otumat))
rownames(taxmat) <- paste0("ASV", 1:nrow(taxmat))
```

We also must ensure R Studio recognizes our OTU data as numeric and not as character data.

```{r}
class(otumat)<-"numeric"
```

Next, we will tell phyloseq where our "OTU" data can be found (i.e., within the otumat object). 

```{r}
OTU = otu_table(otumat, taxa_are_rows = TRUE)
```

We will also tell it where our "Taxa" data can be found (i.e., within the taxmat object).
 
```{r}
TAX = tax_table(taxmat)
```

Finally, we will tell phyloseq to combine our "OTU" and "Taxa" data into an object called physeq. We will also instruct phyloseq to include the names of our samples in the physeq object. 

```{r}
physeq = phyloseq(OTU, TAX)
physeq
sample_names(physeq)
samplenames<-sample_names(physeq)
```
## Step Three: Graph Production

In the following script, we will create graphs showcasing the absolute and relative abundance of our ASVs by taxa. 

We will begin by creating a bar graph of our sample's absolute abundance by Phylum using the plot_bar() function.

```{r}
phylum_barplot <- plot_bar(physeq, fill = "Phylum")
phylum_barplot
```

Here, the space between each dark line represents the absolute abundance of a particular ASV. We can use the script below to remove these dark lines from our graph.

```{r}
phylum_barstacked <- phylum_barplot + geom_bar(aes(fill=Phylum), stat="identity", position="stack")
phylum_barstacked
```

Now, we will use the tax_glom() phyloseq function to "glom" together ASVs based on our taxonomic assignment of choice. In this example, we will combine our taxonomic data by Phylum to visualize relative abundance by Phylum.

```{r}
g_phylum <- tax_glom(physeq, "Phylum")
```

Next we will use plot_bar() to plot our "glommed" g_phylum graph. This graph will allow us to observe the absolute abundance distribution of our samples more easily.

```{r}
plot_bar(g_phylum, fill="Phylum")
```

Now that we've "glommed" together our taxa by Phylum, we can make a relative abundance graph from our absolute abundance data. We can do this by tallying up the ASVs within each taxa in one sample, and dividing by its total number of ASVs. Then we can use psmelt() to remove phyloseq's formatting and make the data easier to plot.

```{r}
ps_phylum_relabun <- transform_sample_counts(g_phylum, function(ASV) ASV/sum(ASV))
taxa_abundance_table_phylum <- psmelt(ps_phylum_relabun)
```

We will then factor our Phylum values for downstream graphing.

```{r} 
taxa_abundance_table_phylum$Phylum<-factor(taxa_abundance_table_phylum$Phylum)
```

Here we will save our relative abundance table as a .csv file for easy access downstream.  

```{r}
write.csv(taxa_abundance_table_phylum, file = "/state_working_directory_here/insert_your_file_name_here_RelativeAbundance.csv")
```

The following command uses plot_bar() to create a relative abundance table of our samples by Phylum. 

```{r}
p_realabun<-plot_bar(ps_phylum_relabun, fill = "Phylum") + labs(y="Relative Abundance (%)", x= "Sample", title = "3-4cm Core DNA and cDNA & Cyanobacterial Viability Assay Relative Abundance") + theme(plot.title = element_text(size = 12))
p_realabun
```

Should we desire to eliminate the dark lines present between each Phylum, we can assign the position argument within ggplot2's geom_bar() function as "stack".

```{r}
p_abun_stacked<- p_realabun + geom_bar(aes(fill=Phylum), stat="identity", position="stack")
p_abun_stacked
```
* Note: The script above can be adapted to graph the relative abundance of other taxonomic assignments (i.e., Order, Family, etc.)  

## Other Graphing Options 

The 16S Amplicon Sequence Variant Processing Pipeline features code to produce relative abundance graphs using ggplot2. This code will be included below. More detailed instructions can be found in the attached script. 

To produce a geom_point graph, where the size of points reflect taxonomic relative abundance, the following script can be employed. 

```{r}
taxa_abundance_table_order$Abundance[taxa_abundance_table_order$Abundance == 0] <- NA
ggplot(data = taxa_abundance_table_order, aes(x = Sample, y = Order, size = Abundance)) +
  geom_point(aes(colour=Order)) +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "3-4cm Core DNA and cDNA & Cyanobacterial Viability Assay Relative Abundance by Phylum", x = "Sample ", y = "Order", size = "Relative Abundance (%)")
```

To produce a heatmap, in which colour is used to display relative abundance, the following script can be used. 

```{r}
ggplot(taxa_abundance_table_phylum, aes(x = Sample, y = Phylum, fill = Abundance)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("darkgreen", "yellow", "red")) +
  labs(title = "Heatmap of Relative Abundance by Phylum in 3-4cm Core DNA and cDNA & Cyanobacterial Viability Assay",
       x = "Sample",
       y = "Phylum",
       fill = "Relative Abundance (%)") +
  theme(plot.title = element_text(hjust = 0.5))
```

# Technologies Used:

This pipeline utilizes R Studio to process and visualize 16S amplicon sequence variant data. dada2, phyloseq, Biostrings, ggplot2, RColorBrewer and tidyverse are used in this pipeline. The most recent version of SILVA, 138.2 is employed for taxonomic assignment.

# Known Issues or Limitations:

## Concerning filtering and trimming reads in Usage Instructions - dada2 Processing

### Step Two: Quality Control 

Recall that when truncating reads you must leave at least a 50 base pair overlap between forward and reverse reads, as this is necessary for the merging. If you have no sequences left after trimming, you were too stringent.  

If you want to speed up downstream computation, consider tightening maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (i.e., maxEE=c(2,5)), and reducing the truncLen to remove low quality tails. Remember though, when choosing truncLen for paired-end reads you must maintain overlap after truncation in order to merge them later.

### Step Three: 

If your the error rates of your data are poor, your data may not be useable. This pipeline cannot help correct the error rates of your samples. 

### Step Four: Merge Paired-End Reads

Sequences that are much longer or shorter than expected may be the result of non-specific priming. You can remove non-target-length sequences from your sequence table (i.e., seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). This is analogous to “cutting a band” in-silico to get amplicons of the targeted length.  

A low number of merged sequences may be the result of stringent truncation. 

### Step Five: Evaluate and Identify ASVs

Oftentimes, large removals of chimeric reads are a result of failure to remove ambiguous nucleotide primer sequences prior to dada2 pipeline processing. In addition, there is a higher risk for chimeras if reads lack necessary overlap for merging.

## Concerning Issues with Installation  

The phyloseq link stated above can also be used to troubleshoot its installation. Issues with installing packages may be attributed to the version of R Studio you are working with. Ensure installation of packages compatible with your system's version of R Studio. You may need to install newer versions of some packages for compatibility purposes. 

Ensure the most recent version of SILVA is used to assign taxonomy. Be sure it is accessible on your system using the path you input in the assignTaxonomy() function. 

## Concerning Script Compatibility 

Portions of this script will need to be edited such that they can be used by your system. For instance, you must manually denote where your working directory is, rather than using the filler name displayed in this file. In addition, since this script was adapted for a Mac, it may need to be altered such that it is compatible with your system. You may use the link below for troubleshooting purposes. 

https://benjjneb.github.io/dada2/dada-installation.html

# License:

This script is not licensed. It can be freely accessed and distributed via github.

# Contact Information and Acknowledgments:

The dada2 portion of this script was adapted from the following pipeline: 

https://benjjneb.github.io/dada2/tutorial.html

The phyloseq portion of this script can be attributed to Lecture 6 from Jacqueline Goordial's Bioinformatics for Environmental Sciences graduate course (ENVS*6452) taught at the University of Guelph. 

## Bibliography

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon     data. Nature methods, 13(7), 581–583. https://doi.org/10.1038/nmeth.3869

Chuvochina M, Gerken J, Frentrup M, Sandikci Y, Goldmann R, Freese HM, Göker M, Sikorski J, Yarza P, Quast C, Peplies J, Glöckner FO, Reimer LC (2026) SILVA in         2026: a global core biodata resource for rRNA within the DSMZ digital diversity. Nucleic Acids Research, gkaf1247.

McMurdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 8(4):e61217

Neuwirth E (2022). RColorBrewer: ColorBrewer Palettes. R package version 1.1-3, https://CRAN.R-project.org/package=RColorBrewer.

Pagès H, Aboyoun P, Gentleman R, DebRoy S (2025). Biostrings: Efficient manipulation of biological strings. R package version 2.76.0,                                   https://github.com/bioconductor/biostrings

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J,     Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43),        1686. doi:10.21105/joss.01686.
