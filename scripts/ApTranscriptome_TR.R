Thermal reactionome of the common ant species *Aphaenogaster*
================================================================
  
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**Technical Report No. 2**

**Department of Biology**

**University of Vermont**

```{r setup, results='hide'}
# Global settings
options(stringAsFactors=FALSE)

# Load libraries
library(R.utils)
library(ggplot2)
library(knitr)
library(knitcitations)
library(pander)
library(stringr)
#library(plyr)
#library(reshape2)
library(qvalue)

# knitr options
opts_chunk$set(cache=TRUE)
```
    
## Summary ##
  
In this technical report, which accompanies the manuscript "...." (in press), we describe the *de novo* assembly of the transcriptome for two ant colonies with in the *Aphaenogaster rudis-picea-fulva* species complex `r citep("10.1155/2012/752815")`.

This script is completely reproducible assuming that R, `knitr` and the other required libraries (listed below) are installed on a standard linux system using the following:
    
    Rscript -e "library(knitr); knit('ApTran_assemble.Rmd')"

The assembled transcriptome, annotation and expression values are downloaded rather than re-run due to the computational demands, but the exact commands for each of these steps are documented below.

```{r download, echo = FALSE}
# download files
#    * Trimmomatic filtered reads ...
# downloadFile(url = "")
#    * Post-processed Trinity assembly ...
#    * Annotation file ...
# unzip into proper directories
```

## Sample description ##

Two ant colonies were used for the transcriptome sequencing. The first, designated A22, was collected at Molly Bog, Vermont in August 2012 by Nick Gotelli and Andrew Nguyen. This colony was putatively identifed as *A. picea*. The second colony, designated Ar, was collected by Lauren Nichols in Raleigh, North Carolina. These colonies were maintained in the lab for 6 months prior to sample collection. 

For each colony, three ants were exposed to one of 12 temperature treatments, every 3.5C ranging from 0C to 38.5C, for one hour in glass tubes in a water bath. The ants were flash frozen and stored at -80C until RNA was extracted using a two step extraction; [RNAzol RT](http://www.mrcgene.com/rnazol.htm) (Molecular Research Center, Inc) followed by an [RNeasy Micro](http://www.qiagen.com/products/catalog/sample-technologies/rna-sample-technologies/total-rna/rneasy-micro-kit) column (Qiagen). Samples from each colony were pooled and sequenced in separate lanes on a 100bp paired-end run of an Illumina HiSeq at the University of Minnesota Genomics Center, yielding 20e6 and 16e6 reads for the A22 and Ar samples, respectively.

## Transcriptome assembly ##

The Illumina reads were filtered using the program [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) `r citep("10.1093/nar/gks540")` to remove Ilumina adapter sequences and filter out bases with quality scores less than ??. 

~~~
TRIMMOMATIC CODE
~~~

This filtering yielded...

These reads were combined and used in *de novo* transcriptome assembly using the program [Trinity](http://trinityrnaseq.sourceforge.net/) `r citep("10.1038/nbt.1883")`. Note that this required ??? GB RAM and ??? hours run-time and was run on the [Vermont Genetics Network](http://vgn.uvm.edu/) computing cluster. 

~~~
TRINITY CODE
~~~

This assembly contained 100,381 unique components (roughly genes) in 126,172 total transcripts (Table 1). 

As we were assembling two divergent colonies into a single transcriptome, we suspected that this assembly would be susceptible to known problems of errors during assembly (e.g. chimeric transcripts that are fusions of two transcripts) and redundancy `r citep("10.1186/1471-2164-14-328")`. To account for this, we performed two post-assembly processing steps.

First, we ran the program [cap3](http://seq.cs.iastate.edu/) `r citep("10.1101/gr.9.9.868")` setting the maximum gap length and band expansion size to 50 `-f 50 -a 50`, no end clipping as the reads were already filtered `k 0`, requiring 90% identity for assembly, and a minimum overlap length of 100 bp `-o 100`. The percent identity threshold of 90% was chosen to liberally collapse orthologous contigs from the two colonies that may have been assembled separately. 

    cap3 Trinity.fasta -f 50 -a 50 -k 0 -p 90 -o 100 > Trinity_cap3.out

The output of `cap3` gives assembled "contigs" and unassembled "singlets" that were concatenated into a single file.

    # check the number of contigs clustered
    grep -c "Contig" Trinity.fasta.cap.contigs
    grep -c "comp" Trinity.fasta.cap.singlets
    # compare to contigs from Trinity output
    grep -c "comp" Trinity.fasta

    # Combine contigs and singlets from CAP3
    cat Trinity.fasta.cap.contigs Trinity.fasta.cap.singlets > Trinity_cap3.fasta

The output file "Trinity.fasta.cap.info" gives specific information on which contigs were collapsed.

Subsequent to running `cap3`, we ran [uclust](http://drive5.com/usearch/manual/uclust_algo.html) to cluster sequences completely contained within longer sequences, again specificing a 90% identity cutoff for clustering. 

    # sort
    uclust --sort Trinity_cap3.fasta --output Trinity_cap3_sorted.fasta
    # cluster by 90% similarity threshold
    uclust --input Trinity_cap3_sorted.fasta --uc Trinity_cap3_uclust.out --id 0.90
    # convert uclust to fasta format
    uclust --uc2fasta Trinity_cap3_uclust.out --input Trinity_cap3_uclust.fasta


These post-processing step removed `r round((126172-105536)/126172 * 100, 0)`% of the initial reads (Table 1).


```{r assemstats, results = "asis", echo=FALSE}
# make table
trinity <- c("126,172", "100,389,539", "358", "795", "16,201", "1,631")
reduced <- c("105,536", "62,648,997", "320", "593", "15,491", "895")

assemstats <- rbind(trinity, reduced)
colnames(assemstats) <- c("Total contigs", "Total length", "Median contig size", "Mean contig size", "N50 contig", "N50 Length")

pandoc.table(assemstats, style="rmarkdown", caption = "Table 1: Statistics for Trinity and cap3+uclust reduced transcriptome assemblies")
```


## Annotation ##

Annotation was performed by uploading the reduced assembly "Trinity_cap3_uclust.fasta" to the web-based annotation program [FastAnnotator](http://fastannotator.cgu.edu.tw/index.php) `r citep("10.1186/1471-2164-13-S7-S9")`.

Results are available as job ID [13894410176993](http://fastannotator.cgu.edu.tw/job.php?jobid=13894410176993#page=basicinfo).


## Identification of thermally-responsive genes ##

Expression was quantified using the program [sailfish]()                                              
    export PYTHONPATH=/opt/software/khmer/python
    export LD_LIBRARY_PATH=/opt/software/Sailfish-0.6.2-Linux_x86-64/lib:$LD_LIBRARY_PATH
    export PATH=/opt/software/Sailfish-0.6.2-Linux_x86-64/bin:$PATH



## Quantify gene expression                                          

Quantify gene expression using [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html). In shell, set up paths to software libraries:
                                                 
    export LD_LIBRARY_PATH=/opt/software/Sailfish-0.6.2-Linux_x86-64/lib:$LD_LIBRARY_PATH
    export PATH=/opt/software/Sailfish-0.6.2-Linux_x86-64/bin:$PATH

Then build the index of the assembly:

    sailfish index -t Trinity_cap3_uclust.fasta -o sailfish-index -k 20 -p 4

Once this is done, quantify expression for the Trimmomatic filtered reads from each colony-treatment sample separately. Note that for each sample, there are three four filtered read files:

- paired.left.fastq
- paired.right.fastq
- unpaired.left.fastq
- unpaired.right.fastq
                                                 
Make a directory for the expression values

    mkdir sailfish-expression                                                 

Then, for each sample, run the following command:
                                                 
    sailfish -i sailfish-index -o sailfish-expression/A22-0 --reads A22-0_ATCACG.paired.left.fastq A22-0_ATCACG.paired.right.fastq A22-0_ATCACG.unpaired.left.fastq A22-0_ATCACG.unpaired.right.fastq
-p 4

which I looped in an R script:                                                 
                                                 
```{r A22expression, eval=FALSE}
# directory containing trimmed reads
readdir <- "../data/ind_files/" 

# list of reads in each of four trimmed classes
readlist <- list.files(readdir)
(paired.left <- readlist[grep("\\.paired.left.fastq$", readlist)])
(paired.right <- readlist[grep("\\.paired.right.fastq$", readlist)])
(unpaired.left <- readlist[grep("unpaired.left.fastq$", readlist)])
(unpaired.right <- readlist[grep("unpaired.right.fastq$", readlist)])

# Loop across each sample and quantify expression

samples <- c("A22-0", "A22-3", "A22-7", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", "A22-31", "A22-35", "A22-38", "Ar-0", "Ar-3", "Ar-7", "Ar-10", "Ar-14", "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38")

for (j in 1:length(samples)) {
    message("Start expression quantification for sample ", samples[j], ": ", Sys.time())
    quantdir <- paste(samples[j], "_quant", sep="")
    samp.paired.l <- paste(readdir, paired.left[j], sep="")
v    samp.paired.r <- paste(readdir, paired.right[j], sep="")
    samp.unpaired.l <- paste(readdir, unpaired.left[j], sep="")
    samp.unpaired.r <- paste(readdir, unpaired.right[j], sep="")
        
    system(paste("sailfish quant -i ../results/trinity-full/sailfish-index -o ../results/trinity-full/sailfish-expression/", quantdir, " --reads ", samp.paired.l, " ", samp.paired.r, " ", samp.unpaired.l, " ", samp.unpaired.r, " -p 4", sep=""))

    message("Done with expression quantification for sample ", samples[j], ": ", Sys.time())
}
```

This generated a directory for each sample `r list.files("../results/trinity-full/sailfish-expression")` and within each directory there are the following files `r list.files("../results/trinity-full/sailfish-expression/A22-0_quant")`

The file *quant_bias_corrected.sf* contains the following columns, following a number of header lines:

1. Transcript ID
2. Transcript Length
3. Transcripts per Million (TPM): computed as described in `r citep("10.1093/bioinformatics/btp692")`, and is meant as an estimate of the number of transcripts, per million observed transcripts, originating from each isoform.
4. Reads Per Kilobase per Million mapped reads (RPKM): classic measure of relative transcript abundance, and is an estimate of the number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

The TPM column for each sample was extracted and combined into a matrix for each colony.

```{r expression_matrix, eval=FALSE, echo=FALSE}
# function to read sailfish-expression file

read.sample <- function(filein, outname) {
    file.df <- read.table(filein, header=FALSE, sep="\t")
    colnames(file.df) <- c("Transcript", "Length", "TPM", "RPKM", "KPKM", "EstimatedNumReads")
    #head(file.df)
    assign(outname, file.df, envir = .GlobalEnv)
}
    
#  read in each file using loop
samples <- c("A22-0", "A22-3", "A22-7", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", "A22-31", "A22-35", "A22-38", "Ar-0", "Ar-3", "Ar-7", "Ar-10", "Ar-14", "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38")

fileinpath <- "../results/trinity-full/sailfish-expression/"
    
for (j in 1:length(samples)) {
    samp <- samples[j]
    outpre <- gsub("-", "_", samp)
    outname <- paste(outpre, "_quant", sep="")
    read.sample(filein=paste(fileinpath, samp, "_quant/quant_bias_corrected.sf", sep=""), outname=outname)
}

# combine TPM into single matrix
A22.TPM <- cbind(A22_0_TPM=A22_0_quant$TPM, A22_3_TPM=A22_3_quant$TPM, A22_7_TPM=A22_7_quant$TPM, A22_10_TPM=A22_10_quant$TPM, A22_14_TPM=A22_14_quant$TPM, A22_17_TPM=A22_17_quant$TPM, A22_21_TPM=A22_21_quant$TPM, A22_24_TPM=A22_24_quant$TPM, A22_28_TPM=A22_28_quant$TPM, A22_31_TPM=A22_31_quant$TPM, A22_35_TPM=A22_35_quant$TPM, A22_38_TPM=A22_38_quant$TPM)
rownames(A22.TPM) <- A22_0_quant[,"Transcript"]
head(A22.TPM)


Ar.TPM <- cbind(Ar_0_TPM=Ar_0_quant$TPM, Ar_3_TPM=Ar_3_quant$TPM, Ar_7_TPM=Ar_7_quant$TPM, Ar_10_TPM=Ar_10_quant$TPM, Ar_14_TPM=Ar_14_quant$TPM, Ar_17_TPM=Ar_17_quant$TPM, Ar_21_TPM=Ar_21_quant$TPM, Ar_24_TPM=Ar_24_quant$TPM, Ar_28_TPM=Ar_28_quant$TPM, Ar_31_TPM=Ar_31_quant$TPM, Ar_35_TPM=Ar_35_quant$TPM, Ar_38_TPM=Ar_38_quant$TPM)
rownames(Ar.TPM) <- Ar_0_quant[,"Transcript"]
head(Ar.TPM)
```


## Identification of thermally-responsive genes

For each colony, identify genes that have a significant linear or quadratic regression.

```{r}
# Define function to report overall model P-value

lmp <- function(modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
            f <- summary(modelobject)$fstatistic
            p <- unname(pf(f[1],f[2],f[3],lower.tail=F))
            attributes(p) <- NULL
            return(p)
    }

# previous results suggest that A22_7 and Ar_7 samples were switched. conservatively remove from further analyses

temps <- c(0, 3.5, 10.5, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)

A22.TPM.sub <- subset(A22.TPM, , select=-A22_7_TPM)
head(A22.TPM.sub)

# Fit quadratic regression to each transcript, retaining only transcripts with significant model

signif <- vector(length=0)
pvals <- vector(length=0)

#for(i in 1:nrow(A22.TPM.sub)) {
for(i in 1:1000) {       
    if(i%%1000 == 0) message(i, " regressions done!")
    out <- lm(unlist(A22.TPM.sub[i,]) ~ temps + I(temps^2))
    #print(summary(out))
    pvals <- c(pvals, lmp(out))
    if(is.na(lmp(out))) print("NA") else {
        if(lmp(out) < 0.05) signif <- c(signif, i)
    }
}

length(signif)


# False discovery rate correction with q-value
length(pvals)
# Remove NA pvals
pvals.complete <- pvals[which(!is.na(pvals))]
length(pvals.complete)
# Calculate q-values
qvals <- qvalue(pvals.complete)
length(qvals$qvalues)
qsummary(qvals)
# q-value diagnostics
png("qval-diag.png")
qplot(qvals)
dev.off()

png("qval-hist.png")
par(mfrow = c(2,1))
  hist(pvals)
  hist(qvals$qvalues)
dev.off()


# Make table with responsive genes reporting both p-val and q-val
responsive.transcripts <- as.character(data[,1])
# Add pvals, which include 'NA's
responsive.transcripts <- cbind(responsive.transcripts, pvals = as.numeric(pvals))
# Drop NA rows, then add qvalues
responsive.transcripts <- as.data.frame(responsive.transcripts, stringsAsFactors = FALSE)
responsive.transcripts <- responsive.transcripts[which(responsive.transcripts$pvals!="NaN"), ]
#dim(responsive.transcripts)
# add qvalues
responsive.transcripts <- cbind(responsive.transcripts, qvals=qvals$qvalues)
responsive.transcripts$pvals <- as.numeric(responsive.transcripts$pvals)
#head(responsive.transcripts)
#str(responsive.transcripts)

# Order by p-value
responsive.transcripts.by.p <- responsive.transcripts[order(responsive.transcripts$pvals), ]
head(responsive.transcripts.by.p)
tail(responsive.transcripts.by.p)

# Take only significant hits, FDR < 0.05
responsive.transcripts.FDR <- responsive.transcripts.by.p[which(responsive.transcripts.by.p$qvals < 0.05),]

# Write all results to file
write.table(file="responsive_transcripts_by_pval.txt", responsive.transcripts.by.p, row.names=FALSE, quote=FALSE, sep="\t")
# Write only significant results to file
write.table(file="responsive_transcripts_by_pval.txt", responsive.transcripts.FDR, row.names=FALSE, quote=FALSE, sep="\t")

### Plot significant transcripts by q-value
# scale values
combined.data.t <- cbind(Transcript_ID = data$Transcript_ID, combined.data)
head(combined.data.t)
signif.transcripts.FDR <- merge(combined.data.t, responsive.transcripts.FDR, by.x = "Transcript_ID", by.y = "responsive.transcripts", all.x = FALSE, all.y = TRUE)
dim(signif.transcripts.FDR)
head(signif.transcripts.FDR)

m <- signif.transcripts.FDR[,2:12]
head(m)

signif.transcripts.FDR.scaled <- data.frame(m[0,])

for (j in 1:nrow(m)) {
    signif.transcripts.FDR.scaled[j,] <- scale(unlist(m[j,]))
}


# 

signif.transcripts <- combined.data[signif, ]
dim(signif.transcripts)

signif.transcripts.scaled <- data.frame(signif.transcripts[0,])
signif.transcripts.scaled

for (j in 1:nrow(signif.transcripts)) {
    signif.transcripts.scaled[j,] <- scale(unlist(signif.transcripts[j,]))
}

head(signif.transcripts.FDR.scaled)


## Plot

# determine min and max values


(ex.min <- min(unlist(signif.transcripts.scaled)))
(ex.max <- max(unlist(signif.transcripts.scaled)))

# Plot loess smooth for each significant transcript


pdf("expression-plot.pdf")
  plot(temps, signif.transcripts.scaled[1,], ylim=c(ex.min, ex.max), ylab="Expected count", xlab="Temp C", pch=16, col="white")
  for(i in 1:nrow(signif.transcripts.scaled)) {
      lines(temps, signif.transcripts.scaled[i,], type="l")
  }
dev.off()


pdf("loess-expression-plot.pdf")
  plot(1, xlim=c(0, 40), ylim=c(ex.min, ex.max), ylab="Expected count", xlab="Temp C", pch=16, col="white")
  for(i in 1:nrow(signif.transcripts.scaled)) {
      loo <- loess(unlist(signif.transcripts.scaled[i,]) ~ temps)      
      lines(predict(loo))
  }
dev.off()

# Make standardized version of all data
data.scaled <- data.frame(data[0,])

for (j in 1:nrow(data)) {
    data.scaled[j,] <- scale(unlist(data[j,]))
}
```

Get annotation of responsive genes. 

```{r }
#### Addendum ####
## compare sailfish to RSEM results

# Read data file for 'expected_counts' from RSEM
#data <- read.table("../data/combined.exp_count.matrix", header=T, sep="\t")
#dim(data)
#data[1:10,1:4]
#str(data)

```


```{r session}
sessionInfo()
```
