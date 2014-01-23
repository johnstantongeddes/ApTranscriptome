Thermal reactionome of the common ant species *Aphaenogaster*
================================================================
  
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**Technical Report No. 3**

**Department of Biology**

**University of Vermont**

```{r setup, results='hide', message = FALSE}
# Global settings
options(stringsAsFactors=FALSE)

# Load libraries
library(R.utils)
library(ggplot2)
library(knitr)
library(knitcitations)
library(pander)
library(stringr)
library(data.table)
library(RCurl)

# Load Bioconductor libraries
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
#biocLite("qvalue")
#biocLite("Rgraphviz")
library(qvalue) 
library(topGO) 
library(Rgraphviz) 


# knitr options
opts_chunk$set(cache=TRUE)

# load custom functions
source("RxNseq.R")
```
    
## Summary ##
  
In this technical report, which accompanies the manuscript "...." (in press), we describe the *de novo* assembly of the transcriptome for two ant colonies with in the *Aphaenogaster rudis-picea-fulva* species complex `r citep("10.1155/2012/752815")`.

This script is completely reproducible assuming that R, `knitr` and the other required libraries (listed below) are installed on a standard linux system using the following:
    
    Rscript -e "library(knitr); knit('ApTran_assemble.Rmd')"

The assembled transcriptome, annotation and expression values are downloaded rather than re-run due to the computational demands, but the exact commands for each of these steps are documented below.

```{r download, echo = TRUE, results = "hide"}
### Transcriptome assembly

### Annotation file
# from either AWS or GoogleDrive
annotationURL <- getURL("http://johnstantongeddes.org/assets/files/Aphaeno_transcriptome_AnnotationTable.txt")
#a2 <- getURL("https://googledrive.com/host/0B75IymziRJ_9Tlg1U1Vxbjk1bzg") # GoogleDrive link

annotationfile <- read.csv(textConnection(annotationURL), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(annotationfile)
str(annotationfile)

# Convert to data.table
annotationtable <- data.table(annotationfile)
head(annotationtable)
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

which I looped in an R script.                                                 
                                                 
```{r sailfish, eval=FALSE, echo=FALSE}
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
    samp.paired.r <- paste(readdir, paired.right[j], sep="")
    samp.unpaired.l <- paste(readdir, unpaired.left[j], sep="")
    samp.unpaired.r <- paste(readdir, unpaired.right[j], sep="")
        
    system(paste("sailfish quant -i ../results/trinity-full/sailfish-index -o ../results/trinity-full/sailfish-expression/", quantdir, " --reads ", samp.paired.l, " ", samp.paired.r, " ", samp.unpaired.l, " ", samp.unpaired.r, " -p 4", sep=""))

    message("Done with expression quantification for sample ", samples[j], ": ", Sys.time())
}
```

This generated a directory for each sample

`r list.files("../results/trinity-full/sailfish-expression")`

and within each directory there are the following files:

`r list.files("../results/trinity-full/sailfish-expression/A22-0_quant")`

The file *quant_bias_corrected.sf* contains the following columns, following a number of header lines:

1. Transcript ID
2. Transcript Length
3. Transcripts per Million (TPM): computed as described in `r citep("10.1093/bioinformatics/btp692")`, and is meant as an estimate of the number of transcripts, per million observed transcripts, originating from each isoform.
4. Reads Per Kilobase per Million mapped reads (RPKM): classic measure of relative transcript abundance, and is an estimate of the number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

The TPM column for each sample was extracted and combined into a matrix for each colony.

```{r expression_matrix, eval=TRUE, echo=TRUE}
    
#  read in each file using loop
samples <- c("A22-0", "A22-3", "A22-7", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", "A22-31", "A22-35", "A22-38", "Ar-0", "Ar-3", "Ar-7", "Ar-10", "Ar-14", "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38")

fileinpath <- "../results/trinity-full/sailfish-expression/"
    
for (j in 1:length(samples)) {
    samp <- samples[j]
    outpre <- gsub("-", "_", samp)
    outname <- paste(outpre, "_quant", sep="")
    read.sailfish.quant(filein=paste(fileinpath, samp, "_quant/quant_bias_corrected.sf", sep=""), outname=outname)
}

# combine TPM into single data.frame
A22.TPM <- data.frame(Transcript = A22_0_quant$Transcript, A22_0_TPM=A22_0_quant$TPM, A22_3_TPM=A22_3_quant$TPM, A22_7_TPM=A22_7_quant$TPM, A22_10_TPM=A22_10_quant$TPM, A22_14_TPM=A22_14_quant$TPM, A22_17_TPM=A22_17_quant$TPM, A22_21_TPM=A22_21_quant$TPM, A22_24_TPM=A22_24_quant$TPM, A22_28_TPM=A22_28_quant$TPM, A22_31_TPM=A22_31_quant$TPM, A22_35_TPM=A22_35_quant$TPM, A22_38_TPM=A22_38_quant$TPM)
head(A22.TPM)
str(A22.TPM)

Ar.TPM <- data.frame(Transcript = Ar_0_quant$Transcript, Ar_0_TPM=Ar_0_quant$TPM, Ar_3_TPM=Ar_3_quant$TPM, Ar_7_TPM=Ar_7_quant$TPM, Ar_10_TPM=Ar_10_quant$TPM, Ar_14_TPM=Ar_14_quant$TPM, Ar_17_TPM=Ar_17_quant$TPM, Ar_21_TPM=Ar_21_quant$TPM, Ar_24_TPM=Ar_24_quant$TPM, Ar_28_TPM=Ar_28_quant$TPM, Ar_31_TPM=Ar_31_quant$TPM, Ar_35_TPM=Ar_35_quant$TPM, Ar_38_TPM=Ar_38_quant$TPM)
head(Ar.TPM)
str(Ar.TPM)
```

## Identification of thermally-responsive genes

For each colony, identify genes that have a significant linear or quadratic regression by fitting the linear model to the expression levels at each temperature for each transcript

$$ TPM ~ temp + temp^2 $$

For this list of P-values, False Discovery Rate (FDR) is applied and q-values are calculated using the [qvalue]() package.

```{r RxN, echo=TRUE, eval=TRUE}
# Previous results suggest that A22_7 and Ar_7 samples were switched. conservatively remove from further analyses

temps <- c(0, 3.5, 10.5, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)

# Identify thermally-responsive genes for A22
A22.TPM.sub <- subset(A22.TPM, , select=-A22_7_TPM)
head(A22.TPM.sub)

RxNseq(mat = A22.TPM.sub, vals = temps, qcrit = 0.05, makeplots = TRUE, prefix = "A22")

# Identify thermally-responsive genes for Ar
Ar.TPM.sub <- subset(Ar.TPM, , select=-Ar_7_TPM)
head(Ar.TPM.sub)

RxNseq(mat = Ar.TPM.sub, vals = temps, qcrit = 0.05, makeplots = TRUE, prefix = "Ar")

save.image("RxN_results.RData")
```

Examine

1. relationships among responsive transcripts between A22 and Ar
2. patterns of responsiveness with temperature

```{r RxN_qcrit}

## Convert RxN df to data.tables for fast sort and merge
A22_RxN_table <- data.table(A22_RxN)
head(A22_RxN_table)
Ar_RxN_table <- data.table(Ar_RxN)
head(Ar_RxN_table)

## Set keys
setkey(annotationtable, Sequence.Name)
setkey(A22_RxN_table, Transcript)
setkey(Ar_RxN_table, Transcript)
# check tables and keys
tables()

# Join tables
A22_DT <- annotationtable[A22_RxN_table]
str(A22_DT)
Ar_DT <- annotationtable[Ar_RxN_table]
str(Ar_DT)

# Pull out transcripts with Q < 0.05 for each colony, order by qval
A22_DT_qcrit <- A22_DT[A22_DT$qval < 0.05]
A22_DT_qcrit <- A22_DT_qcrit[order(qval)]
dim(A22_DT_qcrit)
Ar_DT_qcrit <- Ar_DT[Ar_DT$qval < 0.05]
Ar_DT_qcrit <- Ar_DT_qcrit[order(qval)]
dim(Ar_DT_qcrit)

setkey(A22_DT_qcrit, Sequence.Name)
setkey(Ar_DT_qcrit, Sequence.Name)

# Save significant transcripts with best hits and GO to file
write.table(A22_DT_qcrit[,list(Sequence.Name, best.hit.to.nr, GO.Biological.Process, GO.Cellular.Component, GO.Molecular.Function, Enzyme, Domain, annotation.type, pval, qval, lin.coef, quad.coef)], file = "A22_signif_transcripts_GO.txt", quote = FALSE, sep = "\t", row.names = FALSE)

write.table(Ar_DT_qcrit[,list(Sequence.Name, best.hit.to.nr, GO.Biological.Process, GO.Cellular.Component, GO.Molecular.Function, Enzyme, Domain, annotation.type, pval, qval, lin.coef, quad.coef)], file = "Ar_signif_transcripts_GO.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```

**1. relationships among responsive transcripts between A22 and Ar**

Is there overlap in signficant transcripts between the colonies?

`r round(length(which(A22_DT_qcrit$Sequence.Name %in% Ar_DT_qcrit$Sequence.Name))/length(A22_DT_qcrit$Sequence.Name), 2) * 100`% of responsive transcripts are shared among the two colonies!

Are the responses in the same direction?

Correlation among intercepts is: `r cor.test(A22_DT_qcrit[which(A22_DT_qcrit$Sequence.Name %in% Ar_DT_qcrit$Sequence.Name), intercept], Ar_DT_qcrit[which(Ar_DT_qcrit$Sequence.Name %in% A22_DT_qcrit$Sequence.Name), intercept])`

Correlation among linear coefficients is: `r cor.test(A22_DT_qcrit[which(A22_DT_qcrit$Sequence.Name %in% Ar_DT_qcrit$Sequence.Name), lin.coef], Ar_DT_qcrit[which(Ar_DT_qcrit$Sequence.Name %in% A22_DT_qcrit$Sequence.Name), lin.coef])`

Correlation among quadratic coefficients is: `r cor.test(A22_DT_qcrit[which(A22_DT_qcrit$Sequence.Name %in% Ar_DT_qcrit$Sequence.Name), quad.coef], Ar_DT_qcrit[which(Ar_DT_qcrit$Sequence.Name %in% A22_DT_qcrit$Sequence.Name), quad.coef])`

To summarize:

- 98% of significant transcripts the same in Ar and A22
- 99% correlation in linear models for the significant transcripts for Ar and A22

**2. patterns of responsiveness with temperature**

Distribution of intercepts, linear and quadratic coefficients for each colony.

```{r lm_coefficients, echo = FALSE}
# quantiles
quantile(A22_DT_qcrit$intercept)
quantile(Ar_DT_qcrit$intercept)

quantile(A22_DT_qcrit$lin.coef, na.rm = TRUE)
quantile(Ar_DT_qcrit$lin.coef, na.rm = TRUE)

quantile(A22_DT_qcrit$quad.coef, na.rm = TRUE)
quantile(Ar_DT_qcrit$quad.coef, na.rm = TRUE)

# percent of linear coefficients that are positive


```

Number of A22 transcripts with significant linear responses:
`r length(which(!is.na(A22_DT_qcrit$lin.coef)))`

Number of A22 transcripts with significant quadratic responses:
`r length(which(!is.na(A22_DT_qcrit$quad.coef)))`

Percent of A22 transcripts that significantly increase expression with temperature:
`r round(length(which(A22_DT_qcrit$lin.coef > 0))/length(which(!is.na(A22_DT_qcrit$lin.coef)))*100, 2)`

Percent of A22 transcripts that have significant negative curvature:
`r
round(length(which(A22_DT_qcrit$quad.coef > 0))/length(which(!is.na(A22_DT_qcrit$quad.coef)))*100, 2)`


Number of Ar transcripts with significant linear responses:
`rlength(which(!is.na(Ar_DT_qcrit$lin.coef)))`

Number of Ar transcripts with significant quadratic responses:
`r length(which(!is.na(Ar_DT_qcrit$quad.coef)))`

Percent of Ar transcripts that significantly increase expression with temperature:
`r round(length(which(Ar_DT_qcrit$lin.coef > 0))/length(which(!is.na(Ar_DT_qcrit$lin.coef)))*100, 2)`

Percent of Ar transcripts that have significant negative curvature:
`r round(length(which(Ar_DT_qcrit$quad.coef > 0))/length(which(!is.na(Ar_DT_qcrit$quad.coef)))*100, 2)`


```{r hists, echo=FALSE, eval=FALSE}

## IN PROGRESS ##

# combine data
# reshape long
# make histogram

png("hist_intercept.png")
  h <- ggplot(A22_DT_qcrit, aes(intercept)) + geom_histogram(binwidth = 0.5) +
#   scale_x_log10()
    xlim(0, 100)
  h
dev.off()

png("hist_linear.png")
  h <- ggplot(A22_DT_qcrit, aes(lin.coef)) + geom_histogram(binwidth = 0.5) +
    xlim(-5,5)
  h
dev.off()

png("hist_quad.png")
  h <- ggplot(A22_DT_qcrit, aes(quad.coef)) + geom_histogram(binwidth = 0.1) +
    xlim(-1,1)
  h
dev.off()
```


## Functional annotation and gene set enrichment analysis ##

Get annotation of responsive genes. 

Use [topGO](http://www.bioconductor.org/packages/2.12/bioc/html/topGO.html) to perform gene set enrichment analysis

First need to create gene ID to GO term map file

```{r geneid2go, echo = TRUE, eval = FALSE}

# create geneid2go.map file from FastAnnotator AnnotationTable.txt
geneid2GOmap(annotationfile)
```

Using this gene2GO map file, perform gene set enrichment analysis.

Provide qvalues as the gene score, and select genes with q < 0.05 using custom `selectFDR` function.

```{r gene_enrichment_A22}

# read mappings file
geneID2GO <- readMappings(file = "geneid2go.map")
str(head(geneID2GO))

# create geneList. note that NA values cause problems with topGO
# need to retain these for full ontology, so set NA to 1
A22geneList <- A22_RxN$qval
A22geneList[which(is.na(A22geneList))] <- 1
names(A22geneList) <- A22_RxN$Transcript
str(A22geneList)

# create function to select significant genes for topGO
selectFDR <- function(qvalue) {
    return(qvalue < 0.05)
}

# create topGOdata object
BP_GOdata <- new("topGOdata",
                 description = "BP gene set analysis", ontology = "BP",
                 allGenes = A22geneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

BP_GOdata

# perform enrichment analysis using multiple methods
BP.resultParentChild <- runTest(BP_GOdata, statistic = 'fisher', algorithm = 'parentchild')
BP.resultParentChild

BP.ResTable <- GenTable(BP_GOdata, parentchild = BP.resultParentChild, topNodes = 40)
BP.ResTable

# graph significant nodes

pdf("BP_topGO_sig_nodes.pdf")
showSigOfNodes(BP_GOdata, score(BP.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()

## Cellular Component

# create topGOdata object
CC_GOdata <- new("topGOdata",
                 description = "CC gene set analysis", ontology = "CC",
                 allGenes = A22geneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

CC_GOdata

# perform enrichment analysis using multiple methods
CC.resultParentChild <- runTest(CC_GOdata, statistic = 'fisher', algorithm = 'parentchild')
CC.resultParentChild

CC.ResTable <- GenTable(CC_GOdata, parentchild = CC.resultParentChild, topNodes = 40)
CC.ResTable

# graph significant nodes

pdf("CC_topGO_sig_nodes.pdf")
showSigOfNodes(CC_GOdata, score(CC.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()


#### Molecular Function

# create topGOdata object
MF_GOdata <- new("topGOdata",
                 description = "MF gene set analysis", ontology = "MF",
                 allGenes = A22geneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

MF_GOdata

# perform enrichment analysis using multiple methods
MF.resultParentChild <- runTest(MF_GOdata, statistic = 'fisher', algorithm = 'parentchild')
MF.resultParentChild

MF.ResTable <- GenTable(MF_GOdata, parentchild = MF.resultParentChild, topNodes = 40)
MF.ResTable

# graph significant nodes

pdf("MF_topGO_sig_nodes.pdf")
showSigOfNodes(MF_GOdata, score(MF.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()
```

Perform gene enrichment analysis for Ar

```{r gene_enrichment_Ar}

# create geneList. note that NA values cause problems with topGO
# need to retain these for full ontology, so set NA to 1
A22geneList <- Ar_RxN$qval
ArgeneList[which(is.na(ArgeneList))] <- 1
names(ArgeneList) <- Ar_RxN$Transcript
str(ArgeneList)

# create topGOdata object
Ar.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis", ontology = "BP",
                 allGenes = ArgeneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

Ar.BP.GOdata

# perform enrichment analysis using multiple methods
Ar.BP.resultParentChild <- runTest(Ar.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ar.BP.resultParentChild

Ar.BP.ResTable <- GenTable(Ar.BP.GOdata, parentchild = Ar.BP.resultParentChild, topNodes = 40)
Ar.BP.ResTable

# graph significant nodes

pdf("Ar.BP_topGO_sig_nodes.pdf")
showSigOfNodes(Ar.BP_GOdata, score(Ar.BP.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()

## Cellular Component

# create topGOdata object
CC_GOdata <- new("topGOdata",
                 description = "CC gene set analysis", ontology = "CC",
                 allGenes = ArgeneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

CC_GOdata

# perform enrichment analysis using multiple methods
Ar.CC.resultParentChild <- runTest(Ar.CC_GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ar.CC.resultParentChild

Ar.CC.ResTable <- GenTable(CC_GOdata, parentchild = Ar.CC.resultParentChild, topNodes = 40)
Ar.CC.ResTable

# graph significant nodes

pdf("Ar.CC_topGO_sig_nodes.pdf")
showSigOfNodes(Ar.CC_GOdata, score(Ar.CC.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()


#### Molecular Function

# create topGOdata object
Ar.MF_GOdata <- new("topGOdata",
                 description = "MF gene set analysis", ontology = "MF",
                 allGenes = ArgeneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

Ar.MF_GOdata

# perform enrichment analysis using multiple methods
Ar.MF.resultParentChild <- runTest(Ar.MF_GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ar.MF.resultParentChild

Ar.MF.ResTable <- GenTable(MF_GOdata, parentchild = Ar.MF.resultParentChild, topNodes = 40)
Ar.MF.ResTable

# graph significant nodes

pdf("Ar.MF_topGO_sig_nodes.pdf")
showSigOfNodes(Ar.MF_GOdata, score(Ar.MF.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()
```


Finish ...

Make plots of thermally-responsive genes.

```{r plots, echo=FALSE, eval=FALSE}
### Plot scaled expression values of significant transcripts if makeplots==TRUE

# make new dataframe with scaled values
        m <- subset(responsive.transcripts.FDR, , select = -c(pvals, qvals))
        m.scaled <- data.frame(m[0,-1])
    
        for (j in 1:nrow(m)) {
             m.scaled[j,] <- scale(unlist(m[j,2:ncol(m)])) # skip first column which is transcript ID
        }

        # determine min and max values
        (ex.min <- min(unlist(m.scaled)))
        (ex.max <- max(unlist(m.scaled)))

    # Plot loess smooth for each significant transcript


    pdf(paste(prefix, "_expression_plot.pdf", sep = ""))
      plot(vals, m.scaled[1,], ylim=c(ex.min, ex.max), ylab="Expected count", xlab="Temp C", pch=16, col="white")
      for(i in 1:nrow(m.scaled)) {
          lines(vals, m.scaled[i,], type="l")
      }
    dev.off()


    pdf(paste(prefix, "_loess_expression_plot.pdf", sep = ""))
      plot(1, xlim=c(0, 40), ylim=c(ex.min, ex.max), ylab="Expected count", xlab="Temp C", pch=16, col="white")
      for(i in 1:nrow(m.scaled)) {
        loo <- loess(unlist(m.scaled[i,]) ~ vals)      
        lines(predict(loo))
      }
    dev.off()

```

## Session information ##

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
save.image()
sessionInfo()
```
