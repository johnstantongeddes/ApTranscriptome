Thermal reactionome of the common ant species *Aphaenogaster*
================================================================
  
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**Technical Report No. 3**

**Department of Biology**

**University of Vermont**

```{r setup, echo=FALSE, results='hide', message = FALSE}
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
library(plyr)

# Load Bioconductor libraries
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
#biocLite("qvalue")
#biocLite("Rgraphviz")
suppressMessages(library(qvalue))
library(topGO) 
suppressMessages(library(Rgraphviz))

# load custom functions
source("RxNseq.R")
```

## Summary ##
  
In this technical report, which accompanies the manuscript **Thermal reactionome of a common ant species** (Stanton-Geddes et al., in press), we:

1. describe the *de novo* assembly of the transcriptome for two ant colonies with in the *Aphaenogaster rudis-picea-fulva* species complex `r citep("10.1155/2012/752815")`
2. identify thermally-responsive genes
3. perform gene set enrichment analysis

This script is completely reproducible assuming that R, `knitr` and the other required libraries (listed within the source document) are installed on a standard linux system using the following:
    
    Rscript -e "library(knitr); knit('ApTranscriptome_TR.Rmd')"

The assembled transcriptome, annotation and expression values are downloaded rather than re-run due to the computational demands, but the exact commands for each of these steps are documented below.

```{r download, echo = TRUE, results = "hide"}
### Transcriptome assembly

### Expression data

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


## Transcriptome annotation ##

Annotation was performed by uploading the reduced assembly "Trinity_cap3_uclust.fasta" to the web-based annotation program [FastAnnotator](http://fastannotator.cgu.edu.tw/index.php) `r citep("10.1186/1471-2164-13-S7-S9")`.

Results are available as job ID [13894410176993](http://fastannotator.cgu.edu.tw/job.php?jobid=13894410176993#page=basicinfo).


## Identification of thermally-responsive genes ##

### Quantify gene expression ###

Quantify gene expression using [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html). Make sure that PATHs to the software libraries are set up correctly: 
                                                 
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
                                                 
    sailfish -i sailfish-index -o sailfish-expression/A22-0 --reads A22-0_ATCACG.paired.left.fastq A22-0_ATCACG.paired.right.fastq A22-0_ATCACG.unpaired.left.fastq A22-0_ATCACG.unpaired.right.fastq -p 4

Or, with a loop:                                                 
                                                 
```{r sailfish, eval=TRUE, echo=TRUE, cache=TRUE}
# directory containing trimmed reads
readdir <- "../data/ind_files/" 

# list of reads in each of four trimmed classes
readlist <- list.files(readdir)
(paired.left <- readlist[grep(".\\.paired.left.fastq$", readlist)])
(paired.right <- readlist[grep("\\.paired.right.fastq$", readlist)])
(unpaired.left <- readlist[grep("unpaired.left.fastq$", readlist)])
(unpaired.right <- readlist[grep("unpaired.right.fastq$", readlist)])

# Loop across each sample and quantify expression

# NOTE - samples listed in same order as given by the above lists
samples <- c("A22-0", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", "A22-31", "A22-35", "A22-38", "A22-3", "A22-7", "Ar-0", "Ar-10", "Ar-14", "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38", "Ar-3", "Ar-7")

for (j in 1:length(samples)) {
    message("Start expression quantification for sample ", samples[j], ": ", Sys.time())
    quantdir <- paste(samples[j], "_quant", sep="")
    samp.pos <- grep(paste(paste(samples[j], "_", sep="")), paired.left)
    samp.paired.l <- paste(readdir, paired.left[samp.pos], sep="")
    samp.paired.r <- paste(readdir, paired.right[samp.pos], sep="")
    samp.unpaired.l <- paste(readdir, unpaired.left[samp.pos], sep="")
    samp.unpaired.r <- paste(readdir, unpaired.right[samp.pos], sep="")
    sailfishcmd <- paste("sailfish quant -i ../results/trinity-full/sailfish-index -o ../results/trinity-full/sailfish-expression/", quantdir, " --reads ", samp.paired.l, " ", samp.paired.r, " ", samp.unpaired.l, " ", samp.unpaired.r, " -p 4", sep="")
    print(sailfishcmd)
#    system(sailfishcmd)
    message("Done with expression quantification for sample ", samples[j], ": ", Sys.time(), "\n")
}
```

This generated a directory for each sample

`r list.files("../results/trinity-full/sailfish-expression")`

and within each directory there are the following r:

`r list.files("../results/trinity-full/sailfish-expression/A22-0_quant")`

The file *quant_bias_corrected.sf* contains the following columns, following a number of header lines:

1. Transcript ID
2. Transcript Length
3. Transcripts per Million (TPM): computed as described in `r citep("10.1093/bioinformatics/btp692")`, and is meant as an estimate of the number of transcripts, per million observed transcripts, originating from each isoform.
4. Reads Per Kilobase per Million mapped reads (RPKM): classic measure of relative transcript abundance, and is an estimate of the number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

The TPM column for each sample was extracted and combined into a matrix for each colony.

```{r load_expression_data, eval=TRUE, echo=FALSE, results='hide'}
    
#  read in each file using loop
samples <- c("A22-0", "A22-3", "A22-7", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", "A22-31", "A22-35", "A22-38", "Ar-0", "Ar-3", "Ar-7", "Ar-10", "Ar-14", "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38")

fileinpath <- "../results/trinity-full/sailfish-expression/"

for (j in 1:length(samples)) {
    samp <- samples[j]
    trtval <- as.numeric(str_split_fixed(samp, "-", 2)[2])
    outpre <- gsub("-", "_", samp)
    outname <- paste(outpre, "_quant", sep="")
    read.sailfish.quant(filein=paste(fileinpath, samp, "_quant/quant_bias_corrected.sf", sep=""), outname=outname, samp = samp, trtval = trtval)
}

# combine into long format

A22.TPM <- rbind(A22_0_quant, A22_3_quant, A22_7_quant, A22_10_quant, A22_14_quant, A22_17_quant, A22_21_quant, A22_24_quant, A22_28_quant, A22_31_quant, A22_35_quant, A22_38_quant)
str(A22.TPM)

# convert to data.table

A22.TPM.dt <- data.table(A22.TPM)
setkey(A22.TPM.dt, Transcript, val)

# set "trt" to true values - truncated in file names for convenience

A22.TPM.dt[val==3,val:=3.5]
A22.TPM.dt[val==10,val:=10.5]
A22.TPM.dt[val==17,val:=17.5]
A22.TPM.dt[val==24,val:=24.5]
A22.TPM.dt[val==31,val:=31.5]
A22.TPM.dt[val==38,val:=38.5]
head(A22.TPM.dt)
str(A22.TPM.dt)


# Repeat for Ar
Ar.TPM <- rbind(Ar_0_quant, Ar_3_quant, Ar_7_quant, Ar_10_quant, Ar_14_quant, Ar_17_quant, Ar_21_quant, Ar_24_quant, Ar_28_quant, Ar_31_quant, Ar_35_quant, Ar_38_quant)
Ar.TPM.dt <- data.table(Ar.TPM)
setkey(Ar.TPM.dt, Transcript, val)
Ar.TPM.dt[val==3,val:=3.5]
Ar.TPM.dt[val==10,val:=10.5]
Ar.TPM.dt[val==17,val:=17.5]
Ar.TPM.dt[val==24,val:=24.5]
Ar.TPM.dt[val==31,val:=31.5]
Ar.TPM.dt[val==38,val:=38.5]
head(Ar.TPM.dt)
str(Ar.TPM.dt)

tables()
```

Note that expression levels at each temperature treatment are highly correlated between the two colonies.

```{r exp_correlations, echo=FALSE, results='asis'}
temp <- c(0, 3.5, 7, 10.5, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)
cors <- c(round(cor(Ar_0_quant$TPM, A22_0_quant$TPM), 2), round(cor(Ar_3_quant$TPM, A22_3_quant$TPM), 2), round(cor(Ar_7_quant$TPM, A22_7_quant$TPM), 2), round(cor(Ar_10_quant$TPM, A22_10_quant$TPM), 2), round(cor(Ar_14_quant$TPM, A22_14_quant$TPM), 2), round(cor(Ar_17_quant$TPM, A22_17_quant$TPM), 2), round(cor(Ar_21_quant$TPM, A22_21_quant$TPM), 2), round(cor(Ar_24_quant$TPM, A22_24_quant$TPM), 2), round(cor(Ar_28_quant$TPM, A22_28_quant$TPM), 2), round(cor(Ar_31_quant$TPM, A22_31_quant$TPM), 2), round(cor(Ar_35_quant$TPM, A22_35_quant$TPM), 2), round(cor(Ar_38_quant$TPM, A22_38_quant$TPM), 2))

cortable <- cbind(temp, cors)

pandoc.table(cortable, style="rmarkdown", caption = "correlations between colonies at each temperature treatment")
```

Expression levels are *highly* correlated between colonies. 


## Identification of thermally-responsive genes

To identify transcripts (roughly equivalent to genes) that show thermal responsiveness, I fit the following linear model to each transcript:

$$ TPM = \beta_0 + \beta_1(colony) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(colony * temp) + \\beta_5(colony * temp^2) + \epsilon $$

where TPM is transcripts per million. 

For this list of P-values, False Discovery Rate (FDR) is applied and q-values are calculated using the [qvalue]() package.

Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the analysis to be conservative). 

```{r RxN, echo=TRUE, eval=TRUE}
A22.TPM.dt[,colony:="A22"]
Ar.TPM.dt[,colony:="Ar"]
TPM.dt <- rbind(A22.TPM.dt, Ar.TPM.dt)
TPM.dt$colony <- as.factor(TPM.dt$colony)
str(TPM.dt)

setkey(TPM.dt, val)
TPM.dt.sub <- TPM.dt[val != 7] 
unique(TPM.dt.sub$val)

# define model for RxN function
model <- "TPM ~ colony + val + I(val^2) + colony:val + colony:I(val^2)"

# identify responsive transcripts
RxNout <- RxNseq(f = TPM.dt.sub, model = model)

save.image("RxN_combined_results.RData")
```

Of the `r nrow(RxNout)` transcripts, `r length(which(RxNout$pval < 0.05))` have models with P < 0.05.

Many of these are likely false positives, so I adjusted P-values using false discovery rate (FDR) to identify only those transcripts with less than 5% FDR as significant. 

```{r fdr}
RxNout$qval <- p.adjust(RxNout$pval, method = "fdr")

# subset to significant transcripts
signif.transcripts <- RxNout[which(RxNout$qval < 0.05), ]
```

## Functional annotation

```{r annotation}
# add annotation information
setkey(annotationtable, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
```

```{r table_results, echo=FALSE, results='asis'}
clist <- c('Total', 'Colony', 'Temp.lin', 'Temp.quad', 'Colony:Temp.lin', 'Colony:Temp.quad')
siglist <- c(nrow(signif.transcripts), length(which(!is.na(signif.transcripts$coef.colony))), length(which(!is.na(signif.transcripts$coef.val))), length(which(!is.na(signif.transcripts$'coef.I(val^2)'))), length(which(!is.na(signif.transcripts$'coef.colony:val'))), length(which(!is.na(signif.transcripts$'coef.colony:I(val^2)'))))

sigtable <- data.frame(Coefficient = clist, Number_signficant = siglist)

pandoc.table(sigtable, style="rmarkdown", caption = "Number of transcripts for which each term is significant")
```

Of these, subset to those that have significant responses to temperature, either through a direct effect or interaction with colony. Add annotation information and write results to file. Do the same for transcripts that differ in expression between the colonies.

```{r responsive_transcripts}
responsive.transcripts <- signif.transcripts[!is.na(signif.transcripts$'coef.val') | !is.na(signif.transcripts$'coef.I(val^2)') | !is.na(signif.transcripts$'coef.colony:val') | !is.na(signif.transcripts$'coef.colony:I(val^2)')]
dim(responsive.transcripts)

# join signif transcripts with annotation
responsive.transcripts.ann <- annotationtable[responsive.transcripts]
str(responsive.transcripts.ann)

write.table(responsive.transcripts.ann, file = "Ap_responsive_transcripts_GO.txt", quote = FALSE, sep = "\t", row.names = FALSE)

colony.transcripts <- signif.transcripts[!is.na(signif.transcripts$'coef.colony')]
colony.transcripts.ann <- annotationtable[colony.transcripts]
str(colony.transcripts.ann)
write.table(colony.transcripts.ann, file = "Ap_colony_transcripts_GO.txt", quote = FALSE, sep = "\t", row.names = FALSE)
```

Encouragingly, 2 known candidates (e.g. heat shock proteins?) occur in this list:

`r responsive.transcripts.ann[grep("heat shock protein", responsive.transcripts.ann$best.hit.to.nr), list(best.hit.to.nr)]`


The following table shows the number of responsive transcripts that have maxmimum expression at low, intermediate and high temperatures. 

```{r maxmin, eval=TRUE, echo=FALSE, results = 'asis'}
# merge RxN results with expression values
Ap.dt <- TPM.dt.sub[responsive.transcripts.ann]

# for each transcript, identify temperature of max expression
Ap.max.min.exp <- ddply(Ap.dt, .(Transcript), function(df) {
    lmout <- lm(TPM ~ val + I(val^2), data = df)
    vals <- c(0, 3.5, 10, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)
    newdf <- data.frame(colony = rep("A22", 11), val = vals)
    pout <- predict(lmout, newdata=newdf)
    return(c(max.val = vals[which(pout == max(pout))],
             min.val = vals[which(pout == min(pout))]))
    }
)

# Table maximum and minium expression by temperature
maxt <- table(Ap.max.min.exp$max.val)
mint <- table(Ap.max.min.exp$min.val)

expt <- rbind(Max = maxt, Min = mint)

pandoc.table(expt, style="rmarkdown", caption = "Number of transcripts that have maximum and minimum expression at each temperature level")
```

**Most** transcripts have maximum expression at high (>=31.5) or low (<=10) temperatures. As the genes invovled in each type of thermal response may differ, designate responsive transcripts by max expression at high, low or intermediate (14-28) temperatures.

```{r, eval=FALSE}
# set expression type
responsive.transcripts <- responsive.transcripts[Ap.max.min.exp]
responsive.transcripts[ , exp_type:="intermediate"]
responsive.transcripts[max.val>31, exp_type:="high"]
responsive.transcripts[max.val<11, exp_type:="low"]
```

## Functional annotation ##

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  


## Gene set enrichment analysis ##

Get annotation of responsive genes. 

Use [topGO](http://www.bioconductor.org/packages/2.12/bioc/html/topGO.html) to perform gene set enrichment analysis

First need to create gene ID to GO term map file

```{r geneid2go, echo=TRUE, eval=FALSE}

# create geneid2go.map file from FastAnnotator AnnotationTable.txt
geneid2GOmap(annotationfile)
```

Using this gene2GO map file, perform gene set enrichment analysis.

Provide qvalues as the gene score, and select genes with q < 0.05 using custom `selectFDR` function.

```{r gene_enrichment_A22, echo=FALSE, eval=FALSE}

# read mappings file
geneID2GO <- readMappings(file = "geneid2go.map")
str(head(geneID2GO))

# create geneList. note that NA values cause problems with topGO
# need to retain these for full ontology, so set NA to 1
A22geneList <- A22_RxN$pval
A22geneList[which(is.na(A22geneList))] <- 1
names(A22geneList) <- A22_RxN$transcript
str(A22geneList)

# Function to select top genes (defined above)
selectFDR <- function(pvalue) {
    return(pvalue < 0.033196)
}

# create topGOdata object
A22.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis", ontology = "BP",
                 allGenes = A22geneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

A22.BP.GOdata

# perform enrichment analysis using multiple methods
A22.BP.resultParentChild <- runTest(A22.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
A22.BP.resultParentChild

# 8 significant GO terms

A22.BP.ResTable <- GenTable(A22.BP.GOdata, parentchild = A22.BP.resultParentChild, topNodes = 16)
A22.BP.ResTable

# graph significant nodes

pdf("A22.BP_topGO_sig_nodes.pdf")
showSigOfNodes(A22.BP.GOdata, score(A22.BP.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()


## Cellular Component

# create topGOdata object
A22.CC.GOdata <- new("topGOdata",
                 description = "CC gene set analysis", ontology = "CC",
                 allGenes = A22geneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

A22.CC.GOdata

# perform enrichment analysis using multiple methods
A22.CC.resultParentChild <- runTest(A22.CC.GOdata, statistic = 'fisher', algorithm = 'parentchild')
A22.CC.resultParentChild

# 2 significant GO terms

CC.ResTable <- GenTable(A22.CC.GOdata, parentchild = A22.CC.resultParentChild, topNodes = 10)
CC.ResTable

# graph significant nodes

pdf("A22.CC_topGO_sig_nodes.pdf")
showSigOfNodes(A22.CC.GOdata, score(A22.CC.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()


#### Molecular Function

# create topGOdata object
A22.MF.GOdata <- new("topGOdata",
                 description = "MF gene set analysis", ontology = "MF",
                 allGenes = A22geneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

A22.MF.GOdata

# perform enrichment analysis using multiple methods
A22.MF.resultParentChild <- runTest(A22.MF.GOdata, statistic = 'fisher', algorithm = 'parentchild')
A22.MF.resultParentChild

# 6 significant GO terms

A22.MF.ResTable <- GenTable(A22.MF.GOdata, parentchild = A22.MF.resultParentChild, topNodes = 12)
A22.MF.ResTable

# graph significant nodes

pdf("A22.MF_topGO_sig_nodes.pdf")
showSigOfNodes(A22.MF.GOdata, score(A22.MF.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()

```


```{r}
# Plot
Ap.high <- data.table(Ap.high)
setkey(Ap.high, Transcript)
Ap.high.dt <- Ap.plot.dt[Ap.high]
Ap.low <- data.table(Ap.low)
setkey(Ap.low, Transcript)
Ap.low.dt <- Ap.plot.dt[Ap.low]


# combine
Ap.plot.dt <- rbind(Ap.high.dt, Ap.low.dt)


# scale expression values 
Ap.plot.dt[,exp.scaled:=scale(TPM), by = Transcript]
dim(Ap.plot.dt)

# plot

png("Ap_expression_ggplot.png")
  p <- ggplot(Ap.plot.dt, aes(x=val, y=exp.scaled, group=Transcript)) + geom_line()
  p + facet_grid(. ~ colony)
#  p + geom_line(aes(colour = max.val)) + scale_colour_gradient(low="red")
dev.off()

pdf("Ap_high_expression_smooth_ggplot.pdf")
  p <- ggplot(Ap.high.dt, aes(x=val, y=exp.scaled, group=transcript))
  p + geom_smooth(aes(colour = colony))
  #p + geom_smooth(aes(colour = lin.coef)) + scale_colour_gradient(low="red")
dev.off()
```






```{r hists, echo=FALSE, eval=FALSE}

## IN PROGRESS ##

# combine data
# reshape long
# make histogram

png("hist_'(Intercept)'.png")
  h <- ggplot(A22.RxN.G.qcrit, aes('(Intercept)')) + geom_histogram(binwidth = 0.5) +
#   scale_x_log10()
    xlim(0, 100)
  h
dev.off()

png("hist_linear.png")
  h <- ggplot(A22.RxN.G.qcrit, aes(lin.coef)) + geom_histogram(binwidth = 0.5) +
    xlim(-5,5)
  h
dev.off()

png("hist_quad.png")
  h <- ggplot(A22.RxN.G.qcrit, aes(quad.coef)) + geom_histogram(binwidth = 0.1) +
    xlim(-1,1)
  h
dev.off()
```



Perform gene enrichment analysis for Ar. As there were no genes significant at q < 0.05, use 0.1 as the critical threshold.

```{r gene_enrichment_Ar, echo=FALSE, eval=FALSE}

# create geneList. note that NA values cause problems with topGO
# need to retain these for full ontology, so set NA to 1
ArgeneList <- Ar_RxN$pval
ArgeneList[which(is.na(ArgeneList))] <- 1
names(ArgeneList) <- Ar_RxN$transcript
str(ArgeneList)

# create function to select significant genes for topGO
selectFDR1 <- function(pvalue) {
    return(pvalue < 0.010232)
}


# create topGOdata object
Ar.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis", ontology = "BP",
                 allGenes = ArgeneList, geneSel = selectFDR1,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

Ar.BP.GOdata

# perform enrichment analysis using multiple methods
Ar.BP.resultParentChild <- runTest(Ar.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ar.BP.resultParentChild

Ar.BP.ResTable <- GenTable(Ar.BP.GOdata, parentchild = Ar.BP.resultParentChild, topNodes = 10)
Ar.BP.ResTable

# graph significant nodes

pdf("Ar.BP_topGO_sig_nodes.pdf")
showSigOfNodes(Ar.BP.GOdata, score(Ar.BP.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()

## Cellular Component

# create topGOdata object
Ar.CC.GOdata <- new("topGOdata",
                 description = "CC gene set analysis", ontology = "CC",
                 allGenes = ArgeneList, geneSel = selectFDR1,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

Ar.CC.GOdata

# perform enrichment analysis using multiple methods
Ar.CC.resultParentChild <- runTest(Ar.CC.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ar.CC.resultParentChild

Ar.CC.ResTable <- GenTable(Ar.CC.GOdata, parentchild = Ar.CC.resultParentChild, topNodes = 10)
Ar.CC.ResTable

# graph significant nodes

pdf("Ar.CC_topGO_sig_nodes.pdf")
showSigOfNodes(Ar.CC.GOdata, score(Ar.CC.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()


#### Molecular Function

# create topGOdata object
Ar.MF.GOdata <- new("topGOdata",
                 description = "MF gene set analysis", ontology = "MF",
                 allGenes = ArgeneList, geneSel = selectFDR1,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

Ar.MF.GOdata

# perform enrichment analysis using multiple methods
Ar.MF.resultParentChild <- runTest(Ar.MF.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ar.MF.resultParentChild

Ar.MF.ResTable <- GenTable(Ar.MF.GOdata, parentchild = Ar.MF.resultParentChild, topNodes = 10)
Ar.MF.ResTable

# graph significant nodes

pdf("Ar.MF_topGO_sig_nodes.pdf")
showSigOfNodes(Ar.MF.GOdata, score(Ar.MF.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()
```


## Session information ##

```{r session}
save.image()
sessionInfo()
```
