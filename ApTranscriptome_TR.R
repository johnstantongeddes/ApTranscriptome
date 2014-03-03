Thermal reactionome of the common ant species *Aphaenogaster*
================================================================
  
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**February 28, 2014**

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
#biocLite("Rgraphviz")
library(topGO) 
suppressMessages(library(Rgraphviz))

# load custom functions
source("scripts/RxNseq.R")

# knitr options
opts_chunk$set(cache=TRUE)

# directory to save results
resultsdir <- "results/"
```


## Summary ##
  
In this technical report, which accompanies the manuscript **Thermal reactionome of a common ant species** (Stanton-Geddes et al., in press), we:

1. describe the *de novo* assembly of the transcriptome for two ant colonies with in the *Aphaenogaster rudis-picea-fulva* species complex `r citep("10.1155/2012/752815")`
2. identify thermally-responsive genes
3. perform gene set enrichment analysis

This script is completely reproducible assuming that R, `knitr` and the other required libraries (listed within the source document) are installed on a standard linux system using the following:
    
    Rscript -e "library(knitr); knit('ApTranscriptome_TR.Rmd')"

The assembled transcriptome, annotation and expression values are downloaded rather than re-run due to the computational demands, but the exact commands for each of these steps are documented below.



## Data ##

The Illumina fastq files are available from [https://minilims1.uvm.edu/BCProject-26-Cahan/downloads.html]. The Trimmomatic filtered fastq files should be downloaded, uncompressed and moved to the appropriate directory using the following code.

~~~
# download filtered fastq files, uncompress and move
wget --no-check-certificate https://minilims1.uvm.edu/BCProject-26-Cahan/_downloads/trimmomatic_output.tar.gz
tar -zxvf trimmomatic_output.tar.gz
mkdir -p data/
mv ind_files data/.
~~~



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

Running Trinity is time and computationally-intensive, so the final filtered assembly can be downloaded from [http://johnstantongeddes.org/assets/files/Trinity_filtered_assembly.tgz] 

~~~
# download filtered Trinity assembly, uncompress and move
wget http://johnstantongeddes.org/assets/files/Trinity_filtered_assembly.tgz
tar -zxvf trimmomatic_output.tar.gz
mkdir -p results/
mkdir -p results/trinity-full/
mv Trinity_cap3_uclust.fasta /results/trinity-full/.
~~~



## Transcriptome annotation ##

Annotation was performed by uploading the reduced assembly "Trinity_cap3_uclust.fasta" to the web-based annotation program [FastAnnotator](http://fastannotator.cgu.edu.tw/index.php) `r citep("10.1186/1471-2164-13-S7-S9")`.

#Results are available as job ID [13894410176993](http://fastannotator.cgu.edu.tw/job.php?jobid=13894410176993#page=basicinfo).

This annotation file can be read directly to R:

```{r download, echo = TRUE, results = "hide"}
### Annotation file
# from either AWS or GoogleDrive
annotationURL <- getURL("http://johnstantongeddes.org/assets/files/Aphaeno_transcriptome_AnnotationTable.txt")
#a2 <- getURL("https://googledrive.com/host/0B75IymziRJ_9Tlg1U1Vxbjk1bzg") # GoogleDrive link

annotationfile <- read.csv(textConnection(annotationURL), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(annotationfile)

# Convert to data.table
annotationtable <- data.table(annotationfile)
head(annotationtable)
```

Transcriptome annotations are nearly impossible to visualize in a meaningful way. For lack of better ideas, I created a word cloud.

```{r wordcloud, cache=TRUE}
library(tm)
library(wordcloud)
library(RColorBrewer)
library(stringr)

# subset and extract unique best hits to NCBI nr database
ann <- unique(annotationfile$best.hit.to.nr)
annsplit <- str_split_fixed(ann, pattern = " ", n=2)
annvec <- annsplit[,2]

# create Corpus of terms. convert to dataframe for plotting wordcloud
ap.corpus <- Corpus(VectorSource(annvec))
ap.corpus <- tm_map(ap.corpus, removePunctuation)
ap.corpus <- tm_map(ap.corpus, tolower)
ap.corpus <- tm_map(ap.corpus, function(x) removeWords(x, stopwords("english")))
tdm <- TermDocumentMatrix(ap.corpus)
m <- as.matrix(tdm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
pal <- brewer.pal(4, "Dark2")
pal <- pal[-(1:2)]

#png("wordcloud.png", width=1280,height=800)
wordcloud(d$word,d$freq, scale=c(8,.3),min.freq=2,max.words=100, random.order=T, rot.per=.15, colors=pal, vfont=c("sans serif","plain"))
#dev.off()
```

This wordcloud shows that we do not know much about the transcriptome - mostly hypothetical or predicted proteins. I removed these uninformative words and generated a new wordcloud. 

```{r wordcloud2, cache=TRUE}
# remove 'protein', 'predicted' and 'hypothetical'
rmwords <- c("hypothetical", "protein", "isoform", "subunit", "family", "putative", "partial", "domain", "containing", "predicted", "domaincontaining", "uncharacterized")
ap.corpus.rmwords <- tm_map(ap.corpus, removeWords, rmwords)
tdm <- TermDocumentMatrix(ap.corpus.rmwords)
m <- as.matrix(tdm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
pal <- brewer.pal(8, "Dark2")
pal <- pal[-(1:2)]
#png("wordcloud2.png", width=1280,height=800)
wordcloud(d$word,d$freq, scale=c(8,.3),min.freq=2,max.words=100, random.order=T, rot.per=.15, colors=pal, vfont=c("sans serif","plain"))
#dev.off()
```
                                                 
                                                
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
readdir <- "data/ind_files/" 
# sailfish index directory
sfindex <- "results/trinity-full/sailfish-index"
# sailfish expression directory
sfexpressionroot <- "results/trinity-full/sailfish-expression/"

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
    quantdir <- paste(sfexpressionroot, samples[j], "_quant", sep="")
    samp.pos <- grep(paste(paste(samples[j], "_", sep="")), paired.left)
    samp.paired.l <- paste(readdir, paired.left[samp.pos], sep="")
    samp.paired.r <- paste(readdir, paired.right[samp.pos], sep="")
    samp.unpaired.l <- paste(readdir, unpaired.left[samp.pos], sep="")
    samp.unpaired.r <- paste(readdir, unpaired.right[samp.pos], sep="")
    sailfishcmd <- paste("sailfish quant -i ", sfindex, " -o ", quantdir, " --reads ", samp.paired.l, " ", samp.paired.r, " ", samp.unpaired.l, " ", samp.unpaired.r, " -p 4", sep="")
#   print(sailfishcmd)
    system(sailfishcmd)
    message("Done with expression quantification for sample ", samples[j], ": ", Sys.time(), "\n")
}
```

This generated a directory for each sample

`r list.files(sfexpressionroot)`

and within each directory there are the following r:

`r list.files(quantdir)`

The file *quant_bias_corrected.sf* contains the following columns, following a number of header lines:

1. Transcript ID
2. Transcript Length
3. Transcripts per Million (TPM): computed as described in `r citep("10.1093/bioinformatics/btp692")`, and is meant as an estimate of the number of transcripts, per million observed transcripts, originating from each isoform.
4. Reads Per Kilobase per Million mapped reads (RPKM): classic measure of relative transcript abundance, and is an estimate of the number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

The TPM column for each sample was extracted and combined into a matrix for each colony.

```{r load_expression_data, eval=TRUE, echo=FALSE, results='hide'}
# read in each file using loop
samples <- c("A22-0", "A22-3", "A22-7", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", "A22-31", "A22-35", "A22-38", "Ar-0", "Ar-3", "Ar-7", "Ar-10", "Ar-14", "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38")

for (j in 1:length(samples)) {
    samp <- samples[j]
    trtval <- as.numeric(str_split_fixed(samp, "-", 2)[2])
    outpre <- gsub("-", "_", samp)
    outname <- paste(outpre, "_quant", sep="")
    read.sailfish.quant(filein=paste(sfexpressionroot, samp, "_quant/quant_bias_corrected.sf", sep=""), outname=outname, samp = samp, trtval = trtval)
}

# combine into long format
A22.TPM <- rbind(A22_0_quant, A22_3_quant, A22_7_quant, A22_10_quant, A22_14_quant, A22_17_quant, A22_21_quant, A22_24_quant, A22_28_quant, A22_31_quant, A22_35_quant, A22_38_quant)
str(A22.TPM)

# convert to data.table
A22.TPM <- data.table(A22.TPM)
setkey(A22.TPM, Transcript, val)

# set "trt" to true values - truncated in file names for convenience

A22.TPM[val==3,val:=3.5]
A22.TPM[val==10,val:=10.5]
A22.TPM[val==17,val:=17.5]
A22.TPM[val==24,val:=24.5]
A22.TPM[val==31,val:=31.5]
A22.TPM[val==38,val:=38.5]
head(A22.TPM)
str(A22.TPM)


# Repeat for Ar
Ar.TPM <- rbind(Ar_0_quant, Ar_3_quant, Ar_7_quant, Ar_10_quant, Ar_14_quant, Ar_17_quant, Ar_21_quant, Ar_24_quant, Ar_28_quant, Ar_31_quant, Ar_35_quant, Ar_38_quant)
Ar.TPM <- data.table(Ar.TPM)
setkey(Ar.TPM, Transcript, val)
Ar.TPM[val==3,val:=3.5]
Ar.TPM[val==10,val:=10.5]
Ar.TPM[val==17,val:=17.5]
Ar.TPM[val==24,val:=24.5]
Ar.TPM[val==31,val:=31.5]
Ar.TPM[val==38,val:=38.5]
head(Ar.TPM)
str(Ar.TPM)

tables()
```

Note that expression levels at each temperature treatment are highly correlated between the two colonies.

```{r exp_correlations, echo=FALSE, results='asis'}
# Make Rmarkdown pandoc table of correlations among expression levels between colonies at each temp

temp <- c(0, 3.5, 7, 10.5, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)
cors <- c(round(cor(Ar_0_quant$TPM, A22_0_quant$TPM), 2), round(cor(Ar_3_quant$TPM, A22_3_quant$TPM), 2), round(cor(Ar_7_quant$TPM, A22_7_quant$TPM), 2), round(cor(Ar_10_quant$TPM, A22_10_quant$TPM), 2), round(cor(Ar_14_quant$TPM, A22_14_quant$TPM), 2), round(cor(Ar_17_quant$TPM, A22_17_quant$TPM), 2), round(cor(Ar_21_quant$TPM, A22_21_quant$TPM), 2), round(cor(Ar_24_quant$TPM, A22_24_quant$TPM), 2), round(cor(Ar_28_quant$TPM, A22_28_quant$TPM), 2), round(cor(Ar_31_quant$TPM, A22_31_quant$TPM), 2), round(cor(Ar_35_quant$TPM, A22_35_quant$TPM), 2), round(cor(Ar_38_quant$TPM, A22_38_quant$TPM), 2))

cortable <- cbind(temp, cors)

pandoc.table(cortable, style="rmarkdown", caption = "correlations between colonies at each temperature treatment")
```

## Identification of thermally-responsive genes

To identify transcripts (roughly equivalent to genes) that show thermal responsiveness, I fit the following linear model to each transcript:

$$ TPM = \beta_0 + \beta_1(colony) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(colony * temp) + \\beta_5(colony * temp^2) + \epsilon $$

where TPM is transcripts per million. 

For this list of P-values, correct for multiple testing using False Discovery Rate (FDR).

#Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the analysis to be conservative). 

```{r RxN, echo=TRUE, eval=TRUE, cache=TRUE}
A22.TPM[,colony:="A22"]
Ar.TPM[,colony:="Ar"]
TPM.dt <- rbind(A22.TPM, Ar.TPM)
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
RxNout$padj <- p.adjust(RxNout$pval, method = "fdr")
# Plot FDR values against initial pvalues
par(mfrow = c(2,1))
hist(RxNout$pval)
hist(RxNout$padj)

# subset to significant transcripts
signif.transcripts <- RxNout[which(RxNout$padj < 0.05), ]
```

## Functional annotation ##

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  


```{r annotation}
# add annotation information
setkey(annotationtable, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
```

```{r table_results, echo=FALSE, results='asis'}
clist <- c('Total', 'Colony', 'Temp.lin', 'Temp.quad', 'Colony:Temp.lin', 'Colony:Temp.quad')
siglist <- c(nrow(signif.transcripts), length(which(!is.na(signif.transcripts$coef.colony))), length(which(!is.na(signif.transcripts$coef.val))), length(which(!is.na(signif.transcripts$'coef.I(val^2)'))), length(which(!is.na(signif.transcripts$'coef.colony:val'))), length(which(!is.na(signif.transcripts$'coef.colony:I(val^2)'))))

sigtable <- data.frame(Coefficient = clist, Number_significant = siglist)

pandoc.table(sigtable, style="rmarkdown", caption = "Number of transcripts for which each term is significant")
```

Of these, subset to those that have significant responses to temperature, either through a direct effect or interaction with colony. Add annotation information and write results to file. Do the same for transcripts that differ in expression between the colonies.

```{r responsive_transcripts}
responsive.transcripts <- signif.transcripts[!is.na(signif.transcripts$'coef.val') | !is.na(signif.transcripts$'coef.I(val^2)') | !is.na(signif.transcripts$'coef.colony:val') | !is.na(signif.transcripts$'coef.colony:I(val^2)')]
dim(responsive.transcripts)
# join signif transcripts with annotation
responsive.transcripts <- annotationtable[responsive.transcripts]
str(responsive.transcripts)

# transcripts that differ in expression by colony
colony.transcripts <- signif.transcripts[!is.na(signif.transcripts$'coef.colony')]
colony.transcripts <- annotationtable[colony.transcripts]
str(colony.transcripts)
colony.transcripts <- colony.transcripts[order(colony.transcripts$padj), ]
write.table(colony.transcripts, file = paste(resultsdir, "Ap_colony_transcripts_GO.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

# transcripts that have colony by temperature interactions
interaction.transcripts <- signif.transcripts[!is.na(signif.transcripts$'coef.colony:val') | !is.na(signif.transcripts$'coef.colony:I(val^2)')]
interaction.transcripts <- annotationtable[interaction.transcripts]
str(interaction.transcripts)
```

For responsive transcripts, identify by shape of response:

* High - increase expression with temperature
* Low - decrease expression with temperature
* Intermediate - maximum expression at intermediate temperatures (14 - 28C)
* Bimodal - expressed greater than two standard deviations of expression at both low and high temperatures

```{r expression_shape, eval=TRUE, echo=TRUE}
# merge RxN results with expression values
setkey(TPM.dt.sub, Transcript)
Ap.dt <- TPM.dt.sub[responsive.transcripts]
setkey(Ap.dt, Transcript)

Ap.exp_type <- ddply(Ap.dt, .(Transcript), function(df1) {
    lmout <- lm(TPM ~ val + I(val^2), data = df1)
    vals <- c(0, 3.5, 10, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)
    newdf <- data.frame(val = vals)
    pout <- predict(lmout, newdata=newdf)
    pout <- data.frame(val = vals, exp = pout)

    # get vals of max and min expression
    max.val = vals[which(pout$exp == max(pout$exp))]
    min.val = vals[which(pout$exp == min(pout$exp))]

    # report coefficients
    #coef(lmout)
    exp_type = if(coef(lmout)['val'] > 0 & coef(lmout)['I(val^2)'] > 0) "High" else {
        if(coef(lmout)['val'] < 0 & coef(lmout)['I(val^2)'] < 0) "Low" else {
            if(coef(lmout)['val'] > 0 & coef(lmout)['I(val^2)'] < 0) "Intermediate" else {
                "convex"}}}

    # for transcripts with convex exp_type, check if expression is truly bimodal
    if(exp_type == "convex") {
        if(max(pout[pout$val <= 10, "exp"]) > 2*sd(pout$exp) &
           max(pout[pout$val >= 31.5, "exp"]) > 2*sd(pout$exp)) exp_type = "Bimodal" else {
           # linear increase?
               if(max.val > min.val) exp_type = "High" else exp_type = "Low"
           }
    }
           
    # return values
    return(c(max.val = vals[which(pout$exp == max(pout$exp))],
             min.val = vals[which(pout$exp == min(pout$exp))],
             exp_type = exp_type))
    }
 )

# merge 'exp_type' information with Ap.dt
Ap.exp_type <- data.table(Ap.exp_type)
setkey(Ap.exp_type, Transcript)
Ap.dt <- Ap.dt[Ap.exp_type]

# merge 'exp_type' with responsive.transcripts
setkey(responsive.transcripts, Sequence.Name)
responsive.transcripts <- Ap.exp_type[responsive.transcripts]
dim(responsive.transcripts)
# order by padj
responsive.transcripts <- responsive.transcripts[order(responsive.transcripts$padj), ]
# save results to file
write.table(responsive.transcripts, file = paste(resultsdir, "Ap_responsive_transcripts_GO.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

# merge 'exp_type' with interaction transcripts
setkey(interaction.transcripts, Sequence.Name)
interaction.transcripts <- Ap.exp_type[interaction.transcripts]
dim(interaction.transcripts)
# order by padj
interaction.transcripts <- interaction.transcripts[order(interaction.transcripts$padj), ]
# save results to file
write.table(interaction.transcripts, file = paste(resultsdir, "Ap_interaction_transcripts_GO.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
```

```{r expression_exp_type_table, eval=TRUE, echo=FALSE, results='asis'}
# Table type of response
exp_type.table <- table(Ap.exp_type$exp_type)
pandoc.table(exp_type.table, style="rmarkdown", caption = "Number of transcripts with maximum expression at high, low, intermediate or both high and low (bimodal) temperatures.")
```

Note that among responsive transcripts, there are 25 transcripts with GO term "response to stress" and various heat shock related proteins:

```{r}
unique(Ap.dt[grep("GO:0006950", Ap.dt$GO.Biological.Process), list(Transcript, best.hit.to.nr)])
unique(Ap.dt[grep("shock", Ap.dt$best.hit.to.nr), list(Transcript, best.hit.to.nr)])

unique(Ap.dt.interaction[grep("shock", Ap.dt.interaction$best.hit.to.nr), list(Transcript, best.hit.to.nr)])
```                                                    


Export data for interactive shiny app. 

```{r shiny_file}
# scale expression values 
Ap.dt[,exp.scaled:=scale(TPM), by = Transcript]
str(Ap.dt)
write.csv(Ap.dt, file = paste(resultsdir, "Ap.dt.csv", sep=""), quote = TRUE, row.names = FALSE)

# subset to genes with significant interaction
Ap.dt.interaction <- Ap.dt[!is.na(Ap.dt$'coef.colony:val') | !is.na(Ap.dt$'coef.colony:I(val^2)')]
str(Ap.dt.interaction)
write.csv(Ap.dt.interaction, file = paste(resultsdir, "Ap.dt.interaction.csv", sep=""), quote = TRUE, row.names = FALSE)
```


## Colony-level comparison ##

Compare the expression levels at optimum (21C) between the two colonies for genes in each expression group. Specifically, are genes that are upregulated at high temperatures in A22 more highly expressed at 21C in Ar? Conversely, are genes that are upregulated at low temps more highly expressed at 21C in A22? 

```{r, cache=TRUE, eval=FALSE}
# for each transcript in interaction group
#     predict 
#     predict expression at 21C for each colony

exp_by_colony <- ddply(Ap.dt.interaction, .(Transcript, colony), RxNply)

# check that correct number of rows are output - 2 times the number of transcripts
stopifnot(all.equal(2*length(unique(Ap.dt.interaction$Transcript)), nrow(exp_by_colony)))                
# change 'opt.exp' to numeric
exp_by_colony$opt.exp <- as.numeric(exp_by_colony$opt.exp)
# scale expression values by Transcript so that large expression transcripts don't drive overall pattern
exp_by_colony <- ddply(exp_by_colony, .(Transcript), transform, 
  opt.exp.scaled = scale(opt.exp))

# list of transcripts that are 'high' expressed in A22
A22_high_transcripts <- exp_by_colony[which(exp_by_colony$colony == "A22" & exp_by_colony$exp_type == "High"), "Transcript"]
# dataframe of transcripts from both colonies that are 'high' expressed in A22
A22_high_df <- exp_by_colony[which(exp_by_colony$Transcript %in% A22_high_transcripts), ]

# Compare expression at optimum temp (19.25C) between colonies using t-test
A22_high_df$colony <- as.factor(A22_high_df$colony)
boxplot(data = A22_high_df, log(opt.exp+1) ~ colony)
t.test(opt.exp ~ colony, data = A22_high_df)

# repeate using scaled expression values so outliers don't drive results
boxplot(data = A22_high_df, opt.exp.scaled ~ colony)
t.test(opt.exp.scaled ~ colony, data = A22_high_df)

# this seems wrong...uses absolute values when I actually just want to compare signs. use Wilcoxon signed rank-test

w1 <- wilcox.test(opt.exp ~ colony, data = A22_high_df, alternative = "two.sided", conf.int = TRUE)
w1
```


Compare the width of the 'stable' region for intermediate-expressed transcripts between the colonies. 

Compare the width of the inverse of the bimodally-expressed transcripts between the colonies.


## Gene set enrichment analysis ##

I use [topGO](http://www.bioconductor.org/packages/2.12/bioc/html/topGO.html) to perform gene set enrichment analysis (GSEA) seperately for each expression type (bimodal, intermediate, high, low).

First need to create gene ID to GO term map file

```{r geneid2go, echo=TRUE, eval=TRUE, cache=TRUE}
# create geneid2go.map file from FastAnnotator AnnotationTable.txt
geneid2GOmap(annotationfile)
```

then read map file.

```{r readmap, echo=TRUE, eval=TRUE}
# read mappings file
geneID2GO <- readMappings(file = "geneid2go.map")
str(head(geneID2GO))
```

### GSEA for thermally-responsive transcripts ###

Using this gene2GO map file, perform GSEA for:

**1) all responsive transcripts**

Use `selectFDR` function to select transcripts with adjusted P < 0.05.

```{r gsea_all}
# create geneList. note that NA values cause problems with topGO
# so set any NA to 1 as need to retain for GO analysis
Ap.geneList <- RxNout$padj
Ap.geneList[which(is.na(Ap.geneList))] <- 1
stopifnot(length(which(is.na(Ap.geneList))) == 0)
names(Ap.geneList) <- RxNout$Transcript
str(Ap.geneList)

# Function to select top genes (defined above)
selectFDR <- function(padj) {
    return(padj < 0.05)
}

# create topGOdata object
Ap.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis", ontology = "BP",
                 allGenes = Ap.geneList, geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

Ap.BP.GOdata

# perform enrichment analysis using parentchild method
Ap.BP.resultParentChild <- runTest(Ap.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ap.BP.resultParentChild

# table results
Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, topNodes = 118)
#Ap.BP.ResTable
write.table(Ap.BP.ResTable, file = paste(resultsdir, "Ap_GO.BP_results.txt", sep=""), quote=FALSE, row.names=FALSE, sep = "\t")
pandoc.table(Ap.BP.ResTable)

# graph significant nodes

pdf(paste(resultsdir, "Ap.BP_topGO_sig_nodes.pdf", sep=""))
showSigOfNodes(Ap.BP.GOdata, score(Ap.BP.resultParentChild), firstSigNodes = 10, useInfo = 'all')
dev.off()
```

**2) High and low expressed transcripts**

Use `selectTranscript` function to select transcripts based on 'exp_type'. As only significant transcripts have an 'exp_type' assigned, this is a small subset of the above.

```{r exp_type_function}
selectTranscript <- function(score) {
    return(score == 1)
}

# To select genes, add 'exp_type' to RxN df
RxNout.t <- merge(RxNout, Ap.exp_type, by="Transcript", all=TRUE)
length(which(!is.na(RxNout.t$exp_type)))
```


```{r bim.gsea, echo=TRUE, eval=TRUE}
# create geneList
Ap.bim.geneList <- rep(0, length = nrow(RxNout.t))
Ap.bim.geneList[which(RxNout.t$exp_type=="Bimodal")] <- 1
names(Ap.bim.geneList) <- RxNout.t$Transcript
str(Ap.bim.geneList)
table(Ap.bim.geneList)

# create topGOdata object
Ap.bim.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis",
                 ontology = "BP",
                 allGenes = Ap.bim.geneList,
                 geneSel = selectTranscript,
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GO)

Ap.bim.BP.GOdata

# perform enrichment analysis using parentchild method
Ap.bim.BP.resultParentChild <- runTest(Ap.bim.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ap.bim.BP.resultParentChild
```

Enriched gene sets for genes expressed at both high and low temperatures.

```{r, echo=FALSE, results="asis"}
Ap.bim.BP.ResTable <- GenTable(Ap.bim.BP.GOdata, parentchild = Ap.bim.BP.resultParentChild, topNodes = 52)
#Ap.bim.BP.ResTable
write.table(Ap.bim.BP.ResTable, file = paste(resultsdir, "Ap_bimodal_GO.BP_results.txt", sep=""), quote=FALSE, row.names=FALSE, sep = "\t")
pandoc.table(Ap.bim.BP.ResTable)
```


**3) Intermediate-temperature expressed transcripts**

```{r int.gsea, echo=TRUE, eval=TRUE}
# create geneList
Ap.int.geneList <- rep(0, length = nrow(RxNout.t))
Ap.int.geneList[which(RxNout.t$exp_type=="Intermediate")] <- 1
names(Ap.int.geneList) <- RxNout.t$Transcript
str(Ap.int.geneList)
table(Ap.int.geneList)

# create topGOdata object
Ap.int.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis",
                 ontology = "BP",
                 allGenes = Ap.int.geneList,
                 geneSel = selectTranscript,
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GO)

Ap.int.BP.GOdata

# perform enrichment analysis using parentchild method
Ap.int.BP.resultParentChild <- runTest(Ap.int.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ap.int.BP.resultParentChild
```

Enriched gene sets for genes expressed at intermediate temperatures.

```{r, echo=FALSE, results="asis"}
Ap.int.BP.ResTable <- GenTable(Ap.int.BP.GOdata, parentchild = Ap.int.BP.resultParentChild, topNodes=16)
Ap.int.BP.ResTable
write.table(Ap.int.BP.ResTable, file = paste(resultsdir, "Ap_intermediate_GO.BP_results.txt", sep=""), quote=FALSE, row.names=FALSE, sep = "\t")
pandoc.table(Ap.int.BP.ResTable)
```


**4) High temperature expressed transcripts**


```{r hig.gsea, echo=TRUE, eval=TRUE}
# create geneList
Ap.hig.geneList <- rep(0, length = nrow(RxNout.t))
Ap.hig.geneList[which(RxNout.t$exp_type=="High")] <- 1
names(Ap.hig.geneList) <- RxNout.t$Transcript
str(Ap.hig.geneList)
table(Ap.hig.geneList)

# create topGOdata object
Ap.hig.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis",
                 ontology = "BP",
                 allGenes = Ap.hig.geneList,
                 geneSel = selectTranscript,
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GO)

Ap.hig.BP.GOdata

# perform enrichment analysis using parentchild method
Ap.hig.BP.resultParentChild <- runTest(Ap.hig.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ap.hig.BP.resultParentChild
```

Enriched gene sets for genes expressed at both high temperatures only.

```{r, echo=FALSE, results="asis"}
Ap.hig.BP.ResTable <- GenTable(Ap.hig.BP.GOdata, parentchild = Ap.hig.BP.resultParentChild, topNodes=12)
#Ap.hig.BP.ResTable
write.table(Ap.hig.BP.ResTable, file = paste(resultsdir, "Ap_high_GO.BP_results.txt", sep=""), quote=FALSE, row.names=FALSE, sep = "\t")
pandoc.table(Ap.hig.BP.ResTable)
```


**5) Low temperature expressed transcripts**

```{r low.gsea, echo=TRUE, eval=TRUE}
# create geneList
Ap.low.geneList <- rep(0, length = nrow(RxNout.t))
Ap.low.geneList[which(RxNout.t$exp_type=="Low")] <- 1
names(Ap.low.geneList) <- RxNout.t$Transcript
str(Ap.low.geneList)
table(Ap.low.geneList)

# create topGOdata object
Ap.low.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis",
                 ontology = "BP",
                 allGenes = Ap.low.geneList,
                 geneSel = selectTranscript,
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GO)

Ap.low.BP.GOdata

# perform enrichment analysis using parentchild method
Ap.low.BP.resultParentChild <- runTest(Ap.low.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ap.low.BP.resultParentChild
```

Enriched gene sets for genes expressed at low temperatures only.

```{r, echo=FALSE, results="asis"}
Ap.low.BP.ResTable <- GenTable(Ap.low.BP.GOdata, parentchild = Ap.low.BP.resultParentChild, topNodes=26)
Ap.low.BP.ResTable
write.table(Ap.low.BP.ResTable, file = paste(resultsdir, "Ap_low_GO.BP_results.txt", sep=""), quote=FALSE, row.names=FALSE, sep = "\t")
pandoc.table(Ap.low.BP.ResTable)
```

```{r, echo=FALSE}
# graph significant nodes
pdf(paste(resultsdir, "Ap_bimodal_BP_topGO_sig_nodes.pdf", sep=""))
showSigOfNodes(Ap.bim.BP.GOdata, score(Ap.bim.BP.resultParentChild), firstSigNodes = 52, useInfo = 'all')
dev.off()

pdf(paste(resultsdir, "Ap_intermediate_BP_topGO_sig_nodes.pdf", sep=""))
showSigOfNodes(Ap.int.BP.GOdata, score(Ap.int.BP.resultParentChild), firstSigNodes = 16, useInfo = 'all')
dev.off()

pdf(paste(resultsdir, "Ap_high_BP_topGO_sig_nodes.pdf", sep=""))
showSigOfNodes(Ap.hig.BP.GOdata, score(Ap.hig.BP.resultParentChild), firstSigNodes = 12, useInfo = 'all')
dev.off()

pdf(paste(resultsdir, "Ap_low_BP_topGO_sig_nodes.pdf", sep=""))
showSigOfNodes(Ap.low.BP.GOdata, score(Ap.low.BP.resultParentChild), firstSigNodes = 26, useInfo = 'all')
dev.off()
```

Visualize differens in GO terms among expression types using a wordcloud.

```{gsea_wordcloud, cache=TRUE}
# create Corpus for text mining
GO.dfs <- c(paste(Ap.low.BP.ResTable$Term, collapse=" "), paste(Ap.hig.BP.ResTable$Term, collapse=" "),
            paste(Ap.int.BP.ResTable$Term, collapse=" "), paste(Ap.bim.BP.ResTable$Term, collapse=" "))

vs <- VectorSource(GO.dfs)
GO.corp <- Corpus(vs)
GO.corp

GO.corp <- tm_map(GO.corp, removePunctuation)

# create matrix of terms
term.matrix <- TermDocumentMatrix(GO.corp)
term.matrix <- as.matrix(term.matrix)
colnames(term.matrix) <- c('low', 'high', 'intermediate', 'bimodal')

# comparison cloud
comparison.cloud(term.matrix, max.words=200, random.order=FALSE, scale=c(0.25, 0.5))
```


## Comparison of expression patterns between colonies

A set of `r nrow(interaction.transcripts)` transcripts have expression patterns that depend on colony. In this section, I examine how these transcripts differ between colonies. Specifically, I expect genes that are up-regulated at high temperatures in the northern *A22* colony to be at constitutively higher expression in the more southern *Ar* colony. Conversely, I will examine the extent to which genes upregulated at colder temperatures in the *Ar* colony are constitutively expressed in the *A22* colony.








## Visualize responsive transcripts

Make plots for all significant transcripts

```{r plot_responsive, echo=FALSE, eval=FALSE, cache=TRUE}
# Line plot, expression against temp, faceted by colony
p1 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_line() +
  facet_grid(. ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p1


# Smooth plot, expression against temp, faceted by colony
p2 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_smooth() +
  facet_grid(. ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p2

# Same as p1, faceted by expression type
p3 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_line() +
  facet_grid(exp_type ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p3

# Smooth plot, expression against temp, faceted by colony
p4 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_smooth() +
  facet_grid(exp_type ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p4
```

Make plots for transcripts with significant temperature by colony interaction

```{r plot_interaction_responsive, echo=FALSE, eval=FALSE, cache=TRUE}
# Line plot, expression against temp, faceted by colony
p1 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_line() +
  facet_grid(. ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p1


# Smooth plot, expression against temp, faceted by colony
p2 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_smooth() +
  facet_grid(. ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p2

# Same as p1, faceted by expression type
p3 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_line() +
  facet_grid(exp_type ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p3

# Smooth plot, expression against temp, faceted by colony
p4 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
  geom_smooth() +
  facet_grid(exp_type ~ colony) +
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
p4
```



```{r plot, echo=FALSE, eval=FALSE, cache=TRUE}
# Line plot, expression against temp, faceted by colony
png(paste(resultsdir, "Ap_expression_by_colony_line.png", sep=""))
  p1 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_line() +
    facet_grid(. ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p1
dev.off()

# Smooth plot, expression against temp, faceted by colony
png(paste(resultsdir, "Ap_expression_by_colony_smooth.png", sep=""))
  p2 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_smooth() +
    facet_grid(. ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p2
dev.off()

# Same as p1, faceted by expression type
png(paste(resultsdir, "Ap_expression_by_colony_exp_line.png", sep=""))
  p3 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_line() +
    facet_grid(exp_type ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p3
dev.off()

# Smooth plot, expression against temp, faceted by colony
png(paste(resultsdir, "Ap_expression_by_colony_exp_smooth.png", sep=""))
  p4 <- ggplot(Ap.dt, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_smooth() +
    facet_grid(exp_type ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p4
dev.off()

```


```{r plot_interaction_responsive_to_file, echo=FALSE, eval=FALSE, cache=TRUE}
# Line plot, expression against temp, faceted by colony
png(paste(resultsdir, "Ap_expression_interaction_by_colony_line.png", sep=""))
  p1 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_line() +
    facet_grid(. ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p1
dev.off()

# Smooth plot, expression against temp, faceted by colony
png(paste(resultsdir, "Ap_expression_interaction_by_colony_smooth.png", sep=""))
  p2 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_smooth() +
    facet_grid(. ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p2
dev.off()

# Same as p1, faceted by expression type
png(paste(resultsdir, "Ap_expression_interaction_by_colony_exp_line.png", sep=""))
  p3 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_line() +
    facet_grid(exp_type ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p3
dev.off()

# Smooth plot, expression against temp, faceted by colony
png(paste(resultsdir, "Ap_expression_interaction_by_colony_exp_smooth.png", sep=""))
  p4 <- ggplot(Ap.dt.interaction, aes(x=val, y=exp.scaled, group=Transcript)) +
    geom_smooth() +
    facet_grid(exp_type ~ colony) +
    scale_y_continuous(name="Expression (scaled)") +
    scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
  p4
dev.off()
```


## Session information ##

```{r session}
save.image()
sessionInfo()
```

## References

```{r references, results='asis'}
bibliography()
```
