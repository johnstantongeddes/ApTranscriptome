Thermal reactionome of the common ant species *Aphaenogaster*
================================================================
  
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**February 28, 2014**

**Technical Report No. 3**

**Department of Biology**

**University of Vermont**





## Summary ##
  
In this technical report, which accompanies the manuscript **Thermal reactionome of a common ant species** (Stanton-Geddes et al., in press), we:

1. describe the *de novo* assembly of the transcriptome for two ant colonies with in the *Aphaenogaster rudis-picea-fulva* species complex (<a href="http://dx.doi.org/10.1155/2012/752815">Lubertazzi, 2012</a>)
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

The Illumina reads were filtered using the program [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (<a href="http://dx.doi.org/10.1093/nar/gks540">Lohse et al. 2012</a>) to remove Ilumina adapter sequences and filter out bases with quality scores less than ??. 

~~~
TRIMMOMATIC CODE
~~~

This filtering yielded...

These reads were combined and used in *de novo* transcriptome assembly using the program [Trinity](http://trinityrnaseq.sourceforge.net/) (<a href="http://dx.doi.org/10.1038/nbt.1883">Grabherr et al. 2011</a>). Note that this required ??? GB RAM and ??? hours run-time and was run on the [Vermont Genetics Network](http://vgn.uvm.edu/) computing cluster. 

~~~
TRINITY CODE
~~~

This assembly contained 100,381 unique components (roughly genes) in 126,172 total transcripts (Table 1). 

As we were assembling two divergent colonies into a single transcriptome, we suspected that this assembly would be susceptible to known problems of errors during assembly (e.g. chimeric transcripts that are fusions of two transcripts) and redundancy (<a href="http://dx.doi.org/10.1186/1471-2164-14-328">Yang & Smith, 2013</a>). To account for this, we performed two post-assembly processing steps.

First, we ran the program [cap3](http://seq.cs.iastate.edu/) (<a href="http://dx.doi.org/10.1101/gr.9.9.868">Huang, 1999</a>) setting the maximum gap length and band expansion size to 50 `-f 50 -a 50`, no end clipping as the reads were already filtered `k 0`, requiring 90% identity for assembly, and a minimum overlap length of 100 bp `-o 100`. The percent identity threshold of 90% was chosen to liberally collapse orthologous contigs from the two colonies that may have been assembled separately. 

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

These post-processing step removed 16% of the initial reads (Table 1).


|    &nbsp;     |  Total contigs  |  Total length  |  Median contig size  |
|:-------------:|:---------------:|:--------------:|:--------------------:|
|  **trinity**  |     126,172     |  100,389,539   |         358          |
|  **reduced**  |     105,536     |   62,648,997   |         320          |

Table: Table 1: Statistics for Trinity and cap3+uclust reduced transcriptome assemblies (continued below)

 

|    &nbsp;     |  Mean contig size  |  N50 contig  |  N50 Length  |
|:-------------:|:------------------:|:------------:|:------------:|
|  **trinity**  |        795         |    16,201    |    1,631     |
|  **reduced**  |        593         |    15,491    |     895      |


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

Annotation was performed by uploading the reduced assembly "Trinity_cap3_uclust.fasta" to the web-based annotation program [FastAnnotator](http://fastannotator.cgu.edu.tw/index.php) (<a href="">unknown, unknown</a>).

Results are available as job ID [13894410176993](http://fastannotator.cgu.edu.tw/job.php?jobid=13894410176993#page=basicinfo).

This annotation file can be read directly to R:


```r
### Annotation file from either AWS or GoogleDrive
annotationURL <- getURL("http://johnstantongeddes.org/assets/files/Aphaeno_transcriptome_AnnotationTable.txt")
# a2 <- getURL('https://googledrive.com/host/0B75IymziRJ_9Tlg1U1Vxbjk1bzg')
# # GoogleDrive link

annotationfile <- read.csv(textConnection(annotationURL), header = TRUE, sep = "\t", 
    stringsAsFactors = FALSE)
str(annotationfile)

# Convert to data.table
annotationtable <- data.table(annotationfile)
head(annotationtable)
```


Transcriptome annotations are nearly impossible to visualize in a meaningful way. For lack of better ideas, I created a word cloud.


```r
library(tm)
library(wordcloud)
library(RColorBrewer)
library(stringr)

# subset and extract unique best hits to NCBI nr database
ann <- unique(annotationfile$best.hit.to.nr)
annsplit <- str_split_fixed(ann, pattern = " ", n = 2)
annvec <- annsplit[, 2]

# create Corpus of terms. convert to dataframe for plotting wordcloud
ap.corpus <- Corpus(VectorSource(annvec))
ap.corpus <- tm_map(ap.corpus, removePunctuation)
ap.corpus <- tm_map(ap.corpus, tolower)
ap.corpus <- tm_map(ap.corpus, function(x) removeWords(x, stopwords("english")))
tdm <- TermDocumentMatrix(ap.corpus)
m <- as.matrix(tdm)
```

```
## Error: cannot allocate vector of length 1013898960
```

```r
v <- sort(rowSums(m), decreasing = TRUE)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'sort': Error in is.data.frame(x) : object 'm' not found
## Calls: rowSums -> is.data.frame
```

```r
d <- data.frame(word = names(v), freq = v)
```

```
## Error: object 'v' not found
```

```r
pal <- brewer.pal(4, "Dark2")
pal <- pal[-(1:2)]

# png('wordcloud.png', width=1280,height=800)
wordcloud(d$word, d$freq, scale = c(8, 0.3), min.freq = 2, max.words = 100, 
    random.order = T, rot.per = 0.15, colors = pal, vfont = c("sans serif", 
        "plain"))
```

```
## Error: object 'd' not found
```

```r
# dev.off()
```


This wordcloud shows that we do not know much about the transcriptome - mostly hypothetical or predicted proteins. I removed these uninformative words and generated a new wordcloud. 


```r
# remove 'protein', 'predicted' and 'hypothetical'
rmwords <- c("hypothetical", "protein", "isoform", "subunit", "family", "putative", 
    "partial", "domain", "containing", "predicted", "domaincontaining", "uncharacterized")
ap.corpus.rmwords <- tm_map(ap.corpus, removeWords, rmwords)
```

```
## Error: object 'ap.corpus' not found
```

```r
tdm <- TermDocumentMatrix(ap.corpus.rmwords)
```

```
## Error: object 'ap.corpus.rmwords' not found
```

```r
m <- as.matrix(tdm)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'as.matrix': Error: object 'tdm' not found
```

```r
v <- sort(rowSums(m), decreasing = TRUE)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'sort': Error in is.data.frame(x) : object 'm' not found
## Calls: rowSums -> is.data.frame
```

```r
d <- data.frame(word = names(v), freq = v)
```

```
## Error: object 'v' not found
```

```r
pal <- brewer.pal(8, "Dark2")
pal <- pal[-(1:2)]
# png('wordcloud2.png', width=1280,height=800)
wordcloud(d$word, d$freq, scale = c(8, 0.3), min.freq = 2, max.words = 100, 
    random.order = T, rot.per = 0.15, colors = pal, vfont = c("sans serif", 
        "plain"))
```

```
## Error: object 'd' not found
```

```r
# dev.off()
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
                                                 

```r
# directory containing trimmed reads
readdir <- "data/ind_files/"
# sailfish index directory
sfindex <- "results/trinity-full/sailfish-index"
# sailfish expression directory
sfexpressionroot <- "results/trinity-full/sailfish-expression/"

# list of reads in each of four trimmed classes
readlist <- list.files(readdir)
(paired.left <- readlist[grep(".\\.paired.left.fastq$", readlist)])
```

```
##  [1] "A22-0_ATCACG.paired.left.fastq"  "A22-10_TGACCA.paired.left.fastq"
##  [3] "A22-14_ACAGTG.paired.left.fastq" "A22-17_GCCAAT.paired.left.fastq"
##  [5] "A22-21_CAGATC.paired.left.fastq" "A22-24_ACTTGA.paired.left.fastq"
##  [7] "A22-28_GATCAG.paired.left.fastq" "A22-31_TAGCTT.paired.left.fastq"
##  [9] "A22-35_GGCTAC.paired.left.fastq" "A22-38_CTTGTA.paired.left.fastq"
## [11] "A22-3_CGATGT.paired.left.fastq"  "A22-7_TTAGGC.paired.left.fastq" 
## [13] "Ar-0_AGTCAA.paired.left.fastq"   "Ar-10_CCGTCC.paired.left.fastq" 
## [15] "Ar-14_GTCCGC.paired.left.fastq"  "Ar-17_GTGAAA.paired.left.fastq" 
## [17] "Ar-21_GTGGCC.paired.left.fastq"  "Ar-24_GTTTCG.paired.left.fastq" 
## [19] "Ar-28_CGTACG.paired.left.fastq"  "Ar-31_GAGTGG.paired.left.fastq" 
## [21] "Ar-35_ACTGAT.paired.left.fastq"  "Ar-38_ATTCCT.paired.left.fastq" 
## [23] "Ar-3_AGTTCC.paired.left.fastq"   "Ar-7_ATGTCA.paired.left.fastq"
```

```r
(paired.right <- readlist[grep("\\.paired.right.fastq$", readlist)])
```

```
##  [1] "A22-0_ATCACG.paired.right.fastq"  "A22-10_TGACCA.paired.right.fastq"
##  [3] "A22-14_ACAGTG.paired.right.fastq" "A22-17_GCCAAT.paired.right.fastq"
##  [5] "A22-21_CAGATC.paired.right.fastq" "A22-24_ACTTGA.paired.right.fastq"
##  [7] "A22-28_GATCAG.paired.right.fastq" "A22-31_TAGCTT.paired.right.fastq"
##  [9] "A22-35_GGCTAC.paired.right.fastq" "A22-38_CTTGTA.paired.right.fastq"
## [11] "A22-3_CGATGT.paired.right.fastq"  "A22-7_TTAGGC.paired.right.fastq" 
## [13] "Ar-0_AGTCAA.paired.right.fastq"   "Ar-10_CCGTCC.paired.right.fastq" 
## [15] "Ar-14_GTCCGC.paired.right.fastq"  "Ar-17_GTGAAA.paired.right.fastq" 
## [17] "Ar-21_GTGGCC.paired.right.fastq"  "Ar-24_GTTTCG.paired.right.fastq" 
## [19] "Ar-28_CGTACG.paired.right.fastq"  "Ar-31_GAGTGG.paired.right.fastq" 
## [21] "Ar-35_ACTGAT.paired.right.fastq"  "Ar-38_ATTCCT.paired.right.fastq" 
## [23] "Ar-3_AGTTCC.paired.right.fastq"   "Ar-7_ATGTCA.paired.right.fastq"
```

```r
(unpaired.left <- readlist[grep("unpaired.left.fastq$", readlist)])
```

```
##  [1] "A22-0_ATCACG.unpaired.left.fastq" 
##  [2] "A22-10_TGACCA.unpaired.left.fastq"
##  [3] "A22-14_ACAGTG.unpaired.left.fastq"
##  [4] "A22-17_GCCAAT.unpaired.left.fastq"
##  [5] "A22-21_CAGATC.unpaired.left.fastq"
##  [6] "A22-24_ACTTGA.unpaired.left.fastq"
##  [7] "A22-28_GATCAG.unpaired.left.fastq"
##  [8] "A22-31_TAGCTT.unpaired.left.fastq"
##  [9] "A22-35_GGCTAC.unpaired.left.fastq"
## [10] "A22-38_CTTGTA.unpaired.left.fastq"
## [11] "A22-3_CGATGT.unpaired.left.fastq" 
## [12] "A22-7_TTAGGC.unpaired.left.fastq" 
## [13] "Ar-0_AGTCAA.unpaired.left.fastq"  
## [14] "Ar-10_CCGTCC.unpaired.left.fastq" 
## [15] "Ar-14_GTCCGC.unpaired.left.fastq" 
## [16] "Ar-17_GTGAAA.unpaired.left.fastq" 
## [17] "Ar-21_GTGGCC.unpaired.left.fastq" 
## [18] "Ar-24_GTTTCG.unpaired.left.fastq" 
## [19] "Ar-28_CGTACG.unpaired.left.fastq" 
## [20] "Ar-31_GAGTGG.unpaired.left.fastq" 
## [21] "Ar-35_ACTGAT.unpaired.left.fastq" 
## [22] "Ar-38_ATTCCT.unpaired.left.fastq" 
## [23] "Ar-3_AGTTCC.unpaired.left.fastq"  
## [24] "Ar-7_ATGTCA.unpaired.left.fastq"
```

```r
(unpaired.right <- readlist[grep("unpaired.right.fastq$", readlist)])
```

```
##  [1] "A22-0_ATCACG.unpaired.right.fastq" 
##  [2] "A22-10_TGACCA.unpaired.right.fastq"
##  [3] "A22-14_ACAGTG.unpaired.right.fastq"
##  [4] "A22-17_GCCAAT.unpaired.right.fastq"
##  [5] "A22-21_CAGATC.unpaired.right.fastq"
##  [6] "A22-24_ACTTGA.unpaired.right.fastq"
##  [7] "A22-28_GATCAG.unpaired.right.fastq"
##  [8] "A22-31_TAGCTT.unpaired.right.fastq"
##  [9] "A22-35_GGCTAC.unpaired.right.fastq"
## [10] "A22-38_CTTGTA.unpaired.right.fastq"
## [11] "A22-3_CGATGT.unpaired.right.fastq" 
## [12] "A22-7_TTAGGC.unpaired.right.fastq" 
## [13] "Ar-0_AGTCAA.unpaired.right.fastq"  
## [14] "Ar-10_CCGTCC.unpaired.right.fastq" 
## [15] "Ar-14_GTCCGC.unpaired.right.fastq" 
## [16] "Ar-17_GTGAAA.unpaired.right.fastq" 
## [17] "Ar-21_GTGGCC.unpaired.right.fastq" 
## [18] "Ar-24_GTTTCG.unpaired.right.fastq" 
## [19] "Ar-28_CGTACG.unpaired.right.fastq" 
## [20] "Ar-31_GAGTGG.unpaired.right.fastq" 
## [21] "Ar-35_ACTGAT.unpaired.right.fastq" 
## [22] "Ar-38_ATTCCT.unpaired.right.fastq" 
## [23] "Ar-3_AGTTCC.unpaired.right.fastq"  
## [24] "Ar-7_ATGTCA.unpaired.right.fastq"
```

```r

# Loop across each sample and quantify expression

# NOTE - samples listed in same order as given by the above lists
samples <- c("A22-0", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", 
    "A22-31", "A22-35", "A22-38", "A22-3", "A22-7", "Ar-0", "Ar-10", "Ar-14", 
    "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38", "Ar-3", "Ar-7")

for (j in 1:length(samples)) {
    message("Start expression quantification for sample ", samples[j], ": ", 
        Sys.time())
    quantdir <- paste(sfexpressionroot, samples[j], "_quant", sep = "")
    samp.pos <- grep(paste(paste(samples[j], "_", sep = "")), paired.left)
    samp.paired.l <- paste(readdir, paired.left[samp.pos], sep = "")
    samp.paired.r <- paste(readdir, paired.right[samp.pos], sep = "")
    samp.unpaired.l <- paste(readdir, unpaired.left[samp.pos], sep = "")
    samp.unpaired.r <- paste(readdir, unpaired.right[samp.pos], sep = "")
    sailfishcmd <- paste("sailfish quant -i ", sfindex, " -o ", quantdir, " --reads ", 
        samp.paired.l, " ", samp.paired.r, " ", samp.unpaired.l, " ", samp.unpaired.r, 
        " -p 4", sep = "")
    # print(sailfishcmd)
    system(sailfishcmd)
    message("Done with expression quantification for sample ", samples[j], ": ", 
        Sys.time(), "\n")
}
```

```
## Start expression quantification for sample A22-0: 2014-02-17 17:09:30
## Done with expression quantification for sample A22-0: 2014-02-17 17:12:42
## 
## Start expression quantification for sample A22-10: 2014-02-17 17:12:42
## Done with expression quantification for sample A22-10: 2014-02-17 17:15:54
## 
## Start expression quantification for sample A22-14: 2014-02-17 17:15:54
## Done with expression quantification for sample A22-14: 2014-02-17 17:19:06
## 
## Start expression quantification for sample A22-17: 2014-02-17 17:19:06
## Done with expression quantification for sample A22-17: 2014-02-17 17:22:18
## 
## Start expression quantification for sample A22-21: 2014-02-17 17:22:18
## Done with expression quantification for sample A22-21: 2014-02-17 17:25:29
## 
## Start expression quantification for sample A22-24: 2014-02-17 17:25:29
## Done with expression quantification for sample A22-24: 2014-02-17 17:28:41
## 
## Start expression quantification for sample A22-28: 2014-02-17 17:28:41
## Done with expression quantification for sample A22-28: 2014-02-17 17:31:52
## 
## Start expression quantification for sample A22-31: 2014-02-17 17:31:52
## Done with expression quantification for sample A22-31: 2014-02-17 17:35:04
## 
## Start expression quantification for sample A22-35: 2014-02-17 17:35:04
## Done with expression quantification for sample A22-35: 2014-02-17 17:38:15
## 
## Start expression quantification for sample A22-38: 2014-02-17 17:38:15
## Done with expression quantification for sample A22-38: 2014-02-17 17:41:28
## 
## Start expression quantification for sample A22-3: 2014-02-17 17:41:28
## Done with expression quantification for sample A22-3: 2014-02-17 17:44:41
## 
## Start expression quantification for sample A22-7: 2014-02-17 17:44:41
## Done with expression quantification for sample A22-7: 2014-02-17 17:47:53
## 
## Start expression quantification for sample Ar-0: 2014-02-17 17:47:53
## Done with expression quantification for sample Ar-0: 2014-02-17 17:51:05
## 
## Start expression quantification for sample Ar-10: 2014-02-17 17:51:05
## Done with expression quantification for sample Ar-10: 2014-02-17 17:54:16
## 
## Start expression quantification for sample Ar-14: 2014-02-17 17:54:16
## Done with expression quantification for sample Ar-14: 2014-02-17 17:57:28
## 
## Start expression quantification for sample Ar-17: 2014-02-17 17:57:28
## Done with expression quantification for sample Ar-17: 2014-02-17 18:00:40
## 
## Start expression quantification for sample Ar-21: 2014-02-17 18:00:40
## Done with expression quantification for sample Ar-21: 2014-02-17 18:03:52
## 
## Start expression quantification for sample Ar-24: 2014-02-17 18:03:52
## Done with expression quantification for sample Ar-24: 2014-02-17 18:07:04
## 
## Start expression quantification for sample Ar-28: 2014-02-17 18:07:04
## Done with expression quantification for sample Ar-28: 2014-02-17 18:10:15
## 
## Start expression quantification for sample Ar-31: 2014-02-17 18:10:15
## Done with expression quantification for sample Ar-31: 2014-02-17 18:13:26
## 
## Start expression quantification for sample Ar-35: 2014-02-17 18:13:26
## Done with expression quantification for sample Ar-35: 2014-02-17 18:16:38
## 
## Start expression quantification for sample Ar-38: 2014-02-17 18:16:38
## Done with expression quantification for sample Ar-38: 2014-02-17 18:19:50
## 
## Start expression quantification for sample Ar-3: 2014-02-17 18:19:50
## Done with expression quantification for sample Ar-3: 2014-02-17 18:23:02
## 
## Start expression quantification for sample Ar-7: 2014-02-17 18:23:02
## Done with expression quantification for sample Ar-7: 2014-02-17 18:26:13
```


This generated a directory for each sample



and within each directory there are the following r:



The file *quant_bias_corrected.sf* contains the following columns, following a number of header lines:

1. Transcript ID
2. Transcript Length
3. Transcripts per Million (TPM): computed as described in (<a href="http://dx.doi.org/10.1093/bioinformatics/btp692">Li et al. 2009</a>), and is meant as an estimate of the number of transcripts, per million observed transcripts, originating from each isoform.
4. Reads Per Kilobase per Million mapped reads (RPKM): classic measure of relative transcript abundance, and is an estimate of the number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

The TPM column for each sample was extracted and combined into a matrix for each colony.


```
## Warning: cannot open file
## 'results/trinity-full/sailfish-expression/A22-0_quant/quant_bias_corrected.sf':
## No such file or directory
```

```
## Error: cannot open the connection
```

```
## Error: object 'A22_0_quant' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'A22.TPM' not found
```

```
## Error: object 'A22.TPM' not found
```

```
## Error: object 'Ar_0_quant' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'Ar.TPM' not found
```

```
## Error: object 'Ar.TPM' not found
```


Note that expression levels at each temperature treatment are highly correlated between the two colonies.


```
## Error: object 'A22_0_quant' not found
```

```
## Error: object 'cors' not found
```

```
## Error: object 'cortable' not found
```


## Identification of thermally-responsive genes

To identify transcripts (roughly equivalent to genes) that show thermal responsiveness, I fit the following linear model to each transcript:

$$ TPM = \beta_0 + \beta_1(colony) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(colony * temp) + \\beta_5(colony * temp^2) + \epsilon $$

where TPM is transcripts per million. 

For this list of P-values, correct for multiple testing using False Discovery Rate (FDR).

Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the analysis to be conservative). 


```r
A22.TPM[, `:=`(colony, "A22")]
```

```
##                      Transcript Length     TPM    RPKM    KPKM
##       1:         0|*|Contig6267   9990 0.08478 0.10175 0.10175
##       2:         0|*|Contig6267   9990 0.03762 0.05613 0.05613
##       3:         0|*|Contig6267   9990 0.07277 0.09098 0.09098
##       4:         0|*|Contig6267   9990 0.04037 0.05347 0.05347
##       5:         0|*|Contig6267   9990 0.02216 0.03387 0.03387
##      ---                                                      
## 1266428: 9|*|comp147140_c0_seq1   9030 0.77718 1.22699 1.22699
## 1266429: 9|*|comp147140_c0_seq1   9030 0.75068 1.19098 1.19098
## 1266430: 9|*|comp147140_c0_seq1   9030 1.03635 1.54603 1.54603
## 1266431: 9|*|comp147140_c0_seq1   9030 0.56836 0.84550 0.84550
## 1266432: 9|*|comp147140_c0_seq1   9030 0.43509 0.64956 0.64956
##          EstimatedNumReads sample  val colony
##       1:            0.5124  A22-0  0.0    A22
##       2:            0.2326  A22-3  3.5    A22
##       3:            0.4274  A22-7  7.0    A22
##       4:            0.2918 A22-10 10.5    A22
##       5:            0.1676 A22-14 14.0    A22
##      ---                                     
## 1266428:            4.9758 A22-24 24.5    A22
## 1266429:            3.6599 A22-28 28.0    A22
## 1266430:            5.3611 A22-31 31.5    A22
## 1266431:            3.1595 A22-35 35.0    A22
## 1266432:            3.0201 A22-38 38.5    A22
```

```r
Ar.TPM[, `:=`(colony, "Ar")]
```

```
##                      Transcript Length     TPM    RPKM    KPKM
##       1:         0|*|Contig6267   9990 0.04547 0.09485 0.09485
##       2:         0|*|Contig6267   9990 0.09140 0.15360 0.15360
##       3:         0|*|Contig6267   9990 0.13774 0.16723 0.16723
##       4:         0|*|Contig6267   9990 0.16978 0.31125 0.31125
##       5:         0|*|Contig6267   9990 0.14226 0.27254 0.27254
##      ---                                                      
## 1266428: 9|*|comp147140_c0_seq1   9030 0.51914 0.94175 0.94175
## 1266429: 9|*|comp147140_c0_seq1   9030 0.41781 0.82656 0.82656
## 1266430: 9|*|comp147140_c0_seq1   9030 0.55720 1.11393 1.11393
## 1266431: 9|*|comp147140_c0_seq1   9030 0.29183 0.64831 0.64831
## 1266432: 9|*|comp147140_c0_seq1   9030 0.31166 0.64685 0.64685
##          EstimatedNumReads sample  val colony
##       1:            0.1878   Ar-0  0.0     Ar
##       2:            0.3936   Ar-3  3.5     Ar
##       3:            0.6490   Ar-7  7.0     Ar
##       4:            0.7581  Ar-10 10.5     Ar
##       5:            0.8265  Ar-14 14.0     Ar
##      ---                                     
## 1266428:            2.2410  Ar-24 24.5     Ar
## 1266429:            1.7599  Ar-28 28.0     Ar
## 1266430:            2.3732  Ar-31 31.5     Ar
## 1266431:            1.0365  Ar-35 35.0     Ar
## 1266432:            1.4125  Ar-38 38.5     Ar
```

```r
TPM.dt <- rbind(A22.TPM, Ar.TPM)
TPM.dt$colony <- as.factor(TPM.dt$colony)
str(TPM.dt)
```

```
## Classes 'data.table' and 'data.frame':	2532864 obs. of  9 variables:
##  $ Transcript       : chr  "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" ...
##  $ Length           : int  9990 9990 9990 9990 9990 9990 9990 9990 9990 9990 ...
##  $ TPM              : num  0.0848 0.0376 0.0728 0.0404 0.0222 ...
##  $ RPKM             : num  0.1017 0.0561 0.091 0.0535 0.0339 ...
##  $ KPKM             : num  0.1017 0.0561 0.091 0.0535 0.0339 ...
##  $ EstimatedNumReads: num  0.512 0.233 0.427 0.292 0.168 ...
##  $ sample           : chr  "A22-0" "A22-3" "A22-7" "A22-10" ...
##  $ val              : num  0 3.5 7 10.5 14 17.5 21 24.5 28 31.5 ...
##  $ colony           : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r

setkey(TPM.dt, val)
TPM.dt.sub <- TPM.dt[val != 7]
unique(TPM.dt.sub$val)
```

```
##  [1]  0.0  3.5 10.5 14.0 17.5 21.0 24.5 28.0 31.5 35.0 38.5
```

```r

# define model for RxN function
model <- "TPM ~ colony + val + I(val^2) + colony:val + colony:I(val^2)"

# identify responsive transcripts
RxNout <- RxNseq(f = TPM.dt.sub, model = model)

save.image("RxN_combined_results.RData")
```

```
## Warning: 'package:R.utils' may not be available when loading
```


Of the 105536 transcripts, 22627 have models with P < 0.05.

Many of these are likely false positives, so I adjusted P-values using false discovery rate (FDR) to identify only those transcripts with less than 5% FDR as significant. 


```r
RxNout$padj <- p.adjust(RxNout$pval, method = "fdr")
# Plot FDR values against initial pvalues
par(mfrow = c(2, 1))
hist(RxNout$pval)
hist(RxNout$padj)
```

![plot of chunk fdr](figure/fdr.png) 

```r

# subset to significant transcripts
signif.transcripts <- RxNout[which(RxNout$padj < 0.05), ]
```


## Functional annotation ##

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  



```r
# add annotation information
setkey(annotationtable, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
```



|   Coefficient    |  Number_significant  |
|:----------------:|:--------------------:|
|      Total       |         8817         |
|      Colony      |         7937         |
|     Temp.lin     |         4213         |
|    Temp.quad     |         2477         |
| Colony:Temp.lin  |         2644         |
| Colony:Temp.quad |         2383         |

Table: Number of transcripts for which each term is significant


Of these, subset to those that have significant responses to temperature, either through a direct effect or interaction with colony. Add annotation information and write results to file. Do the same for transcripts that differ in expression between the colonies.


```r
responsive.transcripts <- signif.transcripts[!is.na(signif.transcripts$coef.val) | 
    !is.na(signif.transcripts$"coef.I(val^2)") | !is.na(signif.transcripts$"coef.colony:val") | 
    !is.na(signif.transcripts$"coef.colony:I(val^2)")]
dim(responsive.transcripts)
```

```
## [1] 5580    8
```

```r
# join signif transcripts with annotation
responsive.transcripts <- annotationtable[responsive.transcripts]
str(responsive.transcripts)
```

```
## Classes 'data.table' and 'data.frame':	5580 obs. of  19 variables:
##  $ Sequence.Name        : chr  "100015|*|comp3543055_c0_seq1" "100067|*|comp3557646_c0_seq1" "100089|*|comp11313_c1_seq1" "10016|*|comp130697_c0_seq1" ...
##  $ sequence.length      : int  208 208 208 1320 207 207 207 207 207 207 ...
##  $ best.hit.to.nr       : chr  "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|493136460|ref|WP_006154899.1| tyrosyl-tRNA synthetase " "gi|497544620|ref|WP_009858818.1| lipid-A-disaccharide synthase " "gi|490412587|ref|WP_004285230.1| diaminopimelate decarboxylase " ...
##  $ hit.length           : chr  "25" "68" "67" "299" ...
##  $ E.value              : chr  "7.22e-06" "1.5e-39" "3.37e-23" "1.13e-71" ...
##  $ Bit.score            : chr  "57.514374" "159.812994" "111.355753" "266.598396" ...
##  $ GO.Biological.Process: chr  "-" "-" "GO:0009245 lipid A biosynthetic process" "GO:0009089 lysine biosynthetic process via diaminopimelate" ...
##  $ GO.Cellular.Component: chr  "-" "-" "GO:0009276 Gram-negative-bacterium-type cell wall" "-" ...
##  $ GO.Molecular.Function: chr  "-" "-" "GO:0008915 lipid-A-disaccharide synthase activity" "GO:0030170 pyridoxal phosphate binding | GO:0008836 diaminopimelate decarboxylase activity" ...
##  $ Enzyme               : chr  "-" "-" "-" "-" ...
##  $ Domain               : chr  "-" "-" "-" "pfam02784 Orn_Arg_deC_N | pfam00278 Orn_DAP_Arg_deC" ...
##  $ annotation.type      : chr  "" "" "GO only" "GO & Domain" ...
##  $ pval                 : num  0.00152 0.002161 0.002113 0.000641 0.002516 ...
##  $ coef.colony          : num  0.045717 NA NA 0.000471 0.001497 ...
##  $ coef.val             : num  NA NA 0.018 0.00602 0.00457 ...
##  $ coef.I(val^2)        : num  0.00258 0.00244 NA NA NA ...
##  $ coef.colony:val      : num  NA NA 0.01101 0.00602 NA ...
##  $ coef.colony:I(val^2) : num  0.00261 0.00293 0.00382 NA NA ...
##  $ padj                 : num  0.0254 0.0322 0.0318 0.0138 0.0356 ...
##  - attr(*, "sorted")= chr "Sequence.Name"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r

# transcripts that differ in expression by colony
colony.transcripts <- signif.transcripts[!is.na(signif.transcripts$coef.colony)]
colony.transcripts <- annotationtable[colony.transcripts]
str(colony.transcripts)
```

```
## Classes 'data.table' and 'data.frame':	7937 obs. of  19 variables:
##  $ Sequence.Name        : chr  "100015|*|comp3543055_c0_seq1" "100129|*|comp116318_c0_seq1" "100148|*|comp125464_c0_seq1" "10015|*|comp139203_c0_seq1" ...
##  $ sequence.length      : int  208 207 207 1320 1320 207 207 207 207 207 ...
##  $ best.hit.to.nr       : chr  "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|307171928|gb|EFN63559.1| UPF0439 protein C9orf30-like protein " "gi|322798083|gb|EFZ19922.1| hypothetical protein SINV_07083 " "-" ...
##  $ hit.length           : chr  "25" "62" "69" "-" ...
##  $ E.value              : chr  "7.22e-06" "5.49e-12" "8.02e-20" "-" ...
##  $ Bit.score            : chr  "57.514374" "76.807535" "101.036156" "-" ...
##  $ GO.Biological.Process: chr  "-" "-" "GO:0006508 proteolysis" "-" ...
##  $ GO.Cellular.Component: chr  "-" "-" "-" "-" ...
##  $ GO.Molecular.Function: chr  "-" "-" "GO:0004252 serine-type endopeptidase activity" "-" ...
##  $ Enzyme               : chr  "-" "-" "-" "-" ...
##  $ Domain               : chr  "-" "pfam13873 Myb_DNA-bind_5" "-" "-" ...
##  $ annotation.type      : chr  "" "Domain only" "GO only" "" ...
##  $ pval                 : num  0.00152 0.000354 0.000993 0.000243 0.000641 ...
##  $ coef.colony          : num  4.57e-02 7.92e-06 7.43e-05 7.56e-06 4.71e-04 ...
##  $ coef.val             : num  NA NA NA NA 0.00602 ...
##  $ coef.I(val^2)        : num  0.00258 NA NA NA NA ...
##  $ coef.colony:val      : num  NA NA NA NA 0.00602 ...
##  $ coef.colony:I(val^2) : num  0.00261 NA NA NA NA ...
##  $ padj                 : num  0.02545 0.00921 0.01881 0.00701 0.01385 ...
##  - attr(*, "sorted")= chr "Sequence.Name"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
colony.transcripts <- colony.transcripts[order(colony.transcripts$padj), ]
write.table(colony.transcripts, file = paste(resultsdir, "Ap_colony_transcripts_GO.txt", 
    sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)


# transcripts that have colony by temperature interactions
interaction.transcripts <- signif.transcripts[!is.na(signif.transcripts$"coef.colony:val") | 
    !is.na(signif.transcripts$"coef.colony:I(val^2)")]
interaction.transcripts <- annotationtable[interaction.transcripts]
str(interaction.transcripts)
```

```
## Classes 'data.table' and 'data.frame':	3760 obs. of  19 variables:
##  $ Sequence.Name        : chr  "100015|*|comp3543055_c0_seq1" "100067|*|comp3557646_c0_seq1" "100089|*|comp11313_c1_seq1" "10016|*|comp130697_c0_seq1" ...
##  $ sequence.length      : int  208 208 208 1320 207 207 207 207 207 207 ...
##  $ best.hit.to.nr       : chr  "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|493136460|ref|WP_006154899.1| tyrosyl-tRNA synthetase " "gi|497544620|ref|WP_009858818.1| lipid-A-disaccharide synthase " "gi|490412587|ref|WP_004285230.1| diaminopimelate decarboxylase " ...
##  $ hit.length           : chr  "25" "68" "67" "299" ...
##  $ E.value              : chr  "7.22e-06" "1.5e-39" "3.37e-23" "1.13e-71" ...
##  $ Bit.score            : chr  "57.514374" "159.812994" "111.355753" "266.598396" ...
##  $ GO.Biological.Process: chr  "-" "-" "GO:0009245 lipid A biosynthetic process" "GO:0009089 lysine biosynthetic process via diaminopimelate" ...
##  $ GO.Cellular.Component: chr  "-" "-" "GO:0009276 Gram-negative-bacterium-type cell wall" "-" ...
##  $ GO.Molecular.Function: chr  "-" "-" "GO:0008915 lipid-A-disaccharide synthase activity" "GO:0030170 pyridoxal phosphate binding | GO:0008836 diaminopimelate decarboxylase activity" ...
##  $ Enzyme               : chr  "-" "-" "-" "-" ...
##  $ Domain               : chr  "-" "-" "-" "pfam02784 Orn_Arg_deC_N | pfam00278 Orn_DAP_Arg_deC" ...
##  $ annotation.type      : chr  "" "" "GO only" "GO & Domain" ...
##  $ pval                 : num  0.00152 0.002161 0.002113 0.000641 0.001524 ...
##  $ coef.colony          : num  0.045717 NA NA 0.000471 0.006795 ...
##  $ coef.val             : num  NA NA 0.018 0.00602 0.01652 ...
##  $ coef.I(val^2)        : num  0.00258 0.00244 NA NA 0.04473 ...
##  $ coef.colony:val      : num  NA NA 0.01101 0.00602 0.01652 ...
##  $ coef.colony:I(val^2) : num  0.00261 0.00293 0.00382 NA 0.04473 ...
##  $ padj                 : num  0.0254 0.0322 0.0318 0.0138 0.0255 ...
##  - attr(*, "sorted")= chr "Sequence.Name"
##  - attr(*, ".internal.selfref")=<externalptr>
```


For responsive transcripts, identify by shape of response:

* High - increase expression with temperature
* Low - decrease expression with temperature
* Intermediate - maximum expression at intermediate temperatures (14 - 28C)
* Bimodal - expressed greater than two standard deviations of expression at both low and high temperatures


```r
# merge RxN results with expression values
setkey(TPM.dt.sub, Transcript)
Ap.dt <- TPM.dt.sub[responsive.transcripts]
setkey(Ap.dt, Transcript)

Ap.exp_type <- ddply(Ap.dt, .(Transcript), function(df1) {
    lmout <- lm(TPM ~ val + I(val^2), data = df1)
    vals <- c(0, 3.5, 10, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)
    newdf <- data.frame(val = vals)
    pout <- predict(lmout, newdata = newdf)
    pout <- data.frame(val = vals, exp = pout)
    
    # get vals of max and min expression
    max.val = vals[which(pout$exp == max(pout$exp))]
    min.val = vals[which(pout$exp == min(pout$exp))]
    
    # report coefficients coef(lmout)
    exp_type = if (coef(lmout)["val"] > 0 & coef(lmout)["I(val^2)"] > 0) 
        "High" else {
        if (coef(lmout)["val"] < 0 & coef(lmout)["I(val^2)"] < 0) 
            "Low" else {
            if (coef(lmout)["val"] > 0 & coef(lmout)["I(val^2)"] < 0) 
                "Intermediate" else {
                "convex"
            }
        }
    }
    
    # for transcripts with convex exp_type, check if expression is truly bimodal
    if (exp_type == "convex") {
        if (max(pout[pout$val <= 10, "exp"]) > 2 * sd(pout$exp) & max(pout[pout$val >= 
            31.5, "exp"]) > 2 * sd(pout$exp)) 
            exp_type = "Bimodal" else {
            # linear increase?
            if (max.val > min.val) 
                exp_type = "High" else exp_type = "Low"
        }
    }
    
    # return values
    return(c(max.val = vals[which(pout$exp == max(pout$exp))], min.val = vals[which(pout$exp == 
        min(pout$exp))], exp_type = exp_type))
})

# merge 'exp_type' information with Ap.dt
Ap.exp_type <- data.table(Ap.exp_type)
setkey(Ap.exp_type, Transcript)
Ap.dt <- Ap.dt[Ap.exp_type]

# merge 'exp_type' with responsive.transcripts
setkey(responsive.transcripts, Sequence.Name)
responsive.transcripts <- Ap.exp_type[responsive.transcripts]
dim(responsive.transcripts)
```

```
## [1] 5580   22
```

```r
# order by padj
responsive.transcripts <- responsive.transcripts[order(responsive.transcripts$padj), 
    ]
# save results to file
write.table(responsive.transcripts, file = paste(resultsdir, "Ap_responsive_transcripts_GO.txt", 
    sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

# merge 'exp_type' with interaction transcripts
setkey(interaction.transcripts, Sequence.Name)
interaction.transcripts <- Ap.exp_type[interaction.transcripts]
dim(interaction.transcripts)
```

```
## [1] 3760   22
```

```r
# order by padj
interaction.transcripts <- interaction.transcripts[order(interaction.transcripts$padj), 
    ]
# save results to file
write.table(interaction.transcripts, file = paste(resultsdir, "Ap_interaction_transcripts_GO.txt", 
    sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
```



|  Bimodal  |  High  |  Intermediate  |  Low  |
|:---------:|:------:|:--------------:|:-----:|
|   1990    |  822   |      710       | 2058  |

Table: Number of transcripts with maximum expression at high, low, intermediate or both high and low (bimodal) temperatures.


Note that among responsive transcripts, there are 25 transcripts with GO term "response to stress" and various heat shock related proteins:


```r
unique(Ap.dt[grep("GO:0006950", Ap.dt$GO.Biological.Process), list(Transcript, 
    best.hit.to.nr)])
```

```
##                      Transcript
##  1:   1038|*|comp150483_c5_seq3
##  2:  11281|*|comp146961_c0_seq1
##  3:  12704|*|comp144775_c1_seq1
##  4:     14|*|comp150262_c0_seq1
##  5:           1504|*|Contig2729
##  6:  15115|*|comp132715_c0_seq1
##  7:  17710|*|comp150271_c3_seq3
##  8:          19475|*|Contig1438
##  9:   2087|*|comp150483_c5_seq1
## 10:  21384|*|comp149042_c0_seq3
## 11:  21598|*|comp142101_c0_seq1
## 12:  23441|*|comp114823_c1_seq1
## 13:   2604|*|comp148324_c0_seq4
## 14:   2832|*|comp150483_c5_seq2
## 15:  32312|*|comp933733_c0_seq1
## 16: 37752|*|comp1620595_c0_seq1
## 17:   3995|*|comp145243_c0_seq1
## 18: 47691|*|comp1460938_c0_seq1
## 19: 51985|*|comp1012776_c0_seq1
## 20:    552|*|comp147487_c0_seq1
## 21:  58246|*|comp109744_c0_seq1
## 22:    6075|*|comp92770_c0_seq1
## 23: 75624|*|comp2836178_c0_seq1
## 24:  80544|*|comp132706_c0_seq1
## 25:   9316|*|comp147545_c4_seq2
## 26:           9372|*|Contig4757
##                      Transcript
##                                                                                    best.hit.to.nr
##  1:                              gi|332022897|gb|EGI63169.1| Protein lethal(2)essential for life 
##  2:                                  gi|332029692|gb|EGI69571.1| G-protein coupled receptor Mth2 
##  3:         gi|380028536|ref|XP_003697954.1| PREDICTED: protein lethal(2)essential for life-like 
##  4: gi|332019420|gb|EGI59904.1| Putative fat-like cadherin-related tumor suppressor-like protein 
##  5:                                    gi|332023134|gb|EGI63390.1| Sugar transporter ERD6-like 6 
##  6:                                            gi|194716766|gb|ACF93232.1| heat shock protein 90 
##  7:                                  gi|322799248|gb|EFZ20646.1| hypothetical protein SINV_03807 
##  8:                                  gi|322799248|gb|EFZ20646.1| hypothetical protein SINV_03807 
##  9:                              gi|332022897|gb|EGI63169.1| Protein lethal(2)essential for life 
## 10:                         gi|396467618|ref|XP_003837992.1| hypothetical protein LEMA_P120390.1 
## 11:                              gi|332018201|gb|EGI58806.1| Protein lethal(2)essential for life 
## 12:                           gi|443696809|gb|ELT97425.1| hypothetical protein CAPTEDRAFT_194915 
## 13:                                   gi|307176228|gb|EFN65864.1| hypothetical protein EAG_10145 
## 14:                              gi|332022897|gb|EGI63169.1| Protein lethal(2)essential for life 
## 15:                          gi|367054010|ref|XP_003657383.1| hypothetical protein THITE_2156506 
## 16:                                       gi|295131654|ref|YP_003582317.1| ferritin-like protein 
## 17:                              gi|307211659|gb|EFN87680.1| Heat shock 70 kDa protein cognate 5 
## 18:                      gi|302922354|ref|XP_003053448.1| hypothetical protein NECHADRAFT_102357 
## 19:                                            gi|227018528|gb|ACP18866.1| heat shock protein 30 
## 20:                             gi|332026309|gb|EGI66443.1| RhoA activator C11orf59-like protein 
## 21:                                    gi|493322437|ref|WP_006279741.1| molecular chaperone DnaK 
## 22:             gi|332030037|gb|EGI69862.1| Serine/threonine-protein kinase PINK1, mitochondrial 
## 23:                                                    gi|50418863|ref|XP_457952.1| DEHA2C06072p 
## 24:                                            gi|121605727|ref|YP_983056.1| OsmC family protein 
## 25:                              gi|332020093|gb|EGI60539.1| Heat shock factor-binding protein 1 
## 26:                                  gi|332029691|gb|EGI69570.1| G-protein coupled receptor Mth2 
##                                                                                    best.hit.to.nr
```

```r
unique(Ap.dt[grep("heat shock", Ap.dt$best.hit.to.nr), ])
```

```
##                     Transcript Length    TPM   RPKM   KPKM
## 1:  15115|*|comp132715_c0_seq1    915 0.4329 0.5195 0.5195
## 2: 51985|*|comp1012776_c0_seq1    323 0.0000 0.0000 0.0000
##    EstimatedNumReads sample val colony sequence.length
## 1:            0.2351  A22-0   0    A22             915
## 2:            0.0000  A22-0   0    A22             323
##                                        best.hit.to.nr hit.length   E.value
## 1: gi|194716766|gb|ACF93232.1| heat shock protein 90         305 3.34e-186
## 2: gi|227018528|gb|ACP18866.1| heat shock protein 30          63  1.79e-25
##     Bit.score
## 1: 641.693336
## 2: 121.226673
##                                                                                                   GO.Biological.Process
## 1: GO:0006457 protein folding | GO:0043581 mycelium development | GO:0006950 response to stress | GO:0007049 cell cycle
## 2:                                                                                        GO:0006950 response to stress
##    GO.Cellular.Component
## 1:                     -
## 2:  GO:0005737 cytoplasm
##                                                                    GO.Molecular.Function
## 1: GO:0005525 GTP binding | GO:0005524 ATP binding | GO:0051082 unfolded protein binding
## 2:                                                                                     -
##    Enzyme          Domain annotation.type      pval coef.colony coef.val
## 1:      - pfam00183 HSP90     GO & Domain 4.643e-05   0.0002053       NA
## 2:      -               -         GO only 7.279e-05   0.0188945 0.002078
##    coef.I(val^2) coef.colony:val coef.colony:I(val^2)     padj max.val
## 1:     0.0002269              NA             0.001224 0.002203    38.5
## 2:     0.0030336        0.002078             0.003034 0.002984    38.5
##    min.val exp_type
## 1:    17.5  Bimodal
## 2:      14     High
```

```r
unique(Ap.dt[grep("Heat shock", Ap.dt$best.hit.to.nr), list(Transcript, best.hit.to.nr)])
```

```
##                   Transcript
## 1: 3995|*|comp145243_c0_seq1
## 2: 9316|*|comp147545_c4_seq2
##                                                      best.hit.to.nr
## 1: gi|307211659|gb|EFN87680.1| Heat shock 70 kDa protein cognate 5 
## 2: gi|332020093|gb|EGI60539.1| Heat shock factor-binding protein 1
```


Export data for interactive shiny app. 


```r
# scale expression values
Ap.dt[, `:=`(exp.scaled, scale(TPM)), by = Transcript]
```

```
##                           Transcript Length    TPM   RPKM   KPKM
##      1: 100015|*|comp3543055_c0_seq1    208 0.8024 0.9629 0.9629
##      2: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      3: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      4: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      5: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##     ---                                                         
## 122756:       9|*|comp147140_c0_seq1   9030 0.5572 1.1139 1.1139
## 122757:       9|*|comp147140_c0_seq1   9030 0.5684 0.8455 0.8455
## 122758:       9|*|comp147140_c0_seq1   9030 0.2918 0.6483 0.6483
## 122759:       9|*|comp147140_c0_seq1   9030 0.4351 0.6496 0.6496
## 122760:       9|*|comp147140_c0_seq1   9030 0.3117 0.6469 0.6469
##         EstimatedNumReads sample  val colony sequence.length
##      1:           0.09192  A22-0  0.0    A22             208
##      2:           0.00000   Ar-0  0.0     Ar             208
##      3:           0.00000  A22-3  3.5    A22             208
##      4:           0.00000   Ar-3  3.5     Ar             208
##      5:           0.00000 A22-10 10.5    A22             208
##     ---                                                     
## 122756:           2.37323  Ar-31 31.5     Ar            9030
## 122757:           3.15951 A22-35 35.0    A22            9030
## 122758:           1.03648  Ar-35 35.0     Ar            9030
## 122759:           3.02014 A22-38 38.5    A22            9030
## 122760:           1.41250  Ar-38 38.5     Ar            9030
##                                                         best.hit.to.nr
##      1: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      2: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      3: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      4: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      5: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##     ---                                                               
## 122756:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 122757:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 122758:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 122759:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 122760:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
##         hit.length  E.value   Bit.score
##      1:         25 7.22e-06   57.514374
##      2:         25 7.22e-06   57.514374
##      3:         25 7.22e-06   57.514374
##      4:         25 7.22e-06   57.514374
##      5:         25 7.22e-06   57.514374
##     ---                                
## 122756:       2101      0.0 4740.368306
## 122757:       2101      0.0 4740.368306
## 122758:       2101      0.0 4740.368306
## 122759:       2101      0.0 4740.368306
## 122760:       2101      0.0 4740.368306
##                                                                                                                                                                                                          GO.Biological.Process
##      1:                                                                                                                                                                                                                      -
##      2:                                                                                                                                                                                                                      -
##      3:                                                                                                                                                                                                                      -
##      4:                                                                                                                                                                                                                      -
##      5:                                                                                                                                                                                                                      -
##     ---                                                                                                                                                                                                                       
## 122756: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 122757: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 122758: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 122759: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 122760: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
##                                                     GO.Cellular.Component
##      1:                                                                 -
##      2:                                                                 -
##      3:                                                                 -
##      4:                                                                 -
##      5:                                                                 -
##     ---                                                                  
## 122756: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 122757: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 122758: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 122759: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 122760: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
##                                                                                                                         GO.Molecular.Function
##      1:                                                                                                                                     -
##      2:                                                                                                                                     -
##      3:                                                                                                                                     -
##      4:                                                                                                                                     -
##      5:                                                                                                                                     -
##     ---                                                                                                                                      
## 122756: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 122757: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 122758: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 122759: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 122760: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
##         Enzyme
##      1:      -
##      2:      -
##      3:      -
##      4:      -
##      5:      -
##     ---       
## 122756:      -
## 122757:      -
## 122758:      -
## 122759:      -
## 122760:      -
##                                                                                                                                                     Domain
##      1:                                                                                                                                                  -
##      2:                                                                                                                                                  -
##      3:                                                                                                                                                  -
##      4:                                                                                                                                                  -
##      5:                                                                                                                                                  -
##     ---                                                                                                                                                   
## 122756: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 122757: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 122758: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 122759: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 122760: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
##         annotation.type      pval coef.colony coef.val coef.I(val^2)
##      1:                 0.0015199   0.0457166       NA      0.002578
##      2:                 0.0015199   0.0457166       NA      0.002578
##      3:                 0.0015199   0.0457166       NA      0.002578
##      4:                 0.0015199   0.0457166       NA      0.002578
##      5:                 0.0015199   0.0457166       NA      0.002578
##     ---                                                             
## 122756:     GO & Domain 0.0006629   0.0001827  0.01015            NA
## 122757:     GO & Domain 0.0006629   0.0001827  0.01015            NA
## 122758:     GO & Domain 0.0006629   0.0001827  0.01015            NA
## 122759:     GO & Domain 0.0006629   0.0001827  0.01015            NA
## 122760:     GO & Domain 0.0006629   0.0001827  0.01015            NA
##         coef.colony:val coef.colony:I(val^2)    padj max.val min.val
##      1:              NA             0.002611 0.02545       0      21
##      2:              NA             0.002611 0.02545       0      21
##      3:              NA             0.002611 0.02545       0      21
##      4:              NA             0.002611 0.02545       0      21
##      5:              NA             0.002611 0.02545       0      21
##     ---                                                             
## 122756:         0.02664                   NA 0.01415       0    38.5
## 122757:         0.02664                   NA 0.01415       0    38.5
## 122758:         0.02664                   NA 0.01415       0    38.5
## 122759:         0.02664                   NA 0.01415       0    38.5
## 122760:         0.02664                   NA 0.01415       0    38.5
##         exp_type exp.scaled
##      1:      Low     4.1064
##      2:      Low    -0.3283
##      3:      Low    -0.3283
##      4:      Low    -0.3283
##      5:      Low    -0.3283
##     ---                    
## 122756:      Low    -0.3151
## 122757:      Low    -0.2822
## 122758:      Low    -1.0957
## 122759:      Low    -0.6743
## 122760:      Low    -1.0374
```

```r
str(Ap.dt)
```

```
## Classes 'data.table' and 'data.frame':	122760 obs. of  31 variables:
##  $ Transcript           : chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ Length               : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM                  : num  0.802 0 0 0 0 ...
##  $ RPKM                 : num  0.963 0 0 0 0 ...
##  $ KPKM                 : num  0.963 0 0 0 0 ...
##  $ EstimatedNumReads    : num  0.0919 0 0 0 0 ...
##  $ sample               : chr  "A22-0" "Ar-0" "A22-3" "Ar-3" ...
##  $ val                  : num  0 0 3.5 3.5 10.5 10.5 14 14 17.5 17.5 ...
##  $ colony               : Factor w/ 2 levels "A22","Ar": 1 2 1 2 1 2 1 2 1 2 ...
##  $ sequence.length      : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ best.hit.to.nr       : chr  "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " ...
##  $ hit.length           : chr  "25" "25" "25" "25" ...
##  $ E.value              : chr  "7.22e-06" "7.22e-06" "7.22e-06" "7.22e-06" ...
##  $ Bit.score            : chr  "57.514374" "57.514374" "57.514374" "57.514374" ...
##  $ GO.Biological.Process: chr  "-" "-" "-" "-" ...
##  $ GO.Cellular.Component: chr  "-" "-" "-" "-" ...
##  $ GO.Molecular.Function: chr  "-" "-" "-" "-" ...
##  $ Enzyme               : chr  "-" "-" "-" "-" ...
##  $ Domain               : chr  "-" "-" "-" "-" ...
##  $ annotation.type      : chr  "" "" "" "" ...
##  $ pval                 : num  0.00152 0.00152 0.00152 0.00152 0.00152 ...
##  $ coef.colony          : num  0.0457 0.0457 0.0457 0.0457 0.0457 ...
##  $ coef.val             : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.I(val^2)        : num  0.00258 0.00258 0.00258 0.00258 0.00258 ...
##  $ coef.colony:val      : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.colony:I(val^2) : num  0.00261 0.00261 0.00261 0.00261 0.00261 ...
##  $ padj                 : num  0.0254 0.0254 0.0254 0.0254 0.0254 ...
##  $ max.val              : chr  "0" "0" "0" "0" ...
##  $ min.val              : chr  "21" "21" "21" "21" ...
##  $ exp_type             : chr  "Low" "Low" "Low" "Low" ...
##  $ exp.scaled           : num  4.106 -0.328 -0.328 -0.328 -0.328 ...
##  - attr(*, "sorted")= chr "Transcript"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
write.csv(Ap.dt, file = paste(resultsdir, "Ap.dt.csv", sep = ""), quote = TRUE, 
    row.names = FALSE)

# subset to genes with significant interaction
Ap.dt.interaction <- Ap.dt[!is.na(Ap.dt$"coef.colony:val") | !is.na(Ap.dt$"coef.colony:I(val^2)")]
str(Ap.dt.interaction)
```

```
## Classes 'data.table' and 'data.frame':	82720 obs. of  31 variables:
##  $ Transcript           : chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ Length               : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM                  : num  0.802 0 0 0 0 ...
##  $ RPKM                 : num  0.963 0 0 0 0 ...
##  $ KPKM                 : num  0.963 0 0 0 0 ...
##  $ EstimatedNumReads    : num  0.0919 0 0 0 0 ...
##  $ sample               : chr  "A22-0" "Ar-0" "A22-3" "Ar-3" ...
##  $ val                  : num  0 0 3.5 3.5 10.5 10.5 14 14 17.5 17.5 ...
##  $ colony               : Factor w/ 2 levels "A22","Ar": 1 2 1 2 1 2 1 2 1 2 ...
##  $ sequence.length      : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ best.hit.to.nr       : chr  "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " ...
##  $ hit.length           : chr  "25" "25" "25" "25" ...
##  $ E.value              : chr  "7.22e-06" "7.22e-06" "7.22e-06" "7.22e-06" ...
##  $ Bit.score            : chr  "57.514374" "57.514374" "57.514374" "57.514374" ...
##  $ GO.Biological.Process: chr  "-" "-" "-" "-" ...
##  $ GO.Cellular.Component: chr  "-" "-" "-" "-" ...
##  $ GO.Molecular.Function: chr  "-" "-" "-" "-" ...
##  $ Enzyme               : chr  "-" "-" "-" "-" ...
##  $ Domain               : chr  "-" "-" "-" "-" ...
##  $ annotation.type      : chr  "" "" "" "" ...
##  $ pval                 : num  0.00152 0.00152 0.00152 0.00152 0.00152 ...
##  $ coef.colony          : num  0.0457 0.0457 0.0457 0.0457 0.0457 ...
##  $ coef.val             : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.I(val^2)        : num  0.00258 0.00258 0.00258 0.00258 0.00258 ...
##  $ coef.colony:val      : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.colony:I(val^2) : num  0.00261 0.00261 0.00261 0.00261 0.00261 ...
##  $ padj                 : num  0.0254 0.0254 0.0254 0.0254 0.0254 ...
##  $ max.val              : chr  "0" "0" "0" "0" ...
##  $ min.val              : chr  "21" "21" "21" "21" ...
##  $ exp_type             : chr  "Low" "Low" "Low" "Low" ...
##  $ exp.scaled           : num  4.106 -0.328 -0.328 -0.328 -0.328 ...
##  - attr(*, "sorted")= chr "Transcript"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
write.csv(Ap.dt.interaction, file = paste(resultsdir, "Ap.dt.interaction.csv", 
    sep = ""), quote = TRUE, row.names = FALSE)
```



## Gene set enrichment analysis ##

I use [topGO](http://www.bioconductor.org/packages/2.12/bioc/html/topGO.html) to perform gene set enrichment analysis (GSEA) seperately for each expression type (bimodal, intermediate, high, low).

First need to create gene ID to GO term map file


```r
# create geneid2go.map file from FastAnnotator AnnotationTable.txt
geneid2GOmap(annotationfile)
```


then read map file.


```r
# read mappings file
geneID2GO <- readMappings(file = "geneid2go.map")
```

```
## Warning: cannot open file 'geneid2go.map': No such file or directory
```

```
## Error: cannot open the connection
```

```r
str(head(geneID2GO))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'geneID2GO' not found
```


### GSEA for thermally-responsive transcripts ###

Using this gene2GO map file, perform GSEA for:

**1) all responsive transcripts**

Use `selectFDR` function to select transcripts with adjusted P < 0.05.


```r
# create geneList. note that NA values cause problems with topGO so set any
# NA to 1 as need to retain for GO analysis
Ap.geneList <- RxNout$padj
Ap.geneList[which(is.na(Ap.geneList))] <- 1
stopifnot(length(which(is.na(Ap.geneList))) == 0)
names(Ap.geneList) <- RxNout$Transcript
str(Ap.geneList)
```

```
##  Named num [1:105536] 0.836 0.195 0.93 0.921 0.9 ...
##  - attr(*, "names")= chr [1:105536] "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
```

```r

# Function to select top genes (defined above)
selectFDR <- function(padj) {
    return(padj < 0.05)
}

# create topGOdata object
Ap.BP.GOdata <- new("topGOdata", description = "BP gene set analysis", ontology = "BP", 
    allGenes = Ap.geneList, geneSel = selectFDR, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....
```

```
## Error: object 'geneID2GO' not found
```

```r

Ap.BP.GOdata
```

```
## Error: object 'Ap.BP.GOdata' not found
```

```r

# perform enrichment analysis using parentchild method
Ap.BP.resultParentChild <- runTest(Ap.BP.GOdata, statistic = "fisher", algorithm = "parentchild")
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'runTest': Error: object 'Ap.BP.GOdata' not found
```

```r
Ap.BP.resultParentChild
```

```
## Error: object 'Ap.BP.resultParentChild' not found
```

```r

# table results
Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, 
    topNodes = 118)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'GenTable': Error: object 'Ap.BP.GOdata' not found
```

```r
# Ap.BP.ResTable
write.table(Ap.BP.ResTable, file = paste(resultsdir, "Ap_GO.BP_results.txt", 
    sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
```

```
## Error: object 'Ap.BP.ResTable' not found
```

```r
pandoc.table(Ap.BP.ResTable)
```

```
## Error: object 'Ap.BP.ResTable' not found
```

```r

# graph significant nodes

pdf(paste(resultsdir, "Ap.BP_topGO_sig_nodes.pdf", sep = ""))
showSigOfNodes(Ap.BP.GOdata, score(Ap.BP.resultParentChild), firstSigNodes = 10, 
    useInfo = "all")
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'score': Error: object 'Ap.BP.resultParentChild' not found
```

```r
dev.off()
```

```
## pdf 
##   2
```


**2) High and low expressed transcripts**

Use `selectTranscript` function to select transcripts based on 'exp_type'. As only significant transcripts have an 'exp_type' assigned, this is a small subset of the above.


```r
selectTranscript <- function(score) {
    return(score == 1)
}

# To select genes, add 'exp_type' to RxN df
RxNout.t <- merge(RxNout, Ap.exp_type, by = "Transcript", all = TRUE)
length(which(!is.na(RxNout.t$exp_type)))
```

```
## [1] 5580
```




```r
# create geneList
Ap.bim.geneList <- rep(0, length = nrow(RxNout.t))
Ap.bim.geneList[which(RxNout.t$exp_type == "Bimodal")] <- 1
names(Ap.bim.geneList) <- RxNout.t$Transcript
str(Ap.bim.geneList)
```

```
##  Named num [1:105536] 0 0 0 0 0 0 0 0 0 0 ...
##  - attr(*, "names")= chr [1:105536] "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
```

```r
table(Ap.bim.geneList)
```

```
## Ap.bim.geneList
##      0      1 
## 103546   1990
```

```r

# create topGOdata object
Ap.bim.BP.GOdata <- new("topGOdata", description = "BP gene set analysis", ontology = "BP", 
    allGenes = Ap.bim.geneList, geneSel = selectTranscript, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....
```

```
## Error: object 'geneID2GO' not found
```

```r

Ap.bim.BP.GOdata
```

```
## Error: object 'Ap.bim.BP.GOdata' not found
```

```r

# perform enrichment analysis using parentchild method
Ap.bim.BP.resultParentChild <- runTest(Ap.bim.BP.GOdata, statistic = "fisher", 
    algorithm = "parentchild")
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'runTest': Error: object 'Ap.bim.BP.GOdata' not found
```

```r
Ap.bim.BP.resultParentChild
```

```
## Error: object 'Ap.bim.BP.resultParentChild' not found
```


Enriched gene sets for genes expressed at both high and low temperatures.


```
## Error: error in evaluating the argument 'object' in selecting a method for function 'GenTable': Error: object 'Ap.bim.BP.GOdata' not found
```

```
## Error: object 'Ap.bim.BP.ResTable' not found
```

```
## Error: object 'Ap.bim.BP.ResTable' not found
```



**3) Intermediate-temperature expressed transcripts**


```r
# create geneList
Ap.int.geneList <- rep(0, length = nrow(RxNout.t))
Ap.int.geneList[which(RxNout.t$exp_type == "Intermediate")] <- 1
names(Ap.int.geneList) <- RxNout.t$Transcript
str(Ap.int.geneList)
```

```
##  Named num [1:105536] 0 0 0 0 0 0 0 0 0 0 ...
##  - attr(*, "names")= chr [1:105536] "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
```

```r
table(Ap.int.geneList)
```

```
## Ap.int.geneList
##      0      1 
## 104826    710
```

```r

# create topGOdata object
Ap.int.BP.GOdata <- new("topGOdata", description = "BP gene set analysis", ontology = "BP", 
    allGenes = Ap.int.geneList, geneSel = selectTranscript, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....
```

```
## Error: object 'geneID2GO' not found
```

```r

Ap.int.BP.GOdata
```

```
## Error: object 'Ap.int.BP.GOdata' not found
```

```r

# perform enrichment analysis using parentchild method
Ap.int.BP.resultParentChild <- runTest(Ap.int.BP.GOdata, statistic = "fisher", 
    algorithm = "parentchild")
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'runTest': Error: object 'Ap.int.BP.GOdata' not found
```

```r
Ap.int.BP.resultParentChild
```

```
## Error: object 'Ap.int.BP.resultParentChild' not found
```


Enriched gene sets for genes expressed at intermediate temperatures.


```
## Error: error in evaluating the argument 'object' in selecting a method for function 'GenTable': Error: object 'Ap.int.BP.GOdata' not found
```

```
## Error: object 'Ap.int.BP.ResTable' not found
```

```
## Error: object 'Ap.int.BP.ResTable' not found
```

```
## Error: object 'Ap.int.BP.ResTable' not found
```



**4) High temperature expressed transcripts**



```r
# create geneList
Ap.hig.geneList <- rep(0, length = nrow(RxNout.t))
Ap.hig.geneList[which(RxNout.t$exp_type == "High")] <- 1
names(Ap.hig.geneList) <- RxNout.t$Transcript
str(Ap.hig.geneList)
```

```
##  Named num [1:105536] 0 0 0 0 0 0 0 0 0 0 ...
##  - attr(*, "names")= chr [1:105536] "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
```

```r
table(Ap.hig.geneList)
```

```
## Ap.hig.geneList
##      0      1 
## 104714    822
```

```r

# create topGOdata object
Ap.hig.BP.GOdata <- new("topGOdata", description = "BP gene set analysis", ontology = "BP", 
    allGenes = Ap.hig.geneList, geneSel = selectTranscript, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....
```

```
## Error: object 'geneID2GO' not found
```

```r

Ap.hig.BP.GOdata
```

```
## Error: object 'Ap.hig.BP.GOdata' not found
```

```r

# perform enrichment analysis using parentchild method
Ap.hig.BP.resultParentChild <- runTest(Ap.hig.BP.GOdata, statistic = "fisher", 
    algorithm = "parentchild")
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'runTest': Error: object 'Ap.hig.BP.GOdata' not found
```

```r
Ap.hig.BP.resultParentChild
```

```
## Error: object 'Ap.hig.BP.resultParentChild' not found
```


Enriched gene sets for genes expressed at both high temperatures only.


```
## Error: error in evaluating the argument 'object' in selecting a method for function 'GenTable': Error: object 'Ap.hig.BP.GOdata' not found
```

```
## Error: object 'Ap.hig.BP.ResTable' not found
```

```
## Error: object 'Ap.hig.BP.ResTable' not found
```



**5) Low temperature expressed transcripts**


```r
# create geneList
Ap.low.geneList <- rep(0, length = nrow(RxNout.t))
Ap.low.geneList[which(RxNout.t$exp_type == "Low")] <- 1
names(Ap.low.geneList) <- RxNout.t$Transcript
str(Ap.low.geneList)
```

```
##  Named num [1:105536] 0 0 0 0 0 0 0 0 0 0 ...
##  - attr(*, "names")= chr [1:105536] "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
```

```r
table(Ap.low.geneList)
```

```
## Ap.low.geneList
##      0      1 
## 103478   2058
```

```r

# create topGOdata object
Ap.low.BP.GOdata <- new("topGOdata", description = "BP gene set analysis", ontology = "BP", 
    allGenes = Ap.low.geneList, geneSel = selectTranscript, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....
```

```
## Error: object 'geneID2GO' not found
```

```r

Ap.low.BP.GOdata
```

```
## Error: object 'Ap.low.BP.GOdata' not found
```

```r

# perform enrichment analysis using parentchild method
Ap.low.BP.resultParentChild <- runTest(Ap.low.BP.GOdata, statistic = "fisher", 
    algorithm = "parentchild")
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'runTest': Error: object 'Ap.low.BP.GOdata' not found
```

```r
Ap.low.BP.resultParentChild
```

```
## Error: object 'Ap.low.BP.resultParentChild' not found
```


Enriched gene sets for genes expressed at low temperatures only.


```
## Error: error in evaluating the argument 'object' in selecting a method for function 'GenTable': Error: object 'Ap.low.BP.GOdata' not found
```

```
## Error: object 'Ap.low.BP.ResTable' not found
```

```
## Error: object 'Ap.low.BP.ResTable' not found
```

```
## Error: object 'Ap.low.BP.ResTable' not found
```



```
## Error: error in evaluating the argument 'x' in selecting a method for function 'score': Error: object 'Ap.bim.BP.resultParentChild' not found
```

```
## pdf 
##   2
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'score': Error: object 'Ap.int.BP.resultParentChild' not found
```

```
## pdf 
##   2
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'score': Error: object 'Ap.hig.BP.resultParentChild' not found
```

```
## pdf 
##   2
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'score': Error: object 'Ap.low.BP.resultParentChild' not found
```

```
## pdf 
##   2
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

A set of 3760 transcripts have expression patterns that depend on colony. In this section, I examine how these transcripts differ between colonies. Specifically, I expect genes that are up-regulated at high temperatures in the northern *A22* colony to be at constitutively higher expression in the more southern *Ar* colony. Conversely, I will examine the extent to which genes upregulated at colder temperatures in the *Ar* colony are constitutively expressed in the *A22* colony.








## Visualize responsive transcripts

Make plots for all significant transcripts




Make plots for transcripts with significant temperature by colony interaction














## Session information ##


```r
save.image()
```

```
## Warning: 'package:R.utils' may not be available when loading
```

```r
sessionInfo()
```

```
## R version 3.0.2 (2013-09-25)
## Platform: i686-pc-linux-gnu (32-bit)
## 
## locale:
##  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
##  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
##  [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] wordcloud_2.4        tm_0.5-10            Rcpp_0.10.6         
##  [4] RColorBrewer_1.0-5   Rgraphviz_2.6.0      topGO_2.14.0        
##  [7] SparseM_1.03         GO.db_2.10.1         RSQLite_0.11.4      
## [10] DBI_0.2-7            AnnotationDbi_1.24.0 Biobase_2.22.0      
## [13] BiocGenerics_0.8.0   graph_1.40.1         plyr_1.8            
## [16] RCurl_1.95-4.1       bitops_1.0-6         data.table_1.8.10   
## [19] stringr_0.6.2        pander_0.3.8         knitcitations_0.5-0 
## [22] bibtex_0.3-6         ggplot2_0.9.3.1      R.utils_1.28.4      
## [25] R.oo_1.17.0          R.methodsS3_1.6.1    knitr_1.5           
## 
## loaded via a namespace (and not attached):
##  [1] codetools_0.2-8  colorspace_1.2-4 dichromat_2.0-0  digest_0.6.4    
##  [5] evaluate_0.5.1   formatR_0.10     gtable_0.1.2     httr_0.2        
##  [9] IRanges_1.20.6   labeling_0.2     lattice_0.20-24  MASS_7.3-29     
## [13] munsell_0.4.2    proto_0.3-10     reshape2_1.2.2   scales_0.2.3    
## [17] slam_0.1-31      stats4_3.0.2     tools_3.0.2      XML_3.98-1.1    
## [21] xtable_1.7-1
```


## References


```r
bibliography()
```


- Manfred G Grabherr, Brian J Haas, Moran Yassour, Joshua Z Levin, Dawn A Thompson, Ido Amit, Xian Adiconis, Lin Fan, Raktima Raychowdhury, Qiandong Zeng, Zehua Chen, Evan Mauceli, Nir Hacohen, Andreas Gnirke, Nicholas Rhind, Federica di Palma, Bruce W Birren, Chad Nusbaum, Kerstin Lindblad-Toh, Nir Friedman, Aviv Regev,   (2011) Full-Length Transcriptome Assembly From Rna-Seq Data Without A Reference Genome.  *Nature Biotechnology*  **29**  644-652  [10.1038/nbt.1883](http://dx.doi.org/10.1038/nbt.1883)
- X. Huang,   (1999) Cap3: A Dna Sequence Assembly Program.  *Genome Research*  **9**  868-877  [10.1101/gr.9.9.868](http://dx.doi.org/10.1101/gr.9.9.868)
- B. Li, V. Ruotti, R. M. Stewart, J. A. Thomson, C. N. Dewey,   (2009) Rna-Seq Gene Expression Estimation With Read Mapping Uncertainty.  *Bioinformatics*  **26**  493-500  [10.1093/bioinformatics/btp692](http://dx.doi.org/10.1093/bioinformatics/btp692)
- M. Lohse, A. M. Bolger, A. Nagel, A. R. Fernie, J. E. Lunn, M. Stitt, B. Usadel,   (2012) Robina: A User-Friendly, Integrated Software Solution For Rna-Seq-Based Transcriptomics.  *Nucleic Acids Research*  **40**  W622-W627  [10.1093/nar/gks540](http://dx.doi.org/10.1093/nar/gks540)
- David Lubertazzi,   (2012) The Biology And Natural History of Aphaenogaster Rudis.  *Psyche: A Journal of Entomology*  **2012**  1-11  [10.1155/2012/752815](http://dx.doi.org/10.1155/2012/752815)
- unknown unknown,   (unknown) Unknown.  *Unknown*
- Ya Yang, Stephen A Smith,   (2013) Optimizing de Novo Assembly of Short-Read Rna-Seq Data For Phylogenomics.  *Bmc Genomics*  **14**  328-NA  [10.1186/1471-2164-14-328](http://dx.doi.org/10.1186/1471-2164-14-328)

