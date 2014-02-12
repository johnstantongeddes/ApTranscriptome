Thermal reactionome of the common ant species *Aphaenogaster*
================================================================
  
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**February 10, 2014**

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


```r
### Transcriptome assembly

### Expression data

### Annotation file from either AWS or GoogleDrive
annotationURL <- getURL("http://johnstantongeddes.org/assets/files/Aphaeno_transcriptome_AnnotationTable.txt")
# a2 <- getURL('https://googledrive.com/host/0B75IymziRJ_9Tlg1U1Vxbjk1bzg')
# # GoogleDrive link

annotationfile <- read.csv(textConnection(annotationURL), header = TRUE, sep = "\t", 
    stringsAsFactors = FALSE)
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



## Transcriptome annotation ##

Annotation was performed by uploading the reduced assembly "Trinity_cap3_uclust.fasta" to the web-based annotation program [FastAnnotator](http://fastannotator.cgu.edu.tw/index.php) (<a href="">unknown, unknown</a>).

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
    print(sailfishcmd)
    system(sailfishcmd)
    message("Done with expression quantification for sample ", samples[j], ": ", 
        Sys.time(), "\n")
}
```

```
## Start expression quantification for sample A22-0: 2014-02-11 14:49:51
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-0_quant --reads data/ind_files/A22-0_ATCACG.paired.left.fastq data/ind_files/A22-0_ATCACG.paired.right.fastq data/ind_files/A22-0_ATCACG.unpaired.left.fastq data/ind_files/A22-0_ATCACG.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-0: 2014-02-11 14:53:04
## 
## Start expression quantification for sample A22-10: 2014-02-11 14:53:04
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-10_quant --reads data/ind_files/A22-10_TGACCA.paired.left.fastq data/ind_files/A22-10_TGACCA.paired.right.fastq data/ind_files/A22-10_TGACCA.unpaired.left.fastq data/ind_files/A22-10_TGACCA.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-10: 2014-02-11 14:56:17
## 
## Start expression quantification for sample A22-14: 2014-02-11 14:56:17
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-14_quant --reads data/ind_files/A22-14_ACAGTG.paired.left.fastq data/ind_files/A22-14_ACAGTG.paired.right.fastq data/ind_files/A22-14_ACAGTG.unpaired.left.fastq data/ind_files/A22-14_ACAGTG.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-14: 2014-02-11 14:59:30
## 
## Start expression quantification for sample A22-17: 2014-02-11 14:59:30
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-17_quant --reads data/ind_files/A22-17_GCCAAT.paired.left.fastq data/ind_files/A22-17_GCCAAT.paired.right.fastq data/ind_files/A22-17_GCCAAT.unpaired.left.fastq data/ind_files/A22-17_GCCAAT.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-17: 2014-02-11 15:02:43
## 
## Start expression quantification for sample A22-21: 2014-02-11 15:02:43
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-21_quant --reads data/ind_files/A22-21_CAGATC.paired.left.fastq data/ind_files/A22-21_CAGATC.paired.right.fastq data/ind_files/A22-21_CAGATC.unpaired.left.fastq data/ind_files/A22-21_CAGATC.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-21: 2014-02-11 15:05:55
## 
## Start expression quantification for sample A22-24: 2014-02-11 15:05:55
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-24_quant --reads data/ind_files/A22-24_ACTTGA.paired.left.fastq data/ind_files/A22-24_ACTTGA.paired.right.fastq data/ind_files/A22-24_ACTTGA.unpaired.left.fastq data/ind_files/A22-24_ACTTGA.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-24: 2014-02-11 15:09:08
## 
## Start expression quantification for sample A22-28: 2014-02-11 15:09:08
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-28_quant --reads data/ind_files/A22-28_GATCAG.paired.left.fastq data/ind_files/A22-28_GATCAG.paired.right.fastq data/ind_files/A22-28_GATCAG.unpaired.left.fastq data/ind_files/A22-28_GATCAG.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-28: 2014-02-11 15:12:20
## 
## Start expression quantification for sample A22-31: 2014-02-11 15:12:20
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-31_quant --reads data/ind_files/A22-31_TAGCTT.paired.left.fastq data/ind_files/A22-31_TAGCTT.paired.right.fastq data/ind_files/A22-31_TAGCTT.unpaired.left.fastq data/ind_files/A22-31_TAGCTT.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-31: 2014-02-11 15:15:32
## 
## Start expression quantification for sample A22-35: 2014-02-11 15:15:32
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-35_quant --reads data/ind_files/A22-35_GGCTAC.paired.left.fastq data/ind_files/A22-35_GGCTAC.paired.right.fastq data/ind_files/A22-35_GGCTAC.unpaired.left.fastq data/ind_files/A22-35_GGCTAC.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-35: 2014-02-11 15:18:45
## 
## Start expression quantification for sample A22-38: 2014-02-11 15:18:45
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-38_quant --reads data/ind_files/A22-38_CTTGTA.paired.left.fastq data/ind_files/A22-38_CTTGTA.paired.right.fastq data/ind_files/A22-38_CTTGTA.unpaired.left.fastq data/ind_files/A22-38_CTTGTA.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-38: 2014-02-11 15:21:59
## 
## Start expression quantification for sample A22-3: 2014-02-11 15:21:59
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-3_quant --reads data/ind_files/A22-3_CGATGT.paired.left.fastq data/ind_files/A22-3_CGATGT.paired.right.fastq data/ind_files/A22-3_CGATGT.unpaired.left.fastq data/ind_files/A22-3_CGATGT.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-3: 2014-02-11 15:25:12
## 
## Start expression quantification for sample A22-7: 2014-02-11 15:25:12
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/A22-7_quant --reads data/ind_files/A22-7_TTAGGC.paired.left.fastq data/ind_files/A22-7_TTAGGC.paired.right.fastq data/ind_files/A22-7_TTAGGC.unpaired.left.fastq data/ind_files/A22-7_TTAGGC.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample A22-7: 2014-02-11 15:28:25
## 
## Start expression quantification for sample Ar-0: 2014-02-11 15:28:25
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-0_quant --reads data/ind_files/Ar-0_AGTCAA.paired.left.fastq data/ind_files/Ar-0_AGTCAA.paired.right.fastq data/ind_files/Ar-0_AGTCAA.unpaired.left.fastq data/ind_files/Ar-0_AGTCAA.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-0: 2014-02-11 15:31:38
## 
## Start expression quantification for sample Ar-10: 2014-02-11 15:31:38
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-10_quant --reads data/ind_files/Ar-10_CCGTCC.paired.left.fastq data/ind_files/Ar-10_CCGTCC.paired.right.fastq data/ind_files/Ar-10_CCGTCC.unpaired.left.fastq data/ind_files/Ar-10_CCGTCC.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-10: 2014-02-11 15:34:51
## 
## Start expression quantification for sample Ar-14: 2014-02-11 15:34:51
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-14_quant --reads data/ind_files/Ar-14_GTCCGC.paired.left.fastq data/ind_files/Ar-14_GTCCGC.paired.right.fastq data/ind_files/Ar-14_GTCCGC.unpaired.left.fastq data/ind_files/Ar-14_GTCCGC.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-14: 2014-02-11 15:38:04
## 
## Start expression quantification for sample Ar-17: 2014-02-11 15:38:04
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-17_quant --reads data/ind_files/Ar-17_GTGAAA.paired.left.fastq data/ind_files/Ar-17_GTGAAA.paired.right.fastq data/ind_files/Ar-17_GTGAAA.unpaired.left.fastq data/ind_files/Ar-17_GTGAAA.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-17: 2014-02-11 15:41:17
## 
## Start expression quantification for sample Ar-21: 2014-02-11 15:41:17
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-21_quant --reads data/ind_files/Ar-21_GTGGCC.paired.left.fastq data/ind_files/Ar-21_GTGGCC.paired.right.fastq data/ind_files/Ar-21_GTGGCC.unpaired.left.fastq data/ind_files/Ar-21_GTGGCC.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-21: 2014-02-11 15:44:30
## 
## Start expression quantification for sample Ar-24: 2014-02-11 15:44:30
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-24_quant --reads data/ind_files/Ar-24_GTTTCG.paired.left.fastq data/ind_files/Ar-24_GTTTCG.paired.right.fastq data/ind_files/Ar-24_GTTTCG.unpaired.left.fastq data/ind_files/Ar-24_GTTTCG.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-24: 2014-02-11 15:47:43
## 
## Start expression quantification for sample Ar-28: 2014-02-11 15:47:43
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-28_quant --reads data/ind_files/Ar-28_CGTACG.paired.left.fastq data/ind_files/Ar-28_CGTACG.paired.right.fastq data/ind_files/Ar-28_CGTACG.unpaired.left.fastq data/ind_files/Ar-28_CGTACG.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-28: 2014-02-11 15:50:55
## 
## Start expression quantification for sample Ar-31: 2014-02-11 15:50:55
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-31_quant --reads data/ind_files/Ar-31_GAGTGG.paired.left.fastq data/ind_files/Ar-31_GAGTGG.paired.right.fastq data/ind_files/Ar-31_GAGTGG.unpaired.left.fastq data/ind_files/Ar-31_GAGTGG.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-31: 2014-02-11 15:54:08
## 
## Start expression quantification for sample Ar-35: 2014-02-11 15:54:08
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-35_quant --reads data/ind_files/Ar-35_ACTGAT.paired.left.fastq data/ind_files/Ar-35_ACTGAT.paired.right.fastq data/ind_files/Ar-35_ACTGAT.unpaired.left.fastq data/ind_files/Ar-35_ACTGAT.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-35: 2014-02-11 15:57:20
## 
## Start expression quantification for sample Ar-38: 2014-02-11 15:57:20
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-38_quant --reads data/ind_files/Ar-38_ATTCCT.paired.left.fastq data/ind_files/Ar-38_ATTCCT.paired.right.fastq data/ind_files/Ar-38_ATTCCT.unpaired.left.fastq data/ind_files/Ar-38_ATTCCT.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-38: 2014-02-11 16:00:33
## 
## Start expression quantification for sample Ar-3: 2014-02-11 16:00:33
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-3_quant --reads data/ind_files/Ar-3_AGTTCC.paired.left.fastq data/ind_files/Ar-3_AGTTCC.paired.right.fastq data/ind_files/Ar-3_AGTTCC.unpaired.left.fastq data/ind_files/Ar-3_AGTTCC.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-3: 2014-02-11 16:03:46
## 
## Start expression quantification for sample Ar-7: 2014-02-11 16:03:46
```

```
## [1] "sailfish quant -i results/trinity-full/sailfish-index -o results/trinity-full/sailfish-expression/Ar-7_quant --reads data/ind_files/Ar-7_ATGTCA.paired.left.fastq data/ind_files/Ar-7_ATGTCA.paired.right.fastq data/ind_files/Ar-7_ATGTCA.unpaired.left.fastq data/ind_files/Ar-7_ATGTCA.unpaired.right.fastq -p 4"
```

```
## Done with expression quantification for sample Ar-7: 2014-02-11 16:06:58
```


This generated a directory for each sample

A22-00_trimclip_quant, A22-03_trimclip_quant, A22-07_trimclip_quant, A22-0_quant, A22-10_quant, A22-10_trimclip_quant, A22-14_quant, A22-14_trimclip_quant, A22-17_quant, A22-17_trimclip_quant, A22-21_quant, A22-21_trimclip_quant, A22-24_quant, A22-24_trimclip_quant, A22-28_quant, A22-28_trimclip_quant, A22-31_quant, A22-31_trimclip_quant, A22-35_quant, A22-35_trimclip_quant, A22-38_quant, A22-38_trimclip_quant, A22-3_quant, A22-7_quant, Ar-00_trimclip_quant, Ar-03_trimclip_quant, Ar-07_trimclip_quant, Ar-0_quant, Ar-10_quant, Ar-10_trimclip_quant, Ar-14_quant, Ar-14_trimclip_quant, Ar-17_quant, Ar-17_trimclip_quant, Ar-21_quant, Ar-21_trimclip_quant, Ar-24_quant, Ar-24_trimclip_quant, Ar-28_quant, Ar-28_trimclip_quant, Ar-31_quant, Ar-31_trimclip_quant, Ar-35_quant, Ar-35_trimclip_quant, Ar-38_quant, Ar-38_trimclip_quant, Ar-3_quant, Ar-7_quant

and within each directory there are the following r:

quant_bias_corrected.sf, quant.sf, reads.count_info, reads.sfc

The file *quant_bias_corrected.sf* contains the following columns, following a number of header lines:

1. Transcript ID
2. Transcript Length
3. Transcripts per Million (TPM): computed as described in (<a href="http://dx.doi.org/10.1093/bioinformatics/btp692">Li et al. 2009</a>), and is meant as an estimate of the number of transcripts, per million observed transcripts, originating from each isoform.
4. Reads Per Kilobase per Million mapped reads (RPKM): classic measure of relative transcript abundance, and is an estimate of the number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

The TPM column for each sample was extracted and combined into a matrix for each colony.




Note that expression levels at each temperature treatment are highly correlated between the two colonies.


|  temp  |  cors  |
|:------:|:------:|
|   0    |  0.99  |
|  3.5   |  0.98  |
|   7    |  0.99  |
|  10.5  |   1    |
|   14   |  0.99  |
|  17.5  |  0.98  |
|   21   |  0.98  |
|  24.5  |  0.99  |
|   28   |  0.99  |
|  31.5  |   1    |
|   35   |  0.99  |
|  38.5  |  0.99  |

Table: correlations between colonies at each temperature treatment


## Identification of thermally-responsive genes

To identify transcripts (roughly equivalent to genes) that show thermal responsiveness, I fit the following linear model to each transcript:

$$ TPM = \beta_0 + \beta_1(colony) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(colony * temp) + \\beta_5(colony * temp^2) + \epsilon $$

where TPM is transcripts per million. 

For this list of P-values, False Discovery Rate (FDR) is applied and q-values are calculated using the [qvalue]() package.

Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the analysis to be conservative). 


```r
A22.TPM[, `:=`(colony, "A22")]
```

```
##                      Transcript Length     TPM    RPKM    KPKM
##       1:         0|*|Contig6267   9990 0.08467 0.10151 0.10151
##       2:         0|*|Contig6267   9990 0.03685 0.05494 0.05494
##       3:         0|*|Contig6267   9990 0.07277 0.09098 0.09098
##       4:         0|*|Contig6267   9990 0.04037 0.05347 0.05347
##       5:         0|*|Contig6267   9990 0.02208 0.03365 0.03365
##      ---                                                      
## 1266428: 9|*|comp147140_c0_seq1   9030 0.77842 1.22784 1.22784
## 1266429: 9|*|comp147140_c0_seq1   9030 0.75617 1.19962 1.19962
## 1266430: 9|*|comp147140_c0_seq1   9030 1.03055 1.54239 1.54239
## 1266431: 9|*|comp147140_c0_seq1   9030 0.57364 0.85598 0.85598
## 1266432: 9|*|comp147140_c0_seq1   9030 0.41257 0.61694 0.61694
##          EstimatedNumReads sample  val colony
##       1:            0.5097  A22-0  0.0    A22
##       2:            0.2279  A22-3  3.5    A22
##       3:            0.4274  A22-7  7.0    A22
##       4:            0.2918 A22-10 10.5    A22
##       5:            0.1667 A22-14 14.0    A22
##      ---                                     
## 1266428:            4.9717 A22-24 24.5    A22
## 1266429:            3.6815 A22-28 28.0    A22
## 1266430:            5.3435 A22-31 31.5    A22
## 1266431:            3.1857 A22-35 35.0    A22
## 1266432:            2.8690 A22-38 38.5    A22
```

```r
Ar.TPM[, `:=`(colony, "Ar")]
```

```
##                      Transcript Length     TPM    RPKM    KPKM
##       1:         0|*|Contig6267   9990 0.04653 0.09752 0.09752
##       2:         0|*|Contig6267   9990 0.09165 0.15292 0.15292
##       3:         0|*|Contig6267   9990 0.13826 0.16710 0.16710
##       4:         0|*|Contig6267   9990 0.16976 0.31147 0.31147
##       5:         0|*|Contig6267   9990 0.14294 0.27476 0.27476
##      ---                                                      
## 1266428: 9|*|comp147140_c0_seq1   9030 0.52133 0.94870 0.94870
## 1266429: 9|*|comp147140_c0_seq1   9030 0.42225 0.83172 0.83172
## 1266430: 9|*|comp147140_c0_seq1   9030 0.55667 1.10908 1.10908
## 1266431: 9|*|comp147140_c0_seq1   9030 0.29284 0.64732 0.64732
## 1266432: 9|*|comp147140_c0_seq1   9030 0.30791 0.63980 0.63980
##          EstimatedNumReads sample  val colony
##       1:            0.1926   Ar-0  0.0     Ar
##       2:            0.3926   Ar-3  3.5     Ar
##       3:            0.6505   Ar-7  7.0     Ar
##       4:            0.7585  Ar-10 10.5     Ar
##       5:            0.8320  Ar-14 14.0     Ar
##      ---                                     
## 1266428:            2.2545  Ar-24 24.5     Ar
## 1266429:            1.7762  Ar-28 28.0     Ar
## 1266430:            2.3711  Ar-31 31.5     Ar
## 1266431:            1.0372  Ar-35 35.0     Ar
## 1266432:            1.3942  Ar-38 38.5     Ar
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
##  $ TPM              : num  0.0847 0.0368 0.0728 0.0404 0.0221 ...
##  $ RPKM             : num  0.1015 0.0549 0.091 0.0535 0.0337 ...
##  $ KPKM             : num  0.1015 0.0549 0.091 0.0535 0.0337 ...
##  $ EstimatedNumReads: num  0.51 0.228 0.427 0.292 0.167 ...
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


Of the 105536 transcripts, 22582 have models with P < 0.05.

Many of these are likely false positives, so I adjusted P-values using false discovery rate (FDR) to identify only those transcripts with less than 5% FDR as significant. 


```r
RxNout$qval <- p.adjust(RxNout$pval, method = "fdr")

# subset to significant transcripts
signif.transcripts <- RxNout[which(RxNout$qval < 0.05), ]
```


## Functional annotation


```r
# add annotation information
setkey(annotationtable, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
```



|   Coefficient    |  Number_significant  |
|:----------------:|:--------------------:|
|      Total       |         8753         |
|      Colony      |         7895         |
|     Temp.lin     |         4186         |
|    Temp.quad     |         2459         |
| Colony:Temp.lin  |         2617         |
| Colony:Temp.quad |         2364         |

Table: Number of transcripts for which each term is significant


Of these, subset to those that have significant responses to temperature, either through a direct effect or interaction with colony. Add annotation information and write results to file. Do the same for transcripts that differ in expression between the colonies.


```r
responsive.transcripts <- signif.transcripts[!is.na(signif.transcripts$coef.val) | 
    !is.na(signif.transcripts$"coef.I(val^2)") | !is.na(signif.transcripts$"coef.colony:val") | 
    !is.na(signif.transcripts$"coef.colony:I(val^2)")]
dim(responsive.transcripts)
```

```
## [1] 5540    8
```

```r

# join signif transcripts with annotation
responsive.transcripts.ann <- annotationtable[responsive.transcripts]
str(responsive.transcripts.ann)
```

```
## Classes 'data.table' and 'data.frame':	5540 obs. of  19 variables:
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
##  $ pval                 : num  0.001478 0.002161 0.002116 0.000632 0.002495 ...
##  $ coef.colony          : num  0.044884 NA NA 0.000452 0.001476 ...
##  $ coef.val             : num  NA NA 0.01803 0.00614 0.00457 ...
##  $ coef.I(val^2)        : num  0.00249 0.00244 NA NA NA ...
##  $ coef.colony:val      : num  NA NA 0.01101 0.00614 NA ...
##  $ coef.colony:I(val^2) : num  0.00252 0.00293 0.00383 NA NA ...
##  $ qval                 : num  0.025 0.0323 0.0319 0.0138 0.0355 ...
##  - attr(*, "sorted")= chr "Sequence.Name"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r

write.table(responsive.transcripts.ann, file = paste(resultsdir, "Ap_responsive_transcripts_GO.txt", 
    sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

colony.transcripts <- signif.transcripts[!is.na(signif.transcripts$coef.colony)]
colony.transcripts.ann <- annotationtable[colony.transcripts]
str(colony.transcripts.ann)
```

```
## Classes 'data.table' and 'data.frame':	7895 obs. of  19 variables:
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
##  $ pval                 : num  0.001478 0.000349 0.000972 0.000405 0.000632 ...
##  $ coef.colony          : num  4.49e-02 7.77e-06 7.28e-05 1.29e-05 4.52e-04 ...
##  $ coef.val             : num  NA NA NA NA 0.00614 ...
##  $ coef.I(val^2)        : num  0.00249 NA NA NA NA ...
##  $ coef.colony:val      : num  NA NA NA NA 0.00614 ...
##  $ coef.colony:I(val^2) : num  0.00252 NA NA NA NA ...
##  $ qval                 : num  0.02497 0.00915 0.01859 0.01015 0.01378 ...
##  - attr(*, "sorted")= chr "Sequence.Name"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
write.table(colony.transcripts.ann, file = paste(resultsdir, "Ap_colony_transcripts_GO.txt", 
    sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
```


Identify shape of curve:

* increase - exrpression increases with temperature
* decrease - expression decreases with temperature
* concave - expression greatest at intermediate temperatures
* convex - expression lowest at intermediate temperatures


|  concave  |  convex  |  decrease  |  increase  |
|:---------:|:--------:|:----------:|:----------:|
|    715    |   4272   |    511     |     42     |

Table: Number of transcripts that have maximum and minimum expression at each temperature level


While most transcripts have a convex shape, this summary masks that most of these have maximum expression at low or high, with zero to low elsewhere. Next, identify where maximum expression occurs.



|  &nbsp;   |  0   |  3.5  |  10  |  14  |  17.5  |  21  |  24.5  |
|:---------:|:----:|:-----:|:----:|:----:|:------:|:----:|:------:|
|  **Max**  | 3818 |  101  | 145  | 113  |  169   |  77  |   18   |
|  **Min**  | 196  |  22   | 155  | 557  |  264   | 461  |  655   |

Table: Number of transcripts that have maximum and minimum expression at each temperature level (continued below)

 

|  &nbsp;   |  28  |  31.5  |  35  |  38.5  |
|:---------:|:----:|:------:|:----:|:------:|
|  **Max**  |  7   |   9    |  9   |  1074  |
|  **Min**  | 675  |  379   | 244  |  1932  |


**Most** transcripts have maximum expression at high (>=31.5) or low (<=10) temperatures. As the genes invovled in each type of thermal response may differ, designate responsive transcripts by max expression at high, low or intermediate (14-28) temperatures.



```r
# set expression type
Ap.max.min.exp <- data.table(Ap.max.min.exp)
setkey(Ap.max.min.exp, Transcript)
setkey(Ap.dt, Transcript)
Ap.dt <- Ap.dt[Ap.max.min.exp]
Ap.dt[, `:=`(exp_type, "intermediate")]
```

```
##                           Transcript Length    TPM   RPKM   KPKM
##      1: 100015|*|comp3543055_c0_seq1    208 0.7914 0.9487 0.9487
##      2: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      3: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      4: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      5: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##     ---                                                         
## 121876:       9|*|comp147140_c0_seq1   9030 0.5567 1.1091 1.1091
## 121877:       9|*|comp147140_c0_seq1   9030 0.5736 0.8560 0.8560
## 121878:       9|*|comp147140_c0_seq1   9030 0.2928 0.6473 0.6473
## 121879:       9|*|comp147140_c0_seq1   9030 0.4126 0.6169 0.6169
## 121880:       9|*|comp147140_c0_seq1   9030 0.3079 0.6398 0.6398
##         EstimatedNumReads sample  val colony sequence.length
##      1:           0.09029  A22-0  0.0    A22             208
##      2:           0.00000   Ar-0  0.0     Ar             208
##      3:           0.00000  A22-3  3.5    A22             208
##      4:           0.00000   Ar-3  3.5     Ar             208
##      5:           0.00000 A22-10 10.5    A22             208
##     ---                                                     
## 121876:           2.37111  Ar-31 31.5     Ar            9030
## 121877:           3.18572 A22-35 35.0    A22            9030
## 121878:           1.03722  Ar-35 35.0     Ar            9030
## 121879:           2.86903 A22-38 38.5    A22            9030
## 121880:           1.39419  Ar-38 38.5     Ar            9030
##                                                         best.hit.to.nr
##      1: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      2: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      3: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      4: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      5: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##     ---                                                               
## 121876:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121877:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121878:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121879:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121880:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
##         hit.length  E.value   Bit.score
##      1:         25 7.22e-06   57.514374
##      2:         25 7.22e-06   57.514374
##      3:         25 7.22e-06   57.514374
##      4:         25 7.22e-06   57.514374
##      5:         25 7.22e-06   57.514374
##     ---                                
## 121876:       2101      0.0 4740.368306
## 121877:       2101      0.0 4740.368306
## 121878:       2101      0.0 4740.368306
## 121879:       2101      0.0 4740.368306
## 121880:       2101      0.0 4740.368306
##                                                                                                                                                                                                          GO.Biological.Process
##      1:                                                                                                                                                                                                                      -
##      2:                                                                                                                                                                                                                      -
##      3:                                                                                                                                                                                                                      -
##      4:                                                                                                                                                                                                                      -
##      5:                                                                                                                                                                                                                      -
##     ---                                                                                                                                                                                                                       
## 121876: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121877: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121878: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121879: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121880: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
##                                                     GO.Cellular.Component
##      1:                                                                 -
##      2:                                                                 -
##      3:                                                                 -
##      4:                                                                 -
##      5:                                                                 -
##     ---                                                                  
## 121876: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121877: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121878: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121879: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121880: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
##                                                                                                                         GO.Molecular.Function
##      1:                                                                                                                                     -
##      2:                                                                                                                                     -
##      3:                                                                                                                                     -
##      4:                                                                                                                                     -
##      5:                                                                                                                                     -
##     ---                                                                                                                                      
## 121876: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121877: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121878: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121879: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121880: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
##         Enzyme
##      1:      -
##      2:      -
##      3:      -
##      4:      -
##      5:      -
##     ---       
## 121876:      -
## 121877:      -
## 121878:      -
## 121879:      -
## 121880:      -
##                                                                                                                                                     Domain
##      1:                                                                                                                                                  -
##      2:                                                                                                                                                  -
##      3:                                                                                                                                                  -
##      4:                                                                                                                                                  -
##      5:                                                                                                                                                  -
##     ---                                                                                                                                                   
## 121876: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121877: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121878: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121879: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121880: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
##         annotation.type      pval coef.colony coef.val coef.I(val^2)
##      1:                 0.0014775   0.0448835       NA      0.002491
##      2:                 0.0014775   0.0448835       NA      0.002491
##      3:                 0.0014775   0.0448835       NA      0.002491
##      4:                 0.0014775   0.0448835       NA      0.002491
##      5:                 0.0014775   0.0448835       NA      0.002491
##     ---                                                             
## 121876:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121877:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121878:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121879:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121880:     GO & Domain 0.0006044   0.0001822 0.008766            NA
##         coef.colony:val coef.colony:I(val^2)    qval max.val min.val
##      1:              NA             0.002523 0.02497       0    21.0
##      2:              NA             0.002523 0.02497       0    21.0
##      3:              NA             0.002523 0.02497       0    21.0
##      4:              NA             0.002523 0.02497       0    21.0
##      5:              NA             0.002523 0.02497       0    21.0
##     ---                                                             
## 121876:         0.02394                   NA 0.01335       0    38.5
## 121877:         0.02394                   NA 0.01335       0    38.5
## 121878:         0.02394                   NA 0.01335       0    38.5
## 121879:         0.02394                   NA 0.01335       0    38.5
## 121880:         0.02394                   NA 0.01335       0    38.5
##             exp_type
##      1: intermediate
##      2: intermediate
##      3: intermediate
##      4: intermediate
##      5: intermediate
##     ---             
## 121876: intermediate
## 121877: intermediate
## 121878: intermediate
## 121879: intermediate
## 121880: intermediate
```

```r
Ap.dt[max.val > 31, `:=`(exp_type, "high")]
```

```
##                           Transcript Length    TPM   RPKM   KPKM
##      1: 100015|*|comp3543055_c0_seq1    208 0.7914 0.9487 0.9487
##      2: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      3: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      4: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      5: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##     ---                                                         
## 121876:       9|*|comp147140_c0_seq1   9030 0.5567 1.1091 1.1091
## 121877:       9|*|comp147140_c0_seq1   9030 0.5736 0.8560 0.8560
## 121878:       9|*|comp147140_c0_seq1   9030 0.2928 0.6473 0.6473
## 121879:       9|*|comp147140_c0_seq1   9030 0.4126 0.6169 0.6169
## 121880:       9|*|comp147140_c0_seq1   9030 0.3079 0.6398 0.6398
##         EstimatedNumReads sample  val colony sequence.length
##      1:           0.09029  A22-0  0.0    A22             208
##      2:           0.00000   Ar-0  0.0     Ar             208
##      3:           0.00000  A22-3  3.5    A22             208
##      4:           0.00000   Ar-3  3.5     Ar             208
##      5:           0.00000 A22-10 10.5    A22             208
##     ---                                                     
## 121876:           2.37111  Ar-31 31.5     Ar            9030
## 121877:           3.18572 A22-35 35.0    A22            9030
## 121878:           1.03722  Ar-35 35.0     Ar            9030
## 121879:           2.86903 A22-38 38.5    A22            9030
## 121880:           1.39419  Ar-38 38.5     Ar            9030
##                                                         best.hit.to.nr
##      1: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      2: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      3: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      4: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      5: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##     ---                                                               
## 121876:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121877:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121878:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121879:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121880:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
##         hit.length  E.value   Bit.score
##      1:         25 7.22e-06   57.514374
##      2:         25 7.22e-06   57.514374
##      3:         25 7.22e-06   57.514374
##      4:         25 7.22e-06   57.514374
##      5:         25 7.22e-06   57.514374
##     ---                                
## 121876:       2101      0.0 4740.368306
## 121877:       2101      0.0 4740.368306
## 121878:       2101      0.0 4740.368306
## 121879:       2101      0.0 4740.368306
## 121880:       2101      0.0 4740.368306
##                                                                                                                                                                                                          GO.Biological.Process
##      1:                                                                                                                                                                                                                      -
##      2:                                                                                                                                                                                                                      -
##      3:                                                                                                                                                                                                                      -
##      4:                                                                                                                                                                                                                      -
##      5:                                                                                                                                                                                                                      -
##     ---                                                                                                                                                                                                                       
## 121876: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121877: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121878: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121879: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121880: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
##                                                     GO.Cellular.Component
##      1:                                                                 -
##      2:                                                                 -
##      3:                                                                 -
##      4:                                                                 -
##      5:                                                                 -
##     ---                                                                  
## 121876: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121877: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121878: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121879: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121880: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
##                                                                                                                         GO.Molecular.Function
##      1:                                                                                                                                     -
##      2:                                                                                                                                     -
##      3:                                                                                                                                     -
##      4:                                                                                                                                     -
##      5:                                                                                                                                     -
##     ---                                                                                                                                      
## 121876: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121877: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121878: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121879: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121880: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
##         Enzyme
##      1:      -
##      2:      -
##      3:      -
##      4:      -
##      5:      -
##     ---       
## 121876:      -
## 121877:      -
## 121878:      -
## 121879:      -
## 121880:      -
##                                                                                                                                                     Domain
##      1:                                                                                                                                                  -
##      2:                                                                                                                                                  -
##      3:                                                                                                                                                  -
##      4:                                                                                                                                                  -
##      5:                                                                                                                                                  -
##     ---                                                                                                                                                   
## 121876: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121877: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121878: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121879: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121880: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
##         annotation.type      pval coef.colony coef.val coef.I(val^2)
##      1:                 0.0014775   0.0448835       NA      0.002491
##      2:                 0.0014775   0.0448835       NA      0.002491
##      3:                 0.0014775   0.0448835       NA      0.002491
##      4:                 0.0014775   0.0448835       NA      0.002491
##      5:                 0.0014775   0.0448835       NA      0.002491
##     ---                                                             
## 121876:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121877:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121878:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121879:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121880:     GO & Domain 0.0006044   0.0001822 0.008766            NA
##         coef.colony:val coef.colony:I(val^2)    qval max.val min.val
##      1:              NA             0.002523 0.02497       0    21.0
##      2:              NA             0.002523 0.02497       0    21.0
##      3:              NA             0.002523 0.02497       0    21.0
##      4:              NA             0.002523 0.02497       0    21.0
##      5:              NA             0.002523 0.02497       0    21.0
##     ---                                                             
## 121876:         0.02394                   NA 0.01335       0    38.5
## 121877:         0.02394                   NA 0.01335       0    38.5
## 121878:         0.02394                   NA 0.01335       0    38.5
## 121879:         0.02394                   NA 0.01335       0    38.5
## 121880:         0.02394                   NA 0.01335       0    38.5
##             exp_type
##      1: intermediate
##      2: intermediate
##      3: intermediate
##      4: intermediate
##      5: intermediate
##     ---             
## 121876: intermediate
## 121877: intermediate
## 121878: intermediate
## 121879: intermediate
## 121880: intermediate
```

```r
Ap.dt[max.val < 11, `:=`(exp_type, "low")]
```

```
##                           Transcript Length    TPM   RPKM   KPKM
##      1: 100015|*|comp3543055_c0_seq1    208 0.7914 0.9487 0.9487
##      2: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      3: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      4: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      5: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##     ---                                                         
## 121876:       9|*|comp147140_c0_seq1   9030 0.5567 1.1091 1.1091
## 121877:       9|*|comp147140_c0_seq1   9030 0.5736 0.8560 0.8560
## 121878:       9|*|comp147140_c0_seq1   9030 0.2928 0.6473 0.6473
## 121879:       9|*|comp147140_c0_seq1   9030 0.4126 0.6169 0.6169
## 121880:       9|*|comp147140_c0_seq1   9030 0.3079 0.6398 0.6398
##         EstimatedNumReads sample  val colony sequence.length
##      1:           0.09029  A22-0  0.0    A22             208
##      2:           0.00000   Ar-0  0.0     Ar             208
##      3:           0.00000  A22-3  3.5    A22             208
##      4:           0.00000   Ar-3  3.5     Ar             208
##      5:           0.00000 A22-10 10.5    A22             208
##     ---                                                     
## 121876:           2.37111  Ar-31 31.5     Ar            9030
## 121877:           3.18572 A22-35 35.0    A22            9030
## 121878:           1.03722  Ar-35 35.0     Ar            9030
## 121879:           2.86903 A22-38 38.5    A22            9030
## 121880:           1.39419  Ar-38 38.5     Ar            9030
##                                                         best.hit.to.nr
##      1: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      2: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      3: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      4: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      5: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##     ---                                                               
## 121876:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121877:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121878:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121879:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121880:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
##         hit.length  E.value   Bit.score
##      1:         25 7.22e-06   57.514374
##      2:         25 7.22e-06   57.514374
##      3:         25 7.22e-06   57.514374
##      4:         25 7.22e-06   57.514374
##      5:         25 7.22e-06   57.514374
##     ---                                
## 121876:       2101      0.0 4740.368306
## 121877:       2101      0.0 4740.368306
## 121878:       2101      0.0 4740.368306
## 121879:       2101      0.0 4740.368306
## 121880:       2101      0.0 4740.368306
##                                                                                                                                                                                                          GO.Biological.Process
##      1:                                                                                                                                                                                                                      -
##      2:                                                                                                                                                                                                                      -
##      3:                                                                                                                                                                                                                      -
##      4:                                                                                                                                                                                                                      -
##      5:                                                                                                                                                                                                                      -
##     ---                                                                                                                                                                                                                       
## 121876: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121877: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121878: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121879: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121880: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
##                                                     GO.Cellular.Component
##      1:                                                                 -
##      2:                                                                 -
##      3:                                                                 -
##      4:                                                                 -
##      5:                                                                 -
##     ---                                                                  
## 121876: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121877: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121878: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121879: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121880: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
##                                                                                                                         GO.Molecular.Function
##      1:                                                                                                                                     -
##      2:                                                                                                                                     -
##      3:                                                                                                                                     -
##      4:                                                                                                                                     -
##      5:                                                                                                                                     -
##     ---                                                                                                                                      
## 121876: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121877: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121878: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121879: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121880: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
##         Enzyme
##      1:      -
##      2:      -
##      3:      -
##      4:      -
##      5:      -
##     ---       
## 121876:      -
## 121877:      -
## 121878:      -
## 121879:      -
## 121880:      -
##                                                                                                                                                     Domain
##      1:                                                                                                                                                  -
##      2:                                                                                                                                                  -
##      3:                                                                                                                                                  -
##      4:                                                                                                                                                  -
##      5:                                                                                                                                                  -
##     ---                                                                                                                                                   
## 121876: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121877: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121878: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121879: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121880: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
##         annotation.type      pval coef.colony coef.val coef.I(val^2)
##      1:                 0.0014775   0.0448835       NA      0.002491
##      2:                 0.0014775   0.0448835       NA      0.002491
##      3:                 0.0014775   0.0448835       NA      0.002491
##      4:                 0.0014775   0.0448835       NA      0.002491
##      5:                 0.0014775   0.0448835       NA      0.002491
##     ---                                                             
## 121876:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121877:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121878:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121879:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121880:     GO & Domain 0.0006044   0.0001822 0.008766            NA
##         coef.colony:val coef.colony:I(val^2)    qval max.val min.val
##      1:              NA             0.002523 0.02497       0    21.0
##      2:              NA             0.002523 0.02497       0    21.0
##      3:              NA             0.002523 0.02497       0    21.0
##      4:              NA             0.002523 0.02497       0    21.0
##      5:              NA             0.002523 0.02497       0    21.0
##     ---                                                             
## 121876:         0.02394                   NA 0.01335       0    38.5
## 121877:         0.02394                   NA 0.01335       0    38.5
## 121878:         0.02394                   NA 0.01335       0    38.5
## 121879:         0.02394                   NA 0.01335       0    38.5
## 121880:         0.02394                   NA 0.01335       0    38.5
##         exp_type
##      1:      low
##      2:      low
##      3:      low
##      4:      low
##      5:      low
##     ---         
## 121876:      low
## 121877:      low
## 121878:      low
## 121879:      low
## 121880:      low
```

```r
head(Ap.dt)
```

```
##                      Transcript Length    TPM   RPKM   KPKM
## 1: 100015|*|comp3543055_c0_seq1    208 0.7914 0.9487 0.9487
## 2: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
## 3: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
## 4: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
## 5: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
## 6: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##    EstimatedNumReads sample  val colony sequence.length
## 1:           0.09029  A22-0  0.0    A22             208
## 2:           0.00000   Ar-0  0.0     Ar             208
## 3:           0.00000  A22-3  3.5    A22             208
## 4:           0.00000   Ar-3  3.5     Ar             208
## 5:           0.00000 A22-10 10.5    A22             208
## 6:           0.00000  Ar-10 10.5     Ar             208
##                                                    best.hit.to.nr
## 1: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
## 2: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
## 3: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
## 4: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
## 5: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
## 6: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##    hit.length  E.value Bit.score GO.Biological.Process
## 1:         25 7.22e-06 57.514374                     -
## 2:         25 7.22e-06 57.514374                     -
## 3:         25 7.22e-06 57.514374                     -
## 4:         25 7.22e-06 57.514374                     -
## 5:         25 7.22e-06 57.514374                     -
## 6:         25 7.22e-06 57.514374                     -
##    GO.Cellular.Component GO.Molecular.Function Enzyme Domain
## 1:                     -                     -      -      -
## 2:                     -                     -      -      -
## 3:                     -                     -      -      -
## 4:                     -                     -      -      -
## 5:                     -                     -      -      -
## 6:                     -                     -      -      -
##    annotation.type     pval coef.colony coef.val coef.I(val^2)
## 1:                 0.001478     0.04488       NA      0.002491
## 2:                 0.001478     0.04488       NA      0.002491
## 3:                 0.001478     0.04488       NA      0.002491
## 4:                 0.001478     0.04488       NA      0.002491
## 5:                 0.001478     0.04488       NA      0.002491
## 6:                 0.001478     0.04488       NA      0.002491
##    coef.colony:val coef.colony:I(val^2)    qval max.val min.val exp_type
## 1:              NA             0.002523 0.02497       0      21      low
## 2:              NA             0.002523 0.02497       0      21      low
## 3:              NA             0.002523 0.02497       0      21      low
## 4:              NA             0.002523 0.02497       0      21      low
## 5:              NA             0.002523 0.02497       0      21      low
## 6:              NA             0.002523 0.02497       0      21      low
```



## Functional annotation ##

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  


## Gene set enrichment analysis ##

I use [topGO](http://www.bioconductor.org/packages/2.12/bioc/html/topGO.html) to perform gene set enrichment analysis

First need to create gene ID to GO term map file


```r
# create geneid2go.map file from FastAnnotator AnnotationTable.txt
geneid2GOmap(annotationfile)
```


Using this gene2GO map file, perform gene set enrichment analysis.

Provide qvalues as the gene score, and select genes with q < 0.05 using custom `selectFDR` function.


```r
# read mappings file
geneID2GO <- readMappings(file = "geneid2go.map")
str(head(geneID2GO))
```

```
## List of 6
##  $ 0|*|Contig6267        : chr [1:6] "GO:0035335" "GO:0000188" "GO:0006570" "GO:0017017" ...
##  $ 1|*|comp150820_c2_seq6: chr [1:6] "GO:0030036" "GO:0015074" "GO:0003676" "GO:0003779" ...
##  $ 2|*|Contig6262        : chr [1:6] "GO:0035335" "GO:0000188" "GO:0006570" "GO:0017017" ...
##  $ 3|*|comp149397_c1_seq2: chr [1:4] "GO:0006508" "GO:0005634" "GO:0003677" "GO:0004252"
##  $ 4|*|Contig4755        : chr [1:10] "GO:0055114" "GO:0006355" "GO:0009395" "GO:0005634" ...
##  $ 5|*|Contig3727        : chr [1:7] "GO:0007269" "GO:0050803" "GO:0048488" "GO:0042967" ...
```

```r

# create geneList. note that NA values cause problems with topGO so set any
# NA to 1 as need to retain for GO analysis
Ap.geneList <- RxNout$pval
Ap.geneList[which(is.na(Ap.geneList))] <- 1
stopifnot(length(which(is.na(Ap.geneList))) == 0)
names(Ap.geneList) <- RxNout$Transcript
str(Ap.geneList)
```

```
##  Named num [1:105536] 0.6045 0.0368 0.9005 0.8826 0.8194 ...
##  - attr(*, "names")= chr [1:105536] "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
```

```r

# Function to select top genes (defined above)
selectFDR <- function(qvalue) {
    return(qvalue < 0.5)
}

# create topGOdata object
Ap.BP.GOdata <- new("topGOdata", description = "BP gene set analysis", ontology = "BP", 
    allGenes = Ap.geneList, geneSel = selectFDR, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5600 GO terms found. )
## 
## Build GO DAG topology ..........	( 9054 GO terms and 19829 relations. )
## 
## Annotating nodes ...............	( 33173 genes annotated to the GO terms. )
```

```r

Ap.BP.GOdata
```

```
## 
## ------------------------- topGOdata object -------------------------
## 
##  Description:
##    -  BP gene set analysis 
## 
##  Ontology:
##    -  BP 
## 
##  105536 available genes (all genes from the array):
##    - symbol:  0|*|Contig6267 100000|*|comp2663136_c0_seq1 100001|*|comp3439067_c0_seq1 100002|*|comp2050457_c0_seq1 100003|*|comp2029723_c0_seq1  ...
```

```
## Error: invalid 'digits' argument
```

```r

# perform enrichment analysis using multiple methods
Ap.BP.resultParentChild <- runTest(Ap.BP.GOdata, statistic = "fisher", algorithm = "parentchild")
```

```
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 3585 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	3 nodes to be scored.
## 
## 	 Level 16:	7 nodes to be scored.
## 
## 	 Level 15:	14 nodes to be scored.
## 
## 	 Level 14:	33 nodes to be scored.
## 
## 	 Level 13:	68 nodes to be scored.
## 
## 	 Level 12:	171 nodes to be scored.
## 
## 	 Level 11:	252 nodes to be scored.
## 
## 	 Level 10:	383 nodes to be scored.
## 
## 	 Level 9:	502 nodes to be scored.
## 
## 	 Level 8:	519 nodes to be scored.
## 
## 	 Level 7:	523 nodes to be scored.
## 
## 	 Level 6:	488 nodes to be scored.
## 
## 	 Level 5:	365 nodes to be scored.
## 
## 	 Level 4:	185 nodes to be scored.
## 
## 	 Level 3:	53 nodes to be scored.
## 
## 	 Level 2:	18 nodes to be scored.
```

```r
Ap.BP.resultParentChild
```

```
## 
## Description: BP gene set analysis 
## Ontology: BP 
## 'parentchild' algorithm with the 'fisher : joinFun = union' test
## 3585 GO terms scored: 415 terms with p < 0.01
## Annotation data:
##     Annotated genes: 33173 
##     Significant genes: 21172 
##     Min. no. of genes annotated to a GO: 10 
##     Nontrivial nodes: 3585
```

```r

Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, 
    topNodes = 10)
Ap.BP.ResTable
```

```
##         GO.ID                                        Term Annotated
## 1  GO:0032501            multicellular organismal process      3708
## 2  GO:0032502                       developmental process      3067
## 3  GO:0044707       single-multicellular organism process      3554
## 4  GO:0019538                   protein metabolic process      7991
## 5  GO:0000003                                reproduction      1049
## 6  GO:0044767       single-organism developmental process      2490
## 7  GO:0071840 cellular component organization or bioge...      4671
## 8  GO:0022414                        reproductive process       832
## 9  GO:0006996                      organelle organization      1983
## 10 GO:0048519 negative regulation of biological proces...      1250
##    Significant Expected parentchild
## 1         2993   2366.6      <1e-30
## 2         2472   1957.5      <1e-30
## 3         2855   2268.3      <1e-30
## 4         5695   5100.1      <1e-30
## 5          904    669.5      <1e-30
## 6         1988   1589.2      <1e-30
## 7         3431   2981.2      <1e-30
## 8          720    531.0      <1e-30
## 9         1612   1265.6      <1e-30
## 10        1022    797.8      <1e-30
```

```r
write.table(Ap.BP.ResTable, file = paste(resultsdir, "Ap_GO.BP_results.txt", 
    sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
pandoc.table(Ap.BP.ResTable)
```

```
## 
## ------------------------------------------------------------------
##   GO.ID                Term               Annotated   Significant 
## ---------- ----------------------------- ----------- -------------
## GO:0032501   multicellular organismal       3708         2993     
##                  process                                          
## 
## GO:0032502     developmental process        3067         2472     
## 
## GO:0044707 single-multicellular organism    3554         2855     
##                    process                                        
## 
## GO:0019538   protein metabolic process      7991         5695     
## 
## GO:0000003         reproduction             1049          904     
## 
## GO:0044767 single-organism developmental    2490         1988     
##                    process                                        
## 
## GO:0071840      cellular component          4671         3431     
##              organization or bioge...                             
## 
## GO:0022414     reproductive process          832          720     
## 
## GO:0006996    organelle organization        1983         1612     
## 
## GO:0048519    negative regulation of        1250         1022     
##                biological proces...                               
## ------------------------------------------------------------------
## 
## Table: Table continues below
## 
##  
## ------------------------
##  Expected   parentchild 
## ---------- -------------
##    2367       <1e-30    
## 
##    1957       <1e-30    
## 
##    2268       <1e-30    
## 
##    5100       <1e-30    
## 
##   669.5       <1e-30    
## 
##    1589       <1e-30    
## 
##    2981       <1e-30    
## 
##    531        <1e-30    
## 
##    1266       <1e-30    
## 
##   797.8       <1e-30    
## ------------------------
```

```r

# graph significant nodes

pdf(paste(resultsdir, "Ap.BP_topGO_sig_nodes.pdf", sep = ""))
showSigOfNodes(Ap.BP.GOdata, score(Ap.BP.resultParentChild), firstSigNodes = 10, 
    useInfo = "all")
```

```
## $dag
## A graphNEL graph with directed edges
## Number of Nodes = 21 
## Number of Edges = 28 
## 
## $complete.dag
## [1] "A graph with 21 nodes."
```

```r
dev.off()
```

```
## pdf 
##   2
```


GO analysis for cellular component


```r
# create topGOdata object
Ap.CC.GOdata <- new("topGOdata", description = "CC gene set analysis", ontology = "CC", 
    allGenes = Ap.geneList, geneSel = selectFDR, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 1039 GO terms found. )
## 
## Build GO DAG topology ..........	( 1279 GO terms and 2490 relations. )
## 
## Annotating nodes ...............	( 22645 genes annotated to the GO terms. )
```

```r

Ap.CC.GOdata
```

```
## 
## ------------------------- topGOdata object -------------------------
## 
##  Description:
##    -  CC gene set analysis 
## 
##  Ontology:
##    -  CC 
## 
##  105536 available genes (all genes from the array):
##    - symbol:  0|*|Contig6267 100000|*|comp2663136_c0_seq1 100001|*|comp3439067_c0_seq1 100002|*|comp2050457_c0_seq1 100003|*|comp2029723_c0_seq1  ...
```

```
## Error: invalid 'digits' argument
```

```r

# perform enrichment analysis using multiple methods
Ap.CC.resultParentChild <- runTest(Ap.CC.GOdata, statistic = "fisher", algorithm = "parentchild")
```

```
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 606 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 15:	2 nodes to be scored.
## 
## 	 Level 14:	17 nodes to be scored.
## 
## 	 Level 13:	27 nodes to be scored.
## 
## 	 Level 12:	41 nodes to be scored.
## 
## 	 Level 11:	70 nodes to be scored.
## 
## 	 Level 10:	82 nodes to be scored.
## 
## 	 Level 9:	66 nodes to be scored.
## 
## 	 Level 8:	85 nodes to be scored.
## 
## 	 Level 7:	45 nodes to be scored.
## 
## 	 Level 6:	46 nodes to be scored.
## 
## 	 Level 5:	38 nodes to be scored.
## 
## 	 Level 4:	55 nodes to be scored.
## 
## 	 Level 3:	19 nodes to be scored.
## 
## 	 Level 2:	12 nodes to be scored.
```

```r
Ap.CC.resultParentChild
```

```
## 
## Description: CC gene set analysis 
## Ontology: CC 
## 'parentchild' algorithm with the 'fisher : joinFun = union' test
## 606 GO terms scored: 99 terms with p < 0.01
## Annotation data:
##     Annotated genes: 22645 
##     Significant genes: 14905 
##     Min. no. of genes annotated to a GO: 10 
##     Nontrivial nodes: 606
```

```r

Ap.CC.ResTable <- GenTable(Ap.CC.GOdata, parentchild = Ap.CC.resultParentChild, 
    topNodes = 40)
Ap.CC.ResTable
```

```
##         GO.ID                                        Term Annotated
## 1  GO:0043226                                   organelle      9400
## 2  GO:0043229                     intracellular organelle      9392
## 3  GO:0044422                              organelle part      4810
## 4  GO:0005622                               intracellular     14281
## 5  GO:0044424                          intracellular part     13829
## 6  GO:0044446                intracellular organelle part      4703
## 7  GO:0005811                              lipid particle       277
## 8  GO:0031974                     membrane-enclosed lumen      1363
## 9  GO:0044444                            cytoplasmic part      5992
## 10 GO:0005875              microtubule associated complex       509
## 11 GO:0012505                         endomembrane system       676
## 12 GO:0005789              endoplasmic reticulum membrane       308
## 13 GO:0031090                          organelle membrane      1225
## 14 GO:0043231    intracellular membrane-bounded organelle      6728
## 15 GO:0043227                  membrane-bounded organelle      6740
## 16 GO:0042175 nuclear outer membrane-endoplasmic retic...       314
## 17 GO:0005856                                cytoskeleton      1319
## 18 GO:0005576                        extracellular region      1170
## 19 GO:0044421                   extracellular region part       515
## 20 GO:0032991                      macromolecular complex      8767
## 21 GO:0005700                         polytene chromosome        94
## 22 GO:0000502                          proteasome complex       205
## 23 GO:0005783                       endoplasmic reticulum       650
## 24 GO:0045298                             tubulin complex       241
## 25 GO:0030529                   ribonucleoprotein complex      2123
## 26 GO:0005838              proteasome regulatory particle        53
## 27 GO:0045202                                     synapse       318
## 28 GO:0030880                      RNA polymerase complex       130
## 29 GO:0044432                  endoplasmic reticulum part       366
## 30 GO:0031301              integral to organelle membrane       111
## 31 GO:0022624                proteasome accessory complex        61
## 32 GO:0000123           histone acetyltransferase complex       146
## 33 GO:0070603            SWI/SNF superfamily-type complex        53
## 34 GO:0005819                                     spindle       149
## 35 GO:0030964                  NADH dehydrogenase complex        51
## 36 GO:0000313                         organellar ribosome        60
## 37 GO:0045271                 respiratory chain complex I        51
## 38 GO:0005761                      mitochondrial ribosome        57
## 39 GO:0005746             mitochondrial respiratory chain       107
## 40 GO:0030054                               cell junction       310
##    Significant Expected parentchild
## 1         7073  6187.11     < 1e-30
## 2         7067  6181.84     < 1e-30
## 3         3672  3165.95     < 1e-30
## 4         9905  9399.79     < 1e-30
## 5         9585  9102.29     < 1e-30
## 6         3597  3095.53     < 1e-30
## 7          269   182.32     3.1e-27
## 8         1055   897.13     5.2e-22
## 9         4403  3943.95     6.1e-21
## 10         431   335.03     3.2e-20
## 11         545   444.95     9.4e-17
## 12         258   202.73     7.4e-16
## 13         953   806.30     3.6e-15
## 14        5199  4428.39     4.6e-13
## 15        5206  4436.29     1.1e-12
## 16         261   206.68     2.4e-12
## 17        1054   868.17     2.8e-12
## 18         869   770.10     9.8e-11
## 19         403   338.97     3.1e-10
## 20        5977  5770.46     1.4e-09
## 21          89    61.87     1.8e-08
## 22         173   134.93     1.0e-07
## 23         538   427.83     1.1e-07
## 24         198   158.63     1.2e-07
## 25        1548  1397.36     4.8e-07
## 26          51    34.88     5.5e-07
## 27         249   209.31     6.7e-07
## 28         113    85.57     7.5e-07
## 29         308   240.90     1.1e-06
## 30          91    73.06     1.7e-06
## 31          57    40.15     2.3e-06
## 32         123    96.10     2.3e-06
## 33          50    34.88     3.0e-06
## 34         125    98.07     5.5e-06
## 35          47    33.57     5.7e-06
## 36          58    39.49     6.2e-06
## 37          47    33.57     7.0e-06
## 38          55    37.52     1.0e-05
## 39          92    70.43     1.3e-05
## 40         238   204.04     1.6e-05
```

```r

# graph significant nodes

pdf(paste(resultsdir, "Ap.CC_topGO_sig_nodes.pdf", sep = ""))
showSigOfNodes(Ap.CC.GOdata, score(Ap.CC.resultParentChild), firstSigNodes = 10, 
    useInfo = "all")
```

```
## $dag
## A graphNEL graph with directed edges
## Number of Nodes = 21 
## Number of Edges = 31 
## 
## $complete.dag
## [1] "A graph with 21 nodes."
```

```r
dev.off()
```

```
## pdf 
##   2
```


GO analysis for molecular function


```r
# create topGOdata object
Ap.MF.GOdata <- new("topGOdata", description = "MF gene set analysis", ontology = "MF", 
    allGenes = Ap.geneList, geneSel = selectFDR, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 3151 GO terms found. )
## 
## Build GO DAG topology ..........	( 3739 GO terms and 4617 relations. )
## 
## Annotating nodes ...............	( 34916 genes annotated to the GO terms. )
```

```r

Ap.MF.GOdata
```

```
## 
## ------------------------- topGOdata object -------------------------
## 
##  Description:
##    -  MF gene set analysis 
## 
##  Ontology:
##    -  MF 
## 
##  105536 available genes (all genes from the array):
##    - symbol:  0|*|Contig6267 100000|*|comp2663136_c0_seq1 100001|*|comp3439067_c0_seq1 100002|*|comp2050457_c0_seq1 100003|*|comp2029723_c0_seq1  ...
```

```
## Error: invalid 'digits' argument
```

```r

# perform enrichment analysis using multiple methods
Ap.MF.resultParentChild <- runTest(Ap.MF.GOdata, statistic = "fisher", algorithm = "parentchild")
```

```
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1391 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 15:	2 nodes to be scored.
## 
## 	 Level 14:	14 nodes to be scored.
## 
## 	 Level 13:	8 nodes to be scored.
## 
## 	 Level 12:	14 nodes to be scored.
## 
## 	 Level 11:	21 nodes to be scored.
## 
## 	 Level 10:	33 nodes to be scored.
## 
## 	 Level 9:	63 nodes to be scored.
## 
## 	 Level 8:	113 nodes to be scored.
## 
## 	 Level 7:	229 nodes to be scored.
## 
## 	 Level 6:	398 nodes to be scored.
## 
## 	 Level 5:	269 nodes to be scored.
## 
## 	 Level 4:	146 nodes to be scored.
## 
## 	 Level 3:	66 nodes to be scored.
## 
## 	 Level 2:	14 nodes to be scored.
```

```r
Ap.MF.resultParentChild
```

```
## 
## Description: MF gene set analysis 
## Ontology: MF 
## 'parentchild' algorithm with the 'fisher : joinFun = union' test
## 1391 GO terms scored: 112 terms with p < 0.01
## Annotation data:
##     Annotated genes: 34916 
##     Significant genes: 22273 
##     Min. no. of genes annotated to a GO: 10 
##     Nontrivial nodes: 1391
```

```r

Ap.MF.ResTable <- GenTable(Ap.MF.GOdata, parentchild = Ap.MF.resultParentChild, 
    topNodes = 20)
Ap.MF.ResTable
```

```
##         GO.ID                                        Term Annotated
## 1  GO:0005515                             protein binding      3185
## 2  GO:0046914                transition metal ion binding      4166
## 3  GO:0008233                          peptidase activity      2231
## 4  GO:0016705 oxidoreductase activity, acting on paire...       636
## 5  GO:0004888 transmembrane signaling receptor activit...      1034
## 6  GO:0003729                                mRNA binding       174
## 7  GO:0005198                structural molecule activity      1414
## 8  GO:0043169                              cation binding      8093
## 9  GO:0005488                                     binding     22779
## 10 GO:0016881             acid-amino acid ligase activity       506
## 11 GO:0004497                      monooxygenase activity       462
## 12 GO:0003676                        nucleic acid binding      8699
## 13 GO:0004674    protein serine/threonine kinase activity       759
## 14 GO:0003964        RNA-directed DNA polymerase activity       761
## 15 GO:0004767    sphingomyelin phosphodiesterase activity        43
## 16 GO:0004402          histone acetyltransferase activity       114
## 17 GO:0019787 small conjugating protein ligase activit...       333
## 18 GO:0008270                            zinc ion binding      3017
## 19 GO:0000975               regulatory region DNA binding       131
## 20 GO:0042623                    ATPase activity, coupled      1528
##    Significant Expected parentchild
## 1         2437  2031.72     < 1e-30
## 2         2929  2657.50     2.0e-19
## 3         1580  1423.16     4.7e-19
## 4          467   405.71     1.2e-12
## 5          730   659.59     1.3e-11
## 6          155   111.00     1.3e-10
## 7         1012   901.99     1.6e-10
## 8         5318  5162.54     1.0e-09
## 9        14774 14530.78     7.3e-09
## 10         369   322.78     9.4e-09
## 11         336   294.71     1.9e-08
## 12        5694  5549.11     3.2e-08
## 13         546   484.17     4.9e-08
## 14         524   485.44     6.9e-08
## 15          42    27.43     8.5e-08
## 16          97    72.72     1.2e-07
## 17         268   212.42     1.5e-07
## 18        2189  1924.55     2.1e-07
## 19         106    83.57     4.0e-07
## 20         931   974.71     7.1e-07
```

```r

# graph significant nodes

pdf(paste(resultsdir, "Ap.MF_topGO_sig_nodes.pdf", sep = ""))
showSigOfNodes(Ap.MF.GOdata, score(Ap.MF.resultParentChild), firstSigNodes = 10, 
    useInfo = "all")
```

```
## $dag
## A graphNEL graph with directed edges
## Number of Nodes = 26 
## Number of Edges = 27 
## 
## $complete.dag
## [1] "A graph with 26 nodes."
```

```r
dev.off()
```

```
## pdf 
##   2
```


Note that among responsive transcripts, there are 26 transcripts with GO term "response to stress":

1038|*|comp150483_c5_seq3, 11281|*|comp146961_c0_seq1, 12704|*|comp144775_c1_seq1, 14|*|comp150262_c0_seq1, 1504|*|Contig2729, 15115|*|comp132715_c0_seq1, 17710|*|comp150271_c3_seq3, 19475|*|Contig1438, 2087|*|comp150483_c5_seq1, 21384|*|comp149042_c0_seq3, 21598|*|comp142101_c0_seq1, 23441|*|comp114823_c1_seq1, 2604|*|comp148324_c0_seq4, 2832|*|comp150483_c5_seq2, 32312|*|comp933733_c0_seq1, 37752|*|comp1620595_c0_seq1, 3995|*|comp145243_c0_seq1, 47691|*|comp1460938_c0_seq1, 51985|*|comp1012776_c0_seq1, 552|*|comp147487_c0_seq1, 58246|*|comp109744_c0_seq1, 6075|*|comp92770_c0_seq1, 75624|*|comp2836178_c0_seq1, 80544|*|comp132706_c0_seq1, 9316|*|comp147545_c4_seq2, 9372|*|Contig4757

and this term is included in the list of enriched GO terms:


```r
# significant GO terms
Ap.BP.signif.table <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, 
    topNodes = 415)
Ap.BP.signif <- Ap.BP.signif.table$GO.ID
length(Ap.BP.signif)
```

```
## [1] 415
```

```r

Ap.BP.signif.term <- Ap.BP.signif.table$Term
Ap.BP.signif.term[grep("stress", Ap.BP.signif.term)]
```

```
## [1] "response to stress"               "regulation of response to stress"
```

```r
Ap.BP.signif.term[grep("immune", Ap.BP.signif.term)]
```

```
## [1] "immune system process"               
## [2] "immune response"                     
## [3] "regulation of immune system process" 
## [4] "regulation of innate immune response"
## [5] "regulation of immune response"       
## [6] "innate immune response"
```


Export data for interactive shiny app. 


```r
# scale expression values
Ap.dt[, `:=`(exp.scaled, scale(TPM)), by = Transcript]
```

```
##                           Transcript Length    TPM   RPKM   KPKM
##      1: 100015|*|comp3543055_c0_seq1    208 0.7914 0.9487 0.9487
##      2: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      3: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      4: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##      5: 100015|*|comp3543055_c0_seq1    208 0.0000 0.0000 0.0000
##     ---                                                         
## 121876:       9|*|comp147140_c0_seq1   9030 0.5567 1.1091 1.1091
## 121877:       9|*|comp147140_c0_seq1   9030 0.5736 0.8560 0.8560
## 121878:       9|*|comp147140_c0_seq1   9030 0.2928 0.6473 0.6473
## 121879:       9|*|comp147140_c0_seq1   9030 0.4126 0.6169 0.6169
## 121880:       9|*|comp147140_c0_seq1   9030 0.3079 0.6398 0.6398
##         EstimatedNumReads sample  val colony sequence.length
##      1:           0.09029  A22-0  0.0    A22             208
##      2:           0.00000   Ar-0  0.0     Ar             208
##      3:           0.00000  A22-3  3.5    A22             208
##      4:           0.00000   Ar-3  3.5     Ar             208
##      5:           0.00000 A22-10 10.5    A22             208
##     ---                                                     
## 121876:           2.37111  Ar-31 31.5     Ar            9030
## 121877:           3.18572 A22-35 35.0    A22            9030
## 121878:           1.03722  Ar-35 35.0     Ar            9030
## 121879:           2.86903 A22-38 38.5    A22            9030
## 121880:           1.39419  Ar-38 38.5     Ar            9030
##                                                         best.hit.to.nr
##      1: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      2: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      3: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      4: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##      5: gi|121608385|ref|YP_996192.1| transposase, IS4 family protein 
##     ---                                                               
## 121876:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121877:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121878:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121879:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
## 121880:   gi|322801453|gb|EFZ22114.1| hypothetical protein SINV_07423 
##         hit.length  E.value   Bit.score
##      1:         25 7.22e-06   57.514374
##      2:         25 7.22e-06   57.514374
##      3:         25 7.22e-06   57.514374
##      4:         25 7.22e-06   57.514374
##      5:         25 7.22e-06   57.514374
##     ---                                
## 121876:       2101      0.0 4740.368306
## 121877:       2101      0.0 4740.368306
## 121878:       2101      0.0 4740.368306
## 121879:       2101      0.0 4740.368306
## 121880:       2101      0.0 4740.368306
##                                                                                                                                                                                                          GO.Biological.Process
##      1:                                                                                                                                                                                                                      -
##      2:                                                                                                                                                                                                                      -
##      3:                                                                                                                                                                                                                      -
##      4:                                                                                                                                                                                                                      -
##      5:                                                                                                                                                                                                                      -
##     ---                                                                                                                                                                                                                       
## 121876: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121877: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121878: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121879: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
## 121880: GO:0032418 lysosome localization | GO:0006468 protein phosphorylation | GO:0007264 small GTPase mediated signal transduction | GO:0050808 synapse organization | GO:0009069 serine family amino acid metabolic process
##                                                     GO.Cellular.Component
##      1:                                                                 -
##      2:                                                                 -
##      3:                                                                 -
##      4:                                                                 -
##      5:                                                                 -
##     ---                                                                  
## 121876: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121877: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121878: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121879: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
## 121880: GO:0005765 lysosomal membrane | GO:0031902 late endosome membrane
##                                                                                                                         GO.Molecular.Function
##      1:                                                                                                                                     -
##      2:                                                                                                                                     -
##      3:                                                                                                                                     -
##      4:                                                                                                                                     -
##      5:                                                                                                                                     -
##     ---                                                                                                                                      
## 121876: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121877: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121878: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121879: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
## 121880: GO:0017137 Rab GTPase binding | GO:0004674 protein serine/threonine kinase activity | GO:0005524 ATP binding | GO:0005525 GTP binding
##         Enzyme
##      1:      -
##      2:      -
##      3:      -
##      4:      -
##      5:      -
##     ---       
## 121876:      -
## 121877:      -
## 121878:      -
## 121879:      -
## 121880:      -
##                                                                                                                                                     Domain
##      1:                                                                                                                                                  -
##      2:                                                                                                                                                  -
##      3:                                                                                                                                                  -
##      4:                                                                                                                                                  -
##      5:                                                                                                                                                  -
##     ---                                                                                                                                                   
## 121876: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121877: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121878: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121879: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
## 121880: pfam07714 Pkinase_Tyr | pfam00069 Pkinase | pfam08477 Miro | pfam12796 Ank_2 | pfam13855 LRR_8 | pfam13637 Ank_4 | pfam12799 LRR_4 | pfam00023 Ank
##         annotation.type      pval coef.colony coef.val coef.I(val^2)
##      1:                 0.0014775   0.0448835       NA      0.002491
##      2:                 0.0014775   0.0448835       NA      0.002491
##      3:                 0.0014775   0.0448835       NA      0.002491
##      4:                 0.0014775   0.0448835       NA      0.002491
##      5:                 0.0014775   0.0448835       NA      0.002491
##     ---                                                             
## 121876:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121877:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121878:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121879:     GO & Domain 0.0006044   0.0001822 0.008766            NA
## 121880:     GO & Domain 0.0006044   0.0001822 0.008766            NA
##         coef.colony:val coef.colony:I(val^2)    qval max.val min.val
##      1:              NA             0.002523 0.02497       0    21.0
##      2:              NA             0.002523 0.02497       0    21.0
##      3:              NA             0.002523 0.02497       0    21.0
##      4:              NA             0.002523 0.02497       0    21.0
##      5:              NA             0.002523 0.02497       0    21.0
##     ---                                                             
## 121876:         0.02394                   NA 0.01335       0    38.5
## 121877:         0.02394                   NA 0.01335       0    38.5
## 121878:         0.02394                   NA 0.01335       0    38.5
## 121879:         0.02394                   NA 0.01335       0    38.5
## 121880:         0.02394                   NA 0.01335       0    38.5
##         exp_type exp.scaled
##      1:      low     4.0965
##      2:      low    -0.3295
##      3:      low    -0.3295
##      4:      low    -0.3295
##      5:      low    -0.3295
##     ---                    
## 121876:      low    -0.3185
## 121877:      low    -0.2689
## 121878:      low    -1.0902
## 121879:      low    -0.7400
## 121880:      low    -1.0461
```

```r
str(Ap.dt)
```

```
## Classes 'data.table' and 'data.frame':	121880 obs. of  31 variables:
##  $ Transcript           : chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ Length               : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM                  : num  0.791 0 0 0 0 ...
##  $ RPKM                 : num  0.949 0 0 0 0 ...
##  $ KPKM                 : num  0.949 0 0 0 0 ...
##  $ EstimatedNumReads    : num  0.0903 0 0 0 0 ...
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
##  $ pval                 : num  0.00148 0.00148 0.00148 0.00148 0.00148 ...
##  $ coef.colony          : num  0.0449 0.0449 0.0449 0.0449 0.0449 ...
##  $ coef.val             : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.I(val^2)        : num  0.00249 0.00249 0.00249 0.00249 0.00249 ...
##  $ coef.colony:val      : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.colony:I(val^2) : num  0.00252 0.00252 0.00252 0.00252 0.00252 ...
##  $ qval                 : num  0.025 0.025 0.025 0.025 0.025 ...
##  $ max.val              : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ min.val              : num  21 21 21 21 21 21 21 21 21 21 ...
##  $ exp_type             : chr  "low" "low" "low" "low" ...
##  $ exp.scaled           : num  4.096 -0.329 -0.329 -0.329 -0.329 ...
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
## Classes 'data.table' and 'data.frame':	82016 obs. of  31 variables:
##  $ Transcript           : chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ Length               : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM                  : num  0.791 0 0 0 0 ...
##  $ RPKM                 : num  0.949 0 0 0 0 ...
##  $ KPKM                 : num  0.949 0 0 0 0 ...
##  $ EstimatedNumReads    : num  0.0903 0 0 0 0 ...
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
##  $ pval                 : num  0.00148 0.00148 0.00148 0.00148 0.00148 ...
##  $ coef.colony          : num  0.0449 0.0449 0.0449 0.0449 0.0449 ...
##  $ coef.val             : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.I(val^2)        : num  0.00249 0.00249 0.00249 0.00249 0.00249 ...
##  $ coef.colony:val      : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ coef.colony:I(val^2) : num  0.00252 0.00252 0.00252 0.00252 0.00252 ...
##  $ qval                 : num  0.025 0.025 0.025 0.025 0.025 ...
##  $ max.val              : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ min.val              : num  21 21 21 21 21 21 21 21 21 21 ...
##  $ exp_type             : chr  "low" "low" "low" "low" ...
##  $ exp.scaled           : num  4.096 -0.329 -0.329 -0.329 -0.329 ...
##  - attr(*, "sorted")= chr "Transcript"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
write.csv(Ap.dt.interaction, file = paste(resultsdir, "Ap.dt.interaction.csv", 
    sep = ""), quote = TRUE, row.names = FALSE)
```


# Visualize responsive transcripts

Make plots for all significant transcripts

![plot of chunk plot_responsive](figure/plot_responsive1.png) 

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk plot_responsive](figure/plot_responsive2.png) ![plot of chunk plot_responsive](figure/plot_responsive3.png) 

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk plot_responsive](figure/plot_responsive4.png) 


Make plots for transcripts with significant temperature by colony interaction

![plot of chunk plot_interaction_responsive](figure/plot_interaction_responsive1.png) 

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk plot_interaction_responsive](figure/plot_interaction_responsive2.png) ![plot of chunk plot_interaction_responsive](figure/plot_interaction_responsive3.png) 

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

![plot of chunk plot_interaction_responsive](figure/plot_interaction_responsive4.png) 





```
## pdf 
##   2
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## pdf 
##   2
```

```
## pdf 
##   2
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## pdf 
##   2
```




```
## pdf 
##   2
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## pdf 
##   2
```

```
## pdf 
##   2
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## pdf 
##   2
```



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
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] Rgraphviz_2.6.0      topGO_2.14.0         SparseM_1.03        
##  [4] GO.db_2.10.1         RSQLite_0.11.4       DBI_0.2-7           
##  [7] AnnotationDbi_1.24.0 Biobase_2.22.0       BiocGenerics_0.8.0  
## [10] graph_1.40.1         plyr_1.8             RCurl_1.95-4.1      
## [13] bitops_1.0-6         data.table_1.8.10    stringr_0.6.2       
## [16] pander_0.3.8         knitcitations_0.5-0  bibtex_0.3-6        
## [19] ggplot2_0.9.3.1      R.utils_1.28.4       R.oo_1.17.0         
## [22] R.methodsS3_1.6.1    knitr_1.5           
## 
## loaded via a namespace (and not attached):
##  [1] codetools_0.2-8    colorspace_1.2-4   dichromat_2.0-0   
##  [4] digest_0.6.4       evaluate_0.5.1     formatR_0.10      
##  [7] gtable_0.1.2       httr_0.2           IRanges_1.20.6    
## [10] labeling_0.2       lattice_0.20-24    MASS_7.3-29       
## [13] munsell_0.4.2      proto_0.3-10       RColorBrewer_1.0-5
## [16] reshape2_1.2.2     scales_0.2.3       stats4_3.0.2      
## [19] tools_3.0.2        XML_3.98-1.1       xtable_1.7-1
```

