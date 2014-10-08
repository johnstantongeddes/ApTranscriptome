Thermal reactionome of the common ant species *Aphaenogaster picea* and *A. carolinensis*
========================================================================================
   
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**July 17, 2014**

**Technical Report No. 3**

**Department of Biology**

**University of Vermont**



## Summary ##
  
In this technical report, which accompanies the manuscript **Patterns of stress tolerance versus resistance in genome-wide expression data of parapatric ant species** (Stanton-Geddes et al., in press), we:

1. describe the *de novo* assembly of the transcriptome for two ant species within the *Aphaenogaster rudis-picea-fulva* species complex (Lubertazzi, 2012)
2. identify thermally-responsive genes
3. evaluate differences in reaction norms of expression for each thermally-responsive gene between the two species
4. perform gene set enrichment analysis of thermally-responsive genes for the two species

This script is completely reproducible assuming that R, `knitr` and the other required libraries (listed within the source document) are installed on a standard linux system using the following:
    
    Rscript -e "library(knitr); knit('ApTranscriptome_TR.Rmd')"


## Data ##

The raw Illumina fastq files are available from the NCBI short read archive [link tbd]. The assembled transcriptome, annotation and expression values are downloaded rather than re-run due to the computational demands, but the exact commands for each of these steps are documented below.

## Sample description ##

Two ant colonies were used for the transcriptome sequencing. The first, designated *A22*, was collected at Molly Bog, Vermont in August 2012 by Nick Gotelli and Andrew Nguyen. The second colony, designated *Ar*, was collected by Lauren Nichols in Raleigh, North Carolina. These colonies were maintained in the lab for 6 months prior to sample collection. Bernice Bacon DeMarco (Michigan State University) identified colony *A22* as *A. picea* and *Ar* as *A. carolinensis*. For historical reasons, I refer to these species as colonies at times throughout this technical report.

For each colony, three ants were exposed to one of 12 temperature treatments, every 3.5C ranging from 0C to 38.5C, for one hour in glass tubes in a water bath. The ants were flash frozen and stored at -80C until RNA was extracted using a two step extraction; [RNAzol RT](http://www.mrcgene.com/rnazol.htm) (Molecular Research Center, Inc) followed by an [RNeasy Micro](http://www.qiagen.com/products/catalog/sample-technologies/rna-sample-technologies/total-rna/rneasy-micro-kit) column (Qiagen). Samples from each colony were pooled and sequenced in separate lanes on a 100bp paired-end run of an Illumina HiSeq at the University of Minnesota Genomics Center, yielding 20e6 and 16e6 reads for the A22 and Ar samples, respectively.



## Transcriptome assembly ##

The Illumina reads were filtered using the program [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (Lohse, Bolger, Nagel, Fernie, Lunn, Stitt, and Usadel, 2012) to remove Ilumina adapter sequences, trim bases with PHRED quality scores less than 15 in a 4 bp sliding window and remove final trimmed sequences with length less than 36 bp. The code used was 

~~~
java -jar trimmomatic-0.30.jar PE -threads 40 -phred33 -trimlog trimmomatic.log sample.R1.fastq sample.R2.fastq sample.R1.trimmed.paired.fastq sample.R1.trimmed.unpaired.fastq sample.R2.trimmed.paired.fastq sample.R2.trimmed.unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
~~~

where "sample" was replaced by each sample. This filtering yielded 339,845,787 properly paired reads and 16,860,885 unpaired reads. 

Properly paired and unpaired reads passing the Trimmomatic filter were combined and used in de novo transcriptome assembly using the program [Trinity](http://trinityrnaseq.sourceforge.net/) (Grabherr, Haas, Yassour, Levin, Thompson, Amit, Adiconis, Fan, Raychowdhury, Zeng, Chen, Mauceli, Hacohen, Gnirke, Rhind, di
Palma, Birren, Nusbaum, Lindblad-Toh, Friedman, and Regev, 2011). Unpaired reads resulting from Trimmomatic were concatenated to the set of left (R1) reads per Trinity usage documentation. Assembly on a single 40 CPU machine with 1TB of memory required 20 hours of compute time.

~~~
Trinity.pl --seqType=fq --JM=90G —left=all.R1.trimmed.fastq —right=all.R2.trimmed.fastq --output=trinity --CPU=40 --inchworm_cpu=40 --bflyCPU=5
~~~

This assembly contained 100,381 unique components (roughly genes) in 126,172 total transcripts (Table 1). 

As we were assembling two divergent colonies into a single transcriptome, we suspected that this assembly would be susceptible to known problems of errors during assembly (e.g. chimeric transcripts that are fusions of two transcripts) and redundancy (Yang and Smith, 2013). To account for this, we performed two post-assembly processing steps.

First, we ran the program [cap3](http://seq.cs.iastate.edu/) (Huang, 1999) setting the maximum gap length and band expansion size to 50 `-f 50 -a 50`, no end clipping as the reads were already filtered `k 0`, requiring 90% identity for assembly, and a minimum overlap length of 100 bp `-o 100`. The percent identity threshold of 90% was chosen to liberally collapse orthologous contigs from the two colonies that may have been assembled separately. 

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

Table: Statistics for Trinity and cap3+uclust reduced transcriptome assemblies (continued below)

 

|    &nbsp;     |  Mean contig size  |  N50 contig  |  N50 Length  |
|:-------------:|:------------------:|:------------:|:------------:|
|  **trinity**  |        795         |    16,201    |    1,631     |
|  **reduced**  |        593         |    15,491    |     895      |

To remove contigs that are likely contaminants from bacterial, archael, virus or human sources, we used the program [DeconSeq](http://deconseq.sourceforge.net/) r citep("10.1371/journal.pone.0017288"). We downloaded the bacteria, virus, archae and human [databases of contaminants](ftp://edwards.sdsu.edu:7009/deconseq/db), modified the `DeconSeqConfig.pm` file as described [here](http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/install/deconseq.htm) to point to the databases, and ran DeconSeq specifiying 95% identity over 50% the length of contig

    deconseq.pl -c 50 -i 95 -f Trinity_cap3_uclust.fasta -d Trinity_cap3_uclust -dbs hsref,bast,vir,arch
	
This resulted in removing 5,675 contigs as contaminants, leaving 99,861 "clean" contigs. We spot-checked the contaminants by BLAST and confirmed that they matched bacteria, human or viral sources by greater than 95%. For expression quantification, we use the full assembly to ensure that "contaminant" reads are assigned to the contaminants. After quantification, these transcripts will then be removed from further analyses.

Running Trinity and subsequent programs is time and memory-intensive so the final assembly can be downloaded and used for all further analyses. In addition, the compressed archive contains the "clean" and "contaminant" sequences after filtering with DeconSeq.

~~~
# download filtered Trinity assembly, uncompress and move
wget http://johnstantongeddes.org/assets/files/ApTranscriptome/Aphaenogaster_transcriptome.tar
# check md5sum
md5sum Aphaenogaster_transcriptome.tar
# fd6dbb0b3e88e1c797c9e74611b245b2

# move and extract
mkdir -p results/
mkdir -p results/trinity-full/
mv Aphaenogaster_transcriptome.tar results/trinity-full/.
tar -xvf Aphaenogaster_transcriptome.tar
~~~

To examine the species distribution of BLAST hits in the transcriptome assembly, I used the program [Krona](http://sourceforge.net/p/krona/home/krona/) r citep("doi:10.1186/1471-2105-12-385"). I ... 

~~~
KRONA code
~~~

The interactive visualization is available [here]().


## Transcriptome annotation ##

Annotation was performed by uploading the reduced assembly "Trinity_cap3_uclust.fasta" to the web-based annotation program [FastAnnotator](http://fastannotator.cgu.edu.tw/index.php).

Results are available as job ID [13894410176993](http://fastannotator.cgu.edu.tw/job.php?jobid=13894410176993#page=basicinfo).

This annotation file can be read directly to R:


```r
# URL for annotation file
annotation.URL <- getURL("http://johnstantongeddes.org/assets/files/ApTranscriptome/ApTranscriptome_AnnotationTable_20140113.txt")
# load 
annotation.file <- read.csv(textConnection(annotation.URL), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(annotation.file)
```

```
## 'data.frame':	105536 obs. of  12 variables:
##  $ Sequence.Name        : chr  "0|*|Contig6267" "1|*|comp150820_c2_seq6" "2|*|Contig6262" "3|*|comp149397_c1_seq2" ...
##  $ sequence.length      : int  9990 9944 9711 9639 9558 9436 9410 9396 9103 9030 ...
##  $ best.hit.to.nr       : chr  "gi|110756860|ref|XP_392375.3| PREDICTED: hypothetical protein LOC408844 " "gi|307181425|gb|EFN69020.1| FH2 domain-containing protein 1 " "gi|110756860|ref|XP_392375.3| PREDICTED: hypothetical protein LOC408844 " "gi|332024185|gb|EGI64399.1| AT-rich interactive domain-containing protein 5B " ...
##  $ hit.length           : chr  "598" "1777" "511" "1077" ...
##  $ E.value              : chr  "2.1e-265" "0.0" "1.29e-246" "0.0" ...
##  $ Bit.score            : chr  "904.618736" "3455.802741" "842.252472" "2272.638439" ...
##  $ GO.Biological.Process: chr  "GO:0035335 peptidyl-tyrosine dephosphorylation | GO:0000188 inactivation of MAPK activity | GO:0006570 tyrosine metabolic proce"| __truncated__ "GO:0030036 actin cytoskeleton organization | GO:0015074 DNA integration" "GO:0035335 peptidyl-tyrosine dephosphorylation | GO:0000188 inactivation of MAPK activity | GO:0006570 tyrosine metabolic proce"| __truncated__ "GO:0006508 proteolysis" ...
##  $ GO.Cellular.Component: chr  "-" "-" "-" "GO:0005634 nucleus" ...
##  $ GO.Molecular.Function: chr  "GO:0017017 MAP kinase tyrosine/serine/threonine phosphatase activity | GO:0004725 protein tyrosine phosphatase activity | GO:00"| __truncated__ "GO:0003676 nucleic acid binding | GO:0003779 actin binding | GO:0008270 zinc ion binding | GO:0017048 Rho GTPase binding" "GO:0017017 MAP kinase tyrosine/serine/threonine phosphatase activity | GO:0004725 protein tyrosine phosphatase activity | GO:00"| __truncated__ "GO:0003677 DNA binding | GO:0004252 serine-type endopeptidase activity" ...
##  $ Enzyme               : chr  "3.1.3.16  | 3.1.3.48 " "-" "3.1.3.16  | 3.1.3.48 " "-" ...
##  $ Domain               : chr  "pfam00782 DSPc | pfam00581 Rhodanese" "pfam02181 FH2 | pfam00067 p450 | pfam06367 Drf_FH3 | pfam12795 MscS_porin | pfam01749 IBB | pfam07926 TPR_MLP1_2" "pfam00782 DSPc | pfam00581 Rhodanese" "pfam01388 ARID" ...
##  $ annotation.type      : chr  "GO & Enzyme & Domain" "GO & Domain" "GO & Enzyme & Domain" "GO & Domain" ...
```

```r
# Convert to data.table
annotation.table <- data.table(annotation.file)
str(annotation.table)
```

```
## Classes 'data.table' and 'data.frame':	105536 obs. of  12 variables:
##  $ Sequence.Name        : chr  "0|*|Contig6267" "1|*|comp150820_c2_seq6" "2|*|Contig6262" "3|*|comp149397_c1_seq2" ...
##  $ sequence.length      : int  9990 9944 9711 9639 9558 9436 9410 9396 9103 9030 ...
##  $ best.hit.to.nr       : chr  "gi|110756860|ref|XP_392375.3| PREDICTED: hypothetical protein LOC408844 " "gi|307181425|gb|EFN69020.1| FH2 domain-containing protein 1 " "gi|110756860|ref|XP_392375.3| PREDICTED: hypothetical protein LOC408844 " "gi|332024185|gb|EGI64399.1| AT-rich interactive domain-containing protein 5B " ...
##  $ hit.length           : chr  "598" "1777" "511" "1077" ...
##  $ E.value              : chr  "2.1e-265" "0.0" "1.29e-246" "0.0" ...
##  $ Bit.score            : chr  "904.618736" "3455.802741" "842.252472" "2272.638439" ...
##  $ GO.Biological.Process: chr  "GO:0035335 peptidyl-tyrosine dephosphorylation | GO:0000188 inactivation of MAPK activity | GO:0006570 tyrosine metabolic proce"| __truncated__ "GO:0030036 actin cytoskeleton organization | GO:0015074 DNA integration" "GO:0035335 peptidyl-tyrosine dephosphorylation | GO:0000188 inactivation of MAPK activity | GO:0006570 tyrosine metabolic proce"| __truncated__ "GO:0006508 proteolysis" ...
##  $ GO.Cellular.Component: chr  "-" "-" "-" "GO:0005634 nucleus" ...
##  $ GO.Molecular.Function: chr  "GO:0017017 MAP kinase tyrosine/serine/threonine phosphatase activity | GO:0004725 protein tyrosine phosphatase activity | GO:00"| __truncated__ "GO:0003676 nucleic acid binding | GO:0003779 actin binding | GO:0008270 zinc ion binding | GO:0017048 Rho GTPase binding" "GO:0017017 MAP kinase tyrosine/serine/threonine phosphatase activity | GO:0004725 protein tyrosine phosphatase activity | GO:00"| __truncated__ "GO:0003677 DNA binding | GO:0004252 serine-type endopeptidase activity" ...
##  $ Enzyme               : chr  "3.1.3.16  | 3.1.3.48 " "-" "3.1.3.16  | 3.1.3.48 " "-" ...
##  $ Domain               : chr  "pfam00782 DSPc | pfam00581 Rhodanese" "pfam02181 FH2 | pfam00067 p450 | pfam06367 Drf_FH3 | pfam12795 MscS_porin | pfam01749 IBB | pfam07926 TPR_MLP1_2" "pfam00782 DSPc | pfam00581 Rhodanese" "pfam01388 ARID" ...
##  $ annotation.type      : chr  "GO & Enzyme & Domain" "GO & Domain" "GO & Enzyme & Domain" "GO & Domain" ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

Note that I used the reduced assembly, prior to cleaning out contaminants. As I use the "clean" assembly for read mapping and identification of responsive genes, the contaminants are dis-regarded in downstream analyses. 


## Identification of thermally-responsive genes ##

### Quantify gene expression ###

I quantified gene expression using [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html). To run this program, first make sure that PATHs to the software libraries are set up correctly as described on the sailfish website. 
                                                 
An index of the assembly is built with the command:                                                 

~~~
sailfish index -t results/trinity-full/Trinity_cap3_uclust.fasta -o results/trinity-full/sailfish-index-Trinity-cap3-uclust -k 20 -p 4
~~~

Once this is done, expression is quantified for the Trimmomatic filtered reads from each species-treatment sample separately. Note that for each sample, there are four filtered read files:

- paired.left.fastq
- paired.right.fastq
- unpaired.left.fastq
- unpaired.right.fastq



~~~
# make a directory for the expression values
mkdir -p results/trinity-full/sailfish-expression-Trinity-cap3-uclust
# change to "trinity-full" directory
cd results/trinity-full
# for each sample, run the following command
sailfish quant -i sailfish-index-Trinity-cap3-uclust -o sailfish-expression-Trinity-cap3-uclust/A22-0 -l "T=SE:S=U" -r A22-0_ATCACG.unpaired.left.fastq A22-0_ATCACG.unpaired.right.fastq A22-0_ATCACG.paired.left.fastq A22-0_ATCACG.paired.right.fastq -p 4
~~~
                                           
While it is possible to separately specify the paired-end and orphaned single-end reads in Sailfish v0.6.3, the results are exactly the same as if they are all entered as SE.	

These files are downloaded for convenience:

~~~
# download gene expression quantification
wget http://johnstantongeddes.org/assets/files/ApTranscriptome/sailfish_quant_20140916.tar
# check md5sum
md5sum sailfish_quant_20140916.tar
# 26102c7ef86cf30b5f8e923640378185

# move and extract
mkdir -p /results/trinity-full/sailfish-expression-Trinity-cap3-uclust
mv sailfish_quant_20140916.tar results/trinity-full/.
tar -xvf sailfish_quant_20140916.tar
~~~

For each sample, there is a directory containting a file *quant_bias_corrected.sf*. This file has the following columns, following a number of header lines:

1. Transcript ID
2. Transcript Length
3. Transcripts per Million (TPM): computed as described in (Li, Ruotti, Stewart, Thomson, and Dewey, 2009), and is meant as an estimate of the number of transcripts, per million observed transcripts, originating from each isoform.
4. Reads Per Kilobase per Million mapped reads (RPKM): classic measure of relative transcript abundance, and is an estimate of the number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

The TPM column for each sample was extracted and combined into a matrix for each species.

**Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the combined expression data set for the two species.** 




Note that expression levels at each temperature treatment are highly correlated between the two colonies.




|  Temperature  |  r   |
|:-------------:|:----:|
|       0       | 0.99 |
|      3.5      | 0.98 |
|     10.5      |  1   |
|      14       | 0.99 |
|     17.5      | 0.98 |
|      21       | 0.98 |
|     24.5      | 0.99 |
|      28       | 0.99 |
|     31.5      | 0.99 |
|      35       | 0.99 |
|     38.5      | 0.99 |

Table: Correlations between species for gene expression at temperature treatment


### Remove *contaminant* transcripts

With gene expression quantified using all reads, I now remove the transcripts identified as contaminants by DeconSeq.


```r
# extract transcript names from fasta file rather than loading whole file
system("grep '^>' results/trinity-full/Trinity_cap3_uclust_cont.fa | cut -f 1 -d ' ' | sed -r 's/^.{1}//' > results/trinity-full/cont.list")

# read file
cont.list <- read.table("results/trinity-full/cont.list")
str(cont.list)
```

```
## 'data.frame':	5675 obs. of  1 variable:
##  $ V1: chr  "225|*|Contig6496" "441|*|comp151639_c0_seq1" "2950|*|comp140219_c0_seq8" "3486|*|comp150480_c0_seq1" ...
```

```r
# remove from TPM.dt.sub
setkey(TPM.dt.sub, Transcript)
TPM.dt.sub <- TPM.dt.sub[!cont.list]
str(TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	2160092 obs. of  10 variables:
##  $ Transcript       : chr  "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" ...
##  $ Length           : int  9990 9990 9990 9990 9990 9990 9990 9990 9990 9990 ...
##  $ TPM              : num  0.079 0.0643 0.0357 0.093 0.0395 ...
##  $ RPKM             : num  0.0926 0.1078 0.0532 0.1418 0.0494 ...
##  $ KPKM             : num  0.0926 0.1078 0.0532 0.1418 0.0494 ...
##  $ EstimatedNumKmers: num  2974 1500 1373 2433 1713 ...
##  $ EstimatedNumReads: num  37.1 18.6 17.1 30.1 21.4 ...
##  $ sample           : chr  "A22-0" "Ar-0" "A22-3" "Ar-3" ...
##  $ val              : num  0 0 3.5 3.5 10.5 10.5 14 14 17.5 17.5 ...
##  $ colony           : Factor w/ 2 levels "A22","Ar": 1 2 1 2 1 2 1 2 1 2 ...
##  - attr(*, "sorted")= chr "Transcript"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
length(unique(TPM.dt.sub$Transcript))
```

```
## [1] 98186
```

```r
# remove from annotation.table
setkey(annotation.table, Sequence.Name)
annotation.table <- annotation.table[!cont.list]
str(annotation.table)
```

```
## Classes 'data.table' and 'data.frame':	99861 obs. of  12 variables:
##  $ Sequence.Name        : chr  "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
##  $ sequence.length      : int  9990 208 208 208 208 208 208 208 208 208 ...
##  $ best.hit.to.nr       : chr  "gi|110756860|ref|XP_392375.3| PREDICTED: hypothetical protein LOC408844 " "-" "-" "gi|189208225|ref|XP_001940446.1| 60S ribosomal protein L23 " ...
##  $ hit.length           : chr  "598" "-" "-" "69" ...
##  $ E.value              : chr  "2.1e-265" "-" "-" "5.11e-40" ...
##  $ Bit.score            : chr  "904.618736" "-" "-" "161.159029" ...
##  $ GO.Biological.Process: chr  "GO:0035335 peptidyl-tyrosine dephosphorylation | GO:0000188 inactivation of MAPK activity | GO:0006570 tyrosine metabolic proce"| __truncated__ "-" "-" "GO:0006508 proteolysis | GO:0006412 translation | GO:0042254 ribosome biogenesis" ...
##  $ GO.Cellular.Component: chr  "-" "-" "-" "GO:0005840 ribosome" ...
##  $ GO.Molecular.Function: chr  "GO:0017017 MAP kinase tyrosine/serine/threonine phosphatase activity | GO:0004725 protein tyrosine phosphatase activity | GO:00"| __truncated__ "-" "-" "GO:0008233 peptidase activity | GO:0003735 structural constituent of ribosome" ...
##  $ Enzyme               : chr  "3.1.3.16  | 3.1.3.48 " "-" "-" "-" ...
##  $ Domain               : chr  "pfam00782 DSPc | pfam00581 Rhodanese" "pfam03993 DUF349" "-" "pfam00238 Ribosomal_L14" ...
##  $ annotation.type      : chr  "GO & Enzyme & Domain" "Domain only" "" "GO & Domain" ...
##  - attr(*, "sorted")= chr "Sequence.Name"
##  - attr(*, ".internal.selfref")=<externalptr>
```

### Regression-model to identify thermally-responsive genes

To identify transcripts (roughly equivalent to genes) that show thermal responsiveness, I fit the following linear model to each transcript:

$$ log(TPM + 1) = \beta_0 + \beta_1(species) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(species * temp) + \beta_5(species * temp^2) + \epsilon $$

where TPM is transcripts per million. 


(1) Identify transcripts with overall significant model fit. Adjust *P* values for multiple testing using FDR and retain transcripts with FDR < 0.05. Use log-transformed response to account for outliers.



```r
# define model for RxN function
model <-  "log(TPM+1) ~ colony + val + I(val^2) + colony:val + colony:I(val^2)"

# calculate overall P value and R^2 for each transcript
RxNpval <- ddply(TPM.dt.sub, .(Transcript), .inform="TRUE", modpFunc)
```

Of the 98186 transcripts, 22089 have models with P < 0.05.

Many of these are likely false positives, so I adjust P-values using false discovery rate (FDR). Only those transcripts with less than 5% FDR are retained as significant. 


```r
RxNpval$padj <- p.adjust(RxNpval$pval, method = "fdr")
# Plot FDR values against initial pvalues
par(mfrow = c(2,1))
hist(RxNpval$pval)
hist(RxNpval$padj)
```

![plot of chunk padjust](figure/padjust.png) 

```r
# subset to significant transcripts
signif.transcripts <- RxNpval[which(RxNpval$padj < 0.05), ]

# extract significant transcripts
sig.TPM.dt.sub <- TPM.dt.sub[signif.transcripts$Transcript]
```

At the 5% FDR significance threshold, there are 10525 transcripts with an overall significant model.


(2) Fit linear model to overall significant transcripts; perform stepAIC to retain only significant terms, and save `lm` output to list


```r
# perform model selection for responsive transcripts
# need to use `try` to avoid stopping on error for AIC at Infinity
RxNlmAIC <- try(dlply(sig.TPM.dt.sub, .(Transcript), lmFunc))
```


### Grouping of thermally-responsive transcripts

The set of transcripts with significant expression patterns include those with expression that differs by species, temperature and the interaction of species and temperature. In this section, I am specifically interested in the thermally-responsive transcripts (temperature and species x temperature) so I subset the significant transcripts to examine these. 


```r
interaction.lms <- RxNlmAIC[which(Map(grepFunc, RxNlmAIC, term = "colonyAr:") == TRUE)]
other.lms <- RxNlmAIC[setdiff(names(RxNlmAIC), names(interaction.lms))]
temperature.lms <- other.lms[which(Map(grepFunc, other.lms, term = "val") == TRUE)]
colony.lms <- other.lms[setdiff(names(other.lms), names(temperature.lms))]
responsive.lms <- c(temperature.lms, interaction.lms)
rm(other.lms)

quadratic.lms <- RxNlmAIC[which(Map(grepFunc, RxNlmAIC, term = "2") == TRUE)]
```



|    Coefficient     |  Number.significant  |
|:------------------:|:--------------------:|
|       Total        |        10,525        |
|       Colony       |        1,473         |
|    Temperature     |        2,260         |
| Temperature:Colony |        6,792         |

Table: Number of transcripts with expression that depends on species, temperature or their interaction at 5% FDR  out of 98,186 total transcripts.


### Thermal-response functional types ###

The previous section identified the transcripts with thermally-responsive expression. In this section, I determine the shape of the expression response to temperature for each transcript. Catego
ries of expression response are:

* High - increase expression with temperature
* Low - decrease expression with temperature
* Intermediate - maximum expression at intermediate temperatures (14 - 28C)
* Bimodal - expressed greater than one standard deviation of expression at both low and high temperatures

I do this first for the thermally-responsive transcripts where there is no interaction with species. For the transcripts where thermal-responsive expression depends on species, I determine the functional type of the expression response separately for each species. 



```r
# calculate response type for responsive transcripts
interaction.response.type <- ldply(interaction.lms, .progress = "none", .inform = TRUE, RxNtype)
stopifnot(nrow(interaction.response.type) == length(interaction.lms))

# calculate response types for transcripts without interactions
temperature.response.type <- ldply(temperature.lms, .progress = "none", .inform = TRUE, RxNtype)
stopifnot(nrow(temperature.response.type) == length(temperature.lms))

# merge results
Ap.response.type <- rbind(interaction.response.type, temperature.response.type)
colnames(Ap.response.type)[which(colnames(Ap.response.type) == ".id")] <- "Transcript"
str(Ap.response.type)
```

```
## 'data.frame':	9052 obs. of  9 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100015|*|comp3543055_c0_seq1" "100067|*|comp3557646_c0_seq1" "100089|*|comp11313_c1_seq1" ...
##  $ A22.max   : num  38.5 0 NA 18.5 NA 0 NA 0 0 38.5 ...
##  $ A22.min   : num  0 20.5 NA 38.5 NA 23 NA 28.5 26 0 ...
##  $ A22.opt   : num  1.025 0.945 1 1.057 1 ...
##  $ A22.type  : chr  "High" "Bimodal" "NotResp" "Intermediate" ...
##  $ Ar.max    : num  0 38.5 0 38.5 38.5 18 0 0 NA NA ...
##  $ Ar.min    : num  38.5 0 20.5 13.5 18.5 38.5 38.5 38.5 NA NA ...
##  $ Ar.opt    : num  1.199 1.005 0.935 0.962 0.898 ...
##  $ Ar.type   : chr  "Low" "High" "Bimodal" "High" ...
```

```r
# save results to file
write.table(file = paste(resultsdir, "Ap_responsive_transcripts_", Sys.Date(), ".txt", sep = ""), Ap.response.type, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
```

Next, I compare the number of thermally-responsive in each response category between the two colonies. 


```r
A22.type.table <- table(Ap.response.type[ , 'A22.type'])
Ar.type.table <- table(Ap.response.type[ , 'Ar.type'])

# table
(Ap.type.table <- rbind(A22.type.table, Ar.type.table))
```

```
##                Bimodal High Intermediate  Low NotResp
## A22.type.table    1499 1238          909 4904     502
## Ar.type.table      868 1089         2606 3758     731
```

```r
# Pearson Chi-square test
chi1 <- chisq.test(Ap.type.table)
chi1
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  Ap.type.table
## X-squared = 1191, df = 4, p-value < 2.2e-16
```

```r
# Reorganize for plotting to show overlap among categories between colonies
type.table <- table(Acar = Ap.response.type[ , 'Ar.type'], Apic = Ap.response.type[ , 'A22.type'])
# reorder
tt2 <- type.table[c("Low", "Intermediate", "High", "Bimodal", "NotResp"), c("Low", "Intermediate", "High", "Bimodal", "NotResp")]

# plot
png("mosaic_plot.png")
mosaic.colors <- brewer.pal(5, "Blues")
par(mar = c(5.1, 3.1, 2.1, 7.5))
mosaicplot(tt2, ylab = expression(italic("A. picea")), xlab = "", main = "", color = mosaic.colors, cex.axis = 0.01)
labels <- c("Low", "Intermediate", "High", "Bimodal", "NotResp")
text(c(.1,.5,.75,.85,.99), par("usr")[3] - 0.05, srt = 45, adj = 1, labels = labels, xpd = TRUE, font = 1)
text(.5, par("usr")[4], "A. carolinensis", xpd = TRUE, font = 3)
legend(1.05, .75, fill = mosaic.colors, c("Low", "Intermediate", "High", "Bimodal", "NotResp"), xpd=TRUE)
dev.off()
```

```
## pdf 
##   2
```

Table with 'Total Transcripts' category



The question of biological interest is whether the marginal frequencies differ between the two species. Statistically, this can be addressed using the generalized [McNemar's test](http://en.wikipedia.org/wiki/McNemar's_test) of marginal homogeneity. 


```r
mh.test <- mh_test(as.table(tt2))

# calculate contribution of each off-diagonal to Z~0~
d1 <- sum(tt2[1,2:5]) - sum(tt2[2:5,1])
d2 <- sum(tt2[2,c(1,3:5)]) - sum(tt2[c(1,3:5),2])
d3 <- sum(tt2[3,c(1:2,4:5)]) - sum(tt2[c(1:2,4:5),3])
d4 <- sum(tt2[4,c(1:3,5)]) - sum(tt2[c(1:3,5),4])
d5 <- sum(tt2[5,1:4]) - sum(tt2[1:4,5])

# proportion of test statistic due to off-diagonals
dsum <- sum(abs(d1), abs(d2), abs(d3), abs(d4), abs(d5))
abs(d1)/dsum
```

```
## [1] 0.298
```

```r
abs(d2)/dsum
```

```
## [1] 0.441
```

```r
abs(d3)/dsum
```

```
## [1] 0.0387
```

```r
abs(d4)/dsum
```

```
## [1] 0.164
```

```r
abs(d5)/dsum
```

```
## [1] 0.0594
```

```r
# e2 (Acar Intermediate) contributes ~44% to Z~0~
```


```r
# overall mean for each class
rs <- rowSums(tt2)
cs <- colSums(tt2)
(gm <- (rs + cs) / sum(tt2 * 2))
```

```
##          Low Intermediate         High      Bimodal      NotResp 
##       0.4785       0.1942       0.1285       0.1307       0.0681
```

```r
# calculate observed values for each cell using overall mean for each expression type
Ec <- outer(gm, gm, "*") * sum(tt2)

# get deviations of observed from expected
Ec.dev <- tt2 - Ec

# calculate chi-squared deviation 
Ec.cells <- sign(tt2 - Ec) * (tt2 - Ec)^2 / Ec

md <- melt(Ec.cells)

mh_plot <- qplot(x=Apic, y=Acar, data=md, fill=value, geom="tile", ylim = rev(levels(md$Acar))) + 
  scale_fill_gradient2(limits=c(-600, 600)) 
mh_plot
```

![plot of chunk mh_plot](figure/mh_plot.png) 

The number of thermally-responsive in each response category differs between the colonies, with the msot transcripts expressed at *Low* temperatures in both colonies. For *ApVT*, an equal number of transcripts are expressed at *High* and *Bimodal*, followed by *Intermediate* transcripts. For *AcNC*, transcripts expressed at *Intermediate* temperatures are next most common, followed by *High* and *Bimodal*. 

Interestingly, nearly half of the *High* genes in *AcNC* are *Low* in *ApVT*, and vice versa. In contrast, most of the *Low* genes in one species are also *Low* in the other species.



|   &nbsp;   |   **ApVT**   |  &nbsp;  |    &nbsp;    |  &nbsp;  |  &nbsp;  |  &nbsp;  |  &nbsp;  |
|:----------:|:------------:|:--------:|:------------:|:--------:|:--------:|:--------:|:--------:|
|  **AcNC**  |              |   Low    | Intermediate |   High   | Bimodal  | NotResp  |  Total   |
|            |     Low      |   2654   |     160      |   405    |   299    |   240    |   3758   |
|            | Intermediate |   1314   |     542      |   147    |   536    |    67    |   2606   |
|            |     High     |   466    |      69      |   312    |   153    |    89    |   1089   |
|            |   Bimodal    |   185    |     120      |   160    |   297    |   106    |   868    |
|            |   NotResp    |   285    |      18      |   214    |   214    |    0     |   731    |
|            |    Total     |   4904   |     909      |   1238   |   1499   |   502    |   9052   |

Table: Number of transcripts with maximum expression at high, low, intermediate, both high and low (bimodal) temperatures or are not thermally-responsivefor each species and their overlap.

Table 4 shows the number of transcripts that fall into each expression type for each each species. The totals for each species include the 2260 transcripts that have consistent temperature responses between the two colonies. 

An interesting observation from the matched observations plot (Fig. ...) is that for **NotResp** genes in *A. carolinensis* the **High** and **Bimodal** categories in *A. picea* are over-represented. Finding that the mean expression level of these genes in *A. carolinensis* is greater than the expression level at the optimum temperature (19.5°C) in *A. picea* would be consistent with genetic assimilation.


```r
Ap.high.Ac.nr <- Ap.response.type[which(Ap.response.type$A22.type == "High" & Ap.response.type$Ar.type == "NotResp"), ]

# get mean expression level for Ar transcripts
Ap.high.Ac.nr.TPM <- TPM.dt.sub[Ap.high.Ac.nr$Transcript]
Ap.high.Ac.nr.mean.exp.NR <- ddply(Ap.high.Ac.nr.TPM, .(Transcript), summarize, meanTPM = mean(TPM))

# get expression at optimum for A22 transcripts

Ap.high.Ac.nr.opt.exp <- Ap.high.Ac.nr$A22.opt

t.test(Ap.high.Ac.nr.mean.exp.NR$meanTPM, Ap.high.Ac.nr.opt.exp)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ap.high.Ac.nr.mean.exp.NR$meanTPM and Ap.high.Ac.nr.opt.exp
## t = -55, df = 340, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.975 -0.907
## sample estimates:
## mean of x mean of y 
##    0.0824    1.0234
```

Contrary to expectations, mean expression of the non-responsive *A. carolinensis* genes is less than expression of the **High** genes in *A. picea*.

Next, I check the same for the **Bimodal** genes in *A. picea*.


```r
Ap.bim.Ac.nr <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal" & Ap.response.type$Ar.type == "NotResp"), ]

# get mean expression level for Ar transcripts
Ap.bim.Ac.nr.TPM <- TPM.dt.sub[Ap.bim.Ac.nr$Transcript]
Ap.bim.Ac.nr.mean.exp.NR <- ddply(Ap.bim.Ac.nr.TPM, .(Transcript), summarize, meanTPM = mean(TPM))

# get expression at optimum for A22 transcripts
Ap.bim.Ac.nr.opt.exp <- Ap.bim.Ac.nr$A22.opt

t.test(Ap.bim.Ac.nr.mean.exp.NR$meanTPM, Ap.bim.Ac.nr.opt.exp)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ap.bim.Ac.nr.mean.exp.NR$meanTPM and Ap.bim.Ac.nr.opt.exp
## t = -55.2, df = 375, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.966 -0.900
## sample estimates:
## mean of x mean of y 
##     0.105     1.038
```

As with the **High** transcripts, expression much higher for the responsive *A. picea* transcripts than for the non-responsive *A. carolinensis* transcripts.

## Species comparisons

In this section, I perform a number of comparisons of the thermal reactionomes between the *Ar* and *A22* species. Specifically, I:

- compare profiles of the temperature of maximum expression for thermally-responsive transcripts
- compare basal expression at optimal temperature for thermally-responsive genes
- compare thermal stability of *Intermediate* genes
- compare thermal sensitivity of *Bimodal* genes
- compare the temperature at which gene expression increases, T~on~, for *High* and *Low* genes
- evaluate the extent to which thermal reaction norms differ between species by mean and shape


### Compare temperature of maximum expression between species

Probability density function of peak expression for transcripts that differ in expression between *A22* and *Ar*


```r
# reshape data
Ap.df <- data.frame(Transcript = rep(Ap.response.type$Transcript, times = 2), colony = rep(c("ApVT", "AcNC"), each = length(Ap.response.type$Transcript)), max.val = c(Ap.response.type$A22.max, Ap.response.type$Ar.max))

maxexplot <- ggplot(Ap.df, aes(x=max.val, fill=colony)) + 
  geom_density(alpha=0.2, position="identity") + 
  labs(x = "Temperature of maximum expression", y = "Density") +
  theme(axis.title = element_text(size = rel(2))) +
  theme(axis.text = element_text(size = rel(1.2)))
print(maxexplot)
```

![plot of chunk max_exp_PDF](figure/max_exp_PDF.png) 

```r
# Figure for presentation
png("results/temp_max_expression.png")
maxexplot + theme(legend.position = "none")
```

```
## Error: could not open file 'results/temp_max_expression.png'
```

```r
dev.off()
```

```
## pdf 
##   2
```


### Compare basal expression at optimal temperature for thermally-responsive genes

Genes upregulated in response to thermal stress in one species may have greater basal levels of expression in the other species that experiences those stressful conditions more often. To test this hypothesis, we compared expression levels near the optimal temperature (19.5C) between the two species for genes in that are either in the 'High' or 'Low' expression group in the other species. Specifically, do genes upregulated at high temperatures in *A22* have greater basal expression at optimal temperatures in *Ar*?


```r
# list of transcripts that are 'high' expressed in A22
A22.high.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "High"), ]

# Compare expression at optimum temp (A22.opt) between colonies using t-test
t.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  A22.high.transcripts$A22.opt and A22.high.transcripts$Ar.opt
## t = -0.153, df = 2337, p-value = 0.878
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -1159   991
## sample estimates:
## mean of x mean of y 
##       385       469
```

```r
boxplot(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)
```

![plot of chunk optimum_expression_comparison](figure/optimum_expression_comparison1.png) 

```r
# T test on log-transformed values to control for outliers
t.test(log(A22.high.transcripts$A22.opt+1), log(A22.high.transcripts$Ar.opt+1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(A22.high.transcripts$A22.opt + 1) and log(A22.high.transcripts$Ar.opt + 1)
## t = -3.66, df = 2459, p-value = 0.0002603
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2697 -0.0814
## sample estimates:
## mean of x mean of y 
##      1.11      1.28
```

```r
boxplot(log(A22.high.transcripts$A22.opt+1), log(A22.high.transcripts$Ar.opt+1))
```

![plot of chunk optimum_expression_comparison](figure/optimum_expression_comparison2.png) 

The `t.test` fails to account for the many orders of magnitude difference in expression among transcripts, e.g. non-equal variances. This problem is the key issue in the analysis of differential expression (Bullard, Purdom, Hansen, and Dudoit, 2010; Anders, McCarthy, Chen, Okoniewski, Smyth, Huber, and Robinson, 2013). As my goal is simply to determine if expression is typically greater at optimal temperatures (19.5 C) in *Ar* than *A22* for genes that are up-regulated at high temperatures in *A22*, I use a non-parametric Wilcoxon signed rank-test


```r
w1 <- wilcox.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt, alternative = "two.sided", paired = TRUE, conf.int = TRUE)
w1
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  A22.high.transcripts$A22.opt and A22.high.transcripts$Ar.opt
## V = 198641, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.296 -0.171
## sample estimates:
## (pseudo)median 
##         -0.231
```

Consistent with expectation, there is greater expression at 19.5C for *Ar* than *A22* transcripts for the set of transcripts that are transcripts that are up-regulated at high temperatures in *A22*. Note that A22 had the larger library size so if this was due to TPM not correctly accounting for differences in reads between samples, we would expect to see a positive instead of negative value here.

Next I test the converse: do genes up-regulated at low temperatures in *Ar* have greater basal expression near optimal temperatures in *A22*?


```r
# list of transcripts that are 'high' expressed in Ar
Ar.low.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Low"), ]

# t-test with log values
t.test(log(Ar.low.transcripts$A22.opt+1), log(Ar.low.transcripts$Ar.opt+1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(Ar.low.transcripts$A22.opt + 1) and log(Ar.low.transcripts$Ar.opt + 1)
## t = -6.57, df = 7504, p-value = 5.21e-11
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1727 -0.0934
## sample estimates:
## mean of x mean of y 
##      1.26      1.39
```

```r
boxplot(log(Ar.low.transcripts$A22.opt+1), log(Ar.low.transcripts$Ar.opt+1))
```

![plot of chunk Ar_low_wilcoxon](figure/Ar_low_wilcoxon.png) 

```r
# Wilcoxon signed rank-test
w2 <- wilcox.test(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt, alternative = "two.sided", paired = TRUE, conf.int = TRUE)
w2
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  Ar.low.transcripts$A22.opt and Ar.low.transcripts$Ar.opt
## V = 2198150, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.455 -0.362
## sample estimates:
## (pseudo)median 
##         -0.406
```

Counter to expectations, expression at optimal temperatures is also greater in *Ar* than *A22* for transcripts upregulated at low temperatures in *Ar*. 

To confirm that there are not sample-level issues, I performed the same comparison using *Intermediate* expressed-genes where I do not expect to see a difference in expression.


```r
# list of transcripts that are 'Intermediate' expressed in either colony
Ap.int.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate" | Ap.response.type$A22.type == "Intermediate"), ]
# T test
t.test(log(Ap.int.transcripts$A22.opt+1), log(Ap.int.transcripts$Ar.opt+1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(Ap.int.transcripts$A22.opt + 1) and log(Ap.int.transcripts$Ar.opt + 1)
## t = -6.23, df = 5939, p-value = 4.935e-10
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.237 -0.123
## sample estimates:
## mean of x mean of y 
##      1.33      1.51
```

```r
# Wilcoxon signed rank-test
w3 <- wilcox.test(Ap.int.transcripts$A22.opt, Ap.int.transcripts$Ar.opt, alternative = "two.sided", paired = TRUE, conf.int = TRUE)
w3
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  Ap.int.transcripts$A22.opt and Ap.int.transcripts$Ar.opt
## V = 1433484, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.485 -0.365
## sample estimates:
## (pseudo)median 
##         -0.423
```

The non-parametric test for both comparisions also finds greater expression in *Ar* than *A22* at the optimal temperature.


### Compare thermal stability of *Intermediate* genes

*Intermediate* genes are those that have expression that is shut-down as conditions become stressful, likely non-essential molecular processes. We hypothesized that if the more southern *Ar* species was more thermally-tolerant than *A22*, transcripts with *Intermediate* expression (10-30C) would be active across a wider range of temperatures. To test this with our data, we calculated the standard deviation of the expression function for each temperature transcript that was 'Intermediate' expressed in each species.


```r
# extract 'Intermediate' expressed transcripts for A22
A22.int.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate"), ]
A22.int.lm <- RxNlmAIC[A22.int.transcripts$Transcript]
length(A22.int.lm)
```

```
## [1] 909
```

```r
# apply `transcriptSD` function to all transcripts
A22.int.sd <- unlist(Map(transcriptSD, A22.int.lm, colony = "A22"))
A22.int.sd <- data.frame(colony = "A22", exp_sd = A22.int.sd)

# repeat for Ar
Ar.int.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate"), ]
Ar.int.lm <- RxNlmAIC[Ar.int.transcripts$Transcript]
Ar.int.sd <- unlist(Map(transcriptSD, Ar.int.lm, colony = "Ar"))
Ar.int.sd <- data.frame(colony = "Ar", exp_sd = Ar.int.sd)

# T-test comparing standard deviation of expression between colonies
(t.varint <- t.test(Ar.int.sd$exp_sd, A22.int.sd$exp_sd, alternative = "two.sided"))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.int.sd$exp_sd and A22.int.sd$exp_sd
## t = 9.55, df = 1299, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.360 0.546
## sample estimates:
## mean of x mean of y 
##     10.36      9.91
```

Consistent with our hypothesis, *Intermediate* transcripts in *Ar* are expressed over a significantly wider range of temperatures than in *A22*. 

![plot of chunk plot_thermal_breadth](figure/plot_thermal_breadth.png) 

### Compare thermal sensitivity of *Bimodal* genes

As the converse of the above hypothesis, a species that is especially thermally-sensitive is likely to activate expression of molecular processes for thermal tolerance more quickly. We tested this using the same approach as for the *Intermediate* transcripts, but using the inverse of the *Bimodal* expressed transcripts. 


```r
# extract 'Bimodal' expressed transcripts for A22 species
A22.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$A22.type == "Bimodal"), "Transcript"]]
length(A22.bim.lm)
```

```
## [1] 1499
```

```r
# apply `transcriptSD` function to all transcripts
A22.bim.sd <- unlist(Map(transcriptSD, A22.bim.lm, colony = "A22"))

# repeat for Ar 
Ar.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$Ar.type == "Bimodal"), "Transcript"]]
length(Ar.bim.lm)
```

```
## [1] 868
```

```r
Ar.bim.sd <- unlist(Map(transcriptSD, Ar.bim.lm, colony = "Ar"))

# t-test to compare standard deviation of 'Bimodal' transcripts between colonies
(t.varbim <- t.test(Ar.bim.sd, A22.bim.sd))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.bim.sd and A22.bim.sd
## t = -0.855, df = 1528, p-value = 0.3925
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2126  0.0835
## sample estimates:
## mean of x mean of y 
##      12.6      12.7
```

No difference in the standard deviation of expression for bimodally-expressed transcripts between colonies.


### Compare the temperature at which gene expression increases, T~on~, for *High* and *Low* genes

Thermally-responsive genes could also differ in the temperatures at which they have increased or decreased expression in response to temperature changes. To examine this, I determine the temperature at which each responsive gene has the greated increase or decrease in expression.


```r
# extract TPM data for thermally-responsive transcripts
resp.TPM.dt.sub <- TPM.dt.sub[names(responsive.lms)]
setkey(resp.TPM.dt.sub, Transcript)
str(resp.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	199144 obs. of  10 variables:
##  $ Transcript       : chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ Length           : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM              : num  0 0.3566 0 0.7387 0.0741 ...
##  $ RPKM             : num  0 0.598 0 1.1258 0.0926 ...
##  $ KPKM             : num  0 0.598 0 1.1258 0.0926 ...
##  $ EstimatedNumKmers: num  0 157.7 0 366.2 60.9 ...
##  $ EstimatedNumReads: num  0 1.953 0 4.537 0.759 ...
##  $ sample           : chr  "A22-0" "Ar-0" "A22-3" "Ar-3" ...
##  $ val              : num  0 0 3.5 3.5 10.5 10.5 14 14 17.5 17.5 ...
##  $ colony           : Factor w/ 2 levels "A22","Ar": 1 2 1 2 1 2 1 2 1 2 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r
length(unique(resp.TPM.dt.sub$Transcript))
```

```
## [1] 9052
```

```r
# scale transcripts so can compare
resp.TPM.dt.sub[,TPM.scaled:=scale(TPM), by = Transcript]
```

```
##                          Transcript Length    TPM   RPKM   KPKM
##      1: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      2: 100008|*|comp137625_c0_seq2    208 0.3566 0.5980 0.5980
##      3: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      4: 100008|*|comp137625_c0_seq2    208 0.7387 1.1258 1.1258
##      5: 100008|*|comp137625_c0_seq2    208 0.0741 0.0926 0.0926
##     ---                                                        
## 199140:      9|*|comp147140_c0_seq1   9030 0.7302 1.2322 1.2322
## 199141:      9|*|comp147140_c0_seq1   9030 0.7258 1.0074 1.0074
## 199142:      9|*|comp147140_c0_seq1   9030 0.4534 0.8078 0.8078
## 199143:      9|*|comp147140_c0_seq1   9030 0.5649 0.7758 0.7758
## 199144:      9|*|comp147140_c0_seq1   9030 0.4528 0.8333 0.8333
##         EstimatedNumKmers EstimatedNumReads sample  val colony TPM.scaled
##      1:               0.0             0.000  A22-0  0.0    A22     -0.566
##      2:             157.7             1.953   Ar-0  0.0     Ar      1.131
##      3:               0.0             0.000  A22-3  3.5    A22     -0.566
##      4:             366.2             4.537   Ar-3  3.5     Ar      2.950
##      5:              60.9             0.759 A22-10 10.5    A22     -0.214
##     ---                                                                  
## 199140:           18024.7           222.869  Ar-31 31.5     Ar     -0.243
## 199141:           23968.5           299.433 A22-35 35.0    A22     -0.258
## 199142:            9055.4           112.144  Ar-35 35.0     Ar     -1.130
## 199143:           23690.1           295.836 A22-38 38.5    A22     -0.773
## 199144:           12462.4           153.987  Ar-38 38.5     Ar     -1.132
```

```r
# rename colonies
resp.TPM.dt.sub$colony2 <- ifelse(resp.TPM.dt.sub$colony == "A22", "ApVT", "AcNC")
```

Predict expression for responsive transcripts.


```r
# apply predFunc to all responsive transcripts
resp.TPM.dt.sub.pred <- ddply(resp.TPM.dt.sub, .(Transcript), .inform="TRUE", predFunc)

# setkey to Transcript and colony 
resp.TPM.dt.sub.pred <- data.table(resp.TPM.dt.sub.pred)
setkey(resp.TPM.dt.sub.pred, Transcript, colony)
```


Does the mean expression level for all responsive genes differ among colonies?


```r
t.test(log(Ap.response.type$A22.opt), log(Ap.response.type$Ar.opt))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(Ap.response.type$A22.opt) and log(Ap.response.type$Ar.opt)
## t = -9.09, df = 18041, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.185 -0.120
## sample estimates:
## mean of x mean of y 
##     0.708     0.861
```

```r
# expression level at optimum temp greater in Ar than A22

# non-parametric test
w4 = wilcox.test(Ap.response.type$A22.opt, 
                 Ap.response.type$Ar.opt,
                 alternative = "two.sided", paired = TRUE, conf.int = TRUE)


mean.exp <- ddply(resp.TPM.dt.sub.pred, .(colony, Transcript), summarise,
                  mean.TPM = mean(pTPM))

t.test(log(mean.exp[which(mean.exp$colony == "A22"), "mean.TPM"]), 
       log(mean.exp[which(mean.exp$colony == "Ar"), "mean.TPM"]))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(mean.exp[which(mean.exp$colony == "A22"), "mean.TPM"]) and log(mean.exp[which(mean.exp$colony == "Ar"), "mean.TPM"])
## t = -6.53, df = 18077, p-value = 6.938e-11
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1399 -0.0753
## sample estimates:
## mean of x mean of y 
##     0.747     0.855
```

```r
# non-parametric test
w5 = wilcox.test(mean.exp[which(mean.exp$colony == "A22"), "mean.TPM"],
                 mean.exp[which(mean.exp$colony == "Ar"), "mean.TPM"],
                 alternative = "two.sided", paired = TRUE, conf.int = TRUE)


plot(log(mean.exp[which(mean.exp$colony == "A22"), "mean.TPM"]), 
     log(mean.exp[which(mean.exp$colony == "Ar"), "mean.TPM"]))
```

![plot of chunk mean_exp](figure/mean_exp.png) 




For next analyses, extract list of gene names by response type.


```r
A22.high.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "High"), "Transcript"]
Ar.high.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "High"), "Transcript"]

A22.low.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "Low"), "Transcript"]
Ar.low.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Low"), "Transcript"]

A22.bim.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal"), "Transcript"]
Ar.bim.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Bimodal"), "Transcript"]

A22.int.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate"), "Transcript"]
Ar.int.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate"), "Transcript"]
```


Calculate T~on~ for *High* genes in each species. I use the observed (rather than predicted) values for each species as by nature of fitting a quadratic function, the maximum change tends to be at the extremes (38.5 or 0) for the predicted values.



```r
# transcripts expressed at *High* and *Bimodal* temperatures in A22
A22.high.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(A22.high.transcripts, A22.bim.transcripts), "A22")]
setkey(A22.high.TPM.dt.sub, Transcript)
str(A22.high.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	30107 obs. of  5 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0 0 0.0741 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.01 1.01 1.01 1.01 1.01 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r
# make data.frame for results
l1 <- length(unique(A22.high.TPM.dt.sub$Transcript))
A22.high.T_on <- data.frame(Transcript = unique(A22.high.TPM.dt.sub$Transcript), colony = rep("ApVT", length = l1), type = rep("High", length = l1), T_on = NA)

# loop across transcripts, calculating T_on

for(i in unique(A22.high.TPM.dt.sub$Transcript)) {
    subdf <- A22.high.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val > 14), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM))))+1, val]
    A22.high.T_on[which(A22.high.T_on$Transcript == i), "T_on"] <- T_on
}   

# repeat for Ar
Ar.high.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(Ar.high.transcripts, Ar.bim.transcripts), "Ar")]
setkey(Ar.high.TPM.dt.sub, Transcript)
str(Ar.high.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	21527 obs. of  5 variables:
##  $ Transcript: chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 2 2 2 2 2 2 2 2 2 2 ...
##  $ TPM       : num  0 0 0 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  0.997 0.999 1.001 1.002 1.004 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r
l2 <- length(unique(Ar.high.TPM.dt.sub$Transcript))
Ar.high.T_on <- data.frame(Transcript = unique(Ar.high.TPM.dt.sub$Transcript), colony = rep("AcNC", length = l2), type = rep("High", length = l2), T_on = NA)

for(i in unique(Ar.high.TPM.dt.sub$Transcript)) {
    subdf <- Ar.high.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val > 14), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM))))+1, val]
    Ar.high.T_on[which(Ar.high.T_on$Transcript == i), "T_on"] <- T_on
}   

# determine if T~on~ is greater in *A22* or *Ar* for *High* genes.

(T_on.high.ttest <- t.test(Ar.high.T_on$T_on, A22.high.T_on$T_on))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.high.T_on$T_on and A22.high.T_on$T_on
## t = 5.82, df = 4128, p-value = 6.166e-09
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.578 1.164
## sample estimates:
## mean of x mean of y 
##      33.8      32.9
```

Genes with increased expression at *High* temperatures are on average turned on 1°C higher in *AcNC* than *ApVT*.

Repeat analysis for *Low* genes.


```r
# transcripts expressed at *Low* temperatures in A22
A22.low.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(A22.low.transcripts, A22.bim.transcripts), "A22")]
setkey(A22.low.TPM.dt.sub, Transcript)
str(A22.low.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	70433 obs. of  5 variables:
##  $ Transcript: chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0.751 0 0 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.451 1.27 1.048 0.988 0.954 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r
# make data.frame for results
l3 <- length(unique(A22.low.TPM.dt.sub$Transcript))
A22.low.T_on <- data.frame(Transcript = unique(A22.low.TPM.dt.sub$Transcript), colony = rep("ApVT", length = l3), type = rep("Low", length = l3), T_on = NA)

# loop across transcripts, calculating T_on
for(i in unique(A22.low.TPM.dt.sub$Transcript)) {
    subdf <- A22.low.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val < 21), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM))))+1, val]
    A22.low.T_on[which(A22.low.T_on$Transcript == i), "T_on"] <- T_on
}   

# repeat for Ar
Ar.low.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(Ar.low.transcripts, Ar.bim.transcripts), "A22")]
setkey(Ar.low.TPM.dt.sub, Transcript)
str(Ar.low.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	50886 obs. of  5 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0 0 0.0741 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.01 1.01 1.01 1.01 1.01 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r
l4 <- length(unique(Ar.low.TPM.dt.sub$Transcript))
Ar.low.T_on <- data.frame(Transcript = unique(Ar.low.TPM.dt.sub$Transcript), colony = rep("AcNC", length = l4), type = rep("Low", length = length(unique(Ar.low.TPM.dt.sub$Transcript))), T_on = NA)

for(i in unique(Ar.low.TPM.dt.sub$Transcript)) {
    subdf <- Ar.low.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val < 21), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM))))+1, val]
    Ar.low.T_on[which(Ar.low.T_on$Transcript == i), "T_on"] <- T_on
}   

(T_on.low.ttest <- t.test(Ar.low.T_on$T_on, A22.low.T_on$T_on))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.low.T_on$T_on and A22.low.T_on$T_on
## t = -1.9, df = 10077, p-value = 0.05718
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.30369  0.00456
## sample estimates:
## mean of x mean of y 
##      12.8      12.9
```

Genes with increased expression at *Low* temperatures are on average turned on 0.2°C higher in *ApVT* than *AcNC*.

Visualize T~on~ for both *Low* and *High* genes on the same plot.

![plot of chunk plot_T_on](figure/plot_T_on.png) 


Repeat analysis for when **Intermediate** genes are turned *off*.


```r
# transcripts expressed at *Intermediate* temperatures in A22
A22.int.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(A22.int.transcripts, "A22")]
setkey(A22.int.TPM.dt.sub, Transcript)
str(A22.int.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	9999 obs. of  5 variables:
##  $ Transcript: chr  "100089|*|comp11313_c1_seq1" "100089|*|comp11313_c1_seq1" "100089|*|comp11313_c1_seq1" "100089|*|comp11313_c1_seq1" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0 0 0 0 0.264 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  0.987 1.011 1.044 1.053 1.057 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r
# make data.frame for results
l5 <- length(unique(A22.int.TPM.dt.sub$Transcript))
A22.int.T_off <- data.frame(Transcript = unique(A22.int.TPM.dt.sub$Transcript), colony = rep("ApVT", length = l5), type = rep("Low", length = l5), T_off_high = NA, T_off_low = NA)

# loop across transcripts, calculating T_off at both high and low temperatures
for(i in unique(A22.int.TPM.dt.sub$Transcript)) {
    subdf <- A22.int.TPM.dt.sub[i]
    
    subdf.high <- subdf[which(subdf$val >= 21), ]
    T_off_high <- subdf.high[median(which(diff(subdf.high$TPM) == min(diff(subdf.high$TPM)))), val]
    A22.int.T_off[which(A22.int.T_off$Transcript == i), "T_off_high"] <- T_off_high
    
    subdf.low <- subdf[which(subdf$val <= 21), ]
    T_off_low <- subdf.low[median(which(diff(subdf.low$TPM) == min(diff(subdf.low$TPM)))), val]
    A22.int.T_off[which(A22.int.T_off$Transcript == i), "T_off_low"] <- T_off_low
}   

# repeat for Ar
Ar.int.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(Ar.int.transcripts, "Ar")]
setkey(Ar.int.TPM.dt.sub, Transcript)
str(Ar.int.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	28666 obs. of  5 variables:
##  $ Transcript: chr  "100148|*|comp125464_c0_seq1" "100148|*|comp125464_c0_seq1" "100148|*|comp125464_c0_seq1" "100148|*|comp125464_c0_seq1" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 2 2 2 2 2 2 2 2 2 2 ...
##  $ TPM       : num  0 0 0 0 0.0644 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  0.998 1.003 1.01 1.012 1.012 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r
# make data.frame for results
l6 <- length(unique(Ar.int.TPM.dt.sub$Transcript))
Ar.int.T_off <- data.frame(Transcript = unique(Ar.int.TPM.dt.sub$Transcript), colony = rep("ApVT", length = l6), type = rep("Low", length = l6), T_off_high = NA, T_off_low = NA)

# loop across transcripts, calculating T_off at both high and low temperatures
for(i in unique(Ar.int.TPM.dt.sub$Transcript)) {
    subdf <- Ar.int.TPM.dt.sub[i]
    
    subdf.high <- subdf[which(subdf$val >= 21), ]
    T_off_high <- subdf.high[median(which(diff(subdf.high$TPM) == min(diff(subdf.high$TPM)))), val]
    Ar.int.T_off[which(Ar.int.T_off$Transcript == i), "T_off_high"] <- T_off_high
    
    subdf.low <- subdf[which(subdf$val <= 21), ]
    T_off_low <- subdf.low[median(which(diff(subdf.low$TPM) == min(diff(subdf.low$TPM)))), val]
    Ar.int.T_off[which(Ar.int.T_off$Transcript == i), "T_off_low"] <- T_off_low
}

# compare T_off among species
(T_off.int.high.ttest <- t.test(Ar.int.T_off$T_off_high, A22.int.T_off$T_off_high))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.int.T_off$T_off_high and A22.int.T_off$T_off_high
## t = 7.02, df = 1655, p-value = 3.325e-12
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.871 1.547
## sample estimates:
## mean of x mean of y 
##      28.3      27.1
```

```r
(T_off.int.low.ttest <- t.test(Ar.int.T_off$T_off_low, A22.int.T_off$T_off_low))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.int.T_off$T_off_low and A22.int.T_off$T_off_low
## t = 2.17, df = 1574, p-value = 0.02979
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.0501 0.9722
## sample estimates:
## mean of x mean of y 
##      11.1      10.6
```

As with genes turned on at high temperatures, *A. carolinensis* **Intermediate** genes are down-regulated on average at 1.2C higher temperatures than *A. picea* **Intermediate** genes. 



### Evaluate the extent to which differences in thermal reaction norms are due to mean or shape

For the transcripts that differed in thermal responsiveness due to temperature, was the difference primarily due to differences in the mean value of expression, slope, curvature of a higher order effect? To test this, I rougly follow Murren, Maclean, Diamond, Steiner, Heskel, Handelsman, Ghalambor, Auld, Callahan, Pfennig, Relyea, Schlichting, and Kingsolver (2014) by defining differences among reation norms for individual genes due to changes in the overall mean, slope, curvature and all higher-order shape differences (i.e. wiggle).

- *Mean, M*: overall difference in the mean expression value across all temperatures
- *Slope, S*: difference in overall slope
- *Curvature, C*: average difference in curvature of the reaction norm
- *Wiggle, W*: variability in shape not captured by the previous three measures


```r
varshape.out <- ldply(interaction.lms, .progress = "none", RxNvarshape)

# delta mean
boxplot(varshape.out$A22.mean, varshape.out$Ar.mean)
```

![plot of chunk varshape](figure/varshape1.png) 

```r
boxplot(log(varshape.out$A22.mean), log(varshape.out$Ar.mean))
```

![plot of chunk varshape](figure/varshape2.png) 

```r
t.mean <- t.test(log(varshape.out$A22.mean), log(varshape.out$Ar.mean), paired = TRUE)
t.mean
```

```
## 
## 	Paired t-test
## 
## data:  log(varshape.out$A22.mean) and log(varshape.out$Ar.mean)
## t = -10.7, df = 6791, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1069 -0.0737
## sample estimates:
## mean of the differences 
##                 -0.0903
```

```r
# delta slope
boxplot(log(varshape.out$A22.slope), log(varshape.out$Ar.slope))
```

```
## Warning: Outlier (-Inf) in boxplot 1 is not drawn
## Warning: Outlier (-Inf) in boxplot 2 is not drawn
```

![plot of chunk varshape](figure/varshape3.png) 

```r
t.slope <- t.test(log(varshape.out$A22.slope + 0.1), log(varshape.out$Ar.slope + 0.1), paired = TRUE)
t.slope
```

```
## 
## 	Paired t-test
## 
## data:  log(varshape.out$A22.slope + 0.1) and log(varshape.out$Ar.slope + 0.1)
## t = 7.71, df = 6791, p-value = 1.496e-14
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.0214 0.0360
## sample estimates:
## mean of the differences 
##                  0.0287
```

```r
# delta curvature
boxplot(log(varshape.out$A22.curve), log(varshape.out$Ar.curve))
```

![plot of chunk varshape](figure/varshape4.png) 

```r
t.curvature <- t.test(log(varshape.out$A22.curve + 0.1), log(varshape.out$Ar.curve + 0.1), paired = TRUE)
t.curvature
```

```
## 
## 	Paired t-test
## 
## data:  log(varshape.out$A22.curve + 0.1) and log(varshape.out$Ar.curve + 0.1)
## t = 0.763, df = 6791, p-value = 0.4454
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.00168  0.00381
## sample estimates:
## mean of the differences 
##                 0.00107
```

```r
# delta wiggle - small values so use un-transformed
boxplot(varshape.out$A22.wiggle, varshape.out$Ar.wiggle)
```

![plot of chunk varshape](figure/varshape5.png) 

```r
t.wiggle <- t.test(varshape.out$A22.wiggle, varshape.out$Ar.wiggle, paired = TRUE)
t.wiggle
```

```
## 
## 	Paired t-test
## 
## data:  varshape.out$A22.wiggle and varshape.out$Ar.wiggle
## t = -0.782, df = 6791, p-value = 0.4344
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.000375  0.000161
## sample estimates:
## mean of the differences 
##               -0.000107
```

Reaction norms differ between species by mean and slope, but not by curvature or wiggle.

Next, I partition the differences in the reaction norms into the variation explained by differences in the trait means and slope, as curvature and wiggle did not differ between species.


```r
# calculate differences in mean, slope, curvature and wiggle between colonies for each transcript
# mean = mean1 - mean2
varshape.out$lmean <- log(varshape.out$A22.mean) - log(varshape.out$Ar.mean)

# slope = slope1 - slope2
varshape.out$lslope <- log(varshape.out$A22.slope + 0.1) - log(varshape.out$Ar.slope + 0.1)

# take absolute value of each value and sum to get total differences between reaction norms
varshape.out$ltotal <- abs(varshape.out$lmean) + abs(varshape.out$lslope)

# variation in reaction norms due to differences in mean
varshape.out$prop.lmean <- abs(varshape.out$lmean) / varshape.out$ltotal
varshape.out$prop.lslope <- abs(varshape.out$lslope) / varshape.out$ltotal

# Mean proportion and 95% CI of total variation of each measure
mean(varshape.out$prop.lmean)
```

```
## [1] 0.773
```

```r
(qlmean <- quantile(varshape.out$prop.lmean, probs = c(0.05, 0.5, 0.95)))
```

```
##    5%   50%   95% 
## 0.448 0.817 0.985
```

```r
mean(varshape.out$prop.lslope)
```

```
## [1] 0.227
```

```r
(qlslope <- quantile(varshape.out$prop.lslope, probs = c(0.05, 0.5, 0.95)))
```

```
##     5%    50%    95% 
## 0.0154 0.1832 0.5517
```

From this analysis, about 80% of the differences in reaction norms between species are due to changes in the mean, with the remainder being due to changes in slope.



## Gene set enrichment analysis for thermal-responsive genes

### Functional annotation

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  


```r
setkey(annotation.table, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
signif.transcripts <- annotation.table[signif.transcripts]
setnames(signif.transcripts, "Sequence.Name", "Transcript")
```


```r
responsive.lms.ann.type <- signif.transcripts[names(responsive.lms)]
# combine responsive.lms.ann with "type" information
responsive.lms.ann.type <- merge(responsive.lms.ann.type, Ap.response.type[,c("Transcript", "A22.type", "Ar.type")], by = "Transcript", all = TRUE)
str(responsive.lms.ann.type)
```

```
## Classes 'data.table' and 'data.frame':	9052 obs. of  17 variables:
##  $ Transcript           : chr  "100008|*|comp137625_c0_seq2" "100015|*|comp3543055_c0_seq1" "100067|*|comp3557646_c0_seq1" "100089|*|comp11313_c1_seq1" ...
##  $ sequence.length      : int  208 208 208 208 1321 208 207 1320 1320 207 ...
##  $ best.hit.to.nr       : chr  "-" "gi|121608385|ref|YP_996192.1| transposase, IS4 family protein " "gi|493136460|ref|WP_006154899.1| tyrosyl-tRNA synthetase " "gi|497544620|ref|WP_009858818.1| lipid-A-disaccharide synthase " ...
##  $ hit.length           : chr  "-" "25" "68" "67" ...
##  $ E.value              : chr  "-" "7.22e-06" "1.5e-39" "3.37e-23" ...
##  $ Bit.score            : chr  "-" "57.514374" "159.812994" "111.355753" ...
##  $ GO.Biological.Process: chr  "-" "-" "-" "GO:0009245 lipid A biosynthetic process" ...
##  $ GO.Cellular.Component: chr  "-" "-" "-" "GO:0009276 Gram-negative-bacterium-type cell wall" ...
##  $ GO.Molecular.Function: chr  "-" "-" "-" "GO:0008915 lipid-A-disaccharide synthase activity" ...
##  $ Enzyme               : chr  "-" "-" "-" "-" ...
##  $ Domain               : chr  "-" "-" "-" "-" ...
##  $ annotation.type      : chr  "" "" "" "GO only" ...
##  $ pval                 : num  0.003598 0.000917 0.0021 0.001571 0.001918 ...
##  $ adj.r.squared        : num  0.523 0.604 0.557 0.574 0.562 ...
##  $ padj                 : num  0.0374 0.0135 0.0252 0.0203 0.0236 ...
##  $ A22.type             : chr  "High" "Bimodal" "NotResp" "Intermediate" ...
##  $ Ar.type              : chr  "Low" "High" "Bimodal" "High" ...
##  - attr(*, "sorted")= chr "Transcript"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
write.csv(responsive.lms.ann.type, file = paste(resultsdir, "Ap_responsive_genes_", Sys.Date(), ".csv", sep = ""), row.names = FALSE)
```


### Proportion of responsive transcripts annotated

Of the responsive transcripts, 51% are annotated, while 51% of all transcripts are annotated.


length(which(responsive.lms.ann.type$Transcript %in% annotation.table[which(annotation.table$best.hit.to.nr != "-"), "Sequence.Name"]))

best.hit.to.nr != "-")

### Candidate gene enrichment

First, I explore of candidate genes with the GO term "response to stress" are enriched in the responsive data set. 


```r
# GO 'response to stress' hits in responsive transcripts
GO0006950.responsive <- responsive.lms.ann.type[grep("GO:0006950", responsive.lms.ann.type$GO.Biological.Process), list(Transcript, best.hit.to.nr, A22.type, Ar.type)]

# in high category
GO0006950.responsive[union(with(GO0006950.responsive, grep("High", Ar.type)), with(GO0006950.responsive, grep("High", A22.type))), ]
```

```
##                      Transcript
##  1:           1504|*|Contig2729
##  2:  15115|*|comp132715_c0_seq1
##  3:  20095|*|comp137823_c1_seq1
##  4:   2087|*|comp150483_c5_seq1
##  5:   3269|*|comp141800_c0_seq1
##  6: 75624|*|comp2836178_c0_seq1
##  7:   7632|*|comp148048_c0_seq1
##  8:   8886|*|comp150172_c1_seq4
##  9:  21598|*|comp142101_c0_seq1
## 10:  23441|*|comp114823_c1_seq1
## 11:  32312|*|comp933733_c0_seq1
## 12: 37154|*|comp1975086_c0_seq1
## 13: 47691|*|comp1460938_c0_seq1
## 14: 51985|*|comp1012776_c0_seq1
##                                                               best.hit.to.nr
##  1:               gi|332023134|gb|EGI63390.1| Sugar transporter ERD6-like 6 
##  2:                       gi|194716766|gb|ACF93232.1| heat shock protein 90 
##  3:                          gi|337757286|emb|CBZ98843.1| 60 kDa chaperonin 
##  4:         gi|332022897|gb|EGI63169.1| Protein lethal(2)essential for life 
##  5:         gi|332021988|gb|EGI62314.1| Heat shock 70 kDa protein cognate 4 
##  6:                               gi|50418863|ref|XP_457952.1| DEHA2C06072p 
##  7:             gi|322789999|gb|EFZ15075.1| hypothetical protein SINV_03446 
##  8:             gi|332029691|gb|EGI69570.1| G-protein coupled receptor Mth2 
##  9:         gi|332018201|gb|EGI58806.1| Protein lethal(2)essential for life 
## 10:      gi|443696809|gb|ELT97425.1| hypothetical protein CAPTEDRAFT_194915 
## 11:     gi|367054010|ref|XP_003657383.1| hypothetical protein THITE_2156506 
## 12:             gi|46115086|ref|XP_383561.1| hypothetical protein FG03385.1 
## 13: gi|302922354|ref|XP_003053448.1| hypothetical protein NECHADRAFT_102357 
## 14:                       gi|227018528|gb|ACP18866.1| heat shock protein 30 
##     A22.type Ar.type
##  1:     High    High
##  2:  Bimodal    High
##  3:  NotResp    High
##  4:      Low    High
##  5:     High    High
##  6:  NotResp    High
##  7:      Low    High
##  8:      Low    High
##  9:     High Bimodal
## 10:     High     Low
## 11:     High Bimodal
## 12:     High NotResp
## 13:     High NotResp
## 14:     High NotResp
```

```r
# in low category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Low", Ar.type)), with(GO0006950.responsive, grep("Low", A22.type))), ]
```

```
##                      Transcript
##  1:   1038|*|comp150483_c5_seq3
##  2:  11281|*|comp146961_c0_seq1
##  3:     14|*|comp150262_c0_seq1
##  4:          19475|*|Contig1438
##  5:  20215|*|comp145360_c0_seq1
##  6:  23441|*|comp114823_c1_seq1
##  7:   2604|*|comp148324_c0_seq4
##  8:   4273|*|comp150636_c5_seq1
##  9: 50934|*|comp3428507_c0_seq1
## 10:    6075|*|comp92770_c0_seq1
## 11:   6438|*|comp150878_c2_seq2
## 12:   6689|*|comp141130_c0_seq2
## 13:  80544|*|comp132706_c0_seq1
## 14:           9372|*|Contig4757
## 15:  12704|*|comp144775_c1_seq1
## 16:   1656|*|comp147700_c0_seq1
## 17:  17710|*|comp150271_c3_seq3
## 18:   2087|*|comp150483_c5_seq1
## 19:   3995|*|comp145243_c0_seq1
## 20:    552|*|comp147487_c0_seq1
## 21:   7632|*|comp148048_c0_seq1
## 22:   8886|*|comp150172_c1_seq4
## 23:   9316|*|comp147545_c4_seq2
##                      Transcript
##                                                                                    best.hit.to.nr
##  1:                              gi|332022897|gb|EGI63169.1| Protein lethal(2)essential for life 
##  2:                                  gi|332029692|gb|EGI69571.1| G-protein coupled receptor Mth2 
##  3: gi|332019420|gb|EGI59904.1| Putative fat-like cadherin-related tumor suppressor-like protein 
##  4:                                  gi|322799248|gb|EFZ20646.1| hypothetical protein SINV_03807 
##  5:                                  gi|332020393|gb|EGI60813.1| G-protein coupled receptor Mth2 
##  6:                           gi|443696809|gb|ELT97425.1| hypothetical protein CAPTEDRAFT_194915 
##  7:                                   gi|307176228|gb|EFN65864.1| hypothetical protein EAG_10145 
##  8:                    gi|332026123|gb|EGI66271.1| Multiple inositol polyphosphate phosphatase 1 
##  9:                                                          gi|15010456|gb|AAK77276.1| GH05807p 
## 10:             gi|332030037|gb|EGI69862.1| Serine/threonine-protein kinase PINK1, mitochondrial 
## 11:                            gi|307188496|gb|EFN73233.1| Muscarinic acetylcholine receptor DM1 
## 12:                      gi|332016397|gb|EGI57310.1| Mitochondrial import receptor subunit TOM70 
## 13:                                            gi|121605727|ref|YP_983056.1| OsmC family protein 
## 14:                                  gi|332029691|gb|EGI69570.1| G-protein coupled receptor Mth2 
## 15:         gi|380028536|ref|XP_003697954.1| PREDICTED: protein lethal(2)essential for life-like 
## 16:                              gi|332030522|gb|EGI70210.1| Heat shock 70 kDa protein cognate 3 
## 17:                                  gi|322799248|gb|EFZ20646.1| hypothetical protein SINV_03807 
## 18:                              gi|332022897|gb|EGI63169.1| Protein lethal(2)essential for life 
## 19:                              gi|307211659|gb|EFN87680.1| Heat shock 70 kDa protein cognate 5 
## 20:                             gi|332026309|gb|EGI66443.1| RhoA activator C11orf59-like protein 
## 21:                                  gi|322789999|gb|EFZ15075.1| hypothetical protein SINV_03446 
## 22:                                  gi|332029691|gb|EGI69570.1| G-protein coupled receptor Mth2 
## 23:                              gi|332020093|gb|EGI60539.1| Heat shock factor-binding protein 1 
##                                                                                    best.hit.to.nr
##         A22.type      Ar.type
##  1:          Low          Low
##  2:          Low          Low
##  3:          Low          Low
##  4:          Low          Low
##  5:          Low          Low
##  6:         High          Low
##  7:          Low          Low
##  8:          Low          Low
##  9:          Low          Low
## 10:          Low          Low
## 11:          Low          Low
## 12:          Low          Low
## 13: Intermediate          Low
## 14:          Low          Low
## 15:          Low Intermediate
## 16:          Low Intermediate
## 17:          Low Intermediate
## 18:          Low         High
## 19:          Low Intermediate
## 20:          Low Intermediate
## 21:          Low         High
## 22:          Low         High
## 23:          Low Intermediate
##         A22.type      Ar.type
```

```r
# in bimodal category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Bimodal", Ar.type)), with(GO0006950.responsive, grep("Bimodal", A22.type))), ]
```

```
##                    Transcript
## 1:  19778|*|comp97601_c0_seq1
## 2: 21598|*|comp142101_c0_seq1
## 3: 32312|*|comp933733_c0_seq1
## 4: 58246|*|comp109744_c0_seq1
## 5: 15115|*|comp132715_c0_seq1
## 6: 21384|*|comp149042_c0_seq3
##                                                           best.hit.to.nr
## 1:          gi|322796169|gb|EFZ18745.1| hypothetical protein SINV_07491 
## 2:      gi|332018201|gb|EGI58806.1| Protein lethal(2)essential for life 
## 3:  gi|367054010|ref|XP_003657383.1| hypothetical protein THITE_2156506 
## 4:            gi|493322437|ref|WP_006279741.1| molecular chaperone DnaK 
## 5:                    gi|194716766|gb|ACF93232.1| heat shock protein 90 
## 6: gi|396467618|ref|XP_003837992.1| hypothetical protein LEMA_P120390.1 
##    A22.type Ar.type
## 1:  Bimodal Bimodal
## 2:     High Bimodal
## 3:     High Bimodal
## 4:  NotResp Bimodal
## 5:  Bimodal    High
## 6:  Bimodal NotResp
```

```r
# in intermediate category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Intermediate", Ar.type)), with(GO0006950.responsive, grep("Intermediate", A22.type))), ]
```

```
##                    Transcript
## 1: 11087|*|comp141315_c0_seq1
## 2: 12704|*|comp144775_c1_seq1
## 3:  1656|*|comp147700_c0_seq1
## 4: 17710|*|comp150271_c3_seq3
## 5:  3995|*|comp145243_c0_seq1
## 6:   552|*|comp147487_c0_seq1
## 7:  9316|*|comp147545_c4_seq2
## 8:    97|*|comp150840_c0_seq3
## 9: 80544|*|comp132706_c0_seq1
##                                                                                      best.hit.to.nr
## 1:                                     gi|322792301|gb|EFZ16285.1| hypothetical protein SINV_03698 
## 2:            gi|380028536|ref|XP_003697954.1| PREDICTED: protein lethal(2)essential for life-like 
## 3:                                 gi|332030522|gb|EGI70210.1| Heat shock 70 kDa protein cognate 3 
## 4:                                     gi|322799248|gb|EFZ20646.1| hypothetical protein SINV_03807 
## 5:                                 gi|307211659|gb|EFN87680.1| Heat shock 70 kDa protein cognate 5 
## 6:                                gi|332026309|gb|EGI66443.1| RhoA activator C11orf59-like protein 
## 7:                                 gi|332020093|gb|EGI60539.1| Heat shock factor-binding protein 1 
## 8: gi|350402309|ref|XP_003486440.1| PREDICTED: probable G-protein coupled receptor Mth-like 1-like 
## 9:                                               gi|121605727|ref|YP_983056.1| OsmC family protein 
##        A22.type      Ar.type
## 1: Intermediate Intermediate
## 2:          Low Intermediate
## 3:          Low Intermediate
## 4:          Low Intermediate
## 5:          Low Intermediate
## 6:          Low Intermediate
## 7:          Low Intermediate
## 8: Intermediate Intermediate
## 9: Intermediate          Low
```

```r
# Chi-square test to see if 'response to stress' related genes overrepresented in responsive.lms compared to full list
resp.stress.responsive.count <- nrow(responsive.lms.ann.type[grep("GO:0006950", responsive.lms.ann.type$GO.Biological.Process), list(Transcript, best.hit.to.nr)])
# GO 'response to stress' hits in all transcripts
resp.stress.all.count <- nrow(annotation.table[grep("GO:0006950", annotation.table$GO.Biological.Process), list(Sequence.Name, best.hit.to.nr)])

GO.stress.table <- matrix(rbind(resp.stress.responsive.count, nrow(responsive.lms.ann.type) - resp.stress.responsive.count, resp.stress.all.count, nrow(annotation.table) - resp.stress.all.count), nrow = 2)

GO.stress.Xsq <- chisq.test(GO.stress.table)
GO.stress.Xsq
```

```
## 
## 	Pearson's Chi-squared test with Yates' continuity correction
## 
## data:  GO.stress.table
## X-squared = 2.97, df = 1, p-value = 0.08469
```

While 38 "response to stress" genes are in the responsive set, this is not enriched compared to the whole transcriptome. 


```r
hsp_responsive <- responsive.lms.ann.type[grep("shock", responsive.lms.ann.type$best.hit.to.nr), list(Transcript, best.hit.to.nr, A22.type, Ar.type)]
hsp_responsive
```

```
##                     Transcript
## 1:  15115|*|comp132715_c0_seq1
## 2:   1656|*|comp147700_c0_seq1
## 3:  20675|*|comp147923_c0_seq1
## 4:   3269|*|comp141800_c0_seq1
## 5:   3995|*|comp145243_c0_seq1
## 6: 51985|*|comp1012776_c0_seq1
## 7:   9316|*|comp147545_c4_seq2
##                                                                           best.hit.to.nr
## 1:                                    gi|194716766|gb|ACF93232.1| heat shock protein 90 
## 2:                      gi|332030522|gb|EGI70210.1| Heat shock 70 kDa protein cognate 3 
## 3: gi|340729370|ref|XP_003402977.1| PREDICTED: heat shock protein beta-1-like isoform 2 
## 4:                      gi|332021988|gb|EGI62314.1| Heat shock 70 kDa protein cognate 4 
## 5:                      gi|307211659|gb|EFN87680.1| Heat shock 70 kDa protein cognate 5 
## 6:                                    gi|227018528|gb|ACP18866.1| heat shock protein 30 
## 7:                      gi|332020093|gb|EGI60539.1| Heat shock factor-binding protein 1 
##    A22.type      Ar.type
## 1:  Bimodal         High
## 2:      Low Intermediate
## 3:  Bimodal      Bimodal
## 4:     High         High
## 5:      Low Intermediate
## 6:     High      NotResp
## 7:      Low Intermediate
```

```r
hsp_all <- annotation.table[grep("shock", annotation.table$best.hit.to.nr), list(Sequence.Name, best.hit.to.nr)]
nrow(hsp_all)
```

```
## [1] 122
```

There are also 7 heat shock related genes in the responsive transcripts, out of 122 total heat shock related genes in the transcriptome annotation.


### Gene set enrichment analysis

I use [topGO](http://www.bioconductor.org/packages/2.12/bioc/html/topGO.html) to perform gene set enrichment analysis (GSEA) seperately for each expression type (bimodal, intermediate, high, low).

First need to create gene ID to GO term map file


```r
# create geneid2go.map file from FastAnnotator annotationtable.txt
mapname = "results/geneid2go.map"
if(file.exists(mapname)) print("Gene ID to GO term file exists") else{
    geneid2GOmap(as.data.frame(annotation.table), ontology = c("BP", "CC", "MF"), mapname = mapname)
}
```

```
## [1] "Gene ID to GO term file exists"
```

then read map file.


```r
# read mappings file
geneID2GO <- readMappings(file = mapname)
str(head(geneID2GO))
```

```
## List of 6
##  $ 0|*|Contig6267              : chr [1:6] "GO:0035335" "GO:0000188" "GO:0006570" "GO:0017017" ...
##  $ 100000|*|comp2663136_c0_seq1: chr ""
##  $ 100001|*|comp3439067_c0_seq1: chr ""
##  $ 100002|*|comp2050457_c0_seq1: chr [1:6] "GO:0006508" "GO:0006412" "GO:0042254" "GO:0005840" ...
##  $ 100004|*|comp131141_c1_seq1 : chr ""
##  $ 100005|*|comp3923213_c0_seq1: chr ""
```

### GSEA for thermally-responsive transcripts

Using this gene2GO map file, perform GSEA for:

**1) all responsive transcripts**

Use `selectFDR` function to select transcripts with adjusted P < 0.05.

Perform GSEA for all thermally-responsive transcripts.

For many of the responsive categories, it appears that there is considerable redundancy in the enriched GO terms. I use the [GOSemSim](http://www.bioconductor.org/packages/release/bioc/html/GOSemSim.html) package to determine semantic similarity of enriched GO terms (Yu, Li, Qin, Bo, Wu, and Wang, 2010) and then perform hierarhichal clustering based on this information distance.


```r
# create geneList. note that NA values cause problems with topGO
# so set any NA to 1 as need to retain for GO analysis
Ap.geneList <- RxNpval$padj
Ap.geneList[which(is.na(Ap.geneList))] <- 1
stopifnot(length(which(is.na(Ap.geneList))) == 0)
names(Ap.geneList) <- RxNpval$Transcript
str(Ap.geneList)

# Function to select genes with adjusted P-value < 0.05
selectFDR <- function(padj) {
    return(padj < 0.05)
}

# create topGOdata object
Ap.BP.GOdata <- new("topGOdata",
                 description = "BP gene set analysis", ontology = "BP",
                 allGenes = Ap.geneList, 
                 geneSel = selectFDR,
                 nodeSize = 10,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Ap.BP.GOdata

# perform enrichment analysis using parentchild method
Ap.BP.resultParentChild <- runTest(Ap.BP.GOdata, statistic = 'fisher', algorithm = 'parentchild')
Ap.BP.resultParentChild

# table results
Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, topNodes = 131)
Ap.BP.ResTable
```

As the molecular processes involved in cold and hot temperature tolerance are different, I perform GSEA separately for each response type and compile results in a single table.


```r
## Genes with *High* expression in both colonies
Ap.high <- Ap.response.type[which(Ap.response.type$Ar.type == "High" & Ap.response.type$A22.type == "High"), "Transcript"]
Ap.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.high) <- names(Ap.geneList)
Ap.geneList.high[(which(names(Ap.geneList.high) %in% Ap.high))] <- 1
# check correct number of values set to 1
table(Ap.geneList.high)
# run GSEA
Ap.high.gsea <- gsea(genelist = Ap.geneList.high, geneID2GO = geneID2GO, plotpath = NA)
```


```r
# similarity among enriched GO terms
Ap.high.GO.term.sim <- termSim(Ap.high.gsea$GO.ID, Ap.high.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ap.high.GO.term.sim) <- Ap.high.gsea$Term
colnames(Ap.high.GO.term.sim) <- Ap.high.gsea$Term

# distance matrix
Ap.high.hclust <- hclust(dist(Ap.high.GO.term.sim))
plot(Ap.high.hclust)
```

![plot of chunk gsea_Ap_high_cluster](figure/gsea_Ap_high_cluster.png) 

```r
# report items in each cluster
Ap.high.hclust.report <- GSEAReportClusters(Ap.high.hclust, h = 1)

# add cluster information
Ap.high.gsea$Cluster <- NA

for(i in 1:length(Ap.high.hclust.report)) {
    Ap.high.gsea[Ap.high.hclust.report[[i]], "Cluster"] <- i
}
```


```r
## Genes with *Low* expression in both colonies
Ap.low <- Ap.response.type[which(Ap.response.type$Ar.type == "Low" & Ap.response.type$A22.type == "Low"), "Transcript"]
Ap.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.low) <- names(Ap.geneList)
Ap.geneList.low[(which(names(Ap.geneList.low) %in% Ap.low))] <- 1
# check correct number of values set to 1
table(Ap.geneList.low)
# Run GSEA
Ap.low.gsea <- gsea(genelist = Ap.geneList.low, geneID2GO = geneID2GO)
```


```r
# similarity among enriched GO terms
Ap.low.GO.term.sim <- termSim(Ap.low.gsea$GO.ID, Ap.low.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ap.low.GO.term.sim) <- Ap.low.gsea$Term
colnames(Ap.low.GO.term.sim) <- Ap.low.gsea$Term

# distance matrix
Ap.low.hclust <- hclust(dist(Ap.low.GO.term.sim))
plot(Ap.low.hclust)
```

![plot of chunk gsea_Ap_low_cluster](figure/gsea_Ap_low_cluster.png) 

```r
# report items in each cluster
Ap.low.gsea.cluster <- Ap.low.hclust.report <- GSEAReportClusters(Ap.low.hclust, h = 1.5)

# add cluster information
Ap.low.gsea$Cluster <- NA

for(i in 1:length(Ap.low.hclust.report)) {
    Ap.low.gsea[Ap.low.hclust.report[[i]], "Cluster"] <- i
}
```



```r
# Genes with *Bimodal* expression in both colonies
Ap.bim <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal" & Ap.response.type$Ar.type == "Bimodal"), "Transcript"]
# create gene list, setting value to 1 for "bim" transcripts
Ap.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.bim) <- names(Ap.geneList)
Ap.geneList.bim[(which(names(Ap.geneList.bim) %in% Ap.bim))] <- 1
# check correct number of values set to 1
table(Ap.geneList.bim)
```

```
## Ap.geneList.bim
##     0     1 
## 97889   297
```

```r
# Run GSEA
Ap.bim.gsea <- gsea(genelist = Ap.geneList.bim, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5452 GO terms found. )
## 
## Build GO DAG topology ..........	( 8935 GO terms and 19895 relations. )
## 
## Annotating nodes ...............	( 30467 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 645 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 16:	1 nodes to be scored.
## 
## 	 Level 15:	2 nodes to be scored.
## 
## 	 Level 14:	2 nodes to be scored.
## 
## 	 Level 13:	4 nodes to be scored.
## 
## 	 Level 12:	8 nodes to be scored.
## 
## 	 Level 11:	19 nodes to be scored.
## 
## 	 Level 10:	49 nodes to be scored.
## 
## 	 Level 9:	71 nodes to be scored.
## 
## 	 Level 8:	79 nodes to be scored.
## 
## 	 Level 7:	84 nodes to be scored.
## 
## 	 Level 6:	96 nodes to be scored.
## 
## 	 Level 5:	108 nodes to be scored.
## 
## 	 Level 4:	76 nodes to be scored.
## 
## 	 Level 3:	30 nodes to be scored.
## 
## 	 Level 2:	15 nodes to be scored.
```


```r
# similarity among enriched GO terms
Ap.bim.GO.term.sim <- termSim(Ap.bim.gsea$GO.ID, Ap.bim.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ap.bim.GO.term.sim) <- Ap.bim.gsea$Term
colnames(Ap.bim.GO.term.sim) <- Ap.bim.gsea$Term

# distance matrix
Ap.bim.hclust <- hclust(dist(Ap.bim.GO.term.sim))
plot(Ap.bim.hclust)
```

![plot of chunk gsea_Ap_bim_cluster](figure/gsea_Ap_bim_cluster.png) 

```r
# report items in each cluster
Ap.bim.hclust.report <- GSEAReportClusters(Ap.bim.hclust, h = 1.5)

# add cluster information
Ap.bim.gsea$Cluster <- NA

for(i in 1:length(Ap.bim.hclust.report)) {
    Ap.bim.gsea[Ap.bim.hclust.report[[i]], "Cluster"] <- i
}
```


```r
# genes with *Intermediate* expression in both colonies
Ap.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & Ap.response.type$Ar.type == "Intermediate"), "Transcript"]
# create gene list, setting value to 1 for "int" transcripts
Ap.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.int) <- names(Ap.geneList)
Ap.geneList.int[(which(names(Ap.geneList.int) %in% Ap.int))] <- 1
# check correct number of values set to 1
table(Ap.geneList.int)
# run GSEA
Ap.int.gsea <- gsea(genelist = Ap.geneList.int, geneID2GO = geneID2GO)
```


```r
# similarity among enriched GO terms
Ap.int.GO.term.sim <- termSim(Ap.int.gsea$GO.ID, Ap.int.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ap.int.GO.term.sim) <- Ap.int.gsea$Term
colnames(Ap.int.GO.term.sim) <- Ap.int.gsea$Term

# distance matrix
Ap.int.hclust <- hclust(dist(Ap.int.GO.term.sim))
plot(Ap.int.hclust)
```

![plot of chunk gsea_Ap_int_cluster](figure/gsea_Ap_int_cluster.png) 

```r
# report items in each cluster
(Ap.int.hclust.report <- GSEAReportClusters(Ap.int.hclust, h = 1.5))
```

```
## $`Cluster 1`
##                glycolipid metabolic process 
##                                           1 
##         membrane lipid biosynthetic process 
##                                           2 
##            membrane lipid metabolic process 
##                                           6 
##             glycolipid biosynthetic process 
##                                           7 
##                 macromolecule glycosylation 
##                                          11 
##              sphingolipid metabolic process 
##                                          12 
##              carbohydrate metabolic process 
##                                          14 
## single-organism carbohydrate metabolic p... 
##                                          15 
##                               glycosylation 
##                                          18 
## 
## $`Cluster 2`
##     glycoprotein metabolic process  glycoprotein biosynthetic process 
##                                  3                                  4 
##              protein glycosylation            DNA conformation change 
##                                  5                                  8 
##                      DNA packaging cellular protein metabolic process 
##                                  9                                 16 
##   lipoprotein biosynthetic process 
##                                 17 
## 
## $`Cluster 3`
##                  nucleosome organization 
##                                       10 
## protein-DNA complex subunit organization 
##                                       13
```

```r
# add cluster information
Ap.int.gsea$Cluster <- NA

for(i in 1:length(Ap.int.hclust.report)) {
    Ap.int.gsea[Ap.int.hclust.report[[i]], "Cluster"] <- i
}
```



GSEA for genes that are *Intermediate* in A. carolinensis and *Bimodal* in A. picea


```r
# genes with *Intermediate* expression in Ar and *Bimodal* in A22
Ar.int.A22.bim <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate" & Ap.response.type$A22.type == "Bimodal"), "Transcript"]
Ar.int.A22.bim.geneList <- rep(0, length = length(Ap.geneList))
names(Ar.int.A22.bim.geneList) <- names(Ap.geneList)
Ar.int.A22.bim.geneList[(which(names(Ar.int.A22.bim.geneList) %in% Ar.int.A22.bim))] <- 1
# check correct number of values set to 1
table(Ar.int.A22.bim.geneList)
# run GSEA
Ar.int.A22.bim.gsea <- gsea(genelist = Ar.int.A22.bim.geneList, geneID2GO = geneID2GO, plotpath = NA)
```

Cluster by semantic similarity


```r
# similarity among enriched GO terms
Ar.int.A22.bim.gsea.term.sim <- termSim(Ar.int.A22.bim.gsea$GO.ID, Ar.int.A22.bim.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ar.int.A22.bim.gsea.term.sim) <- Ar.int.A22.bim.gsea$Term
colnames(Ar.int.A22.bim.gsea.term.sim) <- Ar.int.A22.bim.gsea$Term
# distance matrix
Ar.int.A22.bim.gsea.term.sim.hclust <- hclust(dist(Ar.int.A22.bim.gsea.term.sim))
plot(Ar.int.A22.bim.gsea.term.sim.hclust)
```

![plot of chunk gsea_Ar_int_A22_bim_cluster](figure/gsea_Ar_int_A22_bim_cluster.png) 

```r
# report items in each cluster at height 1.5
(Ar.int.A22.bim.hclust.report <- GSEAReportClusters(Ar.int.A22.bim.gsea.term.sim.hclust, h = 1))
```

```
## $`Cluster 1`
## phospholipid catabolic process           histone modification 
##                              1                              7 
## 
## $`Cluster 2`
## L-amino acid transport 
##                      2 
## 
## $`Cluster 3`
## macromolecule metabolic process 
##                               3 
## 
## $`Cluster 4`
## response to organic substance 
##                             4 
## 
## $`Cluster 5`
## glutamate receptor signaling pathway 
##                                    5 
## 
## $`Cluster 6`
## regulation of nitrogen compound metaboli... 
##                                           6
```

GSEA for genes that are *Intermediate* in A. carolinensis and *Low* in A. picea


```r
# genes with *Intermediate* expression in Ar and *Low* in A22
Ar.int.A22.low <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate" & Ap.response.type$A22.type == "Low"), "Transcript"]
Ar.int.A22.low.geneList <- rep(0, length = length(Ap.geneList))
names(Ar.int.A22.low.geneList) <- names(Ap.geneList)
Ar.int.A22.low.geneList[(which(names(Ar.int.A22.low.geneList) %in% Ar.int.A22.low))] <- 1
# check correct number of values set to 1
table(Ar.int.A22.low.geneList)
# run GSEA
Ar.int.A22.low.gsea <- gsea(genelist = Ar.int.A22.low.geneList, geneID2GO = geneID2GO, plotpath = NA)
```


```r
# similarity among enriched GO terms
Ar.int.A22.low.gsea.term.sim <- termSim(Ar.int.A22.low.gsea$GO.ID, Ar.int.A22.low.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ar.int.A22.low.gsea.term.sim) <- Ar.int.A22.low.gsea$Term
colnames(Ar.int.A22.low.gsea.term.sim) <- Ar.int.A22.low.gsea$Term
# distance matrix
Ar.int.A22.low.gsea.term.sim.hclust <- hclust(dist(Ar.int.A22.low.gsea.term.sim))
plot(Ar.int.A22.low.gsea.term.sim.hclust)
```

![plot of chunk gsea_Ar_int_A22_low_cluster](figure/gsea_Ar_int_A22_low_cluster.png) 

```r
# report items in each cluster at height 1
(Ar.int.A22.low.hclust.report <- GSEAReportClusters(Ar.int.A22.low.gsea.term.sim.hclust, h = 1.2))
```

```
## $`Cluster 1`
## cellular protein metabolic process          protein metabolic process 
##                                  1                                  2 
## 
## $`Cluster 2`
## single-organism organelle organization 
##                                      3 
##                 organelle organization 
##                                      5 
## 
## $`Cluster 3`
## membrane docking  vesicle docking 
##                4               14 
## 
## $`Cluster 4`
## respiratory electron transport chain               recombinational repair 
##                                    6                                   11 
##                        glycosylation 
##                                   13 
## 
## $`Cluster 5`
## small GTPase mediated signal transductio... 
##                                           7 
## regulation of actin filament-based proce... 
##                                          23 
## positive regulation of biological proces... 
##                                          34 
## 
## $`Cluster 6`
##            vesicle-mediated transport 
##                                     8 
##                 cellular localization 
##                                     9 
## establishment of localization in cell 
##                                    12 
##               intracellular transport 
##                                    29 
## 
## $`Cluster 7`
## cellular component organization or bioge... 
##                                          10 
##                             membrane fusion 
##                                          19 
##             single-organism membrane fusion 
##                                          26 
## 
## $`Cluster 8`
## ethanolamine-containing compound metabol... 
##                                          15 
##              phospholipid metabolic process 
##                                          16 
##           phospholipid biosynthetic process 
##                                          21 
##    phosphatidylcholine biosynthetic process 
##                                          25 
## 
## $`Cluster 9`
##                       GTP catabolic process 
##                                          17 
##                       GTP metabolic process 
##                                          20 
## guanosine-containing compound catabolic ... 
##                                          22 
## guanosine-containing compound metabolic ... 
##                                          31 
## 
## $`Cluster 10`
##                  macromolecule modification 
##                                          18 
## posttranscriptional regulation of gene e... 
##                                          24 
##             macromolecule metabolic process 
##                                          27 
## 
## $`Cluster 11`
## endoderm development 
##                   28 
## 
## $`Cluster 12`
##             protein glycosylation                       translation 
##                                30                                33 
## glycoprotein biosynthetic process 
##                                35 
## 
## $`Cluster 13`
##      sphingolipid catabolic process sulfur amino acid catabolic process 
##                                  32                                  36
```



GSEA for genes that are *Intermediate* in A. carolinensis and *Low* or *Bimodal* in A. picea


```r
# genes with *Intermediate* expression in Ar and *Low* in A22
Ar.int.A22.lowbim <- c(Ar.int.A22.low, Ar.int.A22.bim)
Ar.int.A22.lowbim.geneList <- rep(0, length = length(Ap.geneList))
names(Ar.int.A22.lowbim.geneList) <- names(Ap.geneList)
Ar.int.A22.lowbim.geneList[(which(names(Ar.int.A22.lowbim.geneList) %in% Ar.int.A22.lowbim))] <- 1
# check correct number of values set to 1
table(Ar.int.A22.lowbim.geneList)
# run GSEA
Ar.int.A22.lowbim.gsea <- gsea(genelist = Ar.int.A22.lowbim.geneList, geneID2GO = geneID2GO, plotpath = NA)
```


```r
# similarity among enriched GO terms
Ar.int.A22.lowbim.gsea.term.sim <- termSim(Ar.int.A22.lowbim.gsea$GO.ID, Ar.int.A22.lowbim.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ar.int.A22.lowbim.gsea.term.sim) <- Ar.int.A22.lowbim.gsea$Term
colnames(Ar.int.A22.lowbim.gsea.term.sim) <- Ar.int.A22.lowbim.gsea$Term
# distance matrix
Ar.int.A22.lowbim.gsea.term.sim.hclust <- hclust(dist(Ar.int.A22.lowbim.gsea.term.sim))
plot(Ar.int.A22.lowbim.gsea.term.sim.hclust)
```

![plot of chunk gsea_Ar_int_A22_lowbim_cluster](figure/gsea_Ar_int_A22_lowbim_cluster.png) 

```r
# report items in each cluster at height 1
(Ar.int.A22.lowbim.hclust.report <- GSEAReportClusters(Ar.int.A22.lowbim.gsea.term.sim.hclust, h = 1.2))
```

```
## $`Cluster 1`
##       cellular protein metabolic process 
##                                        1 
##                protein metabolic process 
##                                        2 
##          macromolecule metabolic process 
##                                        9 
##               macromolecule modification 
##                                       16 
## cellular macromolecule metabolic process 
##                                       38 
## 
## $`Cluster 2`
##            vesicle-mediated transport 
##                                     3 
## establishment of localization in cell 
##                                    18 
##                 cellular localization 
##                                    21 
##               intracellular transport 
##                                    25 
## 
## $`Cluster 3`
## small GTPase mediated signal transductio... 
##                                           4 
## regulation of actin filament-based proce... 
##                                          27 
## regulation of actin cytoskeleton organiz... 
##                                          42 
## 
## $`Cluster 4`
## single-organism organelle organization 
##                                      5 
##                 organelle organization 
##                                      7 
##                   vacuole organization 
##                                     44 
## 
## $`Cluster 5`
##              phospholipid metabolic process 
##                                           6 
## ethanolamine-containing compound metabol... 
##                                          30 
##    phosphatidylcholine biosynthetic process 
##                                          37 
##           phospholipid biosynthetic process 
##                                          41 
## 
## $`Cluster 6`
##                        glycosylation respiratory electron transport chain 
##                                    8                                   12 
##               recombinational repair                 histone modification 
##                                   32                                   39 
## 
## $`Cluster 7`
##                       GTP catabolic process 
##                                          10 
## guanosine-containing compound catabolic ... 
##                                          14 
##                       GTP metabolic process 
##                                          17 
## guanosine-containing compound metabolic ... 
##                                          22 
## 
## $`Cluster 8`
## cellular component organization or bioge... 
##                                          11 
##                             membrane fusion 
##                                          15 
##             single-organism membrane fusion 
##                                          19 
## 
## $`Cluster 9`
## membrane docking  vesicle docking 
##               13               35 
## 
## $`Cluster 10`
## response to organic substance 
##                            20 
## 
## $`Cluster 11`
## regulation of purine nucleotide metaboli... 
##                                          23 
## regulation of purine nucleotide cataboli... 
##                                          26 
## regulation of nucleotide catabolic proce... 
##                                          36 
## regulation of nucleoside metabolic proce... 
##                                          43 
## 
## $`Cluster 12`
##             protein glycosylation glycoprotein biosynthetic process 
##                                24                                28 
## 
## $`Cluster 13`
## L-amino acid transport 
##                     29 
## 
## $`Cluster 14`
##               regulation of mRNA processing 
##                                          31 
## tRNA splicing, via endonucleolytic cleav... 
##                                          34 
## 
## $`Cluster 15`
## positive regulation of hydrolase activit... 
##                                          33 
## 
## $`Cluster 16`
## regulation of protein polymerization 
##                                   40 
## 
## $`Cluster 17`
## regulation of cellular catabolic process 
##                                       45
```


### GSEA for each functional type in each species

GO enrichment for *A. picea* **High** genes.


```r
# A22 'High' genes
A22.geneList.high <- rep(0, length = length(Ap.geneList))
names(A22.geneList.high) <- names(Ap.geneList)
A22.geneList.high[(which(names(A22.geneList.high) %in% A22.high.transcripts))] <- 1
A22.high.gsea <- gsea(genelist = A22.geneList.high, geneID2GO = geneID2GO)
```
Clustering of enriched GO terms for *A. picea* **High** genes.


```r
# similarity among enriched GO terms
A22.high.gsea.term.sim <- termSim(A22.high.gsea$GO.ID, A22.high.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(A22.high.gsea.term.sim) <- A22.high.gsea$Term
colnames(A22.high.gsea.term.sim) <- A22.high.gsea$Term
# distance matrix
A22.high.gsea.term.sim.hclust <- hclust(dist((A22.high.gsea.term.sim)))
plot((A22.high.gsea.term.sim.hclust))
```

![plot of chunk gsea_A22_high_cluster](figure/gsea_A22_high_cluster.png) 

```r
# report items in each cluster at height 1.8
GSEAReportClusters(A22.high.gsea.term.sim.hclust, h = 1.2)
```

```
## $`Cluster 1`
##                                 translation 
##                                           1 
##                             DNA integration 
##                                           3 
## cellular macromolecule biosynthetic proc... 
##                                           8 
##                glutamyl-tRNA aminoacylation 
##                                           9 
## 
## $`Cluster 2`
## erythrose 4-phosphate/phosphoenolpyruvat... 
##                                           2 
##           L-phenylalanine catabolic process 
##                                           4 
## 
## $`Cluster 3`
##    macromolecule metabolic process macromolecule biosynthetic process 
##                                  5                                  7 
##                    gene expression 
##                                 14 
## 
## $`Cluster 4`
##              regulation of transport trivalent inorganic cation transport 
##                                    6                                   10 
## 
## $`Cluster 5`
## sporulation resulting in formation of a ... 
##                                          11 
##           negative regulation of cell death 
##                                          12 
##                  osteoblast differentiation 
##                                          16 
## 
## $`Cluster 6`
## cellular component biogenesis 
##                            13 
## 
## $`Cluster 7`
## cAMP biosynthetic process 
##                        15
```

GO enrichment for *A. picea* **Low** genes.


```r
# A22 'Low' genes
A22.geneList.low <- rep(0, length = length(Ap.geneList))
names(A22.geneList.low) <- names(Ap.geneList)
A22.geneList.low[(which(names(A22.geneList.low) %in% A22.low.transcripts))] <- 1
A22.low.gsea <- gsea(genelist = A22.geneList.low, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. picea* **Low** genes.


```r
# similarity among enriched GO terms
A22.low.gsea.term.sim <- termSim(A22.low.gsea$GO.ID, A22.low.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(A22.low.gsea.term.sim) <- A22.low.gsea$Term
colnames(A22.low.gsea.term.sim) <- A22.low.gsea$Term
# distance matrix
A22.low.gsea.term.sim.hclust <- hclust(dist((A22.low.gsea.term.sim)))
plot((A22.low.gsea.term.sim.hclust))
```

![plot of chunk gsea_A22_low_cluster](figure/gsea_A22_low_cluster.png) 

```r
# report items in each cluster at height 1.8
GSEAReportClusters(A22.low.gsea.term.sim.hclust, h = 2)
```

```
## $`Cluster 1`
##                protein metabolic process 
##                                        1 
##       cellular protein metabolic process 
##                                        2 
##                            glycosylation 
##                                        3 
##               macromolecule modification 
##                                        6 
##          macromolecule metabolic process 
##                                        7 
## cellular macromolecule metabolic process 
##                                       12 
##        glycoprotein biosynthetic process 
##                                       13 
##                    protein glycosylation 
##                                       14 
##                 translational initiation 
##                                       15 
##             protein modification process 
##                                       16 
##           glycoprotein metabolic process 
##                                       26 
##                     histone modification 
##                                       36 
##    cellular protein modification process 
##                                       69 
##              macromolecule glycosylation 
##                                       79 
##                      histone acetylation 
##                                       82 
##                  ether metabolic process 
##                                       90 
##                        protein acylation 
##                                       91 
## 
## $`Cluster 2`
## small GTPase mediated signal transductio... 
##                                           4 
##           regulation of signal transduction 
##                                          41 
##              regulation of cellular process 
##                                          52 
##                       biological regulation 
##                                          61 
##            regulation of biological process 
##                                          63 
##          regulation of response to stimulus 
##                                          66 
## regulation of Rab protein signal transdu... 
##                                          68 
##             Rab protein signal transduction 
##                                          75 
## regulation of intracellular signal trans... 
##                                          93 
##                     regulation of signaling 
##                                          94 
##            regulation of cell communication 
##                                          97 
## 
## $`Cluster 3`
## single-organism organelle organization 
##                                      5 
##                 organelle organization 
##                                      8 
##                       membrane docking 
##                                     10 
##                        vesicle docking 
##                                     39 
## vesicle docking involved in exocytosis 
##                                     59 
## 
## $`Cluster 4`
## nucleoside triphosphate metabolic proces... 
##                                           9 
##      nucleoside phosphate catabolic process 
##                                          17 
##            ribonucleotide metabolic process 
##                                          22 
##                nucleotide catabolic process 
##                                          24 
##                       GTP metabolic process 
##                                          25 
## regulation of purine nucleotide metaboli... 
##                                          27 
##                       GTP catabolic process 
##                                          28 
## regulation of nucleotide metabolic proce... 
##                                          31 
## regulation of nucleoside metabolic proce... 
##                                          32 
## guanosine-containing compound catabolic ... 
##                                          33 
## guanosine-containing compound metabolic ... 
##                                          34 
## regulation of nucleotide catabolic proce... 
##                                          37 
## regulation of purine nucleotide cataboli... 
##                                          40 
## positive regulation of Rab GTPase activi... 
##                                          46 
##           regulation of Ras GTPase activity 
##                                          53 
## purine-containing compound catabolic pro... 
##                                          73 
##         purine nucleotide metabolic process 
##                                          81 
##                nucleoside catabolic process 
##                                          85 
## 
## $`Cluster 5`
##               intracellular transport 
##                                    11 
##                 cellular localization 
##                                    19 
##            vesicle-mediated transport 
##                                    20 
## establishment of localization in cell 
##                                    23 
##            macromolecule localization 
##                                    45 
## establishment of protein localization 
##                                    57 
##                    chloride transport 
##                                    78 
##         regulation of anion transport 
##                                    86 
##                     protein transport 
##                                    87 
## 
## $`Cluster 6`
##                glycolipid metabolic process 
##                                          18 
##            membrane lipid metabolic process 
##                                          21 
## serine family amino acid metabolic proce... 
##                                          30 
##              sphingolipid metabolic process 
##                                          44 
##                      recombinational repair 
##                                          47 
## flavin-containing compound metabolic pro... 
##                                          50 
##                riboflavin metabolic process 
##                                          51 
##                  cell cycle DNA replication 
##                                          62 
##     hexachlorocyclohexane metabolic process 
##                                          64 
##                   styrene metabolic process 
##                                          65 
##           phospholipid biosynthetic process 
##                                          80 
##              phospholipid metabolic process 
##                                          84 
##              sphingolipid catabolic process 
##                                          96 
##                       amino acid activation 
##                                          98 
## 
## $`Cluster 7`
## regulation of phosphate metabolic proces... 
##                                          29 
##      regulation of translational initiation 
##                                          35 
## posttranscriptional regulation of gene e... 
##                                          38 
##             regulation of catabolic process 
##                                          42 
##    regulation of cellular catabolic process 
##                                          48 
## regulation of phosphorus metabolic proce... 
##                                          49 
##                   regulation of translation 
##                                          77 
## regulation of cellular protein metabolic... 
##                                          92 
## 
## $`Cluster 8`
##                  protein complex biogenesis 
##                                          43 
## positive regulation of protein polymeriz... 
##                                          67 
##        regulation of protein polymerization 
##                                          71 
##                    protein complex assembly 
##                                          72 
##             macromolecular complex assembly 
##                                          83 
## 
## $`Cluster 9`
##                     muscle cell development 
##                                          54 
##                regulation of blood pressure 
##                                          55 
##                   hair cell differentiation 
##                                          56 
##                           megagametogenesis 
##                                          58 
##                      embryo sac development 
##                                          60 
##         embryo sac egg cell differentiation 
##                                          70 
##    regulation of unidimensional cell growth 
##                                          74 
##              detection of external stimulus 
##                                          76 
## imaginal disc-derived wing hair organiza... 
##                                          99 
## 
## $`Cluster 10`
## chondroitin sulfate proteoglycan biosynt... 
##                                          88 
##    chondroitin sulfate biosynthetic process 
##                                          89 
##             folic acid biosynthetic process 
##                                          95
```

GO enrichment for *A. picea* **Bimodal** genes.


```r
# A22 'Bim' genes
A22.geneList.bim <- rep(0, length = length(Ap.geneList))
names(A22.geneList.bim) <- names(Ap.geneList)
A22.geneList.bim[(which(names(A22.geneList.bim) %in% A22.bim.transcripts))] <- 1
A22.bim.gsea <- gsea(genelist = A22.geneList.bim, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. picea* **Bimodal** genes.


```r
# similarity among enriched GO terms
A22.bim.gsea.term.sim <- termSim(A22.bim.gsea$GO.ID, A22.bim.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(A22.bim.gsea.term.sim) <- A22.bim.gsea$Term
colnames(A22.bim.gsea.term.sim) <- A22.bim.gsea$Term
# distance matrix
A22.bim.gsea.term.sim.hclust <- hclust(dist((A22.bim.gsea.term.sim)))
plot((A22.bim.gsea.term.sim.hclust))
```

![plot of chunk gsea_A22_bim_cluster](figure/gsea_A22_bim_cluster.png) 

```r
# report items in each cluster at height 1.8
GSEAReportClusters(A22.bim.gsea.term.sim.hclust, h = 1)
```

```
## $`Cluster 1`
##                   glycine catabolic process 
##                                           1 
## serine family amino acid catabolic proce... 
##                                          12 
## 
## $`Cluster 2`
## protein modification process   macromolecule modification 
##                            2                            4 
## 
## $`Cluster 3`
##               protein phosphorylation 
##                                     3 
## cellular protein modification process 
##                                     5 
##                  histone modification 
##                                    13 
## 
## $`Cluster 4`
## regulation of intracellular signal trans... 
##                                           6 
## protein kinase C-activating G-protein co... 
##                                           7 
##         I-kappaB kinase/NF-kappaB signaling 
##                                           9 
##        glutamate receptor signaling pathway 
##                                          10 
## 
## $`Cluster 5`
## porphyrin-containing compound biosynthet... 
##                                           8 
## 
## $`Cluster 6`
## phospholipid catabolic process 
##                             11
```

GO enrichment for *A. picea* **Intermediate** genes.


```r
# A22 'Int' genes
A22.geneList.int <- rep(0, length = length(Ap.geneList))
names(A22.geneList.int) <- names(Ap.geneList)
A22.geneList.int[(which(names(A22.geneList.int) %in% A22.int.transcripts))] <- 1
A22.int.gsea <- gsea(genelist = A22.geneList.int, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. picea* **Intermediate** genes.


```r
# similarity among enriched GO terms
A22.int.gsea.term.sim <- termSim(A22.int.gsea$GO.ID, A22.int.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(A22.int.gsea.term.sim) <- A22.int.gsea$Term
colnames(A22.int.gsea.term.sim) <- A22.int.gsea$Term
# distance matrix
A22.int.gsea.term.sim.hclust <- hclust(dist((A22.int.gsea.term.sim)))
plot((A22.int.gsea.term.sim.hclust))
```

![plot of chunk gsea_A22_int_cluster](figure/gsea_A22_int_cluster.png) 

```r
# report items in each cluster at height 1.8
GSEAReportClusters(A22.int.gsea.term.sim.hclust, h = 1.3)
```

```
## $`Cluster 1`
##        glycolipid metabolic process membrane lipid biosynthetic process 
##                                   1                                   2 
##    membrane lipid metabolic process     glycolipid biosynthetic process 
##                                   3                                   4 
##      sphingolipid metabolic process 
##                                   8 
## 
## $`Cluster 2`
##    glycoprotein metabolic process             protein glycosylation 
##                                 5                                 6 
##       macromolecule glycosylation glycoprotein biosynthetic process 
##                                 7                                 9 
##                     DNA packaging 
##                                13 
## 
## $`Cluster 3`
##    nitric oxide metabolic process nitric oxide biosynthetic process 
##                                10                                11 
## 
## $`Cluster 4`
## nucleosome organization 
##                      12
```

GO enrichment for *A. carolinensis* **High** genes.


```r
# Ar 'High' genes
Ar.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.high) <- names(Ap.geneList)
Ar.geneList.high[(which(names(Ar.geneList.high) %in% Ar.high.transcripts))] <- 1
Ar.high.gsea <- gsea(genelist = Ar.geneList.high, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. carolinensis* **High** genes.


```r
# similarity among enriched GO terms
Ar.high.gsea.term.sim <- termSim(Ar.high.gsea$GO.ID, Ar.high.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ar.high.gsea.term.sim) <- Ar.high.gsea$Term
colnames(Ar.high.gsea.term.sim) <- Ar.high.gsea$Term
# distance matrix
Ar.high.gsea.term.sim.hclust <- hclust(dist((Ar.high.gsea.term.sim)))
plot((Ar.high.gsea.term.sim.hclust))
```

![plot of chunk gsea_Ar_high_cluster](figure/gsea_Ar_high_cluster.png) 

```r
# report items in each cluster
GSEAReportClusters(Ar.high.gsea.term.sim.hclust, h = 1.2)
```

```
## $`Cluster 1`
##              ossification mammary gland development 
##                         1                         6 
##      pancreas development         tissue remodeling 
##                         7                         8 
## 
## $`Cluster 2`
## microtubule anchoring 
##                     2 
## 
## $`Cluster 3`
## macromolecular complex assembly 
##                               3 
## 
## $`Cluster 4`
## regulation of blood pressure            hormone transport 
##                            4                           11 
## 
## $`Cluster 5`
## negative regulation of transmembrane rec... 
##                                           5 
##        cellular response to biotic stimulus 
##                                          12 
## 
## $`Cluster 6`
##          osteoblast differentiation developmental programmed cell death 
##                                   9                                  10
```
GO enrichment for *A. carolinensis* **Low** genes.


```r
# Ar 'Low' genes
Ar.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.low) <- names(Ap.geneList)
Ar.geneList.low[(which(names(Ar.geneList.low) %in% Ar.low.transcripts))] <- 1
Ar.low.gsea <- gsea(genelist = Ar.geneList.low, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. carolinensis* **Low** genes.


```r
# similarity among enriched GO terms
Ar.low.gsea.term.sim <- termSim(Ar.low.gsea$GO.ID, Ar.low.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ar.low.gsea.term.sim) <- Ar.low.gsea$Term
colnames(Ar.low.gsea.term.sim) <- Ar.low.gsea$Term
# distance matrix
Ar.low.gsea.term.sim.hclust <- hclust(dist((Ar.low.gsea.term.sim)))
plot((Ar.low.gsea.term.sim.hclust))
```

![plot of chunk gsea_Ar_low_cluster](figure/gsea_Ar_low_cluster.png) 

```r
# report items in each cluster at height 1.8
GSEAReportClusters(Ar.low.gsea.term.sim.hclust, h = 2)
```

```
## $`Cluster 1`
##                protein metabolic process 
##                                        1 
##                 translational initiation 
##                                        7 
##               macromolecule modification 
##                                       13 
##       cellular protein metabolic process 
##                                       14 
##                        protein acylation 
##                                       18 
##                    protein glycosylation 
##                                       19 
##          macromolecule metabolic process 
##                                       33 
##        glycoprotein biosynthetic process 
##                                       39 
##           glycoprotein metabolic process 
##                                       43 
##             protein modification process 
##                                       59 
## cellular macromolecule metabolic process 
##                                       65 
## 
## $`Cluster 2`
##                               glycosylation 
##                                           2 
##                   styrene metabolic process 
##                                          10 
##           organophosphate catabolic process 
##                                          23 
##            membrane lipid metabolic process 
##                                          25 
##         glycosyl compound catabolic process 
##                                          29 
##                         histone acetylation 
##                                          31 
##             folic acid biosynthetic process 
##                                          44 
## chondroitin sulfate proteoglycan biosynt... 
##                                          61 
##    chondroitin sulfate biosynthetic process 
##                                          62 
##    DNA replication, synthesis of RNA primer 
##                                          63 
## folic acid-containing compound metabolic... 
##                                          64 
##                    malate metabolic process 
##                                          67 
## 
## $`Cluster 3`
## nucleoside triphosphate metabolic proces... 
##                                           3 
##      nucleoside phosphate catabolic process 
##                                           4 
##            ribonucleotide metabolic process 
##                                           5 
##                nucleotide catabolic process 
##                                          11 
## regulation of nucleotide metabolic proce... 
##                                          22 
##                nucleoside catabolic process 
##                                          27 
## purine-containing compound catabolic pro... 
##                                          36 
##                       GTP metabolic process 
##                                          38 
##         purine nucleotide metabolic process 
##                                          41 
## regulation of nucleoside metabolic proce... 
##                                          42 
##           regulation of Ras GTPase activity 
##                                          51 
## guanosine-containing compound metabolic ... 
##                                          53 
## positive regulation of Rab GTPase activi... 
##                                          55 
## regulation of purine nucleotide metaboli... 
##                                          66 
## 
## $`Cluster 4`
## small GTPase mediated signal transductio... 
##                                           6 
##                 regulation of ion transport 
##                                           8 
##           regulation of signal transduction 
##                                          16 
##     regulation of muscle tissue development 
##                                          32 
##               regulation of anion transport 
##                                          50 
##                     regulation of signaling 
##                                          68 
## 
## $`Cluster 5`
## establishment of protein localization 
##                                     9 
##                     protein transport 
##                                    20 
##            macromolecule localization 
##                                    21 
##              protein complex assembly 
##                                    40 
##            protein complex biogenesis 
##                                    45 
##               protein oligomerization 
##                                    46 
##        detection of external stimulus 
##                                    49 
##               intracellular transport 
##                                    57 
## 
## $`Cluster 6`
##                       muscle system process 
##                                          12 
##                          limb morphogenesis 
##                                          17 
##                    renal tubule development 
##                                          24 
##                     muscle cell development 
##                                          26 
##                            limb development 
##                                          28 
##                 skeletal system development 
##                                          34 
##        multicellular organismal homeostasis 
##                                          47 
##                          photomorphogenesis 
##                                          52 
##                          embryo development 
##                                          54 
## sensory perception of mechanical stimulu... 
##                                          56 
## 
## $`Cluster 7`
## regulation of phosphate metabolic proces... 
##                                          15 
## regulation of phosphorus metabolic proce... 
##                                          30 
## posttranscriptional regulation of gene e... 
##                                          35 
##      regulation of translational initiation 
##                                          37 
##             regulation of catabolic process 
##                                          48 
## regulation of cellular protein metabolic... 
##                                          58 
##    regulation of cellular catabolic process 
##                                          60
```

GO enrichment for *A. carolinensis* **Bimodal** genes.


```r
# Ar 'Bim' genes
Ar.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.bim) <- names(Ap.geneList)
Ar.geneList.bim[(which(names(Ar.geneList.bim) %in% Ar.bim.transcripts))] <- 1
Ar.bim.gsea <- gsea(genelist = Ar.geneList.bim, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. carolinensis* **Bimodal** genes.


```r
# similarity among enriched GO terms
Ar.bim.gsea.term.sim <- termSim(Ar.bim.gsea$GO.ID, Ar.bim.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ar.bim.gsea.term.sim) <- Ar.bim.gsea$Term
colnames(Ar.bim.gsea.term.sim) <- Ar.bim.gsea$Term
# distance matrix
Ar.bim.gsea.term.sim.hclust <- hclust(dist((Ar.bim.gsea.term.sim)))
plot((Ar.bim.gsea.term.sim.hclust))
```

![plot of chunk gsea_Ar_bim_cluster](figure/gsea_Ar_bim_cluster.png) 

```r
# report items in each cluster at height 1.8
GSEAReportClusters(Ar.bim.gsea.term.sim.hclust, h = 1)
```

```
## $`Cluster 1`
## hemolymph coagulation 
##                     1 
## 
## $`Cluster 2`
## outflow tract morphogenesis 
##                           2 
## 
## $`Cluster 3`
## cellular aromatic compound metabolic pro... 
##                                           3 
##    cellular macromolecule metabolic process 
##                                           6 
## 
## $`Cluster 4`
## pyrimidine nucleoside monophosphate bios... 
##                                           4 
## nucleobase-containing compound metabolic... 
##                                           5 
## 
## $`Cluster 5`
##      cellular potassium ion transport 
##                                     7 
## potassium ion transmembrane transport 
##                                     9 
## 
## $`Cluster 6`
## organic cyclic compound metabolic proces... 
##                                           8 
##             macromolecule metabolic process 
##                                          10
```

GO enrichment for *A. carolinensis* **Intermediate** genes.


```r
# Ar 'Int' genes
Ar.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.int) <- names(Ap.geneList)
Ar.geneList.int[(which(names(Ar.geneList.int) %in% Ar.int.transcripts))] <- 1
Ar.int.gsea <- gsea(genelist = Ar.geneList.int, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. carolinensis* **Intermediate** genes.


```r
# similarity among enriched GO terms
Ar.int.gsea.term.sim <- termSim(Ar.int.gsea$GO.ID, Ar.int.gsea$GO.ID, ont = "BP", organism = 'fly')
# replace GO IDs with terms
rownames(Ar.int.gsea.term.sim) <- Ar.int.gsea$Term
colnames(Ar.int.gsea.term.sim) <- Ar.int.gsea$Term
# distance matrix
Ar.int.gsea.term.sim.hclust <- hclust(dist((Ar.int.gsea.term.sim)))
plot((Ar.int.gsea.term.sim.hclust))
```

![plot of chunk gsea_Ar_int_cluster](figure/gsea_Ar_int_cluster.png) 

```r
# report items in each cluster
GSEAReportClusters(Ar.int.gsea.term.sim.hclust, h = 1.8)
```

```
## $`Cluster 1`
##       cellular protein metabolic process 
##                                        1 
##                protein metabolic process 
##                                        2 
##                            glycosylation 
##                                        3 
##                    protein glycosylation 
##                                        4 
##        glycoprotein biosynthetic process 
##                                        6 
##          macromolecule metabolic process 
##                                        7 
##               macromolecule modification 
##                                       10 
##           glycoprotein metabolic process 
##                                       11 
##              macromolecule glycosylation 
##                                       19 
## cellular macromolecule metabolic process 
##                                       24 
##             protein modification process 
##                                       40 
##   regulation of translational initiation 
##                                       49 
## 
## $`Cluster 2`
##         glycolipid metabolic process respiratory electron transport chain 
##                                    5                                   15 
##       phospholipid metabolic process             electron transport chain 
##                                   17                                   30 
##       sphingolipid metabolic process     membrane lipid metabolic process 
##                                   32                                   34 
##                     membrane docking      glycolipid biosynthetic process 
##                                   35                                   37 
##  membrane lipid biosynthetic process              lipid metabolic process 
##                                   47                                   51 
## 
## $`Cluster 3`
##                       GTP catabolic process 
##                                           8 
## guanosine-containing compound catabolic ... 
##                                          12 
##                       GTP metabolic process 
##                                          16 
## guanosine-containing compound metabolic ... 
##                                          23 
## regulation of purine nucleotide metaboli... 
##                                          28 
## regulation of purine nucleotide cataboli... 
##                                          33 
## regulation of nucleotide catabolic proce... 
##                                          42 
##      positive regulation of GTPase activity 
##                                          43 
## 
## $`Cluster 4`
## small GTPase mediated signal transductio... 
##                                           9 
##      single-organism organelle organization 
##                                          13 
##                      organelle organization 
##                                          14 
## regulation of actin cytoskeleton organiz... 
##                                          22 
## regulation of actin filament-based proce... 
##                                          27 
## 
## $`Cluster 5`
##               response to organic substance 
##                                          18 
##                  vesicle-mediated transport 
##                                          20 
## cellular component organization or bioge... 
##                                          21 
##                       cellular localization 
##                                          36 
##                     intracellular transport 
##                                          39 
##       establishment of localization in cell 
##                                          44 
## 
## $`Cluster 6`
##                  nucleosome organization 
##                                       25 
## protein-DNA complex subunit organization 
##                                       29 
##     regulation of protein polymerization 
##                                       48 
## 
## $`Cluster 7`
##                                RNA splicing 
##                                          26 
## tRNA splicing, via endonucleolytic cleav... 
##                                          31 
##                          chromatin assembly 
##                                          38 
##               regulation of mRNA processing 
##                                          41 
##                         nucleosome assembly 
##                                          50 
## 
## $`Cluster 8`
##         glycosaminoglycan catabolic process 
##                                          45 
## positive regulation of hydrolase activit... 
##                                          46
```

Combine gene set enrichment analysis for each species and category into single table.


```r
# combine into single table
A22.high.gsea$Type <- "High"
A22.low.gsea$Type <- "Low"
A22.int.gsea$Type <- "Intermediate"
A22.bim.gsea$Type <- "Bimodal"
A22.gsea <- rbind(A22.high.gsea, A22.low.gsea, A22.int.gsea, A22.bim.gsea)
A22.gsea$Species <- "ApVT"

Ar.high.gsea$Type <- "High"
Ar.low.gsea$Type <- "Low"
Ar.bim.gsea$Type <- "Bimodal"
Ar.int.gsea$Type <- "Intermediate"
Ar.gsea <- rbind(Ar.high.gsea, Ar.low.gsea, Ar.int.gsea, Ar.bim.gsea)
Ar.gsea$Species <- "AcNC"

# combine
Ap.gsea.union <- rbind(A22.gsea, Ar.gsea)
# reorder
Ap.gsea.union <- Ap.gsea.union[,c("Species", "Type", "GO.ID", "Term", "Annotated", "Significant", "Expected", "parentchild")]
colnames(Ap.gsea.union)[8] <- "P"

write.csv(Ap.gsea.union, file = paste(resultsdir, "Ap_gsea_union", Sys.Date(), ".csv", sep = ""), row.names = FALSE)
```



## Visualize responsive transcripts

Plots for all genes expressed at *High* temps in GO category "GO:0006950: response to stress"

![plot of chunk plot_GOstress](figure/plot_GOstress1.png) ![plot of chunk plot_GOstress](figure/plot_GOstress2.png) 

Make plots for all genes expressed at *Low* temps in GO category "GO:0006950: response to stress"

![plot of chunk plot_GOstress_low](figure/plot_GOstress_low1.png) ![plot of chunk plot_GOstress_low](figure/plot_GOstress_low2.png) 



```r
trp_A22_high <- ggplot(resp.TPM.dt.sub[A22.high.transcripts][1:220,], aes(x=val, y=TPM.scaled, group=Transcript)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) + 
  facet_grid(. ~ colony2) + 
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
print(trp_A22_high)
```

![plot of chunk ggplot_high](figure/ggplot_high1.png) 

```r
### predicted values
resp.TPM.dt.sub.pred[,pTPM.scaled:=scale(pTPM), by = Transcript]
```

```
##                          Transcript    TPM  val colony pTPM pTPM.scaled
##      1: 100008|*|comp137625_c0_seq2 0.0000  0.0    A22 1.01      -0.609
##      2: 100008|*|comp137625_c0_seq2 0.0000  3.5    A22 1.01      -0.639
##      3: 100008|*|comp137625_c0_seq2 0.0741 10.5    A22 1.01      -0.662
##      4: 100008|*|comp137625_c0_seq2 0.0000 14.0    A22 1.01      -0.655
##      5: 100008|*|comp137625_c0_seq2 0.0000 17.5    A22 1.01      -0.635
##     ---                                                                
## 199140:      9|*|comp147140_c0_seq1 0.6310 24.5     Ar 1.67      -0.510
## 199141:      9|*|comp147140_c0_seq1 0.5952 28.0     Ar 1.63      -0.643
## 199142:      9|*|comp147140_c0_seq1 0.7302 31.5     Ar 1.59      -0.822
## 199143:      9|*|comp147140_c0_seq1 0.4534 35.0     Ar 1.53      -1.042
## 199144:      9|*|comp147140_c0_seq1 0.4528 38.5     Ar 1.47      -1.300
```

```r
pred_A22_high <- ggplot(resp.TPM.dt.sub.pred[A22.high.transcripts][1:220,], aes(x=val, y=pTPM.scaled, group=Transcript)) +
  geom_line() + 
  facet_grid(. ~ colony) + 
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
print(pred_A22_high)
```

![plot of chunk ggplot_high](figure/ggplot_high2.png) 

```r
png("RxN_example.png")
ggplot(resp.TPM.dt.sub.pred["100348|*|comp3450455_c0_seq1"], aes(x=val, y=pTPM.scaled, group=Transcript, color = Transcript)) +
  geom_line() + 
  facet_grid(. ~ colony) +
  labs(y = "Expression (scaled)", x = expression("Temperature ", degree, "C")) +
  theme(legend.position = "none")
dev.off()
```

```
## pdf 
##   2
```

## Shiny interactive web-app

To assist visualization of specific transcripts, I made a interactive web-app using the [shiny](http://www.rstudio.com/shiny/) package. The scripts for this app are in the sub-directory `.\ApRxN-shinyapp`.

Export data for interactive shiny app. 







## Session information


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
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C         LC_TIME=C           
##  [4] LC_COLLATE=C         LC_MONETARY=C        LC_MESSAGES=C       
##  [7] LC_PAPER=C           LC_NAME=C            LC_ADDRESS=C        
## [10] LC_TELEPHONE=C       LC_MEASUREMENT=C     LC_IDENTIFICATION=C 
## 
## attached base packages:
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] GOSemSim_1.22.0      Rcpp_0.11.2          Rgraphviz_2.8.1     
##  [4] topGO_2.16.0         SparseM_1.05         GO.db_2.14.0        
##  [7] RSQLite_0.11.4       DBI_0.3.0            AnnotationDbi_1.26.0
## [10] GenomeInfoDb_1.0.2   Biobase_2.24.0       BiocGenerics_0.10.0 
## [13] graph_1.42.0         RColorBrewer_1.0-5   reshape2_1.4        
## [16] xtable_1.7-3         MASS_7.3-34          plyr_1.8.1          
## [19] RCurl_1.95-4.3       bitops_1.0-6         data.table_1.9.2    
## [22] stringr_0.6.2        pander_0.3.8         knitcitations_1.0-1 
## [25] ggplot2_1.0.0        R.utils_1.33.0       R.oo_1.18.0         
## [28] R.methodsS3_1.6.1    knitr_1.6           
## 
## loaded via a namespace (and not attached):
##  [1] IRanges_1.22.10   RJSONIO_1.3-0     RefManageR_0.8.34
##  [4] XML_3.98-1.1      bibtex_0.3-6      codetools_0.2-8  
##  [7] colorspace_1.2-4  digest_0.6.4      evaluate_0.5.5   
## [10] formatR_1.0       gtable_0.1.2      httr_0.5         
## [13] labeling_0.3      lattice_0.20-29   lubridate_1.3.3  
## [16] memoise_0.2.1     munsell_0.4.2     proto_0.3-10     
## [19] scales_0.2.4      stats4_3.1.1      tools_3.1.1
```

## References

[1] S. Anders, D. J. McCarthy, Y. Chen, et al. "Count-based
differential expression analysis of RNA sequencing data using R
and Bioconductor". In: _Nature Protocols_ 8.9 (Aug. 2013), pp.
1765-1786. DOI: 10.1038/nprot.2013.099. <URL:
http://dx.doi.org/10.1038/nprot.2013.099>.

[2] J. H. Bullard, E. Purdom, K. D. Hansen, et al. "Evaluation of
statistical methods for normalization and differential expression
in mRNA-Seq experiments". In: _BMC Bioinformatics_ 11.1 (2010), p.
94. DOI: 10.1186/1471-2105-11-94. <URL:
http://dx.doi.org/10.1186/1471-2105-11-94>.

[3] M. G. Grabherr, B. J. Haas, M. Yassour, et al. "Full-length
transcriptome assembly from RNA-Seq data without a reference
genome". In: _Nat Biotechnol_ 29.7 (May. 2011), pp. 644-652. DOI:
10.1038/nbt.1883. <URL: http://dx.doi.org/10.1038/nbt.1883>.

[4] X. Huang. "CAP3: A DNA Sequence Assembly Program". In: _Genome
Research_ 9.9 (Sep. 1999), pp. 868-877. DOI: 10.1101/gr.9.9.868.
<URL: http://dx.doi.org/10.1101/gr.9.9.868>.

[5] B. Li, V. Ruotti, R. M. Stewart, et al. "RNA-Seq gene
expression estimation with read mapping uncertainty". In:
_Bioinformatics_ 26.4 (Dec. 2009), pp. 493-500. DOI:
10.1093/bioinformatics/btp692. <URL:
http://dx.doi.org/10.1093/bioinformatics/btp692>.

[6] M. Lohse, A. M. Bolger, A. Nagel, et al. "RobiNA: a
user-friendly, integrated software solution for RNA-Seq-based
transcriptomics". In: _Nucleic Acids Research_ 40.W1 (Jun. 2012),
pp. W622-W627. DOI: 10.1093/nar/gks540. <URL:
http://dx.doi.org/10.1093/nar/gks540>.

[7] D. Lubertazzi. " The Biology and Natural History of
Aphaenogaster rudis ". In: _Psyche: A Journal of Entomology_ 2012
(2012), pp. 1-11. DOI: 10.1155/2012/752815. <URL:
http://dx.doi.org/10.1155/2012/752815>.

[8] C. J. Murren, H. J. Maclean, S. E. Diamond, et al.
"Evolutionary Change in Continuous Reaction Norms". In: _The
American Naturalist_ 183.4 (Apr. 2014), pp. 453-467. DOI:
10.1086/675302. <URL: http://dx.doi.org/10.1086/675302>.

[9] Y. Yang and S. A. Smith. "Optimizing de novo assembly of
short-read RNA-seq data for phylogenomics". In: _BMC Genomics_
14.1 (2013), p. 328. DOI: 10.1186/1471-2164-14-328. <URL:
http://dx.doi.org/10.1186/1471-2164-14-328>.

[10] G. Yu, F. Li, Y. Qin, et al. "GOSemSim: an R package for
measuring semantic similarity among GO terms and gene products".
In: _Bioinformatics_ 26.7 (Feb. 2010), pp. 976-978. DOI:
10.1093/bioinformatics/btq064. <URL:
http://dx.doi.org/10.1093/bioinformatics/btq064>.
