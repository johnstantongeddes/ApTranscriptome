Thermal reactionome of the common ant species *Aphaenogaster picea* and *A. carolinensis*
========================================================================================
   
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**May 6, 2014**

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

Table: Statistics for Trinity and cap3+uclust reduced transcriptome assemblies (continued below)

 

|    &nbsp;     |  Mean contig size  |  N50 contig  |  N50 Length  |
|:-------------:|:------------------:|:------------:|:------------:|
|  **trinity**  |        795         |    16,201    |    1,631     |
|  **reduced**  |        593         |    15,491    |     895      |


To remove contigs that are likely contaminants from bacterial, archael, virus or human sources, we used the program [DeconSeq](http://deconseq.sourceforge.net/) (<a href="http://dx.doi.org/10.1371/journal.pone.0017288">Schmieder et al. 2011</a>). We downloaded the bacteria, virus, archae and human [databases of contaminants](ftp://edwards.sdsu.edu:7009/deconseq/db), modified the `DeconSeqConfig.pm` file as described [here](http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/install/deconseq.htm) to point to the databases, and ran DeconSeq specifiying 95% identity over 50% the length of contig

    deconseq.pl -c 50 -i 95 -f Trinity_cap3_uclust.fasta -d Trinity_cap3_uclust -dbs hsref,bast,vir,arch
	
This resulted in removing 5,675 contigs as contaminants, leaving 99,861 "clean" contigs. We spot-checked the contaminants by BLAST and confirmed that they matched bacteria, human or viral sources by greater than 95%. For expression quantification, we use the full assembly to ensure that "contaminant" reads are assigned to the contaminants. After quantification, these transcripts will then be removed from further analyses.

Running Trinity and subsequent programs is time and memory-intensive so the final assembly can be downloaded and used for all further analyses. In addition, the compressed archive contains the "clean" and "contaminant" sequences after filtering with DeconSeq.

~~~
# download filtered Trinity assembly, uncompress and move
wget http://johnstantongeddes.org/assets/files/Aphaenogaster_transcriptome.tar
# check md5sum
md5sum Aphaenogaster_transcriptome.tar
# fd6dbb0b3e88e1c797c9e74611b245b2

# uncompress and move
tar -xvf Aphaenogaster_transcriptome.tar
mkdir -p results/
mkdir -p results/trinity-full/
mv Trinity_cap3_uclust.fasta results/trinity-full/.
mv Trinity_cap3_uclust_clean.fasta results/trinity-full/.
~~~


## Transcriptome annotation ##

Annotation was performed by uploading the reduced assembly "Trinity_cap3_uclust.fasta" to the web-based annotation program [FastAnnotator](http://fastannotator.cgu.edu.tw/index.php) (<a href="">unknown, unknown</a>).

Results are available as job ID [13894410176993](http://fastannotator.cgu.edu.tw/job.php?jobid=13894410176993#page=basicinfo).

This annotation file can be read directly to R:


```r
### Annotation file

# from either AWS or GoogleDrive
annotationURL <- getURL("http://johnstantongeddes.org/assets/files/ApTranscriptome_AnnotationTable_20140113.txt")
# a2 <- getURL('https://googledrive.com/host/0B75IymziRJ_9Tlg1U1Vxbjk1bzg') #
# GoogleDrive link

annotationfile <- read.csv(textConnection(annotationURL), header = TRUE, sep = "\t", 
    stringsAsFactors = FALSE)
str(annotationfile)
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
annotationtable <- data.table(annotationfile)
head(annotationtable)
```

```
##             Sequence.Name sequence.length
## 1:         0|*|Contig6267            9990
## 2: 1|*|comp150820_c2_seq6            9944
## 3:         2|*|Contig6262            9711
## 4: 3|*|comp149397_c1_seq2            9639
## 5:         4|*|Contig4755            9558
## 6:         5|*|Contig3727            9436
##                                                                   best.hit.to.nr
## 1:      gi|110756860|ref|XP_392375.3| PREDICTED: hypothetical protein LOC408844 
## 2:                  gi|307181425|gb|EFN69020.1| FH2 domain-containing protein 1 
## 3:      gi|110756860|ref|XP_392375.3| PREDICTED: hypothetical protein LOC408844 
## 4: gi|332024185|gb|EGI64399.1| AT-rich interactive domain-containing protein 5B 
## 5:                  gi|322799153|gb|EFZ20592.1| hypothetical protein SINV_04047 
## 6:                                     gi|307168092|gb|EFN61390.1| Endophilin-A 
##    hit.length   E.value   Bit.score
## 1:        598  2.1e-265  904.618736
## 2:       1777       0.0 3455.802741
## 3:        511 1.29e-246  842.252472
## 4:       1077       0.0 2272.638439
## 5:       1210       0.0 2809.706194
## 6:        372 1.55e-245  838.663047
##                                                                                                                                                                               GO.Biological.Process
## 1:                                                                GO:0035335 peptidyl-tyrosine dephosphorylation | GO:0000188 inactivation of MAPK activity | GO:0006570 tyrosine metabolic process
## 2:                                                                                                                          GO:0030036 actin cytoskeleton organization | GO:0015074 DNA integration
## 3:                                                                GO:0035335 peptidyl-tyrosine dephosphorylation | GO:0000188 inactivation of MAPK activity | GO:0006570 tyrosine metabolic process
## 4:                                                                                                                                                                           GO:0006508 proteolysis
## 5:                                                       GO:0055114 oxidation-reduction process | GO:0006355 regulation of transcription, DNA-dependent | GO:0009395 phospholipid catabolic process
## 6: GO:0007269 neurotransmitter secretion | GO:0050803 regulation of synapse structure and activity | GO:0048488 synaptic vesicle endocytosis | GO:0042967 acyl-carrier-protein biosynthetic process
##                                           GO.Cellular.Component
## 1:                                                            -
## 2:                                                            -
## 3:                                                            -
## 4:                                           GO:0005634 nucleus
## 5: GO:0005634 nucleus | GO:0005667 transcription factor complex
## 6:                   GO:0005737 cytoplasm | GO:0016020 membrane
##                                                                                                                                                                                                                                            GO.Molecular.Function
## 1:                                                                                                          GO:0017017 MAP kinase tyrosine/serine/threonine phosphatase activity | GO:0004725 protein tyrosine phosphatase activity | GO:0016301 kinase activity
## 2:                                                                                                                                      GO:0003676 nucleic acid binding | GO:0003779 actin binding | GO:0008270 zinc ion binding | GO:0017048 Rho GTPase binding
## 3:                                                                                                          GO:0017017 MAP kinase tyrosine/serine/threonine phosphatase activity | GO:0004725 protein tyrosine phosphatase activity | GO:0016301 kinase activity
## 4:                                                                                                                                                                                        GO:0003677 DNA binding | GO:0004252 serine-type endopeptidase activity
## 5: GO:0043565 sequence-specific DNA binding | GO:0003700 sequence-specific DNA binding transcription factor activity | GO:0032440 2-alkenal reductase [NAD(P)] activity | GO:0004623 phospholipase A2 activity | GO:0003950 NAD+ ADP-ribosyltransferase activity
## 6:                                                                                                                                                                                                     GO:0042171 lysophosphatidic acid acyltransferase activity
##                   Enzyme
## 1: 3.1.3.16  | 3.1.3.48 
## 2:                     -
## 3: 3.1.3.16  | 3.1.3.48 
## 4:                     -
## 5:             2.4.2.30 
## 6:                     -
##                                                                                                                      Domain
## 1:                                                                                     pfam00782 DSPc | pfam00581 Rhodanese
## 2:         pfam02181 FH2 | pfam00067 p450 | pfam06367 Drf_FH3 | pfam12795 MscS_porin | pfam01749 IBB | pfam07926 TPR_MLP1_2
## 3:                                                                                     pfam00782 DSPc | pfam00581 Rhodanese
## 4:                                                                                                           pfam01388 ARID
## 5: pfam12796 Ank_2 | pfam00644 PARP | pfam13637 Ank_4 | pfam13857 Ank_5 | pfam00023 Ank | pfam13606 Ank_3 | pfam07647 SAM_2
## 6:                                      pfam03114 BAR | pfam00018 SH3_1 | pfam07653 SH3_2 | pfam10455 BAR_2 | pfam08397 IMD
##         annotation.type
## 1: GO & Enzyme & Domain
## 2:          GO & Domain
## 3: GO & Enzyme & Domain
## 4:          GO & Domain
## 5: GO & Enzyme & Domain
## 6:          GO & Domain
```


Note that I used the reduced assembly, prior to cleaning out contaminants. As I use the "clean" assembly for read mapping and identification of responsive genes, the contaminants are dis-regarded in downstream analyses. 


## Identification of thermally-responsive genes ##

### Quantify gene expression ###

I quantified gene expression using [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html). To run this program, first make sure that PATHs to the software libraries are set up correctly as described on the sailfish website. 
                                                 
Then build the index of the assembly:


```r
system("sailfish index -t results/trinity-full/Trinity_cap3_uclust.fasta -o results/trinity-full/sailfish-index-Trinity-cap3-uclust -k 20 -p 4")
```


Once this is done, quantify expression for the Trimmomatic filtered reads from each colony-treatment sample separately. Note that for each sample, there are three four filtered read files:

- paired.left.fastq
- paired.right.fastq
- unpaired.left.fastq
- unpaired.right.fastq
                                                 
Make a directory for the expression values


```r
system("mkdir -p results/trinity-full/sailfish-expression-Trinity-cap3-uclust")
```


Then, for each sample, run the following command in the `results/trinity-full/` directory:
                                                 
    sailfish -i sailfish-index-Trinity-cap3-uclust -o sailfish-expression-Trinity-cap3-uclust/A22-0 -l "T=PE:O=><:S=U" -r A22-0_ATCACG.unpaired.left.fastq A22-0_ATCACG.unpaired.right.fastq -1 A22-0_ATCACG.paired.left.fastq -2 A22-0_ATCACG.paired.right.fastq  -p 4

Or, using a loop in R:                                                 
                                                 

```r
# directory containing trimmed reads
readdir <- "data/ind_files/"
# sailfish index directory
sfindex <- "results/trinity-full/sailfish-index-Trinity-cap3-uclust"
# sailfish expression directory
sfexpressionroot <- "results/trinity-full/sailfish-expression-Trinity-cap3-uclust/"

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
##  [1] "A22-0_ATCACG.unpaired.left.fastq"  "A22-10_TGACCA.unpaired.left.fastq"
##  [3] "A22-14_ACAGTG.unpaired.left.fastq" "A22-17_GCCAAT.unpaired.left.fastq"
##  [5] "A22-21_CAGATC.unpaired.left.fastq" "A22-24_ACTTGA.unpaired.left.fastq"
##  [7] "A22-28_GATCAG.unpaired.left.fastq" "A22-31_TAGCTT.unpaired.left.fastq"
##  [9] "A22-35_GGCTAC.unpaired.left.fastq" "A22-38_CTTGTA.unpaired.left.fastq"
## [11] "A22-3_CGATGT.unpaired.left.fastq"  "A22-7_TTAGGC.unpaired.left.fastq" 
## [13] "Ar-0_AGTCAA.unpaired.left.fastq"   "Ar-10_CCGTCC.unpaired.left.fastq" 
## [15] "Ar-14_GTCCGC.unpaired.left.fastq"  "Ar-17_GTGAAA.unpaired.left.fastq" 
## [17] "Ar-21_GTGGCC.unpaired.left.fastq"  "Ar-24_GTTTCG.unpaired.left.fastq" 
## [19] "Ar-28_CGTACG.unpaired.left.fastq"  "Ar-31_GAGTGG.unpaired.left.fastq" 
## [21] "Ar-35_ACTGAT.unpaired.left.fastq"  "Ar-38_ATTCCT.unpaired.left.fastq" 
## [23] "Ar-3_AGTTCC.unpaired.left.fastq"   "Ar-7_ATGTCA.unpaired.left.fastq"
```

```r
(unpaired.right <- readlist[grep("unpaired.right.fastq$", readlist)])
```

```
##  [1] "A22-0_ATCACG.unpaired.right.fastq"  "A22-10_TGACCA.unpaired.right.fastq"
##  [3] "A22-14_ACAGTG.unpaired.right.fastq" "A22-17_GCCAAT.unpaired.right.fastq"
##  [5] "A22-21_CAGATC.unpaired.right.fastq" "A22-24_ACTTGA.unpaired.right.fastq"
##  [7] "A22-28_GATCAG.unpaired.right.fastq" "A22-31_TAGCTT.unpaired.right.fastq"
##  [9] "A22-35_GGCTAC.unpaired.right.fastq" "A22-38_CTTGTA.unpaired.right.fastq"
## [11] "A22-3_CGATGT.unpaired.right.fastq"  "A22-7_TTAGGC.unpaired.right.fastq" 
## [13] "Ar-0_AGTCAA.unpaired.right.fastq"   "Ar-10_CCGTCC.unpaired.right.fastq" 
## [15] "Ar-14_GTCCGC.unpaired.right.fastq"  "Ar-17_GTGAAA.unpaired.right.fastq" 
## [17] "Ar-21_GTGGCC.unpaired.right.fastq"  "Ar-24_GTTTCG.unpaired.right.fastq" 
## [19] "Ar-28_CGTACG.unpaired.right.fastq"  "Ar-31_GAGTGG.unpaired.right.fastq" 
## [21] "Ar-35_ACTGAT.unpaired.right.fastq"  "Ar-38_ATTCCT.unpaired.right.fastq" 
## [23] "Ar-3_AGTTCC.unpaired.right.fastq"   "Ar-7_ATGTCA.unpaired.right.fastq"
```

```r

# Loop across each sample and quantify expression

# NOTE - samples listed in same order as given by the above lists
samples <- c("A22-0", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", 
    "A22-31", "A22-35", "A22-38", "A22-3", "A22-7", "Ar-0", "Ar-10", "Ar-14", "Ar-17", 
    "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38", "Ar-3", "Ar-7")

for (j in 1:length(samples)) {
    message("Start expression quantification for sample ", samples[j], ": ", Sys.time())
    quantdir <- paste(sfexpressionroot, samples[j], "_quant", sep = "")
    samp.pos <- grep(paste(paste(samples[j], "_", sep = "")), paired.left)
    samp.paired.l <- paste(readdir, paired.left[samp.pos], sep = "")
    samp.paired.r <- paste(readdir, paired.right[samp.pos], sep = "")
    samp.unpaired.l <- paste(readdir, unpaired.left[samp.pos], sep = "")
    samp.unpaired.r <- paste(readdir, unpaired.right[samp.pos], sep = "")
    sailfishcmd <- paste("sailfish quant -i ", sfindex, " -o ", quantdir, " -l 'T=SE:S=U' -r ", 
        samp.paired.l, " ", samp.paired.r, " ", samp.unpaired.l, " ", samp.unpaired.r, 
        " -p 4", sep = "")
    # print(sailfishcmd)
    system(sailfishcmd)
    message("Done with expression quantification for sample ", samples[j], ": ", 
        Sys.time(), "\n")
}
```

```
## Start expression quantification for sample A22-0: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-0: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-10: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-10: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-14: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-14: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-17: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-17: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-21: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-21: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-24: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-24: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-28: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-28: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-31: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-31: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-35: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-35: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-38: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-38: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-3: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-3: 2014-05-08 10:46:02
## 
## Start expression quantification for sample A22-7: 2014-05-08 10:46:02
## Done with expression quantification for sample A22-7: 2014-05-08 10:46:02
## 
## Start expression quantification for sample Ar-0: 2014-05-08 10:46:02
## Done with expression quantification for sample Ar-0: 2014-05-08 10:46:02
## 
## Start expression quantification for sample Ar-10: 2014-05-08 10:46:02
## Done with expression quantification for sample Ar-10: 2014-05-08 10:46:02
## 
## Start expression quantification for sample Ar-14: 2014-05-08 10:46:02
## Done with expression quantification for sample Ar-14: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-17: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-17: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-21: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-21: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-24: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-24: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-28: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-28: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-31: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-31: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-35: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-35: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-38: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-38: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-3: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-3: 2014-05-08 10:46:03
## 
## Start expression quantification for sample Ar-7: 2014-05-08 10:46:03
## Done with expression quantification for sample Ar-7: 2014-05-08 10:46:03
```


This generated a directory for each sample

A22-0_quant, A22-10_quant, A22-14_quant, A22-17_quant, A22-21_quant, A22-24_quant, A22-28_quant, A22-31_quant, A22-35_quant, A22-38_quant, A22-3_quant, A22-7_quant, Ar-0_quant, Ar-10_quant, Ar-14_quant, Ar-17_quant, Ar-21_quant, Ar-24_quant, Ar-28_quant, Ar-31_quant, Ar-35_quant, Ar-38_quant, Ar-3_quant, Ar-7_quant

and within each directory there are the following r:

logs, quant_bias_corrected.sf, quant.sf, reads.count_info, reads.sfc

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
|   7    |  0.98  |
|  10.5  |   1    |
|   14   |  0.99  |
|  17.5  |  0.98  |
|   21   |  0.98  |
|  24.5  |  0.99  |
|   28   |  0.99  |
|  31.5  |  0.99  |
|   35   |  0.99  |
|  38.5  |  0.99  |

Table: correlations between colonies at each temperature treatment



**Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the combined expression data set for the two species.** 


```r
A22.TPM[, `:=`(colony, "A22")]
```

```
##                      Transcript Length    TPM   RPKM   KPKM EstimatedNumReads
##       1:         0|*|Contig6267   9990 0.0786 0.0917 0.0917              2963
##       2:         0|*|Contig6267   9990 0.0356 0.0531 0.0531              1371
##       3:         0|*|Contig6267   9990 0.0706 0.0856 0.0856              2520
##       4:         0|*|Contig6267   9990 0.0395 0.0492 0.0492              1712
##       5:         0|*|Contig6267   9990 0.0238 0.0337 0.0337              1079
##      ---                                                                     
## 1266428: 9|*|comp147140_c0_seq1   9030 0.9003 1.3145 1.3145             35271
## 1266429: 9|*|comp147140_c0_seq1   9030 0.8852 1.2721 1.2721             25430
## 1266430: 9|*|comp147140_c0_seq1   9030 1.0720 1.4999 1.4999             32576
## 1266431: 9|*|comp147140_c0_seq1   9030 0.7239 1.0022 1.0022             23921
## 1266432: 9|*|comp147140_c0_seq1   9030 0.5655 0.7745 0.7745             23814
##             V7 sample  val colony
##       1:  37.0  A22-0  0.0    A22
##       2:  17.1  A22-3  3.5    A22
##       3:  31.5  A22-7  7.0    A22
##       4:  21.4 A22-10 10.5    A22
##       5:  13.4 A22-14 14.0    A22
##      ---                         
## 1266428: 439.2 A22-24 24.5    A22
## 1266429: 317.0 A22-28 28.0    A22
## 1266430: 406.1 A22-31 31.5    A22
## 1266431: 298.8 A22-35 35.0    A22
## 1266432: 297.4 A22-38 38.5    A22
```

```r
Ar.TPM[, `:=`(colony, "Ar")]
```

```
##                      Transcript Length    TPM  RPKM  KPKM EstimatedNumReads
##       1:         0|*|Contig6267   9990 0.0614 0.104 0.104              1438
##       2:         0|*|Contig6267   9990 0.0920 0.140 0.140              2416
##       3:         0|*|Contig6267   9990 0.1615 0.171 0.171              4409
##       4:         0|*|Contig6267   9990 0.2095 0.307 0.307              5095
##       5:         0|*|Contig6267   9990 0.1562 0.235 0.235              5037
##      ---                                                                   
## 1266428: 9|*|comp147140_c0_seq1   9030 0.6406 1.019 1.019             16612
## 1266429: 9|*|comp147140_c0_seq1   9030 0.5848 0.968 0.968             14166
## 1266430: 9|*|comp147140_c0_seq1   9030 0.7167 1.212 1.212             17725
## 1266431: 9|*|comp147140_c0_seq1   9030 0.4328 0.793 0.793              8730
## 1266432: 9|*|comp147140_c0_seq1   9030 0.4626 0.847 0.847             12716
##             V7 sample  val colony
##       1:  17.8   Ar-0  0.0     Ar
##       2:  29.9   Ar-3  3.5     Ar
##       3:  55.1   Ar-7  7.0     Ar
##       4:  63.2  Ar-10 10.5     Ar
##       5:  62.4  Ar-14 14.0     Ar
##      ---                         
## 1266428: 205.6  Ar-24 24.5     Ar
## 1266429: 175.3  Ar-28 28.0     Ar
## 1266430: 219.2  Ar-31 31.5     Ar
## 1266431: 108.1  Ar-35 35.0     Ar
## 1266432: 157.1  Ar-38 38.5     Ar
```

```r
TPM.dt <- rbind(A22.TPM, Ar.TPM)
TPM.dt$colony <- as.factor(TPM.dt$colony)
str(TPM.dt)
```

```
## Classes 'data.table' and 'data.frame':	2532864 obs. of  10 variables:
##  $ Transcript       : chr  "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" ...
##  $ Length           : int  9990 9990 9990 9990 9990 9990 9990 9990 9990 9990 ...
##  $ TPM              : num  0.0786 0.0356 0.0706 0.0395 0.0238 ...
##  $ RPKM             : num  0.0917 0.0531 0.0856 0.0492 0.0337 ...
##  $ KPKM             : num  0.0917 0.0531 0.0856 0.0492 0.0337 ...
##  $ EstimatedNumReads: num  2963 1371 2520 1712 1079 ...
##  $ V7               : num  37 17.1 31.5 21.4 13.4 ...
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
length(unique(TPM.dt.sub$Transcript))
```

```
## [1] 105536
```


### Remove *contaminant* transcripts


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
## Classes 'data.table' and 'data.frame':	2196942 obs. of  10 variables:
##  $ Transcript       : chr  "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" ...
##  $ Length           : int  9990 9990 9990 9990 9990 9990 9990 9990 9990 9990 ...
##  $ TPM              : num  0.0786 0.0614 0.0356 0.092 0.0395 ...
##  $ RPKM             : num  0.0917 0.1036 0.0531 0.1404 0.0492 ...
##  $ KPKM             : num  0.0917 0.1036 0.0531 0.1404 0.0492 ...
##  $ EstimatedNumReads: num  2963 1438 1371 2416 1712 ...
##  $ V7               : num  37 17.8 17.1 29.9 21.4 ...
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
## [1] 99861
```



### Regression-model to identify thermally-responsive genes

To identify transcripts (roughly equivalent to genes) that show thermal responsiveness, I fit the following linear model to each transcript:

$$ log(TPM + 1) = \beta_0 + \beta_1(colony) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(colony * temp) + \beta_5(colony * temp^2) + \epsilon $$

where TPM is transcripts per million. 


(1) Identify transcripts with overall significant model fit. Adjust *P* values for multiple testing using FDR and retain transcripts with FDR < 0.05. Use log-transformed response to account for outliers.



```r
# define model for RxN function
model <- "log(TPM+1) ~ colony + val + I(val^2) + colony:val + colony:I(val^2)"

# define function for `ddply` to calculate overall model P value
modpFunc <- function(foo) {
    # fit lm
    lmout <- eval(parse(text = paste("lm(", model, ", data = foo)", sep = "")))
    
    # calculate overall model significance
    f <- summary(lmout)$fstatistic
    pval <- unname(pf(f[1], f[2], f[3], lower.tail = F))
    attributes(pval) <- NULL
    
    # return pvalue
    c(pval = pval)
}

RxNpval <- ddply(TPM.dt.sub, .(Transcript), .inform = "TRUE", modpFunc)
```


Of the 99861 transcripts, 22339 have models with P < 0.05.

Many of these are likely false positives, so I adjust P-values using false discovery rate (FDR). Only those transcripts with less than 5% FDR are retained as significant. 


```r
RxNpval$padj <- p.adjust(RxNpval$pval, method = "fdr")
# Plot FDR values against initial pvalues
par(mfrow = c(2, 1))
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


At the 5% FDR significance threshold, there are 10597 transcripts with an overall significant model.


(2) Fit linear model to overall significant transcripts; perform stepAIC to retain only significant terms, and save `lm` output to list


```r
# first, define function to use
lmFunc <- function(bar) {
    # fit lm
    lmout <- eval(parse(text = paste("lm(", model, ", data = bar)", sep = "")))
    # stepAIC to drop non-significant terms
    f.lmout <- stepAIC(lmout)
    # return final model
    return(f.lmout)
}

# apply function to all responsive transcripts
RxNlmAIC <- try(dlply(sig.TPM.dt.sub, .(Transcript), lmFunc))
```



### Grouping of thermally-responsive transcripts

The set of transcripts with significant expression patterns include those with expression that differs by colony, temperature and the interaction of colony and temperature. In this section, I am specifically interested in the thermally-responsive transcripts (temperature and colony x temperature) so I subset the significant transcripts to examine these. 


```r
grepFunc <- function(lmitem, term = NA) {
    coefs <- coefficients(lmitem)
    coef.names <- names(coefs)
    keep <- if (length(grep(term, coef.names)) > 0) 
        TRUE else FALSE
}

interaction.lms <- RxNlmAIC[which(Map(grepFunc, RxNlmAIC, term = "colonyAr:") == 
    TRUE)]
other.lms <- RxNlmAIC[setdiff(names(RxNlmAIC), names(interaction.lms))]
temperature.lms <- other.lms[which(Map(grepFunc, other.lms, term = "val") == TRUE)]
colony.lms <- other.lms[setdiff(names(other.lms), names(temperature.lms))]
responsive.lms <- c(temperature.lms, interaction.lms)
rm(other.lms)
```



|    Coefficient     |  Number.significant  |
|:------------------:|:--------------------:|
|       Total        |        10597         |
|       Colony       |         1441         |
|    Temperature     |         2248         |
| Temperature:Colony |         6908         |

Table: Number of transcripts of 99,861 total with expression that depends on colony, temperature or their interaction at 5% FDR.



### Thermal-response functional types ###

The previous section identified the transcripts with thermally-responsive expression. In this section, I determine the shape of the expression response to temperature for each transcript. Catego
ries of expression response are:

* High - increase expression with temperature
* Low - decrease expression with temperature
* Intermediate - maximum expression at intermediate temperatures (14 - 28C)
* Bimodal - expressed greater than two standard deviations of expression at both low and high temperatures

I do this first for the thermally-responsive transcripts where there is no interaction with colony. For the transcripts where thermal-responsive expression depends on colony, I determine the functional type of the expression response separately for each colony. 



```r
# calculate response type for responsive transcripts
interaction.response.type <- ldply(interaction.lms, .progress = "none", .inform = TRUE, 
    RxNtype)
stopifnot(nrow(interaction.response.type) == length(interaction.lms))

# calculate response types for transcripts without interactions
temperature.response.type <- ldply(temperature.lms, .progress = "none", .inform = TRUE, 
    RxNtype)
stopifnot(nrow(temperature.response.type) == length(temperature.lms))

# merge results
Ap.response.type <- rbind(interaction.response.type, temperature.response.type)
colnames(Ap.response.type)[which(colnames(Ap.response.type) == ".id")] <- "Transcript"
str(Ap.response.type)
```

```
## 'data.frame':	9156 obs. of  9 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100015|*|comp3543055_c0_seq1" "100067|*|comp3557646_c0_seq1" "100089|*|comp11313_c1_seq1" ...
##  $ A22.max   : num  38.5 0 NA 18.5 NA 0 NA 0 0 38.5 ...
##  $ A22.min   : num  0 20.5 NA 38.5 NA 23 NA 28.5 26 0 ...
##  $ A22.opt   : num  1.025 0.945 1 1.057 1 ...
##  $ A22.type  : chr  "High" "Bimodal" "NotResp" "Intermediate" ...
##  $ Ar.max    : num  0 38.5 0 38.5 38.5 18 0 0 NA NA ...
##  $ Ar.min    : num  38.5 0 20.5 13.5 18.5 38.5 38.5 38.5 NA NA ...
##  $ Ar.opt    : num  1.198 1.004 0.935 0.962 0.897 ...
##  $ Ar.type   : chr  "Low" "High" "Bimodal" "High" ...
```

```r

# save results to file
write.table(file = paste(resultsdir, "Ap_responsive_transcripts.txt", sep = ""), 
    Ap.response.type, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
```


Next, I compare the number of thermally-responsive in each response category between the two colonies. 


```r
A22.type.table <- table(Ap.response.type[, "A22.type"])
Ar.type.table <- table(Ap.response.type[, "Ar.type"])

# table
Ap.type.table <- rbind(A22.type.table, Ar.type.table)

# Pearson Chi-square test
chisq.test(Ap.type.table)
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  Ap.type.table
## X-squared = 1300, df = 4, p-value < 2.2e-16
```

```r

# Reorganize for plotting to show overlap among categories between colonies
type.table <- table(Ap.response.type[, "A22.type"], Ap.response.type[, "Ar.type"])

# reorder
tt2 <- type.table[c(2, 4, 1, 3, 5), c(2, 4, 1, 3, 5)]

# plot
mosaicplot(tt2, xlab = "AcNC", ylab = "ApVT", main = "Mosaic plot of responsive transcripts")
```

![plot of chunk temperature_response](figure/temperature_response.png) 


The number of thermally-responsive in each response category differs between the colonies, with the msot transcripts expressed at *Low* temperatures in both colonies. For *ApVT*, an equal number of transcripts are expressed at *High* and *Bimodal*, followed by *Intermediate* transcripts. For *AcNC*, transcripts expressed at *Intermediate* temperatures are next most common, followed by *High* and *Bimodal*. 

Interestingly, nearly half of the *High* genes in *AcNC* are *Low* in *ApVT*, and vice versa. In contrast, most of the *Low* genes in one species are also *Low* in the other species.



|       &nbsp;       |  High  |  Low  |  Bimodal  |  Intermediate  |  NotResp  |
|:------------------:|:------:|:-----:|:---------:|:--------------:|:---------:|
|      **High**      |  306   |  420  |    165    |      153       |    221    |
|      **Low**       |  445   | 2632  |    175    |      1372      |    293    |
|    **Bimodal**     |  150   |  302  |    291    |      567       |    222    |
|  **Intermediate**  |   66   |  164  |    111    |      569       |    18     |
|    **NotResp**     |   92   |  245  |    108    |       69       |     0     |

Table: Number of transcripts with maximum expression at high, low, intermediate, both high and low (bimodal) temperatures or are not thermally-responsivefor each colony and their overlap.



Table 4 shows the number of transcripts that fall into each expression type for each each colony. The totals for each colony include the 2248 transcripts that have consistent temperature responses between the two colonies. 



## Colony-level comparison ##

In this section, I compare the thermal reactionome between the *Ar* and *A22* colonies. 

### Plasticity versus constitutive expression

Selection may have acted in response to thermal stress such that some genes are constitutively activated in one colony while plastically-expressed in the other. I evaluate this hypothesis by comparing expression levels at the at optimum (19.5C) temperature between the two colonies for genes in that are either in the 'High' or 'Low' expression group in the other colony. Specifically, I expect that genes upregulated at high temperatures in *A22* are more highly expressed at 19.5C in *Ar*? 


```r
# list of transcripts that are 'high' expressed in A22
A22.high.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "High"), 
    ]

# Compare expression at optimum temp (A22.opt) between colonies using t-test
t.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  A22.high.transcripts$A22.opt and A22.high.transcripts$Ar.opt
## t = -0.184, df = 2342, p-value = 0.8537
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -1177   975
## sample estimates:
## mean of x mean of y 
##       373       475
```

```r
boxplot(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)
```

![plot of chunk optimum_expression_comparison](figure/optimum_expression_comparison1.png) 

```r

# T test on log-transformed values to control for outliers
t.test(log(A22.high.transcripts$A22.opt + 1), log(A22.high.transcripts$Ar.opt + 1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(A22.high.transcripts$A22.opt + 1) and log(A22.high.transcripts$Ar.opt + 1)
## t = -3.62, df = 2514, p-value = 0.0002992
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2630 -0.0782
## sample estimates:
## mean of x mean of y 
##      1.11      1.28
```

```r
boxplot(log(A22.high.transcripts$A22.opt + 1), log(A22.high.transcripts$Ar.opt + 
    1))
```

![plot of chunk optimum_expression_comparison](figure/optimum_expression_comparison2.png) 


The `t.test` fails to account for the many orders of magnitude difference in expression among transcripts, e.g. non-equal variances. This problem is the key issue in the analysis of differential expression (<a href="http://dx.doi.org/10.1186/1471-2105-11-94">Bullard et al. 2010</a>; <a href="http://dx.doi.org/10.1038/nprot.2013.099">Anders et al. 2013</a>). As my goal is simply to determine if expression is typically greater at 19.25C in *Ar* than *A22* for genes that are up-regulated at high temperatures in *A22*, I use a non-parametric Wilcoxon signed rank-test


```r
(w1 <- wilcox.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE))
```

```
## Error: object 'A22.high.transcripts' not found
```


Consistent with expectation, there is greater expression at 19.5C for *Ar* than *A22* transcripts for the set of transcripts that are transcripts that are up-regulated at high temperatures in *A22*. Note that A22 had the larger library size so if this was due to TPM not correctly accounting for differences in reads between samples, we would expect to see a positive instead of negative value here.

Next I test the converse, that transcripts that are up-regulated at low temperatures in *Ar* are more highly-expressed at 19.25C temperatures in *A22*. 


```r
# list of transcripts that are 'high' expressed in Ar
Ar.low.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Low"), 
    ]

# compare expression at optimum temp (Ar.opt) between colonies using t-test
t.test(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.low.transcripts$A22.opt and Ar.low.transcripts$Ar.opt
## t = 0.823, df = 6866, p-value = 0.4107
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -5.42 13.27
## sample estimates:
## mean of x mean of y 
##        16        12
```

```r
boxplot(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt)
```

![plot of chunk Ar_low_wilcoxon](figure/Ar_low_wilcoxon1.png) 

```r
# t-test with log values
t.test(log(Ar.low.transcripts$A22.opt + 1), log(Ar.low.transcripts$Ar.opt + 1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(Ar.low.transcripts$A22.opt + 1) and log(Ar.low.transcripts$Ar.opt + 1)
## t = -6.63, df = 7511, p-value = 3.547e-11
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1737 -0.0944
## sample estimates:
## mean of x mean of y 
##      1.26      1.39
```

```r
boxplot(log(Ar.low.transcripts$A22.opt + 1), log(Ar.low.transcripts$Ar.opt + 1))
```

![plot of chunk Ar_low_wilcoxon](figure/Ar_low_wilcoxon2.png) 

```r

# Wilcoxon signed rank-test
(w2 <- wilcox.test(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE))
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  Ar.low.transcripts$A22.opt and Ar.low.transcripts$Ar.opt
## V = 2210400, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.450 -0.357
## sample estimates:
## (pseudo)median 
##         -0.402
```


Counter to expectations, expression at 19.5C is also greater in *Ar* than *A22* for transcripts upregulated at low temperatures in *Ar*. 

To confirm that there are not sample-level issues, I performed the same comparison using transcripts where I do *not* expect to see a difference in expression.


```r
# list of transcripts that are 'Intermediate' expressed in Ar
Ar.int.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate"), 
    ]
# T test
t.test(Ar.int.transcripts$A22.opt, Ar.int.transcripts$Ar.opt)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.int.transcripts$A22.opt and Ar.int.transcripts$Ar.opt
## t = 0.875, df = 5129, p-value = 0.3815
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -9.76 25.51
## sample estimates:
## mean of x mean of y 
##      40.6      32.7
```

```r
t.test(log(Ar.int.transcripts$A22.opt + 1), log(Ar.int.transcripts$Ar.opt + 1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(Ar.int.transcripts$A22.opt + 1) and log(Ar.int.transcripts$Ar.opt + 1)
## t = -7.01, df = 5454, p-value = 2.717e-12
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.277 -0.156
## sample estimates:
## mean of x mean of y 
##      1.38      1.59
```

```r
# Wilcoxon signed rank-test
wilcox.test(Ar.int.transcripts$A22.opt, Ar.int.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE)
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  Ar.int.transcripts$A22.opt and Ar.int.transcripts$Ar.opt
## V = 1133332, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.607 -0.467
## sample estimates:
## (pseudo)median 
##         -0.535
```

```r

# list of transcripts that are 'Intermediate' expressed in A22
A22.int.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate"), 
    ]
# T test
t.test(A22.int.transcripts$A22.opt, A22.int.transcripts$Ar.opt)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  A22.int.transcripts$A22.opt and A22.int.transcripts$Ar.opt
## t = 0.227, df = 1849, p-value = 0.8203
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -33.6  42.4
## sample estimates:
## mean of x mean of y 
##      50.0      45.6
```

```r
t.test(log(A22.int.transcripts$A22.opt + 1), log(A22.int.transcripts$Ar.opt + 1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(A22.int.transcripts$A22.opt + 1) and log(A22.int.transcripts$Ar.opt + 1)
## t = -1.62, df = 1854, p-value = 0.106
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2086  0.0201
## sample estimates:
## mean of x mean of y 
##      1.44      1.53
```

```r
# Wilcoxon signed rank-test
wilcox.test(A22.int.transcripts$A22.opt, A22.int.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE)
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  A22.int.transcripts$A22.opt and A22.int.transcripts$Ar.opt
## V = 155963, p-value = 8.222e-13
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.412 -0.220
## sample estimates:
## (pseudo)median 
##         -0.309
```


The non-parametric test for both comparisions also finds greater constitutive expression in *Ar* than *A22*.


### Thermal tolerance indicated by region of constant gene expression ###

The 'Intermediate' expressed transcripts are core molecular processes that are expressed at non-stressful temperatures, and shut-off when the organism experiences thermall stress. We hypothesized that if the more southern *Ar* colony was more thermally-tolerant than *A22*, transcripts with 'Intermediate' expression (10-30C) would be active across a wider range of temperatures. To test this with our data, we calculated the standard deviation of the expression function for each temperature transcript that was 'Intermediate' expressed in each colony.


```r
# extract 'Intermediate' expressed transcripts for A22 colony
A22.int.lm <- RxNlmAIC[A22.int.transcripts$Transcript]
length(A22.int.lm)
```

```
## [1] 928
```

```r

# apply `transcriptSD` function to all transcripts
A22.int.sd <- unlist(Map(transcriptSD, A22.int.lm, colony = "A22"))
A22.int.sd <- data.frame(colony = "A22", exp_sd = A22.int.sd)

# repeat for Ar
Ar.int.lm <- RxNlmAIC[Ar.int.transcripts$Transcript]
Ar.int.sd <- unlist(Map(transcriptSD, Ar.int.lm, colony = "Ar"))
Ar.int.sd <- data.frame(colony = "Ar", exp_sd = Ar.int.sd)

# T-test comparing standard deviation of expression between colonies
(t1 <- t.test(Ar.int.sd$exp_sd, A22.int.sd$exp_sd, alternative = "two.sided"))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.int.sd$exp_sd and A22.int.sd$exp_sd
## t = 9.28, df = 1323, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.339 0.520
## sample estimates:
## mean of x mean of y 
##     10.36      9.93
```


Consistent with our hypothesis, 'Intermediate' transcripts in *Ar* are expressed over a significantly wider range of temperatures than in *A22*. 

![plot of chunk plotthermalbreadth](figure/plotthermalbreadth.png) 



### Thermal sensitivity indicated by response of bimodally-expressed transcripts ###

As the converse of the above hypothesis, a colony that is especially thermally-sensitive is likely to activate expression of molecular processes more quickly. We tested this using the same approach as for the 'Intermediate' transcripts, but using the inverse of the 'Bimodal' expressed transcripts. 



```r
# extract 'Bimodal' expressed transcripts for A22 colony
A22.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$A22.type == "Bimodal"), 
    "Transcript"]]
length(A22.bim.lm)
```

```
## [1] 1532
```

```r

# apply `transcriptSD` function to all transcripts
A22.bim.sd <- unlist(Map(transcriptSD, A22.bim.lm, colony = "A22"))

# repeat for Ar
Ar.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$Ar.type == "Bimodal"), 
    "Transcript"]]
length(Ar.bim.lm)
```

```
## [1] 850
```

```r
Ar.bim.sd <- unlist(Map(transcriptSD, Ar.bim.lm, colony = "Ar"))

# t-test to compare standard deviation of 'Bimodal' transcripts between colonies
(t2 <- t.test(Ar.bim.sd, A22.bim.sd))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.bim.sd and A22.bim.sd
## t = -0.43, df = 1465, p-value = 0.6671
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.182  0.116
## sample estimates:
## mean of x mean of y 
##      12.6      12.7
```


No difference in the standard deviation of expression for bimodally-expressed transcripts between colonies.


### Compare peak expression among colonies ###

Probability density function of peak expression for transcripts that differ in expression between *A22* and *Ar*


```r
# reshape data
Ap.df <- data.frame(Transcript = rep(Ap.response.type$Transcript, times = 2), colony = rep(c("ApVT", 
    "ApNC"), each = length(Ap.response.type$Transcript)), max.val = c(Ap.response.type$A22.max, 
    Ap.response.type$Ar.max))

png(paste(resultsdir, "PDF_expression_all.png", sep = ""))
g3 <- ggplot(Ap.df, aes(x = max.val, fill = colony)) + geom_density(alpha = 0.2, 
    position = "identity") + # scale_fill_manual(name = 'Colony', values = c('white', 'black')) +
scale_y_continuous(name = "Density") + scale_x_continuous(name = "Temperature of maximum expression")
suppressWarnings(print(g3))
dev.off()
```

```
## pdf 
##   2
```


For the transcripts that differed in thermal responsiveness due to temperature, was the difference primarily due to differences in the mean value of expression, slope, curvature of a higher order effect? To test this, I use the framework of <a href="http://dx.doi.org/10.1086/675302">Murren et al. (2014)</a>:

- *Offset, O*: overall difference in the mean expression value across all temperatures
- *Slope, S*: difference in overall slope
- *Curvature, C*: average difference in curvature of the reaction norm
- *Wiggle, W*: variability in shape not captured by the previous three measures


```r
varshape.out <- ldply(interaction.lms, .progress = "none", RxNvarshape)

# plot
boxplot(varshape.out$A22.mean, varshape.out$Ar.mean)
```

![plot of chunk offset](figure/offset1.png) 

```r
boxplot(log(varshape.out$A22.mean), log(varshape.out$Ar.mean))
```

![plot of chunk offset](figure/offset2.png) 

```r

# calculate differences in mean, slope, curvature and wiggle between colonies for
# each transcript offset = mean1 - mean2
varshape.out$offset <- varshape.out$A22.mean - varshape.out$Ar.mean
# to compare among transcripts that differ in scale of expression, standardize by
# overall mean for each transcript
varshape.out$s.offset <- NA
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.offset"] <- (varshape.out[i, "offset"]/mean(c(varshape.out[i, 
        "A22.mean"], varshape.out[i, "Ar.mean"])))
}

# t-test
(t3 <- t.test(varshape.out$s.offset, alternative = "two.sided"))
```

```
## 
## 	One Sample t-test
## 
## data:  varshape.out$s.offset
## t = -11.7, df = 6907, p-value < 2.2e-16
## alternative hypothesis: true mean is not equal to 0
## 95 percent confidence interval:
##  -0.0899 -0.0642
## sample estimates:
## mean of x 
##    -0.077
```

```r

# slope = slope1 - slope2
varshape.out$slope <- varshape.out$A22.slope - varshape.out$Ar.slope
# to compare among transcripts that differ in scale of expression, standardize by
# overall mean for each transcript
varshape.out$s.slope <- NA
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.slope"] <- (varshape.out[i, "slope"]/mean(c(varshape.out[i, 
        "A22.slope"], varshape.out[i, "Ar.slope"])))
}
# for transcripts where both samples had slope=0, change s.slope from Inf to 0 as
# there is truly no difference in slope
varshape.out[which(varshape.out$Ar.slope == 0 & varshape.out$A22.slope == 0), "s.slope"] <- 0

# t-test
(t4 <- t.test(varshape.out$s.slope, alternative = "two.sided"))
```

```
## 
## 	One Sample t-test
## 
## data:  varshape.out$s.slope
## t = 12.8, df = 6907, p-value < 2.2e-16
## alternative hypothesis: true mean is not equal to 0
## 95 percent confidence interval:
##  0.192 0.262
## sample estimates:
## mean of x 
##     0.227
```

```r


# curvature = curvature1 - curvature2
varshape.out$curvature <- varshape.out$A22.curve - varshape.out$Ar.curve
# to compare among transcripts that differ in scale of expression, standardize by
# overall mean for each transcript
varshape.out$s.curvature <- NA
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.curvature"] <- (varshape.out[i, "curvature"]/mean(c(varshape.out[i, 
        "A22.curve"], varshape.out[i, "Ar.curve"])))
}
# change 'Inf' values to 0
varshape.out[which(varshape.out$Ar.curve == 0 & varshape.out$A22.curve == 0), "s.curvature"] <- 0
# t-test
(t.test(varshape.out$s.curvature, alternative = "two.sided"))
```

```
## 
## 	One Sample t-test
## 
## data:  varshape.out$s.curvature
## t = 2.45, df = 6907, p-value = 0.01423
## alternative hypothesis: true mean is not equal to 0
## 95 percent confidence interval:
##  0.00898 0.08056
## sample estimates:
## mean of x 
##    0.0448
```

```r


# wiggle
varshape.out$wiggle <- varshape.out$A22.wiggle - varshape.out$Ar.wiggle
varshape.out$s.wiggle <- NA
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.wiggle"] <- (varshape.out[i, "wiggle"]/mean(c(varshape.out[i, 
        "A22.wiggle"], varshape.out[i, "Ar.wiggle"])))
}
# change 'Inf' values to 0
varshape.out[which(varshape.out$Ar.wiggle == 0 & varshape.out$A22.wiggle == 0), "s.wiggle"] <- 0
# note that only 407 transcripts have any 'wiggle'
length(which(varshape.out$s.wiggle != 0))
```

```
## [1] 531
```

```r

# t-test
(t5 <- t.test(varshape.out$s.wiggle, alternative = "two.sided"))
```

```
## 
## 	One Sample t-test
## 
## data:  varshape.out$s.wiggle
## t = -17.4, df = 6907, p-value < 2.2e-16
## alternative hypothesis: true mean is not equal to 0
## 95 percent confidence interval:
##  -0.126 -0.100
## sample estimates:
## mean of x 
##    -0.113
```


Next, I partition the differences in the reaction norms into the variation explained by differences in the trait means (offset), slope, curvature and wiggle.



```r
# weighted mean-standardized values of each measure
mean(abs((varshape.out$s.offset)))
```

```
## [1] 0.394
```

```r
mean(abs((varshape.out$s.slope)))
```

```
## [1] 1.36
```

```r
mean(abs((varshape.out$s.curvature)))
```

```
## [1] 1.32
```

```r
mean(abs((varshape.out$s.wiggle)))
```

```
## [1] 0.152
```

```r

# take absolute value of each value and sum to get total differences between
# reaction norms
varshape.out$s.total <- abs(varshape.out$s.offset) + abs(varshape.out$s.slope) + 
    abs(varshape.out$s.curvature) + abs(varshape.out$s.wiggle)

(t6 <- t.test(varshape.out$s.total, alternative = "two.sided"))
```

```
## 
## 	One Sample t-test
## 
## data:  varshape.out$s.total
## t = 193, df = 6907, p-value < 2.2e-16
## alternative hypothesis: true mean is not equal to 0
## 95 percent confidence interval:
##  3.20 3.26
## sample estimates:
## mean of x 
##      3.23
```

```r

# variation in reaction norms due to differences in mean
varshape.out$prop.mean <- abs(varshape.out$s.offset)/varshape.out$s.total
varshape.out$prop.slope <- abs(varshape.out$s.slope)/varshape.out$s.total
varshape.out$prop.curve <- abs(varshape.out$s.curvature)/varshape.out$s.total
varshape.out$prop.wiggle <- abs(varshape.out$s.wiggle)/varshape.out$s.total

# Mean proportion and 95% CI of total variation of each measure
mean(varshape.out$prop.mean)
```

```
## [1] 0.132
```

```r
quantile(varshape.out$prop.mean, probs = c(0.05, 0.5, 0.95))
```

```
##     5%    50%    95% 
## 0.0152 0.1021 0.3440
```

```r

mean(varshape.out$prop.slope)
```

```
## [1] 0.449
```

```r
quantile(varshape.out$prop.slope, probs = c(0.05, 0.5, 0.95))
```

```
##    5%   50%   95% 
## 0.101 0.443 0.896
```

```r

mean(varshape.out$prop.curve)
```

```
## [1] 0.39
```

```r
quantile(varshape.out$prop.curve, probs = c(0.05, 0.5, 0.95))
```

```
##    5%   50%   95% 
## 0.000 0.447 0.699
```

```r

mean(varshape.out$prop.wiggle)
```

```
## [1] 0.0295
```

```r
quantile(varshape.out$prop.wiggle, probs = c(0.05, 0.5, 0.95))
```

```
##    5%   50%   95% 
## 0.000 0.000 0.312
```


## Temperature of gene activation

Thermally-responsive genes could also differ in the temperatures at which they have increased or decreased expression in response to temperature changes. To examine this, I 








### Identification of biological-functions associated with temperature responses

## Functional annotation

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  


```r
setkey(annotationtable, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
signif.transcripts <- annotationtable[signif.transcripts]
setnames(signif.transcripts, "Sequence.Name", "Transcript")
```






## Gene set enrichment analysis ##

I perform gene set enrichment analysis below, but a quick `grep` shows that there are 26 transcripts with GO term "response to stress", though this is not enriched compared to the frequency of this term in the full dataset.
  

```r
# GO 'response to stress' hits in responsive transcripts
GO0006950.responsive <- responsive.lms.ann.type[grep("GO:0006950", responsive.lms.ann.type$GO.Biological.Process), 
    list(Transcript, best.hit.to.nr, A22.type, Ar.type)]

# in high category
GO0006950.responsive[union(with(GO0006950.responsive, grep("High", Ar.type)), with(GO0006950.responsive, 
    grep("High", A22.type))), ]
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
##  8:  21598|*|comp142101_c0_seq1
##  9:  23441|*|comp114823_c1_seq1
## 10:  32312|*|comp933733_c0_seq1
## 11: 37154|*|comp1975086_c0_seq1
## 12: 47691|*|comp1460938_c0_seq1
## 13: 51985|*|comp1012776_c0_seq1
##                                                               best.hit.to.nr
##  1:               gi|332023134|gb|EGI63390.1| Sugar transporter ERD6-like 6 
##  2:                       gi|194716766|gb|ACF93232.1| heat shock protein 90 
##  3:                          gi|337757286|emb|CBZ98843.1| 60 kDa chaperonin 
##  4:         gi|332022897|gb|EGI63169.1| Protein lethal(2)essential for life 
##  5:         gi|332021988|gb|EGI62314.1| Heat shock 70 kDa protein cognate 4 
##  6:                               gi|50418863|ref|XP_457952.1| DEHA2C06072p 
##  7:             gi|322789999|gb|EFZ15075.1| hypothetical protein SINV_03446 
##  8:         gi|332018201|gb|EGI58806.1| Protein lethal(2)essential for life 
##  9:      gi|443696809|gb|ELT97425.1| hypothetical protein CAPTEDRAFT_194915 
## 10:     gi|367054010|ref|XP_003657383.1| hypothetical protein THITE_2156506 
## 11:             gi|46115086|ref|XP_383561.1| hypothetical protein FG03385.1 
## 12: gi|302922354|ref|XP_003053448.1| hypothetical protein NECHADRAFT_102357 
## 13:                       gi|227018528|gb|ACP18866.1| heat shock protein 30 
##     A22.type Ar.type
##  1:     High    High
##  2:  Bimodal    High
##  3:  NotResp    High
##  4:      Low    High
##  5:     High    High
##  6:  NotResp    High
##  7:      Low    High
##  8:     High Bimodal
##  9:     High     Low
## 10:     High Bimodal
## 11:     High NotResp
## 12:     High NotResp
## 13:     High NotResp
```

```r
# in low category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Low", Ar.type)), with(GO0006950.responsive, 
    grep("Low", A22.type))), ]
```

```
##                      Transcript
##  1:   1038|*|comp150483_c5_seq3
##  2:  11281|*|comp146961_c0_seq1
##  3:          19475|*|Contig1438
##  4:  20215|*|comp145360_c0_seq1
##  5:  23441|*|comp114823_c1_seq1
##  6:   2604|*|comp148324_c0_seq4
##  7:   4273|*|comp150636_c5_seq1
##  8: 50934|*|comp3428507_c0_seq1
##  9:    6075|*|comp92770_c0_seq1
## 10:   6438|*|comp150878_c2_seq2
## 11:   6689|*|comp141130_c0_seq2
## 12:  80544|*|comp132706_c0_seq1
## 13:           9372|*|Contig4757
## 14:  12704|*|comp144775_c1_seq1
## 15:     14|*|comp150262_c0_seq1
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
##  3:                                  gi|322799248|gb|EFZ20646.1| hypothetical protein SINV_03807 
##  4:                                  gi|332020393|gb|EGI60813.1| G-protein coupled receptor Mth2 
##  5:                           gi|443696809|gb|ELT97425.1| hypothetical protein CAPTEDRAFT_194915 
##  6:                                   gi|307176228|gb|EFN65864.1| hypothetical protein EAG_10145 
##  7:                    gi|332026123|gb|EGI66271.1| Multiple inositol polyphosphate phosphatase 1 
##  8:                                                          gi|15010456|gb|AAK77276.1| GH05807p 
##  9:             gi|332030037|gb|EGI69862.1| Serine/threonine-protein kinase PINK1, mitochondrial 
## 10:                            gi|307188496|gb|EFN73233.1| Muscarinic acetylcholine receptor DM1 
## 11:                      gi|332016397|gb|EGI57310.1| Mitochondrial import receptor subunit TOM70 
## 12:                                            gi|121605727|ref|YP_983056.1| OsmC family protein 
## 13:                                  gi|332029691|gb|EGI69570.1| G-protein coupled receptor Mth2 
## 14:         gi|380028536|ref|XP_003697954.1| PREDICTED: protein lethal(2)essential for life-like 
## 15: gi|332019420|gb|EGI59904.1| Putative fat-like cadherin-related tumor suppressor-like protein 
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
##  5:         High          Low
##  6:          Low          Low
##  7:          Low          Low
##  8:          Low          Low
##  9:          Low          Low
## 10:          Low          Low
## 11:          Low          Low
## 12: Intermediate          Low
## 13:          Low          Low
## 14:          Low Intermediate
## 15:          Low      Bimodal
## 16:          Low Intermediate
## 17:          Low Intermediate
## 18:          Low         High
## 19:          Low Intermediate
## 20:          Low Intermediate
## 21:          Low         High
## 22:          Low      Bimodal
## 23:          Low Intermediate
##         A22.type      Ar.type
```

```r
# in bimodal category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Bimodal", Ar.type)), 
    with(GO0006950.responsive, grep("Bimodal", A22.type))), ]
```

```
##                    Transcript
## 1:    14|*|comp150262_c0_seq1
## 2:  19778|*|comp97601_c0_seq1
## 3: 21598|*|comp142101_c0_seq1
## 4: 32312|*|comp933733_c0_seq1
## 5: 58246|*|comp109744_c0_seq1
## 6:  8886|*|comp150172_c1_seq4
## 7: 15115|*|comp132715_c0_seq1
## 8: 21384|*|comp149042_c0_seq3
##                                                                                   best.hit.to.nr
## 1: gi|332019420|gb|EGI59904.1| Putative fat-like cadherin-related tumor suppressor-like protein 
## 2:                                  gi|322796169|gb|EFZ18745.1| hypothetical protein SINV_07491 
## 3:                              gi|332018201|gb|EGI58806.1| Protein lethal(2)essential for life 
## 4:                          gi|367054010|ref|XP_003657383.1| hypothetical protein THITE_2156506 
## 5:                                    gi|493322437|ref|WP_006279741.1| molecular chaperone DnaK 
## 6:                                  gi|332029691|gb|EGI69570.1| G-protein coupled receptor Mth2 
## 7:                                            gi|194716766|gb|ACF93232.1| heat shock protein 90 
## 8:                         gi|396467618|ref|XP_003837992.1| hypothetical protein LEMA_P120390.1 
##    A22.type Ar.type
## 1:      Low Bimodal
## 2:  Bimodal Bimodal
## 3:     High Bimodal
## 4:     High Bimodal
## 5:  NotResp Bimodal
## 6:      Low Bimodal
## 7:  Bimodal    High
## 8:  Bimodal NotResp
```

```r
# in intermediate category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Intermediate", Ar.type)), 
    with(GO0006950.responsive, grep("Intermediate", A22.type))), ]
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

# Chi-square test to see if 'response to stress' related genes overrepresented in
# responsive.lms compared to full list
resp.stress.responsive.count <- nrow(responsive.lms.ann.type[grep("GO:0006950", responsive.lms.ann.type$GO.Biological.Process), 
    list(Transcript, best.hit.to.nr)])
# GO 'response to stress' hits in all transcripts
resp.stress.all.count <- nrow(annotationtable[grep("GO:0006950", annotationtable$GO.Biological.Process), 
    list(Sequence.Name, best.hit.to.nr)])

GO.stress.table <- matrix(rbind(resp.stress.responsive.count, nrow(responsive.lms.ann.type) - 
    resp.stress.responsive.count, resp.stress.all.count, nrow(annotationtable) - 
    resp.stress.all.count), nrow = 2)

GO.stress.Xsq <- chisq.test(GO.stress.table)
GO.stress.Xsq
```

```
## 
## 	Pearson's Chi-squared test with Yates' continuity correction
## 
## data:  GO.stress.table
## X-squared = 2.08, df = 1, p-value = 0.1489
```



There are also 7 heat shock related genes in the responsive transcripts, out of 130 total.


```r
hsp_responsive <- responsive.lms.ann.type[grep("shock", responsive.lms.ann.type$best.hit.to.nr), 
    list(Transcript, best.hit.to.nr, A22.type, Ar.type)]
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

hsp_all <- annotationtable[grep("shock", annotationtable$best.hit.to.nr), list(Sequence.Name, 
    best.hit.to.nr)]
nrow(hsp_all)
```

```
## [1] 130
```




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


### GSEA for thermally-responsive transcripts ###

Using this gene2GO map file, perform GSEA for:

**1) all responsive transcripts**

Use `selectFDR` function to select transcripts with adjusted P < 0.05.


```r
# create geneList. note that NA values cause problems with topGO so set any NA to
# 1 as need to retain for GO analysis
Ap.geneList <- RxNpval$padj
Ap.geneList[which(is.na(Ap.geneList))] <- 1
stopifnot(length(which(is.na(Ap.geneList))) == 0)
names(Ap.geneList) <- RxNpval$Transcript
str(Ap.geneList)
```

```
##  Named num [1:99861] 0.754 0.203 0.924 0.928 0.77 ...
##  - attr(*, "names")= chr [1:99861] "0|*|Contig6267" "100000|*|comp2663136_c0_seq1" "100001|*|comp3439067_c0_seq1" "100002|*|comp2050457_c0_seq1" ...
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
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
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
##  99861 available genes (all genes from the array):
##    - symbol:  0|*|Contig6267 100000|*|comp2663136_c0_seq1 100001|*|comp3439067_c0_seq1 100002|*|comp2050457_c0_seq1 100004|*|comp131141_c1_seq1  ...
```

```
## Error: invalid 'digits' argument
```

```r

# perform enrichment analysis using parentchild method
Ap.BP.resultParentChild <- runTest(Ap.BP.GOdata, statistic = "fisher", algorithm = "parentchild")
```

```
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 3119 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 18:	1 nodes to be scored.
## 
## 	 Level 17:	2 nodes to be scored.
## 
## 	 Level 16:	7 nodes to be scored.
## 
## 	 Level 15:	10 nodes to be scored.
## 
## 	 Level 14:	24 nodes to be scored.
## 
## 	 Level 13:	54 nodes to be scored.
## 
## 	 Level 12:	133 nodes to be scored.
## 
## 	 Level 11:	208 nodes to be scored.
## 
## 	 Level 10:	311 nodes to be scored.
## 
## 	 Level 9:	413 nodes to be scored.
## 
## 	 Level 8:	447 nodes to be scored.
## 
## 	 Level 7:	475 nodes to be scored.
## 
## 	 Level 6:	432 nodes to be scored.
## 
## 	 Level 5:	341 nodes to be scored.
## 
## 	 Level 4:	187 nodes to be scored.
## 
## 	 Level 3:	56 nodes to be scored.
## 
## 	 Level 2:	17 nodes to be scored.
```

```r
Ap.BP.resultParentChild
```

```
## 
## Description: BP gene set analysis 
## Ontology: BP 
## 'parentchild' algorithm with the 'fisher : joinFun = union' test
## 3534 GO terms scored: 125 terms with p < 0.01
## Annotation data:
##     Annotated genes: 30854 
##     Significant genes: 3053 
##     Min. no. of genes annotated to a GO: 10 
##     Nontrivial nodes: 3119
```

```r

# table results
Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, topNodes = 137)
dim(Ap.BP.ResTable)
```

```
## [1] 137   6
```

```r
# pandoc.table(Ap.BP.ResTable)

# graph significant nodes

# pdf(paste(resultsdir, 'Ap.BP_topGO_sig_nodes.pdf', sep=''))
# showSigOfNodes(Ap.BP.GOdata, score(Ap.BP.resultParentChild), firstSigNodes =
# 10, useInfo = 'all') dev.off()
```


**2) High**

Genes with *High* expression in both colonies

Use `selectTranscript` function to select transcripts from 'Ap.response.type' for GSEA.


```r
# Select transcripts
Ap.high <- Ap.response.type[which(Ap.response.type$Ar.type == "High" & Ap.response.type$A22.type == 
    "High"), "Transcript"]
Ap.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.high) <- names(Ap.geneList)
Ap.geneList.high[(which(names(Ap.geneList.high) %in% Ap.high))] <- 1
# check correct number of values set to 1
table(Ap.geneList.high)
```

```
## Ap.geneList.high
##     0     1 
## 99555   306
```

```r

# Run GSEA
Ap.high.gsea <- gsea(genelist = Ap.geneList.high, geneID2GO = geneID2GO, plotpath = NA)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 655 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 15:	1 nodes to be scored.
## 
## 	 Level 14:	2 nodes to be scored.
## 
## 	 Level 13:	3 nodes to be scored.
## 
## 	 Level 12:	7 nodes to be scored.
## 
## 	 Level 11:	22 nodes to be scored.
## 
## 	 Level 10:	46 nodes to be scored.
## 
## 	 Level 9:	66 nodes to be scored.
## 
## 	 Level 8:	74 nodes to be scored.
## 
## 	 Level 7:	89 nodes to be scored.
## 
## 	 Level 6:	103 nodes to be scored.
## 
## 	 Level 5:	110 nodes to be scored.
## 
## 	 Level 4:	84 nodes to be scored.
## 
## 	 Level 3:	32 nodes to be scored.
## 
## 	 Level 2:	15 nodes to be scored.
```

```r
# pandoc.table(Ap.high.gsea) output for Revigo
write.table(Ap.high.gsea[, "GO.ID"], file = "results/Ap_high_gsea.txt", row.names = FALSE, 
    sep = "\t", quote = FALSE)
```



**3) Low**

Genes with *Low* expression in both colonies


```r
# Select transcripts
Ap.low <- Ap.response.type[which(Ap.response.type$Ar.type == "Low" & Ap.response.type$A22.type == 
    "Low"), "Transcript"]
Ap.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.low) <- names(Ap.geneList)
Ap.geneList.low[(which(names(Ap.geneList.low) %in% Ap.low))] <- 1
# check correct number of values set to 1
table(Ap.geneList.low)
```

```
## Ap.geneList.low
##     0     1 
## 97229  2632
```

```r

# Run GSEA
Ap.low.gsea <- gsea(genelist = Ap.geneList.low, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 2338 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	2 nodes to be scored.
## 
## 	 Level 16:	5 nodes to be scored.
## 
## 	 Level 15:	6 nodes to be scored.
## 
## 	 Level 14:	13 nodes to be scored.
## 
## 	 Level 13:	32 nodes to be scored.
## 
## 	 Level 12:	87 nodes to be scored.
## 
## 	 Level 11:	146 nodes to be scored.
## 
## 	 Level 10:	216 nodes to be scored.
## 
## 	 Level 9:	300 nodes to be scored.
## 
## 	 Level 8:	333 nodes to be scored.
## 
## 	 Level 7:	358 nodes to be scored.
## 
## 	 Level 6:	324 nodes to be scored.
## 
## 	 Level 5:	287 nodes to be scored.
## 
## 	 Level 4:	162 nodes to be scored.
## 
## 	 Level 3:	49 nodes to be scored.
## 
## 	 Level 2:	17 nodes to be scored.
```

```r
# pandoc.table(Ap.low.BP.ResTable)
write.table(Ap.low.gsea[, "GO.ID"], file = "results/Ap_low_gsea.txt", row.names = FALSE, 
    sep = "\t", quote = FALSE)
```



**4) Bimodal**

Genes with *Bimodal* expression in both colonies


```r
Ap.bim <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal" & Ap.response.type$Ar.type == 
    "Bimodal"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
Ap.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.bim) <- names(Ap.geneList)
Ap.geneList.bim[(which(names(Ap.geneList.bim) %in% Ap.bim))] <- 1
# check correct number of values set to 1
table(Ap.geneList.bim)
```

```
## Ap.geneList.bim
##     0     1 
## 99570   291
```

```r

# Run GSEA
Ap.bim.gsea <- gsea(genelist = Ap.geneList.bim, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 607 nontrivial nodes
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
## 	 Level 12:	9 nodes to be scored.
## 
## 	 Level 11:	20 nodes to be scored.
## 
## 	 Level 10:	47 nodes to be scored.
## 
## 	 Level 9:	65 nodes to be scored.
## 
## 	 Level 8:	72 nodes to be scored.
## 
## 	 Level 7:	77 nodes to be scored.
## 
## 	 Level 6:	89 nodes to be scored.
## 
## 	 Level 5:	101 nodes to be scored.
## 
## 	 Level 4:	72 nodes to be scored.
## 
## 	 Level 3:	31 nodes to be scored.
## 
## 	 Level 2:	14 nodes to be scored.
```

```r
# pandoc.table(Ap.bim.gsea)
write.table(Ap.bim.gsea[, "GO.ID"], file = "results/Ap_bim_gsea.txt", row.names = FALSE, 
    sep = "\t", quote = FALSE)
```



**5) Intermediate**

Genes with *Intermediate* expression in both colonies


```r
Ap.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & Ap.response.type$Ar.type == 
    "Intermediate"), "Transcript"]
# create gene list, setting value to 1 for 'int' transcripts
Ap.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.int) <- names(Ap.geneList)
Ap.geneList.int[(which(names(Ap.geneList.int) %in% Ap.int))] <- 1
# check correct number of values set to 1
table(Ap.geneList.int)
```

```
## Ap.geneList.int
##     0     1 
## 99292   569
```

```r

# Run GSEA
Ap.int.gsea <- gsea(genelist = Ap.geneList.int, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1341 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 18:	1 nodes to be scored.
## 
## 	 Level 17:	1 nodes to be scored.
## 
## 	 Level 16:	2 nodes to be scored.
## 
## 	 Level 15:	3 nodes to be scored.
## 
## 	 Level 14:	3 nodes to be scored.
## 
## 	 Level 13:	11 nodes to be scored.
## 
## 	 Level 12:	29 nodes to be scored.
## 
## 	 Level 11:	68 nodes to be scored.
## 
## 	 Level 10:	112 nodes to be scored.
## 
## 	 Level 9:	134 nodes to be scored.
## 
## 	 Level 8:	160 nodes to be scored.
## 
## 	 Level 7:	202 nodes to be scored.
## 
## 	 Level 6:	217 nodes to be scored.
## 
## 	 Level 5:	205 nodes to be scored.
## 
## 	 Level 4:	129 nodes to be scored.
## 
## 	 Level 3:	46 nodes to be scored.
## 
## 	 Level 2:	17 nodes to be scored.
```

```r
# pandoc.table(Ap.int.gsea)
write.table(Ap.int.gsea[, "GO.ID"], file = "results/Ap_int_gsea.txt", row.names = FALSE, 
    sep = "\t", quote = FALSE)
```


            Type      GO.ID                                        Term
1           High GO:0009889          regulation of biosynthetic process
2           High GO:0031326 regulation of cellular biosynthetic proc...
3           High GO:0019222             regulation of metabolic process
4           High GO:0031323    regulation of cellular metabolic process
5           High GO:0050794              regulation of cellular process
6           High GO:0050789            regulation of biological process
7           High GO:0080090     regulation of primary metabolic process
8           High GO:0010556 regulation of macromolecule biosynthetic...
9            Low GO:0019538                   protein metabolic process
10           Low GO:0070085                               glycosylation
11           Low GO:0043412                  macromolecule modification
12           Low GO:0044267          cellular protein metabolic process
13           Low GO:0007264 small GTPase mediated signal transductio...
14           Low GO:0009141 nucleoside triphosphate metabolic proces...
15           Low GO:0006643            membrane lipid metabolic process
16           Low GO:0033036                  macromolecule localization
17           Low GO:0045184       establishment of protein localization
18           Low GO:0009101           glycoprotein biosynthetic process
19           Low GO:0006486                       protein glycosylation
20           Low GO:0036211                protein modification process
21           Low GO:0009581              detection of external stimulus
22           Low GO:0009100              glycoprotein metabolic process
23           Low GO:0009259            ribonucleotide metabolic process
24           Low GO:0055001                     muscle cell development
25           Low GO:0032318           regulation of Ras GTPase activity
26           Low GO:0050650 chondroitin sulfate proteoglycan biosynt...
27           Low GO:0030206    chondroitin sulfate biosynthetic process
28           Low GO:0046907                     intracellular transport
29           Low GO:1901292      nucleoside phosphate catabolic process
30           Low GO:0051259                     protein oligomerization
31           Low GO:0070271                  protein complex biogenesis
32           Low GO:0009966           regulation of signal transduction
33           Low GO:0006665              sphingolipid metabolic process
34           Low GO:0006140 regulation of nucleotide metabolic proce...
35           Low GO:0006464       cellular protein modification process
36           Low GO:0061326                    renal tubule development
37           Low GO:0043170             macromolecule metabolic process
38           Low GO:0009118 regulation of nucleoside metabolic proce...
39           Low GO:0009894             regulation of catabolic process
40           Low GO:0010608 posttranscriptional regulation of gene e...
41           Low GO:0046039                       GTP metabolic process
42           Low GO:0009166                nucleotide catabolic process
43           Low GO:0006163         purine nucleotide metabolic process
44           Low GO:0009582               detection of abiotic stimulus
45           Low GO:0031329    regulation of cellular catabolic process
46           Low GO:0006413                    translational initiation
47           Low GO:0042692                 muscle cell differentiation
48           Low GO:0006749               glutathione metabolic process
49           Low GO:0006518                   peptide metabolic process
50           Low GO:0050954 sensory perception of mechanical stimulu...
51           Low GO:1901068 guanosine-containing compound metabolic ...
52           Low GO:0030811 regulation of nucleotide catabolic proce...
53           Low GO:0006417                   regulation of translation
54           Low GO:0051649       establishment of localization in cell
55           Low GO:0001655               urogenital system development
56           Low GO:0033121 regulation of purine nucleotide cataboli...
57           Low GO:0015031                           protein transport
58           Low GO:0048871        multicellular organismal homeostasis
59           Low GO:0006664                glycolipid metabolic process
60           Low GO:0032273 positive regulation of protein polymeriz...
61           Low GO:0019220 regulation of phosphate metabolic proces...
62           Low GO:0019685               photosynthesis, dark reaction
63           Low GO:0065003             macromolecular complex assembly
64           Low GO:0016573                         histone acetylation
65           Low GO:0072002               Malpighian tubule development
66           Low GO:0050654 chondroitin sulfate proteoglycan metabol...
67           Low GO:0030204       chondroitin sulfate metabolic process
68           Low GO:0006108                    malate metabolic process
69           Low GO:0060042     retina morphogenesis in camera-type eye
70           Low GO:0050794              regulation of cellular process
71           Low GO:0033559    unsaturated fatty acid metabolic process
72           Low GO:1900542 regulation of purine nucleotide metaboli...
73       Bimodal GO:0065007                       biological regulation
74       Bimodal GO:0050794              regulation of cellular process
75       Bimodal GO:0050789            regulation of biological process
76       Bimodal GO:0009889          regulation of biosynthetic process
77       Bimodal GO:0031323    regulation of cellular metabolic process
78       Bimodal GO:0080090     regulation of primary metabolic process
79       Bimodal GO:0031326 regulation of cellular biosynthetic proc...
80       Bimodal GO:0019222             regulation of metabolic process
81       Bimodal GO:0010468               regulation of gene expression
82       Bimodal GO:0006351                transcription, DNA-templated
83       Bimodal GO:0010556 regulation of macromolecule biosynthetic...
84       Bimodal GO:0051171 regulation of nitrogen compound metaboli...
85       Bimodal GO:2000112 regulation of cellular macromolecule bio...
86       Bimodal GO:0060255 regulation of macromolecule metabolic pr...
87       Bimodal GO:0032774                    RNA biosynthetic process
88       Bimodal GO:0050896                        response to stimulus
89       Bimodal GO:0044700                   single organism signaling
90       Bimodal GO:0034654 nucleobase-containing compound biosynthe...
91       Bimodal GO:0019219 regulation of nucleobase-containing comp...
92       Bimodal GO:0009059          macromolecule biosynthetic process
93       Bimodal GO:0036211                protein modification process
94       Bimodal GO:0023052                                   signaling
95  Intermediate GO:0006664                glycolipid metabolic process
96  Intermediate GO:0046467         membrane lipid biosynthetic process
97  Intermediate GO:0009100              glycoprotein metabolic process
98  Intermediate GO:0009101           glycoprotein biosynthetic process
99  Intermediate GO:0006486                       protein glycosylation
100 Intermediate GO:0009247             glycolipid biosynthetic process
101 Intermediate GO:0006643            membrane lipid metabolic process
102 Intermediate GO:0006323                               DNA packaging
103 Intermediate GO:0005975              carbohydrate metabolic process
104 Intermediate GO:0071103                     DNA conformation change
105 Intermediate GO:0034728                     nucleosome organization
106 Intermediate GO:0043413                 macromolecule glycosylation
107 Intermediate GO:0044723 single-organism carbohydrate metabolic p...
108 Intermediate GO:0006665              sphingolipid metabolic process
109 Intermediate GO:0071824    protein-DNA complex subunit organization
110 Intermediate GO:0022900                    electron transport chain
111 Intermediate GO:0042158            lipoprotein biosynthetic process
    Annotated Significant Expected        P
1        3430          19     7.89  0.00061
2        3425          19     7.88  0.00116
3        4963          22    11.42  0.00282
4        4500          21    10.36  0.00340
5        8411          31    19.36  0.00462
6        8921          31    20.53  0.00567
7        4480          20    10.31  0.00851
8        3390          19     7.80  0.00858
9        7596         264   203.35  2.0e-07
10        149          12     3.99  1.7e-06
11       3268         136    87.49  3.0e-06
12       5294         196   141.73  2.0e-05
13        659          37    17.64  5.2e-05
14       1980          61    53.01  8.7e-05
15        242          16     6.48  0.00021
16       1636          61    43.80  0.00022
17       1092          46    29.23  0.00023
18        165          13     4.42  0.00025
19        141          12     3.77  0.00026
20       2887         127    77.29  0.00028
21         67           6     1.79  0.00055
22        176          14     4.71  0.00056
23       2256          65    60.40  0.00056
24         51           6     1.37  0.00062
25        149          14     3.99  0.00101
26         16           3     0.43  0.00108
27         16           3     0.43  0.00109
28        960          39    25.70  0.00112
29       1769          53    47.36  0.00116
30         97          11     2.60  0.00129
31        561          25    15.02  0.00133
32        892          39    23.88  0.00151
33        165          12     4.42  0.00173
34        521          24    13.95  0.00187
35       2887         127    77.29  0.00192
36         43           5     1.15  0.00210
37      13944         409   373.30  0.00228
38        500          24    13.39  0.00232
39        615          28    16.46  0.00274
40        862          35    23.08  0.00276
41        919          38    24.60  0.00350
42       1765          53    47.25  0.00352
43       2180          65    58.36  0.00357
44         67           6     1.79  0.00427
45        572          26    15.31  0.00430
46        370          19     9.91  0.00450
47         71           6     1.90  0.00455
48        204          12     5.46  0.00477
49        265          15     7.09  0.00509
50         56           5     1.50  0.00517
51        953          39    25.51  0.00521
52        500          24    13.39  0.00524
53        747          30    20.00  0.00559
54       1323          47    35.42  0.00583
55         88           7     2.36  0.00610
56        500          24    13.39  0.00621
57       1059          43    28.35  0.00634
58         48           5     1.29  0.00637
59        211          13     5.65  0.00672
60         38           4     1.02  0.00714
61        821          34    21.98  0.00717
62         13           2     0.35  0.00727
63        731          30    19.57  0.00772
64        136          12     3.64  0.00810
65         33           4     0.88  0.00827
66         19           3     0.51  0.00858
67         19           3     0.51  0.00859
68         45           4     1.20  0.00907
69         18           3     0.48  0.00909
70       8411         249   225.17  0.00910
71         48           3     1.29  0.00922
72        514          24    13.76  0.00944
73       9263          36    19.21 0.000011
74       8411          34    17.45 0.000035
75       8921          34    18.50 0.000041
76       3430          19     7.11  0.00011
77       4500          23     9.33  0.00011
78       4480          22     9.29  0.00019
79       3425          19     7.10  0.00037
80       4963          23    10.29  0.00042
81       3698          19     7.67  0.00073
82       2978          16     6.18  0.00104
83       3390          19     7.03  0.00106
84       3560          18     7.38  0.00119
85       3388          19     7.03  0.00204
86       4193          20     8.70  0.00262
87       3000          16     6.22  0.00443
88       7343          25    15.23  0.00465
89       4891          18    10.15  0.00483
90       3887          17     8.06  0.00500
91       3544          18     7.35  0.00597
92       7008          24    14.54  0.00637
93       2887          12     5.99  0.00861
94       4892          18    10.15  0.00902
95        211           9     1.53  3.4e-06
96         81           5     0.59  0.00012
97        176           7     1.28  0.00015
98        165           6     1.20  0.00041
99        141           6     1.02  0.00055
100        72           5     0.52  0.00060
101       242           9     1.76  0.00082
102       246           6     1.79  0.00145
103      2084          26    15.13  0.00180
104       431           7     3.13  0.00186
105       151           4     1.10  0.00198
106       141           6     1.02  0.00205
107      1500          21    10.89  0.00238
108       165           5     1.20  0.00354
109       164           4     1.19  0.00510
110       394           8     2.86  0.00898
111        71           3     0.52  0.00928



Next, I perform GSEA for genes in each functional type in one colony but not the other (e.g. the set difference) to gain insight on differences between the colonies.

**A22 'High' genes not in Ar**


```r
A22.high <- Ap.response.type[which(Ap.response.type$A22.type == "High" & Ap.response.type$Ar.type != 
    "High"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.high <- rep(0, length = length(Ap.geneList))
names(A22.geneList.high) <- names(Ap.geneList)
A22.geneList.high[(which(names(A22.geneList.high) %in% A22.high))] <- 1
# check correct number of values set to 1
table(A22.geneList.high)
```

```
## A22.geneList.high
##     0     1 
## 98902   959
```

```r

# Run GSEA
A22.high.gsea <- gsea(genelist = A22.geneList.high, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1095 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	1 nodes to be scored.
## 
## 	 Level 16:	2 nodes to be scored.
## 
## 	 Level 15:	2 nodes to be scored.
## 
## 	 Level 14:	3 nodes to be scored.
## 
## 	 Level 13:	5 nodes to be scored.
## 
## 	 Level 12:	19 nodes to be scored.
## 
## 	 Level 11:	43 nodes to be scored.
## 
## 	 Level 10:	73 nodes to be scored.
## 
## 	 Level 9:	110 nodes to be scored.
## 
## 	 Level 8:	136 nodes to be scored.
## 
## 	 Level 7:	163 nodes to be scored.
## 
## 	 Level 6:	183 nodes to be scored.
## 
## 	 Level 5:	176 nodes to be scored.
## 
## 	 Level 4:	117 nodes to be scored.
## 
## 	 Level 3:	45 nodes to be scored.
## 
## 	 Level 2:	16 nodes to be scored.
```

```r
# pandoc.table(A22.high.gsea, split.table = Inf)
```



**Ar 'High' genes not in A22**


```r
Ar.high <- Ap.response.type[which(Ap.response.type$A22.type != "High" & Ap.response.type$Ar.type == 
    "High"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.high) <- names(Ap.geneList)
Ar.geneList.high[(which(names(Ar.geneList.high) %in% Ar.high))] <- 1
# check correct number of values set to 1
table(Ar.geneList.high)
```

```
## Ar.geneList.high
##     0     1 
## 99108   753
```

```r

# Run GSEA
Ar.high.gsea <- gsea(genelist = Ar.geneList.high, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1068 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	1 nodes to be scored.
## 
## 	 Level 16:	3 nodes to be scored.
## 
## 	 Level 15:	2 nodes to be scored.
## 
## 	 Level 14:	4 nodes to be scored.
## 
## 	 Level 13:	12 nodes to be scored.
## 
## 	 Level 12:	24 nodes to be scored.
## 
## 	 Level 11:	62 nodes to be scored.
## 
## 	 Level 10:	92 nodes to be scored.
## 
## 	 Level 9:	104 nodes to be scored.
## 
## 	 Level 8:	122 nodes to be scored.
## 
## 	 Level 7:	145 nodes to be scored.
## 
## 	 Level 6:	167 nodes to be scored.
## 
## 	 Level 5:	164 nodes to be scored.
## 
## 	 Level 4:	107 nodes to be scored.
## 
## 	 Level 3:	42 nodes to be scored.
## 
## 	 Level 2:	16 nodes to be scored.
```

```r
# pandoc.table(Ar.high.gsea, split.table = Inf)
```



**A22 'Low' genes not in Ar**


```r
A22.low <- Ap.response.type[which(Ap.response.type$A22.type == "Low" & Ap.response.type$Ar.type != 
    "Low"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.low <- rep(0, length = length(Ap.geneList))
names(A22.geneList.low) <- names(Ap.geneList)
A22.geneList.low[(which(names(A22.geneList.low) %in% A22.low))] <- 1
# check correct number of values set to 1
table(A22.geneList.low)
```

```
## A22.geneList.low
##     0     1 
## 97576  2285
```

```r

# Run GSEA
A22.low.gsea <- gsea(genelist = A22.geneList.low, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 2023 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	1 nodes to be scored.
## 
## 	 Level 16:	4 nodes to be scored.
## 
## 	 Level 15:	4 nodes to be scored.
## 
## 	 Level 14:	8 nodes to be scored.
## 
## 	 Level 13:	25 nodes to be scored.
## 
## 	 Level 12:	65 nodes to be scored.
## 
## 	 Level 11:	125 nodes to be scored.
## 
## 	 Level 10:	191 nodes to be scored.
## 
## 	 Level 9:	235 nodes to be scored.
## 
## 	 Level 8:	276 nodes to be scored.
## 
## 	 Level 7:	313 nodes to be scored.
## 
## 	 Level 6:	304 nodes to be scored.
## 
## 	 Level 5:	256 nodes to be scored.
## 
## 	 Level 4:	150 nodes to be scored.
## 
## 	 Level 3:	48 nodes to be scored.
## 
## 	 Level 2:	17 nodes to be scored.
```

```r
# pandoc.table(A22.low.gsea, split.table = Inf)
```



**Ar 'Low' genes not in A22**


```r
Ar.low <- Ap.response.type[which(Ap.response.type$A22.type != "Low" & Ap.response.type$Ar.type == 
    "Low"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.low) <- names(Ap.geneList)
Ar.geneList.low[(which(names(Ar.geneList.low) %in% Ar.low))] <- 1
# check correct number of values set to 1
table(Ar.geneList.low)
```

```
## Ar.geneList.low
##     0     1 
## 98730  1131
```

```r

# Run GSEA
Ar.low.gsea <- gsea(genelist = Ar.geneList.low, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1312 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	1 nodes to be scored.
## 
## 	 Level 16:	2 nodes to be scored.
## 
## 	 Level 15:	2 nodes to be scored.
## 
## 	 Level 14:	5 nodes to be scored.
## 
## 	 Level 13:	12 nodes to be scored.
## 
## 	 Level 12:	40 nodes to be scored.
## 
## 	 Level 11:	72 nodes to be scored.
## 
## 	 Level 10:	98 nodes to be scored.
## 
## 	 Level 9:	144 nodes to be scored.
## 
## 	 Level 8:	162 nodes to be scored.
## 
## 	 Level 7:	187 nodes to be scored.
## 
## 	 Level 6:	210 nodes to be scored.
## 
## 	 Level 5:	185 nodes to be scored.
## 
## 	 Level 4:	129 nodes to be scored.
## 
## 	 Level 3:	46 nodes to be scored.
## 
## 	 Level 2:	16 nodes to be scored.
```

```r
# pandoc.table(Ar.low.gsea, split.table = Inf)
```




**A22 'Bimodal' genes not in Ar**


```r
A22.bim <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal" & Ap.response.type$Ar.type != 
    "Bimodal"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.bim <- rep(0, length = length(Ap.geneList))
names(A22.geneList.bim) <- names(Ap.geneList)
A22.geneList.bim[(which(names(A22.geneList.bim) %in% A22.bim))] <- 1
# check correct number of values set to 1
table(A22.geneList.bim)
```

```
## A22.geneList.bim
##     0     1 
## 98620  1241
```

```r

# Run GSEA
A22.bim.gsea <- gsea(genelist = A22.geneList.bim, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1417 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	1 nodes to be scored.
## 
## 	 Level 16:	2 nodes to be scored.
## 
## 	 Level 15:	2 nodes to be scored.
## 
## 	 Level 14:	7 nodes to be scored.
## 
## 	 Level 13:	21 nodes to be scored.
## 
## 	 Level 12:	47 nodes to be scored.
## 
## 	 Level 11:	85 nodes to be scored.
## 
## 	 Level 10:	122 nodes to be scored.
## 
## 	 Level 9:	153 nodes to be scored.
## 
## 	 Level 8:	170 nodes to be scored.
## 
## 	 Level 7:	201 nodes to be scored.
## 
## 	 Level 6:	210 nodes to be scored.
## 
## 	 Level 5:	205 nodes to be scored.
## 
## 	 Level 4:	128 nodes to be scored.
## 
## 	 Level 3:	45 nodes to be scored.
## 
## 	 Level 2:	17 nodes to be scored.
```

```r
# pandoc.table(A22.bim.gsea, split.table = Inf)
```



**Ar 'Bimodal' genes not in A22**


```r
Ar.bim <- Ap.response.type[which(Ap.response.type$A22.type != "Bimodal" & Ap.response.type$Ar.type == 
    "Bimodal"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.bim) <- names(Ap.geneList)
Ar.geneList.bim[(which(names(Ar.geneList.bim) %in% Ar.bim))] <- 1
# check correct number of values set to 1
table(Ar.geneList.bim)
```

```
## Ar.geneList.bim
##     0     1 
## 99302   559
```

```r

# Run GSEA
Ar.bim.gsea <- gsea(genelist = Ar.geneList.bim, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1073 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 16:	2 nodes to be scored.
## 
## 	 Level 15:	3 nodes to be scored.
## 
## 	 Level 14:	7 nodes to be scored.
## 
## 	 Level 13:	14 nodes to be scored.
## 
## 	 Level 12:	21 nodes to be scored.
## 
## 	 Level 11:	44 nodes to be scored.
## 
## 	 Level 10:	81 nodes to be scored.
## 
## 	 Level 9:	118 nodes to be scored.
## 
## 	 Level 8:	130 nodes to be scored.
## 
## 	 Level 7:	155 nodes to be scored.
## 
## 	 Level 6:	166 nodes to be scored.
## 
## 	 Level 5:	161 nodes to be scored.
## 
## 	 Level 4:	114 nodes to be scored.
## 
## 	 Level 3:	39 nodes to be scored.
## 
## 	 Level 2:	17 nodes to be scored.
```

```r
# pandoc.table(Ar.bim.gsea, split.table = Inf)
```




**A22 'Intermediate' genes not in Ar**


```r
A22.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & Ap.response.type$Ar.type != 
    "Intermediate"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.int <- rep(0, length = length(Ap.geneList))
names(A22.geneList.int) <- names(Ap.geneList)
A22.geneList.int[(which(names(A22.geneList.int) %in% A22.int))] <- 1
# check correct number of values set to 1
table(A22.geneList.int)
```

```
## A22.geneList.int
##     0     1 
## 99502   359
```

```r

# Run GSEA
A22.int.gsea <- gsea(genelist = A22.geneList.int, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 783 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 16:	2 nodes to be scored.
## 
## 	 Level 15:	3 nodes to be scored.
## 
## 	 Level 14:	5 nodes to be scored.
## 
## 	 Level 13:	9 nodes to be scored.
## 
## 	 Level 12:	17 nodes to be scored.
## 
## 	 Level 11:	37 nodes to be scored.
## 
## 	 Level 10:	63 nodes to be scored.
## 
## 	 Level 9:	85 nodes to be scored.
## 
## 	 Level 8:	88 nodes to be scored.
## 
## 	 Level 7:	105 nodes to be scored.
## 
## 	 Level 6:	111 nodes to be scored.
## 
## 	 Level 5:	123 nodes to be scored.
## 
## 	 Level 4:	85 nodes to be scored.
## 
## 	 Level 3:	35 nodes to be scored.
## 
## 	 Level 2:	14 nodes to be scored.
```

```r
# pandoc.table(A22.int.gsea, split.table = Inf)
```



**Ar 'Intermediate' genes not in A22**


```r
Ar.int <- Ap.response.type[which(Ap.response.type$A22.type != "Intermediate" & Ap.response.type$Ar.type == 
    "Intermediate"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.int) <- names(Ap.geneList)
Ar.geneList.int[(which(names(Ar.geneList.int) %in% Ar.int))] <- 1
# check correct number of values set to 1
table(Ar.geneList.int)
```

```
## Ar.geneList.int
##     0     1 
## 97700  2161
```

```r

# Run GSEA
Ar.int.gsea <- gsea(genelist = Ar.geneList.int, geneID2GO = geneID2GO)
```

```
## 
## Building most specific GOs .....	( 5471 GO terms found. )
## 
## Build GO DAG topology ..........	( 8953 GO terms and 19938 relations. )
## 
## Annotating nodes ...............	( 30854 genes annotated to the GO terms. )
## 
## 			 -- Parent-Child Algorithm -- 
## 
## 		 the algorithm is scoring 1940 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 17:	1 nodes to be scored.
## 
## 	 Level 16:	3 nodes to be scored.
## 
## 	 Level 15:	4 nodes to be scored.
## 
## 	 Level 14:	8 nodes to be scored.
## 
## 	 Level 13:	26 nodes to be scored.
## 
## 	 Level 12:	67 nodes to be scored.
## 
## 	 Level 11:	116 nodes to be scored.
## 
## 	 Level 10:	177 nodes to be scored.
## 
## 	 Level 9:	231 nodes to be scored.
## 
## 	 Level 8:	261 nodes to be scored.
## 
## 	 Level 7:	295 nodes to be scored.
## 
## 	 Level 6:	291 nodes to be scored.
## 
## 	 Level 5:	252 nodes to be scored.
## 
## 	 Level 4:	145 nodes to be scored.
## 
## 	 Level 3:	45 nodes to be scored.
## 
## 	 Level 2:	17 nodes to be scored.
```

```r
# pandoc.table(Ar.int.gsea, split.table = Inf)
```



Table of overall GSEA results.

% latex table generated in R 3.1.0 by xtable 1.7-3 package
% Thu May  8 15:52:12 2014
\begin{table}[ht]
\centering
\begin{tabular}{rllllrrrl}
  \hline
 & Colony & Type & GO.ID & Term & Annotated & Significant & Expected & P \\ 
  \hline
1 & ApVT & High & GO:1902222 & erythrose 4-phosphate/phosphoenolpyruvat... &  18 &   2 & 0.13 & 0.00092 \\ 
  2 & ApVT & High & GO:0030435 & sporulation resulting in formation of a ... &  18 &   2 & 0.13 & 0.00293 \\ 
  3 & ApVT & High & GO:0006559 & L-phenylalanine catabolic process &  18 &   2 & 0.13 & 0.00305 \\ 
  4 & ApVT & High & GO:0006412 & translation & 1886 &  25 & 13.39 & 0.00481 \\ 
  5 & ApVT & High & GO:0033002 & muscle cell proliferation &  12 &   2 & 0.09 & 0.00482 \\ 
  6 & ApVT & High & GO:0015074 & DNA integration & 852 &  17 & 6.05 & 0.00503 \\ 
  7 & ApVT & High & GO:0051049 & regulation of transport & 269 &   5 & 1.91 & 0.00653 \\ 
  8 & ApVT & High & GO:0043269 & regulation of ion transport &  81 &   3 & 0.57 & 0.00676 \\ 
  9 & ApVT & High & GO:0006424 & glutamyl-tRNA aminoacylation &  31 &   2 & 0.22 & 0.00857 \\ 
  10 & ApVT & High & GO:0043934 & sporulation &  23 &   2 & 0.16 & 0.00898 \\ 
  11 & ApVT & Low & GO:1902589 & single-organism organelle organization & 1407 &  58 & 32.38 & 7.4e-06 \\ 
  12 & ApVT & Low & GO:0019538 & protein metabolic process & 7596 & 217 & 174.80 & 1.3e-05 \\ 
  13 & ApVT & Low & GO:0006996 & organelle organization & 2042 &  78 & 46.99 & 1.7e-05 \\ 
  14 & ApVT & Low & GO:0043170 & macromolecule metabolic process & 13944 & 350 & 320.87 & 7.7e-05 \\ 
  15 & ApVT & Low & GO:0022406 & membrane docking &  34 &   6 & 0.78 & 0.00010 \\ 
  16 & ApVT & Low & GO:0044260 & cellular macromolecule metabolic process & 12209 & 308 & 280.95 & 0.00039 \\ 
  17 & ApVT & Low & GO:0048278 & vesicle docking &  33 &   6 & 0.76 & 0.00058 \\ 
  18 & ApVT & Low & GO:0044267 & cellular protein metabolic process & 5294 & 161 & 121.82 & 0.00060 \\ 
  19 & ApVT & Low & GO:0051641 & cellular localization & 1501 &  51 & 34.54 & 0.00068 \\ 
  20 & ApVT & Low & GO:0022904 & respiratory electron transport chain & 368 &  12 & 8.47 & 0.00104 \\ 
  21 & ApVT & Low & GO:0000725 & recombinational repair &  48 &   5 & 1.10 & 0.00127 \\ 
  22 & ApVT & Low & GO:0046907 & intracellular transport & 960 &  35 & 22.09 & 0.00205 \\ 
  23 & ApVT & Low & GO:0044707 & single-multicellular organism process & 3354 &  98 & 77.18 & 0.00262 \\ 
  24 & ApVT & Low & GO:0051649 & establishment of localization in cell & 1323 &  45 & 30.44 & 0.00266 \\ 
  25 & ApVT & Low & GO:0007264 & small GTPase mediated signal transductio... & 659 &  25 & 15.16 & 0.00292 \\ 
  26 & ApVT & Low & GO:0009069 & serine family amino acid metabolic proce... & 1251 &  33 & 28.79 & 0.00294 \\ 
  27 & ApVT & Low & GO:0006413 & translational initiation & 370 &  17 & 8.51 & 0.00372 \\ 
  28 & ApVT & Low & GO:0006446 & regulation of translational initiation & 344 &  17 & 7.92 & 0.00374 \\ 
  29 & ApVT & Low & GO:0070085 & glycosylation & 149 &   7 & 3.43 & 0.00392 \\ 
  30 & ApVT & Low & GO:1900542 & regulation of purine nucleotide metaboli... & 514 &  17 & 11.83 & 0.00413 \\ 
  31 & ApVT & Low & GO:0046039 & GTP metabolic process & 919 &  25 & 21.15 & 0.00571 \\ 
  32 & ApVT & Low & GO:0006278 & RNA-dependent DNA replication & 725 &  29 & 16.68 & 0.00599 \\ 
  33 & ApVT & Low & GO:0035315 & hair cell differentiation &  41 &   4 & 0.94 & 0.00716 \\ 
  34 & ApVT & Low & GO:0009166 & nucleotide catabolic process & 1765 &  36 & 40.62 & 0.00806 \\ 
  35 & ApVT & Low & GO:0006184 & GTP catabolic process & 864 &  25 & 19.88 & 0.00832 \\ 
  36 & ApVT & Low & GO:0042439 & ethanolamine-containing compound metabol... &  36 &   4 & 0.83 & 0.00849 \\ 
  37 & ApVT & Low & GO:0006644 & phospholipid metabolic process & 500 &  19 & 11.51 & 0.00857 \\ 
  38 & ApVT & Low & GO:0006904 & vesicle docking involved in exocytosis &  29 &   5 & 0.67 & 0.00865 \\ 
  39 & ApVT & Low & GO:0042726 & flavin-containing compound metabolic pro... &  88 &   6 & 2.03 & 0.00940 \\ 
  40 & ApVT & Low & GO:0035966 & response to topologically incorrect prot... &  58 &   5 & 1.33 & 0.00999 \\ 
  41 & ApVT & Intermediate & GO:0009561 & megagametogenesis &  18 &   2 & 0.06 & 0.00022 \\ 
  42 & ApVT & Intermediate & GO:0009553 & embryo sac development &  20 &   2 & 0.06 & 0.00032 \\ 
  43 & ApVT & Intermediate & GO:0009560 & embryo sac egg cell differentiation &  14 &   2 & 0.04 & 0.00034 \\ 
  44 & ApVT & Intermediate & GO:0048229 & gametophyte development &  28 &   2 & 0.09 & 0.00067 \\ 
  45 & ApVT & Intermediate & GO:0001558 & regulation of cell growth &  74 &   3 & 0.24 & 0.00178 \\ 
  46 & ApVT & Intermediate & GO:0030031 & cell projection assembly & 124 &   3 & 0.39 & 0.00498 \\ 
  47 & ApVT & Intermediate & GO:0015942 & formate metabolic process &  47 &   2 & 0.15 & 0.00499 \\ 
  48 & ApVT & Intermediate & GO:0042126 & nitrate metabolic process &  51 &   2 & 0.16 & 0.00598 \\ 
  49 & ApVT & Intermediate & GO:0006259 & DNA metabolic process & 2944 &  16 & 9.35 & 0.00794 \\ 
  50 & ApVT & Bimodal & GO:0006546 & glycine catabolic process &  53 &   4 & 0.50 & 0.00035 \\ 
  51 & ApVT & Bimodal & GO:0009071 & serine family amino acid catabolic proce... &  56 &   4 & 0.53 & 0.00391 \\ 
  52 & ApVT & Bimodal & GO:0006468 & protein phosphorylation & 1158 &  21 & 11.00 & 0.00469 \\ 
  53 & ApVT & Bimodal & GO:0007215 & glutamate receptor signaling pathway &  88 &   4 & 0.84 & 0.00827 \\ 
  54 & ApVT & Bimodal & GO:0006886 & intracellular protein transport & 667 &  12 & 6.33 & 0.00914 \\ 
  55 & ApNC & High & GO:0001503 & ossification &  37 &   3 & 0.23 & 0.00063 \\ 
  56 & ApNC & High & GO:0090304 & nucleic acid metabolic process & 7050 &  60 & 43.87 & 0.00083 \\ 
  57 & ApNC & High & GO:0006259 & DNA metabolic process & 2944 &  32 & 18.32 & 0.00181 \\ 
  58 & ApNC & High & GO:0006457 & protein folding & 460 &   8 & 2.86 & 0.00199 \\ 
  59 & ApNC & High & GO:0008217 & regulation of blood pressure &  10 &   2 & 0.06 & 0.00242 \\ 
  60 & ApNC & High & GO:0071216 & cellular response to biotic stimulus &  19 &   2 & 0.12 & 0.00424 \\ 
  61 & ApNC & High & GO:0002790 & peptide secretion &  19 &   2 & 0.12 & 0.00508 \\ 
  62 & ApNC & High & GO:0009914 & hormone transport &  23 &   2 & 0.14 & 0.00515 \\ 
  63 & ApNC & High & GO:0030099 & myeloid cell differentiation &  41 &   2 & 0.26 & 0.00623 \\ 
  64 & ApNC & High & GO:0002831 & regulation of response to biotic stimulu... &  24 &   2 & 0.15 & 0.00684 \\ 
  65 & ApNC & High & GO:0018193 & peptidyl-amino acid modification & 451 &   7 & 2.81 & 0.00687 \\ 
  66 & ApNC & High & GO:0035967 & cellular response to topologically incor... &  25 &   2 & 0.16 & 0.00703 \\ 
  67 & ApNC & High & GO:0007249 & I-kappaB kinase/NF-kappaB signaling &  25 &   2 & 0.16 & 0.00777 \\ 
  68 & ApNC & High & GO:0034976 & response to endoplasmic reticulum stress &  30 &   2 & 0.19 & 0.00874 \\ 
  69 & ApNC & High & GO:0001501 & skeletal system development &  60 &   2 & 0.37 & 0.00911 \\ 
  70 & ApNC & Low & GO:0015074 & DNA integration & 852 &  30 & 8.20 & 3.9e-06 \\ 
  71 & ApNC & Low & GO:0009164 & nucleoside catabolic process & 1733 &  24 & 16.68 & 0.00085 \\ 
  72 & ApNC & Low & GO:0006259 & DNA metabolic process & 2944 &  50 & 28.34 & 0.00105 \\ 
  73 & ApNC & Low & GO:1901658 & glycosyl compound catabolic process & 1734 &  24 & 16.69 & 0.00113 \\ 
  74 & ApNC & Low & GO:0046434 & organophosphate catabolic process & 1998 &  26 & 19.23 & 0.00155 \\ 
  75 & ApNC & Low & GO:0046483 & heterocycle metabolic process & 11355 & 141 & 109.30 & 0.00184 \\ 
  76 & ApNC & Low & GO:0008152 & metabolic process & 24338 & 254 & 234.28 & 0.00218 \\ 
  77 & ApNC & Low & GO:0007249 & I-kappaB kinase/NF-kappaB signaling &  25 &   3 & 0.24 & 0.00266 \\ 
  78 & ApNC & Low & GO:0043122 & regulation of I-kappaB kinase/NF-kappaB ... &  21 &   3 & 0.20 & 0.00268 \\ 
  79 & ApNC & Low & GO:0034641 & cellular nitrogen compound metabolic pro... & 11484 & 140 & 110.54 & 0.00325 \\ 
  80 & ApNC & Low & GO:0009059 & macromolecule biosynthetic process & 7008 &  93 & 67.46 & 0.00369 \\ 
  81 & ApNC & Low & GO:1901292 & nucleoside phosphate catabolic process & 1769 &  24 & 17.03 & 0.00385 \\ 
  82 & ApNC & Low & GO:0051049 & regulation of transport & 269 &   6 & 2.59 & 0.00461 \\ 
  83 & ApNC & Low & GO:0033002 & muscle cell proliferation &  12 &   2 & 0.12 & 0.00482 \\ 
  84 & ApNC & Low & GO:0060548 & negative regulation of cell death & 151 &   5 & 1.45 & 0.00532 \\ 
  85 & ApNC & Low & GO:0006760 & folic acid-containing compound metabolic... & 170 &   5 & 1.64 & 0.00580 \\ 
  86 & ApNC & Low & GO:0033013 & tetrapyrrole metabolic process & 198 &   7 & 1.91 & 0.00644 \\ 
  87 & ApNC & Low & GO:0034645 & cellular macromolecule biosynthetic proc... & 6965 &  93 & 67.04 & 0.00644 \\ 
  88 & ApNC & Low & GO:0006424 & glutamyl-tRNA aminoacylation &  31 &   3 & 0.30 & 0.00663 \\ 
  89 & ApNC & Low & GO:0006139 & nucleobase-containing compound metabolic... & 10066 & 126 & 96.90 & 0.00709 \\ 
  90 & ApNC & Low & GO:0072330 & monocarboxylic acid biosynthetic process & 550 &   9 & 5.29 & 0.00830 \\ 
  91 & ApNC & Low & GO:0042126 & nitrate metabolic process &  51 &   3 & 0.49 & 0.00916 \\ 
  92 & ApNC & Low & GO:0051604 & protein maturation & 2323 &  37 & 22.36 & 0.00917 \\ 
  93 & ApNC & Low & GO:0043069 & negative regulation of programmed cell d... & 142 &   5 & 1.37 & 0.00973 \\ 
  94 & ApNC & Low & GO:0003012 & muscle system process &  42 &   3 & 0.40 & 0.00973 \\ 
  95 & ApNC & Intermediate & GO:0044267 & cellular protein metabolic process & 5294 & 182 & 119.08 & 1.0e-08 \\ 
  96 & ApNC & Intermediate & GO:0019538 & protein metabolic process & 7596 & 230 & 170.86 & 1.5e-08 \\ 
  97 & ApNC & Intermediate & GO:0043170 & macromolecule metabolic process & 13944 & 350 & 313.64 & 3.1e-05 \\ 
  98 & ApNC & Intermediate & GO:0006996 & organelle organization & 2042 &  75 & 45.93 & 5.0e-05 \\ 
  99 & ApNC & Intermediate & GO:0071840 & cellular component organization or bioge... & 4430 & 137 & 99.64 & 5.4e-05 \\ 
  100 & ApNC & Intermediate & GO:0006644 & phospholipid metabolic process & 500 &  24 & 11.25 & 8.4e-05 \\ 
  101 & ApNC & Intermediate & GO:0016192 & vesicle-mediated transport & 744 &  31 & 16.73 & 9.5e-05 \\ 
  102 & ApNC & Intermediate & GO:1902589 & single-organism organelle organization & 1407 &  52 & 31.65 & 0.00020 \\ 
  103 & ApNC & Intermediate & GO:0007264 & small GTPase mediated signal transductio... & 659 &  28 & 14.82 & 0.00020 \\ 
  104 & ApNC & Intermediate & GO:0006184 & GTP catabolic process & 864 &  26 & 19.43 & 0.00032 \\ 
  105 & ApNC & Intermediate & GO:0070085 & glycosylation & 149 &   9 & 3.35 & 0.00038 \\ 
  106 & ApNC & Intermediate & GO:0044260 & cellular macromolecule metabolic process & 12209 & 307 & 274.62 & 0.00044 \\ 
  107 & ApNC & Intermediate & GO:0010033 & response to organic substance & 577 &  19 & 12.98 & 0.00045 \\ 
  108 & ApNC & Intermediate & GO:1901069 & guanosine-containing compound catabolic ... & 865 &  26 & 19.46 & 0.00055 \\ 
  109 & ApNC & Intermediate & GO:0061025 & membrane fusion &  44 &   4 & 0.99 & 0.00058 \\ 
  110 & ApNC & Intermediate & GO:0046039 & GTP metabolic process & 919 &  26 & 20.67 & 0.00063 \\ 
  111 & ApNC & Intermediate & GO:1900542 & regulation of purine nucleotide metaboli... & 514 &  18 & 11.56 & 0.00066 \\ 
  112 & ApNC & Intermediate & GO:0051649 & establishment of localization in cell & 1323 &  44 & 29.76 & 0.00075 \\ 
  113 & ApNC & Intermediate & GO:0022406 & membrane docking &  34 &   5 & 0.76 & 0.00079 \\ 
  114 & ApNC & Intermediate & GO:0022904 & respiratory electron transport chain & 368 &  14 & 8.28 & 0.00081 \\ 
  115 & ApNC & Intermediate & GO:0044801 & single-organism membrane fusion &  44 &   4 & 0.99 & 0.00101 \\ 
  116 & ApNC & Intermediate & GO:0043412 & macromolecule modification & 3268 & 107 & 73.51 & 0.00116 \\ 
  117 & ApNC & Intermediate & GO:0046907 & intracellular transport & 960 &  34 & 21.59 & 0.00122 \\ 
  118 & ApNC & Intermediate & GO:0033121 & regulation of purine nucleotide cataboli... & 500 &  18 & 11.25 & 0.00137 \\ 
  119 & ApNC & Intermediate & GO:1901068 & guanosine-containing compound metabolic ... & 953 &  27 & 21.44 & 0.00152 \\ 
  120 & ApNC & Intermediate & GO:0051641 & cellular localization & 1501 &  47 & 33.76 & 0.00178 \\ 
  121 & ApNC & Intermediate & GO:0030811 & regulation of nucleotide catabolic proce... & 500 &  18 & 11.25 & 0.00255 \\ 
  122 & ApNC & Intermediate & GO:0051345 & positive regulation of hydrolase activit... & 221 &  12 & 4.97 & 0.00264 \\ 
  123 & ApNC & Intermediate & GO:0006486 & protein glycosylation & 141 &   9 & 3.17 & 0.00419 \\ 
  124 & ApNC & Intermediate & GO:0043547 & positive regulation of GTPase activity & 182 &  12 & 4.09 & 0.00485 \\ 
  125 & ApNC & Intermediate & GO:0006412 & translation & 1886 &  66 & 42.42 & 0.00518 \\ 
  126 & ApNC & Intermediate & GO:0006446 & regulation of translational initiation & 344 &  15 & 7.74 & 0.00523 \\ 
  127 & ApNC & Intermediate & GO:0032970 & regulation of actin filament-based proce... &  88 &   7 & 1.98 & 0.00588 \\ 
  128 & ApNC & Intermediate & GO:0042439 & ethanolamine-containing compound metabol... &  36 &   4 & 0.81 & 0.00755 \\ 
  129 & ApNC & Intermediate & GO:0000725 & recombinational repair &  48 &   4 & 1.08 & 0.00799 \\ 
  130 & ApNC & Intermediate & GO:0009118 & regulation of nucleoside metabolic proce... & 500 &  18 & 11.25 & 0.00804 \\ 
  131 & ApNC & Intermediate & GO:0032956 & regulation of actin cytoskeleton organiz... &  87 &   7 & 1.96 & 0.00840 \\ 
  132 & ApNC & Intermediate & GO:0006140 & regulation of nucleotide metabolic proce... & 521 &  18 & 11.72 & 0.00853 \\ 
  133 & ApNC & Intermediate & GO:0050684 & regulation of mRNA processing &  59 &   5 & 1.33 & 0.00882 \\ 
  134 & ApNC & Intermediate & GO:0006656 & phosphatidylcholine biosynthetic process &  12 &   2 & 0.27 & 0.00904 \\ 
  135 & ApNC & Intermediate & GO:0031329 & regulation of cellular catabolic process & 572 &  21 & 12.87 & 0.00912 \\ 
  136 & ApNC & Intermediate & GO:0008654 & phospholipid biosynthetic process & 213 &  11 & 4.79 & 0.00921 \\ 
  137 & ApNC & Intermediate & GO:0048278 & vesicle docking &  33 &   5 & 0.74 & 0.00928 \\ 
  138 & ApNC & Bimodal & GO:0006259 & DNA metabolic process & 2944 &  28 & 13.93 & 0.00022 \\ 
  139 & ApNC & Bimodal & GO:0003151 & outflow tract morphogenesis &  10 &   2 & 0.05 & 0.00103 \\ 
  140 & ApNC & Bimodal & GO:0071804 & cellular potassium ion transport &  50 &   3 & 0.24 & 0.00168 \\ 
  141 & ApNC & Bimodal & GO:0071805 & potassium ion transmembrane transport &  50 &   3 & 0.24 & 0.00262 \\ 
  142 & ApNC & Bimodal & GO:0031023 & microtubule organizing center organizati... & 108 &   3 & 0.51 & 0.00431 \\ 
  143 & ApNC & Bimodal & GO:0006725 & cellular aromatic compound metabolic pro... & 11495 &  62 & 54.39 & 0.00521 \\ 
  144 & ApNC & Bimodal & GO:0072528 & pyrimidine-containing compound biosynthe... & 214 &   4 & 1.01 & 0.00649 \\ 
  145 & ApNC & Bimodal & GO:0034641 & cellular nitrogen compound metabolic pro... & 11484 &  61 & 54.34 & 0.00937 \\ 
   \hline
\end{tabular}
\end{table}




## Visualize responsive transcripts

Organize data


```r
# extract TPM data for thermally-responsive transcripts
resp.TPM.dt.sub <- TPM.dt.sub[names(responsive.lms)]
setkey(resp.TPM.dt.sub, Transcript)
str(resp.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	201432 obs. of  10 variables:
##  $ Transcript       : chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ Length           : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM              : num  0 0.3553 0 0.7353 0.0744 ...
##  $ RPKM             : num  0 0.599 0 1.1229 0.0928 ...
##  $ KPKM             : num  0 0.599 0 1.1229 0.0928 ...
##  $ EstimatedNumReads: num  0 157.6 0 366.1 61.2 ...
##  $ V7               : num  0 1.953 0 4.535 0.763 ...
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
## [1] 9156
```

```r
# scale transcripts so can compare
resp.TPM.dt.sub[, `:=`(TPM.scaled, scale(TPM)), by = Transcript]
```

```
##                          Transcript Length    TPM   RPKM   KPKM
##      1: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      2: 100008|*|comp137625_c0_seq2    208 0.3553 0.5990 0.5990
##      3: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      4: 100008|*|comp137625_c0_seq2    208 0.7353 1.1229 1.1229
##      5: 100008|*|comp137625_c0_seq2    208 0.0744 0.0928 0.0928
##     ---                                                        
## 201428:      9|*|comp147140_c0_seq1   9030 0.7167 1.2119 1.2119
## 201429:      9|*|comp147140_c0_seq1   9030 0.7239 1.0022 1.0022
## 201430:      9|*|comp147140_c0_seq1   9030 0.4328 0.7932 0.7932
## 201431:      9|*|comp147140_c0_seq1   9030 0.5655 0.7745 0.7745
## 201432:      9|*|comp147140_c0_seq1   9030 0.4626 0.8465 0.8465
##         EstimatedNumReads      V7 sample  val colony TPM.scaled
##      1:               0.0   0.000  A22-0  0.0    A22     -0.565
##      2:             157.6   1.953   Ar-0  0.0     Ar      1.128
##      3:               0.0   0.000  A22-3  3.5    A22     -0.565
##      4:             366.1   4.535   Ar-3  3.5     Ar      2.938
##      5:              61.2   0.763 A22-10 10.5    A22     -0.211
##     ---                                                        
## 201428:           17725.0 219.162  Ar-31 31.5     Ar     -0.274
## 201429:           23921.0 298.840 A22-35 35.0    A22     -0.251
## 201430:            8729.7 108.110  Ar-35 35.0     Ar     -1.184
## 201431:           23814.1 297.384 A22-38 38.5    A22     -0.759
## 201432:           12716.0 157.120  Ar-38 38.5     Ar     -1.089
```


Make plots for all genes expressed at *High* temps in GO category "GO:0006950: response to stress"

![plot of chunk plot_GOstress](figure/plot_GOstress1.png) ![plot of chunk plot_GOstress](figure/plot_GOstress2.png) 


Make plots for all genes expressed at *Low* temps in GO category "GO:0006950: response to stress"

![plot of chunk plot_GOstress_low](figure/plot_GOstress_low1.png) ![plot of chunk plot_GOstress_low](figure/plot_GOstress_low2.png) 



## Shiny interactive web-app

To assist visualization of specific transcripts, I made a interactive web-app using the [shiny](http://www.rstudio.com/shiny/) package. The scripts for this app are in the sub-directory `.\ApRxN-shinyapp`.

Export data for interactive shiny app. 






```
## Error: object 'w1' not found
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
## R version 3.1.0 (2014-04-10)
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
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] Rgraphviz_2.8.1      topGO_2.16.0         SparseM_1.03        
##  [4] GO.db_2.14.0         RSQLite_0.11.4       DBI_0.2-7           
##  [7] AnnotationDbi_1.26.0 GenomeInfoDb_1.0.2   Biobase_2.24.0      
## [10] BiocGenerics_0.10.0  graph_1.42.0         MASS_7.3-31         
## [13] plyr_1.8.1           RCurl_1.95-4.1       bitops_1.0-6        
## [16] data.table_1.9.2     stringr_0.6.2        pander_0.3.8        
## [19] knitcitations_0.5-0  bibtex_0.3-6         ggplot2_0.9.3.1     
## [22] R.utils_1.29.8       R.oo_1.18.0          R.methodsS3_1.6.1   
## [25] knitr_1.5           
## 
## loaded via a namespace (and not attached):
##  [1] codetools_0.2-8    colorspace_1.2-4   dichromat_2.0-0    digest_0.6.4      
##  [5] evaluate_0.5.3     formatR_0.10       gtable_0.1.2       httr_0.3          
##  [9] IRanges_1.22.3     labeling_0.2       lattice_0.20-29    munsell_0.4.2     
## [13] proto_0.3-10       RColorBrewer_1.0-5 Rcpp_0.11.1        reshape2_1.2.2    
## [17] scales_0.2.3       stats4_3.1.0       tools_3.1.0        XML_3.98-1.1      
## [21] xtable_1.7-3
```


## References


- Simon Anders, Davis J McCarthy, Yunshun Chen, Michal Okoniewski, Gordon K Smyth, Wolfgang Huber, Mark D Robinson,   (2013) Count-Based Differential Expression Analysis of Rna Sequencing Data Using R And Bioconductor.  *Nature Protocols*  **8**  1765-1786  [10.1038/nprot.2013.099](http://dx.doi.org/10.1038/nprot.2013.099)
- James H Bullard, Elizabeth Purdom, Kasper D Hansen, Sandrine Dudoit,   (2010) Evaluation of Statistical Methods For Normalization And Differential Expression in Mrna-Seq Experiments.  *Bmc Bioinformatics*  **11**  94-NA  [10.1186/1471-2105-11-94](http://dx.doi.org/10.1186/1471-2105-11-94)
- Manfred G Grabherr, Brian J Haas, Moran Yassour, Joshua Z Levin, Dawn A Thompson, Ido Amit, Xian Adiconis, Lin Fan, Raktima Raychowdhury, Qiandong Zeng, Zehua Chen, Evan Mauceli, Nir Hacohen, Andreas Gnirke, Nicholas Rhind, Federica di Palma, Bruce W Birren, Chad Nusbaum, Kerstin Lindblad-Toh, Nir Friedman, Aviv Regev,   (2011) Full-Length Transcriptome Assembly From Rna-Seq Data Without A Reference Genome.  *Nature Biotechnology*  **29**  644-652  [10.1038/nbt.1883](http://dx.doi.org/10.1038/nbt.1883)
- X. Huang,   (1999) Cap3: A Dna Sequence Assembly Program.  *Genome Research*  **9**  868-877  [10.1101/gr.9.9.868](http://dx.doi.org/10.1101/gr.9.9.868)
- B. Li, V. Ruotti, R. M. Stewart, J. A. Thomson, C. N. Dewey,   (2009) Rna-Seq Gene Expression Estimation With Read Mapping Uncertainty.  *Bioinformatics*  **26**  493-500  [10.1093/bioinformatics/btp692](http://dx.doi.org/10.1093/bioinformatics/btp692)
- M. Lohse, A. M. Bolger, A. Nagel, A. R. Fernie, J. E. Lunn, M. Stitt, B. Usadel,   (2012) Robina: A User-Friendly, Integrated Software Solution For Rna-Seq-Based Transcriptomics.  *Nucleic Acids Research*  **40**  W622-W627  [10.1093/nar/gks540](http://dx.doi.org/10.1093/nar/gks540)
- David Lubertazzi,   (2012) The Biology And Natural History of Aphaenogaster Rudis.  *Psyche: A Journal of Entomology*  **2012**  1-11  [10.1155/2012/752815](http://dx.doi.org/10.1155/2012/752815)
- Courtney J. Murren, Heidi J. Maclean, Sarah E. Diamond, Ulrich K. Steiner, Mary A. Heskel, Corey A. Handelsman, Cameron K. Ghalambor, Josh R. Auld, Hilary S. Callahan, David W. Pfennig, Rick A. Relyea, Carl D. Schlichting, Joel Kingsolver,   (2014) Evolutionary Change in Continuous Reaction Norms.  *The American Naturalist*  **183**  453-467  [10.1086/675302](http://dx.doi.org/10.1086/675302)
- Robert Schmieder, Robert Edwards, Francisco Rodriguez-Valera,   (2011) Fast Identification And Removal of Sequence Contamination From Genomic And Metagenomic Datasets.  *Plos One*  **6**  e17288-NA  [10.1371/journal.pone.0017288](http://dx.doi.org/10.1371/journal.pone.0017288)
- unknown unknown,   (unknown) Unknown.  *Unknown*
- Ya Yang, Stephen A Smith,   (2013) Optimizing de Novo Assembly of Short-Read Rna-Seq Data For Phylogenomics.  *Bmc Genomics*  **14**  328-NA  [10.1186/1471-2164-14-328](http://dx.doi.org/10.1186/1471-2164-14-328)

