Thermal reactionome of the common ant species *Aphaenogaster picea* and *A. carolinensis*
========================================================================================
   
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**May 14, 2014**

**Technical Report No. 3**

**Department of Biology**

**University of Vermont**





## Summary ##
  
In this technical report, which accompanies the manuscript **Thermal reactionome of a common ant species** (Stanton-Geddes et al., in press), we:

1. describe the *de novo* assembly of the transcriptome for two ant species within the *Aphaenogaster rudis-picea-fulva* species complex (<a href="http://dx.doi.org/10.1155/2012/752815">Lubertazzi, 2012</a>)
2. identify thermally-responsive genes
3. evaluate differences in the expression patterns between the two species
3. perform gene set enrichment analysis of thermally-responsive genes for the two species

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

Two ant colonies were used for the transcriptome sequencing. The first, designated *A22*, was collected at Molly Bog, Vermont in August 2012 by Nick Gotelli and Andrew Nguyen. The second colony, designated *Ar*, was collected by Lauren Nichols in Raleigh, North Carolina. These colonies were maintained in the lab for 6 months prior to sample collection. Bernice Bacon DeMarco (Michigan State University) identified colony *A22* as *A. picea* and *Ar* as *A. carolinensis*.

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
                                                 
    sailfish quant -i sailfish-index-Trinity-cap3-uclust -o sailfish-expression-Trinity-cap3-uclust/A22-0 -l "T=SE:S=U" -r A22-0_ATCACG.unpaired.left.fastq A22-0_ATCACG.unpaired.right.fastq A22-0_ATCACG.paired.left.fastq A22-0_ATCACG.paired.right.fastq -p 4
	
While it is possible to separately specify the paired-end and orphaned single-end reads in Sailfish v0.6.3, the results are exactly the same as if they are all entered as SE.	

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

# calculate overall P value and R^2 for each transcript
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
# perform model selection for responsive transcripts need to use `try` to avoid
# stopping on error for AIC at Infinity
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
write.table(file = paste(resultsdir, "Ap-responsive-transcripts", Sys.Date(), ".txt", 
    sep = ""), Ap.response.type, row.names = FALSE, col.names = TRUE, quote = FALSE, 
    sep = "\t")
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
w1 <- wilcox.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE)
w1
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  A22.high.transcripts$A22.opt and A22.high.transcripts$Ar.opt
## V = 210415, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.294 -0.171
## sample estimates:
## (pseudo)median 
##          -0.23
```


Consistent with expectation, there is greater expression at 19.5C for *Ar* than *A22* transcripts for the set of transcripts that are transcripts that are up-regulated at high temperatures in *A22*. Note that A22 had the larger library size so if this was due to TPM not correctly accounting for differences in reads between samples, we would expect to see a positive instead of negative value here.

Next I test the converse, that transcripts that are up-regulated at low temperatures in *Ar* are more highly-expressed at 19.25C temperatures in *A22*. 


```r
# list of transcripts that are 'high' expressed in Ar
Ar.low.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Low"), 
    ]

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

![plot of chunk Ar_low_wilcoxon](figure/Ar_low_wilcoxon.png) 

```r

# Wilcoxon signed rank-test
w2 <- wilcox.test(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE)
w2
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

![plot of chunk plot_thermal_breadth](figure/plot_thermal_breadth.png) 



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
    "AcNC"), each = length(Ap.response.type$Transcript)), max.val = c(Ap.response.type$A22.max, 
    Ap.response.type$Ar.max))

maxexplot <- ggplot(Ap.df, aes(x = max.val, fill = colony)) + geom_density(alpha = 0.2, 
    position = "identity") + # scale_fill_manual(name = 'Colony', values = c('white', 'black')) +
scale_y_continuous(name = "Density") + scale_x_continuous(name = "Temperature of maximum expression")
suppressWarnings(print(maxexplot))
```

![plot of chunk max_exp_PDF](figure/max_exp_PDF.png) 


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

Thermally-responsive genes could also differ in the temperatures at which they have increased or decreased expression in response to temperature changes. To examine this, I determine the temperature at which each responsive gene has the greated increase or decrease in expression.


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


Predict expression for responsive transcripts.


```r
# apply predFunc to all responsive transcripts
resp.TPM.dt.sub.pred <- ddply(resp.TPM.dt.sub, .(Transcript), .inform = "TRUE", predFunc)

# setkey to Transcript and colony
resp.TPM.dt.sub.pred <- data.table(resp.TPM.dt.sub.pred)
setkey(resp.TPM.dt.sub.pred, Transcript, colony)
```


For next analyses, extract list of gene names by response type.


```r
A22.high <- Ap.response.type[which(Ap.response.type$A22.type == "High" & Ap.response.type$Ar.type != 
    "High"), "Transcript"]
Ar.high <- Ap.response.type[which(Ap.response.type$A22.type != "High" & Ap.response.type$Ar.type == 
    "High"), "Transcript"]

A22.low <- Ap.response.type[which(Ap.response.type$A22.type == "Low" & Ap.response.type$Ar.type != 
    "Low"), "Transcript"]
Ar.low <- Ap.response.type[which(Ap.response.type$A22.type != "Low" & Ap.response.type$Ar.type == 
    "Low"), "Transcript"]

A22.bim <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal" & Ap.response.type$Ar.type != 
    "Bimodal"), "Transcript"]
Ar.bim <- Ap.response.type[which(Ap.response.type$A22.type != "Bimodal" & Ap.response.type$Ar.type == 
    "Bimodal"), "Transcript"]

A22.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & Ap.response.type$Ar.type != 
    "Intermediate"), "Transcript"]
Ar.int <- Ap.response.type[which(Ap.response.type$A22.type != "Intermediate" & Ap.response.type$Ar.type == 
    "Intermediate"), "Transcript"]
```



Calculate `T_on` for *High* genes in each colony.


```r
# transcripts expressed at *High* and *Bimodal* temperatures in A22
A22.high.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(A22.high, A22.bim), "A22")]
setkey(A22.high.TPM.dt.sub, Transcript)
str(A22.high.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	24200 obs. of  5 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0 0 0.0744 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.01 1.01 1.01 1.01 1.01 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r

# make data.frame for results
A22.high.T_on <- data.frame(Transcript = unique(A22.high.TPM.dt.sub$Transcript), 
    colony = rep(A22.high.TPM.dt.sub$colony, length = length(unique(A22.high.TPM.dt.sub$Transcript))), 
    T_on = NA, pT_on = NA)

# loop across transcripts, calculating T_on

for (i in unique(A22.high.TPM.dt.sub$Transcript)) {
    subdf <- A22.high.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val > 14), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM)))) + 1, val]
    pT_on <- subdf[median(which(diff(subdf$pTPM) == max(diff(subdf$pTPM)))) + 1, 
        val]
    A22.high.T_on[which(A22.high.T_on$Transcript == i), "T_on"] <- T_on
    A22.high.T_on[which(A22.high.T_on$Transcript == i), "pT_on"] <- pT_on
}

# repeat for Ar
Ar.high.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(Ar.high, Ar.bim), "Ar")]
str(Ar.high.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	14432 obs. of  5 variables:
##  $ Transcript: chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 2 2 2 2 2 2 2 2 2 2 ...
##  $ TPM       : num  0 0 0 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  0.997 0.999 1.001 1.002 1.004 ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r

Ar.high.T_on <- data.frame(Transcript = unique(Ar.high.TPM.dt.sub$Transcript), colony = rep(Ar.high.TPM.dt.sub$colony, 
    length = length(unique(Ar.high.TPM.dt.sub$Transcript))), T_on = NA, pT_on = NA)

for (i in unique(Ar.high.TPM.dt.sub$Transcript)) {
    subdf <- Ar.high.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val > 14), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM)))) + 1, val]
    pT_on <- subdf[median(which(diff(subdf$pTPM) == max(diff(subdf$pTPM)))) + 1, 
        val]
    Ar.high.T_on[which(Ar.high.T_on$Transcript == i), "T_on"] <- T_on
    Ar.high.T_on[which(Ar.high.T_on$Transcript == i), "pT_on"] <- pT_on
}
```

```
## Error: When i is a data.table (or character vector), x must be keyed (i.e.
## sorted, and, marked as sorted) so data.table knows which columns to join to and
## take advantage of x being sorted. Call setkey(x,...) first, see ?setkey.
```

```r

cor.test(Ar.high.T_on$T_on, Ar.high.T_on$pT_on)
```

```
## Error: 'x' must be a numeric vector
```

```r
cor.test(A22.high.T_on$T_on, A22.high.T_on$pT_on)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  A22.high.T_on$T_on and A22.high.T_on$pT_on
## t = 10.6, df = 2198, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.18 0.26
## sample estimates:
##  cor 
## 0.22
```


Determine if `T_on` is greater in *A22* or *Ar* for *High* genes.


```r
t.test(Ar.high.T_on$T_on, A22.high.T_on$T_on)
```

```
## Error: not enough 'x' observations
```

```r
t.test(Ar.high.T_on$pT_on, A22.high.T_on$pT_on)
```

```
## Error: not enough 'x' observations
```

Genes with increased expression at *High* temperatures are turned on at higher temperatures in *ApVT* than *AcNC*.

Plot `T_on` for *High* genes.


```
## Error: object 'A22.high.T_on' not found
```

```
## Error: object 'Ap.high.T_on' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'print': Error: object 'T_on_plot_high' not found
```

```
## Error: object 'Ap.high.T_on' not found
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'print': Error: object 'pT_on_plot_high' not found
```


Repeat analysis for *Low* genes.



























































