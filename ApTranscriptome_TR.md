Thermal reactionome of the common ant species *Aphaenogaster*
================================================================
   
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**April 29, 2014**

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
	
This resulted in removing 5,675 contigs as contaminants, leaving 99,861 "clean" contigs. We spot-checked the contaminants by BLAST and confirmed that they matched bacteria, human or viral sources by greater than 95%. 


Running Trinity and subsequent programs is time and memory-intensive, so the final filtered clean assembly can be downloaded from [...] 

~~~
# download filtered Trinity assembly, uncompress and move
wget http://johnstantongeddes.org/assets/files/Aphaenogaster_transcriptome.tar
tar -xvf Aphaenogaster_transcriptome.tar
# check md5sum
md5sum Trinity_cap3_uclust_clean.fa
# c604159ee944dddda7b4e2a556bf2f53
mkdir -p results/
mkdir -p results/trinity-full/
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
# a2 <- getURL('https://googledrive.com/host/0B75IymziRJ_9Tlg1U1Vxbjk1bzg')
# # GoogleDrive link

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

I quantified gene expression using [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html). To run this program, first make sure that PATHs to the software libraries are set up correctly: 
                                                 
    export LD_LIBRARY_PATH=/opt/software/Sailfish-0.6.2-Linux_x86-64/lib:$LD_LIBRARY_PATH
    export PATH=/opt/software/Sailfish-0.6.2-Linux_x86-64/bin:$PATH

Then build the index of the assembly:

    sailfish index -t Trinity_cap3_uclust_clean.fasta -o sailfish-index -k 20 -p 4

Once this is done, quantify expression for the Trimmomatic filtered reads from each colony-treatment sample separately. Note that for each sample, there are three four filtered read files:

- paired.left.fastq
- paired.right.fastq
- unpaired.left.fastq
- unpaired.right.fastq
                                                 
Make a directory for the expression values

    mkdir sailfish-expression                                                 

Then, for each sample, run the following command:
                                                 
    sailfish -i sailfish-index -o sailfish-expression/A22-0 --reads A22-0_ATCACG.paired.left.fastq A22-0_ATCACG.paired.right.fastq A22-0_ATCACG.unpaired.left.fastq A22-0_ATCACG.unpaired.right.fastq -p 4

Or, using a loop in R:                                                 
                                                 

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
## Start expression quantification for sample A22-0: 2014-04-29 10:23:06
## Done with expression quantification for sample A22-0: 2014-04-29 10:26:29
## 
## Start expression quantification for sample A22-10: 2014-04-29 10:26:29
## Done with expression quantification for sample A22-10: 2014-04-29 10:29:52
## 
## Start expression quantification for sample A22-14: 2014-04-29 10:29:52
## Done with expression quantification for sample A22-14: 2014-04-29 10:33:14
## 
## Start expression quantification for sample A22-17: 2014-04-29 10:33:14
## Done with expression quantification for sample A22-17: 2014-04-29 10:36:36
## 
## Start expression quantification for sample A22-21: 2014-04-29 10:36:36
## Done with expression quantification for sample A22-21: 2014-04-29 10:39:58
## 
## Start expression quantification for sample A22-24: 2014-04-29 10:39:58
## Done with expression quantification for sample A22-24: 2014-04-29 10:43:20
## 
## Start expression quantification for sample A22-28: 2014-04-29 10:43:20
## Done with expression quantification for sample A22-28: 2014-04-29 10:46:42
## 
## Start expression quantification for sample A22-31: 2014-04-29 10:46:42
## Done with expression quantification for sample A22-31: 2014-04-29 10:50:04
## 
## Start expression quantification for sample A22-35: 2014-04-29 10:50:04
## Done with expression quantification for sample A22-35: 2014-04-29 10:53:26
## 
## Start expression quantification for sample A22-38: 2014-04-29 10:53:26
## Done with expression quantification for sample A22-38: 2014-04-29 10:56:48
## 
## Start expression quantification for sample A22-3: 2014-04-29 10:56:48
## Done with expression quantification for sample A22-3: 2014-04-29 11:00:11
## 
## Start expression quantification for sample A22-7: 2014-04-29 11:00:11
## Done with expression quantification for sample A22-7: 2014-04-29 11:03:33
## 
## Start expression quantification for sample Ar-0: 2014-04-29 11:03:33
## Done with expression quantification for sample Ar-0: 2014-04-29 11:06:55
## 
## Start expression quantification for sample Ar-10: 2014-04-29 11:06:55
## Done with expression quantification for sample Ar-10: 2014-04-29 11:10:17
## 
## Start expression quantification for sample Ar-14: 2014-04-29 11:10:17
## Done with expression quantification for sample Ar-14: 2014-04-29 11:13:39
## 
## Start expression quantification for sample Ar-17: 2014-04-29 11:13:39
## Done with expression quantification for sample Ar-17: 2014-04-29 11:17:01
## 
## Start expression quantification for sample Ar-21: 2014-04-29 11:17:01
## Done with expression quantification for sample Ar-21: 2014-04-29 11:20:24
## 
## Start expression quantification for sample Ar-24: 2014-04-29 11:20:24
## Done with expression quantification for sample Ar-24: 2014-04-29 11:23:46
## 
## Start expression quantification for sample Ar-28: 2014-04-29 11:23:46
## Done with expression quantification for sample Ar-28: 2014-04-29 11:27:08
## 
## Start expression quantification for sample Ar-31: 2014-04-29 11:27:08
## Done with expression quantification for sample Ar-31: 2014-04-29 11:30:30
## 
## Start expression quantification for sample Ar-35: 2014-04-29 11:30:30
## Done with expression quantification for sample Ar-35: 2014-04-29 11:33:53
## 
## Start expression quantification for sample Ar-38: 2014-04-29 11:33:53
## Done with expression quantification for sample Ar-38: 2014-04-29 11:37:15
## 
## Start expression quantification for sample Ar-3: 2014-04-29 11:37:15
## Done with expression quantification for sample Ar-3: 2014-04-29 11:40:37
## 
## Start expression quantification for sample Ar-7: 2014-04-29 11:40:37
## Done with expression quantification for sample Ar-7: 2014-04-29 11:43:59
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
|   0    |  0.87  |
|  3.5   |  0.89  |
|   7    |  0.66  |
|  10.5  |  0.87  |
|   14   |  0.88  |
|  17.5  |  0.87  |
|   21   |  0.79  |
|  24.5  |  0.88  |
|   28   |  0.88  |
|  31.5  |  0.9   |
|   35   |  0.9   |
|  38.5  |  0.91  |

Table: correlations between colonies at each temperature treatment


## Identification of thermally-responsive genes

To identify transcripts (roughly equivalent to genes) that show thermal responsiveness, I fit the following linear model to each transcript:

$$ log(TPM + 1) = \beta_0 + \beta_1(colony) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(colony * temp) + \beta_5(colony * temp^2) + \epsilon $$

where TPM is transcripts per million. 



**Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the analysis.** 


```r
A22.TPM[, `:=`(colony, "A22")]
```

```
##                      Transcript Length  TPM RPKM KPKM EstimatedNumReads
##       1:         0|*|Contig6267   9990 4.09 7.11 7.11              23.7
##       2:         0|*|Contig6267   9990 2.84 4.97 4.97              17.3
##       3:         0|*|Contig6267   9990 5.41 9.39 9.39              28.0
##       4:         0|*|Contig6267   9990 3.81 6.64 6.64              24.8
##       5:         0|*|Contig6267   9990 3.07 5.38 5.38              21.2
##      ---                                                               
## 1198328: 9|*|comp147140_c0_seq1   9030 3.95 6.95 6.95              23.5
## 1198329: 9|*|comp147140_c0_seq1   9030 3.70 6.50 6.50              17.3
## 1198330: 9|*|comp147140_c0_seq1   9030 4.45 7.78 7.78              20.4
## 1198331: 9|*|comp147140_c0_seq1   9030 4.27 7.47 7.47              22.0
## 1198332: 9|*|comp147140_c0_seq1   9030 4.20 7.38 7.38              27.5
##          sample  val colony
##       1:  A22-0  0.0    A22
##       2:  A22-3  3.5    A22
##       3:  A22-7  7.0    A22
##       4: A22-10 10.5    A22
##       5: A22-14 14.0    A22
##      ---                   
## 1198328: A22-24 24.5    A22
## 1198329: A22-28 28.0    A22
## 1198330: A22-31 31.5    A22
## 1198331: A22-35 35.0    A22
## 1198332: A22-38 38.5    A22
```

```r
Ar.TPM[, `:=`(colony, "Ar")]
```

```
##                      Transcript Length  TPM  RPKM  KPKM EstimatedNumReads
##       1:         0|*|Contig6267   9990 4.01  6.91  6.91             12.04
##       2:         0|*|Contig6267   9990 4.72  8.15  8.15             15.93
##       3:         0|*|Contig6267   9990 1.56  2.79  2.79              9.43
##       4:         0|*|Contig6267   9990 3.98  6.89  6.89             14.07
##       5:         0|*|Contig6267   9990 4.84  8.38  8.38             20.46
##      ---                                                                 
## 1198328: 9|*|comp147140_c0_seq1   9030 5.64  9.75  9.75             18.62
## 1198329: 9|*|comp147140_c0_seq1   9030 5.57  9.61  9.61             17.95
## 1198330: 9|*|comp147140_c0_seq1   9030 6.16 10.64 10.64             18.90
## 1198331: 9|*|comp147140_c0_seq1   9030 5.04  8.68  8.68             13.01
## 1198332: 9|*|comp147140_c0_seq1   9030 5.33  9.19  9.19             18.00
##          sample  val colony
##       1:   Ar-0  0.0     Ar
##       2:   Ar-3  3.5     Ar
##       3:   Ar-7  7.0     Ar
##       4:  Ar-10 10.5     Ar
##       5:  Ar-14 14.0     Ar
##      ---                   
## 1198328:  Ar-24 24.5     Ar
## 1198329:  Ar-28 28.0     Ar
## 1198330:  Ar-31 31.5     Ar
## 1198331:  Ar-35 35.0     Ar
## 1198332:  Ar-38 38.5     Ar
```

```r
TPM.dt <- rbind(A22.TPM, Ar.TPM)
TPM.dt$colony <- as.factor(TPM.dt$colony)
str(TPM.dt)
```

```
## Classes 'data.table' and 'data.frame':	2396664 obs. of  9 variables:
##  $ Transcript       : chr  "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" ...
##  $ Length           : int  9990 9990 9990 9990 9990 9990 9990 9990 9990 9990 ...
##  $ TPM              : num  4.09 2.84 5.41 3.81 3.07 ...
##  $ RPKM             : num  7.11 4.97 9.39 6.64 5.38 ...
##  $ KPKM             : num  7.11 4.97 9.39 6.64 5.38 ...
##  $ EstimatedNumReads: num  23.7 17.3 28 24.8 21.2 ...
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

RxNpval <- ddply(TPM.dt.sub, .(Transcript), .progress = "text", .inform = "TRUE", 
    modpFunc)
```

```
##   |                                                                         |                                                                 |   0%  |                                                                         |                                                                 |   1%  |                                                                         |=                                                                |   1%  |                                                                         |=                                                                |   2%  |                                                                         |==                                                               |   2%  |                                                                         |==                                                               |   3%  |                                                                         |==                                                               |   4%  |                                                                         |===                                                              |   4%  |                                                                         |===                                                              |   5%  |                                                                         |====                                                             |   5%  |                                                                         |====                                                             |   6%  |                                                                         |====                                                             |   7%  |                                                                         |=====                                                            |   7%  |                                                                         |=====                                                            |   8%  |                                                                         |======                                                           |   8%  |                                                                         |======                                                           |   9%  |                                                                         |======                                                           |  10%  |                                                                         |=======                                                          |  10%  |                                                                         |=======                                                          |  11%  |                                                                         |=======                                                          |  12%  |                                                                         |========                                                         |  12%  |                                                                         |========                                                         |  13%  |                                                                         |=========                                                        |  13%  |                                                                         |=========                                                        |  14%  |                                                                         |=========                                                        |  15%  |                                                                         |==========                                                       |  15%  |                                                                         |==========                                                       |  16%  |                                                                         |===========                                                      |  16%  |                                                                         |===========                                                      |  17%  |                                                                         |===========                                                      |  18%  |                                                                         |============                                                     |  18%  |                                                                         |============                                                     |  19%  |                                                                         |=============                                                    |  19%  |                                                                         |=============                                                    |  20%  |                                                                         |=============                                                    |  21%  |                                                                         |==============                                                   |  21%  |                                                                         |==============                                                   |  22%  |                                                                         |===============                                                  |  22%  |                                                                         |===============                                                  |  23%  |                                                                         |===============                                                  |  24%  |                                                                         |================                                                 |  24%  |                                                                         |================                                                 |  25%  |                                                                         |=================                                                |  25%  |                                                                         |=================                                                |  26%  |                                                                         |=================                                                |  27%  |                                                                         |==================                                               |  27%  |                                                                         |==================                                               |  28%  |                                                                         |===================                                              |  28%  |                                                                         |===================                                              |  29%  |                                                                         |===================                                              |  30%  |                                                                         |====================                                             |  30%  |                                                                         |====================                                             |  31%  |                                                                         |====================                                             |  32%  |                                                                         |=====================                                            |  32%  |                                                                         |=====================                                            |  33%  |                                                                         |======================                                           |  33%  |                                                                         |======================                                           |  34%  |                                                                         |======================                                           |  35%  |                                                                         |=======================                                          |  35%  |                                                                         |=======================                                          |  36%  |                                                                         |========================                                         |  36%  |                                                                         |========================                                         |  37%  |                                                                         |========================                                         |  38%  |                                                                         |=========================                                        |  38%  |                                                                         |=========================                                        |  39%  |                                                                         |==========================                                       |  39%  |                                                                         |==========================                                       |  40%  |                                                                         |==========================                                       |  41%  |                                                                         |===========================                                      |  41%  |                                                                         |===========================                                      |  42%  |                                                                         |============================                                     |  42%  |                                                                         |============================                                     |  43%  |                                                                         |============================                                     |  44%  |                                                                         |=============================                                    |  44%  |                                                                         |=============================                                    |  45%  |                                                                         |==============================                                   |  45%  |                                                                         |==============================                                   |  46%  |                                                                         |==============================                                   |  47%  |                                                                         |===============================                                  |  47%  |                                                                         |===============================                                  |  48%  |                                                                         |================================                                 |  48%  |                                                                         |================================                                 |  49%  |                                                                         |================================                                 |  50%  |                                                                         |=================================                                |  50%  |                                                                         |=================================                                |  51%  |                                                                         |=================================                                |  52%  |                                                                         |==================================                               |  52%  |                                                                         |==================================                               |  53%  |                                                                         |===================================                              |  53%  |                                                                         |===================================                              |  54%  |                                                                         |===================================                              |  55%  |                                                                         |====================================                             |  55%  |                                                                         |====================================                             |  56%  |                                                                         |=====================================                            |  56%  |                                                                         |=====================================                            |  57%  |                                                                         |=====================================                            |  58%  |                                                                         |======================================                           |  58%  |                                                                         |======================================                           |  59%  |                                                                         |=======================================                          |  59%  |                                                                         |=======================================                          |  60%  |                                                                         |=======================================                          |  61%  |                                                                         |========================================                         |  61%  |                                                                         |========================================                         |  62%  |                                                                         |=========================================                        |  62%  |                                                                         |=========================================                        |  63%  |                                                                         |=========================================                        |  64%  |                                                                         |==========================================                       |  64%  |                                                                         |==========================================                       |  65%  |                                                                         |===========================================                      |  65%  |                                                                         |===========================================                      |  66%  |                                                                         |===========================================                      |  67%  |                                                                         |============================================                     |  67%  |                                                                         |============================================                     |  68%  |                                                                         |=============================================                    |  68%  |                                                                         |=============================================                    |  69%  |                                                                         |=============================================                    |  70%  |                                                                         |==============================================                   |  70%  |                                                                         |==============================================                   |  71%  |                                                                         |==============================================                   |  72%  |                                                                         |===============================================                  |  72%  |                                                                         |===============================================                  |  73%  |                                                                         |================================================                 |  73%  |                                                                         |================================================                 |  74%  |                                                                         |================================================                 |  75%  |                                                                         |=================================================                |  75%  |                                                                         |=================================================                |  76%  |                                                                         |==================================================               |  76%  |                                                                         |==================================================               |  77%  |                                                                         |==================================================               |  78%  |                                                                         |===================================================              |  78%  |                                                                         |===================================================              |  79%  |                                                                         |====================================================             |  79%  |                                                                         |====================================================             |  80%  |                                                                         |====================================================             |  81%  |                                                                         |=====================================================            |  81%  |                                                                         |=====================================================            |  82%  |                                                                         |======================================================           |  82%  |                                                                         |======================================================           |  83%  |                                                                         |======================================================           |  84%  |                                                                         |=======================================================          |  84%  |                                                                         |=======================================================          |  85%  |                                                                         |========================================================         |  85%  |                                                                         |========================================================         |  86%  |                                                                         |========================================================         |  87%  |                                                                         |=========================================================        |  87%  |                                                                         |=========================================================        |  88%  |                                                                         |==========================================================       |  88%  |                                                                         |==========================================================       |  89%  |                                                                         |==========================================================       |  90%  |                                                                         |===========================================================      |  90%  |                                                                         |===========================================================      |  91%  |                                                                         |===========================================================      |  92%  |                                                                         |============================================================     |  92%  |                                                                         |============================================================     |  93%  |                                                                         |=============================================================    |  93%  |                                                                         |=============================================================    |  94%  |                                                                         |=============================================================    |  95%  |                                                                         |==============================================================   |  95%  |                                                                         |==============================================================   |  96%  |                                                                         |===============================================================  |  96%  |                                                                         |===============================================================  |  97%  |                                                                         |===============================================================  |  98%  |                                                                         |================================================================ |  98%  |                                                                         |================================================================ |  99%  |                                                                         |=================================================================|  99%  |                                                                         |=================================================================| 100%
```


Of the 99861 transcripts, 44703 have models with P < 0.05.

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

```
## Error: typeof x.val (double) != typeof i.val (character)
```


At the 5% FDR significance threshold, there are 34366 transcripts with an overall significant model.


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
RxNlmAIC <- dlply(sig.TPM.dt.sub, .(Transcript), .progress = "text", lmFunc)
```

```
## Error: object 'sig.TPM.dt.sub' not found
```



## Thermally-responsive transcripts

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
```

```
## Error: object 'RxNlmAIC' not found
```

```r
other.lms <- RxNlmAIC[setdiff(names(RxNlmAIC), names(interaction.lms))]
```

```
## Error: object 'RxNlmAIC' not found
```

```r
temperature.lms <- other.lms[which(Map(grepFunc, other.lms, term = "val") == 
    TRUE)]
```

```
## Error: object 'other.lms' not found
```

```r
colony.lms <- other.lms[setdiff(names(other.lms), names(temperature.lms))]
```

```
## Error: object 'other.lms' not found
```

```r
responsive.lms <- c(temperature.lms, interaction.lms)
```

```
## Error: object 'temperature.lms' not found
```

```r
rm(other.lms)
```

```
## Warning: object 'other.lms' not found
```



```
## Error: object 'RxNlmAIC' not found
```

```
## Error: object 'siglist' not found
```

```
## Error: object 'sigtable' not found
```



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
interaction.response.type <- ldply(interaction.lms, .progress = "text", .inform = TRUE, 
    RxNtype)
```

```
## Error: object 'interaction.lms' not found
```

```r
stopifnot(nrow(interaction.response.type) == length(interaction.lms))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'interaction.response.type' not found
```

```r

# calculate response types for transcripts without interactions
temperature.response.type <- ldply(temperature.lms, .progress = "text", .inform = TRUE, 
    RxNtype)
```

```
## Error: object 'temperature.lms' not found
```

```r
stopifnot(nrow(temperature.response.type) == length(temperature.lms))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'temperature.response.type' not found
```

```r

# merge results
Ap.response.type <- rbind(interaction.response.type, temperature.response.type)
```

```
## Error: object 'interaction.response.type' not found
```

```r
colnames(Ap.response.type)[which(colnames(Ap.response.type) == ".id")] <- "Transcript"
```

```
## Error: object 'Ap.response.type' not found
```

```r
str(Ap.response.type)
```

```
## Error: object 'Ap.response.type' not found
```

```r

# save results to file
write.table(file = paste(resultsdir, "Ap_responsive_transcripts.txt", sep = ""), 
    Ap.response.type, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
```

```
## Error: object 'Ap.response.type' not found
```


Examine transcripts by temperature response. Compare between colonies.


```
## Error: object 'Ap.response.type' not found
```

```
## Error: object 'Ap.response.type' not found
```

```
## Error: object 'Ap.response.type' not found
```

```
## Error: object 'type.table' not found
```

```
## Error: object 'tt2' not found
```

pdf 
  2 

```
## Error: object 'tt2' not found
```

```
## Error: object 'Ap.type.table' not found
```



Table 4 shows the number of transcripts that fall into each expression type for each each colony. The totals for each colony include the 

```

Error in nrow(temperature.lms) : 
  error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'temperature.lms' not found


```

 transcripts that have consistent temperature responses between the two colonies. 



## Colony-level comparison ##

In this section, I compare the thermal reactionome between the *Ar* and *A22* colonies. 

### Plasticity versus constitutive expression

Selection may have acted in response to thermal stress such that some genes are constitutively activated in one colony while plastically-expressed in the other. I evaluate this hypothesis by comparing expression levels at the at optimum (19.5C) temperature between the two colonies for genes in that are either in the 'High' or 'Low' expression group in the other colony. Specifically, I expect that genes upregulated at high temperatures in *A22* are more highly expressed at 19.5C in *Ar*? Conversely, I expect genes that are upregulated at low temps in *Ar* to have greater constitutive expression in *A22*? 


```r
# list of transcripts that are 'high' expressed in A22
A22.high.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == 
    "High"), ]

# Compare expression at optimum temp (A22.opt) between colonies using t-test
t.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)
boxplot(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)

# T test on log-transformed values to control for outliers
t.test(log(A22.high.transcripts$A22.opt + 1), log(A22.high.transcripts$Ar.opt + 
    1))
boxplot(log(A22.high.transcripts$A22.opt + 1), log(A22.high.transcripts$Ar.opt + 
    1))
```


The `t.test` fails to account for the many orders of magnitude difference in expression among transcripts, e.g. non-equal variances. This problem is the key issue in the analysis of differential expression 

```

Error in base::parse(text = code, srcfile = NULL) : 
  2:0: unexpected end of input
1: citep(c("10.1186/1471-2105-11-94", "10.1038/nprot.2013.099")
   ^

```

. As my goal is simply to determine if expression is typically greater at 19.25C in *Ar* than *A22* for genes that are up-regulated at high temperatures in *A22*, I use a non-parametric Wilcoxon signed rank-test


```r
(w1 <- wilcox.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt, 
    alternative = "two.sided", paired = TRUE, conf.int = TRUE))
```

```
## Error: object 'A22.high.transcripts' not found
```


Yes - transcripts that are up-regulated at high temperatures in *A22* are more highly-expressed at 19.25C in *Ar*. Note that A22 had the larger library size so if this was due to TPM not correctly accounting for differences in reads between samples, we would expect to see a positive instead of negative value here.

Next I test the converse, that transcripts that are up-regulated at low temperatures in *Ar* are more highly-expressed at 19.25C temperatures in *A22*. 


```r
# list of transcripts that are 'high' expressed in Ar
Ar.low.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Low"), 
    ]
```

```
## Error: object 'Ap.response.type' not found
```

```r

# compare expression at optimum temp (Ar.opt) between colonies using t-test
t.test(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt)
```

```
## Error: object 'Ar.low.transcripts' not found
```

```r
boxplot(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'boxplot': Error: object 'Ar.low.transcripts' not found
```

```r
# t-test with log values
t.test(log(Ar.low.transcripts$A22.opt + 1), log(Ar.low.transcripts$Ar.opt + 
    1))
```

```
## Error: object 'Ar.low.transcripts' not found
```

```r
boxplot(log(Ar.low.transcripts$A22.opt + 1), log(Ar.low.transcripts$Ar.opt + 
    1))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'boxplot': Error: object 'Ar.low.transcripts' not found
```

```r

# Wilcoxon signed rank-test
(w2 <- wilcox.test(Ar.low.transcripts$A22.opt, Ar.low.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE))
```

```
## Error: object 'Ar.low.transcripts' not found
```


Counter to expectations, expression at 19.5C is also greater in *Ar* than *A22* for transcripts upregulated at low temperatures in *Ar*. 

To confirm that there are not sample-level issues, I performed the same comparison using transcripts where I do *not* expect to see a difference in expression.


```r
# list of transcripts that are 'Intermediate' expressed in Ar
Ar.int.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate"), 
    ]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# T test
t.test(Ar.int.transcripts$A22.opt, Ar.int.transcripts$Ar.opt)
```

```
## Error: object 'Ar.int.transcripts' not found
```

```r
t.test(log(Ar.int.transcripts$A22.opt + 1), log(Ar.int.transcripts$Ar.opt + 
    1))
```

```
## Error: object 'Ar.int.transcripts' not found
```

```r
# Wilcoxon signed rank-test
wilcox.test(Ar.int.transcripts$A22.opt, Ar.int.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE)
```

```
## Error: object 'Ar.int.transcripts' not found
```

```r

# list of transcripts that are 'Intermediate' expressed in A22
A22.int.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate"), 
    ]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# T test
t.test(A22.int.transcripts$A22.opt, A22.int.transcripts$Ar.opt)
```

```
## Error: object 'A22.int.transcripts' not found
```

```r
t.test(log(A22.int.transcripts$A22.opt + 1), log(A22.int.transcripts$Ar.opt + 
    1))
```

```
## Error: object 'A22.int.transcripts' not found
```

```r
# Wilcoxon signed rank-test
wilcox.test(A22.int.transcripts$A22.opt, A22.int.transcripts$Ar.opt, alternative = "two.sided", 
    paired = TRUE, conf.int = TRUE)
```

```
## Error: object 'A22.int.transcripts' not found
```


The non-parametric test for both comparisions also finds greater constitutive expression in *Ar* than *A22*.


### Thermal tolerance indicated by region of constant gene expression ###

The 'Intermediate' expressed transcripts are core molecular processes that are expressed at non-stressful temperatures, and shut-off when the organism experiences thermall stress. We hypothesized that if the more southern *Ar* colony was more thermally-tolerant than *A22*, transcripts with 'Intermediate' expression (10-30C) would be active across a wider range of temperatures. To test this with our data, we calculated the standard deviation of the expression function for each temperature transcript that was 'Intermediate' expressed in each colony.


```r
# extract 'Intermediate' expressed transcripts for A22 colony
A22.int.lm <- RxNlmAIC[A22.int.transcripts$Transcript]
```

```
## Error: object 'RxNlmAIC' not found
```

```r
length(A22.int.lm)
```

```
## Error: object 'A22.int.lm' not found
```

```r

# apply `transcriptSD` function to all transcripts
A22.int.sd <- unlist(Map(transcriptSD, A22.int.lm, colony = "A22"))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unlist': Error in eval(expr, envir, enclos) : object 'A22.int.lm' not found
## Calls: Map -> standardGeneric -> eval -> eval -> eval
```

```r
A22.int.sd <- data.frame(colony = "A22", exp_sd = A22.int.sd)
```

```
## Error: object 'A22.int.sd' not found
```

```r

# repeat for Ar
Ar.int.lm <- RxNlmAIC[Ar.int.transcripts$Transcript]
```

```
## Error: object 'RxNlmAIC' not found
```

```r
Ar.int.sd <- unlist(Map(transcriptSD, Ar.int.lm, colony = "Ar"))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unlist': Error in eval(expr, envir, enclos) : object 'Ar.int.lm' not found
## Calls: Map -> standardGeneric -> eval -> eval -> eval
```

```r
Ar.int.sd <- data.frame(colony = "Ar", exp_sd = Ar.int.sd)
```

```
## Error: object 'Ar.int.sd' not found
```

```r

# T-test comparing standard deviation of expression between colonies
(t1 <- t.test(Ar.int.sd$exp_sd, A22.int.sd$exp_sd, alternative = "two.sided"))
```

```
## Error: object 'Ar.int.sd' not found
```


Consistent with our hypothesis, 'Intermediate' transcripts in *Ar* are expressed over a significantly wider range of temperatures than in *A22*. 


```
## Error: object 'A22.int.sd' not found
```

```
## Error: object 'Ap.int.sd' not found
```

```
## Error: object 'g1' not found
```



### Thermal sensitivity indicated by response of bimodally-expressed transcripts ###

As the converse of the above hypothesis, a colony that is especially thermally-sensitive is likely to activate expression of molecular processes more quickly. We tested this using the same approach as for the 'Intermediate' transcripts, but using the inverse of the 'Bimodal' expressed transcripts. 



```r
# extract 'Bimodal' expressed transcripts for A22 colony
A22.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$A22.type == "Bimodal"), 
    "Transcript"]]
```

```
## Error: object 'RxNlmAIC' not found
```

```r
length(A22.bim.lm)
```

```
## Error: object 'A22.bim.lm' not found
```

```r

# apply `transcriptSD` function to all transcripts
A22.bim.sd <- unlist(Map(transcriptSD, A22.bim.lm, colony = "A22"))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unlist': Error in eval(expr, envir, enclos) : object 'A22.bim.lm' not found
## Calls: Map -> standardGeneric -> eval -> eval -> eval
```

```r

# repeat for Ar
Ar.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$Ar.type == "Bimodal"), 
    "Transcript"]]
```

```
## Error: object 'RxNlmAIC' not found
```

```r
length(Ar.bim.lm)
```

```
## Error: object 'Ar.bim.lm' not found
```

```r
Ar.bim.sd <- unlist(Map(transcriptSD, Ar.bim.lm, colony = "Ar"))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unlist': Error in eval(expr, envir, enclos) : object 'Ar.bim.lm' not found
## Calls: Map -> standardGeneric -> eval -> eval -> eval
```

```r

# t-test to compare standard deviation of 'Bimodal' transcripts between
# colonies
t2 <- t.test(Ar.bim.sd, A22.bim.sd)
```

```
## Error: object 'Ar.bim.sd' not found
```


No difference in the standard deviation of expression for bimodally-expressed transcripts between colonies.


### Compare peak expression among colonies ###

Probability density function of peak expression for transcripts that differ in expression between *A22* and *Ar*


```r
# reshape data
Ap.df <- data.frame(Transcript = rep(Ap.response.type$Transcript, times = 2), 
    colony = rep(c("ApVT", "ApNC"), each = length(Ap.response.type$Transcript)), 
    max.val = c(Ap.response.type$A22.max, Ap.response.type$Ar.max))
```

```
## Error: object 'Ap.response.type' not found
```

```r

png(paste(resultsdir, "PDF_expression_all.png", sep = ""))
g3 <- ggplot(Ap.df, aes(x = max.val, fill = colony)) + geom_density(alpha = 0.2, 
    position = "identity") + # scale_fill_manual(name = 'Colony', values = c('white', 'black')) +
scale_y_continuous(name = "Density") + scale_x_continuous(name = "Temperature of maximum expression")
```

```
## Error: object 'Ap.df' not found
```

```r
suppressWarnings(print(g3))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'print': Error: object 'g3' not found
```

```r
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
varshape.out <- ldply(interaction.lms, .progress = "text", RxNvarshape)
```

```
## Error: object 'interaction.lms' not found
```

```r

# plot
boxplot(varshape.out$A22.mean, varshape.out$Ar.mean)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'boxplot': Error: object 'varshape.out' not found
```

```r
boxplot(log(varshape.out$A22.mean), log(varshape.out$Ar.mean))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'boxplot': Error: object 'varshape.out' not found
```

```r

# calculate differences in mean, slope, curvature and wiggle between
# colonies for each transcript offset = mean1 - mean2
varshape.out$offset <- varshape.out$A22.mean - varshape.out$Ar.mean
```

```
## Error: object 'varshape.out' not found
```

```r
# to compare among transcripts that differ in scale of expression,
# standardize by overall mean for each transcript
varshape.out$s.offset <- NA
```

```
## Error: object 'varshape.out' not found
```

```r
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.offset"] <- (varshape.out[i, "offset"]/mean(c(varshape.out[i, 
        "A22.mean"], varshape.out[i, "Ar.mean"])))
}
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'varshape.out' not found
```

```r

# t-test
t3 <- t.test(varshape.out$s.offset, alternative = "two.sided")
```

```
## Error: object 'varshape.out' not found
```

```r

# slope = slope1 - slope2
varshape.out$slope <- varshape.out$A22.slope - varshape.out$Ar.slope
```

```
## Error: object 'varshape.out' not found
```

```r
# to compare among transcripts that differ in scale of expression,
# standardize by overall mean for each transcript
varshape.out$s.slope <- NA
```

```
## Error: object 'varshape.out' not found
```

```r
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.slope"] <- (varshape.out[i, "slope"]/mean(c(varshape.out[i, 
        "A22.slope"], varshape.out[i, "Ar.slope"])))
}
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'varshape.out' not found
```

```r
# for transcripts where both samples had slope=0, change s.slope from Inf to
# 0 as there is truly no difference in slope
varshape.out[which(varshape.out$Ar.slope == 0 & varshape.out$A22.slope == 0), 
    "s.slope"] <- 0
```

```
## Error: object 'varshape.out' not found
```

```r

# t-test
t4 <- t.test(varshape.out$s.slope, alternative = "two.sided")
```

```
## Error: object 'varshape.out' not found
```

```r


# curvature = curvature1 - curvature2
varshape.out$curvature <- varshape.out$A22.curve - varshape.out$Ar.curve
```

```
## Error: object 'varshape.out' not found
```

```r
# to compare among transcripts that differ in scale of expression,
# standardize by overall mean for each transcript
varshape.out$s.curvature <- NA
```

```
## Error: object 'varshape.out' not found
```

```r
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.curvature"] <- (varshape.out[i, "curvature"]/mean(c(varshape.out[i, 
        "A22.curve"], varshape.out[i, "Ar.curve"])))
}
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'varshape.out' not found
```

```r
# change 'Inf' values to 0
varshape.out[which(varshape.out$Ar.curve == 0 & varshape.out$A22.curve == 0), 
    "s.curvature"] <- 0
```

```
## Error: object 'varshape.out' not found
```

```r
# t-test
tt.test(varshape.out$s.curvature, alternative = "two.sided")
```

```
## Error: could not find function "tt.test"
```

```r


# wiggle
varshape.out$wiggle <- varshape.out$A22.wiggle - varshape.out$Ar.wiggle
```

```
## Error: object 'varshape.out' not found
```

```r
varshape.out$s.wiggle <- NA
```

```
## Error: object 'varshape.out' not found
```

```r
for (i in 1:nrow(varshape.out)) {
    varshape.out[i, "s.wiggle"] <- (varshape.out[i, "wiggle"]/mean(c(varshape.out[i, 
        "A22.wiggle"], varshape.out[i, "Ar.wiggle"])))
}
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'varshape.out' not found
```

```r
# change 'Inf' values to 0
varshape.out[which(varshape.out$Ar.wiggle == 0 & varshape.out$A22.wiggle == 
    0), "s.wiggle"] <- 0
```

```
## Error: object 'varshape.out' not found
```

```r
# note that only 407 transcripts have any 'wiggle'
length(which(varshape.out$s.wiggle != 0))
```

```
## Error: object 'varshape.out' not found
```

```r

# t-test
t.test(varshape.out$s.wiggle, alternative = "two.sided")
```

```
## Error: object 'varshape.out' not found
```


Next, I partition the differences in the reaction norms into the variation explained by differences in the trait means (offset), slope, curvature and wiggle.



```r
# weighted mean-standardized values of each measure
mean(abs((varshape.out$s.offset)))
```

```
## Error: object 'varshape.out' not found
```

```r
mean(abs((varshape.out$s.slope)))
```

```
## Error: object 'varshape.out' not found
```

```r
mean(abs((varshape.out$s.curvature)))
```

```
## Error: object 'varshape.out' not found
```

```r
mean(abs((varshape.out$s.wiggle)))
```

```
## Error: object 'varshape.out' not found
```

```r

# take absolute value of each value and sum to get total differences between
# reaction norms
varshape.out$s.total <- abs(varshape.out$s.offset) + abs(varshape.out$s.slope) + 
    abs(varshape.out$s.curvature) + abs(varshape.out$s.wiggle)
```

```
## Error: object 'varshape.out' not found
```

```r

t.test(varshape.out$s.total, alternative = "two.sided")
```

```
## Error: object 'varshape.out' not found
```

```r

# variation in reaction norms due to differences in mean
varshape.out$prop.mean <- abs(varshape.out$s.offset)/varshape.out$s.total
```

```
## Error: object 'varshape.out' not found
```

```r
varshape.out$prop.slope <- abs(varshape.out$s.slope)/varshape.out$s.total
```

```
## Error: object 'varshape.out' not found
```

```r
varshape.out$prop.curve <- abs(varshape.out$s.curvature)/varshape.out$s.total
```

```
## Error: object 'varshape.out' not found
```

```r
varshape.out$prop.wiggle <- abs(varshape.out$s.wiggle)/varshape.out$s.total
```

```
## Error: object 'varshape.out' not found
```

```r

# Mean proportion and 95% CI of total variation of each measure
mean(varshape.out$prop.mean)
```

```
## Error: object 'varshape.out' not found
```

```r
quantile(varshape.out$prop.mean, probs = c(0.05, 0.5, 0.95))
```

```
## Error: object 'varshape.out' not found
```

```r

mean(varshape.out$prop.slope)
```

```
## Error: object 'varshape.out' not found
```

```r
quantile(varshape.out$prop.slope, probs = c(0.05, 0.5, 0.95))
```

```
## Error: object 'varshape.out' not found
```

```r

mean(varshape.out$prop.curve)
```

```
## Error: object 'varshape.out' not found
```

```r
quantile(varshape.out$prop.curve, probs = c(0.05, 0.5, 0.95))
```

```
## Error: object 'varshape.out' not found
```

```r

mean(varshape.out$prop.wiggle)
```

```
## Error: object 'varshape.out' not found
```

```r
quantile(varshape.out$prop.wiggle, probs = c(0.05, 0.5, 0.95))
```

```
## Error: object 'varshape.out' not found
```




## Functional annotation ##

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  


```r
setkey(annotationtable, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
signif.transcripts <- annotationtable[signif.transcripts]
setnames(signif.transcripts, "Sequence.Name", "Transcript")
```



```
## Error: object 'responsive.lms' not found
```

```
## Error: object 'responsive.lms.ann.type' not found
```

```
## Error: object 'responsive.lms.ann.type' not found
```

```
## Error: object 'responsive.lms.ann.type' not found
```



## Gene set enrichment analysis ##

I perform gene set enrichment analysis below, but a quick `grep` shows that there are 26 transcripts with GO term "response to stress", though this is not enriched compared to the frequency of this term in the full dataset.
  

```r
# GO 'response to stress' hits in responsive transcripts
GO0006950.responsive <- responsive.lms.ann.type[grep("GO:0006950", responsive.lms.ann.type$GO.Biological.Process), 
    list(Transcript, best.hit.to.nr, A22.type, Ar.type)]
```

```
## Error: object 'responsive.lms.ann.type' not found
```

```r

# in high category
GO0006950.responsive[union(with(GO0006950.responsive, grep("High", Ar.type)), 
    with(GO0006950.responsive, grep("High", A22.type))), ]
```

```
## Error: object 'GO0006950.responsive' not found
```

```r
# in low category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Low", Ar.type)), 
    with(GO0006950.responsive, grep("Low", A22.type))), ]
```

```
## Error: object 'GO0006950.responsive' not found
```

```r
# in bimodal category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Bimodal", Ar.type)), 
    with(GO0006950.responsive, grep("Bimodal", A22.type))), ]
```

```
## Error: object 'GO0006950.responsive' not found
```

```r
# in intermediate category
GO0006950.responsive[union(with(GO0006950.responsive, grep("Intermediate", Ar.type)), 
    with(GO0006950.responsive, grep("Intermediate", A22.type))), ]
```

```
## Error: object 'GO0006950.responsive' not found
```

```r

# Chi-square test to see if 'response to stress' related genes
# overrepresented in responsive.lms compared to full list
resp.stress.responsive.count <- nrow(responsive.lms.ann.type[grep("GO:0006950", 
    responsive.lms.ann.type$GO.Biological.Process), list(Transcript, best.hit.to.nr)])
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'nrow': Error: object 'responsive.lms.ann.type' not found
```

```r
# GO 'response to stress' hits in all transcripts
resp.stress.all.count <- nrow(annotationtable[grep("GO:0006950", annotationtable$GO.Biological.Process), 
    list(Sequence.Name, best.hit.to.nr)])

GO.stress.table <- matrix(rbind(resp.stress.responsive.count, nrow(responsive.lms.ann.type) - 
    resp.stress.responsive.count, resp.stress.all.count, nrow(annotationtable) - 
    resp.stress.all.count), nrow = 2)
```

```
## Error: object 'resp.stress.responsive.count' not found
```

```r

GO.stress.Xsq <- chisq.test(GO.stress.table)
```

```
## Error: object 'GO.stress.table' not found
```

```r
GO.stress.Xsq
```

```
## Error: object 'GO.stress.Xsq' not found
```





There are also 6 heat shock related genes in the responsive transcripts, out of 130 total.


```r
hsp_responsive <- responsive.lms.ann.type[grep("shock", responsive.lms.ann.type$best.hit.to.nr), 
    list(Transcript, best.hit.to.nr, A22.type, Ar.type)]
```

```
## Error: object 'responsive.lms.ann.type' not found
```

```r
hsp_responsive
```

```
## Error: object 'hsp_responsive' not found
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
# create geneList. note that NA values cause problems with topGO so set any
# NA to 1 as need to retain for GO analysis
Ap.geneList <- RxNpval$padj
Ap.geneList[which(is.na(Ap.geneList))] <- 1
stopifnot(length(which(is.na(Ap.geneList))) == 0)
names(Ap.geneList) <- RxNpval$Transcript
str(Ap.geneList)
```

```
##  Named num [1:99861] 0.00525 0.00226 0.15862 0.54076 0.0142 ...
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
## 		 the algorithm is scoring 3531 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union 
## 
## 	 Level 18:	1 nodes to be scored.
## 
## 	 Level 17:	3 nodes to be scored.
## 
## 	 Level 16:	7 nodes to be scored.
## 
## 	 Level 15:	14 nodes to be scored.
## 
## 	 Level 14:	33 nodes to be scored.
## 
## 	 Level 13:	66 nodes to be scored.
## 
## 	 Level 12:	162 nodes to be scored.
## 
## 	 Level 11:	238 nodes to be scored.
## 
## 	 Level 10:	358 nodes to be scored.
## 
## 	 Level 9:	493 nodes to be scored.
## 
## 	 Level 8:	502 nodes to be scored.
## 
## 	 Level 7:	525 nodes to be scored.
## 
## 	 Level 6:	482 nodes to be scored.
## 
## 	 Level 5:	373 nodes to be scored.
## 
## 	 Level 4:	197 nodes to be scored.
## 
## 	 Level 3:	58 nodes to be scored.
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
## 3534 GO terms scored: 17 terms with p < 0.01
## Annotation data:
##     Annotated genes: 30854 
##     Significant genes: 10593 
##     Min. no. of genes annotated to a GO: 10 
##     Nontrivial nodes: 3531
```

```r

# table results
Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, 
    topNodes = 113)
dim(Ap.BP.ResTable)
```

```
## [1] 113   6
```

```r
# pandoc.table(Ap.BP.ResTable)

# graph significant nodes

# pdf(paste(resultsdir, 'Ap.BP_topGO_sig_nodes.pdf', sep=''))
# showSigOfNodes(Ap.BP.GOdata, score(Ap.BP.resultParentChild), firstSigNodes
# = 10, useInfo = 'all') dev.off()
```


**2) High**

Genes with *High* expression in both colonies

Use `selectTranscript` function to select transcripts from 'Ap.response.type' for GSEA.


```r
# Select transcripts
Ap.high <- Ap.response.type[which(Ap.response.type$Ar.type == "High" & Ap.response.type$A22.type == 
    "High"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
Ap.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.high) <- names(Ap.geneList)
Ap.geneList.high[(which(names(Ap.geneList.high) %in% Ap.high))] <- 1
```

```
## Error: object 'Ap.high' not found
```

```r
# check correct number of values set to 1
table(Ap.geneList.high)
```

```
## Ap.geneList.high
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
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
```

```
## Error: object 'Ap.response.type' not found
```

```r
Ap.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.low) <- names(Ap.geneList)
Ap.geneList.low[(which(names(Ap.geneList.low) %in% Ap.low))] <- 1
```

```
## Error: object 'Ap.low' not found
```

```r
# check correct number of values set to 1
table(Ap.geneList.low)
```

```
## Ap.geneList.low
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
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
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
Ap.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.bim) <- names(Ap.geneList)
Ap.geneList.bim[(which(names(Ap.geneList.bim) %in% Ap.bim))] <- 1
```

```
## Error: object 'Ap.bim' not found
```

```r
# check correct number of values set to 1
table(Ap.geneList.bim)
```

```
## Ap.geneList.bim
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```

```r
# pandoc.table(Ap.bim.gsea)
write.table(Ap.bim.gsea[, "GO.ID"], file = "results/Ap_bim_gsea.txt", row.names = FALSE, 
    sep = "\t", quote = FALSE)
```



**5) Intermediate**

Genes with *Intermediate* expression in both colonies


```r
Ap.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & 
    Ap.response.type$Ar.type == "Intermediate"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'int' transcripts
Ap.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.int) <- names(Ap.geneList)
Ap.geneList.int[(which(names(Ap.geneList.int) %in% Ap.int))] <- 1
```

```
## Error: object 'Ap.int' not found
```

```r
# check correct number of values set to 1
table(Ap.geneList.int)
```

```
## Ap.geneList.int
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```

```r
# pandoc.table(Ap.int.gsea)
write.table(Ap.int.gsea[, "GO.ID"], file = "results/Ap_int_gsea.txt", row.names = FALSE, 
    sep = "\t", quote = FALSE)
```



```
## Error: undefined columns selected
```

```
## Error: object 'Ap.gsea1' not found
```

```
## Error: object 'Ap.gsea1' not found
```



Next, I perform GSEA for genes in each functional type in one colony but not the other (e.g. the set difference) to gain insight on differences between the colonies.

**A22 'High' genes not in Ar**


```r
A22.high <- Ap.response.type[which(Ap.response.type$A22.type == "High" & Ap.response.type$Ar.type != 
    "High"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.high <- rep(0, length = length(Ap.geneList))
names(A22.geneList.high) <- names(Ap.geneList)
A22.geneList.high[(which(names(A22.geneList.high) %in% A22.high))] <- 1
```

```
## Error: object 'A22.high' not found
```

```r
# check correct number of values set to 1
table(A22.geneList.high)
```

```
## A22.geneList.high
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```

```r
# pandoc.table(A22.high.gsea)
```



**Ar 'High' genes not in A22**


```r
Ar.high <- Ap.response.type[which(Ap.response.type$A22.type != "High" & Ap.response.type$Ar.type == 
    "High"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.high) <- names(Ap.geneList)
Ar.geneList.high[(which(names(Ar.geneList.high) %in% Ar.high))] <- 1
```

```
## Error: object 'Ar.high' not found
```

```r
# check correct number of values set to 1
table(Ar.geneList.high)
```

```
## Ar.geneList.high
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```



**A22 'Low' genes not in Ar**


```r
A22.low <- Ap.response.type[which(Ap.response.type$A22.type == "Low" & Ap.response.type$Ar.type != 
    "Low"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.low <- rep(0, length = length(Ap.geneList))
names(A22.geneList.low) <- names(Ap.geneList)
A22.geneList.low[(which(names(A22.geneList.low) %in% A22.low))] <- 1
```

```
## Error: object 'A22.low' not found
```

```r
# check correct number of values set to 1
table(A22.geneList.low)
```

```
## A22.geneList.low
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```



**Ar 'Low' genes not in A22**


```r
Ar.low <- Ap.response.type[which(Ap.response.type$A22.type != "Low" & Ap.response.type$Ar.type == 
    "Low"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.low) <- names(Ap.geneList)
Ar.geneList.low[(which(names(Ar.geneList.low) %in% Ar.low))] <- 1
```

```
## Error: object 'Ar.low' not found
```

```r
# check correct number of values set to 1
table(Ar.geneList.low)
```

```
## Ar.geneList.low
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```




**A22 'Bimodal' genes not in Ar**


```r
A22.bim <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal" & Ap.response.type$Ar.type != 
    "Bimodal"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.bim <- rep(0, length = length(Ap.geneList))
names(A22.geneList.bim) <- names(Ap.geneList)
A22.geneList.bim[(which(names(A22.geneList.bim) %in% A22.bim))] <- 1
```

```
## Error: object 'A22.bim' not found
```

```r
# check correct number of values set to 1
table(A22.geneList.bim)
```

```
## A22.geneList.bim
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```



**Ar 'Bimodal' genes not in A22**


```r
Ar.bim <- Ap.response.type[which(Ap.response.type$A22.type != "Bimodal" & Ap.response.type$Ar.type == 
    "Bimodal"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.bim) <- names(Ap.geneList)
Ar.geneList.bim[(which(names(Ar.geneList.bim) %in% Ar.bim))] <- 1
```

```
## Error: object 'Ar.bim' not found
```

```r
# check correct number of values set to 1
table(Ar.geneList.bim)
```

```
## Ar.geneList.bim
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```




**A22 'Intermediate' genes not in Ar**


```r
A22.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & 
    Ap.response.type$Ar.type != "Intermediate"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
A22.geneList.int <- rep(0, length = length(Ap.geneList))
names(A22.geneList.int) <- names(Ap.geneList)
A22.geneList.int[(which(names(A22.geneList.int) %in% A22.int))] <- 1
```

```
## Error: object 'A22.int' not found
```

```r
# check correct number of values set to 1
table(A22.geneList.int)
```

```
## A22.geneList.int
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```



**Ar 'Intermediate' genes not in A22**


```r
Ar.int <- Ap.response.type[which(Ap.response.type$A22.type != "Intermediate" & 
    Ap.response.type$Ar.type == "Intermediate"), "Transcript"]
```

```
## Error: object 'Ap.response.type' not found
```

```r
# create gene list, setting value to 1 for 'bim' transcripts
Ar.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.int) <- names(Ap.geneList)
Ar.geneList.int[(which(names(Ar.geneList.int) %in% Ar.int))] <- 1
```

```
## Error: object 'Ar.int' not found
```

```r
# check correct number of values set to 1
table(Ar.geneList.int)
```

```
## Ar.geneList.int
##     0 
## 99861
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
## 		 the algorithm is scoring 0 nontrivial nodes
## 		 parameters: 
## 			 test statistic:  fisher : joinFun = union
```

```
## Warning: No enrichment can pe performed - there are no feasible GO terms!
```



Make a table with results.


```
## Error: undefined columns selected
```

       GO.ID                             Term Annotated Significant
1 GO:0000002 mitochondrial genome maintenance        22           0
2 GO:0000002 mitochondrial genome maintenance        22           0
3 GO:0000002 mitochondrial genome maintenance        22           0
4 GO:0000002 mitochondrial genome maintenance        22           0
5 GO:0000002 mitochondrial genome maintenance        22           0
6 GO:0000002 mitochondrial genome maintenance        22           0
7 GO:0000002 mitochondrial genome maintenance        22           0
8 GO:0000002 mitochondrial genome maintenance        22           0
  Expected apply(l, 2, format.FUN, dig = 2, eps = 1e-30)         Type    P
1        0                                             1         High ApVT
2        0                                             1          Low ApVT
3        0                                             1 Intermediate ApVT
4        0                                             1      Bimodal ApVT
5        0                                             1         High ApNC
6        0                                             1          Low ApNC
7        0                                             1 Intermediate ApNC
8        0                                             1      Bimodal ApNC




## Visualize responsive transcripts

Organize data


```r
# extract TPM data for thermally-responsive transcripts
resp.TPM.dt.sub <- TPM.dt.sub[names(responsive.lms)]
```

```
## Error: object 'responsive.lms' not found
```

```r
setkey(resp.TPM.dt.sub, Transcript)
```

```
## Error: object 'resp.TPM.dt.sub' not found
```

```r
str(resp.TPM.dt.sub)
```

```
## Error: object 'resp.TPM.dt.sub' not found
```

```r
length(unique(resp.TPM.dt.sub$Transcript))
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'unique': Error: object 'resp.TPM.dt.sub' not found
```

```r
# scale transcripts so can compare
resp.TPM.dt.sub[, `:=`(TPM.scaled, scale(TPM)), by = Transcript]
```

```
## Error: object 'resp.TPM.dt.sub' not found
```


Make plots for all genes expressed at *High* temps in GO category "GO:0006950: response to stress"




Make plots for all genes expressed at *Low* temps in GO category "GO:0006950: response to stress"





## Shiny interactive web-app

To assist visualization of specific transcripts, I made a interactive web-app using the [shiny](http://www.rstudio.com/shiny/) package. The scripts for this app are in the sub-directory `.\ApRxN-shinyapp`.

Export data for interactive shiny app. 








```
## Error: object 'tt2' not found
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
## [1] parallel  methods   stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
##  [1] topGO_2.16.0         SparseM_1.03         GO.db_2.14.0        
##  [4] RSQLite_0.11.4       DBI_0.2-7            AnnotationDbi_1.26.0
##  [7] GenomeInfoDb_1.0.2   Biobase_2.24.0       BiocGenerics_0.10.0 
## [10] graph_1.42.0         MASS_7.3-31          plyr_1.8.1          
## [13] RCurl_1.95-4.1       bitops_1.0-6         data.table_1.9.2    
## [16] stringr_0.6.2        pander_0.3.8         knitcitations_0.5-0 
## [19] bibtex_0.3-6         ggplot2_0.9.3.1      R.utils_1.29.8      
## [22] R.oo_1.18.0          R.methodsS3_1.6.1    knitr_1.5           
## 
## loaded via a namespace (and not attached):
##  [1] codetools_0.2-8    colorspace_1.2-4   dichromat_2.0-0   
##  [4] digest_0.6.4       evaluate_0.5.3     formatR_0.10      
##  [7] grid_3.1.0         gtable_0.1.2       httr_0.3          
## [10] IRanges_1.22.3     labeling_0.2       lattice_0.20-29   
## [13] munsell_0.4.2      proto_0.3-10       RColorBrewer_1.0-5
## [16] Rcpp_0.11.1        reshape2_1.2.2     scales_0.2.3      
## [19] stats4_3.1.0       tools_3.1.0        XML_3.98-1.1      
## [22] xtable_1.7-3
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
- Courtney J. Murren, Heidi J. Maclean, Sarah E. Diamond, Ulrich K. Steiner, Mary A. Heskel, Corey A. Handelsman, Cameron K. Ghalambor, Josh R. Auld, Hilary S. Callahan, David W. Pfennig, Rick A. Relyea, Carl D. Schlichting, Joel Kingsolver,   (2014) Evolutionary Change in Continuous Reaction Norms.  *The American Naturalist*  **183**  453-467  [10.1086/675302](http://dx.doi.org/10.1086/675302)
- Robert Schmieder, Robert Edwards, Francisco Rodriguez-Valera,   (2011) Fast Identification And Removal of Sequence Contamination From Genomic And Metagenomic Datasets.  *Plos One*  **6**  e17288-NA  [10.1371/journal.pone.0017288](http://dx.doi.org/10.1371/journal.pone.0017288)
- unknown unknown,   (unknown) Unknown.  *Unknown*
- Ya Yang, Stephen A Smith,   (2013) Optimizing de Novo Assembly of Short-Read Rna-Seq Data For Phylogenomics.  *Bmc Genomics*  **14**  328-NA  [10.1186/1471-2164-14-328](http://dx.doi.org/10.1186/1471-2164-14-328)

