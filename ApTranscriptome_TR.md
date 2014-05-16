Thermal reactionome of the common ant species *Aphaenogaster picea* and *A. carolinensis*
========================================================================================
   
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**May 16, 2014**

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

Two ant colonies were used for the transcriptome sequencing. The first, designated *A22*, was collected at Molly Bog, Vermont in August 2012 by Nick Gotelli and Andrew Nguyen. The second colony, designated *Ar*, was collected by Lauren Nichols in Raleigh, North Carolina. These colonies were maintained in the lab for 6 months prior to sample collection. Bernice Bacon DeMarco (Michigan State University) identified colony *A22* as *A. picea* and *Ar* as *A. carolinensis*. For historical reasons, I refer to these species as colonies at times throughout this technical report.

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

To examine the species distribution of BLAST hits in the transcriptome assembly, I used the program [Krona](http://sourceforge.net/p/krona/home/krona/) (<a href="">unknown, unknown</a>). I ... 

~~~
KRONA code
~~~

The interactive visualization is available [here]().


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

I quantified gene expression using [sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html). To run this program, first make sure that PATHs to the software libraries are set up correctly as described on the sailfish website. 
                                                 
Then build the index of the assembly:


```r
system("sailfish index -t results/trinity-full/Trinity_cap3_uclust.fasta -o results/trinity-full/sailfish-index-Trinity-cap3-uclust -k 20 -p 4")
```


Once this is done, quantify expression for the Trimmomatic filtered reads from each species-treatment sample separately. Note that for each sample, there are three four filtered read files:

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
paired.left <- readlist[grep(".\\.paired.left.fastq$", readlist)]
paired.right <- readlist[grep("\\.paired.right.fastq$", readlist)]
unpaired.left <- readlist[grep("unpaired.left.fastq$", readlist)]
unpaired.right <- readlist[grep("unpaired.right.fastq$", readlist)]

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
## Start expression quantification for sample A22-0: 2014-05-16 08:09:31
## Done with expression quantification for sample A22-0: 2014-05-16 08:09:52
## 
## Start expression quantification for sample A22-10: 2014-05-16 08:09:52
## Done with expression quantification for sample A22-10: 2014-05-16 08:10:14
## 
## Start expression quantification for sample A22-14: 2014-05-16 08:10:14
## Done with expression quantification for sample A22-14: 2014-05-16 08:10:36
## 
## Start expression quantification for sample A22-17: 2014-05-16 08:10:36
## Done with expression quantification for sample A22-17: 2014-05-16 08:10:59
## 
## Start expression quantification for sample A22-21: 2014-05-16 08:10:59
## Done with expression quantification for sample A22-21: 2014-05-16 08:11:18
## 
## Start expression quantification for sample A22-24: 2014-05-16 08:11:18
## Done with expression quantification for sample A22-24: 2014-05-16 08:11:41
## 
## Start expression quantification for sample A22-28: 2014-05-16 08:11:41
## Done with expression quantification for sample A22-28: 2014-05-16 08:12:02
## 
## Start expression quantification for sample A22-31: 2014-05-16 08:12:02
## Done with expression quantification for sample A22-31: 2014-05-16 08:12:22
## 
## Start expression quantification for sample A22-35: 2014-05-16 08:12:22
## Done with expression quantification for sample A22-35: 2014-05-16 08:12:43
## 
## Start expression quantification for sample A22-38: 2014-05-16 08:12:43
## Done with expression quantification for sample A22-38: 2014-05-16 08:13:06
## 
## Start expression quantification for sample A22-3: 2014-05-16 08:13:06
## Done with expression quantification for sample A22-3: 2014-05-16 08:13:29
## 
## Start expression quantification for sample A22-7: 2014-05-16 08:13:29
## Done with expression quantification for sample A22-7: 2014-05-16 08:13:49
## 
## Start expression quantification for sample Ar-0: 2014-05-16 08:13:49
## Done with expression quantification for sample Ar-0: 2014-05-16 08:14:11
## 
## Start expression quantification for sample Ar-10: 2014-05-16 08:14:11
## Done with expression quantification for sample Ar-10: 2014-05-16 08:14:33
## 
## Start expression quantification for sample Ar-14: 2014-05-16 08:14:33
## Done with expression quantification for sample Ar-14: 2014-05-16 08:14:56
## 
## Start expression quantification for sample Ar-17: 2014-05-16 08:14:56
## Done with expression quantification for sample Ar-17: 2014-05-16 08:15:19
## 
## Start expression quantification for sample Ar-21: 2014-05-16 08:15:19
## Done with expression quantification for sample Ar-21: 2014-05-16 08:15:40
## 
## Start expression quantification for sample Ar-24: 2014-05-16 08:15:40
## Done with expression quantification for sample Ar-24: 2014-05-16 08:16:01
## 
## Start expression quantification for sample Ar-28: 2014-05-16 08:16:01
## Done with expression quantification for sample Ar-28: 2014-05-16 08:16:21
## 
## Start expression quantification for sample Ar-31: 2014-05-16 08:16:21
## Done with expression quantification for sample Ar-31: 2014-05-16 08:16:42
## 
## Start expression quantification for sample Ar-35: 2014-05-16 08:16:42
## Done with expression quantification for sample Ar-35: 2014-05-16 08:17:03
## 
## Start expression quantification for sample Ar-38: 2014-05-16 08:17:03
## Done with expression quantification for sample Ar-38: 2014-05-16 08:17:29
## 
## Start expression quantification for sample Ar-3: 2014-05-16 08:17:29
## Done with expression quantification for sample Ar-3: 2014-05-16 08:17:50
## 
## Start expression quantification for sample Ar-7: 2014-05-16 08:17:50
## Done with expression quantification for sample Ar-7: 2014-05-16 08:18:11
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

The TPM column for each sample was extracted and combined into a matrix for each species.




Note that expression levels at each temperature treatment are highly correlated between the two colonies.





|  Temperature  |  Correlation  |
|:-------------:|:-------------:|
|       0       |     0.99      |
|      3.5      |     0.98      |
|       7       |     0.98      |
|     10.5      |       1       |
|      14       |     0.99      |
|     17.5      |     0.98      |
|      21       |     0.98      |
|     24.5      |     0.99      |
|      28       |     0.99      |
|     31.5      |     0.99      |
|      35       |     0.99      |
|     38.5      |     0.99      |

Table: Correlations between species for gene expression at temperature treatment


**Preliminary [examination](https://minilims1.uvm.edu/BCProject-26-Cahan/methods.html#clustering-of-samples) of the data indicated that the A22_7 and Ar_7 samples may have been switched, so I remove these values from the combined expression data set for the two species.** 


```r
A22.TPM[, `:=`(colony, "A22")]
Ar.TPM[, `:=`(colony, "Ar")]
TPM.dt <- rbind(A22.TPM, Ar.TPM)
TPM.dt$colony <- as.factor(TPM.dt$colony)
str(TPM.dt)

setkey(TPM.dt, val)
TPM.dt.sub <- TPM.dt[val != 7]
unique(TPM.dt.sub$val)
length(unique(TPM.dt.sub$Transcript))
```


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
## Classes 'data.table' and 'data.frame':	2196942 obs. of  10 variables:
##  $ Transcript       : chr  "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" "0|*|Contig6267" ...
##  $ Length           : int  9990 9990 9990 9990 9990 9990 9990 9990 9990 9990 ...
##  $ TPM              : num  0.0785 0.0643 0.0357 0.0932 0.0394 ...
##  $ RPKM             : num  0.0917 0.108 0.053 0.1414 0.0492 ...
##  $ KPKM             : num  0.0917 0.108 0.053 0.1414 0.0492 ...
##  $ EstimatedNumReads: num  2957 1502 1373 2439 1711 ...
##  $ V7               : num  36.9 18.6 17.1 30.2 21.3 ...
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

$$ log(TPM + 1) = \beta_0 + \beta_1(species) + \beta_2(temp) + \beta_3(temp^2) + \beta_4(species * temp) + \beta_5(species * temp^2) + \epsilon $$

where TPM is transcripts per million. 


(1) Identify transcripts with overall significant model fit. Adjust *P* values for multiple testing using FDR and retain transcripts with FDR < 0.05. Use log-transformed response to account for outliers.



```r
# define model for RxN function
model <- "log(TPM+1) ~ colony + val + I(val^2) + colony:val + colony:I(val^2)"

# calculate overall P value and R^2 for each transcript
RxNpval <- ddply(TPM.dt.sub, .(Transcript), .inform = "TRUE", modpFunc)
```


Of the 99861 transcripts, 22343 have models with P < 0.05.

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


At the 5% FDR significance threshold, there are 10605 transcripts with an overall significant model.


(2) Fit linear model to overall significant transcripts; perform stepAIC to retain only significant terms, and save `lm` output to list


```r
# perform model selection for responsive transcripts need to use `try` to
# avoid stopping on error for AIC at Infinity
RxNlmAIC <- try(dlply(sig.TPM.dt.sub, .(Transcript), lmFunc))
```



### Grouping of thermally-responsive transcripts

The set of transcripts with significant expression patterns include those with expression that differs by species, temperature and the interaction of species and temperature. In this section, I am specifically interested in the thermally-responsive transcripts (temperature and species x temperature) so I subset the significant transcripts to examine these. 


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
temperature.lms <- other.lms[which(Map(grepFunc, other.lms, term = "val") == 
    TRUE)]
colony.lms <- other.lms[setdiff(names(other.lms), names(temperature.lms))]
responsive.lms <- c(temperature.lms, interaction.lms)
rm(other.lms)
```



|    Coefficient     |  Number.significant  |
|:------------------:|:--------------------:|
|       Total        |        10605         |
|       Colony       |         1460         |
|    Temperature     |         2264         |
| Temperature:Colony |         6881         |

Table: Number of transcripts of 99,861 total with expression that depends on species, temperature or their interaction at 5% FDR.



### Thermal-response functional types ###

The previous section identified the transcripts with thermally-responsive expression. In this section, I determine the shape of the expression response to temperature for each transcript. Catego
ries of expression response are:

* High - increase expression with temperature
* Low - decrease expression with temperature
* Intermediate - maximum expression at intermediate temperatures (14 - 28C)
* Bimodal - expressed greater than two standard deviations of expression at both low and high temperatures

I do this first for the thermally-responsive transcripts where there is no interaction with species. For the transcripts where thermal-responsive expression depends on species, I determine the functional type of the expression response separately for each species. 



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
## 'data.frame':	9145 obs. of  9 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100015|*|comp3543055_c0_seq1" "100067|*|comp3557646_c0_seq1" "100089|*|comp11313_c1_seq1" ...
##  $ A22.max   : num  38.5 0 NA 18.5 NA 0 NA 0 0 38.5 ...
##  $ A22.min   : num  0 20.5 NA 38.5 NA 23 NA 28.5 26 10.5 ...
##  $ A22.opt   : num  1.025 0.946 1 1.057 1 ...
##  $ A22.type  : chr  "High" "Bimodal" "NotResp" "Intermediate" ...
##  $ Ar.max    : num  0 38.5 0 38.5 38.5 18 0 0 NA NA ...
##  $ Ar.min    : num  38.5 0 20.5 13.5 18.5 38.5 38.5 38.5 NA NA ...
##  $ Ar.opt    : num  1.199 1.004 0.934 0.962 0.897 ...
##  $ Ar.type   : chr  "Low" "High" "Bimodal" "High" ...
```

```r

# save results to file
write.table(file = paste(resultsdir, "Ap-responsive-transcripts", Sys.Date(), 
    ".txt", sep = ""), Ap.response.type, row.names = FALSE, col.names = TRUE, 
    quote = FALSE, sep = "\t")
```


Next, I compare the number of thermally-responsive in each response category between the two colonies. 


```r
A22.type.table <- table(Ap.response.type[, "A22.type"])
Ar.type.table <- table(Ap.response.type[, "Ar.type"])

# table
Ap.type.table <- rbind(A22.type.table, Ar.type.table)

# Pearson Chi-square test
chi1 <- chisq.test(Ap.type.table)
chi1
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  Ap.type.table
## X-squared = 1259, df = 4, p-value < 2.2e-16
```

```r

# Reorganize for plotting to show overlap among categories between colonies
type.table <- table(Ap.response.type[, "Ar.type"], Ap.response.type[, "A22.type"])
# reorder
tt2 <- type.table[c(2, 4, 1, 3, 5), c(2, 4, 1, 3, 5)]

# plot
mosaicplot(t(tt2), ylab = "AcNC", xlab = "ApVT", main = "Mosaic plot of responsive transcripts")
```

![plot of chunk temperature_response](figure/temperature_response.png) 


The number of thermally-responsive in each response category differs between the colonies, with the msot transcripts expressed at *Low* temperatures in both colonies. For *ApVT*, an equal number of transcripts are expressed at *High* and *Bimodal*, followed by *Intermediate* transcripts. For *AcNC*, transcripts expressed at *Intermediate* temperatures are next most common, followed by *High* and *Bimodal*. 

Interestingly, nearly half of the *High* genes in *AcNC* are *Low* in *ApVT*, and vice versa. In contrast, most of the *Low* genes in one species are also *Low* in the other species.



|   &nbsp;   |     ApVT     |  &nbsp;  |  &nbsp;  |  &nbsp;  |    &nbsp;    |  &nbsp;  |  &nbsp;  |
|:----------:|:------------:|:--------:|:--------:|:--------:|:------------:|:--------:|:--------:|
|  **AcNC**  |              |   High   |   Low    | Bimodal  | Intermediate | NotResp  |  Total   |
|            |     High     |   310    |   444    |   149    |      65      |    93    |   1061   |
|            |     Low      |   421    |   2662   |   305    |     164      |   246    |   3798   |
|            |   Bimodal    |   167    |   176    |   296    |     113      |   108    |   860    |
|            | Intermediate |   150    |   1345   |   555    |     551      |    69    |   2670   |
|            |   NotResp    |   221    |   295    |   222    |      18      |    0     |   756    |
|            |    Total     |   1269   |   4922   |   1527   |     911      |   516    |   9145   |

Table: Number of transcripts with maximum expression at high, low, intermediate, both high and low (bimodal) temperatures or are not thermally-responsivefor each species and their overlap.


Table 4 shows the number of transcripts that fall into each expression type for each each species. The totals for each species include the 2264 transcripts that have consistent temperature responses between the two colonies. 



## Species comparison ##

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
Ap.df <- data.frame(Transcript = rep(Ap.response.type$Transcript, times = 2), 
    colony = rep(c("ApVT", "AcNC"), each = length(Ap.response.type$Transcript)), 
    max.val = c(Ap.response.type$A22.max, Ap.response.type$Ar.max))

maxexplot <- ggplot(Ap.df, aes(x = max.val, fill = colony)) + geom_density(alpha = 0.2, 
    position = "identity") + # scale_fill_manual(name = 'Colony', values = c('white', 'black')) +
scale_y_continuous(name = "Density") + scale_x_continuous(name = "Temperature of maximum expression")
suppressWarnings(print(maxexplot))
```

![plot of chunk max_exp_PDF](figure/max_exp_PDF.png) 



### Compare basal expression at optimal temperature for thermally-responsive genes

Genes upregulated in response to thermal stress in one species may have greater basal levels of expression in the other species that experiences those stressful conditions more often. To test this hypothesis, we compared expression levels near the optimal temperature (19.5C) between the two species for genes in that are either in the 'High' or 'Low' expression group in the other species. Specifically, do genes upregulated at high temperatures in *A22* have greater basal expression at optimal temperatures in *Ar*?


```r
# list of transcripts that are 'high' expressed in A22
A22.high.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == 
    "High"), ]

# Compare expression at optimum temp (A22.opt) between colonies using t-test
t.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  A22.high.transcripts$A22.opt and A22.high.transcripts$Ar.opt
## t = -0.155, df = 2396, p-value = 0.8769
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -1129   964
## sample estimates:
## mean of x mean of y 
##       372       455
```

```r
boxplot(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt)
```

![plot of chunk optimum_expression_comparison](figure/optimum_expression_comparison1.png) 

```r

# T test on log-transformed values to control for outliers
t.test(log(A22.high.transcripts$A22.opt + 1), log(A22.high.transcripts$Ar.opt + 
    1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(A22.high.transcripts$A22.opt + 1) and log(A22.high.transcripts$Ar.opt + 1)
## t = -3.69, df = 2520, p-value = 0.000232
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2630 -0.0804
## sample estimates:
## mean of x mean of y 
##      1.10      1.28
```

```r
boxplot(log(A22.high.transcripts$A22.opt + 1), log(A22.high.transcripts$Ar.opt + 
    1))
```

![plot of chunk optimum_expression_comparison](figure/optimum_expression_comparison2.png) 


The `t.test` fails to account for the many orders of magnitude difference in expression among transcripts, e.g. non-equal variances. This problem is the key issue in the analysis of differential expression (<a href="http://dx.doi.org/10.1186/1471-2105-11-94">Bullard et al. 2010</a>; <a href="http://dx.doi.org/10.1038/nprot.2013.099">Anders et al. 2013</a>). As my goal is simply to determine if expression is typically greater at optimal temperatures (19.5 C) in *Ar* than *A22* for genes that are up-regulated at high temperatures in *A22*, I use a non-parametric Wilcoxon signed rank-test


```r
w1 <- wilcox.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt, 
    alternative = "two.sided", paired = TRUE, conf.int = TRUE)
w1
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  A22.high.transcripts$A22.opt and A22.high.transcripts$Ar.opt
## V = 209676, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.300 -0.172
## sample estimates:
## (pseudo)median 
##         -0.234
```


Consistent with expectation, there is greater expression at 19.5C for *Ar* than *A22* transcripts for the set of transcripts that are transcripts that are up-regulated at high temperatures in *A22*. Note that A22 had the larger library size so if this was due to TPM not correctly accounting for differences in reads between samples, we would expect to see a positive instead of negative value here.

Next I test the converse: do genes up-regulated at low temperatures in *Ar* have greater basal expression near optimal temperatures in *A22*?


```r
# list of transcripts that are 'high' expressed in Ar
Ar.low.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Low"), 
    ]

# t-test with log values
t.test(log(Ar.low.transcripts$A22.opt + 1), log(Ar.low.transcripts$Ar.opt + 
    1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(Ar.low.transcripts$A22.opt + 1) and log(Ar.low.transcripts$Ar.opt + 1)
## t = -6.58, df = 7584, p-value = 5.016e-11
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1726 -0.0934
## sample estimates:
## mean of x mean of y 
##      1.26      1.39
```

```r
boxplot(log(Ar.low.transcripts$A22.opt + 1), log(Ar.low.transcripts$Ar.opt + 
    1))
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
## V = 2240946, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.459 -0.365
## sample estimates:
## (pseudo)median 
##         -0.411
```


Counter to expectations, expression at optimal temperatures is also greater in *Ar* than *A22* for transcripts upregulated at low temperatures in *Ar*. 

To confirm that there are not sample-level issues, I performed the same comparison using transcripts where I do *not* expect to see a difference in expression.


```r
# list of transcripts that are 'Intermediate' expressed in Ar
Ar.int.transcripts <- Ap.response.type[which(Ap.response.type$Ar.type == "Intermediate"), 
    ]
# T test
t.test(log(Ar.int.transcripts$A22.opt + 1), log(Ar.int.transcripts$Ar.opt + 
    1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(Ar.int.transcripts$A22.opt + 1) and log(Ar.int.transcripts$Ar.opt + 1)
## t = -7.16, df = 5332, p-value = 9.291e-13
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.281 -0.160
## sample estimates:
## mean of x mean of y 
##      1.37      1.59
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
## V = 1072644, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.619 -0.476
## sample estimates:
## (pseudo)median 
##         -0.545
```

```r

# list of transcripts that are 'Intermediate' expressed in A22
A22.int.transcripts <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate"), 
    ]
# T test
t.test(log(A22.int.transcripts$A22.opt + 1), log(A22.int.transcripts$Ar.opt + 
    1))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  log(A22.int.transcripts$A22.opt + 1) and log(A22.int.transcripts$Ar.opt + 1)
## t = -1.69, df = 1820, p-value = 0.09127
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2085  0.0155
## sample estimates:
## mean of x mean of y 
##      1.41      1.50
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
## V = 149648, p-value = 5.29e-13
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.407 -0.218
## sample estimates:
## (pseudo)median 
##         -0.305
```


The non-parametric test for both comparisions also finds greater expression in *Ar* than *A22* at the optimal temperature.


### Compare thermal stability of *Intermediate* genes

*Intermediate* genes are those that have expression that is shut-down as conditions become stressful, likely non-essential molecular processes. We hypothesized that if the more southern *Ar* species was more thermally-tolerant than *A22*, transcripts with *Intermediate* expression (10-30C) would be active across a wider range of temperatures. To test this with our data, we calculated the standard deviation of the expression function for each temperature transcript that was 'Intermediate' expressed in each species.


```r
# extract 'Intermediate' expressed transcripts for A22
A22.int.lm <- RxNlmAIC[A22.int.transcripts$Transcript]
length(A22.int.lm)
```

```
## [1] 911
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
(t.thermbreadth <- t.test(Ar.int.sd$exp_sd, A22.int.sd$exp_sd, alternative = "two.sided"))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.int.sd$exp_sd and A22.int.sd$exp_sd
## t = 9.38, df = 1301, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.350 0.534
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
A22.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$A22.type == "Bimodal"), 
    "Transcript"]]
length(A22.bim.lm)
```

```
## [1] 1527
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
## [1] 860
```

```r
Ar.bim.sd <- unlist(Map(transcriptSD, Ar.bim.lm, colony = "Ar"))

# t-test to compare standard deviation of 'Bimodal' transcripts between
# colonies
(t.themsens <- t.test(Ar.bim.sd, A22.bim.sd))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.bim.sd and A22.bim.sd
## t = -0.853, df = 1494, p-value = 0.3938
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.213  0.084
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
## Classes 'data.table' and 'data.frame':	201190 obs. of  10 variables:
##  $ Transcript       : chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ Length           : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM              : num  0 0.3557 0 0.7426 0.0742 ...
##  $ RPKM             : num  0 0.5975 0 1.1269 0.0927 ...
##  $ KPKM             : num  0 0.5975 0 1.1269 0.0927 ...
##  $ EstimatedNumReads: num  0 157.5 0 368.6 61.1 ...
##  $ V7               : num  0 1.951 0 4.566 0.763 ...
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
## [1] 9145
```

```r
# scale transcripts so can compare
resp.TPM.dt.sub[, `:=`(TPM.scaled, scale(TPM)), by = Transcript]
```

```
##                          Transcript Length    TPM   RPKM   KPKM
##      1: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      2: 100008|*|comp137625_c0_seq2    208 0.3557 0.5975 0.5975
##      3: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      4: 100008|*|comp137625_c0_seq2    208 0.7426 1.1269 1.1269
##      5: 100008|*|comp137625_c0_seq2    208 0.0742 0.0927 0.0927
##     ---                                                        
## 201186:      9|*|comp147140_c0_seq1   9030 0.7222 1.2184 1.2184
## 201187:      9|*|comp147140_c0_seq1   9030 0.7138 0.9856 0.9856
## 201188:      9|*|comp147140_c0_seq1   9030 0.4376 0.7972 0.7972
## 201189:      9|*|comp147140_c0_seq1   9030 0.5647 0.7729 0.7729
## 201190:      9|*|comp147140_c0_seq1   9030 0.4604 0.8410 0.8410
##         EstimatedNumReads      V7 sample  val colony TPM.scaled
##      1:               0.0   0.000  A22-0  0.0    A22     -0.563
##      2:             157.5   1.951   Ar-0  0.0     Ar      1.117
##      3:               0.0   0.000  A22-3  3.5    A22     -0.563
##      4:             368.6   4.566   Ar-3  3.5     Ar      2.946
##      5:              61.1   0.763 A22-10 10.5    A22     -0.213
##     ---                                                        
## 201186:           17858.3 220.811  Ar-31 31.5     Ar     -0.265
## 201187:           23576.4 294.535 A22-35 35.0    A22     -0.292
## 201188:            8819.6 109.224  Ar-35 35.0     Ar     -1.186
## 201189:           23755.8 296.656 A22-38 38.5    A22     -0.775
## 201190:           12661.3 156.444  Ar-38 38.5     Ar     -1.112
```

```r
# rename colonies
resp.TPM.dt.sub$colony2 <- ifelse(resp.TPM.dt.sub$colony == "A22", "ApVT", "AcNC")
```


Predict expression for responsive transcripts.


```r
# apply predFunc to all responsive transcripts
resp.TPM.dt.sub.pred <- ddply(resp.TPM.dt.sub, .(Transcript), .inform = "TRUE", 
    predFunc)

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

A22.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & 
    Ap.response.type$Ar.type != "Intermediate"), "Transcript"]
Ar.int <- Ap.response.type[which(Ap.response.type$A22.type != "Intermediate" & 
    Ap.response.type$Ar.type == "Intermediate"), "Transcript"]
```



Calculate T~on~ for *High* genes in each species.


```r
# transcripts expressed at *High* and *Bimodal* temperatures in A22
A22.high.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(A22.high, A22.bim), "A22")]
setkey(A22.high.TPM.dt.sub, Transcript)
str(A22.high.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	24090 obs. of  5 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0 0 0.0742 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.01 1.01 1.01 1.01 1.01 ...
##  - attr(*, ".internal.selfref")=<externalptr> 
##  - attr(*, "sorted")= chr "Transcript"
```

```r

# make data.frame for results
l1 <- length(unique(A22.high.TPM.dt.sub$Transcript))
A22.high.T_on <- data.frame(Transcript = unique(A22.high.TPM.dt.sub$Transcript), 
    colony = rep("ApVT", length = l1), type = rep("High", length = l1), T_on = NA)

# loop across transcripts, calculating T_on

for (i in unique(A22.high.TPM.dt.sub$Transcript)) {
    subdf <- A22.high.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val > 14), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM)))) + 1, 
        val]
    A22.high.T_on[which(A22.high.T_on$Transcript == i), "T_on"] <- T_on
}

# repeat for Ar
Ar.high.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(Ar.high, Ar.bim), "Ar")]
setkey(Ar.high.TPM.dt.sub, Transcript)
str(Ar.high.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	14465 obs. of  5 variables:
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
Ar.high.T_on <- data.frame(Transcript = unique(Ar.high.TPM.dt.sub$Transcript), 
    colony = rep("AcNC", length = l2), type = rep("High", length = l2), T_on = NA)

for (i in unique(Ar.high.TPM.dt.sub$Transcript)) {
    subdf <- Ar.high.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val > 14), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM)))) + 1, 
        val]
    Ar.high.T_on[which(Ar.high.T_on$Transcript == i), "T_on"] <- T_on
}

# determine if T~on~ is greater in *A22* or *Ar* for *High* genes.

t.test(Ar.high.T_on$T_on, A22.high.T_on$T_on)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.high.T_on$T_on and A22.high.T_on$T_on
## t = 4.16, df = 2598, p-value = 0.00003236
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.389 1.083
## sample estimates:
## mean of x mean of y 
##      33.6      32.9
```


Genes with increased expression at *High* temperatures are on average turned on 1°C higher in *AcNC* than *ApVT*.

Repeat analysis for *Low* genes.


```r
# transcripts expressed at *Low* temperatures in A22
A22.low.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(A22.low, "A22")]
str(A22.low.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	24860 obs. of  5 variables:
##  $ Transcript: chr  "100191|*|comp2883648_c0_seq1" "100191|*|comp2883648_c0_seq1" "100191|*|comp2883648_c0_seq1" "100191|*|comp2883648_c0_seq1" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0.482 0.298 0 0.418 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.46 1.33 1.15 1.09 1.05 ...
##  - attr(*, "sorted")= chr  "Transcript" "colony"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r

# make data.frame for results
l3 <- length(unique(A22.low.TPM.dt.sub$Transcript))
A22.low.T_on <- data.frame(Transcript = unique(A22.low.TPM.dt.sub$Transcript), 
    colony = rep("ApVT", length = l3), type = rep("Low", length = l3), T_on = NA)

# loop across transcripts, calculating T_on

for (i in unique(A22.low.TPM.dt.sub$Transcript)) {
    subdf <- A22.low.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val < 21), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM)))) + 1, 
        val]
    A22.low.T_on[which(A22.low.T_on$Transcript == i), "T_on"] <- T_on
}

# repeat for Ar
Ar.low.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(Ar.low, "Ar")]
str(Ar.low.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	12496 obs. of  5 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 2 2 2 2 2 2 2 2 2 2 ...
##  $ TPM       : num  0.356 0.743 0.659 0.148 0.094 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.62 1.49 1.3 1.23 1.17 ...
##  - attr(*, "sorted")= chr  "Transcript" "colony"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r

l4 <- length(unique(Ar.low.TPM.dt.sub$Transcript))
Ar.low.T_on <- data.frame(Transcript = unique(Ar.low.TPM.dt.sub$Transcript), 
    colony = rep("AcNC", length = l4), type = rep("Low", length = length(unique(Ar.low.TPM.dt.sub$Transcript))), 
    T_on = NA)

for (i in unique(Ar.low.TPM.dt.sub$Transcript)) {
    subdf <- Ar.low.TPM.dt.sub[i]
    subdf <- subdf[which(subdf$val < 21), ]
    T_on <- subdf[median(which(diff(subdf$TPM) == max(diff(subdf$TPM)))) + 1, 
        val]
    Ar.low.T_on[which(Ar.low.T_on$Transcript == i), "T_on"] <- T_on
}

t.test(Ar.low.T_on$T_on, A22.low.T_on$T_on)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.low.T_on$T_on and A22.low.T_on$T_on
## t = -5.72, df = 1831, p-value = 1.228e-08
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -1.454 -0.712
## sample estimates:
## mean of x mean of y 
##      11.7      12.8
```


Genes with increased expression at *Low* temperatures are on average turned on 1°C higher in *ApVT* than *AcNC*.

Visualize T~on~ for both *Low* and *High* genes on the same plot.

![plot of chunk plot_low_T_on](figure/plot_low_T_on.png) 


### Evaluate the extent to which differences in thermal reaction norms are due to mean or shape

For the transcripts that differed in thermal responsiveness due to temperature, was the difference primarily due to differences in the mean value of expression, slope, curvature of a higher order effect? To test this, I rougly follow <a href="http://dx.doi.org/10.1086/675302">Murren et al. (2014)</a> by defining differences among reation norms for individual genes due to changes in the overall mean, slope, curvature and all higher-order shape differences (i.e. wiggle).

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
## t = -10.9, df = 6880, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1079 -0.0749
## sample estimates:
## mean of the differences 
##                 -0.0914
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
t.slope <- t.test(log(varshape.out$A22.slope + 0.1), log(varshape.out$Ar.slope + 
    0.1), paired = TRUE)
t.slope
```

```
## 
## 	Paired t-test
## 
## data:  log(varshape.out$A22.slope + 0.1) and log(varshape.out$Ar.slope + 0.1)
## t = 7.04, df = 6880, p-value = 2.098e-12
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.0186 0.0330
## sample estimates:
## mean of the differences 
##                  0.0258
```

```r

# delta curvature
boxplot(log(varshape.out$A22.curve), log(varshape.out$Ar.curve))
```

![plot of chunk varshape](figure/varshape4.png) 

```r
t.curvature <- t.test(log(varshape.out$A22.curve + 0.1), log(varshape.out$Ar.curve + 
    0.1), paired = TRUE)
t.curvature
```

```
## 
## 	Paired t-test
## 
## data:  log(varshape.out$A22.curve + 0.1) and log(varshape.out$Ar.curve + 0.1)
## t = 0.379, df = 6880, p-value = 0.7051
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.00202  0.00299
## sample estimates:
## mean of the differences 
##                0.000483
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
## t = -0.984, df = 6880, p-value = 0.3252
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.000815  0.000270
## sample estimates:
## mean of the differences 
##               -0.000272
```


Reaction norms differ between species by mean and slope, but not by curvature or wiggle.

Next, I partition the differences in the reaction norms into the variation explained by differences in the trait means and slope, as curvature and wiggle did not differ between species. To do this, I first standardize each of these parameters so they are directly comparable. 


```r
# calculate differences in mean, slope, curvature and wiggle between
# colonies for each transcript mean = mean1 - mean2
varshape.out$lmean <- log(varshape.out$A22.mean) - log(varshape.out$Ar.mean)

# slope = slope1 - slope2
varshape.out$lslope <- log(varshape.out$A22.slope + 0.1) - log(varshape.out$Ar.slope + 
    0.1)

# take absolute value of each value and sum to get total differences between
# reaction norms
varshape.out$ltotal <- abs(varshape.out$lmean) + abs(varshape.out$lslope)

# variation in reaction norms due to differences in mean
varshape.out$prop.lmean <- abs(varshape.out$lmean)/varshape.out$ltotal
varshape.out$prop.lslope <- abs(varshape.out$lslope)/varshape.out$ltotal


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
## 0.455 0.815 0.985
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
## 0.0152 0.1848 0.5454
```


From this analysis, about 3/4 of the differences in reaction norms between species are due to changes in the mean, with the remainder being due to changes in slope.




## Gene set enrichment analysis for thermal-responsive genes

### Functional annotation

In the previous section, I identified transcripts that show significant responses in expression. Next, I add gene annotation and ontology information to these transcripts.  


```r
setkey(annotationtable, Sequence.Name)
signif.transcripts <- data.table(signif.transcripts)
setkey(signif.transcripts, Transcript)
signif.transcripts <- annotationtable[signif.transcripts]
setnames(signif.transcripts, "Sequence.Name", "Transcript")
```



```r
responsive.lms.ann.type <- signif.transcripts[names(responsive.lms)]
# combine responsive.lms.ann with 'type' information
responsive.lms.ann.type <- merge(responsive.lms.ann.type, Ap.response.type[, 
    c("Transcript", "A22.type", "Ar.type")], by = "Transcript", all = TRUE)
str(responsive.lms.ann.type)
```

```
## Classes 'data.table' and 'data.frame':	9145 obs. of  17 variables:
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
##  $ pval                 : num  0.003629 0.000908 0.0021 0.001646 0.001552 ...
##  $ adj.r.squared        : num  0.522 0.604 0.557 0.571 0.574 ...
##  $ padj                 : num  0.0379 0.0135 0.0254 0.0211 0.0202 ...
##  $ A22.type             : chr  "High" "Bimodal" "NotResp" "Intermediate" ...
##  $ Ar.type              : chr  "Low" "High" "Bimodal" "High" ...
##  - attr(*, "sorted")= chr "Transcript"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
write.csv(responsive.lms.ann.type, file = paste(resultsdir, "Ap_responsive_genes_", 
    Sys.Date(), ".csv", sep = ""), row.names = FALSE)
```


I perform gene set enrichment analysis below, but a quick `grep` shows that there are 26 transcripts with GO term "response to stress", though this is not enriched compared to the frequency of this term in the full dataset.
  

```r
# GO 'response to stress' hits in responsive transcripts
GO0006950.responsive <- responsive.lms.ann.type[grep("GO:0006950", responsive.lms.ann.type$GO.Biological.Process), 
    list(Transcript, best.hit.to.nr, A22.type, Ar.type)]

# in high category
GO0006950.responsive[union(with(GO0006950.responsive, grep("High", Ar.type)), 
    with(GO0006950.responsive, grep("High", A22.type))), ]
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
GO0006950.responsive[union(with(GO0006950.responsive, grep("Low", Ar.type)), 
    with(GO0006950.responsive, grep("Low", A22.type))), ]
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
## 22:          Low         High
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
## 6: 15115|*|comp132715_c0_seq1
## 7: 21384|*|comp149042_c0_seq3
##                                                                                   best.hit.to.nr
## 1: gi|332019420|gb|EGI59904.1| Putative fat-like cadherin-related tumor suppressor-like protein 
## 2:                                  gi|322796169|gb|EFZ18745.1| hypothetical protein SINV_07491 
## 3:                              gi|332018201|gb|EGI58806.1| Protein lethal(2)essential for life 
## 4:                          gi|367054010|ref|XP_003657383.1| hypothetical protein THITE_2156506 
## 5:                                    gi|493322437|ref|WP_006279741.1| molecular chaperone DnaK 
## 6:                                            gi|194716766|gb|ACF93232.1| heat shock protein 90 
## 7:                         gi|396467618|ref|XP_003837992.1| hypothetical protein LEMA_P120390.1 
##    A22.type Ar.type
## 1:      Low Bimodal
## 2:  Bimodal Bimodal
## 3:     High Bimodal
## 4:     High Bimodal
## 5:  NotResp Bimodal
## 6:  Bimodal    High
## 7:  Bimodal NotResp
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

# Chi-square test to see if 'response to stress' related genes
# overrepresented in responsive.lms compared to full list
resp.stress.responsive.count <- nrow(responsive.lms.ann.type[grep("GO:0006950", 
    responsive.lms.ann.type$GO.Biological.Process), list(Transcript, best.hit.to.nr)])
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
## X-squared = 2.1, df = 1, p-value = 0.1469
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

Perform GSEA for all thermally-responsive transcripts.


```r
# create geneList. note that NA values cause problems with topGO so set any
# NA to 1 as need to retain for GO analysis
Ap.geneList <- RxNpval$padj
Ap.geneList[which(is.na(Ap.geneList))] <- 1
stopifnot(length(which(is.na(Ap.geneList))) == 0)
names(Ap.geneList) <- RxNpval$Transcript
str(Ap.geneList)

# Function to select top genes (defined above)
selectFDR <- function(padj) {
    return(padj < 0.05)
}

# create topGOdata object
Ap.BP.GOdata <- new("topGOdata", description = "BP gene set analysis", ontology = "BP", 
    allGenes = Ap.geneList, geneSel = selectFDR, nodeSize = 10, annot = annFUN.gene2GO, 
    gene2GO = geneID2GO)

# Ap.BP.GOdata

# perform enrichment analysis using parentchild method
Ap.BP.resultParentChild <- runTest(Ap.BP.GOdata, statistic = "fisher", algorithm = "parentchild")
Ap.BP.resultParentChild

# table results
Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, 
    topNodes = 10)
Ap.BP.ResTable
```


As the molecular processes involved in cold and hot temperature tolerance are different, I perform GSEA separately for each response type and compile results in a single table.


```r
## Genes with *High* expression in both colonies
Ap.high <- Ap.response.type[which(Ap.response.type$Ar.type == "High" & Ap.response.type$A22.type == 
    "High"), "Transcript"]
Ap.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.high) <- names(Ap.geneList)
Ap.geneList.high[(which(names(Ap.geneList.high) %in% Ap.high))] <- 1
# check correct number of values set to 1
table(Ap.geneList.high)
# run GSEA
Ap.high.gsea <- gsea(genelist = Ap.geneList.high, geneID2GO = geneID2GO, plotpath = NA)

## Genes with *Low* expression in both colonies
Ap.low <- Ap.response.type[which(Ap.response.type$Ar.type == "Low" & Ap.response.type$A22.type == 
    "Low"), "Transcript"]
Ap.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.low) <- names(Ap.geneList)
Ap.geneList.low[(which(names(Ap.geneList.low) %in% Ap.low))] <- 1
# check correct number of values set to 1
table(Ap.geneList.low)
# Run GSEA
Ap.low.gsea <- gsea(genelist = Ap.geneList.low, geneID2GO = geneID2GO)

## Genes with *Bimodal* expression in both colonies
Ap.bim <- Ap.response.type[which(Ap.response.type$A22.type == "Bimodal" & Ap.response.type$Ar.type == 
    "Bimodal"), "Transcript"]
# create gene list, setting value to 1 for 'bim' transcripts
Ap.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.bim) <- names(Ap.geneList)
Ap.geneList.bim[(which(names(Ap.geneList.bim) %in% Ap.bim))] <- 1
# check correct number of values set to 1
table(Ap.geneList.bim)
# Run GSEA
Ap.bim.gsea <- gsea(genelist = Ap.geneList.bim, geneID2GO = geneID2GO)

## Genes with *Intermediate* expression in both colonies
Ap.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & 
    Ap.response.type$Ar.type == "Intermediate"), "Transcript"]
# create gene list, setting value to 1 for 'int' transcripts
Ap.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.int) <- names(Ap.geneList)
Ap.geneList.int[(which(names(Ap.geneList.int) %in% Ap.int))] <- 1
# check correct number of values set to 1
table(Ap.geneList.int)
# Run GSEA
Ap.int.gsea <- gsea(genelist = Ap.geneList.int, geneID2GO = geneID2GO)

## merge results into single table
Ap.high.gsea$Type <- "High"
Ap.low.gsea$Type <- "Low"
Ap.bim.gsea$Type <- "Bimodal"
Ap.int.gsea$Type <- "Intermediate"
# combine
Ap.gsea <- rbind(Ap.high.gsea, Ap.low.gsea, Ap.bim.gsea, Ap.int.gsea)
# reorder
Ap.gsea1 <- Ap.gsea[, c("Type", "GO.ID", "Term", "Annotated", "Significant", 
    "Expected", "parentchild")]
colnames(Ap.gsea1)[7] <- "P"
# write to file
write.csv(Ap.setdiff.gsea, file = paste(resultsdir, "Ap_setdiff_GSEA_", Sys.Date(), 
    ".csv", sep = ""), row.names = FALSE)
```

```
## Error: object 'Ap.setdiff.gsea' not found
```



-----------------------------------------------------
    Type       GO.ID                Term             
------------ ---------- -----------------------------
    High     GO:0009889  regulation of biosynthetic  
                                   process           

    High     GO:0031326    regulation of cellular    
                            biosynthetic proc...     

    High     GO:0019222    regulation of metabolic   
                                   process           

    High     GO:0031323    regulation of cellular    
                              metabolic process      

    High     GO:0050794    regulation of cellular    
                                   process           

    High     GO:0050789   regulation of biological   
                                   process           

    High     GO:0010556  regulation of macromolecule 
                               biosynthetic...       

    Low      GO:0019538   protein metabolic process  

    Low      GO:0043412  macromolecule modification  

    Low      GO:0070085         glycosylation        

    Low      GO:0044267  cellular protein metabolic  
                                   process           

    Low      GO:0007264 small GTPase mediated signal 
                               transductio...        

    Low      GO:0009141    nucleoside triphosphate   
                             metabolic proces...     

    Low      GO:0006643   membrane lipid metabolic   
                                   process           

    Low      GO:0009259   ribonucleotide metabolic   
                                   process           

    Low      GO:0033036  macromolecule localization  

    Low      GO:0036211 protein modification process 

    Low      GO:1901292     nucleoside phosphate     
                              catabolic process      

    Low      GO:0045184   establishment of protein   
                                localization         

    Low      GO:0009101   glycoprotein biosynthetic  
                                   process           

    Low      GO:0006486     protein glycosylation    

    Low      GO:0006163  purine nucleotide metabolic 
                                   process           

    Low      GO:0009166 nucleotide catabolic process 

    Low      GO:0009581     detection of external    
                                  stimulus           

    Low      GO:0009100    glycoprotein metabolic    
                                   process           

    Low      GO:0055001    muscle cell development   

    Low      GO:0009966     regulation of signal     
                                transduction         

    Low      GO:0043170    macromolecule metabolic   
                                   process           

    Low      GO:0006664 glycolipid metabolic process 

    Low      GO:0006413   translational initiation   

    Low      GO:0050650      chondroitin sulfate     
                           proteoglycan biosynt...   

    Low      GO:0030206      chondroitin sulfate     
                            biosynthetic process     

    Low      GO:0051259    protein oligomerization   

    Low      GO:0070271  protein complex biogenesis  

    Low      GO:0006464 cellular protein modification
                                   process           

    Low      GO:0006665    sphingolipid metabolic    
                                   process           

    Low      GO:0010608      posttranscriptional     
                           regulation of gene e...   

    Low      GO:0032318   regulation of Ras GTPase   
                                  activity           

    Low      GO:0006140   regulation of nucleotide   
                             metabolic proce...      

    Low      GO:0061326   renal tubule development   

    Low      GO:0046039     GTP metabolic process    

    Low      GO:0009118   regulation of nucleoside   
                             metabolic proce...      

    Low      GO:0009894    regulation of catabolic   
                                   process           

    Low      GO:0009582 detection of abiotic stimulus

    Low      GO:0046907    intracellular transport   

    Low      GO:0050954     sensory perception of    
                            mechanical stimulu...    

    Low      GO:1901068 guanosine-containing compound
                                metabolic ...        

    Low      GO:0006417   regulation of translation  

    Low      GO:0001894      tissue homeostasis      

    Low      GO:0042692  muscle cell differentiation 

    Low      GO:0006184     GTP catabolic process    

    Low      GO:0031329    regulation of cellular    
                              catabolic process      

    Low      GO:0044260    cellular macromolecule    
                              metabolic process      

    Low      GO:1901069 guanosine-containing compound
                                catabolic ...        

    Low      GO:0048871   multicellular organismal   
                                 homeostasis         

    Low      GO:0072523  purine-containing compound  
                              catabolic pro...       

    Low      GO:0001655 urogenital system development

    Low      GO:0006446  regulation of translational 
                                 initiation          

    Low      GO:0032273    positive regulation of    
                            protein polymeriz...     

    Low      GO:0065003    macromolecular complex    
                                  assembly           

    Low      GO:0019220    regulation of phosphate   
                             metabolic proces...     

    Low      GO:0016573      histone acetylation     

    Low      GO:0010646      regulation of cell      
                                communication        

    Low      GO:0072002 Malpighian tubule development

    Low      GO:0050654      chondroitin sulfate     
                           proteoglycan metabol...   

    Low      GO:0030204 chondroitin sulfate metabolic
                                   process           

    Low      GO:0030811   regulation of nucleotide   
                             catabolic proce...      

    Low      GO:0033559    unsaturated fatty acid    
                              metabolic process      

    Low      GO:0019685 photosynthesis, dark reaction

    Low      GO:0015031       protein transport      

    Low      GO:0023051    regulation of signaling   

    Low      GO:0060042    retina morphogenesis in   
                               camera-type eye       

    Low      GO:0048583   regulation of response to  
                                  stimulus           

  Bimodal    GO:0065007     biological regulation    

  Bimodal    GO:0050789   regulation of biological   
                                   process           

  Bimodal    GO:0050794    regulation of cellular    
                                   process           

  Bimodal    GO:0031323    regulation of cellular    
                              metabolic process      

  Bimodal    GO:0009889  regulation of biosynthetic  
                                   process           

  Bimodal    GO:0080090     regulation of primary    
                              metabolic process      

  Bimodal    GO:0019222    regulation of metabolic   
                                   process           

  Bimodal    GO:0031326    regulation of cellular    
                            biosynthetic proc...     

  Bimodal    GO:0010468 regulation of gene expression

  Bimodal    GO:0044700   single organism signaling  

  Bimodal    GO:0010556  regulation of macromolecule 
                               biosynthetic...       

  Bimodal    GO:0036211 protein modification process 

  Bimodal    GO:0006351 transcription, DNA-templated 

  Bimodal    GO:0051171    regulation of nitrogen    
                            compound metaboli...     

  Bimodal    GO:2000112    regulation of cellular    
                            macromolecule bio...     

  Bimodal    GO:0023052           signaling          

  Bimodal    GO:0050896     response to stimulus     

  Bimodal    GO:0034654     nucleobase-containing    
                            compound biosynthe...    

  Bimodal    GO:0019219         regulation of        
                        nucleobase-containing comp...

  Bimodal    GO:0060255  regulation of macromolecule 
                               metabolic pr...       

  Bimodal    GO:0007154      cell communication      

  Bimodal    GO:0032774   RNA biosynthetic process   

  Bimodal    GO:0009069   serine family amino acid   
                             metabolic proce...      

Intermediate GO:0006664 glycolipid metabolic process 

Intermediate GO:0046467  membrane lipid biosynthetic 
                                   process           

Intermediate GO:0009100    glycoprotein metabolic    
                                   process           

Intermediate GO:0009101   glycoprotein biosynthetic  
                                   process           

Intermediate GO:0006486     protein glycosylation    

Intermediate GO:0009247    glycolipid biosynthetic   
                                   process           

Intermediate GO:0006643   membrane lipid metabolic   
                                   process           

Intermediate GO:0071103    DNA conformation change   

Intermediate GO:0043413  macromolecule glycosylation 

Intermediate GO:0006323         DNA packaging        

Intermediate GO:0034728    nucleosome organization   

Intermediate GO:0031327    negative regulation of    
                            cellular biosynth...     

Intermediate GO:0009890    negative regulation of    
                            biosynthetic proc...     

Intermediate GO:0006665    sphingolipid metabolic    
                                   process           

Intermediate GO:0051172    negative regulation of    
                            nitrogen compound...     

Intermediate GO:0071824  protein-DNA complex subunit 
                                organization         

Intermediate GO:0070085         glycosylation        

Intermediate GO:0044255   cellular lipid metabolic   
                                   process           

Intermediate GO:0044267  cellular protein metabolic  
                                   process           

Intermediate GO:0042158   lipoprotein biosynthetic   
                                   process           

Intermediate GO:0005975    carbohydrate metabolic    
                                   process           
-----------------------------------------------------

Table: Enriched GO terms for each thermal-response category


Next, I perform GSEA for genes in each functional type in one species but not the other (e.g. the set difference) to gain insight on differences between the colonies.



```r
# A22 'High' genes not in Ar**
A22.geneList.high <- rep(0, length = length(Ap.geneList))
names(A22.geneList.high) <- names(Ap.geneList)
A22.geneList.high[(which(names(A22.geneList.high) %in% A22.high))] <- 1
A22.high.gsea <- gsea(genelist = A22.geneList.high, geneID2GO = geneID2GO)

# *Ar 'High' genes not in A22**
Ar.geneList.high <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.high) <- names(Ap.geneList)
Ar.geneList.high[(which(names(Ar.geneList.high) %in% Ar.high))] <- 1
Ar.high.gsea <- gsea(genelist = Ar.geneList.high, geneID2GO = geneID2GO)

# A22 'Low' genes not in Ar**
A22.geneList.low <- rep(0, length = length(Ap.geneList))
names(A22.geneList.low) <- names(Ap.geneList)
A22.geneList.low[(which(names(A22.geneList.low) %in% A22.low))] <- 1
A22.low.gsea <- gsea(genelist = A22.geneList.low, geneID2GO = geneID2GO)
# Ar 'Low' genes not in A22**
Ar.geneList.low <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.low) <- names(Ap.geneList)
Ar.geneList.low[(which(names(Ar.geneList.low) %in% Ar.low))] <- 1
Ar.low.gsea <- gsea(genelist = Ar.geneList.low, geneID2GO = geneID2GO)

# A22 'Bimodal' genes not in Ar**
A22.geneList.bim <- rep(0, length = length(Ap.geneList))
names(A22.geneList.bim) <- names(Ap.geneList)
A22.geneList.bim[(which(names(A22.geneList.bim) %in% A22.bim))] <- 1
A22.bim.gsea <- gsea(genelist = A22.geneList.bim, geneID2GO = geneID2GO)

# Ar 'Bimodal' genes not in A22**
Ar.geneList.bim <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.bim) <- names(Ap.geneList)
Ar.geneList.bim[(which(names(Ar.geneList.bim) %in% Ar.bim))] <- 1
Ar.bim.gsea <- gsea(genelist = Ar.geneList.bim, geneID2GO = geneID2GO)

# A22 'Intermediate' genes not in Ar**
A22.geneList.int <- rep(0, length = length(Ap.geneList))
names(A22.geneList.int) <- names(Ap.geneList)
A22.geneList.int[(which(names(A22.geneList.int) %in% A22.int))] <- 1
table(A22.geneList.int)
A22.int.gsea <- gsea(genelist = A22.geneList.int, geneID2GO = geneID2GO)

# Ar 'Intermediate' genes not in A22**
Ar.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ar.geneList.int) <- names(Ap.geneList)
Ar.geneList.int[(which(names(Ar.geneList.int) %in% Ar.int))] <- 1
Ar.int.gsea <- gsea(genelist = Ar.geneList.int, geneID2GO = geneID2GO)

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
Ap.setdiff.gsea <- rbind(A22.gsea, Ar.gsea)
# reorder
Ap.setdiff.gsea <- Ap.setdiff.gsea[, c("Species", "Type", "GO.ID", "Term", "Annotated", 
    "Significant", "Expected", "parentchild")]
colnames(Ap.setdiff.gsea)[8] <- "P"

write.csv(Ap.setdiff.gsea, file = paste(resultsdir, "Ap_setdiff_GSEA_", Sys.Date(), 
    ".csv", sep = ""), row.names = FALSE)
```



-------------------------------------------------------------------
 Species      Type       GO.ID                  Term               
--------- ------------ ---------- ---------------------------------
  ApVT        High     GO:1902222             erythrose            
                                  4-phosphate/phosphoenolpyruvat...

  ApVT        High     GO:0030435     sporulation resulting in     
                                         formation of a ...        

  ApVT        High     GO:0006559     L-phenylalanine catabolic    
                                               process             

  ApVT        High     GO:0051049      regulation of transport     

  ApVT        High     GO:0043269    regulation of ion transport   

  ApVT        High     GO:0033002     muscle cell proliferation    

  ApVT        High     GO:0006424   glutamyl-tRNA aminoacylation   

  ApVT        High     GO:0043934            sporulation           

  ApVT        High     GO:0015074          DNA integration         

  ApVT        High     GO:0006412            translation           

  ApVT        Low      GO:1902589     single-organism organelle    
                                            organization           

  ApVT        Low      GO:0006996      organelle organization      

  ApVT        Low      GO:0019538     protein metabolic process    

  ApVT        Low      GO:0022406         membrane docking         

  ApVT        Low      GO:0043170      macromolecule metabolic     
                                               process             

  ApVT        Low      GO:0051641       cellular localization      

  ApVT        Low      GO:0048278          vesicle docking         

  ApVT        Low      GO:0046907      intracellular transport     

  ApVT        Low      GO:0000725      recombinational repair      

  ApVT        Low      GO:0022904       respiratory electron       
                                           transport chain         

  ApVT        Low      GO:0051649   establishment of localization  
                                               in cell             

  ApVT        Low      GO:0044260      cellular macromolecule      
                                          metabolic process        

  ApVT        Low      GO:1900542       regulation of purine       
                                       nucleotide metaboli...      

  ApVT        Low      GO:0044707   single-multicellular organism  
                                               process             

  ApVT        Low      GO:0044267    cellular protein metabolic    
                                               process             

  ApVT        Low      GO:0051606       detection of stimulus      

  ApVT        Low      GO:0007264   small GTPase mediated signal   
                                           transductio...          

  ApVT        Low      GO:0009069     serine family amino acid     
                                         metabolic proce...        

  ApVT        Low      GO:0006184       GTP catabolic process      

  ApVT        Low      GO:0033121       regulation of purine       
                                       nucleotide cataboli...      

  ApVT        Low      GO:0030811     regulation of nucleotide     
                                         catabolic proce...        

  ApVT        Low      GO:0046039       GTP metabolic process      

  ApVT        Low      GO:0016192    vesicle-mediated transport    

  ApVT        Low      GO:0019220      regulation of phosphate     
                                         metabolic proces...       

  ApVT        Low      GO:1901069   guanosine-containing compound  
                                            catabolic ...          

  ApVT        Low      GO:0009593       detection of chemical      
                                              stimulus             

  ApVT        Low      GO:0003008          system process          

  ApVT        Low      GO:0015672    monovalent inorganic cation   
                                              transport            

  ApVT        Low      GO:0035315     hair cell differentiation    

  ApVT        Low      GO:0071840        cellular component        
                                      organization or bioge...     

  ApVT        Low      GO:0032501     multicellular organismal     
                                               process             

  ApVT        Low      GO:0032970        regulation of actin       
                                       filament-based proce...     

  ApVT        Low      GO:0042439      ethanolamine-containing     
                                         compound metabol...       

  ApVT        Low      GO:0006904    vesicle docking involved in   
                                             exocytosis            

  ApVT        Low      GO:0042726    flavin-containing compound    
                                          metabolic pro...         

  ApVT    Intermediate GO:0009561         megagametogenesis        

  ApVT    Intermediate GO:0009553      embryo sac development      

  ApVT    Intermediate GO:0009560        embryo sac egg cell       
                                           differentiation         

  ApVT    Intermediate GO:0048229      gametophyte development     

  ApVT    Intermediate GO:0001558     regulation of cell growth    

  ApVT    Intermediate GO:0030031     cell projection assembly     

  ApVT    Intermediate GO:0015942     formate metabolic process    

  ApVT    Intermediate GO:0006259       DNA metabolic process      

  ApVT    Intermediate GO:0042126     nitrate metabolic process    

  ApVT    Intermediate GO:0072528       pyrimidine-containing      
                                        compound biosynthe...      

  ApVT      Bimodal    GO:0006546     glycine catabolic process    

  ApVT      Bimodal    GO:0006468      protein phosphorylation     

  ApVT      Bimodal    GO:0009071     serine family amino acid     
                                         catabolic proce...        

  ApVT      Bimodal    GO:0007215   glutamate receptor signaling   
                                               pathway             

  ApVT      Bimodal    GO:0006886       intracellular protein      
                                              transport            

  AcNC        High     GO:0006259       DNA metabolic process      

  AcNC        High     GO:0001503           ossification           

  AcNC        High     GO:0008217   regulation of blood pressure   

  AcNC        High     GO:0071216    cellular response to biotic   
                                              stimulus             

  AcNC        High     GO:0002790         peptide secretion        

  AcNC        High     GO:0009914         hormone transport        

  AcNC        High     GO:0006457          protein folding         

  AcNC        High     GO:0030099   myeloid cell differentiation   

  AcNC        High     GO:0002831     regulation of response to    
                                          biotic stimulu...        

  AcNC        High     GO:0035967       cellular response to       
                                       topologically incor...      

  AcNC        High     GO:0007249     I-kappaB kinase/NF-kappaB    
                                              signaling            

  AcNC        High     GO:0034976      response to endoplasmic     
                                          reticulum stress         

  AcNC        High     GO:0001501    skeletal system development   

  AcNC        High     GO:0090304      nucleic acid metabolic      
                                               process             

  AcNC        Low      GO:0015074          DNA integration         

  AcNC        Low      GO:0006259       DNA metabolic process      

  AcNC        Low      GO:0051049      regulation of transport     

  AcNC        Low      GO:0010016    shoot system morphogenesis    

  AcNC        Low      GO:0048316         seed development         

  AcNC        Low      GO:0007249     I-kappaB kinase/NF-kappaB    
                                              signaling            

  AcNC        Low      GO:0043122      regulation of I-kappaB      
                                        kinase/NF-kappaB ...       

  AcNC        Low      GO:0008152         metabolic process        

  AcNC        Low      GO:0009164   nucleoside catabolic process   

  AcNC        Low      GO:0072330        monocarboxylic acid       
                                        biosynthetic process       

  AcNC        Low      GO:0048532       anatomical structure       
                                             arrangement           

  AcNC        Low      GO:0043269    regulation of ion transport   

  AcNC        Low      GO:0033002     muscle cell proliferation    

  AcNC        Low      GO:0060548    negative regulation of cell   
                                                death              

  AcNC        Low      GO:0033013      tetrapyrrole metabolic      
                                               process             

  AcNC        Low      GO:0048507       meristem development       

  AcNC        Low      GO:0006424   glutamyl-tRNA aminoacylation   

  AcNC        Low      GO:0048827       phyllome development       

  AcNC        Low      GO:0009793   embryo development ending in   
                                           seed dorman...          

  AcNC        Low      GO:0006633      fatty acid biosynthetic     
                                               process             

  AcNC        Low      GO:1901658    glycosyl compound catabolic   
                                               process             

  AcNC        Low      GO:0006760       folic acid-containing      
                                        compound metabolic...      

  AcNC        Low      GO:0043069      negative regulation of      
                                        programmed cell d...       

  AcNC        Low      GO:0003012       muscle system process      

  AcNC    Intermediate GO:0019538     protein metabolic process    

  AcNC    Intermediate GO:0044267    cellular protein metabolic    
                                               process             

  AcNC    Intermediate GO:0043170      macromolecule metabolic     
                                               process             

  AcNC    Intermediate GO:0006996      organelle organization      

  AcNC    Intermediate GO:0071840        cellular component        
                                      organization or bioge...     

  AcNC    Intermediate GO:1902589     single-organism organelle    
                                            organization           

  AcNC    Intermediate GO:0016192    vesicle-mediated transport    

  AcNC    Intermediate GO:0006184       GTP catabolic process      

  AcNC    Intermediate GO:0006644      phospholipid metabolic      
                                               process             

  AcNC    Intermediate GO:0007264   small GTPase mediated signal   
                                           transductio...          

  AcNC    Intermediate GO:1901069   guanosine-containing compound  
                                            catabolic ...          

  AcNC    Intermediate GO:0051649   establishment of localization  
                                               in cell             

  AcNC    Intermediate GO:0046039       GTP metabolic process      

  AcNC    Intermediate GO:1900542       regulation of purine       
                                       nucleotide metaboli...      

  AcNC    Intermediate GO:0044260      cellular macromolecule      
                                          metabolic process        

  AcNC    Intermediate GO:0061025          membrane fusion         

  AcNC    Intermediate GO:0046907      intracellular transport     

  AcNC    Intermediate GO:0022406         membrane docking         

  AcNC    Intermediate GO:0070085           glycosylation          

  AcNC    Intermediate GO:0051641       cellular localization      

  AcNC    Intermediate GO:1901068   guanosine-containing compound  
                                            metabolic ...          

  AcNC    Intermediate GO:0033121       regulation of purine       
                                       nucleotide cataboli...      

  AcNC    Intermediate GO:0044801     single-organism membrane     
                                               fusion              

  AcNC    Intermediate GO:0010033   response to organic substance  

  AcNC    Intermediate GO:0032970        regulation of actin       
                                       filament-based proce...     

  AcNC    Intermediate GO:0022904       respiratory electron       
                                           transport chain         

  AcNC    Intermediate GO:0030811     regulation of nucleotide     
                                         catabolic proce...        

  AcNC    Intermediate GO:0051345      positive regulation of      
                                        hydrolase activit...       

  AcNC    Intermediate GO:0032956        regulation of actin       
                                       cytoskeleton organiz...     

  AcNC    Intermediate GO:0006486       protein glycosylation      

  AcNC    Intermediate GO:0043412    macromolecule modification    

  AcNC    Intermediate GO:0043547   positive regulation of GTPase  
                                              activity             

  AcNC    Intermediate GO:0006412            translation           

  AcNC    Intermediate GO:0000725      recombinational repair      

  AcNC    Intermediate GO:0006388        tRNA splicing, via        
                                      endonucleolytic cleav...     

  AcNC    Intermediate GO:0009118     regulation of nucleoside     
                                         metabolic proce...        

  AcNC    Intermediate GO:0019220      regulation of phosphate     
                                         metabolic proces...       

  AcNC    Intermediate GO:0031329      regulation of cellular      
                                          catabolic process        

  AcNC    Intermediate GO:0042439      ethanolamine-containing     
                                         compound metabol...       

  AcNC    Intermediate GO:0032271       regulation of protein      
                                           polymerization          

  AcNC    Intermediate GO:0006140     regulation of nucleotide     
                                         metabolic proce...        

  AcNC    Intermediate GO:0048232      male gamete generation      

  AcNC    Intermediate GO:0048278          vesicle docking         

  AcNC    Intermediate GO:0050684   regulation of mRNA processing  

  AcNC    Intermediate GO:0006656        phosphatidylcholine       
                                        biosynthetic process       

  AcNC      Bimodal    GO:0006259       DNA metabolic process      

  AcNC      Bimodal    GO:0003151    outflow tract morphogenesis   

  AcNC      Bimodal    GO:0071804      cellular potassium ion      
                                              transport            

  AcNC      Bimodal    GO:0071805    potassium ion transmembrane   
                                              transport            

  AcNC      Bimodal    GO:0006725    cellular aromatic compound    
                                          metabolic pro...         

  AcNC      Bimodal    GO:0034641    cellular nitrogen compound    
                                          metabolic pro...         

  AcNC      Bimodal    GO:0006139       nucleobase-containing      
                                        compound metabolic...      

  AcNC      Bimodal    GO:0072528       pyrimidine-containing      
                                        compound biosynthe...      

  AcNC      Bimodal    GO:0046483   heterocycle metabolic process  

  AcNC      Bimodal    GO:0051606       detection of stimulus      
-------------------------------------------------------------------

Table: Enriched GO terms found in one species but not the other for each thermal-response category


## Visualize responsive transcripts

Make plots for all genes expressed at *High* temps in GO category "GO:0006950: response to stress"

![plot of chunk plot_GOstress](figure/plot_GOstress1.png) ![plot of chunk plot_GOstress](figure/plot_GOstress2.png) 


Make plots for all genes expressed at *Low* temps in GO category "GO:0006950: response to stress"

![plot of chunk plot_GOstress_low](figure/plot_GOstress_low1.png) ![plot of chunk plot_GOstress_low](figure/plot_GOstress_low2.png) 




```r
trp_A22_high <- ggplot(resp.TPM.dt.sub[A22.high], aes(x = val, y = TPM.scaled, 
    group = Transcript)) + geom_smooth(method = "lm", formula = y ~ poly(x, 
    2)) + facet_grid(. ~ colony2) + scale_y_continuous(name = "Expression (scaled)") + 
    scale_x_continuous(name = expression(paste("Temperature ", degree, "C")))
print(trp_A22_high)
```

![plot of chunk ggplot_high](figure/ggplot_high.png) 


## Shiny interactive web-app

To assist visualization of specific transcripts, I made a interactive web-app using the [shiny](http://www.rstudio.com/shiny/) package. The scripts for this app are in the sub-directory `.\ApRxN-shinyapp`.

Export data for interactive shiny app. 









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
## [1] grid      parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] Rgraphviz_2.8.1      topGO_2.16.0         SparseM_1.03        
##  [4] GO.db_2.14.0         RSQLite_0.11.4       DBI_0.2-7           
##  [7] AnnotationDbi_1.26.0 GenomeInfoDb_1.0.2   Biobase_2.24.0      
## [10] BiocGenerics_0.10.0  graph_1.42.0         xtable_1.7-3        
## [13] MASS_7.3-31          plyr_1.8.1           RCurl_1.95-4.1      
## [16] bitops_1.0-6         data.table_1.9.2     stringr_0.6.2       
## [19] pander_0.3.8         knitcitations_0.5-0  bibtex_0.3-6        
## [22] ggplot2_0.9.3.1      R.utils_1.29.8       R.oo_1.18.0         
## [25] R.methodsS3_1.6.1    knitr_1.5           
## 
## loaded via a namespace (and not attached):
##  [1] codetools_0.2-8    colorspace_1.2-4   dichromat_2.0-0   
##  [4] digest_0.6.4       evaluate_0.5.3     formatR_0.10      
##  [7] gtable_0.1.2       httr_0.3           IRanges_1.22.3    
## [10] labeling_0.2       lattice_0.20-29    munsell_0.4.2     
## [13] proto_0.3-10       RColorBrewer_1.0-5 Rcpp_0.11.1       
## [16] reshape2_1.2.2     scales_0.2.3       stats4_3.1.0      
## [19] tools_3.1.0        XML_3.98-1.1
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

