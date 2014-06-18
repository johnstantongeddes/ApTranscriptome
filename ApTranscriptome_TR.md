Thermal reactionome of the common ant species *Aphaenogaster picea* and *A. carolinensis*
========================================================================================
   
**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

**June 16, 2014**

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
# URL for annotation file
annotation.URL <- getURL("http://johnstantongeddes.org/assets/files/ApTranscriptome_AnnotationTable_20140113.txt")
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
samples <- c("A22-0", "A22-10", "A22-14", "A22-17", "A22-21", "A22-24", "A22-28", "A22-31", "A22-35", "A22-38", "A22-3", "A22-7", "Ar-0", "Ar-10", "Ar-14", "Ar-17", "Ar-21", "Ar-24", "Ar-28", "Ar-31", "Ar-35", "Ar-38", "Ar-3", "Ar-7")

for (j in 1:length(samples)) {
    message("Start expression quantification for sample ", samples[j], ": ", Sys.time())
    quantdir <- paste(sfexpressionroot, samples[j], "_quant", sep="")
    samp.pos <- grep(paste(paste(samples[j], "_", sep="")), paired.left)
    samp.paired.l <- paste(readdir, paired.left[samp.pos], sep="")
    samp.paired.r <- paste(readdir, paired.right[samp.pos], sep="")
    samp.unpaired.l <- paste(readdir, unpaired.left[samp.pos], sep="")
    samp.unpaired.r <- paste(readdir, unpaired.right[samp.pos], sep="")
    sailfishcmd <- paste("sailfish quant -i ", sfindex, " -o ", quantdir, " -l 'T=SE:S=U' -r ", samp.paired.l, " ", samp.paired.r, " ", samp.unpaired.l, " ", samp.unpaired.r, " -p 4", sep="")
    #print(sailfishcmd)
    system(sailfishcmd)
    message("Done with expression quantification for sample ", samples[j], ": ", Sys.time(), "\n")
}
```

```
## Start expression quantification for sample A22-0: 2014-06-05 16:39:19
## Done with expression quantification for sample A22-0: 2014-06-05 16:39:48
## 
## Start expression quantification for sample A22-10: 2014-06-05 16:39:48
## Done with expression quantification for sample A22-10: 2014-06-05 16:40:09
## 
## Start expression quantification for sample A22-14: 2014-06-05 16:40:09
## Done with expression quantification for sample A22-14: 2014-06-05 16:40:28
## 
## Start expression quantification for sample A22-17: 2014-06-05 16:40:28
## Done with expression quantification for sample A22-17: 2014-06-05 16:40:50
## 
## Start expression quantification for sample A22-21: 2014-06-05 16:40:50
## Done with expression quantification for sample A22-21: 2014-06-05 16:41:09
## 
## Start expression quantification for sample A22-24: 2014-06-05 16:41:09
## Done with expression quantification for sample A22-24: 2014-06-05 16:41:30
## 
## Start expression quantification for sample A22-28: 2014-06-05 16:41:30
## Done with expression quantification for sample A22-28: 2014-06-05 16:41:48
## 
## Start expression quantification for sample A22-31: 2014-06-05 16:41:48
## Done with expression quantification for sample A22-31: 2014-06-05 16:42:07
## 
## Start expression quantification for sample A22-35: 2014-06-05 16:42:07
## Done with expression quantification for sample A22-35: 2014-06-05 16:42:27
## 
## Start expression quantification for sample A22-38: 2014-06-05 16:42:27
## Done with expression quantification for sample A22-38: 2014-06-05 16:42:48
## 
## Start expression quantification for sample A22-3: 2014-06-05 16:42:48
## Done with expression quantification for sample A22-3: 2014-06-05 16:43:09
## 
## Start expression quantification for sample A22-7: 2014-06-05 16:43:09
## Done with expression quantification for sample A22-7: 2014-06-05 16:43:28
## 
## Start expression quantification for sample Ar-0: 2014-06-05 16:43:28
## Done with expression quantification for sample Ar-0: 2014-06-05 16:43:49
## 
## Start expression quantification for sample Ar-10: 2014-06-05 16:43:49
## Done with expression quantification for sample Ar-10: 2014-06-05 16:44:10
## 
## Start expression quantification for sample Ar-14: 2014-06-05 16:44:10
## Done with expression quantification for sample Ar-14: 2014-06-05 16:44:30
## 
## Start expression quantification for sample Ar-17: 2014-06-05 16:44:30
## Done with expression quantification for sample Ar-17: 2014-06-05 16:44:51
## 
## Start expression quantification for sample Ar-21: 2014-06-05 16:44:51
## Done with expression quantification for sample Ar-21: 2014-06-05 16:45:12
## 
## Start expression quantification for sample Ar-24: 2014-06-05 16:45:12
## Done with expression quantification for sample Ar-24: 2014-06-05 16:45:31
## 
## Start expression quantification for sample Ar-28: 2014-06-05 16:45:31
## Done with expression quantification for sample Ar-28: 2014-06-05 16:45:51
## 
## Start expression quantification for sample Ar-31: 2014-06-05 16:45:51
## Done with expression quantification for sample Ar-31: 2014-06-05 16:46:11
## 
## Start expression quantification for sample Ar-35: 2014-06-05 16:46:11
## Done with expression quantification for sample Ar-35: 2014-06-05 16:46:32
## 
## Start expression quantification for sample Ar-38: 2014-06-05 16:46:32
## Done with expression quantification for sample Ar-38: 2014-06-05 16:46:55
## 
## Start expression quantification for sample Ar-3: 2014-06-05 16:46:55
## Done with expression quantification for sample Ar-3: 2014-06-05 16:47:14
## 
## Start expression quantification for sample Ar-7: 2014-06-05 16:47:14
## Done with expression quantification for sample Ar-7: 2014-06-05 16:47:35
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
A22.TPM[,colony:="A22"]
Ar.TPM[,colony:="Ar"]
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
##  $ TPM              : num  0.0788 0.0626 0.0355 0.0923 0.0394 ...
##  $ RPKM             : num  0.0918 0.1049 0.0528 0.1411 0.0491 ...
##  $ KPKM             : num  0.0918 0.1049 0.0528 0.1411 0.0491 ...
##  $ EstimatedNumReads: num  2968 1462 1364 2425 1710 ...
##  $ V7               : num  37.1 18.1 17 30 21.3 ...
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

Of the 99861 transcripts, 22306 have models with P < 0.05.

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

At the 5% FDR significance threshold, there are 10579 transcripts with an overall significant model.


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
```


|    Coefficient     |  Number.significant  |
|:------------------:|:--------------------:|
|       Total        |        10579         |
|       Colony       |         1479         |
|    Temperature     |         2239         |
| Temperature:Colony |         6861         |

Table: Number of transcripts of 99,861 total with expression that depends on species, temperature or their interaction at 5% FDR.


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
## 'data.frame':	9100 obs. of  9 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100015|*|comp3543055_c0_seq1" "100067|*|comp3557646_c0_seq1" "100089|*|comp11313_c1_seq1" ...
##  $ A22.max   : num  38.5 0 NA 18.5 NA 0 NA 0 0 38.5 ...
##  $ A22.min   : num  0 20.5 NA 38.5 NA 23 NA 28.5 26 10.5 ...
##  $ A22.opt   : num  1.025 0.945 1 1.057 1 ...
##  $ A22.type  : chr  "High" "Bimodal" "NotResp" "Intermediate" ...
##  $ Ar.max    : num  0 38.5 0 38.5 38.5 18 0 0 NA NA ...
##  $ Ar.min    : num  38.5 0 20.5 13.5 18.5 38.5 38.5 38.5 NA NA ...
##  $ Ar.opt    : num  1.199 1.004 0.934 0.962 0.897 ...
##  $ Ar.type   : chr  "Low" "High" "Bimodal" "High" ...
```

```r
# save results to file
write.table(file = paste(resultsdir, "Ap-responsive-transcripts", Sys.Date(), ".txt", sep = ""), Ap.response.type, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
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
## A22.type.table    1521 1269          918 4886     506
## Ar.type.table      866 1098         2681 3699     756
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
## X-squared = 1269, df = 4, p-value < 2.2e-16
```

```r
# Reorganize for plotting to show overlap among categories between colonies
type.table <- table(Ap.response.type[ , 'Ar.type'], Ap.response.type[ , 'A22.type'])
# reorder 
tt2 <- type.table[c(2,4,1,3,5), c(2,4,1,3,5)]

# plot
mosaic.colors <- brewer.pal(5, "Blues")
par(mar = c(5.1, 3.1, 2.1, 7.5))
mosaicplot(tt2, ylab = expression(italic("A. picea")), xlab = "", main = "", color = mosaic.colors, cex.axis = 0.01)
labels <- c("High", "Low", "Bimodal", "Intermediate", "NotResp")
text(c(.1,.3,.55,.75,.95), par("usr")[3] - 0.05, srt = 45, adj = 1, labels = labels, xpd = TRUE, font = 1)
text(.5, par("usr")[4], "A. carolinensis", xpd = TRUE, font = 3)
legend(1.05, .75, fill = mosaic.colors, c("High", "Low", "Bimodal", "Intermediate", "NotResp"), xpd=TRUE)
```

![plot of chunk temperature_response](figure/temperature_response.png) 

The number of thermally-responsive in each response category differs between the colonies, with the msot transcripts expressed at *Low* temperatures in both colonies. For *ApVT*, an equal number of transcripts are expressed at *High* and *Bimodal*, followed by *Intermediate* transcripts. For *AcNC*, transcripts expressed at *Intermediate* temperatures are next most common, followed by *High* and *Bimodal*. 

Interestingly, nearly half of the *High* genes in *AcNC* are *Low* in *ApVT*, and vice versa. In contrast, most of the *Low* genes in one species are also *Low* in the other species.



|   &nbsp;   |     ApVT     |  &nbsp;  |  &nbsp;  |  &nbsp;  |    &nbsp;    |  &nbsp;  |  &nbsp;  |
|:----------:|:------------:|:--------:|:--------:|:--------:|:------------:|:--------:|:--------:|
|  **AcNC**  |              |   High   |   Low    | Bimodal  | Intermediate | NotResp  |  Total   |
|            |     High     |   322    |   467    |   153    |      65      |    91    |   1098   |
|            |     Low      |   405    |   2594   |   297    |     163      |   240    |   3699   |
|            |   Bimodal    |   166    |   181    |   297    |     114      |   108    |   866    |
|            | Intermediate |   155    |   1350   |   551    |     558      |    67    |   2681   |
|            |   NotResp    |   221    |   294    |   223    |      18      |    0     |   756    |
|            |    Total     |   1269   |   4886   |   1521   |     918      |   506    |   9100   |

Table: Number of transcripts with maximum expression at high, low, intermediate, both high and low (bimodal) temperatures or are not thermally-responsivefor each species and their overlap.

Table 4 shows the number of transcripts that fall into each expression type for each each species. The totals for each species include the 2239 transcripts that have consistent temperature responses between the two colonies. 



## Species comparison

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
  #scale_fill_manual(name = "Colony", values = c("white", "black")) + 
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Temperature of maximum expression")
print(maxexplot)
```

![plot of chunk max_exp_PDF](figure/max_exp_PDF.png) 


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
## t = -0.149, df = 2400, p-value = 0.8814
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -1126   967
## sample estimates:
## mean of x mean of y 
##       375       455
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
## t = -3.65, df = 2523, p-value = 0.0002726
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2633 -0.0791
## sample estimates:
## mean of x mean of y 
##      1.11      1.28
```

```r
boxplot(log(A22.high.transcripts$A22.opt+1), log(A22.high.transcripts$Ar.opt+1))
```

![plot of chunk optimum_expression_comparison](figure/optimum_expression_comparison2.png) 

The `t.test` fails to account for the many orders of magnitude difference in expression among transcripts, e.g. non-equal variances. This problem is the key issue in the analysis of differential expression (<a href="http://dx.doi.org/10.1186/1471-2105-11-94">Bullard et al. 2010</a>; <a href="http://dx.doi.org/10.1038/nprot.2013.099">Anders et al. 2013</a>). As my goal is simply to determine if expression is typically greater at optimal temperatures (19.5 C) in *Ar* than *A22* for genes that are up-regulated at high temperatures in *A22*, I use a non-parametric Wilcoxon signed rank-test


```r
w1 <- wilcox.test(A22.high.transcripts$A22.opt, A22.high.transcripts$Ar.opt, alternative = "two.sided", paired = TRUE, conf.int = TRUE)
w1
```

```
## 
## 	Wilcoxon signed rank test with continuity correction
## 
## data:  A22.high.transcripts$A22.opt and A22.high.transcripts$Ar.opt
## V = 210438, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.293 -0.172
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
## t = -6.51, df = 7387, p-value = 8.071e-11
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1710 -0.0918
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
## V = 2141816, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.451 -0.357
## sample estimates:
## (pseudo)median 
##         -0.403
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
## t = -6.35, df = 6074, p-value = 2.286e-10
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.239 -0.126
## sample estimates:
## mean of x mean of y 
##      1.34      1.53
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
## V = 1489602, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
## 95 percent confidence interval:
##  -0.501 -0.379
## sample estimates:
## (pseudo)median 
##         -0.438
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
## [1] 918
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
(t.thermbreadth <- t.test(Ar.int.sd$exp_sd, A22.int.sd$exp_sd, alternative = "two.sided"))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.int.sd$exp_sd and A22.int.sd$exp_sd
## t = 9.15, df = 1314, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.336 0.519
## sample estimates:
## mean of x mean of y 
##     10.35      9.93
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
## [1] 1521
```

```r
# apply `transcriptSD` function to all transcripts
A22.bim.sd <- unlist(Map(transcriptSD, A22.bim.lm, colony = "A22"))

# repeat for Ar 
Ar.bim.lm <- RxNlmAIC[Ap.response.type[which(Ap.response.type$Ar.type == "Bimodal"), "Transcript"]]
length(Ar.bim.lm)
```

```
## [1] 866
```

```r
Ar.bim.sd <- unlist(Map(transcriptSD, Ar.bim.lm, colony = "Ar"))

# t-test to compare standard deviation of 'Bimodal' transcripts between colonies
(t.themsens <- t.test(Ar.bim.sd, A22.bim.sd))
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  Ar.bim.sd and A22.bim.sd
## t = -0.84, df = 1538, p-value = 0.4012
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.2103  0.0842
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
## Classes 'data.table' and 'data.frame':	200200 obs. of  10 variables:
##  $ Transcript       : chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ Length           : int  208 208 208 208 208 208 208 208 208 208 ...
##  $ TPM              : num  0 0.3555 0 0.7348 0.0743 ...
##  $ RPKM             : num  0 0.5963 0 1.1238 0.0925 ...
##  $ KPKM             : num  0 0.5963 0 1.1238 0.0925 ...
##  $ EstimatedNumReads: num  0 157.5 0 366.1 61.1 ...
##  $ V7               : num  0 1.95 0 4.534 0.762 ...
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
## [1] 9100
```

```r
# scale transcripts so can compare
resp.TPM.dt.sub[,TPM.scaled:=scale(TPM), by = Transcript]
```

```
##                          Transcript Length    TPM   RPKM   KPKM
##      1: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      2: 100008|*|comp137625_c0_seq2    208 0.3555 0.5963 0.5963
##      3: 100008|*|comp137625_c0_seq2    208 0.0000 0.0000 0.0000
##      4: 100008|*|comp137625_c0_seq2    208 0.7348 1.1238 1.1238
##      5: 100008|*|comp137625_c0_seq2    208 0.0743 0.0925 0.0925
##     ---                                                        
## 200196:      9|*|comp147140_c0_seq1   9030 0.7246 1.2147 1.2147
## 200197:      9|*|comp147140_c0_seq1   9030 0.7184 0.9889 0.9889
## 200198:      9|*|comp147140_c0_seq1   9030 0.4516 0.8023 0.8023
## 200199:      9|*|comp147140_c0_seq1   9030 0.5661 0.7738 0.7738
## 200200:      9|*|comp147140_c0_seq1   9030 0.4600 0.8408 0.8408
##         EstimatedNumReads      V7 sample  val colony TPM.scaled
##      1:               0.0   0.000  A22-0  0.0    A22     -0.566
##      2:             157.5   1.950   Ar-0  0.0     Ar      1.124
##      3:               0.0   0.000  A22-3  3.5    A22     -0.566
##      4:             366.1   4.534   Ar-3  3.5     Ar      2.927
##      5:              61.1   0.762 A22-10 10.5    A22     -0.213
##     ---                                                        
## 200196:           17872.1 220.982  Ar-31 31.5     Ar     -0.260
## 200197:           23722.3 296.358 A22-35 35.0    A22     -0.280
## 200198:            8964.3 111.016  Ar-35 35.0     Ar     -1.140
## 200199:           23785.7 297.030 A22-38 38.5    A22     -0.771
## 200200:           12649.4 156.296  Ar-38 38.5     Ar     -1.113
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


Calculate T~on~ for *High* genes in each species.


```r
# transcripts expressed at *High* and *Bimodal* temperatures in A22
A22.high.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(union(A22.high.transcripts, A22.bim.transcripts), "A22")]
setkey(A22.high.TPM.dt.sub, Transcript)
str(A22.high.TPM.dt.sub)
```

```
## Classes 'data.table' and 'data.frame':	30690 obs. of  5 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0 0 0.0743 0 0 ...
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
## Classes 'data.table' and 'data.frame':	21604 obs. of  5 variables:
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
## t = 5.95, df = 4138, p-value = 2.945e-09
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.591 1.171
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
## Classes 'data.table' and 'data.frame':	70477 obs. of  5 variables:
##  $ Transcript: chr  "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" "100015|*|comp3543055_c0_seq1" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0.749 0 0 0 0 ...
##  $ val       : num  0 3.5 10.5 14 17.5 21 24.5 28 31.5 35 ...
##  $ pTPM      : num  1.45 1.27 1.048 0.988 0.954 ...
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
## Classes 'data.table' and 'data.frame':	50215 obs. of  5 variables:
##  $ Transcript: chr  "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" "100008|*|comp137625_c0_seq2" ...
##  $ colony    : Factor w/ 2 levels "A22","Ar": 1 1 1 1 1 1 1 1 1 1 ...
##  $ TPM       : num  0 0 0.0743 0 0 ...
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
## t = -2.17, df = 9939, p-value = 0.02966
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.3263 -0.0169
## sample estimates:
## mean of x mean of y 
##      12.8      12.9
```

Genes with increased expression at *Low* temperatures are on average turned on 0.2°C higher in *ApVT* than *AcNC*.

Visualize T~on~ for both *Low* and *High* genes on the same plot.

![plot of chunk plot_T_on](figure/plot_T_on.png) 

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
## t = -10.7, df = 6860, p-value < 2.2e-16
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.1064 -0.0733
## sample estimates:
## mean of the differences 
##                 -0.0899
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
## t = 7.8, df = 6860, p-value = 7.193e-15
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  0.0215 0.0360
## sample estimates:
## mean of the differences 
##                  0.0288
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
## t = 0.166, df = 6860, p-value = 0.8679
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.00226  0.00268
## sample estimates:
## mean of the differences 
##                 0.00021
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
## t = -0.943, df = 6860, p-value = 0.3457
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -0.000809  0.000284
## sample estimates:
## mean of the differences 
##               -0.000263
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
## 0.450 0.816 0.985
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
## 0.0148 0.1843 0.5505
```

From this analysis, about 3/4 of the differences in reaction norms between species are due to changes in the mean, with the remainder being due to changes in slope.



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
## Classes 'data.table' and 'data.frame':	9100 obs. of  17 variables:
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
##  $ pval                 : num  0.003741 0.000911 0.0021 0.001547 0.001616 ...
##  $ adj.r.squared        : num  0.52 0.604 0.557 0.575 0.572 ...
##  $ padj                 : num  0.0388 0.0135 0.0254 0.0202 0.0209 ...
##  $ A22.type             : chr  "High" "Bimodal" "NotResp" "Intermediate" ...
##  $ Ar.type              : chr  "Low" "High" "Bimodal" "High" ...
##  - attr(*, "sorted")= chr "Transcript"
##  - attr(*, ".internal.selfref")=<externalptr>
```

```r
write.csv(responsive.lms.ann.type, file = paste(resultsdir, "Ap_responsive_genes_", Sys.Date(), ".csv", sep = ""), row.names = FALSE)
```


### Proportion of responsive transcripts annotated

Of the responsive transcripts, 49.3% are annotated, while 48.7% of all transcripts are annotated.


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
GO0006950.responsive[union(with(GO0006950.responsive, grep("Bimodal", Ar.type)), with(GO0006950.responsive, grep("Bimodal", A22.type))), ]
```

```
##                    Transcript
## 1:    14|*|comp150262_c0_seq1
## 2: 21598|*|comp142101_c0_seq1
## 3: 32312|*|comp933733_c0_seq1
## 4: 58246|*|comp109744_c0_seq1
## 5: 15115|*|comp132715_c0_seq1
## 6: 21384|*|comp149042_c0_seq3
##                                                                                   best.hit.to.nr
## 1: gi|332019420|gb|EGI59904.1| Putative fat-like cadherin-related tumor suppressor-like protein 
## 2:                              gi|332018201|gb|EGI58806.1| Protein lethal(2)essential for life 
## 3:                          gi|367054010|ref|XP_003657383.1| hypothetical protein THITE_2156506 
## 4:                                    gi|493322437|ref|WP_006279741.1| molecular chaperone DnaK 
## 5:                                            gi|194716766|gb|ACF93232.1| heat shock protein 90 
## 6:                         gi|396467618|ref|XP_003837992.1| hypothetical protein LEMA_P120390.1 
##    A22.type Ar.type
## 1:      Low Bimodal
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
## X-squared = 1.7, df = 1, p-value = 0.1919
```

While 37 "response to stress" genes are in the responsive set, this is not enriched compared to the whole transcriptome. 

There are also 7 heat shock related genes in the responsive transcripts, out of 130 total.


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
## [1] 130
```

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

For many of the responsive categories, it appears that there is considerable redundancy in the enriched GO terms. I use the [GOSemSim](http://www.bioconductor.org/packages/release/bioc/html/GOSemSim.html) package to determine semantic similarity of enriched GO terms (<a href="http://dx.doi.org/10.1093/bioinformatics/btq064">Yu et al. 2010</a>) and then perform hierarhichal clustering based on this information distance.


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
Ap.BP.ResTable <- GenTable(Ap.BP.GOdata, parentchild = Ap.BP.resultParentChild, topNodes = 10)
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
# report items in each cluster at height 0.6 
GSEAReportClusters(Ap.high.hclust, h = 0.6)
```

```
## $`Cluster 1`
## regulation of biosynthetic process    regulation of metabolic process 
##                                  1                                  5 
## 
## $`Cluster 2`
## regulation of cellular biosynthetic proc... 
##                                           2 
##    regulation of cellular metabolic process 
##                                           6 
## 
## $`Cluster 3`
## cytokinesis 
##           3 
## 
## $`Cluster 4`
## hindbrain development 
##                     4 
## 
## $`Cluster 5`
## regulation of macromolecule biosynthetic... 
##                                           7
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
# report items in each cluster at height 2
GSEAReportClusters(Ap.low.hclust, h = 2)
```

```
## $`Cluster 1`
##             protein metabolic process 
##                                     1 
##            macromolecule modification 
##                                     3 
##    cellular protein metabolic process 
##                                     5 
##                 protein glycosylation 
##                                     9 
##     glycoprotein biosynthetic process 
##                                    10 
##        glycoprotein metabolic process 
##                                    12 
##              translational initiation 
##                                    13 
##          protein modification process 
##                                    16 
##       macromolecule metabolic process 
##                                    32 
## cellular protein modification process 
##                                    38 
##                   histone acetylation 
##                                    59 
##                  histone modification 
##                                    60 
## 
## $`Cluster 2`
##                               glycosylation 
##                                           2 
##            membrane lipid metabolic process 
##                                          11 
## chondroitin sulfate proteoglycan biosynt... 
##                                          21 
##    chondroitin sulfate biosynthetic process 
##                                          22 
##                glycolipid metabolic process 
##                                          24 
##              sphingolipid metabolic process 
##                                          36 
##                   peptide metabolic process 
##                                          58 
## chondroitin sulfate proteoglycan metabol... 
##                                          62 
##       chondroitin sulfate metabolic process 
##                                          63 
##               glutathione metabolic process 
##                                          64 
##             folic acid biosynthetic process 
##                                          66 
##               photosynthesis, dark reaction 
##                                          68 
## 
## $`Cluster 3`
## small GTPase mediated signal transductio... 
##                                           4 
## posttranscriptional regulation of gene e... 
##                                          18 
##                   regulation of translation 
##                                          26 
##             regulation of catabolic process 
##                                          35 
##           regulation of signal transduction 
##                                          41 
##              regulation of cellular process 
##                                          49 
##      regulation of translational initiation 
##                                          50 
##    regulation of cellular catabolic process 
##                                          51 
## regulation of phosphate metabolic proces... 
##                                          52 
## 
## $`Cluster 4`
## nucleoside triphosphate metabolic proces... 
##                                           6 
##            ribonucleotide metabolic process 
##                                          14 
##           regulation of Ras GTPase activity 
##                                          20 
## regulation of nucleotide metabolic proce... 
##                                          25 
##         purine nucleotide metabolic process 
##                                          27 
## regulation of nucleoside metabolic proce... 
##                                          29 
##                nucleotide catabolic process 
##                                          31 
##                       GTP metabolic process 
##                                          33 
##      nucleoside phosphate catabolic process 
##                                          34 
## guanosine-containing compound metabolic ... 
##                                          40 
## regulation of nucleotide catabolic proce... 
##                                          43 
## regulation of purine nucleotide cataboli... 
##                                          45 
##                       GTP catabolic process 
##                                          46 
## guanosine-containing compound catabolic ... 
##                                          54 
## regulation of purine nucleotide metaboli... 
##                                          61 
## 
## $`Cluster 5`
##                  macromolecule localization 
##                                           7 
##       establishment of protein localization 
##                                           8 
##              detection of external stimulus 
##                                          17 
##                  protein complex biogenesis 
##                                          19 
##                     protein oligomerization 
##                                          23 
##                     intracellular transport 
##                                          30 
##                           protein transport 
##                                          42 
##               detection of abiotic stimulus 
##                                          44 
## positive regulation of protein polymeriz... 
##                                          53 
##             macromolecular complex assembly 
##                                          56 
##       establishment of localization in cell 
##                                          69 
## 
## $`Cluster 6`
##                     muscle cell development 
##                                          15 
##                    renal tubule development 
##                                          28 
##                          tissue homeostasis 
##                                          37 
##                 muscle cell differentiation 
##                                          39 
##        multicellular organismal homeostasis 
##                                          47 
##               urogenital system development 
##                                          48 
## sensory perception of mechanical stimulu... 
##                                          55 
##               Malpighian tubule development 
##                                          57 
##     retina morphogenesis in camera-type eye 
##                                          65 
##                     embryonic morphogenesis 
##                                          67
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
## 99564   297
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
## 		 the algorithm is scoring 656 nontrivial nodes
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
## 	 Level 12:	7 nodes to be scored.
## 
## 	 Level 11:	18 nodes to be scored.
## 
## 	 Level 10:	47 nodes to be scored.
## 
## 	 Level 9:	70 nodes to be scored.
## 
## 	 Level 8:	78 nodes to be scored.
## 
## 	 Level 7:	87 nodes to be scored.
## 
## 	 Level 6:	100 nodes to be scored.
## 
## 	 Level 5:	112 nodes to be scored.
## 
## 	 Level 4:	79 nodes to be scored.
## 
## 	 Level 3:	33 nodes to be scored.
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
# report items in each cluster at height 1.5
GSEAReportClusters(Ap.bim.hclust, h = 1.5)
```

```
## $`Cluster 1`
##                       biological regulation 
##                                           1 
##              regulation of cellular process 
##                                           2 
##            regulation of biological process 
##                                           3 
##     regulation of primary metabolic process 
##                                           5 
##          regulation of biosynthetic process 
##                                           6 
##             regulation of metabolic process 
##                                           7 
##               regulation of gene expression 
##                                           9 
## regulation of nitrogen compound metaboli... 
##                                          12 
## regulation of macromolecule metabolic pr... 
##                                          15 
## 
## $`Cluster 2`
##    regulation of cellular metabolic process 
##                                           4 
## regulation of cellular biosynthetic proc... 
##                                           8 
## regulation of macromolecule biosynthetic... 
##                                          11 
## regulation of cellular macromolecule bio... 
##                                          13 
## regulation of nucleobase-containing comp... 
##                                          21 
## 
## $`Cluster 3`
##                transcription, DNA-templated 
##                                          10 
##                    RNA biosynthetic process 
##                                          18 
## nucleobase-containing compound biosynthe... 
##                                          19 
## purine-containing compound catabolic pro... 
##                                          22 
## 
## $`Cluster 4`
## single organism signaling      response to stimulus 
##                        14                        16 
##                 signaling        cell communication 
##                        20                        23 
## 
## $`Cluster 5`
## protein modification process 
##                           17
```


```r
## Genes with *Intermediate* expression in both colonies
Ap.int <- Ap.response.type[which(Ap.response.type$A22.type == "Intermediate" & Ap.response.type$Ar.type == "Intermediate"), "Transcript"]
# create gene list, setting value to 1 for "int" transcripts
Ap.geneList.int <- rep(0, length = length(Ap.geneList))
names(Ap.geneList.int) <- names(Ap.geneList)
Ap.geneList.int[(which(names(Ap.geneList.int) %in% Ap.int))] <- 1
# check correct number of values set to 1
table(Ap.geneList.int)
# Run GSEA
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
# report items in each cluster at height 1.5
GSEAReportClusters(Ap.int.hclust, h = 1.5)
```

```
## $`Cluster 1`
##        glycolipid metabolic process membrane lipid biosynthetic process 
##                                   1                                   2 
##     glycolipid biosynthetic process    membrane lipid metabolic process 
##                                   6                                   7 
##         macromolecule glycosylation      sphingolipid metabolic process 
##                                  11                                  12 
##                       glycosylation 
##                                  14 
## 
## $`Cluster 2`
##     glycoprotein metabolic process  glycoprotein biosynthetic process 
##                                  3                                  4 
##              protein glycosylation            DNA conformation change 
##                                  5                                  8 
##                      DNA packaging cellular protein metabolic process 
##                                  9                                 15 
##   lipoprotein biosynthetic process 
##                                 16 
## 
## $`Cluster 3`
##                  nucleosome organization 
##                                       10 
## protein-DNA complex subunit organization 
##                                       13
```

## merge results into single table
Ap.high.gsea$Type <- "High"
Ap.low.gsea$Type <- "Low"
Ap.bim.gsea$Type <- "Bimodal"
Ap.int.gsea$Type <- "Intermediate"
# combine
Ap.gsea <- rbind(Ap.high.gsea, Ap.low.gsea, Ap.bim.gsea, Ap.int.gsea)
# reorder
Ap.gsea1 <- Ap.gsea[,c("Type", "GO.ID", "Term", "Annotated", "Significant", "Expected", "parentchild")]
colnames(Ap.gsea1)[7] <- "P"
# write to file
write.csv(Ap.setdiff.gsea, file = paste(resultsdir, "Ap_setdiff_GSEA_", Sys.Date(), ".csv", sep = ""), row.names = FALSE)
```

% latex table generated in R 3.1.0 by xtable 1.7-3 package
% Thu Jun  5 17:36:22 2014
\begin{table}[ht]
\centering
\begin{tabular}{rlll}
  \hline
 & Type & GO.ID & Term \\ 
  \hline
1 & High & GO:0009889 & regulation of biosynthetic process \\ 
  2 & High & GO:0031326 & regulation of cellular biosynthetic proc... \\ 
  3 & High & GO:0000910 & cytokinesis \\ 
  4 & High & GO:0030902 & hindbrain development \\ 
  5 & High & GO:0019222 & regulation of metabolic process \\ 
  6 & High & GO:0031323 & regulation of cellular metabolic process \\ 
  7 & High & GO:0010556 & regulation of macromolecule biosynthetic... \\ 
  8 & Low & GO:0019538 & protein metabolic process \\ 
  9 & Low & GO:0070085 & glycosylation \\ 
  10 & Low & GO:0043412 & macromolecule modification \\ 
  11 & Low & GO:0007264 & small GTPase mediated signal transductio... \\ 
  12 & Low & GO:0044267 & cellular protein metabolic process \\ 
  13 & Low & GO:0009141 & nucleoside triphosphate metabolic proces... \\ 
  14 & Low & GO:0033036 & macromolecule localization \\ 
  15 & Low & GO:0045184 & establishment of protein localization \\ 
  16 & Low & GO:0006486 & protein glycosylation \\ 
  17 & Low & GO:0009101 & glycoprotein biosynthetic process \\ 
  18 & Low & GO:0006643 & membrane lipid metabolic process \\ 
  19 & Low & GO:0009100 & glycoprotein metabolic process \\ 
  20 & Low & GO:0006413 & translational initiation \\ 
  21 & Low & GO:0009259 & ribonucleotide metabolic process \\ 
  22 & Low & GO:0055001 & muscle cell development \\ 
  23 & Low & GO:0036211 & protein modification process \\ 
  24 & Low & GO:0009581 & detection of external stimulus \\ 
  25 & Low & GO:0010608 & posttranscriptional regulation of gene e... \\ 
  26 & Low & GO:0070271 & protein complex biogenesis \\ 
  27 & Low & GO:0032318 & regulation of Ras GTPase activity \\ 
  28 & Low & GO:0050650 & chondroitin sulfate proteoglycan biosynt... \\ 
  29 & Low & GO:0030206 & chondroitin sulfate biosynthetic process \\ 
  30 & Low & GO:0051259 & protein oligomerization \\ 
  31 & Low & GO:0006664 & glycolipid metabolic process \\ 
  32 & Low & GO:0006140 & regulation of nucleotide metabolic proce... \\ 
  33 & Low & GO:0006417 & regulation of translation \\ 
  34 & Low & GO:0006163 & purine nucleotide metabolic process \\ 
  35 & Low & GO:0061326 & renal tubule development \\ 
  36 & Low & GO:0009118 & regulation of nucleoside metabolic proce... \\ 
  37 & Low & GO:0046907 & intracellular transport \\ 
  38 & Low & GO:0009166 & nucleotide catabolic process \\ 
  39 & Low & GO:0043170 & macromolecule metabolic process \\ 
  40 & Low & GO:0046039 & GTP metabolic process \\ 
  41 & Low & GO:1901292 & nucleoside phosphate catabolic process \\ 
  42 & Low & GO:0009894 & regulation of catabolic process \\ 
  43 & Low & GO:0006665 & sphingolipid metabolic process \\ 
  44 & Low & GO:0001894 & tissue homeostasis \\ 
  45 & Low & GO:0006464 & cellular protein modification process \\ 
  46 & Low & GO:0042692 & muscle cell differentiation \\ 
  47 & Low & GO:1901068 & guanosine-containing compound metabolic ... \\ 
  48 & Low & GO:0009966 & regulation of signal transduction \\ 
  49 & Low & GO:0015031 & protein transport \\ 
  50 & Low & GO:0030811 & regulation of nucleotide catabolic proce... \\ 
  51 & Low & GO:0009582 & detection of abiotic stimulus \\ 
  52 & Low & GO:0033121 & regulation of purine nucleotide cataboli... \\ 
  53 & Low & GO:0006184 & GTP catabolic process \\ 
  54 & Low & GO:0048871 & multicellular organismal homeostasis \\ 
  55 & Low & GO:0001655 & urogenital system development \\ 
  56 & Low & GO:0050794 & regulation of cellular process \\ 
  57 & Low & GO:0006446 & regulation of translational initiation \\ 
  58 & Low & GO:0031329 & regulation of cellular catabolic process \\ 
  59 & Low & GO:0019220 & regulation of phosphate metabolic proces... \\ 
  60 & Low & GO:0032273 & positive regulation of protein polymeriz... \\ 
  61 & Low & GO:1901069 & guanosine-containing compound catabolic ... \\ 
  62 & Low & GO:0050954 & sensory perception of mechanical stimulu... \\ 
  63 & Low & GO:0065003 & macromolecular complex assembly \\ 
  64 & Low & GO:0072002 & Malpighian tubule development \\ 
  65 & Low & GO:0006518 & peptide metabolic process \\ 
  66 & Low & GO:0016573 & histone acetylation \\ 
  67 & Low & GO:0016570 & histone modification \\ 
  68 & Low & GO:1900542 & regulation of purine nucleotide metaboli... \\ 
  69 & Low & GO:0050654 & chondroitin sulfate proteoglycan metabol... \\ 
  70 & Low & GO:0030204 & chondroitin sulfate metabolic process \\ 
  71 & Low & GO:0006749 & glutathione metabolic process \\ 
  72 & Low & GO:0060042 & retina morphogenesis in camera-type eye \\ 
  73 & Low & GO:0046656 & folic acid biosynthetic process \\ 
  74 & Low & GO:0048598 & embryonic morphogenesis \\ 
  75 & Low & GO:0019685 & photosynthesis, dark reaction \\ 
  76 & Low & GO:0051649 & establishment of localization in cell \\ 
  77 & Bimodal & GO:0065007 & biological regulation \\ 
  78 & Bimodal & GO:0050794 & regulation of cellular process \\ 
  79 & Bimodal & GO:0050789 & regulation of biological process \\ 
  80 & Bimodal & GO:0031323 & regulation of cellular metabolic process \\ 
  81 & Bimodal & GO:0080090 & regulation of primary metabolic process \\ 
  82 & Bimodal & GO:0009889 & regulation of biosynthetic process \\ 
  83 & Bimodal & GO:0019222 & regulation of metabolic process \\ 
  84 & Bimodal & GO:0031326 & regulation of cellular biosynthetic proc... \\ 
  85 & Bimodal & GO:0010468 & regulation of gene expression \\ 
  86 & Bimodal & GO:0006351 & transcription, DNA-templated \\ 
  87 & Bimodal & GO:0010556 & regulation of macromolecule biosynthetic... \\ 
  88 & Bimodal & GO:0051171 & regulation of nitrogen compound metaboli... \\ 
  89 & Bimodal & GO:2000112 & regulation of cellular macromolecule bio... \\ 
  90 & Bimodal & GO:0044700 & single organism signaling \\ 
  91 & Bimodal & GO:0060255 & regulation of macromolecule metabolic pr... \\ 
  92 & Bimodal & GO:0050896 & response to stimulus \\ 
  93 & Bimodal & GO:0036211 & protein modification process \\ 
  94 & Bimodal & GO:0032774 & RNA biosynthetic process \\ 
  95 & Bimodal & GO:0034654 & nucleobase-containing compound biosynthe... \\ 
  96 & Bimodal & GO:0023052 & signaling \\ 
  97 & Bimodal & GO:0019219 & regulation of nucleobase-containing comp... \\ 
  98 & Bimodal & GO:0072523 & purine-containing compound catabolic pro... \\ 
  99 & Bimodal & GO:0007154 & cell communication \\ 
  100 & Intermediate & GO:0006664 & glycolipid metabolic process \\ 
  101 & Intermediate & GO:0046467 & membrane lipid biosynthetic process \\ 
  102 & Intermediate & GO:0009100 & glycoprotein metabolic process \\ 
  103 & Intermediate & GO:0009101 & glycoprotein biosynthetic process \\ 
  104 & Intermediate & GO:0006486 & protein glycosylation \\ 
  105 & Intermediate & GO:0009247 & glycolipid biosynthetic process \\ 
  106 & Intermediate & GO:0006643 & membrane lipid metabolic process \\ 
  107 & Intermediate & GO:0071103 & DNA conformation change \\ 
  108 & Intermediate & GO:0006323 & DNA packaging \\ 
  109 & Intermediate & GO:0034728 & nucleosome organization \\ 
  110 & Intermediate & GO:0043413 & macromolecule glycosylation \\ 
  111 & Intermediate & GO:0006665 & sphingolipid metabolic process \\ 
  112 & Intermediate & GO:0071824 & protein-DNA complex subunit organization \\ 
  113 & Intermediate & GO:0070085 & glycosylation \\ 
  114 & Intermediate & GO:0044267 & cellular protein metabolic process \\ 
  115 & Intermediate & GO:0042158 & lipoprotein biosynthetic process \\ 
   \hline
\end{tabular}
\end{table}


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
GSEAReportClusters(Ar.int.A22.bim.gsea.term.sim.hclust, h = 1)
```

```
## $`Cluster 1`
## phospholipid catabolic process           histone modification 
##                              1                              5 
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
## glutamate receptor signaling pathway 
##                                    4
```



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
# report items in each cluster at height 1.8
GSEAReportClusters(Ar.int.A22.low.gsea.term.sim.hclust, h = 1.6)
```

```
## $`Cluster 1`
##          cellular protein metabolic process 
##                                           1 
##                   protein metabolic process 
##                                           2 
##                               glycosylation 
##                                          16 
## posttranscriptional regulation of gene e... 
##                                          22 
##                                 translation 
##                                          25 
##                       protein glycosylation 
##                                          30 
##             macromolecule metabolic process 
##                                          32 
## 
## $`Cluster 2`
##      single-organism organelle organization 
##                                           3 
##                      organelle organization 
##                                           5 
## cellular component organization or bioge... 
##                                          13 
##                             membrane fusion 
##                                          18 
##             single-organism membrane fusion 
##                                          26 
## 
## $`Cluster 3`
## membrane docking  vesicle docking 
##                4               12 
## 
## $`Cluster 4`
##        respiratory electron transport chain 
##                                           6 
##                      recombinational repair 
##                                           7 
##              phospholipid metabolic process 
##                                           9 
## ethanolamine-containing compound metabol... 
##                                          14 
##           phospholipid biosynthetic process 
##                                          24 
##    phosphatidylcholine biosynthetic process 
##                                          27 
##              sphingolipid catabolic process 
##                                          33 
## 
## $`Cluster 5`
##            vesicle-mediated transport 
##                                     8 
##                 cellular localization 
##                                    10 
## establishment of localization in cell 
##                                    11 
##               intracellular transport 
##                                    31 
##                  endoderm development 
##                                    34 
## 
## $`Cluster 6`
##                       GTP catabolic process 
##                                          15 
##                       GTP metabolic process 
##                                          19 
## guanosine-containing compound catabolic ... 
##                                          21 
## guanosine-containing compound metabolic ... 
##                                          28 
## 
## $`Cluster 7`
## positive regulation of biological proces... 
##                                          17 
## small GTPase mediated signal transductio... 
##                                          20 
## regulation of actin filament-based proce... 
##                                          23 
##     positive regulation of cellular process 
##                                          29
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
##                                          12 
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
##                                 13 
## 
## $`Cluster 4`
##              regulation of transport trivalent inorganic cation transport 
##                                    6                                   10 
## 
## $`Cluster 5`
## organophosphate catabolic process 
##                                 9 
## 
## $`Cluster 6`
##           negative regulation of cell death 
##                                          11 
## sporulation resulting in formation of a ... 
##                                          14 
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
##                   protein metabolic process 
##                                           1 
##          cellular protein metabolic process 
##                                           2 
##             macromolecule metabolic process 
##                                           5 
##                  macromolecule modification 
##                                           6 
##           glycoprotein biosynthetic process 
##                                          10 
##                       protein glycosylation 
##                                          11 
##                    translational initiation 
##                                          12 
##    cellular macromolecule metabolic process 
##                                          20 
##              glycoprotein metabolic process 
##                                          25 
##                protein modification process 
##                                          35 
##      regulation of translational initiation 
##                                          37 
##                   regulation of translation 
##                                          62 
##       cellular protein modification process 
##                                          77 
## regulation of cellular protein metabolic... 
##                                          80 
## 
## $`Cluster 2`
## small GTPase mediated signal transductio... 
##                                           3 
##           regulation of signal transduction 
##                                          46 
## regulation of Rab protein signal transdu... 
##                                          66 
##             Rab protein signal transduction 
##                                          69 
## 
## $`Cluster 3`
##                           glycosylation 
##                                       4 
##            glycolipid metabolic process 
##                                      38 
##                    histone modification 
##                                      43 
##        membrane lipid metabolic process 
##                                      52 
##          sphingolipid metabolic process 
##                                      55 
##              cell cycle DNA replication 
##                                      63 
## hexachlorocyclohexane metabolic process 
##                                      64 
##          phospholipid metabolic process 
##                                      70 
##             macromolecule glycosylation 
##                                      71 
##                 ether metabolic process 
##                                      82 
##                     histone acetylation 
##                                      92 
## 
## $`Cluster 4`
##      single-organism organelle organization 
##                                           7 
##                      organelle organization 
##                                          13 
##                            membrane docking 
##                                          24 
##                   hair cell differentiation 
##                                          31 
##                     muscle cell development 
##                                          48 
##                             vesicle docking 
##                                          56 
##             mechanoreceptor differentiation 
##                                          60 
##                       inner ear development 
##                                          72 
##             epithelial cell differentiation 
##                                          76 
## regulation of actin cytoskeleton organiz... 
##                                          81 
##      vesicle docking involved in exocytosis 
##                                          85 
##                 muscle cell differentiation 
##                                          94 
## imaginal disc-derived wing hair organiza... 
##                                          98 
## 
## $`Cluster 5`
## nucleoside triphosphate metabolic proces... 
##                                           8 
##                       GTP metabolic process 
##                                          14 
##            ribonucleotide metabolic process 
##                                          16 
## regulation of purine nucleotide metaboli... 
##                                          18 
##                       GTP catabolic process 
##                                          19 
## regulation of nucleotide metabolic proce... 
##                                          21 
## regulation of nucleoside metabolic proce... 
##                                          22 
##                nucleotide catabolic process 
##                                          23 
## regulation of nucleotide catabolic proce... 
##                                          26 
## guanosine-containing compound catabolic ... 
##                                          27 
## regulation of purine nucleotide cataboli... 
##                                          28 
## guanosine-containing compound metabolic ... 
##                                          29 
##      nucleoside phosphate catabolic process 
##                                          36 
## positive regulation of Rab GTPase activi... 
##                                          45 
##           regulation of Ras GTPase activity 
##                                          54 
##         purine nucleotide metabolic process 
##                                          65 
## purine-containing compound catabolic pro... 
##                                          73 
## 
## $`Cluster 6`
##                     intracellular transport 
##                                           9 
##       establishment of localization in cell 
##                                          15 
##                       cellular localization 
##                                          17 
##                  vesicle-mediated transport 
##                                          41 
##                  macromolecule localization 
##                                          44 
##              detection of external stimulus 
##                                          51 
##                regulation of blood pressure 
##                                          53 
##       establishment of protein localization 
##                                          58 
##                          chloride transport 
##                                          75 
## sensory perception of mechanical stimulu... 
##                                          84 
##               regulation of anion transport 
##                                          88 
##                           protein transport 
##                                          96 
##               detection of abiotic stimulus 
##                                          97 
## 
## $`Cluster 7`
## posttranscriptional regulation of gene e... 
##                                          30 
## regulation of phosphate metabolic proces... 
##                                          32 
##             regulation of catabolic process 
##                                          39 
##    regulation of cellular catabolic process 
##                                          40 
## regulation of phosphorus metabolic proce... 
##                                          49 
##              regulation of cellular process 
##                                          59 
##          regulation of response to stimulus 
##                                          74 
##                       biological regulation 
##                                          78 
##            regulation of biological process 
##                                          89 
## regulation of actin filament-based proce... 
##                                          90 
##                     regulation of signaling 
##                                          99 
## 
## $`Cluster 8`
## serine family amino acid metabolic proce... 
##                                          33 
##                      recombinational repair 
##                                          42 
## flavin-containing compound metabolic pro... 
##                                          47 
##                riboflavin metabolic process 
##                                          50 
##               glutathione metabolic process 
##                                          79 
##                       amino acid activation 
##                                          83 
## chondroitin sulfate proteoglycan biosynt... 
##                                          86 
##    chondroitin sulfate biosynthetic process 
##                                          87 
##                   peptide metabolic process 
##                                          91 
##             folic acid biosynthetic process 
##                                          95 
##                      mRNA metabolic process 
##                                         100 
## 
## $`Cluster 9`
##                  protein complex biogenesis 
##                                          34 
##             macromolecular complex assembly 
##                                          57 
##                    protein complex assembly 
##                                          61 
## positive regulation of protein polymeriz... 
##                                          67 
##        regulation of protein polymerization 
##                                          68 
##      regulation of protein complex assembly 
##                                          93
```

GO enrichment for *A. picea* **Bimodal** genes.


```r
# A22 'Bim' genes
A22.geneList.bim <- rep(0, length = length(Ap.geneList))
names(A22.geneList.bim) <- names(Ap.geneList)
A22.geneList.bim[(which(names(A22.geneList.bim) %in% A22.bim.transcripts))] <- 1
A22.bim.gsea <- gsea(genelist = A22.geneList.bim, geneID2GO = geneID2GO)
```

Clustering of enriched GO terms for *A. picea* **Low** genes.


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
##                                          10 
## 
## $`Cluster 2`
## regulation of intracellular signal trans... 
##                                           2 
##         I-kappaB kinase/NF-kappaB signaling 
##                                           5 
##        glutamate receptor signaling pathway 
##                                           9 
## 
## $`Cluster 3`
##                  protein phosphorylation 
##                                        3 
##             protein modification process 
##                                        4 
##               macromolecule modification 
##                                        6 
##    cellular protein modification process 
##                                        8 
## cellular macromolecule metabolic process 
##                                       11 
## 
## $`Cluster 4`
## intracellular protein transport 
##                               7
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
##                                 7                                11 
## 
## $`Cluster 3`
##    nitric oxide metabolic process nitric oxide biosynthetic process 
##                                 9                                10 
## 
## $`Cluster 4`
## embryo sac egg cell differentiation 
##                                  12
```




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
Ap.setdiff.gsea <- rbind(A22.gsea, Ar.gsea)
# reorder
Ap.setdiff.gsea <- Ap.setdiff.gsea[,c("Species", "Type", "GO.ID", "Term", "Annotated", "Significant", "Expected", "parentchild")]
colnames(Ap.setdiff.gsea)[8] <- "P"

write.csv(Ap.setdiff.gsea, file = paste(resultsdir, "Ap_setdiff_GSEA_", Sys.Date(), ".csv", sep = ""), row.names = FALSE)
```







```r
# similarity among enriched GO terms
A22.low.GO.term.sim <- termSim(A22.low.gsea$GO.ID, A22.low.gsea$GO.ID, ont = "BP", organism = "fly")
# distance matrix
dist.A22.low.GO.term.sim <- dist(A22.low.GO.term.sim)
plot(hclust(dist(A22.low.GO.term.sim)))
```

![plot of chunk GOSemSim](figure/GOSemSim1.png) 

```r
A22.high.GO.term.sim <- termSim(A22.high.gsea$GO.ID, A22.high.gsea$GO.ID, ont = "BP", organism = "fly")
# distance matrix
plot(hclust(dist(A22.high.GO.term.sim)))
```

![plot of chunk GOSemSim](figure/GOSemSim2.png) 


## Visualize responsive transcripts

Make plots for all genes expressed at *High* temps in GO category "GO:0006950: response to stress"



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
trp_A22_poslinear <- ggplot(resp.TPM.dt.sub[A22.poslinear.transcripts][1:220,], aes(x=val, y=TPM.scaled, group=Transcript)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) + 
  facet_grid(. ~ colony2) + 
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
```

```
## Error: object 'A22.poslinear.transcripts' not found
```

```r
print(trp_A22_poslinear)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'print': Error: object 'trp_A22_poslinear' not found
```

```r
### predicted values
resp.TPM.dt.sub.pred[,pTPM.scaled:=scale(pTPM), by = Transcript]
```

```
##                          Transcript    TPM  val colony pTPM pTPM.scaled
##      1: 100008|*|comp137625_c0_seq2 0.0000  0.0    A22 1.01      -0.610
##      2: 100008|*|comp137625_c0_seq2 0.0000  3.5    A22 1.01      -0.640
##      3: 100008|*|comp137625_c0_seq2 0.0743 10.5    A22 1.01      -0.663
##      4: 100008|*|comp137625_c0_seq2 0.0000 14.0    A22 1.01      -0.656
##      5: 100008|*|comp137625_c0_seq2 0.0000 17.5    A22 1.01      -0.636
##     ---                                                                
## 200196:      9|*|comp147140_c0_seq1 0.6444 24.5     Ar 1.67      -0.494
## 200197:      9|*|comp147140_c0_seq1 0.5862 28.0     Ar 1.64      -0.629
## 200198:      9|*|comp147140_c0_seq1 0.7246 31.5     Ar 1.59      -0.813
## 200199:      9|*|comp147140_c0_seq1 0.4516 35.0     Ar 1.53      -1.040
## 200200:      9|*|comp147140_c0_seq1 0.4600 38.5     Ar 1.47      -1.306
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
pred_A22_poslinear <- ggplot(resp.TPM.dt.sub.pred[A22.poslinear.transcripts][1:220,], aes(x=val, y=pTPM.scaled, group=Transcript)) +
  geom_line() + 
  facet_grid(. ~ colony) + 
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
```

```
## Error: object 'A22.poslinear.transcripts' not found
```

```r
print(pred_A22_poslinear)
```

```
## Error: error in evaluating the argument 'x' in selecting a method for function 'print': Error: object 'pred_A22_poslinear' not found
```

```r
ggplot(resp.TPM.dt.sub.pred["100348|*|comp3450455_c0_seq1"], aes(x=val, y=pTPM.scaled, group=Transcript)) +
  geom_line() + 
  facet_grid(. ~ colony) + 
  scale_y_continuous(name="Expression (scaled)") +
  scale_x_continuous(name=expression(paste("Temperature ", degree, "C")))
```

![plot of chunk ggplot_high](figure/ggplot_high3.png) 

### Overlap between *Ar* Intermediate and *A22* Low genes

An observation from GO analysis is that there is an overlap between enriched terms in the *Low* group for *A22* and the *Intermediate* group for *Ar*. A potential explanation for this pattern is if the *Low* genes in *A22* shut down at high temperatures, but not low temperatures, rather than are up-regulated at low temperatures.


```r
# transcripts where TPM increases at low temps & decreases at high temps
A22.low.Ar.int.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(intersect(A22.low,Ar.int), "A22")]

for(i in unique(A22.low.Ar.int.TPM.dt.sub$Transcript)[1:10]) {
    s <- A22.low.Ar.int.TPM.dt.sub[i]
    X11()
    plot(s$val, s$pTPM, ylim = c(0, max(s$pTPM)+.1))
    points(s$val, s$TPM, pch = 8, col = "red")
    if(s[which(s$val == 38.5) , pTPM] == min(s$pTPM)) print("OffHigh") else {
        print("OnLow")
    }
}
```

```
## [1] "OnLow"
```

```
## [1] "OffHigh"
```

```
## [1] "OffHigh"
```

```
## [1] "OnLow"
```

```
## [1] "OffHigh"
```

```
## [1] "OffHigh"
```

```
## [1] "OnLow"
```

```
## [1] "OnLow"
```

```
## [1] "OnLow"
```

![plot of chunk gene_overlap](figure/gene_overlap.png) 

```
## [1] "OnLow"
```

```r
# tabulate the number of transcripts that have lowest expression at maximum temperature "OffHigh"
A22.low.Ar.int.low.type <- ddply(A22.low.Ar.int.TPM.dt.sub, .(Transcript), summarise, lt = ifelse(pTPM[which(val == 38.5)] == min(pTPM), "OffHigh", "OnLow"))

table(A22.low.Ar.int.low.type$lt)
```

```
## 
## OffHigh   OnLow 
##     560     790
```

```r
# A22 low Ar low
A22.low.Ar.low.TPM.dt.sub <- resp.TPM.dt.sub.pred[J(intersect(A22.low, Ar.low), "A22")]

A22.low.Ar.low.low.type <- ddply(A22.low.Ar.low.TPM.dt.sub, .(Transcript), summarise, lt = ifelse(pTPM[which(val == 38.5)] == min(pTPM), "OffHigh", "OnLow"))

# table results and perform chi-square test
low.table <- rbind(LowInt = table(A22.low.Ar.int.low.type$lt), LowLow = table(A22.low.Ar.low.low.type$lt))
low.table
```

```
##        OffHigh OnLow
## LowInt     560   790
## LowLow    1540  1054
```

```r
chisq.test(low.table)
```

```
## 
## 	Pearson's Chi-squared test with Yates' continuity correction
## 
## data:  low.table
## X-squared = 113, df = 1, p-value < 2.2e-16
```

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
##  [1] GOSemSim_1.22.0        Rcpp_0.11.1            RDAVIDWebService_1.2.0
##  [4] GOstats_2.30.0         Category_2.30.0        Matrix_1.1-3          
##  [7] Rgraphviz_2.8.1        topGO_2.16.0           SparseM_1.03          
## [10] GO.db_2.14.0           RSQLite_0.11.4         DBI_0.2-7             
## [13] AnnotationDbi_1.26.0   GenomeInfoDb_1.0.2     Biobase_2.24.0        
## [16] BiocGenerics_0.10.0    graph_1.42.0           reshape2_1.4          
## [19] xtable_1.7-3           MASS_7.3-33            plyr_1.8.1            
## [22] RCurl_1.95-4.1         bitops_1.0-6           data.table_1.9.2      
## [25] stringr_0.6.2          pander_0.3.8           knitcitations_0.5-0   
## [28] bibtex_0.3-6           ggplot2_1.0.0          R.utils_1.32.4        
## [31] R.oo_1.18.0            R.methodsS3_1.6.1      knitr_1.6             
## 
## loaded via a namespace (and not attached):
##  [1] annotate_1.42.0       AnnotationForge_1.6.1 codetools_0.2-8      
##  [4] colorspace_1.2-4      digest_0.6.4          evaluate_0.5.5       
##  [7] formatR_0.10          genefilter_1.46.1     GSEABase_1.26.0      
## [10] gtable_0.1.2          httr_0.3              IRanges_1.22.7       
## [13] labeling_0.2          lattice_0.20-29       munsell_0.4.2        
## [16] proto_0.3-10          RBGL_1.40.0           rJava_0.9-6          
## [19] scales_0.2.4          splines_3.1.0         stats4_3.1.0         
## [22] survival_2.37-7       tools_3.1.0           XML_3.98-1.1
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
- G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, S. Wang,   (2010) Gosemsim: an R Package For Measuring Semantic Similarity Among go Terms And Gene Products.  *Bioinformatics*  **26**  976-978  [10.1093/bioinformatics/btq064](http://dx.doi.org/10.1093/bioinformatics/btq064)
