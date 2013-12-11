Aphaenogaster transcriptome assembly
=======================================

**Author:** [John Stanton-Geddes](john.stantongeddes.research@gmail.com)

Requires the following programs installed and available on path:

* [rlsim](https://github.com/sbotond/rlsim) to simulate Illumina RNAseq reads
* [simNGS](https://github.com/timmassingham/simNGS) to simulate reads 
* [Trim_Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim <- galore/)
  - requires [cutadapt](https://code.google.com/p/cutadapt/) and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [khmer](https://github.com/ged-lab/khmer) to perform digital normalization and filtering
* [Trinity](http://trinityrnaseq.sourceforge.net/) for transcriptome assembly


```{r setup}
# Global settings
options(stringAsFactors=FALSE)

# Load libraries
library(ggplot2)
library(knitr)
library(Biostrings)

# knitr options
opts_chunk$set(cache=TRUE)

# Load personal functions
source("assemble_functions.R")
```

## *in silico* spike-in ##

Simulate Illumina reads from known mRNA sequence to "spike-in" during transcriptome assembly.
Use these simluated reads to evaluate assembly at end.

Fasta file of known mRNA transcripts for simulated reads with known expression values.
I use a file with 1000 *Arabidopsis* mRNA transcripts downloaded from [European Nucleotide Archive](http://www.ebi.ac.uk/ena/home) 

```{r simulate, eval=FALSE}
knownfasta <- "../data/ena.fasta"
system("mkdir -p ../data/sim")
simout <- "../data/sim"

# Randomly select X sequences; output "known.fasta"
sample.fasta(knownfasta, 100)

## Move intermediate files
system(paste("mv mean-length-fasta.txt ", simout, sep=""))

## Trim fasta titles to just ID
system("cut -f1 -d' ' known.fasta > known1.fasta")

## Assign expression-level to each sequence using the `sel` tool in [rlsim](https://github.com/sbotond/rlsim)

system("sel -d '1.0:g:(1500, 3)' known1.fasta > known-sim.fasta")
system("rm known.fasta known1.fasta") # Clean-up temp files
message("Done with expression levels")

## Simulate fragments created during library prep using `rlsim`. Based on length of fragments estimated when using FLASH to pair reads, use empirical length distribution of 180 bp (SD: 20bp).

system("rlsim -v -n 150000 -d '1:n:(180, 20, 100, 500)' known-sim.fasta > known-sim-frags.fasta")
message("Done with simulated fragmentation levels")

## Generate simulated Illumina paired-end reads using [simNGS](http://www.ebi.ac.uk/goldman-srv/simNGS/) and a default runfile provided with simNGS
# output is 'reads_end1.fq' and 'reads_end2.fq'
system("cat known-sim-frags.fasta | simNGS -p paired -o fastq -O reads /opt/software/simNGS/data/s_3_4x.runfile")
message("Done with Illumina read simulation")

# move intermediate files
system(paste("mv rlsim_report.json ", simout, sep=""))
system(paste("mv sel_report.pdf ", simout, sep=""))

# move files to data directory
system("mv reads_end1.fq ../data/A22-si-R1.fastq")
system("mv reads_end2.fq ../data/A22-si-R2.fastq")
system("mv known-sim.fasta ../data/sim")
system("mv known-sim-frags.fasta ../data/sim")
```

## QC ##

Clip adapter sequence from reads and trim poor quality reads using Trim_Galore

```{r trimclip, eval=FALSE}
# Data directory 
datadir <- "/home/data/Aphaeno_transcriptome/130509_SN1073_0326_BD25DAACXX/Project_Stanton-Geddes_Project_001"
# Directory to store results
system("mkdir -p ../data/trimclip")
clipdir <- "../data/trimclip" 

# List of samples
samplist <- list.files(datadir)
(reads1 <- samplist[grep("A[r2]{1,2}-..-R1.fastq", samplist)])
(reads2 <- samplist[grep("A[r2]{1,2}-..-R2.fastq", samplist)])

# Loop across all the data files and trimclip

for (i in 1:length(reads1)) {
    # Run `trim_galore`.
    r1 <- paste(datadir, reads1[i], sep="/")
    r2 <- paste(datadir, reads2[i], sep="/")
    system(paste("trim_galore --quality 20 --phred33 --fastqc_args '--o", clipdir, "' --length 20 --paired --output_dir", clipdir, r1, r2, sep=" "))
    message("Done with QC for sample: ", reads1[i])
}

# Trimclip simulated spike-in reads
system(paste("trim_galore --quality 20 --phred33 --fastqc_args '--o", clipdir, "' --length 20 --paired --output_dir", clipdir, "../data/A22-si-R1.fastq ../data/A22-si-R2.fastq", sep=" "))
```


## Merge overlapping paired-end reads

With standard insert size of 200 bp and 2 x 100bp paired-end sequencing, many paired-reads are overlapping. FLASH attempts to merge paired-end reads when possible.

```{r FLASH, eval=FALSE}
# make new directory for merged reads
system("mkdir -p ../data/merged")
mergedir <- "../data/merged/"

# list of trimclipped samples
clippedlist <- list.files(clipdir)
(clippedR1 <- clippedlist[grep("A[r2]{1,2}-..-R1_val_1.fq$", clippedlist)])
(clippedR2 <- clippedlist[grep("A[r2]{1,2}-..-R2_val_2.fq$", clippedlist)])

# loop across data files and merge overlapping reads

for (j in 1:length(clippedR1)) {
    outpre <- strsplit(clippedR1[j], split="-R")[[1]][1]
    sampR1 <- paste(clipdir, clippedR1[j], sep="/")
    sampR2 <- paste(clipdir, clippedR2[j], sep="/")
        
    system(paste("flash --phred-offset 33 --interleaved-output --output-directory", mergedir, "--output-prefix", outpre, sampR1, sampR2, sep=" "))
    message("Done with merging reads for sample ", clippedR1[j], ": ", Sys.time())
}
```

## Digital normalization

Normalize reads separately for samples from each colony using [khmer](https://github.com/ged-lab/khmer) software


```{r A22_diginorm, eval=TRUE}
# set PYTHONPATH to khmer module
system("export PYTHONPATH=/opt/software/khmer/python")

# make new directory for diginorm reads
system("mkdir -p ../data/diginorm")
diginormdir <- "../data/diginorm/"

# run `normalize-by-median` for extendedFragments with coverage threshold and kmer of 20. -N 4 -x 3e9 allocates up to 12GB RAM. 

message("Start diginorm on A22 extendedFrags: ", Sys.time())
system("/opt/software/khmer/scripts/normalize-by-median.py -R diginorm.out -C 20 -k 20 -N 4 -x 3e9 --savehash ../data/diginorm/A22-diginorm-C20k20.kh ../data/merged/A22*.extendedFrags.fastq")
message("Done with diginorm on A22 extendedFrags: ", Sys.time())

# to run `normalize-by-median` for notCombined reads, need to add Illumina 1.3 style tags /1 or /2 to reads

# list of trimclipped samples
mergedir <- "../data/merged/"
mergedlist <- list.files(mergedir)
(notCombined <- mergedlist[grep("A22-...notCombined", mergedlist)])

# loop across files, adding tags
for (s in notCombined) {
    system(paste("python add-illumina-tags-interleaved.py ", mergedir, s, " > ", mergedir, s, ".out", sep=""))
}

# run `normalize-by-median`. load hash of previously saved read filtering
message("Start diginorm on A22 notCombined reads: ", Sys.time())
system("/opt/software/khmer/scripts/normalize-by-median.py -R diginorm-paired.out -p -C 20 -k 20 -N 4 -x 4e9 --loadhash ../data/diginorm/A22-diginorm-C20k20.kh --savehash ../data/diginorm/A22-diginorm-C20k20.kh ../data/merged/A22*.notCombined.fastq.out")

message("Done with diginorm on A22 notCombined reads: ", Sys.time())
```

Trim low abundance parts of high coverage reads - these are likely erroneous.
Note that this will orphan some reads with poor quality partners

```{r A22_diginorm_trim, eval=TRUE}
message("Trim low abundance k-mers: ", Sys.time())
system("/opt/software/khmer/scripts/filter-abund.py -V ../data/diginorm/A22-diginorm-C20k20.kh *.keep")

# Separate orphaned from still-paired reads in .notCombined.fastq.out.keep.abundfilt
message("Separate orphaned from still-paired A22 reads: ", Sys.time())
dnlist <- list.files()
(dn <- dnlist[grep("A22-...notCombined.fastq.out.keep.abundfilt", dnlist)])

for (n in dn) {
    system(paste("/opt/software/khmer/scripts/extract-paired-reads.py ", n, sep=""))
}

# Move final files and cleanup

message("Move final files and clean-up: ", Sys.time())

system("mv *.extendedFrags.fastq.keep.abundfilt ../data/diginorm/")
system("mv *.abundfilt.[ps]e ../data/diginorm/")

system("mv *.keep ../data/diginorm/")
system("mv *notCombined.fastq.out.keep.abundfilt ../data/diginorm/")
```

**Repeat diginorm for Ar colony**


Normalize reads separately for samples from each colony using [khmer](https://github.com/ged-lab/khmer) software


```{r Ar_diginorm, eval=TRUE}

# run `normalize-by-median` for extendedFragments with coverage threshold and kmer of 20. -N 4 -x 3e9 allocates up to 12GB RAM. 
message("Start diginorm on Ar extendedFrags: ", Sys.time())
system("/opt/software/khmer/scripts/normalize-by-median.py -R diginorm.out -C 20 -k 20 -N 4 -x 3e9 --savehash ../data/diginorm/Ar-diginorm-C20k20.kh ../data/merged/Ar*.extendedFrags.fastq")
message("Done with diginorm on Ar extendedFrags: ", Sys.time())

# to run `normalize-by-median` for notCombined reads, need to add Illumina 1.3 style tags /1 or /2 to reads

# list of trimclipped samples
mergedir <- "../data/merged/"
mergedlist <- list.files(mergedir)
(notCombined <- mergedlist[grep("Ar-...notCombined", mergedlist)])

# loop across files, adding tags
for (s in notCombined) {
    system(paste("python add-illumina-tags-interleaved.py ", mergedir, s, " > ", mergedir, s, ".out", sep=""))
}

# run `normalize-by-median`. load hash of previously saved read filtering
message("Start diginorm on notCombined reads: ", Sys.time())
system("/opt/software/khmer/scripts/normalize-by-median.py -R diginorm-paired.out -p -C 20 -k 20 -N 4 -x 4e9 --loadhash ../data/diginorm/Ar-diginorm-C20k20.kh --savehash ../data/diginorm/Ar-diginorm-C20k20.kh ../data/merged/Ar*.notCombined.fastq.out")
message("Done with diginorm on Ar notCombined reads: ", Sys.time())
```

Trim low abundance parts of high coverage reads - these are likely erroneous. 
Note that this will orphan some reads with poor quality partners

```{r Ar_diginorm_trim, eval=TRUE}
message("Trim low abundance k-mers: ", Sys.time())
system("/opt/software/khmer/scripts/filter-abund.py -V ../data/diginorm/Ar-diginorm-C20k20.kh *.keep")

# Separate orphaned from still-paired reads in .notCombined.fastq.out.keep.abundfilt
message("Separate orphaned from still-paired reads: ", Sys.time())
dnlist <- list.files()
(dn <- dnlist[grep("Ar-...notCombined.fastq.out.keep.abundfilt", dnlist)])

for (n in dn) {
    system(paste("/opt/software/khmer/scripts/extract-paired-reads.py ", n, sep=""))
}

# Move final files and cleanup

message("Move final files and clean-up: ", Sys.time())

system("mv *.extendedFrags.fastq.keep.abundfilt ../data/diginorm/")
system("mv *.abundfilt.[ps]e ../data/diginorm/")

system("mv *.keep ../data/diginorm")
system("mv *notCombined.fastq.out.keep.abundfilt ../data/diginorm/")
```

## Transcriptome assembly

Use Trinity for transcriptome assembly.
Requires reads split into 'left' and 'right' files.
Combine all *.notCombined into reads 1 and 2

```{r trinity_prep, eval=FALSE}
# split paired-end .notCombined reads that passed through digital normalization
diginormdir <- "../data/diginorm/"
diginormlist <- list.files(diginormdir)
(pe <- diginormlist[grep("A[r2]{1,2}-...notCombined.fastq.out.keep.abundfilt.pe$", diginormlist)])

for (p in pe) {
    system(paste("/opt/software/khmer/scripts/split-paired-reads.py ", diginormdir, p, sep=""))
}

# move split files
system(paste("mv A*.abundfilt.pe.1 ", diginormdir, sep=""))
system(paste("mv A*.abundfilt.pe.2 ", diginormdir, sep=""))

# add paired-end tag to 'extendedFrags.fastq.keep.abundfilt' for Trinity
eflist <- list.files(diginormdir)
(extendedFrags <- eflist[grep("extendedFrags.fastq.keep.abundfilt", eflist)])

# loop across files, adding tags
for (e in extendedFrags) {
    system(paste("python add-illumina-tags.py ", diginormdir, e, " 1 > ", diginormdir, e, ".tag", sep=""))
}
```


Run Trinity

```{r A22-trinity, eval=FALSE}
diginormdir <- "../data/diginorm/"
A22trinitydir <- paste("../results/A22-trinity", Sys.Date(), sep="-")
system(paste("mkdir -p ", A22trinitydir, sep=""))

# concatenate files into left and right
system(paste("cat ", diginormdir, "A22*.1 > ", A22trinitydir, "/A22-r1.fq", sep=""))
system(paste("cat ", diginormdir, "A22*.2 > ", A22trinitydir, "/A22-r2.fq", sep=""))

# add single-end unpaired reads to one file
system(paste("cat ", diginormdir, "A22*.notCombined.fastq.out.keep.abundfilt.se >> ", A22trinitydir, "/A22-r1.fq", sep=""))

# add single-end extendedFrags to one file
system(paste("cat ", diginormdir, "A22*.extendedFrags.fastq.keep.abundfilt.tag >> ", A22trinitydir, "/A22-r1.fq", sep=""))

system(paste("Trinity.pl --seqType fq --JM 50G --left ", A22trinitydir, "/A22-r1.fq --right ", A22trinitydir, "/A22-r2.fq --output ", A22trinitydir, sep=""))

# summary statistics
system(paste("python assemstats2.py 100 ", A22trinitydir, "Trinity.fasta", sep=""))
```

Trinity assembly completed for A22!

```{r Ar-trinity, eval=FALSE}
diginormdir <- "../data/diginorm/"
Artrinitydir <- paste("../results/Ar-trinity", Sys.Date(), sep="-")
system(paste("mkdir -p ", Artrinitydir, sep=""))

# concatenate files into left and right
system(paste("cat ", diginormdir, "Ar*.1 > ", Artrinitydir, "/Ar-r1.fq", sep=""))
system(paste("cat ", diginormdir, "Ar*.2 > ", Artrinitydir, "/Ar-r2.fq", sep=""))

# add single-end unpaired reads to one file
system(paste("cat ", diginormdir, "Ar*.notCombined.fastq.out.keep.abundfilt.se >> ", Artrinitydir, "/Ar-r1.fq", sep=""))

# add single-end extendedFrags to one file
system(paste("cat ", diginormdir, "Ar*.extendedFrags.fastq.keep.abundfilt.tag >> ", Artrinitydir, "/Ar-r1.fq", sep=""))

message("Start Trinity assembly for Ar: ", Sys.time())

system(paste("Trinity.pl --seqType fq --JM 50G --left ", Artrinitydir, "/Ar-r1.fq --right ", Artrinitydir, "/Ar-r2.fq --output ", Artrinitydir, sep=""))

# summary statistics
system(paste("python assemstats2.py 100 ", Artrinitydir, "Trinity.fasta", sep=""))
```

Yay! Trinity assembly completed for Ar!

Remove intermediate files and compress final fastq files to save disk space

```{r assembly_cleanup, eval=FALSE}
# Keep only files directly needed for Trinity assembly
setwd(mergedir)
system("ls | grep -P 'A.*.notCombined.fastq$' | xargs -d'\n' rm")

setwd("../diginorm")
system("ls | grep -P 'A.*.abundfilt$' | xargs -d'\n' rm")
system("ls | grep -P 'A.*.extendedFrags.fastq.keep$' | xargs -d'\n' rm")
system("ls | grep -P 'A.*.extendedFrags.fastq.keep.abundfilt$' | xargs -d'\n' rm")

system("ls | grep -P 'A.*.notCombined.fastq.out.keep.abundfilt.pe$' | xargs -d'\n' rm")


# compress diginorm directory
setwd("../")
system("tar -czf diginorm_fastq.tar.gz diginorm")

# remove directories
system("rm -rf diginorm merged trimclip")

setwd("../scripts")
```

## Evaluate assembly

To evaluate assembly, blast transcripts against known spike-in reads

```{r evaluate_assembly, eval=TRUE}
# load spike-in fasta file
spikein <- "../data/sim/known-sim.fasta"
spikeinfa <- readDNAStringSet(spikein)
length(names(spikeinfa))

# make BLAST database from spike-in fasta
system(paste("makeblastdb -dbtype nucl -in ", spikein, sep=""))

# BLAST assembled transcripts against database
transcripts <- "../results/trinity/Trinity.fasta"
system(paste("blastn -query ", transcripts, " -db ", spikein, " -outfmt 6 -out blast_spikein.txt",sep=""))

# evaluate
blastout <- "blast_spikein.txt"
assembly.spikein.eval(spikein, blastout)
```

## Transcriptome filtering

Reduce redundant contigs by running

* [uclust](http://drive5.com/usearch/manual/uclust_algo.html) to cluster similar sequences with 90% similarity 
* [CAP3](http://genome.cshlp.org/content/9/9/868.long) to merge contigs that overlap by with similarity 

```{r cap3}

## Cluster and assemble for A22

# cap3. allow gap up to 20 bp (-a 20) and require 90% similarity (-p 90)
system(paste("cap3 ", A22trinitydir, "Trinity.fasta -f 20 -a 20 -k 0 -p 90 -o 100 > ", A22trinitydir, "Trinity_cap3.fasta", sep=""))

## uclust. sort then cluster
system(paste("uclust --sort ", A22trinitydir, "Trinity_cap3.fasta --output ", A22trinitydir,"Trinity_cap3_sorted.fasta", sep=""))
system(paste("uclust --input ", A22trinitydir, "Trinity_cap3_sorted.fasta --uc ", A22trinitydir,"Trinity_cap3_uclust.fasta --id 0.90", sep=""))



# cap3. allow gap up to 20 bp (-a 20) and require 90% similarity (-p 90)
system(paste("cap3 ", Artrinitydir, "test.fasta -f 20 -a 20 -k 0 -p 90 -o 100 > ", Artrinitydir, "test_cap3.out", sep=""))

# concatenate merged 'contigs' and unassembled 'singlets'
system(paste("cat ", Artrinitydir, "test.fasta.cap.contigs ", Artrinitydir, "test.fasta.cap.singlets > ", Artrinitydir, "test_cap3.fasta", sep=""))

## uclust. sort then cluster
system(paste("uclust --sort ", A22trinitydir, "test_cap3.fasta --output ", A22trinitydir,"test_cap3_sorted.fasta", sep=""))
system(paste("uclust --input ", A22trinitydir, "test_cap3_sorted.fasta --uc ", A22trinitydir,"test_cap3_uclust.fasta --id 0.90", sep=""))



```


## Identify thermally response genes
