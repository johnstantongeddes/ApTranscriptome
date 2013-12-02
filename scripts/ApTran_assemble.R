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
options(stringAsFactors=FALSE)

library(ggplot2)
library(knitr)
opts_chunk$set(cache=TRUE)

source("assemble_functions.R")
```

## *in silico* spike-in ##

Simulate Illumina reads from known mRNA sequence to "spike-in" during transcriptome assembly.
Use these simluated reads to evaluate assembly at end.

Fasta file of known mRNA transcripts for simulated reads with known expression values.
I use a file with 1000 *Arabidopsis* mRNA transcripts downloaded from [European Nucleotide Archive](http://www.ebi.ac.uk/ena/home) 

```{r simulate}
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

```{r trimclip}
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

```{r FLASH}
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

Normalize reads using khmer


```{r diginorm, eval=FALSE}
# set PYTHONPATH to khmer module
system("export PYTHONPATH=/opt/software/khmer/python")

# make new directory for diginorm reads
system("mkdir -p ../data/diginorm")
diginormdir <- "../data/diginorm/"

# run `normalize-by-median` for extendedFragments with coverage threshold and kmer of 20. -N 4 -x 3e9 allocates up to 12GB RAM. 
message("Start diginorm on extendedFrags: ", Sys.time())
system("/opt/software/khmer/scripts/normalize-by-median.py -R diginorm.out -C 20 -k 20 -N 4 -x 3e9 --savehash ../data/merged/Ap-diginorm-C20k20.kh ../data/merged/*.extendedFrags.fastq")
message("Done with diginorm on extendedFrags: ", Sys.time())

# to run `normalize-by-median` for notCombined reads, need to add Illumina 1.3 style tags /1 or /2 to reads

# list of trimclipped samples
mergedlist <- list.files(mergedir)
(notCombined <- mergedlist[grep("A[r2]{1,2}-...notCombined", mergedlist)])

# loop across files, adding tags
for (s in notCombined) {
    system(paste("python add-illumina-tags-interleaved.py ", mergedir, s, " > ", mergedir, s, ".out", sep=""))
}

# run `normalize-by-median`. load hash of previously saved read filtering
message("Start diginorm on notCombined reads: ", Sys.time())
system("/opt/software/khmer/scripts/normalize-by-median.py -R diginorm-paired.out -p -C 20 -k 20 -N 4 -x 4e9 --loadhash ../data/merged/Ap-diginorm-C20k20.kh --savehash ../data/merged/Ap-diginorm-C20k20.kh ../data/merged/*.notCombined.fastq.out")
message("Done with diginorm on notCombined reads: ", Sys.time())
```

Trim likely erroneous k-mers. Note that this will orphan some reads with poor quality partners

```{r diginorm_trim, eval=FALSE}
message("Trim erroneous k-mers: ", Sys.time())
system("/opt/software/khmer/scripts/filter-abund.py -V ../data/merged/Ap-diginorm-C20k20.kh *.keep")

# Separate orphaned from still-paired reads in .notCombined.fastq.out.keep.abundfilt
message("Separate orphaned from still-paired reads: ", Sys.time())
dnlist <- list.files()
(dn <- dnlist[grep("A[r2]{1,2}-...notCombined.fastq.out.keep.abundfilt", dnlist)])

for (n in dn) {
    system(paste("/opt/software/khmer/scripts/extract-paired-reads.py ", n, sep=""))
}
```

Move final normalized fastq files to diginorm directory and remove intermediate files

```{r diginorm_cleanup, eval=FALSE}
message("Move final files and clean-up: ", Sys.time())
system("mv *.extendedFrags.fastq.out.keep.abundfilt ../data/diginorm/")
system("mv *.abundfilt.pe ../data/diginorm/")
system("mv *.abundfilt.se ../data/diginorm/")

system("rm *.keep")
```


## Transcriptome assemly

Use Trinity for transcriptome assembly.
Requires reads split into 'left' and 'right' files.
Combine all *.notCombined into reads 1 and 2

```{r trinity_prep, eval=FALSE}
# make directory for Trinity output
system("mkdir -p ../results/trinity")
trinitydir <- "../results/trinity/"

# split paired-end .notCombined reads that passed through digital normalization
diginormlist <- list.files(diginormdir)
(pe <- diginormlist[grep("A[r2]{1,2}-...notCombined.fastq.out.keep.abundfilt.pe", diginormlist)])

for (p in pe) {
    system(paste("/opt/software/khmer/scripts/split-paired-reads.py ", diginormdir, p, sep=""))
}

# move split files
system(paste("mv A*.abundfilt.pe.1 ", diginormdir, sep=""))
system(paste("mv A*.abundfilt.pe.2 ", diginormdir, sep=""))

# concatenate files into left and right
system(paste("cat ", diginormdir, "*.1 > ", diginormdir, "Ap-r1.fq", sep=""))
system(paste("cat ", diginormdir, "*.2 > ", diginormdir, "Ap-r2.fq", sep=""))

# add single-end unpaired reads to one file
system(paste("cat ", diginormdir, "*.notCombined.fastq.out.keep.abundfilt.se >> ", diginormdir, "Ap-r1.fq", sep=""))
# add single-end extendedFrags to one file
system(paste("cat ", diginormdir, "*.extendedFrags.fastq.keep.abundfilt >> ", diginormdir, "Ap-r1.fq", sep=""))
```

Run Trinity

```{r trinity}
diginormdir <- "../data/diginorm/"
trinitydir <- "../results/trinity/"

system(paste("Trinity.pl --seqType fq --JM 50G --left ", diginormdir, "Ap-r1.fq --right ", diginormdir, "Ap-r2.fq --output ", trinitydir, sep=""))

# summary statistics
system(paste("python assemstats2.py 100 ", trinitydir, "Trinity.fasta", sep=""))
```
