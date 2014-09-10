#!/bin/Rscript

# list of reads in each of four trimmed classes
readdir <- "../data/ind_files/"

# sailfish index directory
system('mkdir -p ../results/trinity-full/sailfish-index-Trinity-cap3-uclust')
sfindex <- "../results/trinity-full/sailfish-index-Trinity-cap3-uclust"

# make sailfish index
system("sailfish index -t ../results/trinity-full/Trinity_cap3_uclust.fasta -o ../results/trinity-full/sailfish-index-Trinity-cap3-uclust -k 20 -p 4")

# prep for sailfish quantification
readlist <- list.files(readdir)
paired.left <- readlist[grep(".\\.paired.left.fastq$", readlist)]
paired.right <- readlist[grep("\\.paired.right.fastq$", readlist)]
unpaired.left <- readlist[grep("unpaired.left.fastq$", readlist)]
unpaired.right <- readlist[grep("unpaired.right.fastq$", readlist)]

# output directiry
system('mkdir -p ../results/trinity-full/sailfish-expression-Trinity-cap3-uclust')
sfexpressionroot <- "../results/trinity-full/sailfish-expression-Trinity-cap3-uclust/"

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
     print(sailfishcmd)
     system(sailfishcmd)
     message("Done with expression quantification for sample ", samples[j], ": ", Sys.time(), "\n")
}
