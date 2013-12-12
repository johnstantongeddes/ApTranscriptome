## Functions for transcriptome assembly, evaluation and evaluation
## Author: [John Stanton-Geddes](john.stantongeddes.research@uvm.edu)
## Date: 2013-12-11

######################################################
## sample.fasta(filein, X)
##
## Usage:
##
##     sample.fasta(filein, X)
##
##     where filein is a fasta file
##     and X is the number of transcripts to sample
##     output is a fasta file containing the X randomly
##     sampled transcripts. Also outputs mean length of the
##     sampled transcripts


sample.fasta <- function(filein, X) {
    ##-------------Load libraries--------------------------------
    require(Biostrings)

    ## Read known fasta file
    file <- readDNAStringSet(filein)
    length(file)
    mean(width(file))

    # Sample X sequences
    samp.file <- sample(file, X)
    length(samp.file)
    head(samp.file)

    # Write subset of fasta sequences to file
    writeXStringSet(samp.file, file="known.fasta")

    # Write mean length to file
    cat(mean(width(samp.file)), '\n', file="mean-length-fasta.txt")
}


#####################################################
## assembly.spikein.eval(fasta, blastout)
##
## Usage: 
##
##     assembly.spikein.eval(fasta, blastout)
##
##     where fasta is the 'spikein' fasta file and
##     blastout is the result of running blastn of the
##     transcriptome assembly against the known seqeunces
##     specifying `-outfmt 6` to allow reading with read.blast() 
##
## Description:
##     Script to evaluate simulated transcriptome assembly
##     Takes in list of known 'spikein' fasta transcripts
##     genes and results of BLAST output of transcriptome assembly
##     against these spikein transcripts 
##
##     Reports the starting number of genes, number of genes
##     captured in assembly, number of genes missing in assembly
##     and the average length of each of genes in these groups.
##     Also produces histogram showing distribution of captured
##     versus missing genes.
##
#####################################################

assembly.spikein.eval <- function(fasta, blastout) {
    ##-------------Load libraries--------------------------------
    #source("http://bioconductor.org/biocLite.R")
    #biocLite("ShortRead")
    require(Biostrings)
    require(stringr)
    require(RFLPtools)
    require(ggplot2)

    # read fasta file
    known <- readDNAStringSet(fasta)
    # extract names
    known.names <- names(known)
    head(known.names)

    # read blastn file
    suppressWarnings(blast.res <- read.blast(blastout))
    head(blast.res)
    dim(blast.res)
    subject.id.names <- unique(blast.res$subject.id)
    length(unique(subject.id.names))

    cat("Done reading files", '\n')


    ##-----------Directory for results---------------------------------------

    # prefix for output files
    pre <- str_split_fixed(blastout, "\\.", 2)[,1]

    # make output directory
    system(paste("mkdir -p ", pre,"_results", sep=""))
    outdir <- paste(pre,"_results/", sep="")

    ##------------Evaluate assembly against known starting sequences-----------

    # get names of known genes that were mapped in assembly
    captured <- which(known.names %in% subject.id.names)
    length(captured)
    captured.names <- known.names[captured]

    # If some query transcripts are not mapped in assembly, report
    # differences in mean length of missing versus assembled transcripts

    if(length(captured.names) < length(known.names)) {
        known.length <- width(known)
        known.length.captured <- known.length[captured]
        ml.capt <- mean(known.length.captured)
        num.capt <- length(known.length.captured)

        known.length.missing <- known.length[-captured]
        ml.miss <- mean(known.length.missing)
        num.miss <- length(known.length.missing)
        
        # combine two dataframes for histogram
        kcap <- data.frame(status="captured", length=known.length.captured)
        kmis <- data.frame(status="missing", length=known.length.missing)

        kdf <- rbind(kcap, kmis)

                                        # plot
        png(paste(outdir, "hist-", pre, "-missing-vs-captured-length.png", sep=""))
        ggplot(kdf, aes(length, fill=status)) +
            geom_density(alpha=0.2)
        dev.off()

    } else { ml.capt = mean(width(known)); num.capt = length(unique(subject.id.names));
             ml.miss = "NA"; num.miss = 0}


    ##--------Evaluate quality of assembly-------------------------------
    # of known transcripts that are mapped, what proportion of their length is mapped?
    # are more transcripts inferred than original transcripts?

    # Dataframe to collect results
    qual.df <- matrix(ncol=4, nrow=0)
    colnames(qual.df) <- c("gene", "prop.length.mapped", "num.isoforms", "bp.mapped.to.known")

    for(i in subject.id.names) {
        # get number of transcripts
        bin <- blast.res[blast.res$subject.id==i,]
        ni <- nrow(bin)
        # length of known transcript
        kl <- width(known[which(known.names==i)])
        # what proportion of length of known transcript mapped by assembled transcripts
        length.mapped <- round(abs((max(bin$s.end) - min(bin$s.start)))/kl,2)
        # how many base pairs of assembled transcript divided by known base pairs
        bp.mapped <- sum(bin$alignment.length)/kl
        # report
        qual.df <- rbind(qual.df, c(bin$subject.id[1], length.mapped, ni, round(bp.mapped,2)))
    }

    qual.df <- as.data.frame(qual.df)
    qual.df$prop.length.mapped <- as.numeric(as.character(qual.df$prop.length.mapped))
    qual.df$num.isoforms <- as.numeric(as.character(qual.df$num.isoforms))
    qual.df$bp.mapped.to.known <- as.numeric(as.character(qual.df$bp.mapped.to.known))
    str(qual.df)

    # Mean length of original transcript mapped
    mean(qual.df$prop.length.mapped)

    # Mean copies of original transcript mapped. 1=mapped once, <1 mapped to less than full length,
    #  > 1 mapped multiple times (e.g 2=mapped twice)
    mean(qual.df$bp.mapped.to.known)

    ##-------------Report results-------------------------------------------

    ## Report results to file
    cat('Number spike-in', '\t', length(known.names), '\n',
        'Mean length spike-in (bp)', '\t', mean(width(known)), '\n',
        'Number assembled', '\t', num.capt, '\n',
        'Mean length assembled (bp)', '\t', ml.capt, '\n',
        'Mean length missing (bp)', '\t', ml.miss, '\n',
        'Mean proportion of original transcript mapped', '\t', round(mean(qual.df$prop.length.mapped),2), '\n',
        'Mean number assembled transcripts per starting', '\t',    round(mean(qual.df$num.isoforms),2), '\n', 
        'Proportion bp assembled to starting bp', '\t', round(mean(qual.df$bp.mapped),2)
        , '\n',
        file=paste(outdir, pre, "-summary.txt", sep=""))

    ##----------Clean up------------------------------------------------

    # Move BLAST file into out directory
    system(paste("mv", blastout, outdir, sep=" "))

    # Move final folder to results directory
    system(paste("mv ", outdir, " ../results", sep="")) 

    message("Completed evaluation: ", Sys.time())
}
