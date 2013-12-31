## Functions for transcriptome assembly, evaluation and evaluation
## Author: [John Stanton-Geddes](john.stantongeddes.research@uvm.edu)
## Date: 2013-12-11


######################################################################
## sample.fasta
######################################################################

sample.fasta <- function(filein, X) {
  # Sample X transcripts from FASTA file
  #
  # Args:
  #   filein: input FASTA file
  #   X: number of transcripts to sample
  #
  # Returns:
  #   X randomly sampled transcripts
  #   Mean length of the sampled transcripts

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


########################################################################
## assembly.spikein.eval
########################################################################

assembly.spikein.eval <- function(fasta, blastout) {
  # Evalute transcriptome assembly against in silico spike-in reads
  #
  # Args:
  #  fasta: known spike-in fasta file
  #  blastout: result of running blastn of the transcriptome assembly
  #            against the known sequences specifying `-outfmt 6`
  #            to allow reading with read.blast() 
  #
  # Returns:
  #  Reports the starting number of genes, number of genes
  #  captured in assembly, number of genes missing in assembly
  #  and the average length of each of genes in these groups.
  #  Also produces histogram showing distribution of captured
  #  versus missing genes.

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



##########################################################################
## rename.trinity.transcripts
##########################################################################

rename.trinity.transcripts <- function(fasta, prefix, maxdigits=6) {
  # Renames transcripts in fasta format generated by Trinity and optionally
  # post-processed by CAP3 (maybe other program, not tested) to 'useful' names
  #
  # Args:
  #  fasta: FASTA file
  #  prefix: Identifier common to *all* transcripts
  #  maxdigits: Number of digits in *unique* numeric ID for each transcript.
  #             Default 6 (000000)
  #
  # Returns:
  #  Same FASTA sequences with new names
    
  ##--------load libraries------------
  require(Biostrings)
  require(stringr)

  ##--------read assembly file-----------------
  file <- readDNAStringSet(fasta)
  #length(file)
  #head(file)
  #str(file)

  # pull out names
  file.names <- names(file)
  #tail(file.names)
  #length(file.names)

  # strip anything added to sequence (e.g. CAP3)
  file.names.pre <- str_split_fixed(file.names, pattern = " ", n=3)[,1]
  #head(file.names.pre)
  #tail(file.names.pre)

  # split out `seq` ID
  file.names.pre.seq <- str_split_fixed(file.names.pre, pattern="_", n=3)
  colnames(file.names.pre.seq) <- c("comp", "c", "seq")
  #head(file.names.pre.seq)
  #tail(file.names.pre.seq)
  #dim(file.names.pre.seq) 
  # for "Contig" names that do not have "seq", gave value of "seq1"
  file.names.pre.seq[which(file.names.pre.seq[,"seq"]==""),"seq"] <- "seq1"
  #head(file.names.pre.seq)

  # set starting values for component and seq
  compID <- 0
  seqID <- 1

  # vector for new names
  newnames <- vector(length=0)

  # loop across names
  for (i in 1:nrow(file.names.pre.seq)) {
    # check if seq=="seq1"
    if (file.names.pre.seq[i,"seq"] == "seq1") {
      # print status every 10k sequences
      if(i%%1000 == 0) message("Renaming sequence: ", i, "\t", Sys.time())
      # if yes, increment
      compID <- compID+1
      cID <- formatC(compID, width=maxdigits, format="d", flag="0")
      # create new name
      uniqueID <- paste(prefix, "_", cID, ".", seqID, sep="")
    } else { # if seq > 1, keep cID at same and set seqID
        seqID <- substr(file.names.pre.seq[i,"seq"], 4, 6) # select digits, allowing up to 999, following seq
        cID <- formatC(compID, width=maxdigits, format="d", flag="0")
        uniqueID <- paste(prefix, "_", cID, ".", seqID, sep="")
        } # end else
    # add new unique ID to newnames
    newnames <- c(newnames, uniqueID)
    # reset seqID to 1
    seqID <- 1   
    } # end for loop

  # Replace FASTA names with newnames
  names(file) <- newnames
  # head(file)
  # names(file)

  # Write to file
  writeXStringSet(file, file=paste(prefix, ".fasta", sep = ""))
  
} # end rename.trinity.transcripts 


##########################################################################
## split.fasta
##########################################################################

split.fasta <- function(fasta, num_seqs) {
  # Splits fasta file into multiple files with `num_seqs` sequences each
  #
  # Args:
  #  fasta: FASTA file
  #  num.seqs: number of sequences in each output file, except last file which
  #            will contain remainder
  #
  # Returns:
  #  Multiple fasta files with `num_seqs` sequences each

  ##--------load libraries------------
  require(Biostrings)
  require(stringr)
  
  ##--------read assembly file-----------------
  file <- readDNAStringSet(fasta)
  totalseqs <- length(file)
  #head(file)

  ##-------create subsets of length `num_seqs`--------
  # total number of subsets to create
  (num_subsets <- ceiling(totalseqs/num_seqs))
  # starting sequence
  seqs_start <- 1
  # prefix for file outputs
  prefix <- str_split_fixed(fasta, pattern=".fasta", n=2)[1]

  for(i in 1:num_subsets) {
      # if last subset, subset to remaining sequences
      if (i == num_subsets) {
          subfasta <- file[seqs_start:length(width(file)),] # remaining sequences
          # write to file
          (outname <- paste(prefix, "_sub", i, ".fasta", sep=""))
          writeXStringSet(subfasta, filepath=outname, format="fasta")
       } else { # else, extract subset of length `num_seqs`
          seqs_end <- seqs_start + (num_seqs - 1)
          message("Start:", seqs_start, " Finish:", seqs_end)
                                        # extract subset
          subfasta <- file[seqs_start:seqs_end,]
          # write to file
          (outname <- paste(prefix, "_sub", i, ".fasta", sep=""))
          writeXStringSet(subfasta, filepath=outname, format="fasta")
          # increment `seqs_start`
          seqs_start <- seqs_end + 1
      } # end ifelse
  } # end for loop
  
} # end split.fasta function

