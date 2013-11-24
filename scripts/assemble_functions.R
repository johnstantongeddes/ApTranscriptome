## Functions required by assembleR script


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
