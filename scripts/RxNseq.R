############################################################################################
## Scripts used for ApTranscriptome_TR.Rmd
## John Stanton-Geddes
############################################################################################

read.sailfish.quant <- function(filein, outname, samp, trtval) {
    # Read sailfish quantification file
    #
    # Args:
    #  filein: file path to quant_bias_corrected.sf file
    #  outname: name for data.frame in R
    
    file.df <- read.table(filein, header=FALSE, sep="\t", stringsAsFactors = FALSE)
    colnames(file.df) <- c("Transcript", "Length", "TPM", "RPKM", "KPKM", "EstimatedNumReads")
    #head(file.df)
    # add columns with sample ID and treatment
    file.df$sample <- samp
    file.df$val <- trtval
    assign(outname, file.df, envir = .GlobalEnv)
}



RxNseq <- function(df, model = "NA", prefix="RxN") {
    # Identify transcripts with significant reaction norms against a continuous variable
    #
    # Args:
    #  df: data.frame in long format. must contain columns specified in model formula
    #  model: model formula, e.g. "TPM ~ colony + val"
    #  makeplots: should plots be generated? default = TRUE
    #  prefix: prefix for generated files. default="RxNseq"
    #
    # Returns:
    #  Returns dataframe with p-value and regression coefficients
    #  for each transcript
      

    # requires
    require(qvalue)
    require(plyr)
    
    #########################################################################
    # Function to report overall model P-value

    lmp <- function(modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- unname(pf(f[1],f[2],f[3],lower.tail=F))
        attributes(p) <- NULL
        return(p)
    }
    ##########################################################################

    dd <- ddply(df, .(Transcript), function(df) {
        lmout <- eval(parse(text = paste("lm(", model, ", data = df)", sep = "")))
        return(c(pval = lmp(lmout),
                 coef(lmout)))
        } # end function
    ) # end ddply

    # return this dataframe
    assign(paste(prefix, "_RxN", sep=""), dd, envir = .GlobalEnv)
    
} # end function




geneid2GOmap <- function(annotmat) {
    # Create geneid2go.map file from AnnotationTable.txt produced by FastAnnotator 
    #
    # Args:
    #  annotmat: AnnotationTable.txt from FastAnnotator
    #
    # Returns:
    #  geneid2go.map saved in root directory
      

    # requires
    require(stringr)

    # Extract all GO terms and combine, write to file
    for(r in 1:nrow(annotmat)) {
        GO.BP.list <- str_split(annotmat[r,"GO.Biological.Process"], " ")
        GO.BP.terms <- grep('GO', unlist(GO.BP.list), value = TRUE)
        GO.CC.list <- str_split(annotmat[r,"GO.Cellular.Component"], " ")
        GO.CC.terms <- grep('GO', unlist(GO.CC.list), value = TRUE)
        GO.MF.list <- str_split(annotmat[r,"GO.Molecular.Function"], " ")
        GO.MF.terms <- grep('GO', unlist(GO.MF.list), value = TRUE)


        (all.GO.terms <- paste(c(GO.BP.terms, GO.CC.terms, GO.MF.terms), collapse = ", "))
        cat(annotmat[r, "Sequence.Name"], '\t', all.GO.terms, '\n', file = "geneid2go.map", append = TRUE)    
    } # end for loop
} # end function

    

