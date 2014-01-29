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
    file.df$trt <- trtval
    assign(outname, file.df, envir = .GlobalEnv)
}



RxNseq <- function(mat, qcrit = 0.05, makeplots = TRUE, prefix="RxN") {
    # Identify transcripts with significant reaction norms against a continuous variable
    #
    # Args:
    #  mat: matrix of expression values. first column should be transcript IDs and
    #       columns should be ordered by expression values at each assayed level
    #  vals: values that expression were assayed against. same order as columns of 'mat'
    #  qcrit: critical q-value for siginificant hits. default = 0.05
    #  makeplots: should plots be generated? default = TRUE
    #  prefix: prefix for generated files. default="RxNseq"
    #
    # Returns:
    #  1) Returns dataframe with p-value, q-value and regression coefficients
    #     for each transcript
    #  2) Writes tab-delimited text file of only significant transcripts at qcrit,
    #       sorted low to high
    #  3) Saves diagnositic plots for qvalue comparisions if makeplots = TRUE
      

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

    # plyr: split by Transcript, fit lm model to treatment values, return pvals and coefficients
    dd <- ddply(mat, .(Transcript), function(df) {
        lmout <- lm(TPM ~ trt + I(trt^2), data = df)
        return(c(pval = lmp(lmout),
                 coef(lmout)))
    })

    # False discovery rate correction with q-value
    # Remove NA pvals due to all expression values being zero. These have no expression and thus makes no sense to test for differential expression
    pvals.pos <- which(!is.na(dd$pval))
    pvals <- dd[pvals.pos,"pval"]
    
    # Calculate q-values
    qobj <- qvalue(pvals)
    #length(qvals$qvalues)
    qs <- qsummary(qobj)
    # return q summary
    assign(paste(prefix, "_qsummary", sep = ""), qs, envir = .GlobalEnv)
    
    # Add to mat.new
    dd$qval <- NA
    dd$qval[pvals.pos] <- qobj$qvalues

    # Plot qvalue diagnostics if makeplots==TRUE
    if(makeplots==TRUE) {
        message("Plotting q-value diagnostics")
        # q-value diagnostics

        pdf(paste(prefix, "_qval_diag.pdf", sep = ""))
          qplot(qobj)
        dev.off()

        pdf(paste(prefix, "_qval_hist.pdf", sep = ""))
        par(mfrow = c(2,1))
          hist(pvals)
          hist(qobj$qvalues)
        dev.off()

    }
      
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

    
