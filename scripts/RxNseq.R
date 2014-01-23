############################################################################################
## Scripts used for ApTranscriptome_TR.Rmd
## John Stanton-Geddes
############################################################################################

read.sailfish.quant <- function(filein, outname) {
    # Read sailfish quantification file
    #
    # Args:
    #  filein: file path to quant_bias_corrected.sf file
    #  outname: name for data.frame in R
    
    file.df <- read.table(filein, header=FALSE, sep="\t", stringsAsFactors = FALSE)
    colnames(file.df) <- c("Transcript", "Length", "TPM", "RPKM", "KPKM", "EstimatedNumReads")
    #head(file.df)
    assign(outname, file.df, envir = .GlobalEnv)
}



RxNseq <- function(mat, vals, qcrit = 0.05, makeplots = TRUE, prefix="RxN") {
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
    #  2) Writes ab-delimited text file of only significant transcripts at qcrit,
    #       sorted low to high
    #  3) Saves diagnositic plots for qvalue comparisions if makeplots = TRUE
      

    # requires
    require(qvalue)

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

    # Check that length(vals)==ncols(mat)+1
    if(length(vals)+1 != ncol(mat)) stop("Number of martrix columns does not equal number of X values!")

    # Fit quadratic regression to each transcript and record overall model P-value
    # for significant models, report linear and quadratic terms

    # Create new matrix to store results
    # Set columns for p-val, linear and quadratic terms to NA
    mat.new <- mat
    mat.new$pval <- NA
    mat.new$intercept <- NA
    mat.new$lin.coef <- NA
    mat.new$quad.coef <- NA
    
    # loop across transcripts
    for(i in 1:nrow(mat)) {
        if(i%%1000 == 0) message(i, " transcripts examined")
        # fit model - note that first column of transcript IDs is left out
        out <- lm(unlist(mat[i,2:ncol(mat)]) ~ vals + I(vals^2))
        lmp.out <- lmp(out)
        anova.out <- anova(out)
        #pvals <- c(pvals, lmp(out))
        if(is.na(lmp.out)) {
            next } else {
                # record p-value
                mat.new$pval[i] <- lmp.out
                # check each parameter. if significant, record
                # intercept
                if(lmp.out < 0.05) mat.new$intercept[i] <- out$coefficients["(Intercept)"]
                # linear coefficient
                if(anova.out$'Pr(>F)'[1] < 0.05) mat.new$lin.coef[i] <- out$coefficients["vals"]
                # quadratic coefficient
                if(anova.out$'Pr(>F)'[2] < 0.05) mat.new$quad.coef[i] <- out$coefficients["I(vals^2)"]
        } # end if statement
    } # end for loop

    # False discovery rate correction with q-value
    # Remove NA pvals due to all expression values being zero. These have no expression and thus makes no sense to test for differential expression
    pvals.pos <- which(!is.na(mat.new$pval))
    pvals <- mat.new[pvals.pos,"pval"]
    
    # Calculate q-values
    qobj <- qvalue(pvals)
    #length(qvals$qvalues)
    qs <- qsummary(qobj)
    # return q summary
    assign(paste(prefix, "_qsummary", sep = ""), qs, envir = .GlobalEnv)
    
    # Add to mat.new
    mat.new$qval <- NA
    mat.new$qval[pvals.pos] <- qobj$qvalues

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
    assign(paste(prefix, "_RxN", sep=""), mat.new, envir = .GlobalEnv)
    
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

    
