
RxNseq <- function(mat, vals, qcrit = 0.05, makeplots = TRUE, prefix="RxN") {
    # Identify transcripts with significant reaction norms against a continuous variable
    #
    # Args:
    #  mat: matrix of expression values. rownames should be transcript IDs and
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

    # Check that length(vals)==ncols(mat)
    if(length(vals) != ncol(mat)) stop("Number of martrix columns does not equal number of X values!")

    # Fit quadratic regression to each transcript and record overall model P-value
    # for significant models, report linear and quadratic terms

    #signif <- vector(length=0)
    #pvals <- vector(length=0)
    # Set columns for p-val, linear and quadratic terms to NA
    mat.new <- mat
    mat.new$pval <- NA
    mat.new$lin.coef <- NA
    mat.new$quad.coef <- NA
    
    #for(i in 1:nrow(mat)) {
    for(i in 1:nrow(mat)) {
        if(i%%1000 == 0) message(i, " transcripts examined")
        out <- lm(unlist(mat[i,]) ~ vals + I(vals^2))
        #pvals <- c(pvals, lmp(out))
        if(is.na(lmp(out))) {
            next } else {
            # check if linear term, significant. if so, record
                # pvalue
                mat.new$pval[i] <- lmp(out)
                # linear coefficient
                if(anova(out)$'Pr(>F)'[1] < 0.05) mat.new$lin.coef[i] <- out$coefficients["vals"]
                # quadratic coefficient
                if(anova(out)$'Pr(>F)'[2] < 0.05) mat.new$quad.coef[i] <- out$coefficients["I(vals^2)"]
        }
    }

    # False discovery rate correction with q-value
    # Remove NA pvals due to all expression values being zero. These have no expression and thus makes no sense to test for differential expression
    pvals.pos <- which(!is.na(mat.new$pval))
    pvals <- mat.new[pvals.pos,"pval"]
    
    # Calculate q-values
    qobj <- qvalue(pvals)
    #length(qvals$qvalues)
    qsummary(qobj)

    # Add to mat.new
    mat.new$qval <- NA
    mat.new$qval[pvals.pos] <- qobj$qvalues

    # Plot qvalue diagnostics if makeplots==TRUE
    if(makeplots==TRUE) {
        message("Plotting q-value diagnostics")
        # q-value diagnostics

        pdf(paste(prefix, "_qval_diag.pdf", sep = ""))
          qplot(qvals)
        dev.off()

        pdf(paste(prefix, "_qval_hist.pdf", sep = ""))
        par(mfrow = c(2,1))
          hist(pvals)
          hist(qvals$qvalues)
        dev.off()
    }
    
    # Order by p-value
    responsive.transcripts.by.p <- mat.new[order(mat.new$pval), ]
    
    # Return this dataframe
    assign(paste(prefix, "_TPM_stats", sep=""), responsive.transcripts.by.p, envir = .GlobalEnv)
    
    # Dataframe with significant hits, FDR < qcrt
    responsive.transcripts.FDR <- responsive.transcripts.by.p[which(responsive.transcripts.by.p$qval < qcrit),]

    # Write only significant results to file
    write.table(file=paste(prefix, "_FDR_responsive_transcripts.txt", sep = ""), responsive.transcripts.FDR, row.names=FALSE, quote=FALSE, sep="\t")

} # end function
