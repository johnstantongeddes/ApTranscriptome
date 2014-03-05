############################################################################################
## Scripts used for ApTranscriptome_TR.Rmd
## John Stanton-Geddes
############################################################################################


############################################################################################
## read.sailfish.quant
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


############################################################################################
## RxNseq
############################################################################################

RxNseq <- function(f, model = "NA") {
    # Identify transcripts with significant reaction norms against a continuous variable
    #
    # Args:
    #  f: data.frame in long format. must contain columns specified in model formula.
    #     column with name "Transcript" must exist for proper split-apply
    #  model: model formula
    #     NOTE - for function to work correctly, each model term must be separated
    #     by a single space. for example:
    #            TPM ~ colony + val + colony:val
    #     it is important that interactions are grouped together, no spaces
    #
    # Returns:
    #  Returns dataframe with p-value and regression coefficients
    #  for each transcript
      

    # requires
    require(stringr)
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

    dd <- ddply(f, .(Transcript), function(f) {
        # fit model
        lmout <- eval(parse(text = paste("lm(", model, ", data = f)", sep = "")))
        Fvals <- anova(lmout)$'Pr(>F)'

        # assign coefficients only for significant terms, else NA
        coefvec <- vector()
        coefnames <- unlist(str_split(model, pattern = " "))
        coefnames <- paste("coef", coefnames[c(TRUE,FALSE)][-1], sep = ".")
        
        for(i in 1:(length(Fvals)-1)) { # skip last value which is residuals
            coefval <- ifelse(Fvals[i] < 0.05, Fvals[i], NA)
            coefvec <- append(coefvec, coefval)
            names(coefvec)[i] <- coefnames[i]
        }
        
        # return values
        return(c(pval = lmp(lmout),
                 coefvec))
        } # end function
    ) # end ddply
} # end function




############################################################################################
## RxNply
############################################################################################

RxNply <- function(df1) {
  # Function supplied to `ddply` to report the values at which minimum and maximum expression
  # occur, the expression level at the optimum 'temperature' and the overall shape of the 
  # expression function
  #
  # Args:
  #  df1: dataframe in long format containing expression information
  #
  # Returns:
  #  dataframe with maximum value of expression (max.val), minimum value of expression (min.val),
  #  expression at optimum (opt.exp) and the shape of the expression function (exp_type)
  
  lmout <- lm(TPM ~ val + I(val^2), data = df1)
  # NOTE - added a point to predict at 19.25C as this is the mean of the end points (0, 38.5)
  vals <- c(0, 3.5, 10, 14, 17.5, 19.25, 21, 24.5, 28, 31.5, 35, 38.5)
  newdf <- data.frame(val = vals)
  pout <- predict(lmout, newdata=newdf)
  pout <- data.frame(val = vals, exp = pout)
  
  # if all data is zero, set values to NA 
  if(coef(lmout)['(Intercept)'] == 0 & coef(lmout)['val'] == 0 & coef(lmout)['I(val^2)'] == 0) {
    max.val = NA
    min.val = NA
    opt.exp = NA
    exp_type = NA
  } else { # else set values based on predicted expression levels
    
    # get vals of max and min expression, and expression at 
    max.val <- vals[which(pout$exp == max(pout$exp))]
    min.val <- vals[which(pout$exp == min(pout$exp))]
    
    # report coefficients
    #coef(lmout)
    exp_type = if(coef(lmout)['val'] > 0 & coef(lmout)['I(val^2)'] > 0) "High" else {
      if(coef(lmout)['val'] < 0 & coef(lmout)['I(val^2)'] < 0) "Low" else {
        if(coef(lmout)['val'] > 0 & coef(lmout)['I(val^2)'] < 0) "Intermediate" else {
          "convex"} # end if
        } # end if
      } # end if
    
    # for transcripts with convex exp_type, assign to 'bimodal' if expression at both ends greater
    # than 2 SD, otherwise assign to 'low' or 'high' depending on where expression is higher
    if(exp_type == "convex") { 
      if(max(pout[pout$val <= 10, "exp"]) > 2*sd(pout$exp) &
           max(pout[pout$val >= 31.5, "exp"]) > 2*sd(pout$exp)) exp_type = "Bimodal" else {
             # linear increase?
             if(max.val > min.val) exp_type = "High" else exp_type = "Low"
           } # end if
    } # end if
  } # end else
  
  # return values
  return(c(max.val = max.val,
           min.val = min.val,
           # report expression level at optimum, 19.25C
           opt.exp = round(pout$exp[6], 4), 
           exp_type = exp_type))
} # end RxNply




############################################################################################
## geneid2GOmap
############################################################################################

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

