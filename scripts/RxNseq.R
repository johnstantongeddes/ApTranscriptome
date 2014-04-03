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
## RxNtype
############################################################################################

RxNtype <- function(lmitem) {
    # Function passed to `ldply` that defines transcripts into response type
    #
    # Arguments:
    #  lmitem: object of class `lm`
    #
    # Returns
    #   dataframe with values of max expression, min expression and expression type for each colony
    
  # predict 
  pred.vals <- seq(from = 0, to = 38.5, by = 0.5)
  colony.levels <- as.factor(rep(unlist(lmitem$xlevels)))
  # for transcripts with no 'colony' effect in model, no levels, so create dummy
  if(length(colony.levels) == 0) colony.levels <- c("A22", "Ar")
  newdf <- data.frame(val = pred.vals, colony = rep(colony.levels, each = length(pred.vals)))
  pout <- predict(lmitem, newdata=newdf)
  pout <- cbind(newdf, logTPM = pout)

  # back-transform from log-scale
  pout$pTPM <- exp(pout$logTPM)
  
  # prediction can give values very close to zero or negative. biologically, these are meaningless so change to zero
  pout$pTPM <- ifelse(pout$pTPM < 0, 0, round(pout$pTPM, 5))
  
  # calculate for A22
  A22.pout <- pout[pout$colony == "A22", ]
  # if pTPM doesn't vary, set max and min values to NA
  if(sd(A22.pout$pTPM) < 0.01) {
      A22.max.val <- NA
      A22.min.val <- NA
      A22.opt <- A22.pout[A22.pout$val == 19.5, "pTPM"]
      A22.exp.type <- "NotExp"
  } else {
          # else, set values from data
          A22.max.val <- median(A22.pout[which(A22.pout$pTPM == max(A22.pout$pTPM)), "val"])
          # for min.val, take median point as it may be string of points at zero
          A22.min.val <- median(A22.pout[which(A22.pout$pTPM == min(A22.pout$pTPM)), "val"])
          A22.opt <- A22.pout[A22.pout$val == 19.5, "pTPM"]
          # determine expression shape
          A22.exp.type <- NA
          # expression type defined as "Bimodal" if max expression below 10C and above 30C is greater than one standard deviation above expression at the optimum temperature (19.5C)
          if(max(A22.pout[A22.pout$val < 10, "pTPM"])[1] > (A22.pout[A22.pout$val == 19.5, "pTPM"] + sd(A22.pout$pTPM)) & max(A22.pout[A22.pout$val > 30, "pTPM"])[1] > (A22.pout[A22.pout$val == 19.5, "pTPM"] + sd(A22.pout$pTPM))) A22.exp.type <- "Bimodal" else {
              # expression type "High" if max expression at temperature above 30C
              if(A22.max.val > 30) A22.exp.type <- "High" else {
              # expression type "Low" if max expression at temperature below 10C
                  if(A22.max.val < 10) A22.exp.type <- "Low" else
                  A22.exp.type <- "Intermediate"
              } # end else
          } # end else
      } # end if/else

  # calculate for Ar
  Ar.pout <- pout[pout$colony == "Ar", ]
  # if pTPM doesn't vary, set max and min values to NA
  if(sd(Ar.pout$pTPM) < 0.01) {
      Ar.max.val <- NA
      Ar.min.val <- NA
      Ar.opt <- Ar.pout[Ar.pout$val == 19.5, "pTPM"]
      Ar.exp.type <- "NotExp"
  } else {
          # else, set values from data
          Ar.max.val <- median(Ar.pout[which(Ar.pout$pTPM == max(Ar.pout$pTPM)), "val"])
          # for min.val, take median point as it may be string of points at zero
          Ar.min.val <- median(Ar.pout[which(Ar.pout$pTPM == min(Ar.pout$pTPM)), "val"])
          Ar.opt <- Ar.pout[Ar.pout$val == 19.5, "pTPM"]
          # determine expression shape
          Ar.exp.type <- NA
          # expression type defined as "Bimodal" if max expression below 10C and above 30C is greater than one standard deviation above expression at the optimum temperature (19.5C)
          if(max(Ar.pout[Ar.pout$val < 10, "pTPM"])[1] > (Ar.pout[Ar.pout$val == 19.5, "pTPM"] + sd(Ar.pout$pTPM)) & max(Ar.pout[Ar.pout$val > 30, "pTPM"])[1] > (Ar.pout[Ar.pout$val == 19.5, "pTPM"] + sd(Ar.pout$pTPM))) Ar.exp.type <- "Bimodal" else {
              # expression type "High" if max expression at temperature above 30C
              if(Ar.max.val > 30) Ar.exp.type <- "High" else {
              # expression type "Low" if max expression at temperature below 10C
                  if(Ar.max.val < 10) Ar.exp.type <- "Low" else
                  Ar.exp.type <- "Intermediate"
              } # end else
          } # end else
      } # end if/else

  data.frame(A22.max = A22.max.val,
           A22.min = A22.min.val,
           A22.opt = A22.opt,
           A22.type = A22.exp.type,
           Ar.max = Ar.max.val,
           Ar.min = Ar.min.val,
           Ar.opt = Ar.opt,
           Ar.type = Ar.exp.type)
           
} # end RxNtype


############################################################################################
## transcriptSD
############################################################################################

transcriptSD <- function(lmitem, colony) {
    # function passed to Map() that calculates the standard deviation of expression for each 'Intermediate' transcript
    #
    # Arguments:
    #  lmitem: List item of class `lm`
    #  colony: Colony that standard deviation of expression is being calculated for
    #
    # Returns:
    #   list of standard deviation of expression for each transcript
    
    # predict only for A22
    pred.vals <- seq(from = 0, to = 38.5, by = 0.5)
    newdf <- data.frame(val = pred.vals, colony = rep(colony, length = length(pred.vals)))
    pout <- predict(lmitem, newdata=newdf)
    pout <- cbind(newdf, pTPM = pout) # note this is on log scale

    # prediction can give values very close to zero or negative. biologically, these are meaningless so change to zero
    pout$pTPM <- ifelse(pout$pTPM < 0, 0, round(pout$pTPM, 5))

    # draw 1000 'temps' weighted by predicted expression
    random.draw <- sample(size = 1000, pred.vals, prob = pout$pTPM, replace = TRUE)

    # calculate standard deviation of weighted random.draw
    data.frame(exp_sd = sd(random.draw))
    } # end function



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

