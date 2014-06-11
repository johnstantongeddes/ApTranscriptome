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
## modpFunc
############################################################################################

modpFunc <- function(foo) {
  # Function supplied to `ddply` to calculate overall model P value for a linear model
  #
  # Args:
  #  foo: data.frame model is fit to
  #
  #  'model' must be specified in environment  
  #
  # Returns:
  #  data.frame of Transcripts with pval

  lmout <- eval(parse(text = paste("lm(", model, ", data = foo)", sep = "")))
  
  # calculate overall model significance
  f <- summary(lmout)$fstatistic
  pval <- unname(pf(f[1],f[2],f[3],lower.tail=F))
  attributes(pval) <- NULL
  
  # return pvalue
  data.frame(Transcript = unique(foo$Transcript), pval = pval, adj.r.squared = summary(lmout)$adj.r.squared) 
}


############################################################################################
## modpFunc
############################################################################################

lmFunc <- function(bar) {
  # Function supplied to `dlply` to calculate overall model P value for a linear model
  #
  # Args:
  #  bar: data.frame model model is fit to
  #
  #  'model' must be specified in environment  
  #
  # Returns:
  #  list of class lm
  
  
  lmout <- eval(parse(text = paste("lm(", model, ", data = bar)", sep = "")))  
  # stepAIC to drop non-significant terms
  f.lmout <- stepAIC(lmout)
  # return final model
  return(f.lmout) 
}


############################################################################################
## grepFunc
############################################################################################

grepFunc <- function(lmitem, term = NA) {
    # Function to identify models (after stepwise AIC selection) that contain a specific term
    #
    # Arguments:
    #  lmitem: object of class `lm`
    #  term: term to grep if in model
    #
    # Returns
    #   TRUE/FALSE list of linear models that do or do not include term
    
    coefs <- coefficients(lmitem)
    coef.names <- names(coefs)
    keep <- if(length(grep(term, coef.names)) > 0) TRUE else FALSE
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

    # check class lm
    if(!"lm" %in% class(lmitem)) stop("Object not of class lm")
    
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
    if(sd(A22.pout$pTPM) < 0.00001) {
        A22.max.val <- NA
        A22.min.val <- NA
        A22.opt <- A22.pout[A22.pout$val == 19.5, "pTPM"]
        A22.exp.type <- "NotResp"
    } else {
        # else, set values from data
        A22.max.val <- median(A22.pout[which(A22.pout$pTPM == max(A22.pout$pTPM)), "val"])
        # for min.val, take median point as it may be string of points at zero
        A22.min.val <- median(A22.pout[which(A22.pout$pTPM == min(A22.pout$pTPM)), "val"])
        A22.opt <- A22.pout[A22.pout$val == 19.5, "pTPM"]

        # determine expression type
        A22.exp.type <- NA
        # "Bimodal" if max expression below 10C and above 30C is greater than one standard deviation above expression at the optimum temperature (19.5C)
        if(max(A22.pout[A22.pout$val < 10, "pTPM"])[1] > (A22.pout[A22.pout$val == 19.5, "pTPM"] + sd(A22.pout$pTPM)) & max(A22.pout[A22.pout$val > 30, "pTPM"])[1] > (A22.pout[A22.pout$val == 19.5, "pTPM"] + sd(A22.pout$pTPM))) A22.exp.type <- "Bimodal" else {
              # "PositiveLinear" if max expression at temperature above 30C and minimum expression at 0C
              if(A22.max.val > 30 & A22.min.val == 0) A22.exp.type <- "PositiveLinear" else {
                  # "HighOn" if max expression above 30C and minimum not at 0C (due to quadratic curvature)
                  if(A22.max.val > 30 & A22.min.val != 0) A22.exp.type <- "HighOn" else {
                      # "NegativeLinear" if max expression below 10C and minimum expression at 38.5
                      if(A22.max.val < 10 & A22.min.val == 38.5) A22.exp.type <- "NegativeLinear" else {
                          # "LowOn" if max expression below 10C and minimum not at 38.5C (due to quadratic curvature)
                          if(A22.max.val < 10 & A22.min.val != 38.5) A22.exp.type <- "LowOn" else {
                              A22.exp.type <- "Intermediate"
                          }
                      }
                  }
              }
          }
    }
    

        
    # calculate for Ar
    Ar.pout <- pout[pout$colony == "Ar", ]
    # if pTPM doesn't vary, set max and min values to NA
    if(sd(Ar.pout$pTPM) < 0.00001) {
        Ar.max.val <- NA
        Ar.min.val <- NA
        Ar.opt <- Ar.pout[Ar.pout$val == 19.5, "pTPM"]
        Ar.exp.type <- "NotResp"
    } else {
            # else, set values from data
            Ar.max.val <- median(Ar.pout[which(Ar.pout$pTPM == max(Ar.pout$pTPM)), "val"])
            # for min.val, take median point as it may be string of points at zero
            Ar.min.val <- median(Ar.pout[which(Ar.pout$pTPM == min(Ar.pout$pTPM)), "val"])
            Ar.opt <- Ar.pout[Ar.pout$val == 19.5, "pTPM"]

            # determine expression type
            Ar.exp.type <- NA
            # "Bimodal" if max expression below 10C and above 30C is greater than one standard deviation above expression at the optimum temperature (19.5C)
            if(max(Ar.pout[Ar.pout$val < 10, "pTPM"])[1] > (Ar.pout[Ar.pout$val == 19.5, "pTPM"] + sd(Ar.pout$pTPM)) & max(Ar.pout[Ar.pout$val > 30, "pTPM"])[1] > (Ar.pout[Ar.pout$val == 19.5, "pTPM"] + sd(Ar.pout$pTPM))) Ar.exp.type <- "Bimodal" else {
              # "PositiveLinear" if max expression at temperature above 30C and minimum expression at 0C
              if(Ar.max.val > 30 & Ar.min.val == 0) Ar.exp.type <- "PositiveLinear" else {
                  # "HighOn" if max expression above 30C and minimum not at 0C (due to quadratic curvature)
                  if(Ar.max.val > 30 & Ar.min.val != 0) Ar.exp.type <- "HighOn" else {
                      # "NegativeLinear" if max expression below 10C and minimum expression at 38.5
                      if(Ar.max.val < 10 & Ar.min.val == 38.5) Ar.exp.type <- "NegativeLinear" else {
                          # "LowOn" if max expression below 10C and minimum not at 38.5C (due to quadratic curvature)
                          if(Ar.max.val < 10 & Ar.min.val != 38.5) Ar.exp.type <- "LowOn" else {
                              Ar.exp.type <- "Intermediate"
                          }
                      }
                  }
              }
          }
        }


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

    # check class lm
    if(!"lm" %in% class(lmitem)) stop("Object not of class lm")
    
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
## RxNvarshape
############################################################################################

RxNvarshape <- function(lmitem) {
    # check class lm
    if(!"lm" %in% class(lmitem)) stop("Object not of class lm")

    # calculate mean expression across all temps
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
    A22.mean <- round(mean(A22.pout$pTPM), 4)
    Si <- diff(A22.pout$pTPM)
    A22.slope <- round(abs(sum(Si))/(nrow(A22.pout)-1), 4)
    Ci <- diff(Si)
    A22.curve <- round(abs(sum(Ci))/(nrow(A22.pout)-2), 4)
    A22.wiggle <- round((sum(abs(Ci))/(nrow(A22.pout)-2)) - A22.curve, 4)
    
    # calculate for Ar
    Ar.pout <- pout[pout$colony == "Ar", ]
    Ar.mean <- round(mean(Ar.pout$pTPM), 4)
    Si <- diff(Ar.pout$pTPM)
    Ar.slope <- round(abs(sum(Si))/(nrow(Ar.pout)-1), 4)
    Ci <- diff(Si)
    Ar.curve <- round(abs(sum(Ci))/(nrow(Ar.pout)-2), 4)
    Ar.wiggle <- round((sum(abs(Ci))/(nrow(Ar.pout)-2)) - Ar.curve, 4)
    
    
    
    # combine into data.frame
    c(A22.mean = A22.mean, A22.slope = A22.slope, A22.curve = A22.curve, A22.wiggle = A22.wiggle, Ar.mean = Ar.mean, Ar.slope = Ar.slope, Ar.curve = Ar.curve, Ar.wiggle = Ar.wiggle)
}


############################################################################################
## predFunc
############################################################################################

  predFunc <- function(foo) {
  # Function supplied to `ddply` to predict expression at each temperature
  #
  # Args:
  #  foo: data.frame model is fit to
  #
  #  'model' must be specified in environment  
  #
  # Returns:
  #  data.frame of Transcripts with pval

  lmout <- eval(parse(text = paste("lm(", model, ", data = foo)", sep = "")))
  
  # predict values
  pout <- predict(lmout)
  # back-transform from log-scale
  foo$pTPM <- exp(pout)
  
  # prediction can give values very close to zero or negative. biologically, these are meaningless so change to zero
  foo$pTPM <- ifelse(foo$pTPM < 0, 0, round(foo$pTPM, 5))

  pred.out.df <- data.frame(Transcript = foo$Transcript, TPM = foo$TPM, val = foo$val, colony = foo$colony, pTPM = foo$pTPM)
}



############################################################################################
## geneid2GOmap
############################################################################################

geneid2GOmap <- function(annotmat, ontology = c("BP", "CC", "MF"), mapname = "geneid2go.map") {
  # Create geneid2go.map file from AnnotationTable.txt produced by FastAnnotator 
  #
  # Args:
  #  annotmat: AnnotationTable.txt from FastAnnotator
  #
  # Returns:
  #  geneid2go.map saved in root directory
    
  # requires
  require(stringr)

  # check that mapname doesn't already exist
  if(file.exists(mapname)) stop(paste("File ", mapname, " already exists!", sep =""))
  
  # Extract all GO terms and combine, write to file

  
  for(r in 1:nrow(annotmat)) {
     all.GO.terms <- vector(length=0)
     if("BP" %in% ontology) {
          GO.BP.list <- str_split(annotmat[r,"GO.Biological.Process"], " ")
          GO.BP.terms <- grep('GO', unlist(GO.BP.list), value = TRUE)
          all.GO.terms <- paste(c(all.GO.terms, GO.BP.terms), collapse = ", ")
      }
      if("CC" %in% ontology) {
          GO.CC.list <- str_split(annotmat[r,"GO.Cellular.Component"], " ")
          GO.CC.terms <- grep('GO', unlist(GO.CC.list), value = TRUE)
          all.GO.terms <- paste(c(all.GO.terms, GO.CC.terms), collapse = ", ")
      }
      if("MF" %in% ontology) {
          GO.MF.list <- str_split(annotmat[r,"GO.Molecular.Function"], " ")
          GO.MF.terms <- grep('GO', unlist(GO.MF.list), value = TRUE)
          all.GO.terms <- paste(c(all.GO.terms, GO.MF.terms), collapse = ", ")
      }
        
    cat(annotmat[r, "Sequence.Name"], '\t', all.GO.terms, '\n', file = mapname, append = TRUE)    
  } # end for loop
} # end function


############################################################################################
## gsea
############################################################################################

gsea <- function(genelist, geneID2GO, plotpath=NA) {
  # Run topGO using parentChild 
  #
  # Args:
  #  genelist: list of genes with 0 for non-significant and 1 for significant genes. 'names' need to correspond to gene names in geneID2GO file
  # geneID2GO: gene read mappings created using `geneID2GOMap` function from topGO
  # plotpath: file path to save plot of significant nodes, max of 10 due to size limitations. if NA, no plot is made
  #
  # Returns:
  #  data.frame with enriched GO terms at p < 0.01
  
  # create topGOdata object
  GOdata <- new("topGOdata",
                description = "BP gene set analysis",
                ontology = "BP",
                allGenes = genelist,
                geneSel = function(genelist) {return(genelist == 1)},
                nodeSize = 10,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO)
  
  # perform enrichment analysis using parentchild method
  resultParentChild <- runTest(GOdata, statistic = 'fisher', algorithm = 'parentchild')
  
  # get number of significant GO terms p < 0.01
  numsignodes <- length(which(score(resultParentChild) < 0.01))
  
  # if 'plotpath' is set, make plot of top 10 significant nodes 
#  if(!is.na(plotpath)) {
#    # plot nodes
#    pdf(plotpath)
#    showSigOfNodes(GOdata, score(resultParentChild), firstSigNodes = min(numsignodes, 10), useInfo = 'all')
#    dev.off()
#  }  
  
  # result table
  resTable <- GenTable(GOdata, parentchild = resultParentChild, topNodes = numsignodes)
}




############################################################################################
## GSEAReportClusters
############################################################################################

GSEAReportClusters <- function(hclustobj, h) {
    # Identify members of cluster from `hclust` object
    #                                    
    # Args:
    #   hclustobj: hclust() object
    #   h: height, determined visually from looking at plot(hclustobj), to split clusters
    #
    # Returns:
    #   list of terms in each cluster
      
    clusters <- cutree(hclustobj, h = h)
    clustlist <- list()
    for(i in 1:max(clusters)) {
        li <- list(which(clusters == i))
        names(li) <- paste("Cluster", i)
        clustlist <- c(clustlist, li)
    }

    return(clustlist)
}
