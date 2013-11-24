#################################################################
## Script to determine which transcripts show a significant norm
## of reaction across a thermal gradient
##
## John Stanton-Geddes
## 2013-10-24
#################################################################

# Settings

# Load libraries
library(plyr)
library(reshape2)

# Command-line arguments

setwd("../data/A22-merge-diginorm-oases-assembly/sailfish-out")

#read.tables <- function(file.names, ...) {
#        require(plyr)
#        ldply(file.names, function(fn) data.frame(Filename=fn, read.table(fn, sep="\t", ...)))
#    }

# Read all data files
#data <- read.tables(c("A22-00_quant/quant_bias_corrected.sf", "A22-03_quant/quant_bias_corrected.sf"))
#colnames(data) <- c("Filename", "transcript", "length", "TPM", "RPKM")

# Cast data so that each transcript in one row
# test <- head(data)
#colnames(test) <- c("Filename", "transcript", "length", "TPM", "RPKM")
#head(test)
#dcastest <- dcast(melt(test, id.vars="transcript", ID~variable + TPM))

quant00 <- read.table("A22-00_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant00) <- c("transcript", "length", "TPM00", "RPKM")
head(quant00)
dim(quant00)

quant03 <- read.table("A22-03_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant03) <- c("transcript", "length", "TPM03", "RPKM")

quant07 <- read.table("A22-07_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant07) <- c("transcript", "length", "TPM07", "RPKM")

quant10 <- read.table("A22-10_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant10) <- c("transcript", "length", "TPM10", "RPKM")

quant14 <- read.table("A22-14_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant14) <- c("transcript", "length", "TPM14", "RPKM")

quant17 <- read.table("A22-17_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant17) <- c("transcript", "length", "TPM17", "RPKM")

quant21 <- read.table("A22-21_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant21) <- c("transcript", "length", "TPM21", "RPKM")

quant24 <- read.table("A22-24_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant24) <- c("transcript", "length", "TPM24", "RPKM")

quant28 <- read.table("A22-28_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant28) <- c("transcript", "length", "TPM28", "RPKM")

quant31 <- read.table("A22-31_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant31) <- c("transcript", "length", "TPM31", "RPKM")

quant35 <- read.table("A22-35_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant35) <- c("transcript", "length", "TPM35", "RPKM")

quant38 <- read.table("A22-38_quant/quant_bias_corrected.sf", sep="\t")
colnames(quant38) <- c("transcript", "length", "TPM38", "RPKM")

head(quant38)

# Bind TPM columns to one data.frame

quantall <- cbind(quant00[,1:3], quant03$TPM03, quant07$TPM07, quant10$TPM10, quant14$TPM14, quant17$TPM17, quant21$TPM21, quant24$TPM24, quant28$TPM28, quant31$TPM31, quant35$TPM35, quant38$TPM38)
head(quantall)
dim(quantall)

# Fit quadratic regression to each transcript, retaining only transcripts with significant linear or quadratic term

temps <- c(0, 3.5, 7, 10.5, 14, 17.5, 21, 24.5, 28, 31.5, 35, 38.5)

quanttest <- quantall[1:50,]

signif <- vector(length=0)

for(i in 1:nrow(quanttest)) {
    #print(i)
    out <- lm(unlist(quantall[i,3:14]) ~ temps + I(temps^2))
    #print(anova(out))
    if(is.na(anova(out)$'Pr(>F)'[2])) cat(i, " ", "NaN", "\n") else
      if((anova(out)$'Pr(>F)'[2]) < 0.05) signif <- c(signif, i) else
        if(is.na(anova(out)$'Pr(>F)'[1])) cat(i, " ", "NaN", "\n") else
          if((anova(out)$'Pr(>F)'[1]) < 0.05) signif <- c(signif, i)
}

length(signif)

signif.transcripts <- quantall[signif, ]

dim(signif.transcripts)

signif.transcripts[1,3:14]
scale(unlist(signif.transcripts[1,3:14]))



ex.max <- max(unlist(signif.transcripts[,3:14]))
ex.min <- min(unlist(signif.transcripts[,3:14]))

pdf("expression.pdf")
  plot(temps, signif.transcripts[1,3:14], type = "l", ylim=c(ex.min, ex.max), ylab="TPM", xlab="Temp C")
  for(i in 1:nrow(signif.transcripts)) {
      lines(temps, signif.transcripts[i, 3:14], type="l")
  }
dev.off()     



