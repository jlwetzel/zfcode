library(ggplot2)

# Create the list of directories
fings <- c("F2")#c('F1', 'F2', 'F3')
strins <- c('low', 'high')
filts <- c("filt_10e-4_025_0_c")#c('filt_10e-4_025_0_c', 'filt_10e-4_05_0')
inPref <- '../../data/b1hData/antonProcessed'
outDir <- '../../figures/contextAnalysis'
inFiles <- vector()
outTags <- vector()
for (f in fings)
  for (s in strins)
    for (filt in filts) {
      inFiles <- c(inFiles, paste(inPref, f, s, filt,
                                "bindVectsUnique4pos.txt",
                                sep = '/'))
      outTags <- c(outTags, paste(f, s, filt, sep = '_'))
    }

# Import data and create a histogram of the entropy
# of the binding distribution across all of the 64 
# possible triplets
dfList <- list()
for (i in 1:length(inFiles)) {
  dfList[[outTags[i]]] <- read.table(inFiles[i], sep = '\t',
                                     header = TRUE)
#  g <- ggplot(dfList[[outTags[i]]], aes(entropy)) +
#    geom_histogram() +
#    xlab("Entropy") +
#    ylab("Number of Canonical Proteins")
#  plotFile <- paste(outDir, 
#                    paste0("entropy_", outTags[i], '.pdf'),
#                    sep = '/')
#  ggsave(plotFile, plot = g)
}

# Compare high and low stringency overlap
for (f in fings)
  for (filt in filts) {
    lowTag <- paste(f,"low", filt, sep = '_')
    highTag <- paste(f,"high", filt, sep = '_')
    
    # Look at histogram of PCC for intersection of 
    # high vs low stringency
    interProts <- intersect(dfList[[lowTag]]$prot,
                            dfList[[highTag]]$prot)
    
    pccs <- vector()
    freqs <- vector()
    freqsNegCorr <- vector()
    for (p in interProts) {
      v1Vals <- as.numeric(dfList[[lowTag]][dfList[[lowTag]]$prot == p,])
      v1 <- v1Vals[4:length(v1Vals)]
      v2Vals <- as.numeric(dfList[[highTag]][dfList[[highTag]]$prot == p,])
      v2 <- v2Vals[4:length(v2Vals)]
      freqs <- c(freqs, v1Vals[3])
      pcc <- cor(v1, v2)
      pccs <- c(pccs, cor(v1, v2))
      if (pcc < 0) {
        freqsNegCorr <- c(freqsNegCorr, v1Vals[3])
      }
    }
    
    pccFrame <- data.frame(prot = interProts, pcc = pccs)
    plotFile <- paste("intersectHighLow", "PCC",
                     f, s, filt, sep = '_')
    plotFile <- paste0(outDir, '/', plotFile, '.pdf')
    #g <- ggplot(pccFrame, aes(pcc)) + 
    #  geom_histogram() +
    #  xlab("PCC between high and low") +
    #  ylab("Number of Canonical Proteins")
    #ggsave(plotFile, plot = g)
  }