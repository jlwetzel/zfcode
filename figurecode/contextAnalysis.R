library(ggplot2)

# Create the list of directories
fings <- c('F1', 'F2', 'F3')
strins <- c('low', 'high')
filts <- c('filt_10e-4_025_0_c', 'filt_10e-4_05_0', 'filt_10e-4_0_5')
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
      
      # Cosine similarity
      num <- sum(v1*v2)
      den <- sqrt(sum(v1*v1)) * sqrt(sum(v2*v2))
      cosSim <- num/den
      pccs <- c(pccs, cosSim)
      
      # To use Jaccard coefficient 
      #set1 <- which(as.logical(v1))
      #set2 <- which(as.logical(v2))
      #int <- intersect(set1, set2)
      #un <- union(set1, set2)
      #jaccard <- length(int)/length(un)
      #pccs <- c(pccs, jaccard)  
      
      # To use Pearson Correlation
      #pcc <- cor(v1, v2)
      #pccs <- c(pccs, pcc)
    }
    
    pccFrame <- data.frame(prot = interProts, pcc = pccs)
    plotFile <- paste("intersectHighLow", "CosineSim",
                     f, filt, sep = '_')
    plotFile <- paste0(outDir, '/', plotFile, '.pdf')
    print(plotFile)
    g <- ggplot(pccFrame, aes(pcc)) + 
      geom_histogram() +
      xlab("Cosine similarity between high and low") +
      ylab("Number of Canonical Proteins")
    ggsave(plotFile, plot = g)
  }