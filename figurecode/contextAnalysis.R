library(ggplot2)

getCosineSim <- function(v1, v2) {
  # Returns the cosine similarity of vectors v1 and v2
  num <- sum(v1*v2)
  den <- sqrt(sum(v1*v1)) * sqrt(sum(v2*v2))
  cosSim <- num/den
  cosSim
}

getJaccard <- function(v1, v2) {
  # Returns the jaccard coefficient for two sets
  int <- intersect(set1, set2)
  un <- union(set1, set2)
  jaccard <- length(int)/length(un)
  jaccard
}

makeSimPlot <- function(simFrame, plotFile, xlabel) {

    g <- ggplot(simFrame, aes(sim)) + 
    geom_histogram() +
    xlab(xlabel) +
    ylab("Number of Canonical Proteins")

    ggsave(plotFile, plot = g)
}
  
compareHighLowString <- function(dfList, fings, filts, simType, outDir) {
  # Compare high and low stringency overlapping proteins
  # binding preferences by outputting a histogram for each 
  # of the fing/filt combinations using simType metric
  # compare the binding vectors


  for (f in fings)
    for (filt in filts) {
      lowTag <- paste(f,"low", filt, sep = '_')
      highTag <- paste(f,"high", filt, sep = '_')
      
      # Look at histogram of similarity for intersection of 
      # high vs low stringency
      interProts <- intersect(dfList[[lowTag]]$prot,
                              dfList[[highTag]]$prot)
      
      simVect <- vector()
      for (p in interProts) {
        v1Vals <- as.numeric(dfList[[lowTag]][dfList[[lowTag]]$prot == p,])
        v1 <- v1Vals[4:length(v1Vals)]
        v2Vals <- as.numeric(dfList[[highTag]][dfList[[highTag]]$prot == p,])
        v2 <- v2Vals[4:length(v2Vals)]
        
        # Cosine similarity
        if (simType == 'cosine') {
          sim <- getCosineSim(v1, v2)
          simTag <- 'CosineSim'
          xlabel <- "Cosine similarity between high and low"
        }
        else if (simType == 'cosine_bin'){
          v1 <- as.numeric(as.logical(v1))
          v2 <- as.numeric(as.logical(v2))
          sim <- getCosineSim(v1, v2)
          simTag <- 'CosineSim_binary'
          xlabel <- "Cosine similarity between high and low"
        }
        else if (simType == 'pcc') {
          sim <- cor(v1, v2)
          simTag <- 'PCC'
          xlabel <- "PCC between high and low"
        }
        else if (simType == 'jaccard'){
          set1 <- which(as.logical(v1))
          set2 <- which(as.logical(v2))
          sim <- getJaccard(v1, v2)
          simTag <- 'Jaccard'
          xlabel <- "Jaccard between high and low"
        }

        simVect <- c(simVect, sim)
      }
      
      simFrame <- data.frame(prot = interProts, sim = simVect)
      plotFile <- paste("intersectHighLow", simTag, f, 
                        filt, sep = '_')
      plotFile <- paste0(outDir, '/', plotFile, '.pdf')

      makeSimPlot(simFrame, plotFile, xlabel)
    }
}

compareF2vsF3 <- function(dfList, strins, filts, simType, outDir) {
  # Compare F2 vs F3 intersecting proteins
  # binding preferences by outputting a histogram for each 
  # of the strin/filt combinations using simType metric
  # compare the binding vectors

  for (s in strins)
    for (filt in filts) {
      F2Tag <- paste("F2",s, filt, sep = '_')
      F3Tag <- paste("F3",s, filt, sep = '_')
      
      print(F2Tag)
      #print(head(dfList[[F2Tag]]))
      print(F3Tag)
      #print(head(dfList[[F3Tag]]))
      # Look at histogram of similarity for intersection of 
      # high vs low stringency
      interProts <- intersect(dfList[[F2Tag]]$prot,
                              dfList[[F3Tag]]$prot)
      print(interProts)
      print(length(interProts))
      
      simVect <- vector()
      for (p in interProts) {
        v1Vals <- as.numeric(dfList[[F2Tag]][dfList[[F2Tag]]$prot == p,])
        v1 <- v1Vals[4:length(v1Vals)]
        v2Vals <- as.numeric(dfList[[F3Tag]][dfList[[F3Tag]]$prot == p,])
        v2 <- v2Vals[4:length(v2Vals)]
        
        # Cosine similarity
        if (simType == 'cosine') {
          sim <- getCosineSim(v1, v2)
          simTag <- 'CosineSim'
          xlabel <- "Cosine similarity between F2 and F3"
        }
        else if (simType == 'cosine_bin'){
          v1 <- as.numeric(as.logical(v1))
          v2 <- as.numeric(as.logical(v2))
          sim <- getCosineSim(v1, v2)
          simTag <- 'CosineSim_binary'
          xlabel <- "Cosine similarity between F2 and F3"
        }
        else if (simType == 'pcc') {
          sim <- cor(v1, v2)
          simTag <- 'PCC'
          xlabel <- "PCC between F2 and F3"
        }
        else if (simType == 'jaccard'){
          set1 <- which(as.logical(v1))
          set2 <- which(as.logical(v2))
          sim <- getJaccard(v1, v2)
          simTag <- 'Jaccard'
          xlabel <- "Jaccard between F2 and F3"
        }

        simVect <- c(simVect, sim)
      }
      
      simFrame <- data.frame(prot = interProts, sim = simVect)
      plotFile <- paste("F2vsF3", simTag, s, 
                        filt, sep = '_')
      plotFile <- paste0(outDir, '/', plotFile, '.pdf')

      makeSimPlot(simFrame, plotFile, xlabel)
    }
}

runHighVsLowAnalysis <- function(simType) {
  # Perform a similarity analysis between high and low
  # stringency data for each finger/filter of interest

  # Create the list of directories

  fings <- c('F1', 'F2', 'F3')
  strins <- c('low', 'high')
  filts <- c('filt_10e-4_025_0_c')
  inPref <- '../../data/b1hData/antonProcessed'
  outDir <- '../../figures/contextAnalysis/highVsLow'
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

  # Import data sets from binding vector files
  dfList <- list()
  for (i in 1:length(inFiles)) {
    dfList[[outTags[i]]] <- read.table(inFiles[i], sep = '\t',
                                       header = TRUE)
  }
  compareHighLowString(dfList, fings, filts, simType, outDir)
}

runF2vF3SimAnalysis <- function(simType) {
  # Perform a similarity analysis bewteen F2 and F3 
  # for the given combined stringency/filter

  fings <- c('F2', 'F3')
  strins <- c('union', 'inter')
  filts <- c('filt_10e-4_025_0_c')
  inPref <- '../../data/b1hData/antonProcessed'
  outDir <- '../../figures/contextAnalysis/F2vsF3'
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

  print(inFiles)

  # Import data sets from binding vector files
  dfList <- list()
  for (i in 1:length(inFiles)) {
    dfList[[outTags[i]]] <- read.table(inFiles[i], sep = '\t',
                                       header = TRUE)
  }
  print(names(dfList))
  compareF2vsF3(dfList, strins, filts, simType, outDir)
}

main <- function() {
  #runHighVsLowAnalysis('pcc')
  #runF2vF3SimAnalysis('pcc')
  runF2vF3SimAnalysis('cosine')
  #runF2vF3SimAnalysis('cosine_bin')
}

main()