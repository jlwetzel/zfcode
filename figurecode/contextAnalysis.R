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

makeSimHist <- function(simFrame, plotFile, xlabel) {

    g <- ggplot(simFrame, aes(sim)) + 
    geom_histogram(fill = "royalblue") +
    xlab(xlabel) +
    theme_bw() + 
    ylab("Number of Canonical Proteins")

    ggsave(plotFile, plot = g)
}

plotWeightVsSim <- function(df, plotFile, xlabel, ylabel) {
  
  g <- ggplot(df, aes(x = weight, y = sim)) +
  geom_point(color = "royalblue", size = 1) +
  xlab(xlabel) + 
  scale_x_continuous(limits = c(0, 0.01)) + 
  ylab(ylabel) +
  theme_bw()

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
      weights <- vector()
      for (p in interProts) {
        v1Vals <- as.numeric(dfList[[lowTag]][dfList[[lowTag]]$prot == p,])
        v1 <- v1Vals[4:length(v1Vals)]
        v2Vals <- as.numeric(dfList[[highTag]][dfList[[highTag]]$prot == p,])
        v2 <- v2Vals[4:length(v2Vals)]
        weight <- min(v1Vals[3], v2Vals[3])

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

        weights <- c(weights, weight)
        simVect <- c(simVect, sim)
      }
      
      #print(simVect)
      simFrame <- data.frame(prot = interProts, sim = simVect,
                             weight = weights)
      plotFile <- paste("intersectHighLow", simTag, f, 
                        filt, sep = '_')
      plotFile <- paste0(outDir, '/', plotFile, '.eps')

      makeSimHist(simFrame, plotFile, xlabel)

      plotFile <- paste("intersectHighLow", simTag, 'vWeight', f, 
                        filt, sep = '_')
      plotFile <- paste0(outDir, '/', plotFile, '.eps')
      plotWeightVsSim(simFrame, plotFile, "Protein weight",
                      xlabel)
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
      plotFile <- paste0(outDir, '/', plotFile, '.eps')

      makeSimHist(simFrame, plotFile, xlabel)
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
  # for the given combined stringency/filter ...
  # Create a histogram of similarity scores between 
  # binding vectors

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

runF2vF3DiffAnalysis1 <- function() {
  # For each DNA triplet, how many of the canonical proteins
  # that are in F2 high and F2 low are in either F3 low or 
  # F3 high ... and vice-versa

  strins <- c('inter', 'union', 'high', 'low')
  fings <- c('F3', 'F2', 'F3', 'F3') 

  filts <- c('filt_10e-4_025_0_c')
  inPref <- '../../data/b1hData/antonProcessed'
  outDir <- '../../figures/contextAnalysis/F2vsF3/diffAnalysis1'
  dirs <- paste(inPref, fings, strins, filts, sep = '/')
  bases <- c('A', 'C', 'G', 'T')
  triplets <- vector()
  for (b1 in bases)
    for (b2 in bases)
      for (b3 in bases)
        triplets <- c(triplets, paste0(b1, b2, b3))

  # Read in the two datasets to be compared
  faFrame <- read.table(paste(dirs[1], 'all_4pos.txt', sep = '/'),
                           sep = '\t', header = TRUE)
  fbFrame <- read.table(paste(dirs[2], 'all_4pos.txt', sep = '/'),
                           sep = '\t', header = TRUE)
  faLowFrame <- read.table(paste(dirs[3], 'all_4pos.txt', sep = '/'),
                           sep = '\t', header = TRUE)
  faHighFrame <- read.table(paste(dirs[4], 'all_4pos.txt', sep = '/'),
                           sep = '\t', header = TRUE)

  # Create vectors for a dataframe
  targs <- vector()
  cmpType <- vector()
  numProts1 <- vector()
  numProts2 <- vector()
  interSize <- vector()
  unionSize <- vector()

  for (t in triplets) {
    
    # Compare fa and fb
    targs <- c(targs, t)
    cmpType <- c(cmpType, "fafb")
    faProts <- faFrame[faFrame$targ==t,]$prot
    fbProts <- fbFrame[fbFrame$targ==t,]$prot
    numProts1 <- c(numProts1, length(faProts))
    numProts2 <- c(numProts2, length(fbProts))
    interSize <- c(interSize, length(intersect(faProts,
                                               fbProts)))
    unionSize <- c(unionSize, length(union(faProts,fbProts)))
    
    # Compare high and low for one of the the first set
    targs <- c(targs, t)
    cmpType <- c(cmpType, "lowhigh")
    lProts <- faLowFrame[faLowFrame$targ==t,]$prot
    hProts <- faHighFrame[faHighFrame$targ==t,]$prot
    numProts1 <- c(numProts1, length(hProts))
    numProts2 <- c(numProts2, length(lProts))
    interSize <- c(interSize, length(intersect(lProts,
                                               hProts)))
    unionSize <- c(unionSize, length(union(lProts,hProts)))
  }

  print(targs)
  print(interSize)
  print(unionSize)

  # Create the dataframe
  countFrame <- data.frame(targ = targs, cmpType = cmpType,
                           numProts1 = numProts1, numProts2 = numProts2,
                           interSize = interSize, unionSize = unionSize)


  ftag <- paste0(fings[1], strins[1], '_', fings[2], strins[2])
  tableName <- paste(outDir, (paste0(ftag,'.txt')), sep = '/')
  write.table(file = tableName, countFrame, row.names=FALSE, quote=FALSE)

  # Make the graph/s
  ctypeBreaks <- c("fafb", "lowhigh")
  ctypeLabs <- c(paste0(fings[1], " inter vs. ", fings[2], " union"), 
                 paste0(fings[1], " high vs. ", fings[1], " low"))

  # Side by side boxplot
  g <- ggplot(countFrame, aes(x = cmpType, y = interSize/unionSize,
                              fill = "royalblue")) +
    geom_boxplot(fill = "royalblue") +
    ylab("Jaccard index") +
    scale_x_discrete("", breaks = ctypeBreaks, labels = ctypeLabs) +
    theme_bw()
  plotName <- paste(outDir, (paste0(ftag,'_diff_box.pdf')), sep = '/')
  ggsave(plotName, plot = g)

  # Faceted bargraph
  g <- ggplot(countFrame, aes(x = cmpType, y = interSize/unionSize,
                              fill = cmpType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("", breaks = ctypeBreaks,
                        labels = ctypeLabs,
                        values = c("indianred", "royalblue")) +
    facet_wrap(~targ, nrow = 8, ncol = 8) +
    theme_bw() +
    ylab("Jaccard index") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
         axis.title.x = element_blank(), legend.position = "bottom",
         legend.direction = "horizontal")

  plotName <- paste(outDir, paste0(ftag,'_diff_bar.pdf'), sep = '/')
  print(plotName)
  ggsave(plotName, plot = g)
}

runWeightedFractionAnalysis <- function(varPos, refSet) {
  # Find the weighted fraction of a particular high 
  # stringency dataset contained in the other datasets
  #
  # varPos is the number of variable positions to consider
  # 4 corresponds to (-1,2,3,6), 6 to (-1,1,2,3,5,6)
  # 
  # refSet is the reference set (e.g. "F2high")



  # Read in the various datasets
  inPref <- '../../data/b1hData/antonProcessed'
  strins <- c('unionHigh', 'unionLow')#c('low', 'high')
  fings <- c('F2F3')#c('F2', 'F3')
  filt <- 'filt_10e-4_025_0_c'
  dsets <- list()
  for (f in fings)
    for (s in strins) {
      dsetLab <- paste0(f,s)
      path <- paste(inPref, f, s, filt, sep = '/')
      if (varPos == 4) {
        path <- paste(path, 'all_4pos.txt', sep = '/')
        print(path)
        dsets[[dsetLab]] <- read.table(file=path, sep = '\t',
                                       header = TRUE)
      } else {
        path <- paste(path, 'all.txt', sep = '/')
        dsets[[dsetLab]] <- read.table(path, sep = '\t',
                                       header = FALSE)
        names(dsets[[dsetLab]]) <- c('targ', 'prot', 'freq')
      }
    }

  # Set the output directory
  outDir <- '../../figures/contextAnalysis/F2vsF3/weightedFraction'

  # Get the set of triplets
  bases <- c('A', 'C', 'G', 'T')
  triplets <- vector()
  for (b1 in bases)
    for (b2 in bases)
      for (b3 in bases)
        triplets <- c(triplets, paste0(b1, b2, b3))

  # Get the fraction of the weighted fraction of the 
  # reference variables' proteins contained in each 
  # other set
  frac <- vector()
  targ <- vector()
  dsetName <- vector()
  for (n in names(dsets)) {
    if (n != refSet) 
      for (t in triplets) {
        
        tsetRef <- subset(dsets[[refSet]], targ == t)
        tset <- subset(dsets[[n]], targ == t)
        freqsRef <- tsetRef$freq
        names(freqsRef) <- as.character(tsetRef$prot)
        
        seqOlap <- intersect(as.character(tset$prot),
                             names(freqsRef))
        fracCovered <- 0.0
        for (s in seqOlap) {
          fracCovered <- fracCovered + freqsRef[s]
        }
      frac <- c(frac, fracCovered)
      dsetName <- c(dsetName, n)
      targ <- c(targ, t)
      }
    }

  fracFrame <- data.frame(targ = targ, frac = frac, 
                          dsetName = dsetName)

  #print(fracFrame)

  if (refSet == 'F2high') {
    fracFrame$dsetName <- factor(fracFrame$dsetName, 
                                 levels = c("F2low", "F3low", "F3high"),
                                 labels = c("F2 low", "F3 low", "F3 high"))
    llabs <- c("F2 low", "F3 high", "F3 low")
    ylabs <- paste0("Weighted fraction of F2 high stringency found")
  } else if (refSet == 'F3high') {
    fracFrame$dsetName <- factor(fracFrame$dsetName, 
                                 levels = c("F3low", "F2low", "F2high"),
                                 labels = c("F3 low", "F2 low", "F2 high"))
    llabs <- c("F3 low", "F2 high", "F2 low")
    ylabs <- paste0("Weighted fraction of F3 high stringency found")
  } else if (refSet == 'F2F3unionHigh') {
    fracFrame$dsetName <- factor(fracFrame$dsetName, 
                                 levels = c("F2F3unionLow"),
                                 labels = c("F2F3 Union Low"))
    llabs <- c("F2F3 Union Low")
    ylabs <- paste0("Weighted fraction of high stringency found in low stringency")
  }

  #cols <- c("indianred", "royalblue", "gray50")

  # Make the plot
  g <- ggplot(fracFrame, aes(x = dsetName, y = frac)) +
    #scale_fill_manual("", breaks = llabs, values = cols) + 
    geom_boxplot(aes(outlier.size = 1), fill = "royalblue") +
    ylab(ylabs) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=9),
          axis.text.y = element_text(size=9),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  plotName <- paste(outDir, 
                    paste0(refSet, '_', 
                           varPos, 'pos_weightedFrac.eps'),
                    sep = '/')
  print(plotName)
  tableName <- paste(outDir, 
                    paste0(refSet, '_', 
                           varPos, 'pos_weightedFrac.txt'),
                    sep = '/')
  write.table(file = tableName, fracFrame, row.names=FALSE)
  ggsave(plotName, plot = g, width = 4, height = 5.5)

}

parseWeightedJaccard <- function(dframe, triplets) {
  # Compute the weighted Jaccard coefficient 
  # of binding canonical protein for all 
  # pairs of DNA triplets in the dframe

  # Create a jaccard dataframe
  trip1 <- vector()
  trip2 <- vector()
  wjac <- vector()
  for (i in 1:length(triplets)) {
    for (j in 1:length(triplets)) {
      
      # Get the weighted jaccard similarity for this
      # pair of triplets ( (sum of intersection)/2)
      t1 <- triplets[i]
      t2 <- triplets[j]
      if (i <= j){
        tframe1 <- subset(dframe, targ == t1)
        tframe2 <- subset(dframe, targ == t2)
        f1 <- tframe1$freq
        f2 <- tframe2$freq
        names(f1) <- tframe1$prot
        names(f2) <- tframe2$prot
        interProts <- intersect(names(f1), names(f2))
        unionProts <- union(names(f1), names(f2))
        interWeight <- 0.0
      for (p in interProts)
        interWeight <- interWeight + f1[p] + f2[p]
      wjaccard <- interWeight/2
      
      # Correct for rounding errors
      if (wjaccard > 1)
        wjaccard = 1
      
      } else {
        wjaccard <- 0
      }
      
      # Add new entry to each vector for new dframe
      trip1 <- c(trip1, t1)
      trip2 <- c(trip2, t2)
      wjac <- c(wjac, wjaccard)
    }
  }
  jacFrame <- data.frame(trip1 = trip1, trip2 = trip2,
                         wjac = wjac)
  jacFrame
}

makeTripletHeatmap <- function(fing, strin, filt,
                               noParse = FALSE) {
  
  # Get the set of triplets
  bases <- c('A', 'C', 'G', 'T')
  triplets <- vector()
  for (b1 in bases)
    for (b2 in bases)
      for (b3 in bases)
        triplets <- c(triplets, paste0(b1, b2, b3))

  # Set the output directory
  outDir <- '../../figures/contextAnalysis/tripletHeatMap'
  inPref <- '../../data/b1hData/antonProcessed'
  
  # Compute the coefficients if not already done
  if (noParse) {
    path <- paste(outDir, paste(fing, strin, 
                                'wJaccTrip.txt', sep = '_'),
                  sep = '/')
    jacFrame <- read.table(file=path, header = TRUE)
    print(nrow(jacFrame))
  } else {
    path <- paste(inPref, fing, strin, filt,
                  'all_4pos.txt', sep = '/')
    dframe <- read.table(file=path, sep = '\t',
                       header = TRUE)
    jacFrame <- parseWeightedJaccard(dframe, triplets)
    tableName <- paste(outDir, paste(fing, strin, 
                                   'wJaccTrip.txt', sep = '_'),
                        sep = '/')
    write.table(file = tableName, jacFrame, row.names=FALSE, 
                quote=FALSE)
    print(nrow(jacFrame))
  } 

  jacFrame$trip2 <- factor(jacFrame$trip2, levels = rev(triplets))

  # Plot the heatmap and save to file
  br <- seq(0,1, 0.05)
  print(names(jacFrame))
  g <- ggplot(jacFrame, aes(x = trip1, y = trip2)) + 
    geom_tile(aes(fill = wjac)) +
    #scale_fill_gradient2(breaks = br,
    #                     low = 'white', high = 'royalblue',
    #                     limits = c(0,1), guide = 'legend') +
    scale_fill_gradient2("Weighted Jaccard", low = 'white', high = 'royalblue',
                         limits = c(0,1))+#, guide = "legend") +
    geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = 64.5)) + 
    geom_segment(aes(x = 16.5, xend = 16.5, y = 0.5, yend = 65 - 16.5)) + 
    geom_segment(aes(x = 32.5, xend = 32.5, y = 0.5, yend = 65 - 32.5)) + 
    geom_segment(aes(x = 48.5, xend = 48.5, y = 0.5, yend = 65 - 48.5)) +
    geom_segment(aes(y = 0.5, yend = 0.5, x = 0.5, xend = 64.5)) +  
    geom_segment(aes(y = 16.5, yend = 16.5, x = 0.5, xend = 65 - 16.5)) + 
    geom_segment(aes(y = 32.5, yend = 32.5, x = 0.5, xend = 65 - 32.5)) + 
    geom_segment(aes(y = 48.5, yend = 48.5, x = 0.5, xend = 65 - 48.5)) + 
    #geom_segment(aes(y = 49.5, yend = 49.5, x = 0.5, xend = 65 - 49.5),
    #             size = 0.1, color = "gray80") + 
    
    theme_bw() + 

    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(),
          legend.position = c(0.65, 0.65),
          legend.background = element_rect(fill="white", color = "black",
                                           size=.5)) +

    #draws x and y axis line
    theme(axis.line = element_line(color = 'black', size = 0.5)) + 

    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle=90, vjust=0.5, size = 8))

  plotName <- paste(outDir, paste(fing, strin, 
                                   'wJaccTrip.eps', sep = '_'),
                        sep = '/')
  ggsave(plotName, plot = g, width = 7.75, height = 7.5)
}

main <- function() {
  #runHighVsLowAnalysis('cosine')
  #runF2vF3SimAnalysis('pcc')
  #runF2vF3SimAnalysis('cosine')
  #runF2vF3SimAnalysis('cosine_bin')
  #runWeightedFractionAnalysis(6, "F2F3unionHigh")
  #runWeightedFractionAnalysis(4, "F2F3unionHigh")
  makeTripletHeatmap("F3", "union", 
                     'filt_10e-4_025_0_c', noParse = TRUE)
}

main()