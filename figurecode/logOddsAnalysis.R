library(ggplot2)

helixPosNames <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')
basePosNames <- c('n1', 'n2', 'n3')
nucs = c('A', 'C', 'G', 'T')
aminos = c('A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

logOddsAnalysis <- function(data) {
  
  # Returns a list of 20X20 matrices, one for each
  # (apos1, apos2, bpos, b) combination
  #chiVals <- list()
  oddsMats <- list()
  for (i in 1:(length(helixPosNames) - 1)) {
    for (j in (i+1):length(helixPosNames)) {
      for (bpos in basePosNames) {
        for (b in nucs){
          apos1 <- helixPosNames[i]
          apos2 <- helixPosNames[j]
          name <- paste(apos1, apos2, bpos, b, sep = '_')
          oddsMats[[name]] <- getOddsMat(data, apos1, apos2, bpos, b)
        }
        #name <- paste(apos1, apos2, bpos, sep = '.')
        #chiVals[[name]] <- getChiVals(data, apos1, apos2, bpos)
      }
    }
  }
  oddsMatsFrame <- makeOddsMatFrame(oddsMats)
  makeOddsHeatPlots(oddsMatsFrame)
  write.table('../../figures/logOdds/F2_low_cut10bc_0_5.txt',
              quote = FALSE, sep = '\t')
}

getOddsMat <- function(data, apos1, apos2, bpos, b) {
  
  oddsMat <- matrix(nrow = length(aminos), ncol = length(aminos))
  
  # Get the triplet count distribution and the pairwise count 
  # distributiuons using add-one smoothing to correct for zeros
  nucTab = table(data[,bpos]) + 1
  bFreq <- nucTab[b]/sum(nucTab)
  tripleTab <- table(data[,apos1], data[,apos2], data[,bpos]) + 1
  doubleTab1 <- table(data[,apos1], data[,bpos]) + 1
  doubleTab2 <- table(data[,apos2], data[,bpos]) + 1
  tsum <- sum(tripleTab)

  # Get a log odds score
  for (i in 1:length(aminos)) {
    for (j in 1:length(aminos)) {
      a1 <- aminos[i]
      a2 <- aminos[j]
      tCount <- tripleTab[a1, a2, b] 
      dCount1 <- doubleTab1[a1,b]
      dCount2 <- doubleTab2[a2,b]
      tFreq <- tCount/tsum
      dFreq1 <- dCount1/tsum
      dFreq2 <- dCount2/tsum
      oddsMat[i,j] <- log2( (tFreq/bFreq)/((dFreq1/bFreq)*(dFreq2/bFreq)) )
    }
  }
  rownames(oddsMat) <- aminos
  colnames(oddsMat) <- aminos
  oddsMat
}

makeOddsMatFrame <- function(oddsMats){
  # Convert the list of odds matrices to a dataframe
  names <- names(oddsMats)
  f <- data.frame()
  for (n in names)
    f <- rbind(f, mat2Frame(oddsMats[[n]], n))
  f
}

mat2Frame <- function(mat, name){
  # Utility function for mapping a log odds marix to 
  # a data.frame according to its name
  
  vecLen <- nrow(mat)*ncol(mat)
  amino1 <- vector(mode = 'character')
  amino2 <- vector(mode = 'character')
  score <- vector(mode = 'numeric')
  labVals <- strsplit(name, '_')
  apos1 <- rep(labVals[[1]][1], vecLen)
  apos2 <- rep(labVals[[1]][2], vecLen)
  bpos <- rep(labVals[[1]][3], vecLen)
  base <- rep(labVals[[1]][4], vecLen)
  label <- rep(name, vecLen)
  
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)) {
      amino1 <- c(amino1, rownames(mat)[i])
      amino2 <- c(amino2, colnames(mat)[j])
      score <- c(score, mat[i,j])
    }
  }
  dframe <- data.frame(amino1 = amino1, amino2 = amino2, 
             base = base, score = score, apos1 = apos1, 
             apos2 = apos2, bpos = bpos, label = label) 
}

makeOddsHeatPlots <- function(dframe) {
  # Create one faceted heat-plot for each pair 
  # of amino acid positions.
  
  for (i in 1:(length(helixPosNames) - 1)) {
    for (j in (i+1):length(helixPosNames)) {
      for (n in basePosNames) {
        print(i)
        print(helixPosNames[1])
        print(helixPosNames[2])
        h1 <- helixPosNames[i]
        h2 <- helixPosNames[j]
        print(h1)
        fname <- paste('../../figures/logOdds/',h1, '_', h2,
                       '_', n, '.pdf', sep='')
        print(fname)
        subframe <- subset(dframe, apos1 == h1 & apos2 == h2
                           & bpos == n)
        makeHeatPlot(fname, subframe)
      }
    }
  }
}

makeHeatPlot <- function(fname, dframe) {
  # Makes a faceted heatplot
  g <- ggplot(dframe, aes(amino1, amino2)) + 
    geom_tile(aes(fill = score)) +
    scale_fill_gradient2(low='firebrick', mid='white', high='steelblue',
                         midpoint = 0, na.value = "black") +
    facet_wrap(~base) +
    xlab("Amino 1") +
    ylab("Amino 2") + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0))
  ggsave(fname, plot = g, width = 7.0, height = 7.0)
}

getChiVals <- function(data, apos1, apos2, bpos) {
  
  # Returns a matrix of chi-squared scores (one for each of a1, a2 pair)
  # to evaluate the significance of change in nucleotide distribution 
  # for triplet vs pairwise interactions.
  
  chiMat <- matrix(nrow = length(aminos), ncol = length(aminos))
  
  # Get the triplet count distribution and the pairwise count 
  # distributiuons using add-one smoothing to correct for zeros
  nucTab = table(data[,bpos]) + 1
  tripleTab <- table(data[,apos1], data[,apos2], data[,bpos]) + 1
  doubleTab1 <- table(data[,apos1], data[,bpos]) + 1
  doubleTab2 <- table(data[,apos2], data[,bpos]) + 1
  tsum <- sum(tripleTab)
  
  # Get a log odds score or a chi-square value
  for (i in 1:length(aminos)) {
    for (j in 1:length(aminos)) {
      a1 <- aminos[i]
      a2 <- aminos[j]
      chiSum <- 0
      for (b in nucs) {
        tCount <- tripleTab[a1, a2, b]
        dFreq1 <- doubleTab1[a1,b]/tsum
        dFreq2 <- doubleTab2[a2,b]/tsum
        eCount = (dFreq1*dFreq2)*tsum
        chiSum <- chiSum + (tCount - eCount)^2/eCount
      }
      chiMat[i,j] <- chiSum
    }
  }
  rownames(chiMat) <- aminos
  colnames(chiMat) <- aminos
  chiMat
}

#data <- read.csv("../../data/b1hData/newDatabase/6varpos/F2/low/protein_seq_cut10bc_0_5/all.csv")
#oddsMat <- getOddsMat(data, "a0", "a6", "n3", 'T')
#chiMat <- getChiVals(data, "a0", "a6", "n3")