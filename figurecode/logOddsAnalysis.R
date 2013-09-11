helixPosNames <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')
basePosNames <- c('n1', 'n2', 'n3')
nucs = c('A', 'C', 'G', 'T')
aminos = c('A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

logOddsAnalysis <- function(data) {
  
  # Returns a list of 20X20 matrices, one for each
  # (apos1, apos2, bpos, b) combination
  chiVals <- list()
  oddsMats <- list()
  for (i in 1:(length(helixPosNames) - 1)) {
    for (j in (i+1):length(helixPosNames)) {
      for (bpos in basePosNames) {
        for (b in nucs){
          apos1 <- helixPosNames[i]
          apos2 <- helixPosNames[j]
          name <- paste(apos1, apos2, bpos, b, sep = '.')
          oddsMats[[name]] <- getOddsMat(data, apos1, apos2, bpos, b)
        }
        name <- paste(apos1, apos2, bpos, sep = '.')
        chiVals[[name]] <- getChiVals(data, apos1, apos2, bpos)
      }
    }
  }
  oddsMats
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

  # Get a log odds score or a chi-square value
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

makeLevelPlot <- function(fname, mat) {
  # Makes a levelplot for the matrix values
  #bluered = colorRampPalette(c("blue", "red"),
  #                           space = "Lab") 
  pdf(fname, height = 4.55, width = 2.84)
  #print(levelplot(t(mat), col.regions = bluered,
  #                xlab = "", ylab = ""))
  print(levelplot(t(mat), xlab = "", ylab = ""))
  dev.off()
}

#data <- read.csv("../../data/b1hData/newDatabase/6varpos/F2/low/protein_seq_cut10bc_0_5/all.csv")
#oddsMat <- getOddsMat(data, "a0", "a6", "n3", 'T')
#chiMat <- getChiVals(data, "a0", "a6", "n3")