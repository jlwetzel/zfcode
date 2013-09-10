library(infotheo)
library(lattice)

mutInfoAnalysis <- function(data) {
  # Perform analysis of mutual information for different contact 
  # positions along the helix and the dna strand
  
  helixPosNames <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')
  basePosNames <- c('n1', 'n2', 'n3')
  helixMutInfo <- getMutInfo(data, helixPosNames, helixPosNames)
  helixMutInfo <- makeTriangular(helixMutInfo, diag = TRUE)
  baseMutInfo <- getMutInfo(data, basePosNames, basePosNames)
  baseMutInfo <- makeTriangular(baseMutInfo, diag = TRUE)
  baseHelixMutInfo <- getMutInfo(data, helixPosNames, basePosNames)
  x <- list(helixMutInfo, baseMutInfo, baseHelixMutInfo)
  names(x) <- c('helix', 'bases', 'contacts')
  x
}

getNormMutInfo <- function(data, axis1, axis2) {
  # Takes as input a data frame (data) and two vectors
  # which are sets of names for columns in that dataframe
  # Returns a matrix giving normalized mutual information for 
  # the columns in axis1 against columns in axis2
  
  infoMat <- matrix(nrow = length(axis1), ncol = length(axis2))
  rownames(infoMat) <- axis1
  colnames(infoMat) <- axis2
  
  # Compute Shannon entropy for each separate variable
  ent1 <- vector(mode = 'numeric', size = length(axis1))
  ent2 <- vector(mode = 'numeric', size = length(axis1))
  for (i in 1:length(axis1))
    ent[i] <- natstobits(entropy(data[axis1[i]]))
  for (i in 1:length(axis2))
    ent[i] <- natstobits(entropy(data[axis2[i]]))
  
  # Compute the normalized mutual information
  for (i in 1:length(axis1)) {
    for (j in 1:length(axis2)) {
      condEnt <- natstobits(condentropy(data[axis1[i]], 
                                        data[axis2[j]]))
      infoMat[i, j] <- (ent1 - condEnt)/(min(ent1[i], ent2[j]))
    }
  }
  infoMat  
}

makeTriangular <- function(mat, diag = FALSE) {
  # Sets the upper-right elements of mat to 0
  # If diag = TRUE, also set the diagonal elements to 0
  if (diag){
    for (i in 1:min(nrow(mat), ncol(mat)))
      mat[i,i] <- 0
  }
  for (i in 1:nrow(mat)) {
    for (j in 1:i-1) {
      mat[i,j] <- 0
    }
  }
  mat
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


