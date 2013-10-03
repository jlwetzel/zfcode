library(infotheo)
library(ggplot2)

fing <- 'F1'
strin <- 'high'

#filtPrefix <- '../../data/b1hData/newDatabase/6varpos' # My files
filtPrefix <- '../../data/b1hData/antonProcessed' # Anton files
outDirPrefix <- '../../figures/mutInfo'
filtPrefix <- paste(filtPrefix, fing, strin, sep = '/')
outDirPrefix <- paste(outDirPrefix, fing, strin, sep = '/')
#filters <- c('cut10bc_0_5', 'cut10bc_025', 'cut3bc_0_5', 'cut3bc_025')
#filtDirs <- paste(filtPrefix, paste("protein_seq", filters, sep = '_'),sep = '/') # My files
filters <- c('filt_10e-4_025_0_c')
filtDirs <- paste(filtPrefix, filters,sep = '/') # Anton files

outDirs <- paste(outDirPrefix, filters, sep = '/')

runAllAnalysis <- function(filtDirs, outDirs, filters) {
  for (i in 1:length(filters)) {
    print(paste(filtDirs[i], 'all.csv', sep = '/'))
    data <- read.csv(file = paste(filtDirs[i], 'all.csv', sep = '/'))
    mutInfoAnalysis(data, outDirs[i])
  }
}

mutInfoAnalysis <- function(data, outDir) {
  # Perform analysis of mutual information for different contact 
  # positions along the helix and the dna strand
  
  helixPosNames <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')
  basePosNames <- c('b1', 'b2', 'b3')
  nucs = c('A', 'C', 'G', 'T')
  
  # Get mutual information for helix against helix
  helixMutInfo <- getNormMutInfo(data, helixPosNames, helixPosNames)
  helixMutInfo <- makeTriangular(helixMutInfo, diag = TRUE)
  helixFrame <- mat2Frame(helixMutInfo)
  writeFrame(paste(outDir, "data", "helixMutInfo.txt", sep = '/'),
             helixFrame)
  makeHeatPlot(paste(outDir, "helixMutInfo.pdf", sep = '/'),
               helixFrame)

  # Get mutual information for helix against helix in the 
  # context of specific bases in specific positions
  for (bpos in basePosNames) {
    for (b in nucs) {
      label <- paste('helixMutInfo', bpos, b, sep = '_')
      tfile <- paste(paste(paste(outDir, "data", sep = '/'),
                     label, sep = '/'), 'txt', sep = '.')
      pfile <- paste(paste(outDir,label, sep = '/'), 
                     'pdf', sep = '.')
      dsub <- data[data[[bpos]] == b,]
      print(nrow(dsub))
      helixMutInfo <- getNormMutInfo(dsub, helixPosNames, helixPosNames)
      helixMutInfo <- makeTriangular(helixMutInfo, diag = TRUE)
      helixFrame <- mat2Frame(helixMutInfo)
      writeFrame(tfile, helixFrame)
      makeHeatPlot(pfile, helixFrame)
    }
  }

  # Get mutual information for base positions against base positions
  baseMutInfo <- getNormMutInfo(data, basePosNames, basePosNames)
  baseMutInfo <- makeTriangular(baseMutInfo, diag = TRUE)
  baseFrame <- mat2Frame(baseMutInfo)
  writeFrame(paste(outDir, "data", "baseMutInfo.txt", sep = '/'),
             baseFrame)
  makeHeatPlot(paste(outDir, "baseMutInfo.pdf", sep = '/'),
               baseFrame)

  # Get mutual information for base-amino position contacts
  contactMutInfo <- getNormMutInfo(data, helixPosNames, basePosNames)
  contactFrame <- mat2Frame(contactMutInfo)
  writeFrame(paste(outDir, "data", "contactMutInfo.txt", sep = '/'),
             contactFrame)
  makeHeatPlot(paste(outDir, "contactMutInfo.pdf", sep = '/'),
               contactFrame)

  # Get mutual information for base-amino contacts in the context
  # of specific bases in specific positions
  for (bpos in basePosNames){
    for (b in nucs) {
      label <- paste('contactMutInfo', bpos, b, sep = '_')
      tfile <- paste(paste(paste(outDir, "data", sep = '/'),
                     label, sep = '/'), 'txt', sep = '.')
      pfile <- paste(paste(outDir,label, sep = '/'), 
                     'pdf', sep = '.')
      
      dsub <- data[data[[bpos]] == b,]
      contactMutInfo <- getNormMutInfo(dsub, helixPosNames, basePosNames)
      contactFrame <- mat2Frame(contactMutInfo)
      writeFrame(tfile, contactFrame)
      makeHeatPlot(pfile, contactFrame)
    }
  }
}

writeFrame <- function(dframe, fname){
  write.table(fname, dframe, sep = '\t', row.names = FALSE,
              quote = FALSE)
}

getNormMutInfo <- function(data, axis1, axis2) {
  # Takes as input a data frame (data) and two vectors
  # which are sets of names for columns in that dataframe
  # Returns a matrix giving normalized mutual information for 
  # the columns in axis1 against columns in axis2
  
  infoMat <- matrix(nrow = length(axis1), ncol = length(axis2))
  rownames(infoMat) <- axis1
  colnames(infoMat) <- axis2
  
  print(axis1)
  print(axis2)
  # Compute Shannon entropy for each separate variable
  ent1 <- vector(mode = 'numeric', length = length(axis1))
  ent2 <- vector(mode = 'numeric', length = length(axis2))
  for (i in 1:length(axis1))
    ent1[i] <- natstobits(entropy(data[axis1[i]]))
  for (i in 1:length(axis2))
    ent2[i] <- natstobits(entropy(data[axis2[i]]))
  
  # Compute the normalized mutual information
  for (i in 1:length(axis1)) {
    for (j in 1:length(axis2)) {
      condEnt <- natstobits(condentropy(data[axis1[i]], 
                                        data[axis2[j]]))
      infoMat[i, j] <- (ent1[i] - condEnt)/(min(ent1[i], ent2[j]))
    }
  }
  infoMat  
}

makeTriangular <- function(mat, diag = FALSE) {
  # Sets the upper-right elements of mat to 0
  # If diag = TRUE, also set the diagonal elements to 0
  if (diag) {
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

mat2Frame <- function(mat){
  # Utility function for mapping a log odds marix to 
  # a data.frame according to its name
  
  vecLen <- nrow(mat)*ncol(mat)
  var1 <- vector(mode = 'character')
  var2 <- vector(mode = 'character')
  score <- vector(mode = 'numeric')
  
  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)) {
      var1 <- c(var1, rownames(mat)[i])
      var2 <- c(var2, colnames(mat)[j])
      score <- c(score, mat[i,j])
    }
  }
  dframe <- data.frame(var1 = var1, var2 = var2,
                       score = score)
}

makeHeatPlot <- function(fname, dframe) {
  # Makes a faceted heatplot
  xl <- paste("")
  yl <- paste("")
  g <- ggplot(dframe, aes(var1, var2)) + 
    geom_tile(aes(fill = score)) +
    scale_fill_gradient(low = 'white', high = 'darkblue',
                        limits = c(0,1)) +
    xlab(xl) +
    ylab(yl)
  ggsave(fname, plot = g)
}

runAllAnalysis(filtDirs, outDirs, filters)