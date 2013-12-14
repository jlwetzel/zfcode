library(infotheo)
library(ggplot2)
library(grid)
library(scales)

fing <- 'F2F3'
strin <- 'unionUnions'

#filtPrefix <- '../../data/b1hData/newDatabase/6varpos' # My files
filtPrefix <- '../../data/b1hData/antonProcessed' # Anton files
outDirPrefix <- '../../figures/mutInfo'
filtPrefix <- paste(filtPrefix, fing, strin, sep = '/')
outDirPrefix <- paste(outDirPrefix, fing, strin, sep = '/')
#filters <- c('cut10bc_0_5', 'cut10bc_025', 'cut3bc_0_5', 'cut3bc_025')
#filtDirs <- paste(filtPrefix, paste("protein_seq", filters, sep = '_'),sep = '/') # My files
filters <- c('filt_10e-4_025_0_c')#, 'filt_10e-4_05_0', 'filt_10e-4_0_5')
filtDirs <- paste(filtPrefix, filters,sep = '/') # Anton files

outDirs <- paste(outDirPrefix, filters, sep = '/')

runAllAnalysis <- function(filtDirs, outDirs, filters) {
  for (i in 1:length(filters)) {
    print(paste(filtDirs[i], 'all.csv', sep = '/'))
    data <- read.csv(file = paste(filtDirs[i], 'all.csv', sep = '/'))
    mutInfoAnalysis(data, outDirs[i])
  }
}

mats2Vectors <- function(mats, numTimes){
  # Converts each position of a list of x matrices into 
  # a corresponding vector of x values

  vectList <- list()
  for (i in rownames(mats[[1]])) {
    for (j in colnames(mats[[1]])) {
      vectList[[paste(i,j)]] <- vector(length = numTimes)
    }
  }

  for (k in 1:length(mats)){
    for (i in rownames(mats[[1]])) {
      for (j in colnames(mats[[1]])) {
        vectList[[paste(i,j)]][k] <- mats[[k]][i,j]
      }
    }
  }
  vectList
}

makeStatMat <- function(vectList, rnames, cnames, 
                        stat = "mean"){
  # Return a matrix of a summary statistic 
  # for each contact pair across the multiple 
  # randomization experiments

  statMat <- matrix(nrow = 6, ncol = 3)
  rownames(statMat) <- rnames
  colnames(statMat) <- cnames
  for (i in rnames) {
    for (j in cnames) {
      if (stat == "mean"){
        statMat[i,j] <- mean(vectList[[paste(i,j)]])
      } else if (stat == "sd"){
        statMat[i,j] <- sd(vectList[[paste(i,j)]])
      }
    }
  }
  statMat
}

makeEmpPvalCountMat <- function(vectList, rnames,
                                cnames, actualMat) {
  # Count how many times the values for each contact 
  # pair in the randomized experiment are greater than 
  # or equal to the actual values.

  countMat <- matrix(nrow = 6, ncol = 3)
  rownames(countMat) <- rnames
  colnames(countMat) <- cnames

  for (k in length(vectList)){
    for (i in rnames) {
      for (j in cnames) {
        vect <- vectList[[paste(i,j)]]
        #print(paste(i,j))
        #print(vect)
        #print(actualMat[i,j])
        #print(length(vect[vect >= actualMat[i,j]]))
        countMat[i,j] <- length(vect[vect >= actualMat[i,j]])
      }
    }
  }
  countMat
}

randAnalysis <- function(data, helixPosNames, basePosNames, 
                         outDir, actualMat, numTimes = 100) {
  
  # Create numTimes randomized MI matrices by permuting
  # the DNA triplet labels with respect to the helices

  set.seed(972364)
  randMats <- list()
  for (t in 1:numTimes) {
    # Create a randomized dataframe where labels of triplets
    # are shuffled with respect to proteins
    print(t)
    randData <- data
    randOrder <- sample(seq(1,length(data$b1)))
    randData$b1 <- randData$b1[randOrder]
    randData$b2 <- randData$b2[randOrder]
    randData$b3 <- randData$b3[randOrder]
    #write.csv(randData,
    #          paste(outDir, "data", "shuffleTriples.csv", sep = '/'))
    contactMutInfo <- getNormMutInfo(randData, helixPosNames, basePosNames)
    randMats[[t]] <- contactMutInfo
  }

  # Convert the set of randomized matrices to 
  # a list of vectors (one vector per contact pair)
  # and write to file
  #print(randMats)
  vectList <- mats2Vectors(randMats, numTimes)
  vectFrame <- data.frame(matrix(unlist(vectList), 
                          ncol = length(vectList)))
  names(vectFrame) <- names(vectList)
  writeFrame(paste(outDir, "data", 
             paste0("rand_contactMutInfoRawNums_", numTimes ,".txt"),
             sep = '/'), vectFrame)
  #print(vectList)

  # Get a matrix of means
  avgMat <- makeStatMat(vectList, helixPosNames,
                        basePosNames, stat = "mean")
  
  # Get a matrix of standard deviation
  sdMat <- makeStatMat(vectList, helixPosNames,
                        basePosNames, stat = "sd")
  # Get a matrix of counts for empirical p-vals
  pvalCountMat <- makeEmpPvalCountMat(vectList, helixPosNames,
                                      basePosNames, actualMat)
  print(pvalCountMat)
  # Convert each matrix to a dataframe
  avgFrame <- mat2Frame(avgMat)
  sdFrame <- mat2Frame(sdMat)
  pvalCountFrame <- mat2Frame(pvalCountMat)

  # Combine the mats into a single dataframe
  statFrame <- avgFrame[c(1,2)]
  statFrame$nrand <- rep(numTimes, nrow(statFrame))
  statFrame$mean <- avgFrame$score
  statFrame$sd <- sdFrame$score
  statFrame$count <- pvalCountFrame$score
  statFrame$raw.pval <- statFrame$count/numTimes
  statFrame$bonf.pval <- statFrame$raw.pval*nrow(statFrame)*ncol(statFrame)

  # Write the dataframe into a text file
  writeFrame(paste(outDir, "data", 
             paste0("rand_contactMutInfoStats_", numTimes ,".txt"),
             sep = '/'), statFrame)
  #makeHeatPlot(paste(outDir, "rand_contactMutInfo.eps", sep = '/'),
  #             contactFrame, "Helix position", "Base position")
}

mutInfoAnalysis <- function(data, outDir) {
  # Perform analysis of mutual information for different contact 
  # positions along the helix and the dna strand
  
  helixPosNames <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')
  basePosNames <- c('b1', 'b2', 'b3')
  nucs = c('A', 'C', 'G', 'T')
  
  if(FALSE) {
  # Get mutual information for helix against helix
  helixMutInfo <- getNormMutInfo(data, helixPosNames, helixPosNames)
  helixMutInfo <- makeTriangular(helixMutInfo, diag = TRUE)
  helixFrame <- mat2Frame(helixMutInfo)
  writeFrame(paste(outDir, "data", "helixMutInfo.txt", sep = '/'),
             helixFrame)
  makeHeatPlot(paste(outDir, "helixMutInfo.eps", sep = '/'),
               helixFrame, "Helix position", "Helix position")

  
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
  makeHeatPlot(paste(outDir, "baseMutInfo.eps", sep = '/'),
               baseFrame, "Base position", "Base position")
  } #END IF FALSE

  # Get mutual information for base-amino position contacts
  contactMutInfo <- getNormMutInfo(data, helixPosNames, basePosNames)
  contactFrame <- mat2Frame(contactMutInfo)
  writeFrame(paste(outDir, "data", "contactMutInfo.txt", sep = '/'),
             contactFrame)
  makeHeatPlot(paste(outDir, "contactMutInfo.eps", sep = '/'),
               contactFrame, "Helix position", "Base position")
  
  # Get mutual information for contacts with randomly 
  # shuffled triples
  #print(system.time(randAnalysis(data, helixPosNames, 
  #                               basePosNames, outDir, 
  #                               contactMutInfo, 10000)))

  if (FALSE){
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
      makeHeatPlot(pfile, contactFrame, "Helix position",
                   "Base position")
    }
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
  
  #print(axis1)
  #print(axis2)
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
  # Utility function for mapping a log odds matrix to 
  # a data.frame according to its name
  
  hpos <- vector()
  bpos <- vector()
  score <- vector()

  xIsHelix <- FALSE
  yIsHelix <- FALSE

  for (i in 1:nrow(mat)){
    for (j in 1:ncol(mat)) {
      if (rownames(mat)[i] == "a0") {
        yIsHelix <- TRUE
        hpos <- c(hpos, "-1")
      } else {
        hpos <- c(hpos, substr(rownames(mat)[i], 2, 2)) 
      }
      if (colnames(mat)[j] == "a0"){
        xIsHelix <- TRUE
        bpos <- c(bpos, "-1")
      } else { 
        bpos <- c(bpos, substr(colnames(mat)[j], 2, 2))
      }
      score <- c(score, mat[i,j])
    }
  }
  dframe <- data.frame(hpos = hpos, bpos = bpos,
                       score = score)

  if (xIsHelix) {
    dframe$bpos <- factor(as.factor(dframe$bpos),
                          levels = c("-1", "1", "2", 
                                     "3", "5", "6"))
  } else {
    dframe$bpos <- factor(as.factor(dframe$bpos),
                          levels = c("1", "2", "3"))
  }

  if (yIsHelix) {
    dframe$hpos <- factor(as.factor(dframe$hpos),
                          levels = c("-1", "1", "2", 
                                     "3", "5", "6"))
  } else {
    dframe$hpos <- factor(as.factor(dframe$hpos),
                          levels = c("1", "2", "3"))
  }

  dframe

}

makeHeatPlot <- function(fname, dframe, xl, yl) {
  # Makes a heatplot
  print(fname)
  br <- seq(0,1, 0.05)
  g <- ggplot(dframe, aes(hpos, bpos)) + 
    geom_tile(aes(fill = score)) +
    scale_fill_gradient2("Normalized\nmutual\ninformation", breaks = br,
                         low = 'white', high = 'royalblue',
                         limits = c(0,0.4), guide = 'legend') +
    xlab(xl) +
    ylab(yl) +
    theme_bw() +
    #guides(colour = guide_legend(nrow = 3, title.position = "bottom"))+
    theme(legend.key.size = unit(0.4, "cm"))#,
          #legend.title = element_blank())
  ggsave(fname, plot = g, height = 2.7, width = 5.5)
}

runAllAnalysis(filtDirs, outDirs, filters)