library(ggplot2)

nucs = c('A', 'C', 'G', 'T')
aminos = c('A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

getBaseProbs <- function(allData) {
	# Get sums for independent probs for bases
	contacts <- names(allData[2:length(names(allData))])
	baseProbs <- list()
	for (c in contacts) {
		baseProbs[[c]] <- vector(length = length(nucs))
	  	names(baseProbs[[c]]) <- nucs
	  	for (b in nucs) {
	    	bfreqs <- allData$freq[substr(as.character(allData[,c]),
	        	                         1, 1) == b]
	    	baseProbs[[c]][b] <- sum(bfreqs)/sum(allData$freq)
	  	}
	}
	baseProbs
}

getAminoProbs <- function(allData) {
	# Get sums for independent probs for aminos
	contacts <- names(allData[2:length(names(allData))])
	aminoProbs <- list()
	for (c in contacts) {
		aminoProbs[[c]] <- vector(length = length(aminos))
	  	names(aminoProbs[[c]]) <- aminos
	  	for (a in aminos) {
	    	afreqs <- allData$freq[substr(as.character(allData[,c]),
	     	                             2, 2) == a]
	    	if (length(afreqs) == 0)
	      		aminoProbs[[c]][a] <- NA
	    	else
	      		aminoProbs[[c]][a] <- sum(afreqs)/sum(allData$freq)
	  	}
	}
	aminoProbs
}

getJointProbs <- function(allData) {
# Log likelihood ratio analysis
	contacts <- names(allData)[2:length(names(allData))]
	jProbs <- list()
	for (c in contacts){
	  	jProbs[[c]] <- matrix(nrow = 4, ncol = 20)
	  	rownames(jProbs[[c]]) <- nucs
	  	colnames(jProbs[[c]]) <- aminos
	  	#print("Here?")
	  	for (i in 1:length(nucs)) {
	    	for (j in 1:length(aminos)) {
	      		b <- nucs[i]
	      		a <- aminos[j]
	      		pair <- paste0(b,a)
	      		pfreqs <- allData$freq[as.character(allData[,c])
	            		                 == pair]
	      		if (length(pfreqs) == 0)
	        		jProbs[[c]][i,j] <- NA
	      		else
	      			jProbs[[c]][i, j] <- sum(pfreqs)/sum(allData$freq)
	    	}
	  	}
	}
	jProbs
}

getLLmats <- function(bProbs, aProbs, jProbs) {
	contacts <- names(jProbs)
	llmats <- list()
	for (c in contacts) {
		llmats[[c]] <- matrix(nrow = nrow(jProbs[[c]]),
		                      ncol = ncol(jProbs[[c]]))
		rownames(llmats[[c]]) <- rownames(jProbs[[c]])
		colnames(llmats[[c]]) <- colnames(jProbs[[c]])
		for (i in 1:nrow(jProbs[[c]])) {
	    	for (j in 1:ncol(jProbs[[c]])) {
	      		b <- names(jProbs[[c]])[i]
	      		a <- names(aProbs[[c]])[j]
	      		pair <- paste0(b,a)
	      		if (is.na(bProbs[[c]][i]) || is.na(aProbs[[c]][j]) || is.na(jProbs[[c]][i,j]))
	        		llmats[[c]][i,j] <- NA
	      		else {
	      			jointProb <- jProbs[[c]][i,j]
	      			indepProb <- bProbs[[c]][i] * aProbs[[c]][j]
	      			llmats[[c]][i, j] <- log2(jointProb/indepProb)
	      		}
	    	}
	  	}
	}
	llmats
}

matsToFrame <- function(mats) {
  	contacts <- names(mats)
  	con <- vector()
  	base <- vector()
  	score <- vector()
  	amino <- vector()
  	for (c in contacts) {
    	mat <- mats[[c]]
    	for (i in 1:nrow(mat)) {
      		for (j in 1:ncol(mat)) {
       			con <- c(con, c)
        		base <- c(base, rownames(mat)[i])
        		amino <- c(amino, colnames(mat)[j])
        		score <- c(score, mat[i,j])
      		}
    	} 
  	}
  	data.frame(contact = con, base = base, amino = amino,
    	       score = score)
}

makeScorePlot <- function(fname, dframe){
	g <- ggplot(dframe, aes(amino, base)) +
    geom_tile(aes(fill = score)) +
    scale_fill_gradient2(low = "red", mid = "white", 
                        high = "blue", na.value = "gray50") +
    facet_wrap(~contact)
 	ggsave(fname, plot = g)
}

writeFrame <- function(fname, dframe){
	write.table(dframe, fname, sep = '\t', row.names = FALSE,
	            quote = FALSE)
}

runllAnalysis <- function(inDir, outDir, fname, label) {
	allData <- read.csv(paste(inDir, fname, sep = '/'))
	baseProbs <- getBaseProbs(allData)
	aminoProbs <- getAminoProbs(allData)
	jointProbs <- getJointProbs(allData)
	llmats <- getLLmats(baseProbs, aminoProbs, jointProbs)
	llmatsFrame <- matsToFrame(llmats)
	makeScorePlot(paste(outDir, paste0("llRatios_",label,".pdf"), 
	              sep = '/'),llmatsFrame)
	writeFrame(paste(outDir, paste0("llRatios_",label,".txt"),
	           sep = '/'),llmatsFrame)
}

runAll_llAnalysis <- function() {
	fname <- 'all_factor_c_add.csv'
	label <- 'c'
	inDirPref <- '../../data/b1hData/newDatabase/6varpos'
	outDirPref <- '../../figures/llAnalysis'
	fingers <- c('F2','F3')
	strins <- c('low', 'high')
	filts <- c('cut10bc_025', 'cut10bc_0_5',
	           'cut3bc_025', 'cut3bc_0_5')
	inDirs <- vector()
	outDirs <- vector()
	for (f in fingers)
		for (s in strins)
			for (filt in filts){
				protfilt <- paste('protein_seq', filt, sep = '_')
				inDirs <- c(inDirs, 
				            paste(inDirPref, f, s, protfilt, sep = '/'))
				outDirs <- c(outDirs,
				             paste(outDirPref, f, s, filt, sep = '/'))
			}

	for (i in 1:length(inDirs)) {
		runllAnalysis(inDirs[i], outDirs[i], fname, label)
	}
}

runAll_llAnalysis()




