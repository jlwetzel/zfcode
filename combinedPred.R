# Trains one predictor per target 3-mer and combines their 
# outputs to give a general prediction of specificity
# given a protein.

library(caret)
library(doMC)
registerDoMC(cores = 8)

bases <- c('A', 'C', 'G', 'T')
apos <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')

getTargDframe <- function(path, targ, writeToFile = FALSE) {
	# Creates a dataframe from the targetfile where each 
	# where first column is frequency followed by 
	# a0, a1, a2, a3, a5, a6

	targData <- read.table(paste0(path, targ, '.txt'), 
	                       sep = '\t',
	                       stringsAsFactors = FALSE)
	
	aList <- list()
	freqVec <- vector()
	for (ap in apos)
		aList[[ap]] <- vector()

	for (i in 1:nrow(targData)) {
		aseq <- targData[i,1]
		freqVec <- c(freqVec, targData[i,2])
		for (j in 1:length(apos)) {
			ap <- apos[j]
			aList[[ap]] <- c(aList[[ap]], substr(aseq, j, j)) 
		}
	}

	df <- data.frame(freq = freqVec)
	for (ap in names(aList))
		df <- cbind(df, aList[[ap]])
	names(df) <- c('freq', names(aList))

	if (writeToFile)
		write.csv(df, file = paste0(path, targ, '.csv'))

	df
}

addFreqs <- function(dframe, combineOn) {
	# Add the frequencies for each dataframe in 
	# targ data if they share the same values for 
	# all of the variable names in the combineOn vector

	# Get a unique set of keys
	keys <- vector()
	for (i in 1:nrow(dframe)) {
		freq <- dframe[i, 1]
		seq <- ""
		for (col in combineOn)
			seq <- paste0(seq, as.character(dframe[i, col]))
		keys <- c(keys, seq)
	keys <- unique(keys)
	freqs <- vector(mode = "numeric", length = length(keys))
	names(freqs) <- keys
	for (n in names(freqs))
		freqs[n] <- 0
	}

	# Assign values to the keys
	for (i in 1:nrow(dframe)) {
		freq <- dframe[i, 1]
		seq <- ""
		for (col in combineOn)
			seq <- paste0(seq, as.character(dframe[i, col]))
		freqs[seq] <- freqs[seq] + freq
	}

	# Convert key/value pairs back to unique observations
	aList = list()
	freqVec = vector()
	for (ap in combineOn)
		aList[[ap]] <- vector()
	for (n in names(freqs)) {
		freqVec <- c(freqVec, freqs[n])
		for (j in 1:length(combineOn)) {
			ap <- combineOn[j]
			aList[[ap]] <- c(aList[[ap]], substr(n, j, j)) 
		}
	}

	# Return the new dataframe
	df <- data.frame(freq = freqVec)
	for (ap in names(aList))
		df <- cbind(df, aList[[ap]])
	names(df) <- c('freq', names(aList))
	df
}

getDframes <- function(path, combineOn) {
	# Returns a list of dataFrames, one per DNA target

	targData <- list()

	targs <- vector()
	for (b1 in bases)
		for (b2 in bases)
			for (b3 in bases)
				targs <- c(targs, paste0(b1,b2,b3))



	for (targ in targs) {
		print(targ)
		targData[[targ]] <- getTargDframe(path, targ, 
		                                  writeToFile = FALSE)
		targData[[targ]] <- addFreqs(targData[[targ]],
		                                 combineOn)
	}
	targData
}

fitPred <- function(dframe, formStr, fitControl, predType) {
	train(form = as.formula(formStr), data = dframe,
	      method = 'lm',
	      trControl = fitControl)
}

getAllPredictors <- function(targFrames, formStr, predType = 'lm') {
	# Returns one ols predictor for each of the targs

	olsPreds <- list()
	fitControl <- trainControl(method = "repeatedcv", number = 10,
	                           repeats = 10)
	for (targ in names(targFrames))
		olsPreds[[targ]] <- fitPred(targFrames[[targ]],
		                            formStr, fitControl,
		                            predType)
	olsPreds
}

