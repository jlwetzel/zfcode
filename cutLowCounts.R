cutLowCounts <- function(allprefix, rdataprefix){
	
	# Read in the full data file in a frame write frame the 
	# frame to an R object file.
	allfile <- paste0(allprefix,'/all.txt')
	data <- read.delim(allfile, header = FALSE)
	print(paste0("Processing", allfile, '.'))
	save(data, file=paste0(rdataprefix, 'all.rda'))

	# Get the sorted list of targets and make a vector for storing 
	# the largest count(frequency) that we want to keep
	sortedTargs <- sort(unique(data$V1))
	minCounts <- vector(mode = "numeric", 
	                    length = length(sortedTargs))	

	# Get a list of ten smallest frequencies for each target,
	# subset the larger data frame based on these frequencies, 
	# then rbind the subsets together to make cut10data frame.
	cut10data <- data.frame()
	lastrows <- 0
	for (i in 1:length(sortedTargs)) {
		sortedCounts <- sort(unique(subset(data, 
		                     V1 == sortedTargs[i])$V3),
		                     decreasing = TRUE) 
    	minCounts[i] <- tail(sortedCounts, 10)[1]
    	cut10data <- rbind(cut10data, subset(data, 
    	                   V3 > minCounts[i] & V1 == sortedTargs[i]))
    	rows <- nrow(cut10data) - lastrows
    	lastrows <- nrow(cut10data)
    	print(paste0(sortedTargs[i], ":  ", rows, " (", 
    	      	     rows/nrow(subset(data, V1 == sortedTargs[i])), ")"))
	}

	# Print some possibly useful statistics
	print(paste0("Total num rows before removing low counts:", 
	      nrow(data)))
	print(paste0("Total num rows after removing low counts:", 
	      nrow(cut10data)))
	print(paste0("Ratio: ", nrow(cut10data)/nrow(data)))

	# Write the frame to file and write the frame and 
	# other info calulated to an R object file.
	write.table(cut10data, file = paste0(allprefix, '/all_cut10.txt'), 
	            sep = "\t", quote = FALSE, eol = '\n', 
	            row.names = FALSE, col.names = FALSE)
	save(cut10data, sortedTargs, minCounts, 
	     file = paste0(rdataprefix, 'cut10.rda'))
}


main <- function() {
	npos <- '_5pos_'#_6pos_'
	nposdir <- '5varpos'#6varpos'
	fings <- c('F2', 'F3')#c('F1', 'F2', 'F3')
	strins <- c('low', 'high')
	for (f in fings)
		for (s in strins){
			npos <- '_6pos_'
			rdataprefix <- paste0('../rsessions/',f,'_',s,npos)
			allprefix <- paste('../data/b1hdata/newDatabase',
			                   nposdir, f, s, 'protein_seqs_JSD',
			                   sep = '/')
			cutLowCounts(allprefix, rdataprefix)
		}
}

main()