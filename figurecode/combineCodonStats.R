# Combines the codon statistics for the various finger/stringency
# combinations and outputs ggplot heatmaps against expected
# frequencies for those codons.

library(lattice)

# Normalized per amino acid residue
ecolibiasMat1 <- matrix(c(rep(0.04761905, 6), rep(0.04761905, 6),
                       	 rep(0.03091061, 6), rep(0.01670844, 6),
                       	 rep(0.02255639, 6), rep(0.00297619, 6),
                       	 rep(0.04761905, 6), rep(0.04761905, 6),
                       	 rep(0.04761905, 6), rep(0.04761905, 6),
                       	 rep(0.00732601, 6), rep(0.04029304, 6),
                       	 rep(0.03670635, 6), rep(0.00793651, 6),
                       	 rep(0.00626566, 6), rep(0.03446115, 6),
                       	 rep(0.04761905, 6), rep(0.04761905, 6),
                       	 rep(0.02017756, 6), rep(0.02744149, 6),
                       	 rep(0.0359389,  6), rep(0.01168014, 6),
                       	 rep(0.01763668, 6), rep(0.02998236, 6),
                       	 rep(0.04761905, 6), rep(0.04761905, 6),
                       	 rep(0.01420217, 6), rep(0.01086048, 6),
                       	 rep(0.04761905, 6), rep(0.04761905, 6),
                       	 rep(0.04761905, 6), rep(0.00689223, 6)),
						 nrow = 32, ncol=6, byrow = TRUE)

# Not normalized per amino acid residue
ecolibiasMat2 <- matrix(c(rep(2.6, 6), rep(1.2, 6),
                       	 rep(2.4, 6), rep(1.3, 6),
                       	 rep(1.5, 6), rep(0.2, 6),
                       	 rep(2.7, 6), rep(2.6, 6),
                       	 rep(1.1, 6), rep(2.9, 6),
                       	 rep(0.4, 6), rep(2.4, 6),
                       	 rep(2.2, 6), rep(0.5, 6),
                       	 rep(0.9, 6), rep(5.2, 6),
                       	 rep(2.3, 6), rep(1.9, 6),
                       	 rep(2.3, 6), rep(3.2, 6),
                       	 rep(3.0,  6), rep(0.9, 6),
                       	 rep(1.4, 6), rep(2.4, 6),
                       	 rep(1.4, 6), rep(0.000001, 6),
                       	 rep(1.0, 6), rep(0.8, 6),
                       	 rep(0.6, 6), rep(1.4, 6),
                       	 rep(1.8, 6), rep(1.1, 6)),
						 nrow = 32, ncol=6, byrow = TRUE)
for (i in 1:ncol(ecolibiasMat2))
     ecolibiasMat2[,i] = ecolibiasMat2[,i]/sum(ecolibiasMat2[,i])

unifbgMat <- matrix(c(rep(1/32, 32)), nrow = 32, ncol = 6, byrow=TRUE)

getmatchStrs <- function(matchStr, strList){
	# Get the list of strings that match the regex
	matches <- regexpr(matchStr, strList)
	matches <- ifelse(matches == 1, TRUE, FALSE)
	strList[matches]
}

combineMats <- function(path, flist){
	# Returns a matrix that averages all the matrices
	# given by the files in flist
	
	# figure out the size of the matrix
	data <- read.delim(paste0(path, flist[1]))
	mat <- as.matrix(data[,2:ncol(data)])
	totalMat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
	rownames(totalMat) <- data[,1]
	colnames(totalMat) <- colnames(mat)

	for (f in flist) {
  		data <- read.delim(paste0(path, f))
  		mat <- as.matrix(data[,2:ncol(data)])
  		totalMat = totalMat + mat
	}
	totalMat <- totalMat/length(flist)
	totalMat
}

makeLevelPlot <- function(fname, mat) {
	# Makes a levelplot for the matrix values
	bluered = colorRampPalette(c("red", "white", "blue"),
	                           space = "Lab") 
	pdf(fname, height = 4.55, width = 2.84)
	print(levelplot(t(mat), at=seq(-5,5,length = 100), 
	                col.regions = bluered,
	                xlab = "Codon Position",
	                ylab = "Codon"))
	dev.off
}

main <- function() {
	fings <- c('F1', 'F2', 'F3')
	strins <- c('low', 'high')

	for (fing in fings)
		for (strin in strins){
			prefix <- '../../data/b1hdata/newDatabase/6varpos'
			path <- paste(prefix, fing, strin, "combined_nuc_seq",
			              "statistics/", sep = '/')
			print(path)

			# Get the list of matrix files
			flist <- list.files(path)
			matchStr <- "[ACGT]{3}(.*)_codonStats.txt"
			flist <- getmatchStrs(matchStr, flist)
			
			# Get the combined matrix and write to file
			totalMat <- combineMats(path, flist)
			write.table(totalMat, file = paste0(path, 
			            'all_codonStats.txt'))

			# Make a log-odds plot for ecoli codon bias an 
			# uniform codon background frequencies
			plotpath <- paste0('../../figures/codonBias/',
			                   fing, '_', strin, '_ecoli1.pdf')
			makeLevelPlot(plotpath, log2(totalMat/ecolibiasMat1))
			plotpath <- paste0('../../figures/codonBias/',
			                   fing, '_', strin, '_ecoli2.pdf')
			makeLevelPlot(plotpath, log2(totalMat/ecolibiasMat2))
			plotpath <- paste0('../../figures/codonBias/',
			                   fing, '_', strin, '_unif.pdf')
			makeLevelPlot(plotpath, log2(totalMat/unifbgMat))
		}
}

main()