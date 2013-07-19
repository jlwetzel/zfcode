library(ggplot2)

frequencyVsNumCombos <- function(data) {
	# Boxplots of frequencies grouped by number 
	# of combinations
	ggplot(data, aes(x = as.factor(V5), y = V3)) +
	geom_boxplot() + 
	scale_x_continuous(limits = c(0, 0.0015))
	xlab("# Possible Codon Combinations") +
	ylab("Protein Frequency") +
	#opts(axis.text.x=theme_text(angle=90, hjust=1)) + 
	theme_bw()
}

targetVsEntropybp <- function(data){
	# Boxplot the entropy for each target
	ggplot(data, aes(x = as.factor(V1), y = V7)) +
	geom_boxplot() + 
	xlab("Target Sequence") +
	ylab("Normalized Shannon Entropy") +
	#opts(axis.text.x=theme_text(angle=90, hjust=1)) + 
	theme_bw()
}

countVsEntropybp <- function(data) {
	# Boxplot the entropy for each 
	# number of codon combinations
	ggplot(data, aes(x = as.factor(V5), y = V7)) +
	geom_boxplot() + 
	xlab("# Possible Codon Combinations") +
	ylab("Normalized Shannon Entropy") +
	#opts(axis.text.x=theme_text(angle=90, hjust=1)) + 
	theme_bw()
}

entropyHistogram <- function(data) {
	# Histogram of the entropy
	ggplot(data, aes(x = V7)) + 
	geom_histogram(binwidth = 0.01) + 
	xlab("Normalized Shannon Entropy") +
	ylab("# of Proteins") + 
	theme_bw()
}

freqVsEntropyscatter <- function(data, xlim = NULL) {
	# Frequency of protein vs entropy scatterplot
	if (is.null(xlim))
		xlim = c(0, max(data$V3))
	ggplot(data, aes(x = V3, y = V7)) +
	geom_point(size = 0.5) + 
	scale_x_continuous(limits = xlim) +
	scale_y_continuous(limits = c(0,1)) + 
	stat_smooth(method = lm, col = 'blue') + 
	xlab("Protein Frequency") +
	ylab("Normalized Shannon Entropy") +
	theme_bw()
}

main <- function() {
	fings <- c('F1', 'F2', 'F3')
	strins <- c('low', 'high')
	kept <- c('cut10')
	for (f in fings)
		for (s in strins)
			for (k in kept) {
				npos <- '_6pos_'
				pref <- paste0(f,'_',s,npos,k)
				#dfile <- paste0('../../rsessions/', pref, '.rda')
				fdir <- '../../figures/entropy_plots/'
				#load(dfile)
				#if (k == 'cut10')
			#		data <- cut10data

				data <- read.delim(paste0('../../data/b1hdata/',
				                   'newDatabase/6varpos/', f, '/', s, 
				                   '/protein_seqs_JSD/all_cut10.txt'),
									header = FALSE)
				data <- subset(data, V7 > 0.0 | 
				               (V7 == 0.00 & V4 == 1 & V5 == 1))

				#print(dfile)
				#print(nrow(data))
				#print(data)

				# The histogram
				pdf(paste0(fdir, pref,'_entropy_hist.pdf'), height=4.55,
				    width=6.25)
				print(entropyHistogram(data))
				dev.off()

				if (k == 'cut10') {
					# The scatterplots
					#pdf(paste0(fdir, pref,'_entropy_scatter.pdf'), height=4.55,
				    #	width=6.25)
					#print(freqVsEntropyscatter(data))
					#dev.off()
					pdf(paste0(fdir, pref,'_entropyVfreq_scatter_zoom.pdf'), height=4.55,
				    	width=6.25)
					print(freqVsEntropyscatter(data, c(0, 0.0015)))
					dev.off()
				}

				# The boxplots
				pdf(paste0(fdir, pref,'_entropyVnumComb_box.pdf'), height=4.55,
				    width=6.25)
				print(countVsEntropybp(data))
				dev.off()

				#pdf(paste0(fdir, pref,'_freqVnumComb_box.pdf'), height=4.55,
				#    width=6.25)
				#print(frequencyVsNumCombos(data))
				#dev.off()

				#pdf(paste0(fdir, pref,'_freqVnumComb_box_zoom.pdf'), height=4.55,
				#    width=6.25)
				#print(frequencyVsNumCombos(data))
				#dev.off()

				pdf(paste0(fdir, pref,'_entropyVtarget_box.pdf'), height=4.55,
				    width=6.25)
				print(targetVsEntropybp(data))
				dev.off()
			}
}

main()