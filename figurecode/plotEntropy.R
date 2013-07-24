library(ggplot2)

frequencyVsNumCombos <- function(data) {
	# Boxplots of frequencies grouped by number 
	# of combinations
	ggplot(data, aes(x = as.factor(V5), y = V3)) +
	geom_boxplot(fill = 'dodgerblue1') + 
	scale_x_continuous(limits = c(0, 0.0015))
	xlab("# Possible Codon Combinations") +
	ylab("Protein Frequency") +
	theme_bw() + 
	opts(axis.text.x=theme_text(angle=90, hjust=1)) 
}

targetVsEntropybp <- function(data){
	# Boxplot the entropy for each target
	ggplot(data, aes(x = as.factor(V1), y = V7)) +
	geom_boxplot(fill = 'dodgerblue1') + 
	xlab("Target Sequence") +
	ylab("Normalized Shannon Entropy") +
	theme_bw() +
	opts(axis.text.x=theme_text(angle=90, hjust=1))
}

countVsEntropybp <- function(data) {
	# Boxplot the entropy for each 
	# number of codon combinations
	ggplot(data, aes(x = as.factor(V5), y = V7)) +
	geom_boxplot(fill = 'dodgerblue1') + 
	xlab("# Possible Codon Combinations") +
	ylab("Normalized Shannon Entropy") +
	theme_bw() + 
	opts(axis.text.x=theme_text(angle=90, hjust=1))
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

codonComboScatter <- function(data) {
	# Number of possible codons combinations vs.
	# fraction of possible codon combinations observed
	ggplot(data, aes(x = V5, y = V4/V5)) +
	geom_point(position = position_jitter(w = 0.2, h =0.0),
	           size = 0.5) +
	scale_y_continuous(limits = c(0, 1)) + 
	stat_smooth(method = lm, color = 'blue') + 
	xlab('# Possible Codon Combinations') + 
	ylab('# Observed Codon Comb./# Possible Codon Comb.') + 
	theme_bw()
}

codonComboBox <- function(data, type = 'rel', ylims = NULL) {
	# Number of possible codons combinations vs.
	# fraction of possible codon combinations observed
	# If type = 'abs' then will plot actual observed
	# number on y axis instead of frequency

	if (type == 'rel'){
		g <- ggplot(data, aes(x = as.factor(V5), y = V4/V5)) +
		ylab('Fraction Observed')
	}
	else if (type == 'abs'){
		g <- ggplot(data, aes(x = as.factor(V5), y = V4)) +
		ylab('# Codon Combinations Observed') 
	}
	
	if (!is.null(ylims))
		g <- g + scale_y_continuous(limits = ylims)

	g <- g + 
	geom_boxplot(fill = 'dodgerblue1') +
	xlab('# Possible Codon Combinations') + 
	theme_bw() + 
	opts(axis.text.x=theme_text(angle=90, hjust=1))

	g
}

codonComboBar <- function(data) {
	# Plots a bar graph for how many proteins 
	# there are that have x number of codon combinations
	ggplot(data, aes(x = as.factor(V5))) +
	geom_bar() + 
	xlab('# Possible Codon Combinations') + 
	ylab('# of Proteins') +
	theme_bw() + 
	opts(axis.text.x=theme_text(angle=90, hjust=1))
}

printPlot <- function (data, figurePath, plotFunc, ...) {
	# Just calls a plotting function.
	pdf(figurePath, height=4.55,width=6.25)
	print(plotFunc(data, ...))
	dev.off()
}

main <- function() {
	fings <- c('F2')#c('F1', 'F2', 'F3')
	strins <- c('low', 'high')
	kept <- c('cut10')
	for (f in fings)
		for (s in strins)
			for (k in kept) {
				npos <- '_5pos_'
				pref <- paste0(f,'_',s,npos,k)
				
				fdir <- '../../figures/entropy_plots/cut10_filter0s/'
				data <- read.delim(paste0('../../data/b1hdata/',
				                   'newDatabase/5varpos/', f, '/', s, 
				                   '/protein_seqs_JSD/all_cut10.txt'),
									header = FALSE)
				
				# Uncomment if want to remove all proteins with 
				# 0 entropy.
				#data <- subset(data, V7 > 0)

				# Uncomment if want to remove 0 entropy proteins 
				# EXCEPT for proteins that only CAN code one way.
				data <- subset(data, V7 > 0 | 
				               (V4 == 1 & V5 == 1))

				# Print all the plots
				printPlot(data,paste0(fdir, pref,'_entropy_hist.pdf'),
				          entropyHistogram)
                if (k == 'cut10') {
                	printPlot(data, 
                	          paste0(fdir,pref,'_entropyVfreq_scatter.pdf'),
               	          freqVsEntropyscatter)
                	printPlot(data, 
                	          paste0(fdir,pref,'_entropyVfreq_scatter_zoom.pdf'),
                	          freqVsEntropyscatter, c(0, 0.0015))
                }
                printPlot(data, 
                          paste0(fdir,pref,'_entropyVnumComb_box.pdf'),
                          countVsEntropybp)
                printPlot(data, 
                          paste0(fdir,pref,'_entropyVtarget_box.pdf'),
                          targetVsEntropybp)

				#printPlot(data, 
                #          paste0(fdir,pref,'_codonCombo_box_abs.pdf'),
                #          codonComboBox, type = 'abs')				
                #printPlot(data, 
                #          paste0(fdir,pref,'_codonCombo_box_abs_zoom.pdf'),
                #          codonComboBox, type = 'abs', ylims = c(0,10))
                #printPlot(data, 
                #          paste0(fdir,pref,'_codonCombo_box_rel.pdf'),
                #          codonComboBox, 'rel')
                #printPlot(data,
               # 	      paste0(fdir,pref,'_codonCombo_bar.pdf'),
               #           codonComboBar)
			}
}

main()