# A script for plotting some figures for related to analyzing
# the usefullness of the JSD filter for the binding data.

library(ggplot2)

countVsJSDbp <- function(data) {
	ggplot(data, aes(x = as.factor(V5), y = V6)) +
	geom_boxplot() + 
	xlab("# Possible Condon Combinations") +
	ylab("Jensen-Shannon Divergence") +
	theme_bw()
}

jsdHistogram <- function(data) {
	ggplot(data, aes(x = V6)) + 
	geom_histogram(binwidth = 0.01) + 
	xlab("Jensen-Shannon Divergence") +
	ylab("# of Proteins") + 
	theme_bw()
}

freqVsJSDscatter <- function(data, xlim = NULL) {
	if (is.null(xlim))
		xlim = c(0, max(data$V3))
	ggplot(data, aes(x = V3, y = V6)) +
	geom_point(size = 0.5) + 
	scale_x_continuous(limits = xlim) +
	scale_y_continuous(limits = c(0,1)) + 
	xlab("Protein Frequency") +
	ylab("Jensen-Shannon Divergence") +
	theme_bw()
}

main <- function() {
	fings <- c('F1', 'F2', 'F3')
	strins <- c('low', 'high')
	kept <- c('all', 'cut10')
	for (f in fings)
		for (s in strins)
			for (k in kept) {
				npos <- '_6pos_'
				pref <- paste0(f,'_',s,npos,k)
				dfile <- paste0('../../rsessions/', pref, '.rda')
				fdir <- '../../figures/jsd_plots/'
				load(dfile)
				if (k == 'cut10')
					data <- cut10data

				# The histogram
				pdf(paste0(fdir, pref,'_JSDhist.pdf'), height=4.55,
				    width=6.25)
				print(jsdHistogram(data))
				dev.off()

				if (k == 'cut10') {
					# The scatterplots
					pdf(paste0(fdir, pref,'_JSDscatter.pdf'), height=4.55,
				    	width=6.25)
					print(freqVsJSDscatter(data))
					dev.off()
					pdf(paste0(fdir, pref,'_JSDscatter_zoom.pdf'), height=4.55,
				    	width=6.25)
					print(freqVsJSDscatter(data, c(0, 0.0015)))
					dev.off()
				}

				# The boxplot
				pdf(paste0(fdir, pref,'_JSDbox.pdf'), height=4.55,
				    width=6.25)
				print(countVsJSDbp(data))
				dev.off()
			}
}

main()