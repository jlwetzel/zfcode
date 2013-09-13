library(ggplot2)

data <- read.csv("../../figures/aminoFreqs/F2_low_combinedFilters.txt")
outFileStem <- "../../figures/aminoFreqs/F2_low_"
helixPosNames <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')

for (apos in helixPosNames) {
	filter <- "filter"
	g <- ggplot(data, aes_string(x = apos, fill = filter)) +
		geom_bar() +
		theme_bw() +
		xlab(paste("Residue (pos ", apos, ")", sep = '')) + 
		ylab("# Times Observed")
	ggsave(file = paste(outFileStem, apos, '.pdf', sep = ''), plot = g)
}
