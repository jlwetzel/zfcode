library(ggplot2)

# Get the set of triplets
  bases <- c('A', 'C', 'G', 'T')
  triplets <- vector()
  for (b1 in bases)
    for (b2 in bases)
      for (b3 in bases)
        triplets <- c(triplets, paste0(b1, b2, b3))

prefix <- '../../figures/positionVariation/'
fing <- 'F2'
strin <- 'union'
filt <- 'filt_10e-4_025_0_c'
pos <- c('a-1', 'a2', 'a3', 'a6')

tabLab <- paste0(prefix, paste(fing, strin, filt, sep = '_'))
tabLabs <- vector()
for (p in pos) {
	tLab <- paste(paste(tabLab, p, sep = '_'))
	tabLabs <- c(tabLabs, tLab)
}

for (tl in tabLabs){
	tFile <- paste0(tl, '.txt')
	data <- read.table(tFile, header = TRUE)
	print(head(data))
	data$t2 <- factor(data$t2, levels = rev(triplets))
	g <- ggplot(data, aes(x = t1, y = t2)) +
		geom_tile(aes(fill = score)) +
		scale_fill_gradient2(limits = c(-max(data$score), max(data$score)), 
		                     midpoint=0, 
		                     low = "red", mid="white", 
		                     high = "green") + 

		geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = 64.5)) + 
    geom_segment(aes(x = 16.5, xend = 16.5, y = 0.5, yend = 65 - 16.5)) + 
    geom_segment(aes(x = 32.5, xend = 32.5, y = 0.5, yend = 65 - 32.5)) + 
    geom_segment(aes(x = 48.5, xend = 48.5, y = 0.5, yend = 65 - 48.5)) +
    geom_segment(aes(y = 0.5, yend = 0.5, x = 0.5, xend = 64.5)) +  
    geom_segment(aes(y = 16.5, yend = 16.5, x = 0.5, xend = 65 - 16.5)) + 
    geom_segment(aes(y = 32.5, yend = 32.5, x = 0.5, xend = 65 - 32.5)) + 
    geom_segment(aes(y = 48.5, yend = 48.5, x = 0.5, xend = 65 - 48.5)) + 

		theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(),
          legend.position = c(0.65, 0.65),
          legend.background = element_rect(fill="white", color = "black",
                                           size=.5)) +

    #draws x and y axis line
    theme(axis.line = element_line(color = 'black', size = 0.5)) + 

    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle=90, vjust=0.5, size = 8))

	plotName <- paste0(tl, '.pdf') 
	ggsave(plotName, plot = g, width = 7.75, height = 7.5)
} 
