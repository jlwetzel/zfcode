library(infotheo)
library(ggplot2)

timesCanonicalHelixObserved <- function(pathPref,
                                        fing, strin, filt,
                                        skipParse = FALSE){
	# Find out how often we see each distinct canonical
	# helic per triplet.

	outDir <- '../../figures/fig1supp/timesCanonObserved'
	helixPos <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')

	
	# Get the set of triplets
	bases <- c('A', 'C', 'G', 'T')
	triplets <- vector()
		for (b1 in bases)
  		for (b2 in bases)
  			for (b3 in bases)
      		triplets <- c(triplets, paste0(b1, b2, b3))

  if (!skipParse) {
	  fpath <- paste(pathPref, fing, strin, filt, 
	  	             'all.csv', sep = '/')
	  df <- read.csv(file = fpath, header = TRUE)

	  # Find mean and sd of # of times each can. helix observed
	  # per triplet
	  trip <- vector()
	  avg <- vector()
	  med <- vector()
	  stDev <- vector()
	  for (t in triplets) {
	  	dfsub <- subset(df, b1 == substr(t, 1, 1) &
	  	                b2 == substr(t, 2, 2) &
	  	                b3 == substr(t, 3, 3))
			canHelices <- vector()
			for (i in 1:nrow(dfsub)) {
				canHelix <- paste0(as.character(dfsub[i,]$a0),
				                   as.character(dfsub[i,]$a2),
				                   as.character(dfsub[i,]$a3),
				                   as.character(dfsub[i,]$a6))
				canHelices <- c(canHelices, canHelix)
			}
			trip <- c(trip, t)
			avg <- c(avg, mean(as.vector(table(canHelices))))
			med <- c(med, median(as.vector(table(canHelices))))
			stDev <- c(stDev, sd(as.vector(table(canHelices))))
			print(t)
			cFrame <- data.frame(count = as.vector(table(canHelices)))
			#print(cFrame$count)
			g <- ggplot(cFrame, aes(count)) +
				geom_histogram(fill = "royalblue") +
				xlab("Number of times canonical helix observed") +
				ylab("Count") +
				theme_bw()
			plotName <- paste0(paste(paste(outDir, paste0(t, "_timesCanonObserved"), sep = '/'), 
                  fing, strin, filt, "hist", sep = '_'), '.eps')
			ggsave(plotName, plot = g)
			print(as.vector(table(canHelices)))
			#print(sd(table(canHelices)))
  	}
  	df <- data.frame(trip = trip, avg = avg, stDev = stDev, med = med)
  	tableName <- paste0(paste(paste(outDir, "timesCanonObserved", sep = '/'), 
                      fing, strin, filt, sep = '_'), '.txt')
  	write.table(df, file = tableName)
	} else {
		tableName <- paste0(paste(paste(outDir, "timesCanonObserved", sep = '/'), 
                      fing, strin, filt, sep = '_'), '.txt')
		df <- read.table(tableName)
	}

  # Barplot
  g <- ggplot(df, aes(x = trip, y = avg)) +
  geom_bar(stat = "identity", fill = "royalblue") + 
  geom_errorbar(aes(ymin=max(0, avg-stDev), ymax=avg+stDev), colour="black", width=.1) +
  coord_flip() + 
  xlab("DNA triplet") + 
  ylab("Number of times canonical helix observed") +
  theme_bw() + 
  theme(axis.text.y = element_text(size=6))
  plotName <- paste0(paste(paste(outDir, "timesCanonObserved", sep = '/'), 
                      fing, strin, filt, "barPlot", sep = '_'), '.eps')
  ggsave(plotName, plot = g)

  # Boxplot of means
  df$dummyVar <- as.factor(rep(1,length(df$trip)))
  g <- ggplot(df, aes(x = dummyVar, y = avg)) + 
  geom_boxplot(fill = "royalblue") + 
  ylab("Mean number of times canonical helix observed per DNA triplet") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank())#, 
        #axis.title.y = element_text(size = 10))
  plotName <- paste0(paste(paste(outDir, "timesCanonObserved", sep = '/'), 
                      fing, strin, filt, "boxPlot_mean", sep = '_'), '.eps')
  ggsave(plotName, plot = g)

  # Boxplot of medians
  df$dummyVar <- as.factor(rep(1,length(df$trip)))
  g <- ggplot(df, aes(x = dummyVar, y = med)) + 
  geom_boxplot(fill = "royalblue") + 
  ylab("Median number of times canonical helix observed per DNA triplet") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank())#, 
        #axis.title.y = element_text(size = 10))
  plotName <- paste0(paste(paste(outDir, "timesCanonObserved", sep = '/'), 
                      fing, strin, filt, "boxPlot_median", sep = '_'), '.eps')
  ggsave(plotName, plot = g)

}

entropyPerPosPerTriplet <- function(pathPref, 
                                      fing, strin, filt){
	# Find the per-positon entropy for helices 
	# from each triplet

	outDir <- '../../figures/fig1supp/entropyPerPos'
	helixPos <- c('a0', 'a1', 'a2', 'a3', 'a5', 'a6')

	# Get the set of triplets
	bases <- c('A', 'C', 'G', 'T')
	triplets <- vector()
		for (b1 in bases)
  		for (b2 in bases)
  			for (b3 in bases)
      		triplets <- c(triplets, paste0(b1, b2, b3))

  fpath <- paste(pathPref, fing, strin, filt, 
  	             'all.csv', sep = '/')
  df <- read.csv(file = fpath, header = TRUE)
  print(head(df))

  trip <- vector()
  hpos <- vector()
  ent <- vector()
  for (t in triplets) {
  	dfsub <- subset(df, b1 == substr(t, 1, 1) &
  	                b2 == substr(t, 2, 2) &
  	                b3 == substr(t, 3, 3))
		for (h in helixPos)	{
			trip <- c(trip, t)
			hpos <- c(hpos, h)
			ent <- c(ent, natstobits(entropy(dfsub[[h]])))
			#print(paste(trip, h, ent))
		}
  }
  entFrame <- data.frame(trip = trip, hpos = hpos, ent = ent)
  entFrame$hpos <- factor(as.factor(hpos), 
                          levels = helixPos,
                          labels = c('-1', '1', '2', '3', '5', '6'))
  
  tableName <- paste0(paste(paste(outDir, "entropyPerPos", sep = '/'), 
                      fing, strin, filt, sep = '_'), '.txt')
  write.table(entFrame, file = tableName)

  # Faceted barplot
  g <- ggplot(entFrame, aes(x = hpos, y = ent)) +
  geom_bar(stat = "identity", fill = "royalblue") + 
  facet_wrap(~trip) +
  xlab("Helix position") + 
  ylab("entropy") +
  theme_bw()
  plotName <- paste0(paste(paste(outDir, "entropyPerPos", sep = '/'), 
                      fing, strin, filt, "barPlot", sep = '_'), '.eps')
  ggsave(plotName, plot = g)

  # Boxplot
  g <- ggplot(entFrame, aes(x = hpos, y = ent)) +
  geom_boxplot(fill = "royalblue") + 
  xlab("Helix position") + 
  ylab("entropy") +
  theme_bw()
  plotName <- paste0(paste(paste(outDir, "entropyPerPos", sep = '/'), 
                      fing, strin, filt, "boxPlot", sep = '_'), '.eps')
  ggsave(plotName, plot = g)

}

numTargsVsWeightSpearmans <- function(pathPref, fing, strin, 
                                      filt, outDir, countFrame){
	# Compute spearman corelation between canoncial protein 
	# frequencies and number of targets.

	tableName <- paste0(paste(paste(outDir, "ntargsVfreqSpearman", sep = '/'), 
                      fing, strin, filt, sep = '_'), '.txt')

	# Get the set of triplets
	bases <- c('A', 'C', 'G', 'T')
	triplets <- vector()
		for (b1 in bases)
  		for (b2 in bases)
  			for (b3 in bases)
      		triplets <- c(triplets, paste0(b1, b2, b3))

  
  cts <- countFrame$ntargs
  names(cts) <- as.character(countFrame$prot)

  trip <- vector()
  rho <- vector()
  pval <- vector()
  for (t in triplets) {
  	fpath <- paste(pathPref, fing, strin, filt, 
	  	             paste0('4pos_',t,'.txt'), sep = '/')
  	df <- read.table(file = fpath, header = FALSE, sep = '\t')
  	names(df) <- c("prot", "freq")
  	count <- vector()
		for (prot in df$prot)
  		count <- c(count, cts[[prot]])

  	test <- cor.test(df$freq, count,
  	                 method = "spearman")
  	trip <- c(trip, t)
  	rho <- c(rho, test$estimate)
  	pval <- c(pval, test$p.value)
  	if (t == "AGC") {
  		print(cbind(df, count))
  	}

	}
	spearFrame <- data.frame(trip = trip, rho = rho, pval = pval)
	spearFrame$trip <- reorder(spearFrame$trip, spearFrame$rho)
	spearFrame$BHpval <- p.adjust(spearFrame$pval, method = "BH")
	write.table(spearFrame, file = tableName)
	print(spearFrame)

	
	g <- ggplot(spearFrame, aes(x = trip, y = rho)) +
	geom_bar(stat = 'identity', fill = 'royalblue') + 
	ylab("Spearman correlation (# of triplets bound vs. helix frequency)") + 
	xlab("DNA triplet") +
	coord_flip() +
	theme_bw() +
	theme(axis.text.y = element_text(size=6),
	      axis.title.x = element_text(size=10))

	plotName <- paste0(paste(paste(outDir, "ntargsVfreqSpearman", sep = '/'), 
                      fing, strin, filt, sep = '_'), '.pdf')
	ggsave(plotName, plot = g, width = 4.5)
	print(table(cut(spearFrame$BHpval, breaks = c(0, 0.001, 0.01, 0.05, 1))))
}

numTargsVsWeight <- function(pathPref, fing, strin,
                             filt, skipParse = FALSE) {
	# Examine correlation between number of targets bound 
	# and weight of the protein across the dataset
	
	outDir <- '../../figures/fig1supp/numTargsVsWeight'
	tableName <- paste0(paste(paste(outDir, "numTargsVsWeight", sep = '/'), 
                      fing, strin, filt, sep = '_'), '.txt')

	if (!skipParse){
		fpath <- paste(pathPref, fing, strin, filt, 
	  	             'bindVectsUnique4pos.txt', sep = '/')
	  df <- read.table(file = fpath, header = TRUE, sep = '\t')
	  #print(head(df))
	  prot <- vector()
	  wts <- vector()
	  numTargs <- vector()
	  for (i in 1:nrow(df)) {
	  	p <- as.character(df[i,]$prot)
	  	bindVect <- as.logical(df[i,4:ncol(df)])
	  	ntargs <- length(bindVect[bindVect == TRUE])
	  	wt <- (df[i,]$weight)/ntargs
	  	
	  	prot <- c(prot, p)
	  	wts <- c(wts, wt)
	  	numTargs <- c(numTargs, ntargs)
	  }
	  countFrame <- data.frame(prot = prot, weight = wts, 
	                           ntargs = numTargs)
	  write.table(countFrame, file = tableName)
	  print(head(countFrame))
	  } else {
	  	countFrame <- read.table(tableName)
	  }

	  print(head(countFrame))
	  g <- ggplot(countFrame, aes(x = ntargs, y = weight)) + 
	  geom_jitter(color = "royalblue", size = 0.5) +
	  stat_smooth(method = 'lm', color = 'black') + 
	  xlab("Number of distinct DNA triplets bound by protein") +
	  ylab("Weight of the protein") +
	  theme_bw()
	  plotName <- paste0(paste(paste(outDir, "numTargsVsWeight", sep = '/'), 
	                     fing, strin, filt, sep = '_'), '.eps')
	  ggsave(plotName, plot = g)

	  # Compute per-triplet spearman correlations
	  numTargsVsWeightSpearmans(pathPref, fing, strin, filt, 
	                            outDir, countFrame)
}

main <- function(){
	pathPref <- '../../data/b1hData/antonProcessed/'
	fing <- "F2F3"
	strin <- "unionUnions"
	filt <- 'filt_10e-4_025_0_c'
	#variationPerPosPerTriplet(pathPref, fing, strin, filt)
	#timesCanonicalHelixObserved(pathPref, fing, 
	#                            strin, filt, skipParse = FALSE)
	numTargsVsWeight(pathPref, fing, strin, 
	                  filt, skipParse = TRUE)
}

main()