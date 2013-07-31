proteinDir = paste0('../data/b1hData/newDatabase/6varpos/F2/',
                    'low/protein_seq_cut3bc_025/')
data <- read.csv(paste0(proteinDir, 'all.csv'))
minFreq <- min(data$freq)
data$freq <- data$freq*exp(1)/minFreq
data$freq <- log(data$freq)
olmNoIntercept <- lm(freq ~ n1*a6 + n2*a3 + 
                     n3*a2 + n3*a0 + 0, data = data)
