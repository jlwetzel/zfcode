F2 High Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/6varpos/F2/high/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F2_high_6pos_all.rda"
lab2 <- "../rsessions/F2_high_6pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/6varpos/F2/high/protein_seqs_JSD/all.txt"
data <- read.delim(fname, header = FALSE)
lapply(data[1, ], class)
```

```
## $V1
## [1] "factor"
## 
## $V2
## [1] "factor"
## 
## $V3
## [1] "numeric"
## 
## $V4
## [1] "integer"
## 
## $V5
## [1] "integer"
## 
## $V6
## [1] "numeric"
```

```r
head(data)
```

```
##    V1     V2       V3 V4 V5     V6
## 1 AAA SAGSYN 0.007025 30 36 0.2002
## 2 AAA SPGSYN 0.006014 31 36 0.2458
## 3 AAA SPGSHN 0.004414 31 36 0.2265
## 4 AAA SAGSHN 0.004373 29 36 0.2552
## 5 AAA SAGSFN 0.004272 24 36 0.2496
## 6 AAA SAGAYN 0.003932 20 24 0.1994
```

```r
length(unique(data$V2))
```

```
## [1] 1581220
```


## Drop the 10 lowest counts for each target

```r
sortedTargs <- sort(unique(data$V1))
minCounts <- vector(mode = "numeric", length = length(sortedTargs))

cut10data <- data.frame()
lastrows <- 0
for (i in 1:length(sortedTargs)) {
    minCounts[i] <- tail(sort(unique(subset(data, V1 == sortedTargs[i])$V3), 
        decreasing = TRUE), 10)[1]
    cut10data <- rbind(cut10data, subset(data, V3 > minCounts[i] & V1 == sortedTargs[i]))
    rows <- nrow(cut10data) - lastrows
    lastrows <- nrow(cut10data)
    print(paste0(sortedTargs[i], ":  ", rows, " (", rows/nrow(subset(data, V1 == 
        sortedTargs[i])), ")"))
}
```

```
## [1] "AAA:  1221 (0.0799502357255107)"
## [1] "AAC:  2259 (0.135067264573991)"
## [1] "AAG:  1039 (0.0713451898647257)"
## [1] "AAT:  698 (0.0873045653533458)"
## [1] "ACA:  1888 (0.134607158134892)"
## [1] "ACC:  1101 (0.0911423841059603)"
## [1] "ACG:  2861 (0.100590675761198)"
## [1] "ACT:  3673 (0.156945690723412)"
## [1] "AGA:  3310 (0.255776215130206)"
## [1] "AGC:  779 (0.071956401256235)"
## [1] "AGG:  2205 (0.138574660633484)"
## [1] "AGT:  3223 (0.102571446757049)"
## [1] "ATA:  10850 (0.149067128293903)"
## [1] "ATC:  1999 (0.131929778247096)"
## [1] "ATG:  2828 (0.270337443839021)"
## [1] "ATT:  1518 (0.129114570043378)"
## [1] "CAA:  23159 (0.245935412618008)"
## [1] "CAC:  1516 (0.105681422098292)"
## [1] "CAG:  1822 (0.0905656625907148)"
## [1] "CAT:  2040 (0.112427666023698)"
## [1] "CCA:  769 (0.0162452204406701)"
## [1] "CCC:  765 (0.0554387999130372)"
## [1] "CCG:  14140 (0.211246563881917)"
## [1] "CCT:  1207 (0.0860176738882554)"
## [1] "CGA:  3783 (0.252385082393755)"
## [1] "CGC:  1306 (0.0904244270580904)"
## [1] "CGG:  4890 (0.134606914776481)"
## [1] "CGT:  2021 (0.117863183064093)"
## [1] "CTA:  17197 (0.168355409361019)"
## [1] "CTC:  828 (0.0227891999009165)"
## [1] "CTG:  10796 (0.168872203973096)"
## [1] "CTT:  850 (0.0323612274423209)"
## [1] "GAA:  1433 (0.0869644374317271)"
## [1] "GAC:  3701 (0.111368560423688)"
## [1] "GAG:  2992 (0.170864028325053)"
## [1] "GAT:  1610 (0.0407275302926817)"
## [1] "GCA:  1014 (0.023942764043352)"
## [1] "GCC:  2002 (0.0406563502700946)"
## [1] "GCG:  3103 (0.0958574032312873)"
## [1] "GCT:  4291 (0.111495089123318)"
## [1] "GGA:  1361 (0.0347069924006732)"
## [1] "GGC:  1017 (0.0175233040991092)"
## [1] "GGG:  3374 (0.146320308773147)"
## [1] "GGT:  2266 (0.0616011961397309)"
## [1] "GTA:  3242 (0.234961588636034)"
## [1] "GTC:  1935 (0.0898871185023459)"
## [1] "GTG:  3339 (0.102294660090071)"
## [1] "GTT:  2619 (0.108120381455641)"
## [1] "TAA:  4607 (0.289075735709356)"
## [1] "TAC:  776 (0.077121844563705)"
## [1] "TAG:  587 (0.0196979865771812)"
## [1] "TAT:  2277 (0.0654931400465959)"
## [1] "TCA:  1089 (0.0545946758911114)"
## [1] "TCC:  858 (0.112039697048838)"
## [1] "TCG:  10599 (0.248482006798734)"
## [1] "TCT:  1157 (0.0771333333333333)"
## [1] "TGA:  5038 (0.178925311645417)"
## [1] "TGC:  237 (0.0507060333761232)"
## [1] "TGG:  4705 (0.15575344279661)"
## [1] "TGT:  3442 (0.11659496629518)"
## [1] "TTA:  19613 (0.225641674624084)"
## [1] "TTC:  890 (0.0491495471614756)"
## [1] "TTG:  761 (0.0902621278614636)"
## [1] "TTT:  634 (0.0301832896929303)"
```

```r

nrow(cut10data)
```

```
## [1] 225110
```

```r
nrow(data)
```

```
## [1] 1835562
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.1226
```


## Make plot of JSD for all data

```r
JSDbpAll <- ggplot(data, aes(x = V1, y = V6)) + geom_boxplot() + xlab("Target") + 
    ylab("JSD")
print(JSDbpAll)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


## Make plot of JSD for cut10data

```r
JSDbpCut10 <- ggplot(cut10data, aes(x = V1, y = V6)) + geom_boxplot() + xlab("Target") + 
    ylab("JSD")
print(JSDbpCut10)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


## Save the data

```r
write.table(cut10data, file = tableFile, sep = "\t", quote = FALSE, eol = "\n", 
    row.names = FALSE, col.names = FALSE)
save(data, file = lab1)
save(cut10data, JSDbpAll, JSDbpCut10, file = lab2)
```

