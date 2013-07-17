F1 Low Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/6varpos/F1/low/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F1_low_6pos_all.rda"
lab2 <- "../rsessions/F1_low_6pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/6varpos/F1/low/protein_seqs_JSD/all.txt"
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
##    V1     V2       V3 V4 V5      V6
## 1 AAA DLRNAR 0.001553 18 54 0.79510
## 2 AAA YPSWEA 0.001523  8 12 0.27620
## 3 AAA DLRNGR 0.001429 19 54 0.70752
## 4 AAA DARNAR 0.001182  9 36 0.91539
## 5 AAA TKEYLY 0.001135  6  6 0.09209
## 6 AAA YPSWKA 0.001038  5 12 0.36664
```

```r
length(unique(data$V2))
```

```
## [1] 2459511
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
## [1] "AAA:  8820 (0.22388627998477)"
## [1] "AAC:  800 (0.225415610030995)"
## [1] "AAG:  3681 (0.310370994940978)"
## [1] "AAT:  2988 (0.173297761280594)"
## [1] "ACA:  2194 (0.144275662523838)"
## [1] "ACC:  1129 (0.185874217978268)"
## [1] "ACG:  5198 (0.19773280584297)"
## [1] "ACT:  1008 (0.225503355704698)"
## [1] "AGA:  4745 (0.138750804140593)"
## [1] "AGC:  1823 (0.246918596776378)"
## [1] "AGG:  1529 (0.0857495373226403)"
## [1] "AGT:  9923 (0.158032202067175)"
## [1] "ATA:  1586 (0.0618082618862042)"
## [1] "ATC:  2129 (0.159798844104181)"
## [1] "ATG:  2032 (0.111875791444145)"
## [1] "ATT:  4883 (0.162177422033279)"
## [1] "CAA:  1756 (0.168505901544957)"
## [1] "CAC:  1149 (0.213370473537604)"
## [1] "CAG:  1649 (0.23446608844021)"
## [1] "CAT:  2710 (0.169820779546309)"
## [1] "CCA:  2754 (0.112398987837728)"
## [1] "CCC:  2323 (0.145724860422809)"
## [1] "CCG:  2418 (0.134191686553083)"
## [1] "CCT:  2024 (0.145517290962686)"
## [1] "CGA:  7788 (0.155822328931573)"
## [1] "CGC:  1638 (0.171662125340599)"
## [1] "CGG:  2526 (0.172282089755831)"
## [1] "CGT:  4812 (0.0180357116085216)"
## [1] "CTA:  3477 (0.0125694082942912)"
## [1] "CTC:  11549 (0.247630687422274)"
## [1] "CTG:  1638 (0.146459227467811)"
## [1] "CTT:  1919 (0.171768707482993)"
## [1] "GAA:  1624 (0.145454545454545)"
## [1] "GAC:  7570 (0.0534982332155477)"
## [1] "GAG:  2303 (0.149245026245869)"
## [1] "GAT:  7389 (0.0392701878208739)"
## [1] "GCA:  2141 (0.133221330346587)"
## [1] "GCC:  703 (0.120956641431521)"
## [1] "GCG:  1965 (0.1629623486482)"
## [1] "GCT:  2175 (0.123138764649267)"
## [1] "GGA:  4577 (0.240439167892414)"
## [1] "GGC:  10913 (0.257491387853334)"
## [1] "GGG:  2528 (0.195167142746854)"
## [1] "GGT:  10503 (0.164559341950646)"
## [1] "GTA:  1849 (0.144917313269065)"
## [1] "GTC:  1692 (0.086760332273613)"
## [1] "GTG:  1685 (0.144896379740304)"
## [1] "GTT:  1184 (0.154146595495378)"
## [1] "TAA:  13232 (0.245491651205937)"
## [1] "TAC:  3597 (0.00763872666652509)"
## [1] "TAG:  11008 (0.127218934911243)"
## [1] "TAT:  7203 (0.0622773646896075)"
## [1] "TCA:  1481 (0.113713144963145)"
## [1] "TCC:  4652 (0.149721605355476)"
## [1] "TCG:  802 (0.1335998667333)"
## [1] "TCT:  6986 (0.269417662938681)"
## [1] "TGA:  11083 (0.312135635226857)"
## [1] "TGC:  11066 (0.333393588816582)"
## [1] "TGG:  1362 (0.127123389957066)"
## [1] "TGT:  6803 (0.146505868418219)"
## [1] "TTA:  2136 (0.165286698135108)"
## [1] "TTC:  9783 (0.37041384271705)"
## [1] "TTG:  15231 (0.329696733554127)"
## [1] "TTT:  10420 (0.249389689339907)"
```

```r

nrow(cut10data)
```

```
## [1] 284244
```

```r
nrow(data)
```

```
## [1] 2785109
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.1021
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

