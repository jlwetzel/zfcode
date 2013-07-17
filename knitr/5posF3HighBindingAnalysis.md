F3 (5pos) High Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/5varpos/F3/high/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F3_high_5pos_all.rda"
lab2 <- "../rsessions/F3_high_5pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/5varpos/F3/high/protein_seqs_JSD/all.txt"
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
##    V1    V2       V3 V4 V5     V6
## 1 ACG DSNLR 0.005959 25 27 0.2720
## 2 ACG DSNGR 0.003181 17 18 0.3042
## 3 ACG DSNSR 0.003170 19 27 0.5277
## 4 ACG DTNRR 0.002689 18 18 0.2964
## 5 ACG DSDRR 0.002587 22 27 0.4441
## 6 ACG DSNVR 0.002471 14 18 0.2672
```

```r
length(unique(data$V2))
```

```
## [1] 14312
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
## [1] "ACG:  456 (0.0450548364786088)"
## [1] "AGC:  174 (0.0655861289106672)"
## [1] "CGT:  27 (0.0229982964224872)"
## [1] "TGT:  11 (0.0108910891089109)"
```

```r

nrow(cut10data)
```

```
## [1] 668
```

```r
nrow(data)
```

```
## [1] 14958
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.04466
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

