F2 (5pos) High Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/5varpos/F2/high/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F2_high_5pos_all.rda"
lab2 <- "../rsessions/F2_high_5pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/5varpos/F2/high/protein_seqs_JSD/all.txt"
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
##    V1    V2        V3 V4 V5     V6
## 1 AAG ASNRS 0.0011818 50 54 0.3141
## 2 AAG ASNRA 0.0009738 36 36 0.2364
## 3 AAG TWNRS 0.0008412 18 18 0.2103
## 4 AAG RDNRS 0.0008148 27 27 0.3421
## 5 AAG RDNSS 0.0006801 27 27 0.1656
## 6 AAG ASNRN 0.0006675 18 18 0.1585
```

```r
length(unique(data$V2))
```

```
## [1] 791482
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
## [1] "AAG:  47149 (0.504267379679144)"
## [1] "CAA:  56225 (0.299059604480708)"
## [1] "CAC:  37943 (0.300364937501484)"
## [1] "CAG:  11336 (0.113541666666667)"
## [1] "CAT:  34118 (0.226660023251952)"
## [1] "GAG:  33147 (0.210511942791457)"
## [1] "GGA:  48439 (0.274019641119634)"
## [1] "GGG:  50914 (0.456537723498503)"
## [1] "TAG:  53798 (0.342902670660973)"
```

```r

nrow(cut10data)
```

```
## [1] 373069
```

```r
nrow(data)
```

```
## [1] 1260837
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.2959
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

