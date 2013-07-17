F2 (5pos) Low Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/5varpos/F2/low/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F2_low_5pos_all.rda"
lab2 <- "../rsessions/F2_low_5pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/5varpos/F2/low/protein_seqs_JSD/all.txt"
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
## 1 AAG RDNRS 0.0005985 27 27 0.2922
## 2 AAG RYHRR 0.0005406 27 27 0.4881
## 3 AAG RDNRT 0.0005404 18 18 0.3104
## 4 AAG RDNRL 0.0004840 27 27 0.3634
## 5 AAG RDNSS 0.0004752 27 27 0.1713
## 6 AAG RDNTS 0.0004720 18 18 0.1483
```

```r
length(unique(data$V2))
```

```
## [1] 1170937
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
## [1] "AAG:  68368 (0.3457137222579)"
## [1] "CAA:  56632 (0.343786802646755)"
## [1] "CAC:  39187 (0.228886669353473)"
## [1] "CAG:  66430 (0.293789858213114)"
## [1] "CAT:  202183 (0.338301756405601)"
## [1] "GAG:  75274 (0.318035862162208)"
## [1] "GGA:  29592 (0.200601968600016)"
## [1] "GGG:  28087 (0.200953001023117)"
## [1] "TAG:  39612 (0.220725163404156)"
```

```r

nrow(cut10data)
```

```
## [1] 605365
```

```r
nrow(data)
```

```
## [1] 2060883
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.2937
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

