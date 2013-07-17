F1 High Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/6varpos/F1/high/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F1_high_6pos_all.rda"
lab2 <- "../rsessions/F1_high_6pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/6varpos/F1/high/protein_seqs_JSD/all.txt"
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
##    V1     V2        V3 V4 V5     V6
## 1 AAA DLRNAR 0.0011165 14 54 0.8520
## 2 AAA YLSLTQ 0.0011130 13 54 0.7262
## 3 AAA DARNAR 0.0010914 12 36 0.8880
## 4 AAA DLRNGR 0.0010402 15 54 0.6922
## 5 AAA DQRNWR 0.0008602  9  9 0.5023
## 6 AAA HPSHTL 0.0007108  7 36 0.6563
```

```r
length(unique(data$V2))
```

```
## [1] 1920137
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
## [1] "AAA:  2585 (0.0592360044913953)"
## [1] "AAC:  948 (0.0472041029726634)"
## [1] "AAG:  1015 (0.0497134740657295)"
## [1] "AAT:  661 (0.0552398462309878)"
## [1] "ACA:  1139 (0.0794226344048532)"
## [1] "ACC:  2553 (0.0414017903476907)"
## [1] "ACG:  17424 (0.293852770048065)"
## [1] "ACT:  1389 (0.043864081349081)"
## [1] "AGA:  2061 (0.0182434585561024)"
## [1] "AGC:  649 (0.0193529148650663)"
## [1] "AGG:  2122 (0.115728621291449)"
## [1] "AGT:  4301 (0.130337282947968)"
## [1] "ATA:  1563 (0.0938569627094217)"
## [1] "ATC:  1983 (0.139441670768582)"
## [1] "ATG:  2013 (0.133523481029451)"
## [1] "ATT:  2473 (0.0951080686101069)"
## [1] "CAA:  1015 (0.063592506735167)"
## [1] "CAC:  874 (0.0775716694772344)"
## [1] "CAG:  1704 (0.113358169238957)"
## [1] "CAT:  2068 (0.141537198001506)"
## [1] "CCA:  2535 (0.142800811176205)"
## [1] "CCC:  1708 (0.141931194947648)"
## [1] "CCG:  1697 (0.12363397930934)"
## [1] "CCT:  1921 (0.122559652928416)"
## [1] "CGA:  2203 (0.0743980277599541)"
## [1] "CGC:  1473 (0.136972289380696)"
## [1] "CGG:  2496 (0.172960986764604)"
## [1] "CGT:  55266 (0.348463735584714)"
## [1] "CTA:  3231 (0.0227588101460198)"
## [1] "CTC:  2990 (0.141706161137441)"
## [1] "CTG:  961 (0.0800833333333333)"
## [1] "CTT:  10329 (0.327977645826057)"
## [1] "GAA:  2024 (0.121255691349149)"
## [1] "GAC:  8947 (0.133740919011032)"
## [1] "GAG:  1855 (0.116366601844301)"
## [1] "GAT:  7250 (0.104829381145171)"
## [1] "GCA:  1487 (0.122882406412693)"
## [1] "GCC:  10624 (0.173807770961145)"
## [1] "GCG:  1031 (0.0624772754817598)"
## [1] "GCT:  2041 (0.136832931080719)"
## [1] "GGA:  2786 (0.0266246177370031)"
## [1] "GGC:  1779 (0.084252900781435)"
## [1] "GGG:  2164 (0.147281018171919)"
## [1] "GGT:  2149 (0.0758292166549047)"
## [1] "GTA:  1520 (0.109234638878908)"
## [1] "GTC:  678 (0.0346980552712385)"
## [1] "GTG:  2508 (0.178721584835744)"
## [1] "GTT:  745 (0.0495906277041869)"
## [1] "TAA:  2703 (0.0564843064320642)"
## [1] "TAC:  8056 (0.121899919802684)"
## [1] "TAG:  444 (0.00827246981666418)"
## [1] "TAT:  20288 (0.26509519018437)"
## [1] "TCA:  1817 (0.107279919702427)"
## [1] "TCC:  676 (0.0364971385379549)"
## [1] "TCG:  993 (0.0755592756049308)"
## [1] "TCT:  3182 (0.0998368473895582)"
## [1] "TGA:  3528 (0.0954055004191568)"
## [1] "TGC:  1509 (0.106297548605241)"
## [1] "TGG:  1417 (0.0812779626018126)"
## [1] "TGT:  2127 (0.0295498749652681)"
## [1] "TTA:  4255 (0.0923815106710958)"
## [1] "TTC:  3646 (0.118592245641426)"
## [1] "TTG:  3668 (0.068464769015399)"
## [1] "TTT:  3195 (0.112432698736672)"
```

```r

nrow(cut10data)
```

```
## [1] 248442
```

```r
nrow(data)
```

```
## [1] 2197512
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.1131
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

