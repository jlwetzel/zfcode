F3 High Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/6varpos/F3/high/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F3_high_6pos_all.rda"
lab2 <- "../rsessions/F3_high_6pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/6varpos/F3/high/protein_seqs_JSD/all.txt"
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
##    V1     V2        V3 V4  V5     V6
## 1 AAA WRSWLA 0.0020935 23  54 0.6408
## 2 AAA WRSWRA 0.0015897 22  54 0.6809
## 3 AAA WPSWRA 0.0012864 12  36 0.5553
## 4 AAA WRSWRN 0.0012477 13  27 0.7929
## 5 AAA TSVSRV 0.0010848 49 216 0.7172
## 6 AAA WPSWRN 0.0009133  8  18 0.5086
```

```r
length(unique(data$V2))
```

```
## [1] 2290439
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
## [1] "AAA:  20906 (0.287340049754663)"
## [1] "AAC:  25557 (0.37157603954638)"
## [1] "AAG:  5017 (0.259988599264134)"
## [1] "AAT:  2137 (0.0803715822332544)"
## [1] "ACA:  13636 (0.233161773506831)"
## [1] "ACC:  320 (0.0844104457926668)"
## [1] "ACG:  4791 (0.156896777574011)"
## [1] "ACT:  23514 (0.398184681557245)"
## [1] "AGA:  271 (0.088015589477103)"
## [1] "AGC:  9811 (0.180180345631852)"
## [1] "AGG:  19280 (0.339969318121705)"
## [1] "AGT:  5887 (0.24769638574494)"
## [1] "ATA:  20524 (0.38570972167409)"
## [1] "ATC:  5382 (0.225018814282131)"
## [1] "ATG:  12116 (0.28399981247949)"
## [1] "ATT:  18773 (0.415727351241225)"
## [1] "CAA:  25836 (0.35996321787834)"
## [1] "CAC:  23279 (0.384078534895232)"
## [1] "CAG:  2675 (0.112451656297293)"
## [1] "CAT:  2952 (0.164484314927286)"
## [1] "CCA:  26055 (0.347961377689339)"
## [1] "CCC:  2494 (0.0846456692913386)"
## [1] "CCG:  3089 (0.208702114721978)"
## [1] "CCT:  29249 (0.389213429320417)"
## [1] "CGA:  42231 (0.277738683222297)"
## [1] "CGC:  5105 (0.206170994709422)"
## [1] "CGG:  561 (0.183333333333333)"
## [1] "CGT:  2345 (0.0686837326460079)"
## [1] "CTA:  11223 (0.287614361497655)"
## [1] "CTC:  30284 (0.372030171248864)"
## [1] "CTG:  16332 (0.412299303241442)"
## [1] "CTT:  19168 (0.379293968656008)"
## [1] "GAA:  2501 (0.0836175192243397)"
## [1] "GAC:  2129 (0.0632088355798349)"
## [1] "GAG:  4017 (0.147002854424358)"
## [1] "GAT:  3618 (0.156908665105386)"
## [1] "GCA:  3766 (0.136005778259299)"
## [1] "GCC:  2597 (0.0978117585025046)"
## [1] "GCG:  3549 (0.167232117613797)"
## [1] "GCT:  4031 (0.182762060210374)"
## [1] "GGA:  3205 (0.203724891939995)"
## [1] "GGC:  18809 (0.328868917524872)"
## [1] "GGG:  5038 (0.210900870730074)"
## [1] "GGT:  4997 (0.234821428571429)"
## [1] "GTA:  2571 (0.227784176486223)"
## [1] "GTC:  3632 (0.209011912297865)"
## [1] "GTG:  4800 (0.253177910227333)"
## [1] "GTT:  6325 (0.252223152689716)"
## [1] "TAA:  66346 (0.406417308846771)"
## [1] "TAC:  19404 (0.34006309148265)"
## [1] "TAG:  2212 (0.0981845621199343)"
## [1] "TAT:  4122 (0.190507001894902)"
## [1] "TCA:  49122 (0.250791349276043)"
## [1] "TCC:  2687 (0.175094487162779)"
## [1] "TCG:  18277 (0.372748964982767)"
## [1] "TCT:  2542 (0.15017427778106)"
## [1] "TGA:  3634 (0.225728306105969)"
## [1] "TGC:  16164 (0.291316728544137)"
## [1] "TGG:  632 (0.192272588986918)"
## [1] "TGT:  14232 (0.400337552742616)"
## [1] "TTA:  16645 (0.206063682281866)"
## [1] "TTC:  21307 (0.431499220316329)"
## [1] "TTG:  15673 (0.269221519857084)"
## [1] "TTT:  68535 (0.399374147757071)"
```

```r

nrow(cut10data)
```

```
## [1] 827919
```

```r
nrow(data)
```

```
## [1] 2849574
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.2905
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

