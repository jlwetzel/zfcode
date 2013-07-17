F2 Low Stringency Analysis
========================================================

## Load libraries


```r
library(ggplot2)
```


## Import the data and basic stuff

```r
tableFile <- "../data/b1hdata/newDatabase/6varpos/F2/low/protein_seqs_JSD/all_cut10.txt"
lab1 <- "../rsessions/F2_low_6pos_all.rda"
lab2 <- "../rsessions/F2_low_6pos_cut10.rda"
fname <- "../data/b1hData/newDatabase/6varpos/F2/low/protein_seqs_JSD/all.txt"
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
##    V1     V2       V3 V4  V5     V6
## 1 AAA SAGSSN 0.003159 56 108 0.3854
## 2 AAA SPGSYN 0.002467 31  36 0.3231
## 3 AAA SAGSYN 0.002302 29  36 0.1941
## 4 AAA SSGSYN 0.002201 37  54 0.3117
## 5 AAA SPGSHN 0.001985 26  36 0.3360
## 6 AAA SAGSHN 0.001937 25  36 0.3109
```

```r
length(unique(data$V2))
```

```
## [1] 1987116
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
## [1] "AAA:  3421 (0.115863984285037)"
## [1] "AAC:  6154 (0.103751159066004)"
## [1] "AAG:  5881 (0.166501514679652)"
## [1] "AAT:  5309 (0.100081060191905)"
## [1] "ACA:  5603 (0.112824953182578)"
## [1] "ACC:  3934 (0.116500829187396)"
## [1] "ACG:  6985 (0.114765949756009)"
## [1] "ACT:  9640 (0.148861916674388)"
## [1] "AGA:  4961 (0.136261261261261)"
## [1] "AGC:  2742 (0.121531779097598)"
## [1] "AGG:  5535 (0.129364745477493)"
## [1] "AGT:  8297 (0.145812097993041)"
## [1] "ATA:  2636 (0.0684337599626159)"
## [1] "ATC:  4638 (0.103494443700629)"
## [1] "ATG:  7360 (0.152678089864332)"
## [1] "ATT:  7031 (0.111470471660721)"
## [1] "CAA:  2366 (0.0888705254854825)"
## [1] "CAC:  5218 (0.101175010664288)"
## [1] "CAG:  7293 (0.120529516758115)"
## [1] "CAT:  5045 (0.0994794336869503)"
## [1] "CCA:  3199 (0.138359067514381)"
## [1] "CCC:  3718 (0.128477141573655)"
## [1] "CCG:  3001 (0.13077962260862)"
## [1] "CCT:  4319 (0.139655952919873)"
## [1] "CGA:  3871 (0.160769166874325)"
## [1] "CGC:  4684 (0.156655518394649)"
## [1] "CGG:  7554 (0.129974707066536)"
## [1] "CGT:  4743 (0.145889083694749)"
## [1] "CTA:  5245 (0.139039843066564)"
## [1] "CTC:  2427 (0.115897044076214)"
## [1] "CTG:  4787 (0.0960434974519482)"
## [1] "CTT:  3135 (0.0944163353812794)"
## [1] "GAA:  2864 (0.117410732587218)"
## [1] "GAC:  4265 (0.103963533541342)"
## [1] "GAG:  5583 (0.153387548766416)"
## [1] "GAT:  4685 (0.0912543825477211)"
## [1] "GCA:  3854 (0.103443648173498)"
## [1] "GCC:  3827 (0.110047158960202)"
## [1] "GCG:  6997 (0.124728154301401)"
## [1] "GCT:  6329 (0.117593504394196)"
## [1] "GGA:  21057 (0.25192922005671)"
## [1] "GGC:  3811 (0.101012510602205)"
## [1] "GGG:  6689 (0.139674253497599)"
## [1] "GGT:  4897 (0.114461351471379)"
## [1] "GTA:  11768 (0.39195310418332)"
## [1] "GTC:  4106 (0.122570823009642)"
## [1] "GTG:  8606 (0.162389614310514)"
## [1] "GTT:  7212 (0.12097423510467)"
## [1] "TAA:  5236 (0.229719650769973)"
## [1] "TAC:  3665 (0.164180441696905)"
## [1] "TAG:  2227 (0.14867481140263)"
## [1] "TAT:  6263 (0.149417883385819)"
## [1] "TCA:  4579 (0.242596026490066)"
## [1] "TCC:  3136 (0.165820642978003)"
## [1] "TCG:  4642 (0.167357681075819)"
## [1] "TCT:  4619 (0.20134257443006)"
## [1] "TGA:  3278 (0.0943526567267285)"
## [1] "TGC:  1216 (0.155220832269594)"
## [1] "TGG:  7529 (0.144952927359889)"
## [1] "TGT:  9170 (0.129954792171535)"
## [1] "TTA:  6039 (0.135720064724919)"
## [1] "TTC:  2880 (0.13334568015557)"
## [1] "TTG:  5192 (0.123754588358679)"
## [1] "TTT:  5072 (0.233431516936672)"
```

```r

nrow(cut10data)
```

```
## [1] 342025
```

```r
nrow(data)
```

```
## [1] 2529759
```

```r
nrow(cut10data)/nrow(data)
```

```
## [1] 0.1352
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

