

require(data.table)
require(dplyr)


##### TOP 5% FOR EACH SPECIES TO FILTER MERGED MATRIX #####
# 0.903909 is top 5% FOR YRI
# 1.185418 is top 5% for mus
# 1.481161 is top 5% for chimp
# 1.682461 is top 5% for gorilla


############################## IMPORT YRI DATA #######################
### YRI data with elevated CN in YRI removed (see Bal_sel_Feb26.R)
#YRI.bed.old <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_bn10/YRIrmcnv_unlifted.bed")
# all unique, YRI.test <- data.frame(unique(setDT(YRI.bed), by=c("chr", "start", "stop")))

## filtered for HWe and SB
YRI.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/YRIrmcnv_unlifted.bed")


#test <- anti_join(YRI.bed,YRI.bed.old)
#dim(test) #70704 rows

#some TajDs are slightly different
#YRI.bed[5320796,]
#YRI.bed.old[5320796,]

YRI.ordered <- YRI.bed[with(YRI.bed, order(TajD)), ]
(nrow(YRI.ordered) - ceiling(nrow(YRI.ordered)*0.01)) 
YRI.ordered[2637599,] #1.71 is top 1% (used to be with 2kb windows that 1.74 was top 5%) (1.6 is merged chr1)
(nrow(YRI.ordered) - ceiling(nrow(YRI.ordered)*0.02)) #top 2%
YRI.ordered[2610957,] #1.37 is top 2% 
(nrow(YRI.ordered) - ceiling(nrow(YRI.ordered)*0.05)) #top 5% (is 0.86 in merged with mus for chr1)
YRI.ordered[2531029,] #0.903909 is top 5% 

#YRI.mus.ordered <- YRI.mus.overlap[with(YRI.mus.overlap, order(s1.TajD)), ]
#(nrow(YRI.mus.ordered) - ceiling(nrow(YRI.mus.ordered)*0.02))
#YRI.mus.ordered[278975,]
###1.5 for human

#YRI.mus.ordered <- YRI.mus.overlap[with(YRI.mus.overlap, order(s2.TajD)), ]
#(nrow(YRI.mus.ordered) - ceiling(nrow(YRI.mus.ordered)*0.02))
#YRI.mus.ordered[277082,]
### 1.5 for mus

############################## IMPORT MUS DATA #######################
mus.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/mus_filteredliftchr.bed")
mus.bed <- mus.bed[,-4]
colnames(mus.bed) <- c("chr", "start", "stop", "TajD")
mus.bed$chr <- gsub("chr", "", mus.bed$chr)

mus.bed$length <-  mus.bed$stop - mus.bed$start
## drop windows that lifted over to greater than 10kb
mus.bed <- mus.bed[-which(mus.bed$length >8000),] 
mus.bed <- mus.bed[,-5]


mus.ordered <- mus.bed[with(mus.bed, order(TajD)), ]
mus.ordered$length <-  mus.ordered$stop - mus.ordered$start
## drop windows that lifted over to greater than 4kb 
dim(mus.ordered) #1692258

(nrow(mus.ordered) - ceiling(nrow(mus.ordered)*0.01)) 
mus.ordered[1675335,] #(1.7 top 1%, with 2kb win 1.63)

(nrow(mus.ordered) - ceiling(nrow(mus.ordered)*0.02)) 
mus.ordered[1712467,] #1.5 is top 2%

(nrow(mus.ordered) - ceiling(nrow(mus.ordered)*0.05)) 
mus.ordered[1607645,] #1.185 is top 5%


## check unique
#mustest.bed <- data.frame(unique(setDT(mus.bed), by=c("chr", "start", "stop")))




############################## IMPORT CHIMP DATA #######################
chimp.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/chimp_ALLchr_lifted1.bed")# remove fourth column
#chimp.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/chimp_allchr_lifted1.bed")
# remove fourth column
chimp.bed <- chimp.bed[,-4]
colnames(chimp.bed) <- c("chr", "start", "stop", "TajD")
chimp.bed$chr <- gsub("chr", "", chimp.bed$chr)

####chimp problem,some duplicate window with different TajD... must both lift to same spot
2:  16 32916851 32918851  0.651225    32917851
3:  16 32916851 32918851  0.316821    32917851
### get rid of duplicate rows
#chimp.bed <- data.frame(unique(setDT(chimp.bed), by=c("chr", "start", "stop")))

### also some end up with the same midpoint
2:  10 19422317 19424263  0.456862    19423290
3:  10 19422289 19424290 -0.355175    19423290

chimp.bed$length <-  chimp.bed$stop - chimp.bed$start
chimp.bed <- chimp.bed[-which(chimp.bed$length > 8000),1:4] #used 8kb max

chimp.ordered <- chimp.bed[with(chimp.bed, order(TajD)), ]

(nrow(chimp.ordered) - ceiling(nrow(chimp.ordered)*0.01)) 
chimp.ordered[2601224,] #2.1

(nrow(chimp.ordered) - ceiling(nrow(chimp.ordered)*0.02)) 
chimp.ordered[2351580,] #1.85 is top 2%

(nrow(chimp.ordered) - ceiling(nrow(chimp.ordered)*0.05)) 
chimp.ordered[2496124,] #1.48 is top 5% 1.5 is top 5% for 2kb windows


##################### IMPORT GORILLA DATA ###################

gor.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/liftedtohg19_largechr_gorilla.bed")

gor.bed <- gor.bed[,-4]
colnames(gor.bed) <- c("chr", "start", "stop", "TajD")
gor.bed$chr <- gsub("chr", "", gor.bed$chr)

gor.bed$length <-  gor.bed$stop - gor.bed$start
## drop windows that lifted over to greater than 10kb
gor.bed <- gor.bed[-which(gor.bed$length >8000),] 
gor.bed <- gor.bed[,-5]


gor.ordered <- gor.bed[with(gor.bed, order(TajD)), ]

(nrow(gor.ordered) - ceiling(nrow(gor.ordered)*0.01)) 
gor.ordered[1675335,]

(nrow(gor.ordered) - ceiling(nrow(gor.ordered)*0.02)) 
gor.ordered[1712467,] 

(nrow(gor.ordered) - ceiling(nrow(gor.ordered)*0.05)) 
gor.ordered[2234362,] #1.682461 is top 5%









### plot unmerged data, not overlapping okay
yri.1 <- intersect(which(YRI.bed$chr==1), which(YRI.bed$start>=17977000))
yri.2 <- intersect(yri.1, which(YRI.bed$stop<=18153000))
yri.p <- YRI.bed[yri.2,]
yri.p$midpoint <- ((yri.p$stop + yri.p$start)/2)
yri.p$species <- 'Human'


chimp.1 <- intersect(which(chimp.bed$chr==1), which(chimp.bed$start>=17977000))
chimp.2 <- intersect(chimp.1, which(chimp.bed$stop<=18153000))
chimp.p <- chimp.bed[chimp.2,]
chimp.p$midpoint <- ((chimp.p$stop + chimp.p$start)/2)
chimp.p$species <- 'Chimp'

mus.1 <- intersect(which(mus.bed$chr==1), which(mus.bed$start>=17977000))
mus.2 <- intersect(mus.1, which(mus.bed$stop<=18153000))
mus.p <- mus.bed[mus.2,]
mus.p$midpoint <- ((mus.p$stop + mus.p$start)/2)
mus.p$species <- 'Mouse'


toplot <- rbind(yri.p, chimp.p, mus.p)

ggplot(toplot) + geom_line(aes(x=midpoint, y=TajD, color=species)) +
  ggtitle("Windows plotted as midpoints") +
  labs(x= 'Chr 1 Position', y='Tajima\'s D') + theme_bw()






######################### GET DATA OVERLAP FOR EACH PAIR OF SPECIES ##################################

## results from May 3rd
#HUMAN MOUSE OVERLAP
YRI.mus.overlap <- get.overlap.data(YRI.bed, mus.bed)
YRI.mus.bin <- bin.widen(YRI.mus.overlap, 1.7, 1.7) #top 1%
cor.test(as.numeric(YRI.mus.bin$s1.TajD), as.numeric(YRI.mus.bin$s2.TajD)) #c=.001, p=0.008

YRI.mus.bin <- bin.widen(YRI.mus.overlap, 0.9, 1.2) #top 5%
cor.test(as.numeric(YRI.mus.bin$s1.TajD), as.numeric(YRI.mus.bin$s2.TajD)) #c=.0069, p=2x10-16, limit of R

### HUMAN CHIMP OVERLAP
YRI.chimp.overlap <- get.overlap.data(YRI.bed, chimp.bed)
YRI.chimp.bin <- bin.widen(YRI.chimp.overlap, 1.17, 1.46) #top 1%
cor.test(as.numeric(YRI.chimp.bin$s1.TajD), as.numeric(YRI.chimp.bin$s2.TajD)) #c=.0068, p=2.2x10-16, limit of R

YRI.chimp.bin <- bin.widen(YRI.chimp.overlap, 0.90, 1.5) #top 5%
cor.test(as.numeric(YRI.chimp.bin$s1.TajD), as.numeric(YRI.chimp.bin$s2.TajD)) #c=.0075, p=2.2x10-16, limit of R

###CHIMP MOUSE OVERLAP
chimp.mus.overlap <- get.overlap.data(chimp.bed, mus.bed)
chimp.mus.bin <- bin.widen(chimp.mus.overlap, 1.17, 1.46) #top 1%
cor.test(as.numeric(YRI.chimp.bin$s1.TajD), as.numeric(YRI.chimp.bin$s2.TajD)) #c=.0068, p=2.2x10-16, limit of R

chimp.mus.bin <- bin.widen(chimp.mus.overlap, 1.5, 1.2) #top 5%
cor.test(as.numeric(chimp.mus.bin$s1.TajD), as.numeric(chimp.mus.bin$s2.TajD)) #

setkey(chimp.mus.bin, "chr", "s2.midpoint")
chimp.mus.bin <- unique(chimp.mus.bin) 
dim(chimp.mus.bin) # was 1843679


#### Do stats with hypergeometric ###
## chimp-human, genes overlapping top .1% of regions
#7177 windows overlap
#129729 positive human windows
#2459457 = total number of genes - 129729 (human genes)
#394 is number of chimp genes

### chimp-human
phyper(7177, 129729, 2459457, 124954, lower.tail =  FALSE) #3.853582e-33

###human-mus
phyper(4592, 84646, 1830608, 86983, lower.tail =  FALSE) #5.019186e-35

###chimp-mus
phyper(4321, 85422, 1843679, 87444, lower.tail =  FALSE) #5.162557e-14


### Do stats with Fisher's exact test ###
#The Fisher's exact test forms this problem slightly different but its calculation is also based on the hypergeometric distribution. It starts by constructing a contingency table:
matrix(c(n - union(A,B), setdiff(A,B), setdiff(B,A), intersect(A,B)), nrow=2)

#matrix(c(n - union(A,B), setdiff(A,B), setdiff(B,A), intersect(A,B)), nrow=2)
fisher.test(matrix(c(110, 20, 60, 10), nrow=2), alternative="greater")

### chimp-human
fisher.test(matrix(c(2341680, 122552, 117777, 7177), nrow=2)) p-value < 2.2e-16, still same prob
  129729+124954 - 7177 = union




### prep input file for putting three species together 
## I should go back and do this better later
## change species 2 to species 3
setnames(YRI.mus.overlap,2:4,c("s3.start", "s3.stop", "s3.TajD"))
setnames(YRI.mus.overlap,9,"s3.midpoint")

#remove s1 from species one col names
setnames(YRI.mus.overlap,6:8,c("start", "stop", "TajD"))

## remove midpoint/dist columns
YRI.mus.overlap[ ,`:=`(midpoint = NULL, s1.midpoint = NULL, dist = NULL)]


### add chimp data onto YRI and mus merged data
YRI.mus.chimp.overlap <- get.overlap.data.addon(data.frame(YRI.mus.overlap), chimp.bed)


## how often are mus/chimp points reused?
dim(YRI.mus.chimp.overlap)
setkey(YRI.mus.chimp.overlap, "chr", "s3.midpoint")
dim(unique(YRI.mus.chimp.overlap)) #73% of mus points are unique
setkey(YRI.mus.chimp.overlap, "chr", "s2.midpoint")
dim(unique(YRI.mus.chimp.overlap)) #97% of chimp points are unique
setkey(YRI.mus.chimp.overlap, "chr", "s1.midpoint")
dim(unique(YRI.mus.chimp.overlap)) #100% of YRI points are unique


### Get rid of duplicate values for mus and chimp, keep 71% of data ##
dim(YRI.mus.chimp.overlap)
#[1] 1819658      15
setkey(YRI.mus.chimp.overlap, "chr", "s3.midpoint")
YRI.mus.chimp.overlap <- unique(YRI.mus.chimp.overlap) 
setkey(YRI.mus.chimp.overlap, "chr", "s2.midpoint")
YRI.mus.chimp.overlap <- unique(YRI.mus.chimp.overlap) 
dim(YRI.mus.chimp.overlap)
#[1] 1294989      15


### plot merged data with their own midpoints

### get three different dataframes
plot.YRI <- data.frame(YRI.mus.chimp.overlap)[7960:8035, c("s1.midpoint", "s1.TajD" )]
plot.YRI$species <- 'Human'
colnames(plot.YRI)[1:2] <- c("pos", "TajD")

plot.chimp <- data.frame(YRI.mus.chimp.overlap)[7960:8035, c("s2.midpoint", "s2.TajD" )]
plot.chimp$species <- 'Chimp'
colnames(plot.chimp)[1:2] <- c("pos", "TajD")

plot.mus <- data.frame(YRI.mus.chimp.overlap)[7960:8035, c("s3.midpoint", "s3.TajD" )]
plot.mus$species <- 'Mouse'
colnames(plot.mus)[1:2] <- c("pos", "TajD")

toplot <- rbind(plot.YRI, plot.chimp, plot.mus)

ggplot(toplot) + geom_line(aes(x=pos, y=TajD, color=species)) +
  ggtitle("Windows plotted as midpoints") +
  labs(x= 'Chr 1 Position', y='Tajima\'s D') + theme_bw()





### plot merged data with human midpoints
### get three different dataframes
plot.YRI <- data.frame(YRI.mus.chimp.overlap)[7960:8035, c("s1.midpoint", "s1.TajD" )]
plot.YRI$species <- 'Human'
colnames(plot.YRI)[1:2] <- c("pos", "TajD")

plot.chimp <- data.frame(YRI.mus.chimp.overlap)[7960:8035, c("s1.midpoint", "s2.TajD" )]
plot.chimp$species <- 'Chimp'
colnames(plot.chimp)[1:2] <- c("pos", "TajD")

plot.mus <- data.frame(YRI.mus.chimp.overlap)[7960:8035, c("s1.midpoint", "s3.TajD" )]
plot.mus$species <- 'Mouse'
colnames(plot.mus)[1:2] <- c("pos", "TajD")

toplot <- rbind(plot.YRI, plot.chimp, plot.mus)

ggplot(toplot) + geom_line(aes(x=pos, y=TajD, color=species)) +
  ggtitle("Points with data from all species") +
  labs(x= 'Chr 1 Position', y='Tajima\'s D') + theme_bw()




#### find top 5% cutoffs for the merged data table
my.cutoff.row <- round(nrow(YRI.mus.chimp.overlap)*0.05)
YRI.ordered <- YRI.mus.chimp.overlap[order(-rank(s1.TajD))]
YRI.cut <- YRI.ordered[my.cutoff.row,s1.TajD]

chimp.ordered <- YRI.mus.chimp.overlap[order(-rank(s2.TajD))]
chimp.cut <- chimp.ordered[my.cutoff.row,s2.TajD]

mus.ordered <- YRI.mus.chimp.overlap[order(-rank(s3.TajD))]
mus.cut <- mus.ordered[my.cutoff.row,s3.TajD]


### get binary cutoffs for three species ##
#YRI.mus.chimp.overlap <- get.overlap.data(data.frame(YRI.mus.overlap), chimp.bed)

YRI.mus.chimp.bin <- cutoffs.3(YRI.mus.chimp.overlap, YRI.cut, chimp.cut, mus.cut) #top 5% human, chimp, mus

nrow(YRI.mus.chimp.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1,])

#227/293 = 77.4% FDR when keep duplicate values for mus and chimp

############### With three species (no duplicate values), I have a 75% FDR
# > YRI.cut
# [1] 0.848862
# > chimp.cut
# [1] 1.460143
# > mus.cut
# [1] 1.172498
# > YRI.mus.chimp.bin <- cutoffs.3(YRI.mus.chimp.overlap, YRI.cut, chimp.cut, mus.cut) #top 5% human, chimp, mus
# > nrow(YRI.mus.chimp.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1,])
# [1] 215
# > dim(YRI.mus.chimp.bin)
# [1] 1294989      15
# > .05*.05*.05*1294989
# [1] 161.8736
# > 162/215
# [1] 0.7534884

#119100:119500 okay


#data for plot below
data.frame(YRI.mus.chimp.bin)[7960:8035,c("s1.TajD", "s2.TajD", "s3.TajD", "s1.midpoint", "s2.midpoint", "s3.midpoint")]

example <- data.frame(YRI.mus.chimp.bin)[7960:8035,c("chr","s1.TajD", "s2.TajD", "s3.TajD", "s1.midpoint", "s2.midpoint", "s3.midpoint")]
#colnames(example) <- c("chr", "Human", "Chimp", "Mouse", "midpoint")
minus <- rep(2000, nrow(example))
plus <- rep(1900, nrow(example))
example$win.start <- (example$s1.midpoint - minus)
example$win.stop <- (example$s1.midpoint + plus)

example <- example[c(1,5,8,9,6,7,2:4)]

colnames(example) <- c("chr", "midpoint", "win.start", "win.stop", "Chimp.midpoint", "Mouse.midpoint", "Human", "Chimp", "Mouse")

write.table(example, "/Users/boettger/postdoc/bal_sel/TajD/plots/example_data_verbose.txt", quote=F, row.names=F, sep = "\t")

#### plot binary data #### row 20395 has all three pos TajD
### get three different dataframes
plot.YRI <- data.frame(YRI.mus.chimp.bin)[7960:8035, c("s1.midpoint", "s1.TajD" )]
plot.YRI$species <- 'Human'
colnames(plot.YRI)[1:2] <- c("pos", "TajD")

plot.chimp <- data.frame(YRI.mus.chimp.bin)[7960:8035, c("s1.midpoint", "s2.TajD" )]
plot.chimp$species <- 'Chimp'
colnames(plot.chimp)[1:2] <- c("pos", "TajD")

plot.mus <- data.frame(YRI.mus.chimp.bin)[7960:8035, c("s1.midpoint", "s3.TajD" )]
plot.mus$species <- 'Mouse'
colnames(plot.mus)[1:2] <- c("pos", "TajD")

toplot <- rbind(plot.YRI, plot.chimp, plot.mus)

ggplot(toplot) + geom_line(aes(x=pos, y=TajD, color=species)) +
  ggtitle("Binary calls using top 5% Tajima's D cutoff for each species") +
  labs(x= 'Chr 1 Position', y='Tajima\'s D') + theme_bw()

####### import positions of SNPs in trans-species haplotypes #######
ts.haps <- read.table("/Users/boettger/postdoc/bal_sel/molly_positions_hg19.txt", header=T)
colnames(ts.haps) <- c("chr", "midpoint")
ts.haps$SNP_coord <- ts.haps$midpoint

coding <- read.table("/Users/boettger/postdoc/bal_sel/molly_shared_coding.txt", header=T)
colnames(coding) <- c("chr", "midpoint")
coding$SNP_coord <- coding$midpoint
coding$coding_chr <- coding$chr

# remove 'X' chr
coding <- coding[-336,]

cutoff.05.all3 <- data.frame(YRI.chimp.gor.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1,])[c(1,2,3,5,6,7,9,10)]

## convert to dt
cutoff.05.dt <- data.table(cutoff.05.all3)
coding.dt <- data.table(coding)
ts.haps.dt <- data.table(ts.haps)
  
### make sure chromosomes are numeric (species with X/Y etc will be strings)
# cutoff.05.dt[, chr := as.numeric(chr)] <- already numeric
coding.dt[, chr := as.numeric(as.character(chr))]
coding.dt[, chr := as.numeric(as.character(coding_chr))]
  
### set keys as chr and midpoint for both s1 and s2
setkey(cutoff.05.dt, "chr", "midpoint")
setkey(coding.dt, "chr", "midpoint")
setkey(ts.haps.dt, "chr", "midpoint")
  
  
cutoff.05.coding <-  coding.dt[cutoff.05.dt, roll= "nearest"]
cutoff.05.tshaps <-  ts.haps.dt[cutoff.05.dt, roll= "nearest"]


#find distance between midpoints
cutoff.05.coding <- cbind(cutoff.05.coding, cutoff.05.coding[,.(dist = abs(midpoint-SNP_coord))])
cutoff.05.tshaps <- cbind(cutoff.05.tshaps, cutoff.05.tshaps[,.(dist = abs(midpoint-SNP_coord))])

MAXdist = 100000
cutoff.05.coding[dist<=MAXdist]
cutoff.05.tshaps[dist<=MAXdist]

### didn't merge right, these SNPs are not on the right chromosomes!!! 
 chr midpoint SNP_coord s2.start  s2.stop s3.start  s3.stop s1.start  s1.stop dist
1:   6 77737000  77743797 77735219 77739223 77736574 77737569 77735000 77739000 6797
2:   6 77741000  77743797 77739223 77743225 77740589 77741776 77739000 77743000 2797
3:   6 77742000  77743797 77740224 77744226 77740589 77742502 77740000 77744000 1797
4:  11 11512000  11516007 11510011 11514016 11510182 11514002 11510000 11514000 4007
5:  11 11513000  11516007 11511011 11514143 11511232 11516021 11511000 11515000 3007

#MHC Location: chr6:28,477,797-33,448,354







###### Spatial cross-correlation for a series of binary variables ##########

http://stats.stackexchange.com/questions/168157/spatial-cross-correlation-what-is-the-correlation-statistic-used-by-the-correlo


## Mantel test?
http://www.ats.ucla.edu/stat/r/faq/mantel_test.htm

## Mornan's I
http://www.ats.ucla.edu/stat/r/faq/morans_i.htm


#http://gis.stackexchange.com/questions/49430/what-is-an-appropriate-statistic-to-measure-spatial-autocorrelation-of-points-wi
#"Binary data is a normal use case for spatial autocorrelation."
#Point Pattern Analysis" techniques. Specifically "Ripley's K"

#"1 dimensional" "spatial cross correlation" +binary

### I think I want ccm from the MTS package

################ SIMULATIONS ###############

#dim(YRI.chimp.overlap)
#[1] 2589173      11

## how many positive s1?
dim(YRI.chimp.overlap[s1.TajD==1])
#129728

## how many positive s1?
dim(YRI.chimp.overlap[s2.TajD==1])
#124952

test.data <- cbind(rep(0, nrow(YRI.chimp.overlap)), rep(0, nrow(YRI.chimp.overlap)))
s1.pos <- sample(seq(from = 1, to = nrow(YRI.chimp.overlap)), size = dim(YRI.chimp.overlap[s1.TajD==1]), replace = F)
s2.pos <- sample(seq(from = 1, to = nrow(YRI.chimp.overlap)), size = dim(YRI.chimp.overlap[s2.TajD==1]), replace = F)

test.data[s1.pos,1] <- 1
test.data[s2.pos,2] <- 1
#cor.test(test.data[,1], test.data[,2]) #p-value = 0.8047

#dim(YRI.chimp.overlap[s2.TajD==1 & s1.TajD==1])
#[1] 7177   11

length(intersect(which(test.data[,1]==1), which(test.data[,2]==1)))
[1] 6242, 6282, 

#7177-6242
#[1] 935
#> 935/7177
#[1] 0.1302773

## we do 13% better than random chance
#but this has a big difference in p-value
cor.test(as.numeric(YRI.chimp.bin$s1.TajD), as.numeric(YRI.chimp.bin$s2.TajD))
#p-value < 2.2e-16
#0.007570071

#
### should look for overlap of three

#DBG 6/20/16 transition probabilities 3 species

nrow(YRI.mus.chimp.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1,])
#215

#hum:
table(YRI.mus.chimp.bin$s1.TajD[2:nrow(YRI.mus.chimp.bin)] - YRI.mus.chimp.bin$s1.TajD[1:nrow(YRI.mus.chimp.bin)-1])
#      -1       0       1 
#   26373 1242243   26373 
table(YRI.mus.chimp.bin$s1.TajD[2:nrow(YRI.mus.chimp.bin)] + YRI.mus.chimp.bin$s1.TajD[1:nrow(YRI.mus.chimp.bin)-1])
#       0       1       2 
# 1203868   52746   38375 
#mus
table(YRI.mus.chimp.bin$s2.TajD[2:nrow(YRI.mus.chimp.bin)] - YRI.mus.chimp.bin$s2.TajD[1:nrow(YRI.mus.chimp.bin)-1])
#      -1       0       1 
#   28190 1238609   28190 
table(YRI.mus.chimp.bin$s2.TajD[2:nrow(YRI.mus.chimp.bin)] + YRI.mus.chimp.bin$s2.TajD[1:nrow(YRI.mus.chimp.bin)-1])
#       0       1       2 
# 1202050   56380   36559
#chi
table(YRI.mus.chimp.bin$s3.TajD[2:nrow(YRI.mus.chimp.bin)] - YRI.mus.chimp.bin$s3.TajD[1:nrow(YRI.mus.chimp.bin)-1])
#      -1       0       1 
#   25597 1243795   25597 
table(YRI.mus.chimp.bin$s3.TajD[2:nrow(YRI.mus.chimp.bin)] + YRI.mus.chimp.bin$s3.TajD[1:nrow(YRI.mus.chimp.bin)-1])
#       0       1       2 
# 1204643   51194   39152 

####### three species results 10,000 simulations #######

#look at 10,000 tests of chimp-human MCMC
hum_mus_chimp_10000 <- data.frame(read.table("/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/three_species.txt", header=F))
#make numeric
hum_mus_chimp_10000 <- as.numeric(hum_mus_chimp_10000[1:10000,])
hist(hum_mus_chimp_10000)

z <- (nrow(YRI.mus.chimp.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1,]) - mean(hum_mus_chimp_10000)) / sd(hum_mus_chimp_10000)
pvalue2sided=2*pnorm(-abs(z))
8.885272e-18









########################## FUNCTIONS ###################
### this is extremely slow for chimp-human
species1.bed <- data.frame(YRI.mus.overlap)
#species1.bed <- YRI.bed
species2.bed <- chimp.bed 

unique.merge <- function(dist){
  if (.N==2){
    min(dist)
  }
}


YRI.mus.overlap[1:10,print(.SD),by=c("chr","s2.midpoint")]

### fuck this works, but I'm not sure why
### how does it know to take the min of dist and not some other col???
YRI.mus.overlap.out = YRI.mus.overlap[, lapply(.SD,min), by=c("chr","s2.midpoint")]


## fuck, this doesn't work
YRI.mus.overlap[, .SD[lapply(min(dist))], by=c("chr","s2.midpoint")]
YRI.mus.overlap[, .SD[unique.merge(dist)], by=c("chr","s2.midpoint")]


get.overlap.data.addon <- function(species1.bed, species2.bed, MAXdist = 2000){

  ### get midpoints of windows
  species1.bed$s1.midpoint <- round((species1.bed$stop + species1.bed$start)/2)
  colnames(species1.bed)[c(5,6,7,9)] <- c("s1.start", "s1.stop", "s1.TajD", "midpoint")

  species2.bed$s2.midpoint <- round((species2.bed$stop + species2.bed$start)/2)
  colnames(species2.bed) <- c("chr", "s2.start", "s2.stop", "s2.TajD", "midpoint")
  
  #get rid of duplicate midpoints
  species2.bed <- data.frame(unique(setDT(species2.bed), by=c("chr", "midpoint")))

  ### convert to data.table
  s1.dt <- data.table(species1.bed)
  s2.dt <- data.table(species2.bed)
  
  ### keep only chromosomes that s1 has
  my.chr <- intersect(s2.dt$chr, as.character(unique(s1.dt$chr)))
  s2.dt <- s2.dt[chr %in% my.chr,]
  
  ### make sure chromosomes are numeric (species with X/Y etc will be strings)
  s1.dt[, chr := as.numeric(chr)]
  s2.dt[, chr := as.numeric(chr)]
  
  ### set keys as chr and midpoint for both s1 and s2
  setkey(s1.dt, "chr", "midpoint")
  setkey(s2.dt, "chr", "midpoint")
  
  #s1.s2 <- s2.dt[ s1.dt , list( chr, s1.TajD, s1.midpoint, s2.TajD, s2.start, s2.stop, s2.midpoint) , roll = "nearest", by=chr ]
  
  s1.s2 <-  s2.dt[s1.dt, roll= "nearest"]
  ## chr 1 isn't prsent???
  ## s1 midpoint not present?
  #s1.s2[chr==1,]
  
  #recalculate actual midpoints (ereased/unclear because merged on non-exact value), looks like s1.midpoint is maintained
  s1.s2 <- cbind(s1.s2, s1.s2[, .(s2.midpoint = round((s2.start + s2.stop)/2), s1.midpoint = (s1.start + s1.stop)/2)])
  
  #find distance between midpoints
  s1.s2 <- cbind(s1.s2, s1.s2[,.(dist = abs(s1.midpoint-s2.midpoint))])
  
  s1.s2.keep <- s1.s2[dist<=MAXdist]
  return(s1.s2.keep)
}

species1.bed <- chimp.bed
species2.bed <- mus.bed
get.overlap.data <- function(species1.bed, species2.bed, MAXdist = 2000){

  ### get midpoints of windows
  species1.bed$s1.midpoint <- round((species1.bed$stop + species1.bed$start)/2)
  colnames(species1.bed) <- c("chr", "s1.start", "s1.stop", "s1.TajD", "midpoint")

  species2.bed$s2.midpoint <- round((species2.bed$stop + species2.bed$start)/2)
  colnames(species2.bed) <- c("chr", "s2.start", "s2.stop", "s2.TajD", "midpoint")
  
  #get rid of duplicate midpoints
  species1.bed <- data.frame(unique(setDT(species1.bed), by=c("chr", "midpoint")))
  species2.bed <- data.frame(unique(setDT(species2.bed), by=c("chr", "midpoint")))

  ### convert to data.table
  s1.dt <- data.table(species1.bed)
  s2.dt <- data.table(species2.bed)

  
  ### keep only chromosomes that s1 has
  #my.chr <- intersect(s2.dt$chr, as.character(unique(s1.dt$chr)))
  my.chr <- c(1:22)
  s1.dt <- s1.dt[chr %in% my.chr,]
  s2.dt <- s2.dt[chr %in% my.chr,]

  ### make sure chromosomes are numeric (species with X/Y etc will be strings)
  s1.dt[, chr := as.numeric(as.character(chr))]
  s2.dt[, chr := as.numeric(as.character(chr))]
  
  ### set keys as chr and midpoint for both s1 and s2
  setkey(s1.dt, "chr", "midpoint")
  setkey(s2.dt, "chr", "midpoint")
  
  #s1.s2 <- s2.dt[ s1.dt , list( chr, s1.TajD, s1.midpoint, s2.TajD, s2.start, s2.stop, s2.midpoint) , roll = "nearest", by=chr ]
  
  s1.s2 <-  s2.dt[s1.dt, roll= "nearest"]
  ## chr 1 isn't prsent???
  ## s1 midpoint not present?
  #s1.s2[chr==1,]
  
  #recalculate actual midpoints (ereased/unclear because merged on non-exact value), looks like s1.midpoint is maintained
  s1.s2 <- cbind(s1.s2, s1.s2[, .(s2.midpoint = round((s2.start + s2.stop)/2), s1.midpoint = (s1.start + s1.stop)/2)])
  
  #find distance between midpoints, keep only those below cutoff
  s1.s2 <- data.table(cbind(s1.s2, s1.s2[,.(dist = abs(s1.midpoint-s2.midpoint))]))
  
  s1.s2.keep <- s1.s2[dist<=MAXdist,]
  
  ### keep closer of windows from s2, so no duplcate windows 
  s1.s2.return <- s1.s2.keep[s1.s2.keep[, .I[which.min(dist)], by=c("chr","s2.midpoint")]$V1]
  
  return(s1.s2.return)
}


cutoffs.3 <- function(overlap.data, s1.cutoff, s2.cutoff, s3.cutoff){
  overlap.data[s1.TajD <= s1.cutoff, s1.TajD:= 0]
  overlap.data[s2.TajD <= s2.cutoff, s2.TajD:= 0]
  overlap.data[s3.TajD <= s3.cutoff, s3.TajD:= 0]
  overlap.data[s1.TajD > s1.cutoff, s1.TajD:= 1]
  overlap.data[s2.TajD > s2.cutoff, s2.TajD:= 1]
  overlap.data[s3.TajD > s3.cutoff, s3.TajD:= 1]
  
  return(overlap.data)
}



s1.cutoff <- 1
s2.cutoff <- 1.5
## works line by line, but not as function :(
# I think it isn't returnign the right thing...?
#64000 makes 66% of human windows positive
#20000 makes 30% of windows positive
#10kb makes 17% positive
bin.widen <- function(overlap.data, s1.cutoff, s2.cutoff, cutoff_dist=10000){
  overlap.data[s1.TajD <= s1.cutoff, s1.TajD:= 0]
  overlap.data[s2.TajD <= s2.cutoff, s2.TajD:= 0]
  overlap.data[s1.TajD > s1.cutoff, s1.TajD:= 1]
  overlap.data[s2.TajD > s2.cutoff, s2.TajD:= 1]
  
  ## make new column for nearby midpoints
  #s1.nearTajD <- rep(0,nrow(overlap.data))
  #overlap.data <- cbind(overlap.data, s1.nearTajD)
  #overlap.data[s1.TajD==1, s1.nearTajD:=1] #set actual high TajD windows to 1
  
  #overlap.data[s1.TajD == 1,
  #update_nearness(.BY$chr, .BY$s1.midpoint, cutoff_dist),
  #by=c('s1.midpoint','chr')]
  
  return(overlap.data)
}

# update is nearness finds all the points in overlap.data that are on the same
# chr as the input and are less then cutoff_dist away from near_midpoint and
# sets their s1.nearTajD to 1 and sets their near.midpoint to this.
update_nearness <- function(near_chr, near_midpoint, cutoff_dist) {
  cat(near_chr, ', ', near_midpoint, '\n') #print location to know where we are
  overlap.data[
    chr == near_chr &
    s1.TajD == 0 &
    abs(as.numeric(s1.midpoint) - as.numeric(near_midpoint)) < cutoff_dist,
    c('s1.nearTajD', 's1.near.midpoint') := list(1, near_midpoint)]
  # we return NA here because otherwise the code above will return the whole data.table
  # which we don't want to do in that loop above.
  return(NA)
}


##old
bin.widen <- function(overlap.data, s1.cutoff, s2.cutoff, cutoff_dist=64000){
  overlap.data[which(as.numeric(overlap.data$s1.TajD) < s1.cutoff),"s1.TajD"] <- 0
  overlap.data[which(as.numeric(overlap.data$s2.TajD) < s2.cutoff ),"s2.TajD"] <- 0
  overlap.data[which(as.numeric(overlap.data$s1.TajD) >= s1.cutoff),"s1.TajD"] <- 1
  overlap.data[which(as.numeric(overlap.data$s2.TajD) >= s2.cutoff ),"s2.TajD"] <- 1
  
  overlap.data$s1.nearTajD <- 0
  
  overlap.data <- data.table(overlap.data)
  
  ##!!!!WARNING, THIS CURRENTLY LEAVES ACTUAL +TAJD LOCATIONS AS 0
  overlap.data[s1.TajD == 1,
  update_nearness(.BY$chr, .BY$s1.midpoint, cutoff_dist),
  by=c('s1.midpoint','chr')]
  
  return(overlap.data)
}




#Pearson/Phi coefficient
cor.test(as.numeric(overlap.data$s1.nearTajD), as.numeric(overlap.data$s2.TajD))

bin.mat <- data.frame(cbind(as.numeric(overlap.data$s1.nearTajD), as.numeric(overlap.data$s2.TajD)))








### method for binary comp
### I think I actually want resemblance coefficients
https://www.researchgate.net/post/Can_anyone_help_with_a_correlation_coefficient_between_two_binary_variables

### binary pearson is phi coefficient
https://www.biostars.org/p/80580/


library(simba)
abis.jacc <- sim(bin.mat, method="jaccard")

overlap.data[s1.TajD == 1,
  update_nearness(.BY$chr, .BY$s1.midpoint, cutoff_dist),
  by=c('s1.midpoint','chr')]









species1.bed <- YRI.bed
species2.bed <- mus.bed





























species1.bed <- data.frame(YRI.bed)
species2.bed <- chimp.bed
######pld
get.overlap.data <- function(species1.bed, species2.bed, MAXdist = 1000){
  my.chr <- intersect(species1.bed$chr, species1.bed$chr)
  
  ### get midpoints
  species1.bed$s1.midpoint <- round((species1.bed$stop + species1.bed$start)/2)
  colnames(species1.bed) <- c("chr", "s1.start", "s1.stop", "s1.TajD", "s1.midpoint")

  species2.bed$s2.midpoint <- round((species2.bed$stop + species2.bed$start)/2)
  colnames(species2.bed) <- c("chr", "s2.start", "s2.stop", "s2.TajD", "s2.midpoint")

  my.output <- ""
  
  for (i in my.chr) {
    s1.chr <- ""
    s2.chr <- ""
    
    s1.chr <- species1.bed[which(species1.bed[,1] == i),]
    s1.dt <- data.table(s1.chr, key = "s1.midpoint")
    s2.chr <- species2.bed[which(species2.bed[,1] == i),]
    s2.dt <- data.table(s2.chr, key = "s2.midpoint")
    s2.dt <- unique(s2.dt, by=c("chr", "s2.midpoint"))

    ### combine
    s1.s2.chr <- s2.dt[ s1.dt , list( chr, s1.TajD, s1.midpoint, s2.TajD, s2.start, s2.stop, s2.midpoint) , roll = "nearest" ]
    s1.s2.chr <- data.frame(s1.s2.chr)
    
    ###re-calc actual mus midpoint
    s1.s2.chr$s2.midpoint <- round((s1.s2.chr$s2.stop + s1.s2.chr$s2.start)/2)
  
    ## get distance between midpoints
    s1.s2.chr$dist <- abs(s1.s2.chr$s1.midpoint - s1.s2.chr$s2.midpoint)
    # remove those that are greater than MAXdist away and take only TajD values
    s1.s2.chr <- s1.s2.chr[which(s1.s2.chr$dist < MAXdist),]
    
    my.output <- rbind(my.output, s1.s2.chr)
  }
  my.output <- my.output[-1,]
  return(my.output)
}


#get binary for tajima's D cutoffs for both species
#widen signals of species1 by a given distance
#windows within 54kb + 2kb to account for midpoint
#overlap.data <- YRI.mus.overlap



### Dan's code
# set cutoff dist
#cutoff_dist <- 56e3
# make a new column for 'nearness'
#overlap.data$s1.nearTajD <- 0
# if you don't run the above before updated nearTajD, then just do this:
# overlap.data[is.na(s1.nearTajD), s1.nearTajD := 0]









######### old stuff below #######

# this code runs update_nearness() for every set of rows with a unique 
# combination of chromosome and s1.midpoint (which happens to be only one row)
# using the by=... argument.
# Update_nearness is passed the .BY chr and s1.midpoint.
overlap.data[s1.TajD == 1,
  update_nearness(.BY$chr, .BY$s1.midpoint, cutoff_dist),
  by=c('s1.midpoint','chr')]


hist(overlap.data[s1.TajD == 0 & s1.nearTajD == 1, as.numeric(s1.near.midpoint) - as.numeric(s1.midpoint)])

# this:
sum(as.numeric(overlap.data$s1.nearTajD))
# should be way bigger than this:
sum(as.numeric(overlap.data$s1.TajD))


#####
(TajD=1)    (Taj=0)
10          1    <---- 9
26          5    <---- 4
30          9    <---- 1
45          40   <---- 5
53          10   <---- 0


1 R
5 R
9 R
10 R
10 L
26 R
30 L
40 L
45 R
53 L









60          45



old crap:
  
binary.distance <- function(overlap.data, s1.cutoff, s2.cutoff, dist=56000){
  overlap.data[which(overlap.data$s1.TajD < s1.cutoff),"s1.TajD"] <- 0
  overlap.data[which(overlap.data$s2.TajD < s2.cutoff ),"s2.TajD"] <- 0
  overlap.data[which(overlap.data$s1.TajD >= s1.cutoff),"s1.TajD"] <- 1
  overlap.data[which(overlap.data$s2.TajD >= s2.cutoff ),"s2.TajD"] <- 1
  
  my.chr <- unique(overlap.data$chr)
  my.output <- ""
  
  for (j in my.chr){
    this.chr <- overlap.data[which(overlap.data[,1]==j),]
    cat("chr", this.chr[1,1])
    
    for (i in which(this.chr$s1.TajD ==1)){
      


ddply(df,.(Season),myFun,k = 1)

      setDT(df)[s1.TajD == 1, expand, by=chr]
      
      cat("row", i)
      distances <- abs(as.numeric(this.chr$s1.midpoint) - as.numeric(this.chr[i,"s1.midpoint"]))
      near <- which(distances < dist) 
      if (length(near) > 0){
        this.chr[near,1] <- 1
        my.output <- rbind(my.output, this.chr)
      }
    }
    
  }
  my.output <- my.output[-1,]
  return(overlap.data)
}






overlap.data.dt <- setDT(overlap.data)
test <- setDT(overlap.data)[s1.TajD == 1, expand(overlap.data, 56000), by=chr]

expand <- function(this.chr,dist){
  dt[, .(b, c)]
 distances <- abs(as.numeric(this.chr$s1.midpoint) - as.numeric(this.chr[i,"s1.midpoint"]))
 near <- which(distances < dist) 
 if (length(near) > 0){
  this.chr[near,1] <- 1
  }
    }


binary.distance <- function(overlap.data, s1.cutoff, s2.cutoff, dist=56000){
  overlap.data[which(overlap.data$s1.TajD < s1.cutoff),"s1.TajD"] <- 0
  overlap.data[which(overlap.data$s2.TajD < s2.cutoff ),"s2.TajD"] <- 0
  overlap.data[which(overlap.data$s1.TajD >= s1.cutoff),"s1.TajD"] <- 1
  overlap.data[which(overlap.data$s2.TajD >= s2.cutoff ),"s2.TajD"] <- 1
  
  my.chr <- unique(overlap.data$chr)
  my.output <- ""
  
  for (j in my.chr){
    this.chr <- overlap.data[which(overlap.data[,1]==j),]
    cat("chr", this.chr[1,1])
    
    for (i in which(this.chr$s1.TajD ==1)){
      cat("row", i)
      distances <- abs(as.numeric(this.chr$s1.midpoint) - as.numeric(this.chr[i,"s1.midpoint"]))
      near <- which(distances < dist) 
      if (length(near) > 0){
        this.chr[near,1] <- 1
        my.output <- rbind(my.output, this.chr)
      }
    }
    
  }
  my.output <- my.output[-1,]
  return(overlap.data)
}





#######dbg 6/18/16
# extract Interval objects from data.table and collapse 
# Intervals into unions (of overlapping intervals) and
# convert starts and ends back into a data.table. Do
# this separately for every `chr` using 'by'.

union_ivls = YRI.chimp.bin[s2.TajD == 1 & s1.TajD == 1, 
    as.data.table(
        interval_union(
           Intervals(cbind(s1.start, s1.stop)))),
    by='chr']
setnames(union_ivls, names(union_ivls), c('chr', 's1.start', 's1.stop'))
setkey(union_ivls, 'chr', 's1.start', 's1.stop')
setkey(YRI.chimp.bin, 'chr', 's1.start', 's1.stop')

# Use the data.table foverlaps function to match the original intervals to their 
# respective union-ed intervals.
union_ivl_merged <- foverlaps(
    YRI.chimp.bin[s2.TajD == 1 & s1.TajD == 1],
    union_ivls,
    type='within')

# offline, check how many unions and how large they are
#union_ivl_merged[, .N,  by=c('chr', 's1.start', 's1.stop')]
#table(union_ivl_merged[, .N,  by=c('chr', 's1.start', 's1.stop')]$N)

# combine the unions:
# s1.midpoint: mean of s1.start and s1.stop
# s2.start: min of s2.start
# s2.stop: max of s2.stop
# s2.midpoint: mean of (min of s2.start) and (max of s2.stop)
# dist is distance between new s1 and s2 midpoints
unions_w_metadata <- union_ivl_merged[, list(
    's1.midpoint'= mean(s1.start, s1.stop),
    's2.start'= min(s2.start),
    's2.stop'= max(s2.stop),
    's2.midpoint'= mean(min(s2.start), max(s2.stop)),
    'dist'=  abs(mean(s1.start, s1.stop) - mean(min(s2.start), max(s2.stop))),
    's1.TajD'=1,
    's2.TajD'=1),
    by=c('chr', 's1.start', 's1.stop')]

YRI.chimp.bin[s2.TajD == 1 & s1.TajD == 1] <- unions_w_metadata

######
write.csv(YRI.chimp.bin[, list(chr, s1.TajD, s2.TajD)], 
  file='/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/hum_chi_tajD_states.csv',
  quote=F, row.names=F)

########################### get human-chimp transition probabilities #################

#human sum:
table(YRI.chimp.bin$s1.TajD[2:nrow(YRI.chimp.bin)] + YRI.chimp.bin$s1.TajD[1:nrow(YRI.chimp.bin)-1])
#      0       1       2 
#2419372   80168   89645 
# 0->0: 2419372
# 1->1: 89645

#human diff:
table(YRI.chimp.bin$s1.TajD[2:nrow(YRI.chimp.bin)] - YRI.chimp.bin$s1.TajD[1:nrow(YRI.chimp.bin)-1])
#     -1       0       1 
#  40084 2509017   40084
# 0->1: 40084
# 1->0: 40084

#chimp sum:
table(YRI.chimp.bin$s2.TajD[2:nrow(YRI.chimp.bin)] + YRI.chimp.bin$s2.TajD[1:nrow(YRI.chimp.bin)-1])
#      0       1       2 
#2421470   85522   82193 
# 0->0: 2421470
# 1->1: 82193

#diff:
table(YRI.chimp.bin$s2.TajD[2:nrow(YRI.chimp.bin)] - YRI.chimp.bin$s2.TajD[1:nrow(YRI.chimp.bin)-1])
#     -1       0       1 
#  42761 2503663   42761
# 0->1: 42761
# 1->0: 42761


#look at 10,000 tests of chimp-human MCMC
hum_chimp_10000 <- data.frame(read.table("/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/test_full.txt", header=F))
#make numeric
hum_chimp_10000 <- as.numeric(hum_chimp_10000[1:10000,])
hist(hum_chimp_10000)

z <- (nrow(YRI.chimp.bin[s1.TajD==1 & s2.TajD==1,]) - mean(hum_chimp_10000)) / sd(hum_chimp_10000)
pvalue2sided=2*pnorm(-abs(z))
8.885272e-18


########################### get human-mus transition probabilities #################

#human sum:
table(YRI.mus.bin$s1.TajD[2:nrow(YRI.mus.bin)] + YRI.mus.bin$s1.TajD[1:nrow(YRI.mus.bin)-1])
#      0       1       2 
#1716784   58354   55469 
# 0->0: 1716784
# 1->1: 55469

#human diff:
table(YRI.mus.bin$s1.TajD[2:nrow(YRI.mus.bin)] - YRI.mus.bin$s1.TajD[1:nrow(YRI.mus.bin)-1])
#     -1       0       1 
#  29177 1772253   29177 
# 0->1: 29177
# 1->0: 29177

#mus sum:
table(YRI.mus.bin$s2.TajD[2:nrow(YRI.mus.bin)] + YRI.mus.bin$s2.TajD[1:nrow(YRI.mus.bin)-1])
#      0       1       2 
#1719148   48953   62506 
# 0->0: 1719148
# 1->1: 62506

#diff:
table(YRI.mus.bin$s2.TajD[2:nrow(YRI.mus.bin)] - YRI.mus.bin$s2.TajD[1:nrow(YRI.mus.bin)-1])
#     -1       0       1 
#  24477 1781654   24476 
# 0->1: 24476
# 1->0: 24476

#look at 10,000 tests of mus-human MCMC
hum_mus_10000 <- data.frame(read.table("/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/human_mus_full.txt", header=F, stringsAsFactors=F))
#make numeric
hum_mus_10000 <- as.numeric(hum_mus_10000[1:10000,])
hist(hum_mus_10000)
 mean(hum_mus_10000)
[1] 4023.408
> sd(hum_mus_10000)
[1] 107.4135

z <- (nrow(YRI.mus.bin[s1.TajD==1 & s2.TajD==1,]) - mean(hum_mus_10000)) / sd(hum_mus_10000)
pvalue2sided=2*pnorm(-abs(z))
1.200043e-07


########################### get chimp-mus transition probabilities #################

#chimp sum:
table(chimp.mus.bin$s1.TajD[2:nrow(chimp.mus.bin)] + chimp.mus.bin$s1.TajD[1:nrow(chimp.mus.bin)-1])
#      0       1       2 
#1725594   65324   52760 
# 0->0: 1725594
# 1->1: 52760

#chimp diff:
table(chimp.mus.bin$s1.TajD[2:nrow(chimp.mus.bin)] - chimp.mus.bin$s1.TajD[1:nrow(chimp.mus.bin)-1])
#      -1       0       1 
#   32662 1778354   32662 
# 0->1: 32662
# 1->0: 32662

#mus sum:
table(chimp.mus.bin$s2.TajD[2:nrow(chimp.mus.bin)] + chimp.mus.bin$s2.TajD[1:nrow(chimp.mus.bin)-1])
#       0       1       2 
# 1731952   48564   63162 
# 0->0: 1731952
# 1->1: 63162

#diff:
table(chimp.mus.bin$s2.TajD[2:nrow(chimp.mus.bin)] - chimp.mus.bin$s2.TajD[1:nrow(chimp.mus.bin)-1])
#      -1       0       1 
#   24282 1795114   24282
# 0->1: 24282
# 1->0: 24282


#look at 10,000 tests of mus-chimpan MCMC
chimp_mus_10000 <- data.frame(read.table("/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/chimp_mus_full.txt", header=F, stringsAsFactors=F))
#make numeric
chimp_mus_10000 <- as.numeric(chimp_mus_10000[1:10000,])
hist(chimp_mus_10000)
 mean(chimp_mus_10000)
[1] 4051.048
sd(chimp_mus_10000)
[1] 104.5665

z <- (nrow(chimp.mus.bin[s1.TajD==1 & s2.TajD==1,]) - mean(chimp_mus_10000)) / sd(chimp_mus_10000)
pvalue2sided=2*pnorm(-abs(z))
9.8e-03
