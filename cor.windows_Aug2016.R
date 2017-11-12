## filtered for HWe and SB, 
YRI.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/YRIrmcnv_unlifted.bed")


mus.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/mus_filteredliftchr.bed")
#mus.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_10kbwin/mus_allchr_10kb.bed")
mus.bed <- mus.bed[,-4]
colnames(mus.bed) <- c("chr", "start", "stop", "TajD")
mus.bed$chr <- gsub("chr", "", mus.bed$chr)
mus.bed$length <-  mus.bed$stop - mus.bed$start
## drop windows that lifted over to greater than 20kb
mus.bed <- mus.bed[-which(mus.bed$length >8000),] 
mus.bed <- mus.bed[,-5]

chimp.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/chimp_ALLchr_lifted1.bed")
# remove fourth column
chimp.bed <- chimp.bed[,-4]
colnames(chimp.bed) <- c("chr", "start", "stop", "TajD")
chimp.bed$chr <- gsub("chr", "", chimp.bed$chr)
chimp.bed$length <-  chimp.bed$stop - chimp.bed$start
chimp.bed <- chimp.bed[-which(chimp.bed$length > 8000),1:4] #used 20kb max


gor.bed <- read.table("/Users/boettger/postdoc/bal_sel/TajD/liftovers_Ap2016/4kb/liftedtohg19_largechr_gorilla.bed")
gor.bed <- gor.bed[,-4]
colnames(gor.bed) <- c("chr", "start", "stop", "TajD")
gor.bed$chr <- gsub("chr", "", gor.bed$chr)
gor.bed$length <-  gor.bed$stop - gor.bed$start
## drop windows that lifted over to greater than 20kb
gor.bed <- gor.bed[-which(gor.bed$length > 8000),] 
gor.bed <- gor.bed[,-5]


##############COMBINE DATA ###################

##################### ADD DATA TOGETHER IN LARGE DATA TABLE ##############

### start with YRI and mus
YRI.mus.overlap <- get.overlap.data(YRI.bed, mus.bed)
## make mus s4
setnames(YRI.mus.overlap, c("s2.start", "s2.stop", "s2.TajD", "s2.midpoint"), c("s4.start", "s4.stop", "s4.TajD", "s4.midpoint"))
## 
YRI.mus.chimp.overlap <- get.overlap.data.addon(data.frame(YRI.mus.overlap), chimp.bed)
# set chimp to s3
setnames(YRI.mus.chimp.overlap, c("s2.start", "s2.stop", "s2.TajD", "s2.midpoint"), c("s3.start", "s3.stop", "s3.TajD", "s3.midpoint"))

### add on Gorilla ###
YRI.gor.chimp.mus.overlap <- get.overlap.data.addon(data.frame(YRI.mus.chimp.overlap), gor.bed)

################ find top 5% cutoffs for the merged data table ################
my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
YRI.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s1.TajD))]
YRI.cut <- YRI.ordered[my.cutoff.row,s1.TajD] #0.847376 for YRI, 0.69

my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
gor.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s2.TajD))]
gor.cut <- gor.ordered[my.cutoff.row,s2.TajD] # gorilla 1.624389, 1.48

my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
chimp.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s3.TajD))]
chimp.cut <- chimp.ordered[my.cutoff.row,s3.TajD] # chimp 1.459305, 1.37

my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
mus.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s4.TajD))]
mus.cut <- mus.ordered[my.cutoff.row,s4.TajD] # mus 1.170837, 0.99

YRI.gor.chimp.mus.bin<- cutoffs.4(YRI.gor.chimp.mus.overlap, YRI.cut, gor.cut, chimp.cut, mus.cut)




######### look at primate overlap ############

YRI.gor.overlap <- get.overlap.data(YRI.bed, gor.bed)
## make gor s3
setnames(YRI.gor.overlap, c("s2.start", "s2.stop", "s2.TajD", "s2.midpoint"), c("s3.start", "s3.stop", "s3.TajD", "s3.midpoint"))

### add chimp
YRI.chimp.gor.overlap <- get.overlap.data.addon(data.frame(YRI.gor.overlap), chimp.bed)
### remove duplicate s1.midpoint column
YRI.chimp.gor.overlap <- YRI.chimp.gor.overlap[, c(-16), with=F]

################ find top 5% cutoffs for the merged data table of 4 mammals ################
#remake datatable
YRI.gor.chimp.mus.overlap <- get.overlap.data.addon(data.frame(YRI.mus.chimp.overlap), gor.bed)


my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
YRI.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s1.TajD))]
YRI.cut <- YRI.ordered[my.cutoff.row,s1.TajD] #0.847376 for YRI, 0.69

my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
gor.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s2.TajD))]
gor.cut <- gor.ordered[my.cutoff.row,s2.TajD] # gorilla 1.624389, 1.48

my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
chimp.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s3.TajD))]
chimp.cut <- chimp.ordered[my.cutoff.row,s3.TajD] # chimp 1.459305, 1.37

my.cutoff.row <- round(nrow(YRI.gor.chimp.mus.overlap)*0.05)
mus.ordered <- YRI.gor.chimp.mus.overlap[order(-rank(s4.TajD))]
mus.cut <- mus.ordered[my.cutoff.row,s4.TajD] # mus 1.170837, 0.99

## within 1kb (offset=5000, gets around 15% of windows per species)
YRI.gor.chimp.mus.bin<- data.table(cutoffs.4(YRI.gor.chimp.mus.overlap, YRI.cut, gor.cut, chimp.cut, mus.cut))
YRI.gor.chimp.mus.bin[, c("s1.ext","s2.ext","s3.ext", "s4.ext") := bin.widen.4(.SD, offset=5000), by='chr']
mamtot.pos <- YRI.gor.chimp.mus.bin[, list(chr, s1.start, s1.stop)]
# (also can be done like:)
#ape.pos <- YRI.chimp.gor.bin[, c("chr","s1.start","s1.stop"), with=F]
mamtot.ivl.merged <- intervals_merge(mam.pos)
mamtot.ivl.merged[, length := end - start]
sum(mamtot.ivl.merged$length) #1750038000  
setkey(mamtot.ivl.merged, chr,start,end)

### positions of overlapping chimp windows
mam.overlap <- YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1,]
mamov.pos <- mam.overlap[,list(chr,s1.start,s1.stop)]
colnames(mamov.pos) <- c("chr", "start", "stop")
mam.ivl.merged <- intervals_merge(mamov.pos)
mam.ivl.merged[, chr := as.numeric(chr)]
mam.ivl.merged[, length := end - start]
setkey(mam.ivl.merged, chr,start,end)

############## make venn diagram ################
YRI.gor.chimp.mus.bin

### make unique identifier for each row 
YRI.gor.chimp.mus.bin[,identifier := seq(1:nrow(YRI.gor.chimp.mus.bin))]

hum.rows <- YRI.gor.chimp.mus.bin[s1.ext==1,identifier]
chimp.rows <- YRI.gor.chimp.mus.bin[s3.ext==1,identifier]
gor.rows <- YRI.gor.chimp.mus.bin[s2.ext==1,identifier]
mus.rows <- YRI.gor.chimp.mus.bin[s4.ext==1,identifier]

v1 <- venn.diagram(list(Human=hum.rows, Chimpanzee=chimp.rows, Gorilla=gor.rows, Mouse=mus.rows), 
                   cat.cex=1.5,fontfamily = "sans", 
                   cat.fontfamily= "sans", 
                   lty = 'blank',
                   cat.pos = c(-10, 10, -10, 10),
                   filename=NULL, fill=rainbow(4))
grid.newpage()
grid.draw(v1)


### great ape
YRI.chimp.gor.bin[,s3.midpoint := seq(1:nrow(YRI.chimp.gor.bin))]

hum.rows <- YRI.chimp.gor.bin[s1.TajD==1,s3.midpoint]
chimp.rows <- YRI.chimp.gor.bin[s2.TajD==1,s3.midpoint]
gor.rows <- YRI.chimp.gor.bin[s3.TajD==1,s3.midpoint]

v1 <- venn.diagram(list(Human=hum.rows, Chimpanzee=chimp.rows, Gorilla=gor.rows),cex=1.5, 
                   cat.cex=1.5,fontfamily = "sans", 
                   cat.fontfamily= "sans", 
                   lty = 'blank',
                   cat.pos = c(-30, 28, 180),
                   cross.area=6, filename=NULL, fill=rainbow(3))
grid.newpage()
grid.draw(v1)


##### make unique identifier for each row travel
YRI.mus.bin <- data.table(YRI.mus.bin)
YRI.mus.bin[,label := seq(1:nrow(YRI.mus.bin))]
YRI.mus.bin[, s1.TajD := as.numeric(s1.TajD)]
YRI.mus.bin[, s4.TajD := as.numeric(s4.TajD)]
YRI.mus.bin[, label := as.numeric(label)]

hum.rows <- YRI.mus.bin[s1.TajD==1,label]
mus.rows <- YRI.mus.bin[s4.TajD==1,label]

v1 <- venn.diagram(list(Human=hum.rows, Mouse=mus.rows),
                   cat.cex=1.5,fontfamily = "sans", 
                   cat.fontfamily= "sans", 
                   cex=1.5,
                   cat.pos = c(-20, 20),
                   scaled = FALSE,
                   lty = 'blank',filename=NULL,fill=rainbow(2))
grid.newpage()
grid.draw(v1)


#draw.pairwise.venn

venn.plot <- draw.pairwise.venn(62469, 62469, 3803,cat.cex=NA,fill=rainbow(2), cex=2, scaled = FALSE)
grid.newpage()
grid.draw(venn.plot)






# with offset=7000
dim(YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1]) #1847
## what do we expect by chance?
nrow(YRI.gor.chimp.mus.bin[s1.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s2.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s3.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s4.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin)
#1131 61% FDR

# with offset=9000
dim(YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1]) #3292
## what do we expect by chance?
nrow(YRI.gor.chimp.mus.bin[s1.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s2.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s3.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s4.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin)
#2040 61% FDR



# with offset=5000 HERE IS WHERE I CAN GENERATE NUMBERS TO FILL IN SLIDE 49
dim(YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1])
[1] 849   #(0.07% of region)
## what do we expect by chance?
nrow(YRI.gor.chimp.mus.bin[s1.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s2.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s3.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s4.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin)
[1] 518.8748
FDR:0.61116

## get rid of HLA, about 15% of genome for each
YRI.gor.chimp.mus.bin.noHLA <- YRI.gor.chimp.mus.bin[!(chr==6 & s1.midpoint>=28477797 & s1.midpoint<=33448356)]
dim(YRI.gor.chimp.mus.bin.noHLA[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1])
## lost 7
### get transition probabilities #####
trans.prob.mat(YRI.gor.chimp.mus.bin.noHLA$s1.ext)
trans.prob.mat(YRI.gor.chimp.mus.bin.noHLA$s2.ext)
trans.prob.mat(YRI.gor.chimp.mus.bin.noHLA$s3.ext)
trans.prob.mat(YRI.gor.chimp.mus.bin.noHLA$s4.ext)

##### get genes that overlap at least 3 species,  ###
#human, chimp, gorilla already done
Aoverlap <- YRI.gor.chimp.mus.bin.noHLA[s1.TajD==1 & s2.TajD==1 & s4.TajD==1]
A.pos <-Aoverlap[,list(chr,s1.start,s1.stop)]
Boverlap <- YRI.gor.chimp.mus.bin.noHLA[s2.TajD==1 & s3.TajD==1 & s4.TajD==1]
B.pos <-Boverlap[,list(chr,s1.start,s1.stop)]

colnames(A.pos) <- c("chr", "start", "stop")
colnames(B.pos) <- c("chr", "start", "stop")

colnames(HCnobinnoHLA.pos) <- c("chr", "start", "stop")

####### what genes are here?
mam.overlap <- YRI.gor.chimp.mus.bin.noHLA[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1,]
## add others zorbing
mam.overlap <- data.table(rbind(mam.overlap, Aoverlap, Boverlap, HCnobinnoHLA.overlap))
mamov.pos <- mam.overlap[,list(chr,s1.start,s1.stop)]

colnames(mamov.pos) <- c("chr", "start", "stop")
mam.ivl.merged <- intervals_merge(mamov.pos)
mam.ivl.merged[, chr := as.numeric(chr)]
mam.ivl.merged[, length := end - start]

mam.ivl.merged.regions <- paste(mam.ivl.merged$chr, mam.ivl.merged$start, mam.ivl.merged$end, sep = ":" )
mam.ivl.merged.genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), filters = "chromosomal_region", values = mam.ivl.merged.regions, mart = grch37)
mam.ivl.merged.pcgenes <- mam.ivl.merged.genes[which(mam.ivl.merged.genes[,2]=="protein_coding"),]

write.table(mam.ivl.merged.pcgenes[,1], "/Users/boettger/postdoc/bal_sel/gene_lists/mam.overlap.pcgenes_aug8.txt", row.names=F, col.names=F, quote=F)

all.pos<- data.table(rbind(mamov.pos, A.pos, B.pos, HCnobinnoHLA.pos))
all.ivl.merged <- intervals_merge(all.pos)
all.ivl.merged.regions <- paste(all.ivl.merged$chr, all.ivl.merged$start, all.ivl.merged$end, sep = ":" )
all.ivl.merged.genes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"), filters = "chromosomal_region", values = all.ivl.merged.regions, mart = grch37)

all.ivl.merged.pcgenes <- all.ivl.merged.genes[which(all.ivl.merged.genes[,2]=="protein_coding"),]
length(unique(all.ivl.merged.pcgenes$ensembl_gene_id) ## all
write.table(all.ivl.merged.pcgenes[,1], "/Users/boettger/postdoc/bal_sel/gene_lists/all.overlap.pcgenes_aug8.txt", row.names=F, col.names=F, quote=F)

##from grasp
IRgenes <- read.table("/Users/boettger/postdoc/bal_sel/gene_lists/IR_1e5_genelist.txt")
##30, looks like 11% but may be more

Behavgenes <- read.table("/Users/boettger/postdoc/bal_sel/gene_lists/Behav_1e5_genelist.txt")
#27

Behavgenes_prim <- read.table("/Users/boettger/postdoc/bal_sel/gene_lists/Behav_1e5_genelist_prim.txt")
unique(Behavgenes_prim$V1) #26 of 85, 30%
## only 9 immune related though!!!

aug8.genes <- rbind(ape.overlap.pcgenes, mam.ivl.merged.pcgenes)
aug8.genesu <- which(unique(aug8.genes$ensembl_gene_id)


#with offset=4000 (.04%)
 dim(YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1])
[1] 513  25
## what do we expect by chance?
nrow(YRI.gor.chimp.mus.bin[s1.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s2.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s3.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s4.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin)
[1] 309
FDR: 0.6023392
sum(mam.ivl.merged$length) #1018000 (about 2.6%)
[1] 1018000
sum(nearape.ivl.merged$length) #46026000 (same as above)

## with offset=3000
dim(YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1]) #253
## what do we expect by chance?
nrow(YRI.gor.chimp.mus.bin[s1.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s2.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s3.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s4.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin)
#162 (64%)

## with offset=2000... why is this a different answer than no bin.widen at all?
dim(YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1]) #111
## what do we expect by chance?
nrow(YRI.gor.chimp.mus.bin[s1.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s2.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s3.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin[s4.ext==1])/nrow(YRI.gor.chimp.mus.bin)*nrow(YRI.gor.chimp.mus.bin)
#71 (64%)


#YRI.gor.chimp.mus.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1 & s4.TajD==1] #10 found with overlapping windwos


YRI.gor.chimp.mus.bin.noHLA <- YRI.gor.chimp.mus.bin[!(chr==6 & s1.midpoint>=28477797 & s1.midpoint<=33448356)]
dim(YRI.gor.chimp.mus.bin[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1]) #849
dim(YRI.gor.chimp.mus.bin.noHLA[s1.ext==1 & s2.ext==1 & s3.ext==1 & s4.ext==1]) #842








################ find top 5% cutoffs for the merged data table of great apes ################
### add chimp
YRI.chimp.gor.overlap <- get.overlap.data.addon(data.frame(YRI.gor.overlap), chimp.bed)
### remove duplicate s1.midpoint column
YRI.chimp.gor.overlap <- YRI.chimp.gor.overlap[, c(-16), with=F]


my.cutoff.row <- round(nrow(YRI.chimp.gor.overlap)*0.05)
YRI.ordered <- YRI.chimp.gor.overlap[order(-rank(s1.TajD))]
YRI.cut <- YRI.ordered[my.cutoff.row,s1.TajD] #0.899046 for YRI, now  0.89902

my.cutoff.row <- round(nrow(YRI.chimp.gor.overlap)*0.05)
chimp.ordered <- YRI.chimp.gor.overlap[order(-rank(s2.TajD))]
chimp.cut <- chimp.ordered[my.cutoff.row,s2.TajD] # chimp 1.478248

my.cutoff.row <- round(nrow(YRI.chimp.gor.overlap)*0.05)
gor.ordered <- YRI.chimp.gor.overlap[order(-rank(s3.TajD))]
gor.cut <- gor.ordered[my.cutoff.row,s3.TajD] # gorilla 1.675276

YRI.chimp.gor.bin <- cutoffs.3(YRI.chimp.gor.overlap, YRI.cut, chimp.cut, gor.cut)
YRI.chimp.gor.bin <- data.table(YRI.chimp.gor.bin)

### windows that overlap
YRI.chimp.gor.bin[, c("s1.ext","s2.ext","s3.ext") := bin.widen.3(.SD, offset=1001), by='chr']
ape.pos <- YRI.chimp.gor.bin[, list(chr, s1.start, s1.stop)]
# (also can be done like:)
#ape.pos <- YRI.chimp.gor.bin[, c("chr","s1.start","s1.stop"), with=F]
apetot.ivl.merged <- intervals_merge(ape.pos)
apetot.ivl.merged[, length := end - start]
sum(apetot.ivl.merged$length) #2407055000  
setkey(apetot.ivl.merged, chr,start,end)

### positions of overlapping chimp windows
nearape.overlap <- YRI.chimp.gor.bin[s1.ext==1 & s2.ext==1 & s3.ext==1,]
nearape.pos <- nearape.overlap[,list(chr,s1.start,s1.stop)]
colnames(nearape.pos) <- c("chr", "start", "stop")
nearape.ivl.merged <- intervals_merge(nearape.pos)
nearape.ivl.merged[, chr := as.numeric(chr)]
nearape.ivl.merged[, length := end - start]
setkey(nearape.ivl.merged, chr,start,end)

sum(nearape.ivl.merged$length) #5390000 (2kb), 13245000 (4kb), 18425000 (5kb) 18487000 (5.001kb) 46026000 (9kb)

#no offset
dim(YRI.chimp.gor.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1])
507
## what do we expect by chance?
nrow(YRI.chimp.gor.bin[s1.TajD==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s2.TajD==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s3.TajD==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin)
270 #53% FDR


#offset=1001
dim(YRI.chimp.gor.bin[s1.ext==1 & s2.ext==1 & s3.ext==1])
1276
## what do we expect by chance?
nrow(YRI.chimp.gor.bin[s1.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s2.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s3.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin)
824.8651 77% FDR

#offset=2000
dim(YRI.chimp.gor.bin[s1.ext==1 & s2.ext==1 & s3.ext==1])
2744
## what do we expect by chance?
nrow(YRI.chimp.gor.bin[s1.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s2.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s3.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin)
1949 77% FDR


#offset=5000
dim(YRI.chimp.gor.bin[s1.ext==1 & s2.ext==1 & s3.ext==1])
11805
## what do we expect by chance?
nrow(YRI.chimp.gor.bin[s1.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s2.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin[s3.ext==1])/nrow(YRI.chimp.gor.bin)*nrow(YRI.chimp.gor.bin)
9091.953 77% FDR




### positions of overlapping windows (no widen)
apenobin.overlap <- YRI.chimp.gor.bin[s1.TajD==1 & s2.TajD==1 & s3.TajD==1,] #507 rows
apenobin.pos <- apenobin.overlap[,list(chr,s1.start,s1.stop)]
colnames(apenobin.pos) <- c("chr", "start", "stop")
apenobin.ivl.merged <- intervals_merge(apenobin.pos)
apenobin.ivl.merged[, chr := as.numeric(chr)]
apenobin.ivl.merged[, length := end - start]
setkey(apenobin.ivl.merged, chr,start,end)
sum(apenobin.ivl.merged$length) #1431000

foverlaps(molly.haps.dt, apenobin.ivl.merged, type="any", nomatch=0L) #none
foverlaps(coding.dt, apenobin.ivl.merged, type="any", nomatch=0L) #none


# old:
    chr     start      stop   i.start    i.stop
 1:   1   3306000   3325000   3322097   3322097
 2:   1 200142000 200156000 200143281 200143281
 3:   1 201451000 201476000 201459399 201459399
 4:   1 202080000 202102000 202092580 202092580
 5:   1 223809000 223823000 223816383 223816383
 6:   1 238040000 238050000 238048466 238048466
 7:   2   1487000   1512000   1500437   1500437
 8:   4   6294000   6306000   6303247   6303247
 9:   6  76632000  76656000  76640781  76640781
10:   9  32784000  32806000  32784265  32784265
11:  12  49711000  49744000  49723651  49723651
12:  15  63427000  63439000  63433785  63433785
13:  15  73615000  73629000  73615146  73615146
14:  16  58713000  58728000  58713894  58713894
15:  17   1683000   1704000   1690752   1690752
16:  17  76482000  76499000  76486870  76486870
17:  17  76969000  76974000  76972231  76972231
18:  18  72197000  72206000  72201918  72201918
19:  19   8541000   8571000   8567475   8567475
20:  20  17476000  17508000  17477592  17477592

    chr     start       end length   i.start      stop
 1:   2 129385000 129402000  17000 129385155 129385155
 2:   2 129385000 129402000  17000 129388993 129388993
 3:   3 116436000 116456000  20000 116454070 116454070
 4:   3 116436000 116456000  20000 116454385 116454385
 5:   4  57915000  57922000   7000  57918296  57918296
 6:   4  57915000  57922000   7000  57918427  57918427
 7:   4  57915000  57922000   7000  57918492  57918492
 8:   4  57915000  57922000   7000  57920048  57920048
 9:   4  57915000  57922000   7000  57920087  57920087
10:   5 145348000 145364000  16000 145362862 145362862
11:   5 145348000 145364000  16000 145363859 145363859
12:   7  47791000  47802000  11000  47800478  47800478
13:   7 155041000 155062000  21000 155049958 155049958
14:   7 155041000 155062000  21000 155051811 155051811
15:  11  81491000  81509000  18000  81492293  81492293
16:  11 131548000 131552000   4000 131548195 131548195
17:  18  58417000  58441000  24000  58437878  58437878
18:  18  58417000  58441000  24000  58438910  58438910
19:  22  47641000  47650000   9000  47642100  47642100
20:  22  47641000  47650000   9000  47645285  47645285

# check these by hand
YRI.chimp.gor.overlap[
  s1.start > 3306000-5001 &
  s1.stop < 3325000 + 5001 &
  chr == 1 &
  (s1.TajD | s2.TajD | s3.TajD)]

YRI.chimp.gor.bin[
  s1.start > 3306000-5001 &
  s1.stop < 3325000 + 5001 &
  chr == 1 &
  (s1.TajD | s2.TajD | s3.TajD)]

s1.midpoint: 3317000, 3318000 -> c(3317000-5e3, 3318000+5e3)
s2.midpoint: 3314000 -> c(3314468-5e3, 3314468+5e3)
s3.midpoint: 3324000, 3325000 -> c(3324228-5e3, 3325140+5e3)

3306000 to 3325000
   chr s2.start s2.stop s2.TajD midpoint s3.start s3.stop s3.TajD s1.start s1.stop
1:   1  3312486 3316449       1  3314000  3311657 3315659       0  3312000 3316000
2:   1  3315449 3319364       0  3317000  3314659 3318674       0  3315000 3319000
3:   1  3316449 3320276       0  3318000  3315659 3319713       0  3316000 3320000
4:   1  3322272 3325902       0  3324000  3322165 3326290       1  3322000 3326000
5:   1  3323271 3326680       0  3325000  3323165 3327114       1  3323000 3327000
   s1.TajD s3.midpoint s1.midpoint dist s2.midpoint
1:       0     3313658     3314000  468     3314468
2:       1     3316666     3317000  406     3317406
3:       1     3317686     3318000  362     3318362
4:       0     3324228     3324000   87     3324087
5:       0     3325140     3325000   24     3324976




########################## OVERLAP WITH CODING SNPS ################
coding <- read.table("/Users/boettger/postdoc/bal_sel/molly_shared_coding.txt", header=T)
coding.df <- cbind(coding, coding$coord)
colnames(coding.df) <- c("chr", "start", "end")
coding.dt <- data.table(coding.df)
coding.dt[, chr := as.numeric(chr)]
setkey(coding.dt, chr,start,end)
##how many are possible?
nrow(foverlaps(coding.dt, apetot.ivl.merged, type="any", nomatch=0L))  ##320

##### overlapping ape windows ###
foverlaps(molly.haps.dt, nearape.ivl.merged, type="any", nomatch=0L)
   chr     start       end length   i.start      stop
1:   7 155046000 155056000  10000 155049958 155049958
2:   7 155046000 155056000  10000 155051811 155051811
foverlaps(coding.dt, nearape.ivl.merged, type="any", nomatch=0L)
   chr    start      end length  i.start     stop
1:  18 56245000 56249000   4000 56247744 56247744


# 12:   2 101278037 101278037
# 13:   2 124678599 124678599
# 14:   2 124680876 124680876
# 15:   2 129385155 129385155
# 16:   2 129388993 129388993
# 17:   2 168652708 168652708
# 18:   2 168655935 168655935
# 19:   2 199543380 199543380
# 20:   2 199544091 199544091

### what do we expect??
18425000/2407055000*320
[1] 1.265592

 ############################# Human Chimp overlap ######################
dim(YRI.chimp.gor.bin)
[1] 2162824      15

### windows that overlap
YRI.chimp.gor.bin[, c("s1.ext","s2.ext","s3.ext") := bin.widen.3(.SD, offset=4000), by='chr']
ape.pos <- YRI.chimp.gor.bin[, list(chr, s1.start, s1.stop)]
# (also can be done like:)
#ape.pos <- YRI.chimp.gor.bin[, c("chr","s1.start","s1.stop"), with=F]
apetot.ivl.merged <- intervals_merge(ape.pos)
apetot.ivl.merged[, length := end - start]
sum(apetot.ivl.merged$length) #2407055000  
setkey(apetot.ivl.merged, chr,start,end)

### positions of overlapping chimp windows (bin widen)
HC.overlap <- YRI.chimp.gor.bin[s1.ext==1 & s2.ext==1,] #47863 rows with 4kb offset
HC.pos <- HC.overlap[,list(chr,s1.start,s1.stop)]
colnames(HC.pos) <- c("chr", "start", "stop")
HC.ivl.merged <- intervals_merge(HC.pos)
HC.ivl.merged[, chr := as.numeric(chr)]
HC.ivl.merged[, length := end - start]
setkey(HC.ivl.merged, chr,start,end)

sum(HC.ivl.merged$length) #4kb = 70460000, #5kb=87194000, new 4kb = 70347000, new 4kb =70351000


### positions of overlapping chimp windows
HCnobin.overlap <- YRI.chimp.gor.bin[s1.TajD==1 & s2.TajD==1,]
HCnobin.pos <- HCnobin.overlap[,list(chr,s1.start,s1.stop)]
colnames(HCnobin.pos) <- c("chr", "start", "stop")
HCnobin.ivl.merged <- intervals_merge(HCnobin.pos)
HCnobin.ivl.merged[, chr := as.numeric(chr)]
HCnobin.ivl.merged[, length := end - start]
setkey(HCnobin.ivl.merged, chr,start,end)
sum(HCnobin.ivl.merged$length) #15883000, now 15887000


############### haplotypes ###
### molly's haplotypes
haps <- read.table("/Users/boettger/postdoc/bal_sel/TajD/molly_haplotypes_genes.txt", header=T)
haps.dt <- data.table(haps)[, list(start=min(Position), end=max(Position)),by=c('gene', 'chr')]


foverlaps(molly.haps.dt, HC.ivl.merged, type="any", nomatch=0L)
    chr     start       end length   i.start      stop
 1:   1 237790000 237802000  12000 237795496 237795496
 2:   1 237790000 237802000  12000 237796276 237796276
 3:   2 129386000 129397000  11000 129388993 129388993
 4:   3 116449000 116457000   8000 116454070 116454070
 5:   3 116449000 116457000   8000 116454385 116454385
 6:   3 143683000 143694000  11000 143685046 143685046
 7:   3 143683000 143694000  11000 143688035 143688035
 8:   4  33310000  33340000  30000  33323144  33323144
 9:   4  33310000  33340000  30000  33324082  33324082
10:   4  56134000  56160000  26000  56144663  56144663
11:   4  56134000  56160000  26000  56147967  56147967
12:   4  57915000  57930000  15000  57918296  57918296
13:   4  57915000  57930000  15000  57918427  57918427
14:   4  57915000  57930000  15000  57918492  57918492
15:   4  57915000  57930000  15000  57920048  57920048
16:   4  57915000  57930000  15000  57920087  57920087
17:   4  84565000  84577000  12000  84567638  84567638
18:   4  84565000  84577000  12000  84569070  84569070
19:   4 144649000 144661000  12000 144655406 144655406
20:   4 144649000 144661000  12000 144655489 144655489
21:   4 144649000 144661000  12000 144655823 144655823
22:   4 144649000 144661000  12000 144658471 144658471
23:   4 144649000 144661000  12000 144658748 144658748
24:   4 144649000 144661000  12000 144659625 144659625
25:   4 144662000 144670000   8000 144662054 144662054
26:   5 128312000 128321000   9000 128312506 128312506
27:   5 128312000 128321000   9000 128314616 128314616
28:   7 155045000 155056000  11000 155049958 155049958
29:   7 155045000 155056000  11000 155051811 155051811
30:   8  72531000  72538000   7000  72535024  72535024
31:   8  72531000  72538000   7000  72536493  72536493
32:   9  78336000  78349000  13000  78337373  78337373
33:   9 103852000 103864000  12000 103854887 103854887
34:   9 103852000 103864000  12000 103858769 103858769
35:  10 124420000 124428000   8000 124424638 124424638
36:  10 124420000 124428000   8000 124426504 124426504
37:  11  81476000  81504000  28000  81489841  81489841
38:  11  81476000  81504000  28000  81492293  81492293
39:  14  22318000  22331000  13000  22321419  22321419
40:  14  22318000  22331000  13000  22322973  22322973
41:  16    934000    952000  18000    939565    939565
42:  16    934000    952000  18000    941234    941234
43:  18  58435000  58454000  19000  58437878  58437878
44:  18  58435000  58454000  19000  58438910  58438910
45:  19   8666000   8678000  12000   8672644   8672644
46:  19   8666000   8678000  12000   8673855   8673855
47:  22  47641000  47647000   6000  47642100  47642100
48:  22  47641000  47647000   6000  47645285  47645285

### merge haps.dt with ...?
colnames(YRI.chimp.gor.bin)[9:10] <- c("start", "end")
YRI.chimp.gor.bin <- data.frame(YRI.chimp.gor.bin)[,c(1:12,14:19)]
YRI.chimp.gor.bin.dt <- data.table(YRI.chimp.gor.bin)
## get rid of duplicate

##setkeys
setkey(YRI.chimp.gor.bin.dt, chr, start, end)
setkey(haps.dt, chr, start, end)
haps.apes <- foverlaps(haps.dt, YRI.chimp.gor.bin.dt, type="any", nomatch=0L)[, 
  list(chr, s1.midpoint)]

setkey(haps.apes, chr, s1.midpoint)
setkey(YRI.chimp.gor.bin.dt, chr, s1.midpoint)

#### merge with YRI.chimp.gor.bin.dt
YRI.chimp.gor.bin.dt$transspec.hap <- 0
YRI.chimp.gor.bin.dt[haps.apes, transspec.hap := 1]
## calculate transition probabilities from here

##expect 4 with no bin widen
108138/2162777*108138/2162777*1844/2162777*2162777
[1] 4.609928

## get 43
dim(YRI.chimp.gor.bin.dt[s1.TajD==1 & s2.TajD==1 & transspec.hap==1])
[1] 43... now I find only 36, what happened?
### which haplotypes were they?
setkey(YRI.chimp.gor.bin.dt, chr, start, end)
setkey(haps.dt, chr, start, end)
foverlaps(
  haps.dt, 
  YRI.chimp.gor.bin.dt[s1.TajD==1 & s2.TajD==1 & transspec.hap==1,
    list(chr, start, end)], 
  type="any", nomatch=0L)

##MCMC simulations
#look at 10,000 tests of chimp-human MCMC
hum_chimp_trans_1000 <- as.matrix(read.csv("/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/chimp_human_transspec.txt", header=F))
#make numeric
hum_chimp_trans_1000 <- as.numeric(hum_chimp_trans_1000)
hist(hum_chimp_trans_1000)

z <- (nrow(YRI.chimp.gor.bin.dt[s1.TajD==1 & s2.TajD==1 & transspec.hap==1]) - mean(hum_chimp_trans_1000)) / sd(hum_chimp_trans_1000)
pvalue2sided=2*pnorm(-abs(z))
#1.380646e-21, but this is p-value for all three together not compared to trans-species haps



#look at 10,000 tests of chimp-human-gorlla-mus MCMC
mam_1000 <- as.matrix(read.csv("/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/hum_chp_gor_mus_wi5kb.txt", header=F))
#make numeric
mam_1000 <- as.numeric(mam_1000)
hist(mam_1000)

z <- (nrow(YRI.chimp.gor.bin.dt[s1.TajD==1 & s2.TajD==1 & transspec.hap==1]) - mean(mam_1000)) / sd(mam_1000)
pvalue2sided=2*pnorm(-abs(z))
#1.380646e-21, but this is p-value for all three together not compared to trans-species haps





### what about overlapping windows,not the same window???
dim(YRI.chimp.gor.bin.dt[s1.ext==1 & s2.ext==1 & transspec.hap==1]) #105
setkey(YRI.chimp.gor.bin.dt, chr, start, end)
setkey(haps.dt, chr, start, end)
my4kb.overlap<- foverlaps(
  haps.dt, 
  YRI.chimp.gor.bin.dt[s1.ext==1 & s2.ext==1 & transspec.hap==1,
    list(chr, start, end)], 
  type="any", nomatch=0L)



##expect 48 with 4kb bin widen --> not as good
350574/2162777*351240/2162777*1844/2162777*2162777


## for whole genome
trans.prob.mat(YRI.chimp.gor.bin.dt$s1.TajD)
trans.prob.mat(YRI.chimp.gor.bin.dt$s2.TajD)

trans.prob.mat(as.numeric(YRI.chimp.gor.bin.dt[, s2.TajD==1 & s1.TajD==1]))

trans.prob.mat(YRI.chimp.gor.bin.dt$transspec.hap)


### no HLA
#MHC Location: chr6:28,477,797-33,448,354
#+/-2kb: 28477795-33448356

## expect for no HLA (4.57)
107336/2158982*107763/2158982*1844/2158982*2158982

YRI.chimp.gor.bin.noHLA <- YRI.chimp.gor.bin.dt[!(chr==6 & s1.midpoint>=28477797 & s1.midpoint<=33448356)]
dim(YRI.chimp.gor.bin.noHLA[s1.TajD==1 & s2.TajD==1 & transspec.hap==1])
### which haplotypes were they?
setkey(YRI.chimp.gor.bin.noHLA, chr, start, end)
setkey(haps.dt, chr, start, end)
my.overlap<- foverlaps(
  haps.dt, 
  YRI.chimp.gor.bin.noHLA[s1.TajD==1 & s2.TajD==1 & transspec.hap==1,
    list(chr, start, end)], 
  type="any", nomatch=0L)
## 14 haplotypes/genes of 123 have overlap with old human/chimp haplotypes

###### what genes overlapped all three species???
HCnobinnoHLA.overlap <- YRI.chimp.gor.bin.noHLA[s1.TajD==1 & s2.TajD==1 &s3.TajD==1,]
HCnobinnoHLA.pos <- HCnobinnoHLA.overlap[,list(chr,start,end)]
colnames(HCnobinnoHLA.pos) <- c("chr", "start", "stop")
HCnobinnoHLA.ivl.merged <- intervals_merge(HCnobinnoHLA.pos)

ape.overlap.regions <- paste(HCnobinnoHLA.ivl.merged$chr, HCnobinnoHLA.ivl.merged$start, HCnobinnoHLA.ivl.merged$end, sep = ":" )

ape.overlap.genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"), filters = "chromosomal_region", values = ape.overlap.regions, mart = grch37)
ape.overlap.pcgenes <- ape.overlap.genes[which(ape.overlap.genes[,3]=="protein_coding"),]

write.table(ape.overlap.pcgenes[,1], "/Users/boettger/postdoc/bal_sel/gene_lists/ape.overlap.pcgenes_aug8.txt", row.names=F, col.names=F, quote=F)



trans.prob.mat(YRI.chimp.gor.bin.noHLA$s1.TajD)
trans.prob.mat(YRI.chimp.gor.bin.noHLA$s2.TajD)

trans.prob.mat(as.numeric(YRI.chimp.gor.bin.noHLA[, s2.TajD==1 & s1.TajD==1]))

trans.prob.mat(YRI.chimp.gor.bin.noHLA$transspec.hap)

#### get MCMC output 
trans.mcmc <- read.table("/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/comb_chimphuman_transspec_runs.txt")

z <- (nrow(YRI.chimp.gor.bin.noHLA[s1.TajD==1 & s2.TajD==1 & transspec.hap==1]) - mean(trans.mcmc$V1)) / sd(trans.mcmc$V1)
pvalue2sided=2*pnorm(-abs(z))




foverlaps(coding.dt, HCnobin.ivl.merged, type="any", nomatch=0L)
   chr     start       end length   i.start     i.end
1:   1 202089000 202093000   4000 202092580 202092580
2:  12 108617000 108621000   4000 108618616 108618616
3:  16  81279000  81283000   4000  81279120  81279120

phyper(3, 15883000, 2407055000-15883000, 320, lower.tail =  FALSE) #p=.16


foverlaps(coding.dt, HC.ivl.merged, type="any", nomatch=0L)
    chr     start       end length   i.start     i.end
 1:   1 201456000 201471000  15000 201459399 201459399
 2:   1 202085000 202097000  12000 202092580 202092580
 3:   1 223814000 223818000   4000 223816383 223816383
 4:   2 228542000 228562000  20000 228560728 228560728
 5:   7  95493000  95503000  10000  95499234  95499234
 6:   8 146024000 146032000   8000 146028316 146028316
 7:  12 102050000 102061000  11000 102055018 102055018
 8:  12 108605000 108625000  20000 108618616 108618616
 9:  12 124345000 124353000   8000 124352590 124352590
10:  15  22938000  22948000  10000  22939192  22939192
11:  16  81270000  81283000  13000  81279120  81279120
12:  18  56245000  56249000   4000  56247744  56247744

## what did I expect? 
70347000/2407055000*320 = 9
phyper(12, 70347000, 2407055000-70347000, 320, lower.tail =  FALSE) p=0.14







############################# FUNCTIONS ##############################

cutoffs.2 <- function(overlap.data, s1.cutoff, s4.cutoff){
  overlap.data[s1.TajD <= s1.cutoff, s1.TajD:= 0]
  overlap.data[s4.TajD <= s4.cutoff, s4.TajD:= 0]
  overlap.data[s1.TajD > s1.cutoff, s1.TajD:= 1]
  overlap.data[s4.TajD > s4.cutoff, s4.TajD:= 1]
  return(overlap.data)
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

cutoffs.4 <- function(overlap.data, s1.cutoff, s2.cutoff, s3.cutoff, s4.cutoff){
  overlap.data[s1.TajD < s1.cutoff, s1.TajD:= 0]
  overlap.data[s2.TajD < s2.cutoff, s2.TajD:= 0]
  overlap.data[s3.TajD < s3.cutoff, s3.TajD:= 0]
  overlap.data[s4.TajD < s4.cutoff, s4.TajD:= 0]
  overlap.data[s1.TajD >= s1.cutoff, s1.TajD:= 1]
  overlap.data[s2.TajD >= s2.cutoff, s2.TajD:= 1]
  overlap.data[s3.TajD >= s3.cutoff, s3.TajD:= 1]
  overlap.data[s4.TajD >= s4.cutoff, s4.TajD:= 1]
  
  return(overlap.data)
}



bin.widen.3 <- function(overlap.data, offset=2000) {
  #s1
  if (any(overlap.data$s1.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s1.TajD == T, 
          list(s1.midpoint-offset, s1.midpoint+offset)]))
    s1_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s1.midpoint, ivls), length)])
  } else {
    s1_ext <- rep(0,nrow(overlap.data))
  }
  
  #s2
  if (any(overlap.data$s2.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s2.TajD == T, 
          list(s2.midpoint-offset, s2.midpoint+offset)]))
    s2_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s2.midpoint, ivls), length)])
  } else {
    s2_ext <- rep(0,nrow(overlap.data))
  }
  
  #s3
  if (any(overlap.data$s3.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s3.TajD == T, 
          list(s3.midpoint-offset, s3.midpoint+offset)]))
    s3_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s3.midpoint, ivls), length)])
  } else {
    s3_ext <- rep(0,nrow(overlap.data))
  }
  
  return(data.table(s1_ext, s2_ext, s3_ext))
}

bin.widen.2 <- function(overlap.data, offset=4000) {
  #s1
  if (any(overlap.data$s1.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s1.TajD == T, 
          list(s1.midpoint-offset, s1.midpoint+offset)]))
    s1_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s1.midpoint, ivls), length)])
  } else {
    s1_ext <- rep(0,nrow(overlap.data))
  }
  
  #s2
  if (any(overlap.data$s2.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s2.TajD == T, 
          list(s2.midpoint-offset, s2.midpoint+offset)]))
    s2_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s2.midpoint, ivls), length)])
  } else {
    s2_ext <- rep(0,nrow(overlap.data))
  }
  
  return(data.table(s1_ext, s2_ext))
}


bin.widen.4 <- function(overlap.data, offset=5000) {
  #s1
  if (any(overlap.data$s1.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s1.TajD == T, 
          list(s1.midpoint-offset, s1.midpoint+offset)]))
    s1_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s1.midpoint, ivls), length)])
  } else {
    s1_ext <- rep(0,nrow(overlap.data))
  }
  
  #s2
  if (any(overlap.data$s2.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s2.TajD == T, 
          list(s2.midpoint-offset, s2.midpoint+offset)]))
    s2_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s2.midpoint, ivls), length)])
  } else {
    s2_ext <- rep(0,nrow(overlap.data))
  }
  
  #s3
  if (any(overlap.data$s3.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s3.TajD == T, 
          list(s3.midpoint-offset, s3.midpoint+offset)]))
    s3_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s3.midpoint, ivls), length)])
  } else {
    s3_ext <- rep(0,nrow(overlap.data))
  }
  
  #s4
  if (any(overlap.data$s4.TajD == T)) {
    ivls = interval_union(
      Intervals(
        overlap.data[s4.TajD == T, 
          list(s4.midpoint-offset, s4.midpoint+offset)]))
    s4_ext <- as.numeric(overlap.data[, lapply(interval_overlap(s4.midpoint, ivls), length)])
  } else {
    s4_ext <- rep(0,nrow(overlap.data))
  }
 
  return(data.table(s1_ext, s2_ext, s3_ext, s4_ext))
}




get.overlap.data.addon <- function(species1.bed, species2.bed, MAXdist = 2000){

  ### get midpoints of windows
  #species1.bed$s1.midpoint <- round((species1.bed$stop + species1.bed$start)/2)
  #colnames(species1.bed)[c(5,6,7,9)] <- c("s1.start", "s1.stop", "s1.TajD", "midpoint")

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
  ## chr 1 isn't present???
  ## s1 midpoint not present?
  #s1.s2[chr==1,]
  
  #recalculate actual midpoints (ereased/unclear because merged on non-exact value), looks like s1.midpoint is maintained
  s1.s2 <- cbind(s1.s2, s1.s2[, .(s2.midpoint = round((s2.start + s2.stop)/2), s1.midpoint = (s1.start + s1.stop)/2)])
  
  #find distance between midpoints and filter
  s1.s2[, dist := abs(s1.midpoint-s2.midpoint)]
  s1.s2.keep <- s1.s2[dist<=MAXdist]
  
  ### keep closer of windows from s2, so no duplcate windows 
  s1.s2.return <- s1.s2.keep[s1.s2.keep[, .I[which.min(dist)], by=c("chr","s2.midpoint")]$V1]
  
  
  return(s1.s2.return)
}


get.overlap.data <- function(species1.bed, species2.bed, MAXdist = 5000){

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


intervals_merge <- function(ivl.dt) {
  # takes a data table of intervals with 3 columns (chr, start, and end) and
  # merges intervals that overlap.
  require(intervals)
  required_col_names <- c('chr','start','end')

  if (any(names(ivl.dt) != required_col_names)) {
    names(ivl.dt) <- c('chr','start','end')
  }
  
  
  single_chr <- function(start, end) {
    merged_ivls <- as.matrix(interval_union(Intervals(cbind(start, end))))
    return(data.table(start=merged_ivls[,1], end=merged_ivls[,2]))
  }
  as.data.table(ivl.dt)[, single_chr(start, end), by=chr]
}


trans.prob.mat <- function(table.bin.species){
  from1 <-table(table.bin.species[2:length(table.bin.species)] - table.bin.species[1:length(table.bin.species)-1])
  from0 <- table(table.bin.species[2:length(table.bin.species)] + table.bin.species[1:length(table.bin.species)-1])
  output <- matrix(nrow=1,ncol=4)
  colnames(output) <- c("0:0", "0:1", "1:0", "1:1")
  output[1,1] <- from0[1]
  output[1,2] <- from1[1]
  output[1,3] <- from1[3]  # I think this should be [3]
  output[1,4] <- from0[3]
  return(output)
}
 