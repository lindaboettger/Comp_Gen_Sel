######## takes input file ".geno" from angsd 
## and finds positions that lack both homozygous classes


### 0. get arguments #######################################

args = commandArgs(trailingOnly=TRUE)

### 1. Load packages and functions  ########################

require(data.table)
`%ni%` <- Negate(`%in%`)

#geno <- data.table(read.table('/Users/boettger/Desktop/test.txt', header=F))
geno <- data.table(read.table(args[1],header=F))

#get rid of chr and make numeric
ngeno <- geno[,-1]
ngeno[, names(ngeno) := lapply(.SD, as.numeric)]
has0 <- ngeno[apply(ngeno, 1, function(r) any(r == 0))]
has2 <- ngeno[apply(ngeno, 1, function(r) any(r == 2))]

both.homoz.pos <- merge(has0, has2, by='V2')$V2

#output list of snps to remove because a homozygous class is missing
out.geno <- geno[V2 %ni% both.homoz.pos]
have.both <- geno[V2 %in% both.homoz.pos]
write.table(out.geno, paste(args[1], "lacking_homoz.txt", sep="_"), row.names=F, col.names=F, quote=F)
write.table(have.both, paste(args[1], "both_homoz.txt", sep="_"), row.names=F, col.names=F, quote=F)



