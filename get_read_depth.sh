#!/bin/bash

## Use tools
use Samtools
#use Python-2.6
use GCC-5.2
use Python-2.7 
use GridEngine8

## set variables
#win=4000
#step=1000
#nchr=20

# Get species-specific variables #
#species_param=/seq/vgb/linda/bal_sel/scripts_2017/YRI_param_hg38.sh
source $species_param
mkdir -p $bampath/readDepth


#-maxDepth 250 \ 
# Defaults to 100. Sites with more than maxDepth reads are counted as having 100 reads.

#set a reasonable max depth
numindivs=$(wc -l < $bamfilelist)
myMaxDepth=$(($numindivs*50))

/seq/vgb/linda/bin/Feb2017_angsd/angsd/angsd \
    -r $chr \
    -bam $bamfilelist \
    -out $bampath/readDepth/$species'_'$chr'_readDepth' \
    -ref $ref \
    -doCounts 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -uniqueOnly 1 \
    -minMapQ $minMapQ \
    -minQ $minQ \
    -maxDepth $myMaxDepth \
    -doDepth 1

# awk 'BEGIN { OFS = "\n" } { $1=$1; print }' $bampath/$species'_readDepth.depthGlobal' | sort -g \
# 	> $bampath/$species'_readDepth.depthGlobal.ordered'


# #get rid of positions fully removed by other filters
# totpos=$(grep -v "0" $bampath/$species'_readDepth.depthGlobal.ordered' | wc -l)
# grep -v "0" $bampath/$species'_readDepth.depthGlobal.ordered' \
# 	> $bampath/$species'_readDepth.depthGlobal.ordered.no0'
# linemin=$(echo "$totpos*0.05" | bc)
# lineminint=${linemin%.*}
# linemax=$(echo "$totpos*0.95" | bc)
# linemaxint=${linemax%.*}
# linemean=$(echo "$totpos/2" | bc)
# linemeanint=${linemean%.*}

# MinDepth=$(sed "${lineminint}q;d" $bampath/$species'_readDepth.depthGlobal.ordered.no0')
# MaxDepth=$(sed "${linemaxint}q;d" $bampath/$species'_readDepth.depthGlobal.ordered.no0')
# myMean=$(sed "${linemeanint}q;d" $bampath/$species'_readDepth.depthGlobal.ordered.no0')

# #keep track of cutoffs
# echo "mean: $myMean" > $bampath/$species'_readDepth.outstats'
# echo "max cut: $MaxDepth" >> $bampath/$species'_readDepth.outstats'
# echo "min cut: $MinDepth" >> $bampath/$species'_readDepth.outstats'

# rm $bampath/$species'_readDepth.depthGlobal.ordered.no0'
# rm $bampath/$species'_readDepth.depthGlobal.ordered'


#2687
#4251