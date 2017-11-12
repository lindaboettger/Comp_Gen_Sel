#!/bin/bash

## Use the samtools package
use Samtools
use Python-2.6
use GCC-5.2

## set variables
#win=4000
#step=1000
#nchr=20

HWEpval=0.00001
SBpval=0.000001 # this is the smallest number the snp stats keeps
parentpath='/seq/vgb/linda/bal_sel/from_oracle'
bampath=$parentpath'/YRI_bams/YRI_bams_chr'$chr
outpath=$parentpath'/YRI_out_beta/YRI_out_chr'$chr

#### make folders for this chromosome ###
#mkdir -p $bampath
mkdir -p $outpath

#set nchr based on number of individuals in analysis
npeeps=$(wc -l $bampath/YRIrentschr$chr.bamfilelist | grep -Po '^\d+')
let "nchr=$npeeps*2"


# Set human-sepcific variables 
MinDepth=70
MaxDepth=700
minInd=22
minMapQ=30
minQ=20


################## Calculate MAF.##################

## for genotype likelihoods
	# -GL=0: 
	# 1: SAMtools
	# 2: GATK ----- we are using GATK genotype likelihoods
	# 3: SOAPsnp
	# 4: SYK

#snp_pval 1e-2 limits to likely only truly variable sites

#0.917-35-g89d1f50 (htslib: 1.3.2-208-gd8d0323-dirty) build(Feb 13 2017 12:17:15)
#https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
############### get genotype prob MAF, filter on quality etc.  screen -dr 293
## here I fixed the dir name and didn't put in minHWEpval
## domaf is 4 -- caluclated from genotype probabilities
# do counts requried for min and max depth filters
#GL 2 GATK genotype likelihoods
#SNP_pval - only likely real snps
#not using baq because we already did local realignment around indels
#screen -dr 357
# use GCC-5.1 screen -dr 27567 YRIoutFold_vsmalltest2_$chr
#angsd version: 0.917-35-g89d1f50 (htslib: 1.3.2-208-gd8d0323-dirty) build(Feb 13 2017 12:17:15). <- maybe this is the problem??!!
/seq/vgb/linda/bin/Feb2017_angsd/angsd/angsd \
    -nThreads 15 \
    -bam $bampath/YRIrentschr$chr.bamfilelist \
    -out $outpath/YRIoutFold_$chr \
    -ref /seq/vgb/linda/references/hg19/Homo_sapiens_assembly19.fasta \
    -GL 2 \
    -snp_pval 1e-2 \
    -doCounts 1 \
    -minInd $minInd \
    -domaf 1 \
    -doMajorMinor 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -uniqueOnly 1 \
    -setMinDepth $MinDepth \
    -setMaxDepth $MaxDepth \
    -minMapQ $minMapQ \
    -minQ $minQ \
    -doSnpStat 1 \
    -fold 1


#It didnt' like my index file, so I re-made the index in my own folder
 #cp /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta /seq/vgb/linda/references/hg19/Homo_sapiens_assembly19.fasta
 #samtools faidx /seq/vgb/linda/references/hg19/Homo_sapiens_assembly19.fasta

#-ref /seq/vgb/references/human_hg19/hg19.fa

# ##### ################# filter HWE and SB ############################
# ## unzip
gunzip $outpath'/YRIoutFold_'$chr'.snpStat.gz'

# ## get rid of columns after random newline (also gets rid of header)
grep -e "^$chr" $outpath'/YRIoutFold_'$chr'.snpStat' | awk -F":" '$1=$1' OFS="\t" > $outpath'/sep.tmp'

# ### filter on HWE
awk -v HWE=$HWEpval '$11 < HWE' $outpath'/sep.tmp' > $outpath'/badHWE.txt'

# ### filter on SB (I seem to be using SB3)
awk -v SB=$SBpval '$9 < SB' $outpath'/sep.tmp' > $outpath'/badSB.txt'

cat $outpath'/badHWE.txt' $outpath'/badSB.txt' > $outpath'/badSNPs.txt'

# ## check that two are duplicates between badSB and badSNP
#sort $outpath'/badSNPs.txt' | uniq -c | grep -v '^ *1 '

# ### keep MAFs for SNPs that aren't in badSNPs
awk -F"\t" 'FNR==NR{a[$1,$2];next} !(($1,$2) in a)' $outpath'/badSNPs.txt' <(gzip -dc $outpath/'YRIoutFold'_$chr'.mafs') > $outpath'/YRIfiltered_'$chr'.mafs'
# gzip $outpath'/YRIfiltered_'$chr'.thetas'

rm $outpath/sep.tmp
gzip $outpath'/YRIoutFold_'$chr'.snpStat'



# ################ MAF to beta input ################
maf=$outpath'/YRIfiltered_'$chr'.mafs'

# ## they aren't whole numbers because genotype probs used to calculate freq
awk '{if (NR!=1) print $2, $6*$8*2, $8*2}' $maf > $parentpath/'YRI_out_beta/YRI_BetaIn_chr'$chr'.txt'


# CALCULATE BETA ####

# # ### Parameters 
# # * -i: Path of input file
# # * -w: Total window size (default: 1000, corresponding to 500 bp on either side of the center (i.e. core) SNP)
# # * -p: Value of p (default: 20)
# # * -m: Minimum folded frequency of core SNP, exclusive, can range from 0 to .5 (default: 0)
# * -fold: Use folded version (default: false)
## run betascan un with a 2000 basepair window, a p parameter value of 50 and using the folded version of Beta:
# -p: Value of p (default: 20)
betascan=/seq/vgb/linda/bin/BetaScan-master/BetaScan.py
python $betascan -i $parentpath/'YRI_out_beta/YRI_BetaIn_chr'$chr'.txt' -w 1000 -fold 

# # get top 1% of values on this --- do this when all chr are done!!!!
# my_length=$(wc -l < $parentpath/'YRI_out_beta/YRI_BetaIn_chr'$chr'.txt')
# last_val=$((my_length / 100))
# sort --key 2 --numeric-sort $parentpath/'YRI_out_beta/Betas_YRI_BetaIn_chr'$chr'.txt' | tail -$last_val


# # The first column contains the coordinate of each variant, and the second c
# # ontains the frequency of the derived allele, in number of haploid individuals, of the variant. The third column contains the sample si
# # ze, in number of haploid individuals, that were used to calculate the frequency of that variant. The file should be sorted by position
# #  (the unix command sort -g will do this for you). Variants with frequencies of exactly 0% or 100% should not be included. In practice,
# #  for folded Beta, it doesn't matter if the derived, ancestral, or already folded allele frequency is used in the second column, as Bet
# # aScan will fold the frequency anyway. The scan should be run on each chromosome separately. An example of a sample file is below:

# head -10 $parentpath/'YRI_out_beta/YRI_BetaIn_chr'$chr'.txt' > $parentpath/'YRI_out_beta/YRI_BetaIn_test.txt'
# python $betascan -i $parentpath/'YRI_out_beta/YRI_BetaIn_test.txt' -w 2000 -fold



# #b is, in effect, a weighted average of SNP counts based on their frequency similarity to the core SNP. 
# ## saple output
# # 16050252        0.392194257909
# # 16050408        3.25659986413
# # 16050612        1.67836502531
# # 16050678        2.88298962367
# # 16050822        0.320272302691
# # 16050933        0.507666396863
# # 16050967        0.599795052797
# # 16051107        3.24838710462
# # 16051241        2.7521537148
# # 16051255        3.08845756698

# # ```
# # 14  2 99  
# # 25  1 100  
# # 47  99  100
# # 48  82  95
# # 103 10  100
# # 245 93  96
# # ```





