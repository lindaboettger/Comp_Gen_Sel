#!/bin/bash

## Use tools
use Samtools
#use Python-2.6
use GCC-5.2
use Python-2.7 
use GridEngine8
use R-3.3

## set variables
#win=4000
#step=1000
#nchr=20

# Get species-specific variables #
source $species_param

# HWEpval=0.00001
# SBpval=0.000001 # this is the smallest number the snp stats keeps
# parentpath='/seq/vgb/linda/bal_sel/from_oracle'
# bampath=$parentpath'/YRI_bams/YRI_bams_chr'$chr
# outpath=$parentpath'/YRI_out_beta/YRI_out_chr'$chr

#### make folders for this chromosome ###
#mkdir -p $bampath
#mkdir -p $wd  ## this isnt working... can't set working before it's made
##make directories for liftovers
mkdir -p $outpath
mkdir -p $liftpath
mkdir -p $mafpath
mkdir -p $liftout


#set nchr based on number of individuals in analysis
# npeeps=$(wc -l $bams$chr.bamfilelist | grep -Po '^\d+')
# let "nchr=$npeeps*2"




if [ $new_species == T ] 
  then
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
        #    -nThreads 15 \
        #-doGlf 1 log scaled likelihood ratios to the most likely
        # -doGlf 2 Beagle 
    #doGeno 2 is to find sites with both homozygous classes (requires doPost)

    /seq/vgb/linda/bin/Feb2017_angsd/angsd/angsd \
        -bam $bamfilelist \
        -out $maffile$chr \
        -ref $ref \
        -GL 2 \
        -doGlf 1 \
        -r $chr \
        -snp_pval 1e-3 \
        -doCounts 1 \
        -minInd $minInd \
        -doGeno 2 \
        -doPost 1 \
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


    # ##### ################# filter HWE and SB ############################
    # ## unzip
    gunzip $maffile$chr'.snpStat.gz'
    gunzip $maffile$chr'.geno.gz'

    
    # ## get rid of columns after random newline (also gets rid of header)
    grep -e "^$chr" $maffile$chr'.snpStat' | awk -F":" '$1=$1' OFS="\t" > $mafpath/$chr'_sep.tmp'

    # ### filter on HWE
    awk -v HWE=$HWEpval '{ if ($11 < HWE) print $0 }' $mafpath/$chr'_sep.tmp' \
         > $mafpath/$chr'_badHWE.txt'


    ###filter on Strand bias
    awk -v SB=$SB 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($7) > SB) print $0}' \
        $mafpath/$chr'_sep.tmp' \
        > $mafpath/$chr'_badSB.txt'

    # filter for only keep snps with both homozygous classes ###
    # The number of columns will vary and this seemed hard 
    # to work with in awk

    # also get rid of -NaN and NaN (acutally the above SB filter captures nan)
    #awk '$7=="-nan" || $7=="nan"' $mafpath/$chr'_sep.tmp' >> $mafpath/$chr'_badSB.txt'


    Rscript /seq/vgb/linda/bal_sel/scripts_2017/both_homoz.R \
        $maffile$chr'.geno' \

    ## combine
    cat $mafpath/$chr'_badHWE.txt' $mafpath/$chr'_badSB.txt' > $mafpath/$chr'_badSNPs.txt'

    ### Get rid of positions without both homozygotes
    ### if we are filtering for that
    if [$filt_both_homo == T]; then
        echo "getting rid of snps without both homozygous classes"
        cp $mafpath/$chr'_badSNPs.txt' $mafpath/$chr'_badSNPs_tmp.txt'
        cat $mafpath/$chr'_badSNPs_tmp.txt' $maffile$chr'.geno_lacking_homoz.txt' \
            > $mafpath/$chr'_badSNPs.txt'
        rm $mafpath/$chr'_badSNPs_tmp.txt'
    fi

    # ## check that two are duplicates between badSB and badSNP
    #sort $outpath'/badSNPs.txt' | uniq -c | grep -v '^ *1 '
    nlines=$(wc -l < $mafpath/$chr'_badSNPs.txt')


    if [ $nlines > 0 ]; then
        echo "there are bad snps"
        # ### keep MAFs for SNPs that aren't in badSNPs
        awk -F"\t" 'FNR==NR{a[$1,$2];next} !(($1,$2) in a)' $mafpath/$chr'_badSNPs.txt' <(gzip -dc $maffile$chr'.mafs.gz') \
        > $maffile$chr'_filtered.mafs'
    else
        echo "there are no bad snps"
        gunzip $maffile$chr'.mafs.gz' > $maffile$chr'_filtered.mafs'
        cp $maffile$chr'.mafs' $maffile$chr'_filtered.mafs'
        gzip $maffile$chr'.mafs'
    fi
    
    
    # gzip $outpath'/YRIfiltered_'$chr'.thetas'

    rm $mafpath/$chr'_sep.tmp'
    gzip $maffile$chr'.snpStat'
    gzip $maffile$chr'.geno'


    
fi







################ MAF to beta input ################
maf=$maffile$chr'_filtered.mafs'

# ## they aren't whole numbers because genotype probs used to calculate freq
# put win in this file so it will be in the outfile
awk '{if (NR!=1) print $2, $6*$8*2, $8*2}' $maf > $wd/$species'_BetaIn_'$chr'_win'$win'.txt'

gzip $maffile$chr'_filtered.mafs'



# # CALCULATE BETA ####

# # ### Parameters 
# # * -i: Path of input file
# # * -w: Total window size (default: 1000, corresponding to 500 bp on either side of the center (i.e. core) SNP)
# # * -p: Value of p (default: 20)
# # * -m: Minimum folded frequency of core SNP, exclusive, can range from 0 to .5 (default: 0)
# * -fold: Use folded version (default: false)
## run betascan un with a 2000 basepair window, a p parameter value of 50 and using the folded version of Beta:
# # -p: Value of p (default: 20)
# betascan=/seq/vgb/linda/bin/BetaScan-master/BetaScan.py
# python $betascan -i $wd/$species'_BetaIn_'$chr'_win'$win'.txt' -w $win -fold 


# #convert betas into window output for liftover #
# # $wd/'Betas_'$species'_BetaIn_chr'$chr'.txt'


# cat $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'.txt' | python /seq/vgb/linda/bal_sel/scripts_2017/format_beta_for_liftover.py \
#     > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_toliftover.bed'


## call script to submit betascan by chr and format it for liftover, do the liftover for this chr
betascan=/seq/vgb/linda/bin/BetaScan-master/BetaScan.py


# -p: Value of p (default: 20)
#run betascan

echo "I am about to start betascan"
python $betascan -i $wd/$species'_BetaIn_'$chr'_win'$win'.txt' -w $win -fold 

## Smooth windows using smwin and smstep --FIX THIS SO THAT IT ALSO PRINTS OUT THE BP POSITION
python /seq/vgb/linda/bal_sel/scripts_2017/smooth_win.py $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'.txt' $smwin $smstep \
 > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'sm.txt'

#convert smoothed betas into window output for liftover #
cat $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'sm.txt' | python /seq/vgb/linda/bal_sel/scripts_2017/format_smooth_beta_for_liftover.py \
    > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_smoothed_toliftover.bed'

#also convert bases to bed for lift
$wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'.txt'

awk -v chr=$chr '{print chr, $1, $1+1, chr"_"$1, $2}' $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'.txt' \
   > $wd/$species'_BetaIn_'$chr'_win'$win'_bp_toliftover.bed'


#convert betas into window output for liftover #
# cat $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'sm.txt' | python /seq/vgb/linda/bal_sel/scripts_2017/format_beta_for_liftover.py \
#     > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_toliftover.bed'



### do the liftover for this chr
#do the liftover if this isn't human
echo $species

if ! [ $species == 'YRI' ]; then
  /seq/vgb/linda/bin/liftOver \
    -minMatch=$minmatch \
    $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_smoothed_toliftover.bed' \
    $this_chain \
    $liftout/$species'_to_hg38_'${chr}'_win'$win'_smoothed.bed' \
    $liftpath/$species'_to_hg38_'${chr}'_win'$win'_smoothed_Unmapped.bed'

    echo "I am not a human!!!"

fi




# #do the liftover if this isn't human
# ### NOTE GORILLA WILL REQUIRE 2 LIFTOVERS 
# echo $species

# if ! [ $species == 'YRI' ] && ! [ $species == 'gorilla' ]; then
#   /seq/vgb/linda/bin/liftOver \
#     -minMatch=$minmatch \
#     $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_toliftover.bed' \
#     $this_chain \
#     $liftout/$species'_to_hg19_'${chr}'_win'$win'_.bed' \
#     $liftpath/$species'_to_hg19_'${chr}_Unmapped.bed

#     echo "I am not a gorilla!!!"

#     else
#     if [ $species == 'gorilla' ]
#     then

#     echo "I am a fucking gorilla"
#         #### gorgor4 to gorgor3to human ####
#       /seq/vgb/linda/bin/liftOver \
#       -minMatch=0.5 \
#       $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_toliftover.bed' \
#       $chainpath/gorGor4ToGorGor3.over.chain \
#       $liftpath/gorgor4_to_gorgor3_${chr}.bed \
#       $liftpath/gorgor4_to_gorgor3_${chr}_Unmapped.bed


#       #### gorilla to human ####
#       /seq/vgb/linda/bin/liftOver \
#       -minMatch=0.1 \
#       $liftpath/gorgor4_to_gorgor3_${chr}.bed \
#       $chainpath/gorGor3ToHg19.over.chain \
#       $liftout/gorgor3_to_hg19_${chr}'_win'$win'_.bed' \
#       $liftpath/gorgor3_to_hg19_${chr}_Unmapped.bed

#     fi
# fi


############ this needs to go one level up !!!!




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





