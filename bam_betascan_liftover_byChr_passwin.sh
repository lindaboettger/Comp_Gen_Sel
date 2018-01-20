#!/bin/bash

## Use tools
use Samtools
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



# if this is the first time we're running this species
# run bams through angsd commands to get maf etc. 
if [ $new_species == T ] 
  then
    ################## Calculate MAF.##################

    ## for genotype likelihoods
      # -GL=0: 
      # 1: SAMtools
      # 2: GATK ----- we are using GATK genotype likelihoods
      # 3: SOAPsnp
      # 4: SYK
    ############### get genotype prob MAF, filter on quality etc.  screen -dr 293
    ## here I fixed the dir name and didn't put in minHWEpval
    ## domaf is 4 -- caluclated from genotype probabilities
    # do counts requried for min and max depth filters
    #GL 2 GATK genotype likelihoods
    #SNP_pval - only likely real snps
    #not using baq because we already did local realignment around indels
            #    -nThreads 15 \
        #-doGlf 1 log scaled likelihood ratios to the most likely
        # -doGlf 2 Beagle 
    #doGeno 2 is to find sites with both homozygous classes (requires doPost)

    ####### new version gets windows with data to call 0 beta
    /seq/vgb/linda/bin/Dec2017git_angsd/angsd/angsd \
        -r $chr \
        -bam $bamfilelist \
        -out $maffile$chr \
        -ref $ref \
        -GL 2 \
        -doGlf 1 \
        -snp_pval 1e-6 \
        -doCounts 1 \
        -minInd $minInd \
        -doGeno 2 \
        -doPost 1 \
        -domaf 1 \
        -doMajorMinor 1 \
        -doHWE 1 \
        -remove_bads 1 \
        -only_proper_pairs 1 \
        -uniqueOnly 1 \
        -setMinDepth $MinDepth \
        -setMaxDepth $MaxDepth \
        -minMapQ $minMapQ \
        -minQ $minQ \
        -doSnpStat 1 \
        -doSaf 1 \
        -anc $ref \
        -fold 1

fi



###################################### FILTERING ############################################
#unzip
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

#keep track of how many snps total
wc -l < $maffile$chr'.snpStat' > $maffile$chr'.Nsnps'

# filter for only keep snps with both homozygous classes ###
# The number of columns will vary and this seemed hard 
# to work with in awk

# also get rid of -NaN and NaN (acutally the above SB filter captures nan)
#awk '$7=="-nan" || $7=="nan"' $mafpath/$chr'_sep.tmp' >> $mafpath/$chr'_badSB.txt'


Rscript /seq/vgb/linda/bal_sel/scripts_2017/both_homoz.R \
    $maffile$chr'.geno' \

#replace  with tab (it skips these otherwise) and get first two cols
perl -p -i -e 's/ /\t/g' $maffile$chr'.geno_lacking_homoz.txt'
awk '{print $1,$2}' $maffile$chr'.geno_lacking_homoz.txt' \
    > $maffile$chr'.geno_lacking_homoz_tmp.txt'
mv $maffile$chr'.geno_lacking_homoz_tmp.txt' $maffile$chr'.geno_lacking_homoz.txt' 


## combine HWE and SB filts
cat $mafpath/$chr'_badHWE.txt' $mafpath/$chr'_badSB.txt' | awk '{print $1,$2}' \
 > $mafpath/$chr'_badSNPs.txt'

### which snps were only filtered out for 
### lacking homoz but passed HWE and SB??
### just keep position and give score of 0
awk -F, 'FNR==NR {a[$1]=$0; next}; !($1 in a) {print $1,$2}' \
    $mafpath/$chr'_badSNPs.txt' \
    $maffile$chr'.geno_lacking_homoz.txt' \
    > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_nohomo_passFILT.txt'

awk '{print $2,"0"}' $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_nohomo_passFILT.txt' \
    > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_nohomo_passFILT.tmp'
cp $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_nohomo_passFILT.tmp' \
    $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_nohomo_passFILT.txt'
rm $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_nohomo_passFILT.tmp'



### Get rid of positions without both homozygotes
### if we are filtering for that
if [ $filt_both_homo == T ]; then
    echo "getting rid of snps without both homozygous classes"
    cp $mafpath/$chr'_badSNPs.txt' $mafpath/$chr'_badSNPs_tmp.txt'
    cat $mafpath/$chr'_badSNPs_tmp.txt' $maffile$chr'.geno_lacking_homoz.txt' \
        > $mafpath/$chr'_badSNPs.txt'
    rm $mafpath/$chr'_badSNPs_tmp.txt'
fi

# get unique values
sort $mafpath/$chr'_badSNPs.txt' | uniq > $mafpath/$chr'_badSNPs_tmp.txt'
cp $mafpath/$chr'_badSNPs_tmp.txt' $mafpath/$chr'_badSNPs.txt'
rm $mafpath/$chr'_badSNPs_tmp.txt'
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

#431023 /seq/vgb/linda/bal_sel/population_data/mafs/mus/mus_chr2_filtered.mafs

# gzip $outpath'/YRIfiltered_'$chr'.thetas'

rm $mafpath/$chr'_sep.tmp'
gzip $maffile$chr'.snpStat'
gzip $maffile$chr'.geno'


    


################################ MAF to beta input ################################
maf=$maffile$chr'_filtered.mafs'

# ## they aren't whole numbers because genotype probs used to calculate freq
# put win in this file so it will be in the outfile
awk '{if (NR!=1) print $2, $7*$9*2, $9*2}' $maf > $wd/$species'_BetaIn_'$chr'_win'$win'.txt'

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


## call script to submit betascan by chr and format it for liftover, do the liftover for this chr
betascan=/seq/vgb/linda/bin/BetaScan-master/BetaScan.py


# -p: Value of p (default: 20)
#run betascan

echo "I am about to start betascan"
python $betascan -i $wd/$species'_BetaIn_'$chr'_win'$win'.txt' -w $win -fold 

# add in 0 betas every 1kb or so where there is sequence without a snp (had reads, in chr, not in gap)
# the window starts and stops from angsd are different --DOES THIS MATTER????
# I'm taking the midpoint position as an artificial score of 0 for a snp


#get just the passing positions where the pos
#isn't already present in beta output 
#(will only happen if we didn't filter them out)
awk -F, 'FNR==NR {a[$1]=$0; next}; !($1 in a) {print $1,$2}' \
    $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'.txt' \
    $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_nohomo_passFILT.txt' \
    > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'addPos.txt'

#add the pass filt positions and sort
cat $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'addPos.txt' \
    $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'.txt' \
    | sort -s -n -k 1,1 \
    > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'comb.txt'


## Smooth windows using smwin and smstep 
## May want to update this to use previous window starts and stops
## $maffile$chr'pass_windows.txt'
python /seq/vgb/linda/bal_sel/scripts_2017/smooth_win.py $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'comb.txt' $smwin $smstep \
 > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'sm.txt'

#convert smoothed betas into window output for liftover #
cat $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'sm.txt' | python /seq/vgb/linda/bal_sel/scripts_2017/format_smooth_beta_for_liftover.py \
    > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_smoothed_toliftover.bed'


################################### do the liftover for this chr ################################
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


