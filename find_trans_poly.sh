#!/bin/bash


# Get species-specific variables #
source $species_param

tolift=$parentpath/trans_poly/$species'_maf_tolift'
liftout=$parentpath/trans_poly/$species'_maf_liftout'

mkdir -p $tolift
mkdir -p $liftout

################ MAF to lift input ################
maf=$mafpath/$species'_'$chr'_filtered.mafs'


## get 100bp on either side
# ## they aren't whole numbers because genotype probs used to calculate freq
#chr
#snp - 100bp
#snp + 100 bp
#MAF
#snp pos, major, minor, n indivs


if [ $type == "region" ]; then

 #### for regions ####
awk '{if (NR!=1) print $1, $2-100, $2+100, $6, $1":"$2"_"$3"_"$4"_"$8}' $maf \
    > $tolift/$species'_'$chr'.bed'

 
 #### I did this back when using hg19
  #add "chr if it is YRI"
  # if [ $species == 'YRI' ]; then
  #   awk '{if (NR!=1) print "chr"$1, $2-100, $2+100, $6, $1":"$2"_"$3"_"$4"_"$8}' $maf \
  #     > $tolift/$species'_'$chr'.bed'
  #   else 
  #     awk '{if (NR!=1) print $1, $2-100, $2+100, $6, $1":"$2"_"$3"_"$4"_"$8}' $maf \
  #     > $tolift/$species'_'$chr'.bed'
  # fi


  ### for snps
  #awk '{if (NR!=1) print $1, $2, $2+1, $6, $3"_"$4"_"$8}' $maf > $tolift/$species'_SNP_'$chr'.bed'



  ############### for regions ###########
  # liftover
  #-minMatch default of 0.95)
  #"required that at least 95% of the flanking 100 bp in both directions lift over to the other reference genome"
  if ! [ $species == 'YRI' ]; then
    /seq/vgb/linda/bin/liftOver \
      -minMatch=0.95 \
      $tolift/$species'_'$chr'.bed' \
      $this_chain \
      $liftout/$species'_to_hg38_'${chr}.bed \
      $liftpath/$species'_to_hg38_'${chr}_Unmapped.bed

      echo "I am not a gorilla!!!"

      else
      if [ $species == 'YRI' ]
      then

      echo "I am human"
          #### human to chimp ####
        /seq/vgb/linda/bin/liftOver \
        -minMatch=0.95 \
        $tolift/$species'_'$chr'.bed' \
        $chainpath/hg38ToPanTro5.over.chain \
        $liftout/YRI_to_PanTro5_${chr}.bed \
        $liftpath/YRI_to_PanTro5_${chr}_Unmapped.bed
      fi
  fi
fi

if [ $type == "snp" ]; then

# #add "chr if it is YRI" --- had to do this when it was hg19
#   if [ $species == 'YRI' ]; then
#     awk '{if (NR!=1) print "chr"$1, $2, $2+1, $6, $1":"$2"_"$3"_"$4"_"$8}' $maf \
#       > $tolift/$species'_SNP_'$chr'.bed'
#     else 
# awk '{if (NR!=1) print $1, $2, $2+1, $6, $1":"$2"_"$3"_"$4"_"$8}' $maf \
#   > $tolift/$species'_SNP_'$chr'.bed'
#   fi


awk '{if (NR!=1) print $1, $2, $2+1, $6, $1":"$2"_"$3"_"$4"_"$8}' $maf \
   > $tolift/$species'_SNP_'$chr'.bed'


  
  ###### for SNPs ###############
  if ! [ $species == 'YRI' ]; then
    /seq/vgb/linda/bin/liftOver \
      -minMatch=0.95 \
      $tolift/$species'_SNP_'$chr'.bed' \
      $this_chain \
      $liftout/$species'_to_hg38_SNP_'${chr}.bed \
      $liftpath/$species'_to_hg38_'${chr}_SNPUnmapped.bed

      echo "I am a chimp!!"

      else
      if [ $species == 'YRI' ]
      then
    

      echo "I am human"
          #### gorgor4 to gorgor3to human ####
        /seq/vgb/linda/bin/liftOver \
        -minMatch=0.95 \
        $tolift/$species'_SNP_'$chr'.bed' \
        $chainpath/hg38ToPanTro5.over.chain \
        $liftout/YRI_to_PanTro5_SNP_${chr}.bed \
        $liftpath/YRI_to_PanTro5_SNP_${chr}_Unmapped.bed


      fi
  fi
fi








#combine


#sort by human chrom



#!/bin/bash



# ############### run on each chromosome separately ###############
# awk 'BEGIN {
#   FS = OFS = "\t"
#   }
# NR == FNR {
#   # while reading the 1st file
#   # store its records in the array f
#   f[$2] = $0
#   next
#   }
# $1 in f {
#   # when match is found
#   # print all values
#   print f[$2], $0
#   }' file1 file2 > 