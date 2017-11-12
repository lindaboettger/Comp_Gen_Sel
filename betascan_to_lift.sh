use Python-2.7 


source $species_param
betascan=/seq/vgb/linda/bin/BetaScan-master/BetaScan.py


# -p: Value of p (default: 20)
#run betascan

echo "I am about to start betascan"
python $betascan -i $wd/$species'_BetaIn_'$chr'_win'$win'.txt' -w $win -fold 




#convert betas into window output for liftover #
cat $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'.txt' | python /seq/vgb/linda/bal_sel/scripts_2017/format_beta_for_liftover.py \
    > $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_toliftover.bed'



### do the liftover for this chr
#do the liftover if this isn't human
### NOTE GORILLA WILL REQUIRE 2 LIFTOVERS 
echo $species

if ! [ $species == 'YRI' ] && ! [ $species == 'gorilla' ]; then
  /seq/vgb/linda/bin/liftOver \
    -minMatch=$minmatch \
    $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_toliftover.bed' \
    $this_chain \
    $liftout/$species'_to_hg19_'${chr}'_win'$win'.bed' \
    $liftpath/$species'_to_hg19_'${chr}_Unmapped.bed

    echo "I am not a gorilla!!!"

    else
    if [ $species == 'gorilla' ]
    then

    echo "I am a fucking gorilla"
        #### gorgor4 to gorgor3to human ####
      /seq/vgb/linda/bin/liftOver \
      -minMatch=0.5 \
      $wd/'Betas_'$species'_BetaIn_'$chr'_win'$win'_toliftover.bed' \
      $chainpath/gorGor4ToGorGor3.over.chain \
      $liftpath/gorgor4_to_gorgor3_${chr}.bed \
      $liftpath/gorgor4_to_gorgor3_${chr}_Unmapped.bed


      #### gorilla to human ####
      /seq/vgb/linda/bin/liftOver \
      -minMatch=0.1 \
      $liftpath/gorgor4_to_gorgor3_${chr}.bed \
      $chainpath/gorGor3ToHg19.over.chain \
	  $liftout/$species'_to_hg19_'${chr}'_win'$win'_.bed' \
      $liftpath/$species'_to_hg19_'${chr}_Unmapped.bed

    fi
fi
