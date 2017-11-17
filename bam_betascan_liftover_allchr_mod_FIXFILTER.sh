#!/bin/bash


use GridEngine8


source $species_param

rm /seq/vgb/linda/bal_sel/out/$species'beta_out'
rm /seq/vgb/linda/bal_sel/err/$species'beta_err'

while read c; do
qsub -S /bin/bash -q vert -v chr=$c,new_species=$new_species,species_param=$species_param,win=$win,smwin=$smwin,smstep=$smstep,wd=$wd,filt_both_homo=$filt_both_homo \
 -M boettger@broadinstitute.org \
 -o /seq/vgb/linda/bal_sel/out/$species'beta_out' \
 -e /seq/vgb/linda/bal_sel/err/$species'beta_err' \
 -N $species'_'$c'_'$win  \
 -wd $wd \
 -p -10 \
 /seq/vgb/linda/bal_sel/scripts_2017/bam_betascan_liftover_byChr_mod_FIXFILTER.sh
done < $mychrs

### how many chrs are there?
nchr=$(wc -l < $mychrs)

#### what prefix do all chromosomes have?
if [ $species == 'gorilla' ]; then
  chrpre='CYU'
else
  chrpre='chr'
fi


# check that there are still  jobs with this outfile name is `ps -ef`; there will always
# be at least 1 because the original grep command is returned
while ! [ `qstat -r | grep $species'_'$chrpre | wc -l` == '0' ]
do 
  echo "I am waiting for chrs"
    sleep 60
done

echo "I think I am done!!"

# if nchrs not equal out files now, then die with error

if [ $species == 'YRI' ]; then
  nfiles=`ls $wd/*'_smoothed_toliftover.bed' -1 | wc -l`
else
  nfiles=`ls $liftout/$species'_to_hg38_'*'_win'$win'_smoothed.bed' -1 | wc -l`
fi


echo $nfiles
if ! [ $nchr == $nfiles ]; then 
    echo "Number of chroms ($nchr) does not match number of Liftover outfiles ($nfiles)"
    exit 1
else
    echo "Beta and Liftover runs complete, all $nchr outfiles present."
fi


if [ $species == 'YRI' ] 
  then
  #concatinate the chrs from pre-liftover files
    cat $wd/*smoothed_toliftover.bed \
      > $wd'/YRI_all_win'$win'_smoothed.bed'

    #get rid of underscores for betas
    awk '{sub(/_/," ",$5)};1' $wd'/YRI_all_win'$win'_smoothed.bed' > $wd'/YRI_all_win'$win'_smoothed.tmp'
    awk '{sub(/_/," ",$6)};1' $wd'/YRI_all_win'$win'_smoothed.tmp' > $wd'/YRI_all_win'$win'_smoothed.bed'
    rm $wd'/YRI_all_win'$win'_smoothed.tmp'

      
    #get rid of chr files
    #rm $wd/*toliftover.bed 

    else

      #concatinate post-liftover files
      cat $liftout/*'_win'$win'_smoothed.bed' > $liftout/$species'_to_hg38_win'$win'_smoothed_all'.bed

      #rm underscores
      awk '{sub(/_/," ",$5)};1' $liftout/$species'_to_hg38_win'$win'_smoothed_all'.bed > $liftout/$species'_to_hg38_win'$win'_smoothed_all'.tmp
      awk '{sub(/_/," ",$6)};1' $liftout/$species'_to_hg38_win'$win'_smoothed_all'.tmp > $liftout/$species'_to_hg38_win'$win'_smoothed_all'.bed
      rm $liftout/$species'_to_hg38_win'$win'_smoothed_all'.tmp


      #and unlifted stuff
      #cat $liftpath/$species'_to_hg38_'chr*'_win'$win_smoothed_Unmapped.bed > $liftpath/$species'_to_hg38_all_win'$win_smoothed_Unmapped.bed
    


fi

## get rid of files
#rm $liftpath/$species'_to_hg38_'chr*'_win'$win_Unmapped.bed
#rm $liftout/*chr*'_win'$win'.bed'
