use Samtools #samtools 1.3
use Java-1.8
picard='/home/unix/boettger/bin/picard-tools-2.0.1'
GATK='/broad/software/free/Linux/redhat_6_x86_64/pkgs/GATK3-2.2'

#set path and ref
mypath='/seq/vgb/linda/bal_sel/from_oracle/brown_rat_bams'
ref='/seq/references/Rattus_norvegicus/v0/Rattus_norvegicus.fasta'


#identify all bam files for this individual
#indv=sample_56
mybams="$mypath/*$indv*"

### sort (already done)
# SAVEIFS=$IFS
# IFS=$(echo -en "\n\b")
# # set bams and sort them 
# for f in $mybams
# do
#   #echo "$f"
#   ## get name of this bam
#   thisbam=$(perl -e '$ARGV[0] =~ /.+brown_rat_bams\/(.+).bam/ && print $1' $f)
#   # #sort
#   samtools sort -T $mypath/tmp/$thisbam'.sorted' -@ 8 -o $mypath/$thisbam'.sorted.bam' $f

# done
# # restore $IFS
# IFS=$SAVEIFS

## get list of sorted bams
mysorted="$mypath/*$indv.sorted.bam"


### merge sorted bams for this individual
samtools merge $mypath/$indv'_merged.bam' $mysorted

## to fix incorrect mate pair problems, 
# e.g. 'Mapped mate should have mate reference name'
# we will sort by read name and run samtools fixmate
# on the merged bam. This will properly point the 
# unmapped reads to their mapped mates and solve
# Picard/GATK errors.

samtools sort -n -T $mypath/tmp/$indv'.sorted' -@ 8 -o $mypath/$indv'.tofixsorted.bam' $mypath/$indv'_merged.bam'

samtools fixmate $mypath/$indv'.tofixsorted.bam' $mypath/$indv'_fixed.bam'

samtools sort -T $mypath/tmp/$indv'.sorted' -@ 8 -o $mypath/$indv'_fixedsorted.bam' $mypath/$indv'_fixed.bam'

samtools index $mypath/$indv'_fixedsorted.bam'

### Note: if I ever re-merge from scratch, we should run
# the command right from the merge like this:
# samtools merge -n - $mysorted \
#     | samtools fixmate ...
#     | ... > $mypath/$indv'_fixed.bam'

# ## add read group info
# #"The GATK no longer supports SAM files without read groups"
java -Xmx2G -jar $picard/picard.jar AddOrReplaceReadGroups \
    INPUT=$mypath/$indv'_fixedsorted.bam' \
    OUTPUT=$mypath/$indv'_mergedRG.bam' RGLB=Ratlib \
    VALIDATION_STRINGENCY=LENIENT \
    RGPL=Illumina RGPU=1 RGSM=$indv RGDS=AD

#mark duplicates
java -jar $picard/picard.jar MarkDuplicates \
    INPUT=$mypath/$indv'_mergedRG.bam' \
    OUTPUT=$mypath/$indv'_mkdup.bam' \
    VALIDATION_STRINGENCY=LENIENT \
    METRICS_FILE=$mypath/$indv'_mkdup.metrics'

samtools index $mypath/$indv'_mkdup.bam'

#get indel intervals to realign
java -jar $GATK/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -I $mypath/$indv'_mkdup.bam' \
 -R $ref \
 -o $mypath/$indv'_realign.intervals'
   
#indel realign using intervals and recalculate BAQ
java -jar $GATK/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -I $mypath/$indv'_mkdup.bam' \
 -R $ref \
 -o $mypath/$indv'_realign_mkdup.bam' \
 -targetIntervals $mypath/$indv'_realign.intervals'
 #-baq RECALCULATE apparently this is bad to do when you do indel realign
  
#sort and index final bamfile
samtools sort -o $mypath/finished_bams/$indv'_out.sorted.bam' $mypath/$indv'_realign_mkdup.bam'
samtools index $mypath/finished_bams/$indv'_out.sorted.bam'

# # ### remove intermediate files
# rm $indv*_merged.bam
# rm $indv*_merged.bam.bai
# rm $indv*_fixed.bam
# rm $indv*_fixed.bam.bai
# rm $mypath/$indv'_mkdupRG.bam'
# rm $mypath/$indv'_mkdupRG.bam.bai'
# rm $mypath/$indv'_realign_mkdup.bam'
# rm $indv*.sorted.bam





