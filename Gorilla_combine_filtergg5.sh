use Samtools #samtools 1.3
use Java-1.8
picard='/home/unix/boettger/bin/picard-tools-2.0.1'
GATK='/broad/software/free/Linux/redhat_6_x86_64/pkgs/GATK3-2.2'
use Picard-Tools

#set path and ref
mypath='/seq/vgb/linda/bal_sel/population_data/bams/gorilla'
ref='/seq/vgb/references/gorilla/gorGor5/gorGor5.fa'


####NEED TO UNCOMMENT BELOW ###
mybams="$mypath/$indv*"

### merge sorted bams for this individual
samtools merge $mypath/$indv'_merged.bam' $mybams

samtools index $mypath/$indv'_merged.bam'

rm the old bams
rm $indv*.sorted.bam

## add read group info
#"The GATK no longer supports SAM files without read groups"
java -Xmx2G -jar $picard/picard.jar AddOrReplaceReadGroups \
    INPUT=$mypath/$indv'_merged.bam' \
    OUTPUT=$mypath/$indv'_mergedRG.bam' RGLB=Gorlib \
    RGPL=Illumina RGPU=1 RGSM=$indv RGDS=AD

#mark duplicates
java -jar $picard/picard.jar MarkDuplicates \
        INPUT=$mypath/$indv'_mergedRG.bam' \
        OUTPUT=$mypath/$indv'_mkdup.bam' \
        METRICS_FILE=$mypath/$indv'_mkdup.metrics'

samtools index $mypath/$indv'_mkdup.bam'

#get indel intervals to realign
java -jar $GATK/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -I $mypath/$indv'_mkdup.bam' \
 -R $ref \
 -o $mypath/$indv'_realign.intervals'
   
#indel realign using intervals
java -jar $GATK/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -I $mypath/$indv'_mkdup.bam' \
 -R $ref \
 -o $mypath/$indv'_realign_mkdup.bam' \
 -targetIntervals $mypath/$indv'_realign.intervals'
#  #-baq RECALCULATE apparently this is bad to do when you do indel realign
  
#sort and index final bamfile
samtools sort -T $mypath/tmp/$indv'_out.sorted.bam' -o $mypath/$indv'_out.sorted.bam' $mypath/$indv'_realign_mkdup.bam'
samtools index $mypath/$indv'_out.sorted.bam'

# remove intermediate files
rm $indv*_merged.bam
rm $indv*_merged.bam.bai
rm $mypath/$indv'_mkdupRG.bam'
rm $mypath/$indv'_mkdupRG.bam.bai'
rm $mypath/$indv'_realign_mkdup.bam'
rm $mypath/$indv'_realign.intervals'





