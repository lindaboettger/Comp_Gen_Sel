use Samtools

use Java-1.8
picard='/home/unix/boettger/bin/picard-tools-2.0.1'
GATK='/broad/software/free/Linux/redhat_6_x86_64/pkgs/GATK3-2.2'

#set path and ref
mypath='/seq/vgb/linda/bal_sel/population_data/bams/chimp'
ref='/seq/vgb/references/chimp/pantro5/panTro5_EBV.fa'


samtools index $mypath/$indv'_sorted.bam'

### Note: if I ever re-merge from scratch, we should run
# the command right from the merge like this:
# samtools merge -n - $mysorted \
#     | samtools fixmate ...
#     | ... > $mypath/$indv'_fixed.bam'

# ## add read group info
# #"The GATK no longer supports SAM files without read groups"
java -Xmx2G -jar $picard/picard.jar AddOrReplaceReadGroups \
    INPUT=$mypath/$indv'.sorted.bam' \
    OUTPUT=$mypath/$indv'_mergedRG.bam' RGLB=chilib \
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
samtools sort -o $mypath/$indv'_out.sorted.bam' $mypath/$indv'_realign_mkdup.bam'
samtools index $mypath/$indv'_out.sorted.bam'
