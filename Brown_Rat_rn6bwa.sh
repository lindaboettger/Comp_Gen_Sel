#!/bin/bash


#ERR316487_1_rat_51.fastq.gz


## Use the samtools package
use GATK3
use BWA
use Samtools
picard='/home/unix/boettger/bin/picard-tools-2.0.1'

mypath='/seq/vgb/linda/bal_sel/population_data'
ref='/seq/vgb/references/brown_rat/rn6.fa' #this is rn6


bwa mem $ref $mypath/fastq/rat2/$indv_1'.fastq.gz' \
     $mypath/fastq/rat2/$indv_2'.fastq.gz' \
     | samtools view -bS \
    > $mypath/fastq/rat2/$indv'-pe.bam'

# #sort
samtools sort -T $mypath/bams/rat/tmp/$indv'.sorted' -o $mypath/bams/rat/$indv'.sorted.bam' $mypath/fastq/rat2/$indv'-pe.bam'

#rm $mypath/new_bams/$indiv_id'-pe.bam'




