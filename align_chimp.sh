#!/bin/bash

## Use the samtools package
use GATK3
use BWA
use Samtools
picard='/home/unix/boettger/bin/picard-tools-2.0.1'

mypath='/seq/vgb/linda/bal_sel/population_data'
ref='/seq/vgb/references/chimp/pantro5/panTro5_EBV.fa' 
#/seq/vgb/linda/bal_sel/population_data/fastq/chimp

#make read names the same
rm_readname='s/^([@+]\w+\d+\.\d+)\.\d(\s)/$1$2/'
zcat $mypath/fastq/chimp/$sra_id'_1.fastq.gz' | perl -pe $rm_readname | gzip > $mypath/fastq/chimp/$sra_id'_edit_1.fastq.gz'
zcat $mypath/fastq/chimp/$sra_id'_2.fastq.gz' | perl -pe $rm_readname | gzip > $mypath/fastq/chimp/$sra_id'_edit_2.fastq.gz'

# rm $mypath/fastq/chimp/$sra_id'_1.fastq.gz'
# rm $mypath/fastq/chimp/$sra_id'_2.fastq.gz'



bwa mem -t 26 $ref $mypath/fastq/chimp/$sra_id'_edit_1.fastq.gz' \
     $mypath/fastq/chimp/$sra_id'_edit_2.fastq.gz' \
     | samtools view -bS \
    > $mypath/bams/chimp/$sra_id'-pe.bam'

# #sort
samtools sort -T $mypath/bams/chimp/tmp/$sra_id'.sorted' -o $mypath/bams/chimp/$sra_id'.sorted.bam' $mypath/bams/chimp/$sra_id'-pe.bam'

rm $mypath/bams/chimp/$sra_id'-pe.bam'





