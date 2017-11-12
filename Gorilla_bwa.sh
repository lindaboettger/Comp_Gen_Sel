#!/bin/bash

indiv_id=$(perl -e '$ARGV[0] =~ /(.+).sra/ && print $1' $indv)
#indiv_id=Kowali1

## Use the samtools package
use GATK3
use BWA
use Samtools
picard='/home/unix/boettger/bin/picard-tools-2.0.1'

mypath='/seq/vgb/linda/bal_sel/population_data'
ref='/seq/vgb/references/gorilla/gorGor5/gorGor5.fa'

#Get fastq files from .sra
/seq/vgb/software/ncbi/bin/fastq-dump -I --split-3 $mypath/fastq/gorilla/$indv --gzip --outdir $mypath/fastq/gorilla


#make read names the same
rm_readname='s/\w+\d+\.\d+\.\d+\s+HWI/HWI/'
zcat $mypath/fastq/gorilla/$indiv_id'_1.fastq.gz' | perl -pe $rm_readname | gzip > $mypath/fastq/gorilla/$indiv_id'_edit_1.fastq.gz'
zcat $mypath/fastq/gorilla/$indiv_id'_2.fastq.gz' | perl -pe $rm_readname | gzip > $mypath/fastq/gorilla/$indiv_id'_edit_2.fastq.gz'

rm $mypath/fastq/gorilla/$indiv_id'_1.fastq.gz'
rm $mypath/fastq/gorilla/$indiv_id'_2.fastq.gz'

#/seq/vgb/linda/bal_sel/population_data/bams/gorilla

bwa mem $ref $mypath/fastq/gorilla/$indiv_id'_edit_1.fastq.gz' \
     $mypath/fastq/gorilla/$indiv_id'_edit_2.fastq.gz' \
     | samtools view -bS \
    > $mypath/bams/gorilla/$indiv_id'-pe.bam'


# bwa mem -t 3 $ref $mypath/fastq/gorilla/$indiv_id'_edit_1.fastq.gz' \
#      $mypath/fastq/gorilla/$indiv_id'_edit_2.fastq.gz' \
#      | samtools view -@ 3 -bS \
#     > $mypath/bams/gorilla/$indiv_id'-pe.bam'

# #sort
samtools sort -T $mypath/bams/gorilla/tmp/$indiv_id'.sorted' -o $mypath/bams/gorilla/$indiv_id'.sorted.bam' $mypath/bams/gorilla/$indiv_id'-pe.bam'

#samtools sort -T $mypath/bams/gorilla/tmp/$indiv_id'.sorted' -@ 3 -o $mypath/bams/gorilla/$indiv_id'.sorted.bam' $mypath/bams/gorilla/$indiv_id'-pe.bam'

rm $mypath/bams/gorilla/$indiv_id'-pe.bam'




