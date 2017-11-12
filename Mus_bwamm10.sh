#!/bin/bash

indiv_id=$(perl -e '$ARGV[0] =~ /.*_.*_(\w+)/ && print $1' $indv)

## Use the samtools package
use GATK3
use BWA
use Samtools


parentpath='/seq/vgb/linda/bal_sel/population_data'
ref='/seq/references/Mus_musculus_assembly10/v0/Mus_musculus_assembly10.fasta'


bwa mem $ref $parentpath/fastq/mus/$indv'_1_1_sequence.fq.gz' \
     $parentpath/fastq/mus/$indv'_1_2_sequence.fq.gz' \
     | samtools view -bS \
    > $parentpath/bams/mus/$indiv_id'-pe.bam'

#sort
samtools sort -T $parentpath/bams/mus/tmp/$indiv_id'.sorted' -o $parentpath/bams/mus/$indiv_id'.sorted.bam' $parentpath/bams/mus/$indiv_id'-pe.bam'
