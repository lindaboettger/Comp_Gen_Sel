#!/bin/bash

#get sraid without file extension
sra_id=$(echo $myfile | perl -pe 's/(SRR\d+).sra/$1/')



## Use the samtools package
use GATK3
use BWA
use Samtools
picard='/home/unix/boettger/bin/picard-tools-2.0.1'

mypath='/seq/vgb/linda/bal_sel/from_oracle/'
ref='/seq/vgb/references/orangutan/ponAbe2.fa' 


Get fastq files from .sra
/seq/vgb/software/ncbi/bin/fastq-dump -I --split-3 $mypath/orangutan_sra/$myfile --gzip --outdir $mypath/orangutan_fastq

#make read names the same
#ERROR: paired reads have different names: "SRR1518520.1.1", "SRR1518520.1.2"
rm_readname='s/^([@+]\w+\d+\.\d+)\.\d(\s)/$1$2/'
zcat $mypath/orangutan_fastq/$sra_id'_1.fastq.gz' | perl -pe $rm_readname | gzip > $mypath/orangutan_fastq/$sra_id'_edit_1.fastq.gz'
zcat $mypath/orangutan_fastq/$sra_id'_2.fastq.gz' | perl -pe $rm_readname | gzip > $mypath/orangutan_fastq/$sra_id'_edit_2.fastq.gz'

rm $mypath/orangutan_fastq/$sra_id'_1.fastq.gz'
rm $mypath/orangutan_fastq/$sra_id'_2.fastq.gz'



bwa mem -t 26 $ref $mypath/orangutan_fastq/$sra_id'_edit_1.fastq.gz' \
     $mypath/orangutan_fastq/$sra_id'_edit_2.fastq.gz' \
     | samtools view -bS \
    > $mypath/orangutan_bams/$sra_id'-pe.bam'

# #sort
samtools sort -T $mypath/orangutan_bams/tmp/$sra_id'.sorted' -o $mypath/orangutan_bams/$sra_id'.sorted.bam' $mypath/orangutan_bams/$sra_id'-pe.bam'

rm $mypath/orangutan_bams/$sra_id'-pe.bam'





