#!/bin/bash
use Samtools

bampath=/seq/vgb/linda/bal_sel/from_oracle/brown_rat_bams
fastqpath=/seq/vgb/linda/bal_sel/population_data/fastq/rat2
#$ERR316550_rat_59.sorted.bam

samtools bam2fq $bampath/$ERRrat'.sorted.bam' > $fastqpath/$ERRrat.fastq 
zcat $fastqpath/$ERRrat.fastq  | grep '^@.*/1$' -A 3 --no-group-separator | gzip > $fastqpath/$ERRrat_1.fastq 
zcat $fastqpath/$ERRrat.fastq  | grep '^@.*/2$' -A 3 --no-group-separator | gzip > $fastqpath/$ERRrat_2.fastq 
