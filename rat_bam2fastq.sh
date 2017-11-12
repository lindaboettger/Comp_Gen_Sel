#!/bin/bash
use Samtools
ues Perl-5.8

#echo "hello"

bampath=/seq/vgb/linda/bal_sel/from_oracle/brown_rat_bams
fastqpath=/seq/vgb/linda/bal_sel/population_data/fastq/rat2
mkdir -p $fastqpath/tmp
# #$ERR316550_rat_59.sorted.bam


# bam to fastq
#samtools bam2fq $bampath/$ERRrat'.sorted.bam' | gzip > $fastqpath/$ERRrat.fastq.gz


### separate /1 and /2, make the read names the same and sort the reads
rm_readsample='s/(@.+)\/\d/$1/'
zcat $fastqpath/$ERRrat.fastq.gz  | grep '^@.*/1$' -A 3 --no-group-separator | perl -pe $rm_readsample \
	| paste - - - - | sort -T "$fastqpath/tmp" -k1,1 -t " " | tr "\t" "\n" | gzip > $fastqpath/$ERRrat'_1.fastq.gz'


zcat $fastqpath/$ERRrat.fastq.gz | grep '^@.*/2$' -A 3 --no-group-separator | perl -pe $rm_readsample \
	| paste - - - - | sort -T "$fastqpath/tmp" -k1,1 -t " " | tr "\t" "\n" | gzip > $fastqpath/$ERRrat'_2.fastq.gz'



#make read samples the same
#@HJ2CJCCXX160105:7:1212:9973:43167/2 or @HJ2CJCCXX160105:7:1212:9973:43167/1
