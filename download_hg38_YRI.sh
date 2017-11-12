#!/bin/bash


use Samtools

#download from 1kg 
samtools view -bh 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/YRI/'$indiv'/alignment/'$indiv'.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram' \
	 -T /seq/references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta \
	 -o '/seq/vgb/linda/bal_sel/population_data/bams/YRI/'$indiv'.bam'

samtools index '/seq/vgb/linda/bal_sel/population_data/bams/YRI/'$indiv'.bam'

