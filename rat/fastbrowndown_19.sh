#!/bin/bash

#ebi http://www.ebi.ac.uk/ena/data/view/ERP001276&display=html
#these bams were mapped to the rn5 reference genome

outdir='/seq/vgb/linda/bal_sel/population_data/fastq/rat'

wget -O $outdir/ERR319179_1_rat_55.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319179/ERR319179_1.fastq.gz
wget -O $outdir/ERR319179_2_rat_55.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319179/ERR319179_2.fastq.gz
wget -O $outdir/ERR319180_1_rat_56.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319180/ERR319180_1.fastq.gz
wget -O $outdir/ERR319180_2_rat_56.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319180/ERR319180_2.fastq.gz
wget -O $outdir/ERR319181_1_rat_57.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319181/ERR319181_1.fastq.gz
wget -O $outdir/ERR319181_2_rat_57.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319181/ERR319181_2.fastq.gz
wget -O $outdir/ERR319182_1_rat_58.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319182/ERR319182_1.fastq.gz
wget -O $outdir/ERR319182_2_rat_58.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319182/ERR319182_2.fastq.gz
wget -O $outdir/ERR319183_1_rat_59.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319183/ERR319183_1.fastq.gz
wget -O $outdir/ERR319183_2_rat_59.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR319/ERR319183/ERR319183_2.fastq.gz