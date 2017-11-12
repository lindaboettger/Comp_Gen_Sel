outdir=/seq/vgb/linda/bal_sel/population_data/fastq/chimp


#get sra file
#wget -O $outdir/$ERRid'.sra' ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR225/$ERRid/$ERRid'.sra'

#convert to fastq
/seq/vgb/software/ncbi/bin/fastq-dump -I --split-3 $outdir/$ERRid'.sra' --gzip --outdir $outdir

