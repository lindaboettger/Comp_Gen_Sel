use Python-3.4
#input variables
#chr=22

#directories
anc_dir=/seq/vgb/linda/ancestral_bases/human_chimp/human_ancestor_GRCh37_e59
snp_dir=/seq/vgb/linda/bal_sel/trans_poly
script_dir=/seq/vgb/linda/bal_sel/scripts_2017

## snp file
my_snps=$snp_dir/trans_snps_unfiltered.txt
chr_snps=$snp_dir'/trans_snps_unfiltered_'$chr'.txt'
outfile='/seq/vgb/linda/bal_sel/trans_poly/human_chimp/ancestrals_'$chr'.txt'

#make chr-specific query snp file (I think this is faster than finding them separately with python)

grep "^$chr " $my_snps > $chr_snps

#ancestral file for this chr
anc_snps=$anc_dir'/human_ancestor_'$chr'.fa'

###### fix chromosome names
carrot='>'
fa_header=$carrot$chr

sed -i "1s/.*/$fa_header/" $anc_snps 



### call python script that calls samtools
python $script_dir/get_ancestral_base_byChr.py $anc_snps $chr_snps $chr > $outfile


#samtools faidx human_ancestor_22.fa 22:50693534-50693534