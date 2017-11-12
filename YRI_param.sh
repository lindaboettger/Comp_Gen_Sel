#Set human-specific variables for getting MAF for betascan

## set paths ##
parentpath='/seq/vgb/linda/bal_sel/from_oracle'
mychrs=$parentpath'/lists/YRI_chr.txt'
bampath=$parentpath'/YRI_bams/YRI_bams_chr'$chr
outpath=$parentpath'/YRI_out_beta'
bams=$bampath/YRIrentschr
outs=$outpath/YRIoutFold_
betas=$parentpath/'YRI_out_beta/YRI_'
species='YRI'
liftpath=$parentpath/$species'_out_beta'/$species'_lift'
chainpath=$parentpath'/chains'

#reference
ref=/seq/vgb/linda/references/hg19/Homo_sapiens_assembly19.fasta


## filters ##
HWEpval=0.00001
SBpval=0.000001 # this is the smallest number the snp stats keeps
MinDepth=70
MaxDepth=700
minInd=22
minMapQ=30
minQ=20

