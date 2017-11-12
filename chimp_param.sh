#Set human-specific variables for getting MAF for betascan

## set paths ##
parentpath='/seq/vgb/linda/bal_sel/from_oracle'
mychrs=$parentpath'/lists//chimp_chr.txt'
bampath=$parentpath'/chimp_bams/chimp_bams_'$chr
outpath=$parentpath'/chimp_out_beta'
bams=$bampath/chimp
outs=$outpath/chimpoutFold_
betas=$parentpath/'chimp_out_beta/chimp_'
species='chimp'
liftpath=$parentpath/$species'_out_beta'/$species'_lift'
chainpath=$parentpath'/chains'
liftout=$parentpath/$species'_out_beta/'$species'_out_lift'
this_chain=$chainpath'/panTro2ToHg19.over.chain'

#reference
ref=/seq/vgb/references/chimp/panTro2.fasta


## filters ##
HWEpval=0.0217
SBpval=0.0001 # this is the smallest number the snp stats keeps
MinDepth=15
MaxDepth=175
minInd=5
minMapQ=30
minQ=20

#liftover
minmatch=0.1


