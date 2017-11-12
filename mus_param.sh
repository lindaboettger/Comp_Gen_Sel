#Set human-specific variables for getting MAF for betascan

## set paths ##
parentpath='/seq/vgb/linda/bal_sel/from_oracle'
mychrs=$parentpath'/lists/mus_chr.txt'
bampath=$parentpath'/mus_bams/mus_bams_'$chr
outpath=$parentpath'/mus_out_beta'
bams=$bampath/mus
outs=$outpath/musoutFold_
betas=$parentpath/'mus_out_beta/mus_'
species='mus'
liftpath=$parentpath/$species'_out_beta'/$species'_lift'
chainpath=$parentpath'/chains'
liftout=$parentpath/$species'_out_beta/'$species'_out_lift'
this_chain=$chainpath'/mm9ToHg19.over.chain'

#reference
ref=/seq/references/Mus_musculus_assembly9/v1/Mus_musculus_assembly9.fasta


## filters ##
HWEpval=0.02
SBpval=0.0001 # this is the smallest number the snp stats keeps
MinDepth=50
MaxDepth=275
minInd=4
minMapQ=30
minQ=20

#liftover
minmatch=0.1

