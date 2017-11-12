#Set human-specific variables for getting MAF for betascan

## set paths ##
species='YRI'
parentpath='/seq/vgb/linda/bal_sel/population_data'
mychrs='/seq/vgb/linda/bal_sel/from_oracle/lists/YRI_chr_hg38.txt'
bampath=$parentpath'/bams/'$species
betapath=$parentpath'/betas/'$species
bamfilelist=$bampath/$species'.bamfilelist'
mafpath=$parentpath'/mafs/'$species
maffile=$mafpath/$species'_'
mafs=$mafpath/$species'_'
#outs=$outpath/YRIoutFold_ 
#betas=$parentpath/'YRI_out_beta/YRI_' do I use this?
liftpath=$parentpath'/lifts/'$species
liftout=$parentpath'/liftout/'$species
chainpath='/seq/vgb/linda/bal_sel/chains'

#reference
ref=/seq/references/Homo_sapiens_assembly38/v0/Homo_sapiens_assembly38.fasta


## below depth calculated from get_read_depth.sh
HWEpval=0.00001
SBpval=0.000001 # this is the smallest number the snp stats keeps
MinDepth=1
MaxDepth=619
minInd=22
minMapQ=30
minQ=20

