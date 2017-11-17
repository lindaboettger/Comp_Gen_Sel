## set paths ##
species='mus'
parentpath='/seq/vgb/linda/bal_sel/population_data'
mychrs='/seq/vgb/linda/bal_sel/lists/mus_mm10_chain_chrs.txt'
bampath=$parentpath'/bams/'$species
betapath=$parentpath'/betas/'$species
bamfilelist=$bampath/$species'.bamfilelist'
mafpath=$parentpath'/mafs/'$species
maffile=$mafpath/$species'_'
mafs=$mafpath/$species'_'
liftpath=$parentpath'/lifts/'$species
liftout=$parentpath'/liftout/'$species
chainpath='/seq/vgb/linda/bal_sel/chains'
#reference
ref='/seq/references/Mus_musculus_assembly10/v0/Mus_musculus_assembly10.fasta'
this_chain=$chainpath/mm10ToHg38.over.chain


## filters ##
HWEpval=0.05
SB=1.01 # fixed nov 16
MinDepth=130 # fixed nov 16
MaxDepth=230 # fixed nov 16
minInd=5
minMapQ=30
minQ=20

#liftover
minmatch=0.1

