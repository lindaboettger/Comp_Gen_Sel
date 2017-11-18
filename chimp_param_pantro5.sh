#Set human-specific variables for getting MAF for betascan

## set paths ##
species='chimp'
parentpath='/seq/vgb/linda/bal_sel/population_data'
mychrs='/seq/vgb/linda/bal_sel/lists/chimp_pantro5_chain_chrs.txt' ## updated on nov 16
#mychrs='/seq/vgb/linda/bal_sel/scripts_2017/chimp_test_chrs_pantro5.txt'
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
ref='/seq/vgb/references/chimp/pantro5/panTro5_EBV.fa'
this_chain=$chainpath/panTro5ToHg38.over.chain



## filters ## -- need to update read depth!!!
HWEpval=0.0217
SB=2.39 # fixed Nov 16
MinDepth=50 # fixed Nov 16
MaxDepth=125 # fixed Nov 16
minInd=6
minMapQ=30
minQ=20

#These were the cutoffs from get_read_depth.sh, but they look very wrong
#I plotted the data.. there is a very long tail so this 95% thing may not work
# I adjusted the above parameters 
#mean: 1829
#max cut: 2388215
#min cut: 95

#liftover
minmatch=0.1


