#Set human-specific variables for getting MAF for betascan

## set paths ##
species='gorilla'
parentpath='/seq/vgb/linda/bal_sel/population_data'
#mychrs='/seq/vgb/linda/bal_sel/lists/large_gorilla_chr_testers_gorgor5.txt' 
mychrs='/seq/vgb/linda/bal_sel/chains/gorgor5_chain_chrs.txt'
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
ref='/seq/vgb/references/gorilla/gorGor5/gorGor5.fa'
this_chain=$chainpath/gorGor5ToHg38.over.chain



## filters ## -- need to update read depth!!!
HWEpval=0.0217
SBpval=0.0001 # this is the smallest number that snp stats keeps
MinDepth=200
MaxDepth=610
minInd=15
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


