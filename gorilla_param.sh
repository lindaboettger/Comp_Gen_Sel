#Set human-specific variables for getting MAF for betascan

## set paths ##
parentpath='/seq/vgb/linda/bal_sel/from_oracle'
mychrs=$parentpath'/lists/gorilla_largechrs.txt'
bampath=$parentpath'/gorilla_bams/gorilla_bams'$chr
outpath=$parentpath'/gorilla_out_beta'
bams=$bampath/gorilla
outs=$outpath/gorillaoutFold_
betas=$parentpath/'gorilla_out_beta/gorilla_'
species='gorilla'
liftpath=$parentpath/$species'_out_beta'/$species'_lift'
chainpath=$parentpath'/chains'
liftout=$parentpath/$species'_out_beta/'$species'_out_lift'


#reference
ref='/seq/vgb/references/gorilla/gorGor4.fa'

##filters
#HWE and SB are what i had before, but these seem pretty liberal.... maybe because many individuals (27dip)..?
HWEpval=0.00001
SBpval=0.000001 # this is the smallest number the snp stats keeps
MinDepth=100
MaxDepth=550 
minInd=22
minMapQ=30
minQ=20

#liftover, not used here
minmatch_gorgor=0.5
minmatch_gorhum=0.1


