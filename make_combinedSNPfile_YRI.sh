#libraries
use GCC-5.1
use VCFtools

## set variables
parentpath='/seq/vgb/linda/bal_sel'
ancestral_seq='/seq/vgb/linda/bal_sel/ancestral/chimpHg19.fa'
bampath=$parentpath'/from_oracle/YRI_bams/YRI_bams_chr'$chr
outpath=$parentpath'/combinedSNPfiles/YRI'$chr

## make out directory 
mkdir -p $outpath

### generate saf file GATK genotype likelihoods
#P = number of threads
#doMajorMinor 5 forces the major allele to be the ancestral
#similar to example http://www.popgen.dk/angsd/index.php/ANGSD
#tried with two versions of angsd "29 argument  -setMinDepth is unknown will exit"
#"29 argument   -setMaxDepth is unknown will exit"
/seq/vgb/linda/bin/Feb2017_angsd/angsd/angsd \
    -bam $bampath'/YRIrentschr'$chr'.bamfilelist' \
    -anc $ancestral_seq \
    -out $outpath/'YRI_'$chr'.snpfile' \
    -doSnpStat 1 \
    -nThreads 15 \
    -snp_pval 1e-6 \
    -minInd 22 \
    -GL 2 \
    -domaf 1 \
    -doMajorMinor 5 \
    -hwe_pval 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -uniqueOnly 1 \
    -minMapQ 20 \
    -minQ 20

    #-setMinDepth 70 \
    #    -setMaxDepth 700 \

