import sys
import re
import os
import subprocess

from Bio import SeqIO

# #set variables
# anc_snps = '/seq/vgb/linda/ancestral_bases/human_chimp/human_ancestor_GRCh37_e59/human_ancestor_22.fa'
# chr_snps = '/seq/vgb/linda/bal_sel/trans_poly/trans_snps_unfiltered_22.txt'
# outfile = /seq/vgb/linda/bal_sel/trans_poly
# chrnum = '22'

anc_snps = sys.argv[1]
chr_snps = sys.argv[2]
chrnum = sys.argv[3]
#outfile = sys.argv[4]


# Open chromosome
chr_fasta = next(SeqIO.parse(open(anc_snps), 'fasta'))

with open(chr_snps) as f:
    for line in f:
        # get the bp position (last element of split on space)
        # subtract one b/c base needs to be pythonic index
        bp = int(line.split(" ")[-1].rstrip()) - 1
        anc_base = chr_fasta[bp]
        #make output line and write to file
        print(chrnum, bp + 1, anc_base, sep='\t')

## HOW DO I USE SAMTOOLS???
### call samtools and write output base to variable
# Popen('use Samtools', shell=True)

# ## for each snp in chr_snps 
# with open(chr_snps) as f:
#     for line in f:
#       #get the bp position (last element of split on space)
#         bp = s.split(" ")[-1:]
#         pos = chrnum + ":" + bp + "-" + bp
#         #call samtools and save output to variable
#         proc = subprocess.Popen('samtools faidx anc_snps pos', stdout=subprocess.PIPE)
#       anc_base = proc.stdout.read()
#       #make output line and write to file
#       outline = print(chrnum, anc_base, sep='\t')
#       FILE.writelines(outline) 


# #close outfile when all lines are written
# FILE.close()
# chr_snps.close()
