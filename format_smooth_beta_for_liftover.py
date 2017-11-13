import sys
import re
import os

chrnum = str(os.environ["chr"]) #changed int to str
win = int(os.environ["win"])

for line in sys.stdin:
    win_start, win_end, high_score, high_pos, mean_score, nsnps = line.split()
    chr_name = str(chrnum) #don't want to add chr string


    converted_line = "\t".join([
        chr_name,
        str(win_start),
        str(win_end),
        '_'.join([chr_name, str(win_start),str(win_end)]),
        '_'.join([high_score, high_pos, mean_score, nsnps])])

    print converted_line

