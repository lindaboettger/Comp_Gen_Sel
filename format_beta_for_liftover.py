import sys
import re
import os

chrnum = str(os.environ["chr"]) #changed int to str
win = int(os.environ["win"])

for line in sys.stdin:
    pos, score = line.split()
    win_start = int(pos) - (win/2-1)
    win_end = int(pos) + (win/2-1)
    chr_name = str(chrnum) #don't want to add chr string

    converted_line = "\t".join([
            chr_name,
            str(win_start),
            str(win_end),
            '_'.join([chr_name, str(win_start),str(win_end)]),
            score])

    print converted_line
