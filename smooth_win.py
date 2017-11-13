import sys
import numpy as np 
#variables we read in
file = open(sys.argv[1])
win_size = int(sys.argv[2])
step_size = int(sys.argv[3])



start_pos = 0
end_pos = win_size

betas = []
pos = []


for line in file:
	pos_i, beta = line.split("\t")
	pos_i = int(pos_i)
	beta = float(beta)

	if pos_i < start_pos:
		print("SNP occurs before start of window!!")
		exit(1)
	if pos_i > end_pos:
		#calculate max for all snps in current window
		num_snps = len(pos)
		if num_snps > 0:
			max_beta = max(betas)
			#selectstring = df['ColumnD'].where(df['ColumnF'] == '#')
			max_pos = pos[np.argmax(betas)]
			mean_beta = sum(betas)/float(num_snps)
			print("{}\t{}\t{}\t{}\t{}\t{}".format(
					start_pos, end_pos, max_beta, 
					max_pos, mean_beta, num_snps))

		#update new window boundaries
		start_pos = start_pos + step_size
		end_pos = end_pos + step_size

		#remove snps outside of new window
		while len(pos) > 0 and pos[0] < start_pos:
			pos.pop(0)
			betas.pop(0)
	else:
		pos.append(pos_i)
		betas.append(beta)