import time
import os
import sys
from collections import Counter

if len(sys.argv) == 4:
	file1 = sys.argv[1]
	mdist_min = float(sys.argv[2])
	no_of_res = int(sys.argv[3])
else:
	print 'Analyse.py <align_output.txt> <mapp-min> <No.Of residue match>'
	sys.exit()
	

aline = open(file1, 'r').readlines()
arr = []
for line in aline:
	line = line.strip()
	l = line.split('\t')
	if len(l) == 4:
		if l[2] != 'None':
			l1 = l[2].split(' ')
			if len(l1) == 5:
				_, res, minim, _, _ = l1
				res = int(res.split('/')[0])
				minim = float(minim)
				if res >= no_of_res and minim > mdist_min:
					arr.append(l[0])
					
print Counter(arr).most_common(1)[0]






