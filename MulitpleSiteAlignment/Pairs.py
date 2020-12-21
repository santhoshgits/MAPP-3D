import time
import os
import sys

if len(sys.argv) == 2:
	fold = sys.argv[1]
else:
	print 'Pairs.py <Site-Folder>'
	sys.exit()
	
dire = os.getcwd()
out = open('PairList.txt', 'w')
for i in os.listdir(dire+'/'+fold):
	for j in os.listdir(dire+'/'+fold):
		out.write(i+'\t'+j+'\n')
out.close()
	














