import os
import shutil
import sys

if len(sys.argv) == 2:
	fold = sys.argv[1]
else:
	print "pdb_res.py <input_folder>"
	sys.exit()


dire = os.getcwd()

a1 = dire+"/"+fold
def a():
	dic = {}
	out = open("PDBSize.txt",'w')
	for i in os.listdir(a1):
		#print i
		aline = open(a1+"/"+i, 'r').readlines()
		c = 0
		for line in aline:
			line = line.strip()
			if line[:4] == "ATOM":
				if line[13:15] == "CA":
					c += 1
		#print c
		dic[i] = c
		out.write(i+"\t"+str(c)+"\n")
	out.close()
	return dic
	#out.close()
dic = a()


