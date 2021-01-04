import time
import random
import math
import numpy as np
import sys
from collections import Counter, defaultdict
import copy
import re


if len(sys.argv) == 3:
	file1 = sys.argv[1]
	file2 = sys.argv[2]
else:
	print "python2.7 pocket_matrix7.py <site1.pdb> <site2.pdb>" # site1.pdb site2.pdb
	sys.exit()
	
start = time.time()
start1 = time.time()
	
aline = open(file1, 'r').readlines()
bline = open(file2, 'r').readlines()	
	



def PairWise(res, coord):
	arr = []
	for i in range(len(res)):
		x_ca, y_ca, z_ca = coord[i][0]
		x_cb, y_cb, z_cb = coord[i][1]
		x_cn, y_cn, z_cn = coord[i][2]
		
		for j in range(len(res)):
			if i != j:
				x1_ca, y1_ca, z1_ca = coord[j][0]
				x1_cb, y1_cb, z1_cb = coord[j][1]
				x1_cn, y1_cn, z1_cn = coord[j][2]
				ans_ca = math.sqrt(pow(( x_ca - x1_ca ),2) + pow(( y_ca - y1_ca ),2) + pow(( z_ca - z1_ca ),2))
				ans_cb = math.sqrt(pow(( x_cb - x1_cb ),2) + pow(( y_cb - y1_cb ),2) + pow(( z_cb - z1_cb ),2))
				ans_cn = math.sqrt(pow(( x_cn - x1_cn ),2) + pow(( y_cn - y1_cn ),2) + pow(( z_cn - z1_cn ),2))
				#print res[i], res[j], ans_ca, ans_cb, ans_cn
				#time.sleep(1)
				arr.append( [res[i], res[j], ans_ca, ans_cb, ans_cn] )
	return arr			


def center_residues(arr, coord):
	coord -= np.mean(coord, axis=0)
	dic = {1:"   ", 2:"  ", 3:" ", 4:""}
	brr = []
	for i in range(len(arr)):
		#print arr[i]
		j = [ "%.3f"%j for j in coord[i]]
		val = ''.join([ dic[len(k.split(".")[0])]+k for k in j ])
		#print val
		brr.append(arr[i][:30]+val+arr[i][54:])
		#print '\n' # j[:30]+i+j[54:]
		#time.sleep(1)
	return brr

def file_process(arr):

	'''
	
	remove redundancy and h_bond
	store coordinate
	1) From whole coordinate -> to centroid
	2) CB -> Again Simple. If Glycine, then Ca
	3) CA
	
	'''
	whole_dic = {}
	h_dic = {"H":0}
	brr, het_arr, coord = [], [], []
	for i in arr:
		i = i.strip()
		if i[:4] == "ATOM":
			if i[74:78].strip() not in h_dic:
				var = i[13:16].strip()+" "+i[17:26].strip()
				#var = i[13:16].strip()+" "+i[17:26].strip()
				if var not in whole_dic:
					whole_dic[var] = 0
					brr.append(i)
					coord.append([i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])
					
		if i[:6] == "HETATM":
			brr.append(i)	
			coord.append([i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])
	
	coord = np.asarray(coord, dtype='float')
	brr = center_residues(brr, coord)				
	
	dic1, dic2 = defaultdict(list), defaultdict(list)
	
	
	
	for line in brr:
		if line[:4] == "ATOM":
			val = line[17:20]+'-'+line[21:22]+'-'+line[22:26].strip()
			dic1[val].append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
			dic2[val].append(line[13:16].strip())
			
		

	res, coord1 = [], []
	for i,j in dic2.items():
		#print i,' is i'
		coord = np.asarray(dic1[i], dtype='float')
		
		cn = np.mean(coord, axis=0)
		if i[:3] == "GLY":
			for j1 in range(len(j)):
				if j[j1] == 'CA':
					cb = coord[j1]
				if j[j1] == 'CA':
					ca = coord[j1]
		
		else:
			for j1 in range(len(j)):
				if j[j1] == 'CB':
					cb = coord[j1]
				if j[j1] == 'CA':
					ca = coord[j1]	
					#print j1, coord[j1]

		res.append(i)
		if len(coord1) == 0:
			coord1 = [[ca, cb, cn]]
		else:	
			coord1 = np.append(coord1, [[ca, cb, cn]], axis=0)

	
	arr = PairWise(res, coord1)
	
	arr = sorted(arr, key = lambda x:float(x[3]))
	
	#random.shuffle(arr)
	#dic = defaultdict(list)
	dic, dic1_pairs = defaultdict(list), defaultdict(list)
	arr1 = []
	
	for i in arr:
		#print(i)
		dic[i[0]+' '+i[1]].append(i[2])
		dic[i[0]+' '+i[1]].append(i[3])
		dic[i[0]+' '+i[1]].append(i[4])
		arr1.append(i[0]+' '+i[1])
		dic1_pairs[i[0]].append(i[1])

	return arr1, dic, dic1_pairs, brr, het_arr
	
	
	

	


def compare_matrices(mat1, mat2):

	if (mat2[0]-1.1) < mat1[0] < (mat2[0]+1.1):
		if (mat2[1]-1.1) < mat1[1] < (mat2[1]+1.1):
			if (mat2[2]-1.1) < mat1[2] < (mat2[2]+1.1):
				return True
			else:
				return False
		else:
			return False	
	else:
		return False


def CheckDistance(S1, S2, v1, v2):
	#S1, S2 = S1[1:], S2[1:]
	#print 
	#if res_dic1[v1] == res_dic2[v2]::
	if compare_matrices(res_dic1[v1], res_dic2[v2]):
		v1a = v1.split(' ')[1]
		v2a = v2.split(' ')[1]
		#print (S1, S2)
		#print (v1, v2,'---')
		for i in zip(S1, S2):
			p1, p2 = i
			p1 = p1.split(' ')[1]
			p2 = p2.split(' ')[1]
			#print (res_dic1)
			#print (S1,S2)
			#time.sleep(1)
			if p1 == v1a or p2 == v2a:
				return False
			if not compare_matrices(res_dic1[p1+' '+v1a], res_dic2[p2+' '+v2a]):
				return False
			#print (p1, p2, v1a, v2a)
		return True
		
	else:
		return False
	


def Recursion(res_pair1, sequence):
	
	if sequence == 'start':
		for i in res_arr2:
			#if res_dic1[res_pair1] == res_dic2[i]:
			#print(res_pair1)
			if compare_matrices(res_dic1[res_pair1], res_dic2[i]):
				#print (i, res_dic1[res_pair1], res_dic2[i] )
				return res_pair1, i, 'First'
		return None, None, 'First'
		
	else:
		#print('\nelse--', res_pair1, sequence)
		#print (res_pair1, sequence)
		res1a, res1b = res_pair1.split(' ')
		res2a, res2b = sequence.split(' ')
		#print (res_pair1,' <__> ', sequence)
		for i in res_pairs_dic1[res1b]:
			if i != res1a and i not in dic_single1:
				for j in res_pairs_dic2[res2b]:
					if j != res2a and j not in dic_single2:
						Check = CheckDistance(SequenceArrays1, SequenceArrays2, res1b+' '+i, res2b+' '+j)
						#Check = CheckDistance(SequenceArrays1, SequenceArrays2, [[res1b],[i]], [[res2b],[j]])
						if Check:
							if i+'\t'+j not in dic_pair_captch:
								#print i,res1b,'---',dic_single1
								seq1 = ' '.join(SequenceArrays1)+' '+res1b+' '+i
								seq2 = ' '.join(SequenceArrays2)+' '+res2b+' '+j
								#print (seq1,'----')
								if seq1+'\t'+seq2 not in dic_whole:	 
									#print (res1b+' '+i, res2b+' '+j,' ekse out')
									return res1b+' '+i, res2b+' '+j, 'Next'
								#print (i,j,'---')
								
		return None, None, 'Next'
		
		
	
def SortedArr():
	#print SequenceArrays1
	#print SequenceArrays2
	return True
	arr = [ i[0].split(' ')[1]+' '+i[1].split(' ')[1] for i in zip(SequenceArrays1, SequenceArrays2) ]
	arr = sorted(arr)
	print len(arr),arr
	if len(arr) < 3: # return False for sensitive cases. execute quick break statement if same seq is obtained twice. so dont use this
		return True
	arr = '_'.join(arr)
	#print arr
	#time.sleep(11)
	if arr not in SortedArrDic:
		SortedArrDic[arr] = 0
		return True
	return False	



def PairNext(S1, S2):
	dic1, dic2 = {}, {}
	if len(S1) <= 1:
		return None, None, None, None
	for i in zip(S1[:-1], S2[:-1]):
		
		dic1[i[0].split(' ')[0]] = 0
		dic1[i[0].split(' ')[1]] = 0
		dic2[i[1].split(' ')[0]] = 0
		dic2[i[1].split(' ')[1]] = 0

	S1, S2 = S1[:-1], S2[:-1]
	p1 = S1[-1].split(' ')[1]
	p2 = S2[-1].split(' ')[1]
	for i in res_pairs_dic1[p1]:
		#if i+' '+p1 not in dic1 and p1+' '+i not in dic1:
		if i not in dic1:
			for j in res_pairs_dic2[p2]:
				if j not in dic2:
				#if j+' '+p2 not in dic2 and p2+' '+j not in dic2:  # res1b+' '+i, res2b+' '+j
					Check = CheckDistance(S1, S2, p1+' '+i, p2+' '+j)
					#Check = CheckDistance(SequenceArrays1, SequenceArrays2, p1+' '+i, p2+' '+j)
					if Check:
						seq1 = ' '.join(S1)+' '+p1+' '+i
						seq2 = ' '.join(S2)+' '+p2+' '+j	
						#print(seq1, seq2,'--')
						if seq1+'\t'+seq2 not in dic_whole:
							return(S1, S2, p1+' '+i, p2+' '+j)
	#print('Still')						
	return PairNext(S1, S2)			
				
				
		



def run():
	
	global dic_single1, dic_single2, SequenceArrays1, SequenceArrays2, dic_pair_captch, dic_whole

	dic_loop1, dic_loop2 = {}, {}
	Final1, Final2 = [], []
	BreakLoop = False
	#print(res_pairs_dic1)

	dic_whole, SortedArrDic = {}, {}
	for i in res_arr1:
		if BreakLoop:
			break
		#print '\n\n\n\n\n',i,'+++++++'
		#if i.split(' ')[0] != 'ARG-A-97': # 'ARG-A-97 TRP-A-41':
		#if i == 'ARG-A-97 TRP-A-41':
				#pass
		#print (i,' starting i --')
		ans = 'start'
		SequenceArrays1, SequenceArrays2 = [], []
		dic_loop1, dic_loop2 = {}, {}
		dic_single1, dic_single2 = {}, {}
		
		# add two nested while loops
		InitiateFirstBreak = False
		ObtainedCount = 0
		#SortedArrDic = {}
		#print 'NExt Loop\n\n'
		while True:
			dic_pair_captch = {}
			#print 'aa'
			#print ans,'--\n'
			#time.sleep(1)
			#if not ans:
			#	time.sleep(11)
			#	break
			#print InitiateFirstBreak
			if ans != 'start':
				if not SequenceArrays1:
					break
			if InitiateFirstBreak:
				#break
				#start1 = time.time()
				#print 'CAME'
				#print(SequenceArrays1,'hh', ObtainedCount)
				#time.sleep(1)
				if not SequenceArrays1:
					break
				if ObtainedCount <= 1:
					break
				#print('----------\n\n')
				#time.sleep(11)
				
				if len(SequenceArrays1) <= 2:
					seq1 = ' '.join(SequenceArrays1)
					seq2 = ' '.join(SequenceArrays2)
					dic_whole[seq1+'\t'+seq2] = 0
					break
				
				'''
				if len(SequenceArrays1) <= 1:
					if SequenceArrays1:
						
						seq1 = ' '.join(SequenceArrays1)
						seq2 = ' '.join(SequenceArrays2)
						#print (SequenceArrays1, SequenceArrays2)
						dic_whole[seq1+'\t'+seq2] = 0
						for j in zip(SequenceArrays1, SequenceArrays1):
							#print j
							#time.sleep(1)
							dic_whole[j[0]+'\t'+j[1]] = 0
						
					break
				'''	
				ObtainedCount = 0
				#print 'Came End'
				#print(SequenceArrays1, SequenceArrays2, len(SequenceArrays1))
				SequenceArrays1, SequenceArrays2, i1, ans = PairNext(copy.deepcopy(SequenceArrays1), copy.deepcopy(SequenceArrays2))	
				if not i1:
					break
				SequenceArrays1.append(i1)
				SequenceArrays2.append(ans)
				#print (SequenceArrays1, SequenceArrays2,'\n')
				#time.sleep(11)
				#print SortedArr()
				#print '_________'
				for j in zip(SequenceArrays1, SequenceArrays2):
					dic_pair_captch[j[0].split(' ')[0]+'\t'+j[1].split(' ')[0]] = 0

				i = i1
				
				#end1 = time.time()
				#print end1-start1
			
			while True:
				
				if not ans:
					break
				i, ans, CheckPoint = Recursion(i, ans)
				#print (i, ans, 'Obtained Pairs')
				ObtainedCount += 1
				if not ans:
					#print 'break called'
					break
				if i+'\t'+ans in dic_whole:
					break	
					
				InitiateFirstBreak = True	
				dic_single1[i.split(' ')[0]] = 0	
				dic_single2[ans.split(' ')[0]] = 0
				#dic_pair_captch[i.split(' ')[1]+'\t'+ans.split(' ')[1]] = 0
				dic_pair_captch[i.split(' ')[0]+'\t'+ans.split(' ')[0]] = 0
				#print (dic_pair_captch)
				SequenceArrays1.append(i)	
				SequenceArrays2.append(ans)	

			if not SortedArr():
				BreakLoop = True
				#BreakLoop = False
				#print SequenceArrays2
				#print SequenceArrays1
				#print 'False Detected'
				#time.sleep(11)
			seq1 = ' '.join(SequenceArrays1)
			seq2 = ' '.join(SequenceArrays2)
			#print (SequenceArrays1, SequenceArrays2)
			#print '\n'
			#time.sleep(.1)
			#print 'Enter dic', seq1, seq2
			dic_whole[seq1+'\t'+seq2] = 0
			#print dic_whole
			Final1.append(SequenceArrays1)
			Final2.append(SequenceArrays2)
	return Final1, Final2
	
	
	
	
def process_hits(Final1, Final2):
	arr = []
	for i in zip(Final1, Final2):
		#print (i, len(i[0]))
		arr.append([i[0], i[1], len(i[0])])

	arr = sorted(arr, key = lambda x:int(x[2]), reverse=True)
	NewArr = []
	#print len(arr)
	#print arr[14:16]
	for i in arr:
		#print i
		if int(i[2]) > 3:
			#print i
			#time.sleep(11)
			val1 = [ j.split(' ')[1] for j in i[0][:-1] ]
			val2 = [ j.split(' ')[1] for j in i[1][:-1] ]
			NewArr.append([val1, val2, len(val1)])
	'''
	if len(NewArr) < 10:
		for i in arr:
			if int(i[2]) == 3:
				val1 = [ j.split(' ')[1] for j in i[0][:-1] ]
				val2 = [ j.split(' ')[1] for j in i[1][:-1] ]
				if res_dic1[' '.join(val1)][2] < 4:
					#print res_dic1[' '.join(val1)][2],'0000', val1, val2
					NewArr.append([val1, val2, len(val1)])	
	'''

	#print len(NewArr)
	NewArray = []
	if not NewArr:
		return None
	NewArray.append(NewArr[0]) # base condition
	for i in NewArr:
		check = True
		for j in NewArray:
			dic = { j[0][k]+' '+j[1][k]:0 for k in range(len(j[0])) }
			#print j
			#print [ j[0][k]+' '+j[1][k] for k in range(len(j[0])) ]
			#time.sleep(11)
			dic_count = sum([ 1 for k in range(len(i[0])) if i[0][k]+' '+i[1][k] in dic  ])
			if dic_count == len(j[0]):
				check = False
				break
		if check:
			NewArray.append(i)
	#print len(NewArray)		
	return NewArray
	
	


	
# END OF MAIN MAPP CODE
# START OF ALIGNMENT CODE




def rmsd(V, W):
    D = len(V[0])
    N = len(V)
    result = 0.0
    for v, w in zip(V, W):
        result += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(result/N)


def kabsch_rmsd(P, Q, translate=False):
    P = kabsch_rotate(P, Q)
    return rmsd(P, Q)


def kabsch_rotate(P, Q):
    U = kabsch(P, Q)
    P = np.dot(P, U)
    return P	


def kabsch(P, Q, check):
    if check:
        P -= np.mean(P, axis=0)
        Q -= np.mean(Q, axis=0)
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        #print "reflection detected"
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return U


def centroid(X):
    C = X.mean(axis=0)
    return C
   
   

def blosum():
	global residue_pairs_dictionary
	residue_dict_single = {'G':'GLY','A':'ALA','V':'VAL','L':'LEU','I':'ILE','T':'THR','S':'SER','Y':'TYR','W':'TRP',\
	  'P':'PRO','F':'PHE','N':'ASN','Q':'GLN','D':'ASP','E':'GLU','R':'ARG','K':'LYS','H':'HIS',\
	   'C':'CYS','M':'MET'
	  }
	aline = []
	aline.append("A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *")
	aline.append("A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4")
	aline.append("R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4")
	aline.append("N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4")
	aline.append("D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4")
	aline.append("C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4")
	aline.append("Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4")
	aline.append("E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4")
	aline.append("G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4")
	aline.append("H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4")
	aline.append("I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4")
	aline.append("L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4")
	aline.append("K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4")
	aline.append("M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4")
	aline.append("F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4")
	aline.append("P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4")
	aline.append("S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4")
	aline.append("T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4")
	aline.append("W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4")
	aline.append("Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4")
	aline.append("V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4")
	res_info = aline[0].strip()
	res_info = re.sub(" {1,}"," ", res_info)
	res_info_split = res_info.split(" ")[:-4]
	ans = []
	for line in aline[1:]:
		line = line.strip()
		line = re.sub(" {1,}"," ",line)
		residue = line.split(" ")[0]
		min_max = []
		for i in line.split(" ")[1:-4]:
			min_max.append(int(i))
		minimum = min(min_max)
		if minimum < 0:
			new_min = abs(minimum)
			new_min_max = []
			for i in min_max:
				new_min_max.append(i+new_min)
		maxima = float(max(new_min_max))

		for j,k in zip(res_info_split, new_min_max):
			ans.append(residue_dict_single[residue]+" "+residue_dict_single[j]+" "+str(k))

	dic_temp = {}
	residue_pairs_dictionary = {}
	for i in ans:
		i0 = i.split(" ")[0]
		i1 = i.split(" ")[1]
		val = int(i.split(" ")[2])
		if i0 not in dic_temp:
			dic_temp[i0] = 0
			residue_pairs_dictionary[i0] = {}
		residue_pairs_dictionary[i0][i1] = val
blosum()
   
   
   
def dihedral1(aa1, aa2):	
	arr = ["_CA","_CN","_N"]
	arr3 = []
	# [('_CA', '_CN', '_N'), ('_CA', '_N', '_CN'), ('_CN', '_CA', '_N'), ('_CN', '_N', '_CA'), ('_N', '_CA', '_CN'), ('_N', '_CN', '_CA')]
	arr3.append(["_CA", "_CN", "_N"])
	arr1, arr2 = [], []
	ans = []
	for j in arr3:
		arr1, arr2 = [], []
		for i, k in zip(arr, j):
			if aa1+i in arr1_dihedral_dic and aa2+k in arr2_dihedral_dic:
				arr1.append(arr1_dihedral_dic[aa1+i])
				arr2.append(arr2_dihedral_dic[aa2+k])
		if len(arr1) == 3:
			try:
				arr1 = np.asarray(arr1, dtype=float)
				arr2 = np.asarray(arr2, dtype=float)
				ap1, ap2, ap3 = arr1[0], arr1[1], arr1[2]
				bp1, bp2, bp3 = arr2[0], arr2[1], arr2[2]
				av1 = ap3 - ap1
				av2 = ap2 - ap1
				acp = np.cross(av1, av2)
				bv1 = bp3 - bp1
				bv2 = bp2 - bp1
				bcp = np.cross(bv1, bv2)
				a1, b1, c1, a2, b2, c2 = acp[0], acp[1], acp[2], bcp[0], bcp[1], bcp[2]
				d = ( a1 * a2 + b1 * b2 + c1 * c2 ) 
				e1 = math.sqrt( a1 * a1 + b1 * b1 + c1 * c1) 
				e2 = math.sqrt( a2 * a2 + b2 * b2 + c2 * c2) 
				d = d / (e1 * e2) 
				A = math.degrees(math.acos(d))
				ans.append(A)
			except:
				ans.append(10.0)
		else:
			ans.append(10.0)

	for i in range(0, 360, 10):
		if i <= ans[0] <= i+10:
			return 360.0-i
	return 360.0	 
    
    
def SiteGen():
	arr = copy.deepcopy(B_all) 
	minim = []   
	#print len(site2_coord)
	ans_dic = defaultdict(list)
	for i in range(0, len(arr)):
		info1, name1 = pdb1_res_info[i]
		if name1 == "CA":
			#p#rint name1
			x, y, z = arr[i]
			for j in range(0, len(site2_coord)):
				info2, name2 = pdb2_res_info[j]
				if name2 == "CA":
					x1, y1, z1 = site2_coord[j]
					x_ans = pow((x-x1),2)
					y_ans = pow((y-y1),2)
					z_ans = pow((z-z1),2)
					ans = math.sqrt(x_ans + y_ans + z_ans)
					minim.append(ans)
					if ans < 1.25:
						#print info1, info2, ans
						ans_dic[info1].append(info2)
	
	
	global site_check_dic_len, site_check_res_corr, site_check_sum1, arr1_dihedral_dic, arr2_dihedral_dic
	arr1_dihedral_dic, arr2_dihedral_dic = {}, {}
	
		
	dic_cent, dic_ca, dic_n = defaultdict(list), defaultdict(list), defaultdict(list)			
	for i in range(0, len(arr)):
		info, name = pdb1_res_info[i]
		dic_cent[info].append(arr[i])
		if name == "CA":
			dic_ca[info].append(arr[i])
		if name == "N":	
			dic_n[info].append(arr[i])
	for i in dic_cent.items():
		arr1_dihedral_dic[i[0]+"_CN"] = np.mean(i[1], axis=0)
		arr1_dihedral_dic[i[0]+"_CA"] = np.asarray(dic_ca[i[0]])[0]
		arr1_dihedral_dic[i[0]+"_N"] = np.asarray(dic_n[i[0]])[0]		
			
	
	dic_cent, dic_ca, dic_n = defaultdict(list), defaultdict(list), defaultdict(list)
	for i in range(0, len(site2_coord)):
		info, name = pdb2_res_info[i]
		dic_cent[info].append(site2_coord[i])
		if name == "CA":
			dic_ca[info].append(site2_coord[i])
		if name == "N":	
			dic_n[info].append(site2_coord[i])
	for i in dic_cent.items():
		arr2_dihedral_dic[i[0]+"_CN"] = np.mean(i[1], axis=0)
		arr2_dihedral_dic[i[0]+"_CA"] = np.asarray(dic_ca[i[0]])[0]
		arr2_dihedral_dic[i[0]+"_N"] = np.asarray(dic_n[i[0]])[0]		
			
		
	arr, dihed_factor, site_check_sum1 = [], [], []	
	for i in ans_dic.items():
		dihed_val = dihedral1(i[0], Counter(i[1]).most_common(1)[0][0])
		dihed_factor.append(dihed_val)
		#print i, dihed_val,Counter(i[1])
		if dihed_val > 245:
			arr.append([i[0]+" "+Counter(i[1]).most_common(1)[0][0]])
			#print arr[-1]
		site_check_sum1.append(residue_pairs_dictionary[i[0][:3]][Counter(i[1]).most_common(1)[0][0][:3]])

	site_check_sum1 = sum(site_check_sum1)*sum(dihed_factor)		
	return site_check_sum1, arr
    
    
def SiteGen1(arr):
	
	#print len(arr), len(site1a), len(site2a)
	arr1 = []
	dic = {1:"   ", 2:"  ", 3:" ", 4:""}
	for i in arr:
		var = ""
		for j in i:
			j1 = "%.3f"%j
			var += dic[len(j1.split(".")[0])]+j1
		#print var	
		#print i
		#print var
		#time.sleep(1)
		arr1.append(var)
	#for i in zip(arr,arr1):
	#	print i
	
	out = open("frag.pdb", 'w')
	#print len(arr1), len(site1a),' ---'
	for i1 in range(len(arr1)):
		i = arr1[i1]
		j = site1a[i1]
		out.write(j[:30]+i+j[54:]+"\n")		
	out.close()	
	
	out = open('fixed.pdb', 'w')
	for i in site2a:
		out.write(i+"\n")
	out.close()	
	
	 
def site_gen_het(site_gen_het):
	dic1, dic2 = {}, {}
	out = open("align.txt", 'w')
	#print site_gen_het,'-----'
	for i in site_gen_het:
		#print i
		i1 = i[0].split(" ")
		#print i
		out.write(i1[1]+" "+i1[0]+"\n")
		dic1[i1[0]] = 0
		dic2[i1[1]] = 0
	out.close()	
	arr = B_all
	arr1 = []
	dic = {1:"   ", 2:"  ", 3:" ", 4:""}
	for i in arr:
		var = ""
		for j in i:
			j1 = "%.3f"%j
			var += dic[len(j1.split(".")[0])]+j1
		#print var	
		arr1.append(var)


	out = open("site1.pdb", 'w')
	for i1 in range(0, len(arr1)):
		i = arr1[i1]
		j = site1a[i1]
		#print i,j
		if j[:4] == "ATOM":
			res_info = j[17:20]+"-"+j[21:22]+"-"+j[22:26].strip()
			#print res_info
			if res_info in dic1:
				#print j
				out.write(j[:30]+i+j[54:]+"\n")
		else:	
			out.write(j[:30]+i+j[54:]+"\n")
	out.close()	
	#print site1a
	out = open("site2.pdb", 'w')
	for i in site2a:
		if i[:4] == "ATOM":
			res_info = i[17:20]+"-"+i[21:22]+"-"+i[22:26].strip()
			#print i, res_info
			if res_info in dic2:
				#print i, res_info
				out.write(i+"\n")
		else:	
			out.write(i+"\n")
	out.close()	
 
    
def print_scores(arr):
	res1, res2, summer = [], [], []
	control1, control2 = [], []
	min_max, min_max2 = 0, []
	fin_val = ""
	for i in arr:
		#print i,' i1'
		i1, i2 = i[0].split(' ')
		#print i
		res1.append(i2)
		res2.append(i1)
		i1, i2 = i1[:3], i2[:3]
		#print residue_pairs_dictionary[i1][i2],"ass",i1,i2
		summer.append(residue_pairs_dictionary[i1][i2])
		control1.append(residue_pairs_dictionary[i1][i1])
		control2.append(residue_pairs_dictionary[i2][i2])
	res1 = sorted(set(res1))
	res2 = sorted(set(res2))
	min_max = len(arr)
	min_max2.append(pdb2_ln)
	min_max2.append(pdb1_ln)
	summer = float(sum(summer))
	#print summer, control1, control2,"asaa"
	if len(control1) == 0 or len(control2) == 0:
		return 'No Score'
	blosum_score = summer/((sum(control1)/2)+(sum(control2)/2))
	#print sum(control2)/2	
	#print np.mean(control)	
	#print control1, control2
	#print str(len(res1))+"/"+str(pdb1_ln)
	fin_val = str(len(res1))+"/"+str(pdb1_ln)+" "+str(len(res2))+"/"+str(pdb2_ln)+" "+str(float(min_max)/min(min_max2))+" "+str(float(min_max)/max(min_max2))+\
	" "+str(blosum_score)
	return fin_val
	
	

def MainCode(aline, bline):
	global res_dic1, res_dic2, res_arr1, res_arr2, res_pairs_dic1, res_pair2_dic2
	global dic_loop1, dic_loop2, dic_whole, SortedArrDic
	global dic_single1, dic_single2, site1a, site2a
	global pdb1_ln, pdb2_ln

	global res_arr1, res_arr2, res_dic1, res_dic2, res_pairs_dic1, res_pairs_dic2, B_all
	global pdb1_res_info, pdb2_res_info, site1_coord, site2_coord
	
	res_arr1, res_dic1, res_pairs_dic1, pdb1_lines, pdb1_het_lines = file_process(aline)
	res_arr2, res_dic2, res_pairs_dic2, pdb2_lines, pdb2_het_lines = file_process(bline)	
	Final1, Final2 = run()	
	NewArray = process_hits(Final1, Final2)
	if not NewArray:
		print 'None\tNone'
		sys.exit()
	#for i in NewArray[:3]:
	#	print i
	#print pdb1_lines[:3]
	#print pdb1_het_lines[:3]
	#time.sleep(11)
	# Note bring coordinate to center befroe doing any opertation
	
	pdb1_trans_coord, pdb1_res_info, pdb1_generated_coord, pdb1_ca_dic = [], [], [], defaultdict(list)
	pdb2_trans_coord, pdb2_res_info, pdb2_generated_coord, pdb2_ca_dic = [], [], [], defaultdict(list)
	pdb1_ln, pdb2_ln = 0, 0
	for line in pdb1_lines:
		if line[:4] == "ATOM":
			res1 = line[17:20]+"-"+line[21:22]+"-"+line[22:26].strip()
			pdb1_res_info.append([ res1, line[13:16].strip() ])
			pdb1_generated_coord.append(line)
			pdb1_trans_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
			if line[13:16].strip() == "CA":
				pdb1_ln += 1
				pdb1_ca_dic[res1].append(line[28:38].strip())
				pdb1_ca_dic[res1].append(line[38:46].strip())
				pdb1_ca_dic[res1].append(line[46:54].strip())
		
	for line in pdb2_lines:
		if line[:4] == "ATOM":
			res1 = line[17:20]+"-"+line[21:22]+"-"+line[22:26].strip()
			pdb2_res_info.append([ line[17:20]+"-"+line[21:22]+"-"+line[22:26].strip(), line[13:16].strip() ])
			pdb2_generated_coord.append(line)
			pdb2_trans_coord.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])	
			if line[13:16].strip() == "CA":
				pdb2_ln += 1
				pdb2_ca_dic[res1].append(line[28:38].strip())
				pdb2_ca_dic[res1].append(line[38:46].strip())
				pdb2_ca_dic[res1].append(line[46:54].strip())
		
		
	pdb1_trans_coord = np.asarray(pdb1_trans_coord, dtype='float')	
	pdb2_trans_coord = np.asarray(pdb2_trans_coord, dtype='float')
		
	ResLists = []
	maxi, index_ln = 0, 0
	#for i in [NewArray[7]]:
	#print NewArray
	NewCount = 0
	for i in NewArray:
		#print i,'--'
		site1_arr, site2_arr = [], []
		site1_coord, site2_coord = copy.deepcopy(pdb1_trans_coord), copy.deepcopy(pdb2_trans_coord)
		
		#print i
		for j in range(len(i[0])):
			#print i[0][j], pdb1_ca_dic[i[0][j]]
			#print i[0][j], i[1][j]
			site1_arr.append(pdb1_ca_dic[i[0][j]])
			site2_arr.append(pdb2_ca_dic[i[1][j]])
		site1_arr, site2_arr = np.asarray(site1_arr, dtype=float), np.asarray(site2_arr, dtype=float)	
		#print site1_arr	
		#print pdb1_trans_coord
		
		U = kabsch(copy.deepcopy(site1_arr), copy.deepcopy(site2_arr), True)	
		B_all = copy.deepcopy(site1_coord)
		B_all -= site1_arr.mean(axis=0)
		B_all = np.dot(B_all, U)
		B_all += site2_arr.mean(axis=0)
		
		site1a = pdb1_generated_coord
		site2a = pdb2_generated_coord
		
		score, new_res_list = SiteGen()
		#print i,score,'---\n'
		if score > maxi:
			maxi = score
			index_ln = len(i[0])
		else:
			if len(i[0]) < index_ln:
				NewCount += 1
				
		ResLists.append([new_res_list, score])
		
	ResLists = sorted(ResLists, key=lambda x:float(x[1]), reverse=True)	
	#print ResLists[0]
	line1 = ResLists[0][0]
	#print ResLists[0][0]
	site1_arr, site2_arr = [], []
	#print ResLists[0:4]
	MAPP_scores = print_scores(line1)
	MAPP_seqs = '_'.join([i[0] for i in line1])
	#print MAPP_seqs,[i for i in line1]
	if not MAPP_seqs:
		MAPP_seqs = "No Match"
	#if not MAPP_scores:
			
	print MAPP_seqs, MAPP_scores,' --'
	for i in line1:
		#print i,i[0]
		i1 = i[0].split(' ')
		site1_arr.append(pdb1_ca_dic[i1[0]])
		site2_arr.append(pdb2_ca_dic[i1[1]])
	site1_arr, site2_arr = np.asarray(site1_arr, dtype='float'), np.asarray(site2_arr, dtype='float')
	
	
	# Rotation and alignment
	
	
	site1, site2 = [], []
	for line in pdb1_lines:
		if line[:4] == "ATOM" or line[:6] == "HETATM":
			site1.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
	for line in pdb2_lines:
		if line[:4] == "ATOM" or line[:6] == "HETATM":
			site2.append([line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])		
	site1 = np.asarray(site1, dtype='float')
	site2 = np.asarray(site2, dtype='float')
	
	site1a, site2a = pdb1_lines, pdb2_lines
	#print site1_arr, site2_arr
	if len(site1_arr) == 0 or len(site2_arr) == 0:
		print 'No result'
		sys.exit()
	site1_arr_cnt = site1_arr.mean(axis=0)
	site2_arr_cnt = site2_arr.mean(axis=0)
	site2_new = copy.deepcopy(site2)
	site2 -= site2_arr_cnt
	
	
	U = kabsch(site1_arr, site2_arr, False)	
	B_all = copy.deepcopy(site1)
	B_all -= site1_arr.mean(axis=0)
	B_all = np.dot(B_all, U)
	B_all += site2_arr.mean(axis=0)
	
	SiteGen1(B_all)
	site_gen_het(line1)


MainCode(aline, bline)



#print count










