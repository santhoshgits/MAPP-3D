import os
import time
import sys
from collections import Counter, defaultdict
import copy


if len(sys.argv) == 4:
	files = sys.argv[1]
	pdb = sys.argv[2]
	nos = int(sys.argv[3])
else:
	print "Motif.py <align-output> <PDB-Rep> <No. of Residue>"
	sys.exit()

# 6SRB_A_GSH_A_300.pdb
aline = open(files, 'r').readlines()
arr = []
for line in aline:
	line = line.strip()
	l = line.split('\t')
	if len(l[2].split(' ')) == 5: 
		if l[0] == pdb:
			score = l[2].split(' ')[0].split('/')[0]
			if int(score) >= nos:
				arr.append(l[3])
				#print '\n'

res_dict = {'GLY':'G' ,'ALA':'A' ,'LEU':'L' ,'ILE':'I' ,'TRP':'W' ,'SER':'S' ,'THR':'T' ,'TYR':'W', 'PHE':'F' ,'PRO':'P' ,'ASP':'D' ,\
		   'GLU':'E' ,'HIS':'H' ,'CYS':'C' ,'MET':'M' ,'VAL':'V' ,'ASN':'N' ,'GLN':'Q' ,'ARG':'R' ,'LYS':'K' }


def AvgDistance(arr1, arr2):
	count_arr = []
	#print arr2
	#print len(arr2)
	#print len(arr1)
	#print '\n'
	#print arr1[0:6], arr2[0:6]
	for i in zip(arr1, arr2):
		if i[0] != '-' and i[1] != '-': 
			#print i
			count_arr.append(abs(int(i[0].split('-')[-1]) - (int(i[1].split('-')[-1])-1)))
	ln = len(count_arr)
	#print count_arr
	if not count_arr:
		return None
	arr = []
	maxi = Counter(count_arr).most_common(1)[0][1]
	#arr.append(('THR-A-101', 'ARG-A-104'))
	#print maxi,' maxis'
	for i, j in Counter(count_arr).items():
		arr.append([i, j])
		
	'''	
	arr = []	
	arr.append([3,93])	
	arr.append([2,93])	
	
	arr.append([4,93])	
	arr.append([5,93])	
	arr.append([6,93])	
	arr.append([9,93])	
	arr.append([11,93])	
	arr.append([12,93])	
	'''	
	gap_arr = []	
	arr = sorted(arr, key = lambda x:int(x[1]), reverse=True)
	for i in arr:
		if float(i[1])/maxi >= .3:
			#print i,' --', i[1]
			gap_arr.append(i[0])
	gap_arr = map(int, gap_arr)
	gap_arr = sorted(gap_arr)
	#print gap_arr
	
	gap_val = None
	
	if gap_arr[0] == 0:
		#gap_val = "NO"
		return None
	if 3 > gap_arr[0] > 0:
		gap_val = ""
		#print 'aas'
		maxi = max([i for i in gap_arr if i < 3])
		for _ in range(0, maxi):
			#print 'hi'
			gap_val += 'x'
			gap_val += '-'
		#print gap_val,'---'
		return gap_val[:-1]	
		
	else:
		
		arr1, arr2 = [], []
		#gap_arr.append(gap_arr[-1])
		#if len(gap_arr)> 1:
		#	print 'got it', gap_arr
		#	time.sleep(11)
		count = 0
		
		if len(gap_arr) == 1:
			arr1.append(str(gap_arr[0]))
			#break
		for i in range(len(gap_arr)-1):
			count += 1
			#print gap_arr[i] - gap_arr[i+1], gap_arr[i],']]]', count, len(gap_arr)
			
			if abs(gap_arr[i] - gap_arr[i+1]) <= 1:
				#print gap_arr[i],'++'
				arr2.append(gap_arr[i])
			else:
				#print 'hi'
				if arr2:
					arr2.append(gap_arr[i])
					arr1.append(str(arr2[0])+','+str(arr2[-1]))
					#print str(arr2[0])+'-'+str(arr2[-1]),' --', arr2
					arr2 = []
				else:
					arr1.append(str(gap_arr[i]))
			if count == len(gap_arr)-1:
				if arr2:
					arr1.append(str(arr2[0])+','+str(gap_arr[i+1]))
				else:
					#print gap_arr[i+1], arr1
					arr1.append(str(gap_arr[i+1]))
				break	
	#print arr2
	#print arr1,'arr1'
	#time.sleep(1)
	#gap_arr = map(str, gap_arr)
	#if len(arr1) > 1:
	#	print arr1
	arr1 = arr1[0]
	
	#gap_val = 'x'+'('+''.join(arr1)+')'		
	gap_val = 'x'+'('+arr1+')'		
	#print gap_arr
	return gap_val






def SeqWeights(arr):
	arr1 = []
	for i in arr:
		if i[0] != '-':
			arr1.append(i[:3])
	#print arr1
	maxi = Counter(arr1).most_common(1)[0][1]
	arr2 = []
	count1, count2 = 0, 0
	for i in Counter(arr1).items():
		#print i
		count1 += int(i[1])
		if float(i[1])/maxi > .5:
			count2 += int(i[1])
			arr2.append(res_dict[i[0][:3]])
	#print arr2, str(count2)+'/'+str(count1)
	arr2 = "".join(arr2)
	if len(arr2) > 1:
		arr2 = '['+arr2+']'
	return arr2, str(count2)+'/'+str(count1)





def Motifs(seqs):
	arr = []
	#print seqs
	for i in seqs:		
		for j in i.split('_'):
			arr.append(j.split(' ')[0])
			#print j
			#print j.split(' ')[0]
			
	arr = sorted(set(arr))
	arr = sorted(arr, key = lambda x:int(x.split('-')[-1]))
	#print arr
	dic1 = {}
	for j in arr:
		dic1[j] = '-'
	#print dic1
	#motifs = []
	residue_conserv = defaultdict(list)
	for i in seqs:
		#motif = []
		dic = copy.deepcopy(dic1)
		for j in i.split('_'):
			dic[j.split(' ')[0]] = j.split(' ')[1]
		#print residue_conserv,'\n\n';time.sleep(1)	
		for j in arr:
			#print j, dic[j]
			#motif.append([j, dic[j]])
			residue_conserv[j].append(dic[j])
		#motifs.append(motif)	
	MotifRange = []
	for i in range(0, len(arr)-1):
		#print arr[i],' he'
		#print residue_conserv[arr[i]]
		#print residue_conserv[arr[i+1]]
		#print arr[i],residue_conserv[arr[i]][:10]
		ResPattern, ResFreq = SeqWeights(residue_conserv[arr[i]])
		print ResPattern, ResFreq
		#time.sleep(11)
		distance = AvgDistance(residue_conserv[arr[i]], residue_conserv[arr[i+1]])
		#print distance
		#print '\n'
		#time.sleep(1)
		if distance:
			#print ResPattern+'-'+distance
			MotifRange.append(ResPattern)
			MotifRange.append(distance)
		else:
			MotifRange.append(ResPattern)
		#print ResPattern,distance
		
		#print residue_conserv[arr[i]]
		#print MotifRange
		#print '\n\n'
		#time.sleep(1)
	return "-".join(MotifRange)		




motif_seq = Motifs(arr)
print motif_seq



























