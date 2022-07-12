import re

file1 = open('COVID19_Human_WGCNA_Edges1.txt','r')
file2 = open('WGCNA_Reduced_Edges_0.15.txt','w')
line1 = file1.readline()
file2.write(line1)
line1 = file1.readline()
while line1:
	#line1 = line1.rstrip()
	weight = float(line1.split('\t')[2])
	#if (weight >= 0.2):
	if (weight >= 0.19):
		file2.write(line1)
	line1 = file1.readline()





















