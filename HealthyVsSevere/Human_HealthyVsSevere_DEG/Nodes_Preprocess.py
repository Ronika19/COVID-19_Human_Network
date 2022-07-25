import re

class Nodes:
	def data_extract(self, infile, deletes):
		f = open(infile, 'r')
		lines = f.readlines()
		dict_array = {}
		split_l = ((lines[1]).rstrip()).split('\t')
		for i in range(len(split_l)):
			dict_array['arr_'+str(i)] = []
		if (deletes == 0):
			for l in lines:
				split_l = (l.rstrip()).split('\t')
				for x in range(len(split_l)):
					dict_array['arr_'+str(x)].append(split_l[x])
		elif (deletes > 0):
			for l in lines[deletes:]:
				split_l = (l.rstrip()).split('\t')
				for x in range(len(split_l)):
					dict_array['arr_'+str(x)].append(split_l[x])
		f.close(); #print(dict_array);
		return dict_array

	def Nodes_Process(self, infile1, outfile1, infile2, infile3, infile4, outfile2, outfile3, outfile4, outfile5, infile5, outfile6, outfile7):
		file1 = open(infile1,'r'); file2 = open(outfile1,'w');
		for line1 in file1:
			lines = re.sub('\..*','',line1.rstrip()); #print(lines);
			file2.write(lines+'\n')

		dict1 = self.data_extract(infile2,5)
		lines = dict1['arr_8']
		gene_id, transcript_id = [],[];
		for l in lines:
			split_l = l.split(';'); 
			split_l[0] = re.sub('gene_id "','',split_l[0])
			split_l[2] = re.sub(' transcript_id "','',split_l[2]); 
			gene_id.append((split_l[0]).rstrip('"'))
			transcript_id.append((split_l[2]).rstrip('"')); 

		dict2 = self.data_extract(infile3,0)
		entrezid, geneid, transcriptid = dict2['arr_4'],dict2['arr_1'],dict2['arr_2']
		
		dict3 = self.data_extract(infile4,1)
		transcript, logfc = dict3['arr_0'],dict3['arr_2']
		file5 = open(outfile2,'w'); file5_1 = open(outfile4,'w'); file6 = open(outfile3,'w'); file6_1 = open(outfile5,'w');
		for i in range(len(transcript)):
			if (transcript[i] in gene_id):
				indices1 = int(gene_id.index(transcript[i]))
				gene = gene_id[indices1]
				if gene in geneid:
					indices2 = int(geneid.index(gene))
					file5.write(entrezid[indices2]+'\n')
					file6.write(gene+'\t'+entrezid[indices2]+'\n')
					if (float(logfc) >= 1):
						file5_1.write(entrezid[indices2]+'\n')
					elif (float(logfc) <= -1):
						file6_1.write(entrezid[indices2]+'\n')
			#elif (transcript[i] not in transcript_id):
			#	file6.write(transcript[i]+'\n')
		
		dict4 = self.data_extract(infile5,1)
		transcript = dict4['arr_0']
		file8 = open(outfile6,'w'); file9 = open(outfile7,'w');
		for i in range(len(transcript)):
			if (transcript[i] in gene_id):
				indices1 = int(gene_id.index(transcript[i]))
				gene = gene_id[indices1]
				if gene in geneid:
					indices2 = int(geneid.index(gene))
					file8.write(entrezid[indices2]+'\n')
					file9.write(gene+'\t'+entrezid[indices2]+'\n')
			#elif (transcript[i] not in transcript_id):
			#	file9.write(transcript[i]+'\n')
			

		file1.close(); file2.close(); file5.close(); file5_1.close(); file6.close(); file6_1.close(); file8.close(); file9.close();


