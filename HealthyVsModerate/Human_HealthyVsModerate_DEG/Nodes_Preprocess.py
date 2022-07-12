import re

class Nodes:
	def Nodes_Process(infile1, outfile1, infile2, infile3, infile4, outfile2, outfile3, outfile4, outfile5, infile5, outfile6, outfile7):
		

		file1 = open(infile1,'r')
		line1 = file1.readline()
		file2 = open(outfile1,'w')
		while line1:
			line1 = line1.rstrip()
			lines = re.sub('\..*','',line1); #print(lines);
			file2.write(lines+'\n')
			line1 = file1.readline()

		file2 = open(infile2,'r')
		for i in range(6):
			line2 = file2.readline()
		gene_id, transcript_id = [],[];
		while line2:
			line2 = line2.rstrip()
			split_line2 = line2.split('\t'); #print(split_line2[8])
			split_line2_8 = (split_line2[8]).split(';'); #print(split_line2_8);
			split_line2_8[0] = re.sub('gene_id "','',split_line2_8[0])
			split_line2_8[2] = re.sub(' transcript_id "','',split_line2_8[2]); #print(split_line2_8[0], split_line2_8[2]);
			gene_id.append(split_line2_8[0].rstrip('"'))
			transcript_id.append(split_line2_8[2].rstrip('"')); 
			line2 = file2.readline()

		file3 = open(infile3,'r')
		line3 = file3.readline()
		entrezid, geneid, transcriptid = [],[],[];
		while line3:
			line3 = line3.rstrip()
			split_line3 = line3.split('\t')
			#split_line3[4] = re.sub('\..*','',split_line3[4]); #print(split_line3[4]);
			entrezid.append(split_line3[4])
			geneid.append(split_line3[1])
			transcriptid.append(split_line3[2])
			line3 = file3.readline()

		file4 = open(infile4,'r')
		for i in range(2):
			line4 = file4.readline()
		file5 = open(outfile2,'w')
		file6 = open(outfile3,'w')
		file5_1 = open(outfile4,'w')
		file6_1 = open(outfile5,'w')
		while line4:
			line4 = line4.rstrip()
			split_line4 = line4.split('\t')
			transcript = split_line4[0]
			logfc = float(split_line4[2])
			if (transcript in gene_id):
				indices1 = int(gene_id.index(transcript))
				gene = gene_id[indices1]
				#file6.write(transcript+'\n')
				if gene in geneid:
					indices2 = int(geneid.index(gene))
					file5.write(entrezid[indices2]+'\n')
					file6.write(gene+'\t'+entrezid[indices2]+'\n')
					if (logfc >= 1):
						file5_1.write(entrezid[indices2]+'\n')
					elif (logfc <= -1):
						file6_1.write(entrezid[indices2]+'\n')

			#elif (transcript not in transcript_id):
			#	file6.write(transcript+'\n')
			line4 = file4.readline()

		file7 = open(infile5,'r')
		for i in range(2):
			line7 = file7.readline()
		file8 = open(outfile6,'w')
		file9 = open(outfile7,'w')
		while line7:
			line7 = line7.rstrip()
			split_line7 = line7.split('\t')
			transcript = split_line7[0]
			if (transcript in gene_id):
				indices1 = int(gene_id.index(transcript))
				gene = gene_id[indices1]
				#file9.write(transcript+'\n')
				if gene in geneid:
					indices2 = int(geneid.index(gene))
					file8.write(entrezid[indices2]+'\n')
					file9.write(gene+'\t'+entrezid[indices2]+'\n')
			#elif (transcript not in transcript_id):
			#	file9.write(transcript+'\n')
			line7 = file7.readline()

		file1.close(); file2.close(); file3.close(); file4.close(); file5.close(); file5_1.close(); file6.close(); file6_1.close();
		file7.close(); file8.close(); file9.close();


