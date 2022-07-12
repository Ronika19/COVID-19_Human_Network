import re
from .Nodes_Preprocess import Nodes

class Data_Processer:
	def DataPreprocess():
		infile_list = ['Input_Files/COVID19_Human_Nodes_All.txt', 'Input_Files/Homo_sapiens.GRCh38.101.gtf', 'Input_Files/Homo_sapiens.GRCh38.103.entrez.tsv', 'Input_Files/Moderate_Vs_Severe_Significant.txt', 'Input_Files/Moderate_Vs_Severe.txt']
		outfile_entrez_list = ['Output_Files/ModerateVsSevere_Significant_EntrezIDs.txt', 'Output_Files/ModerateVsSevere_Up_Significant_EntrezIDs.txt', 'Output_Files/ModerateVsSevere_Down_Significant_EntrezIDs.txt', 'Output_Files/ModerateVsSevere_EntrezIDs.txt']
		outfile_genehandles_list = ['Output_Files/COVID19_Human_NodesAll.txt', 'Output_Files/ModerateVsSevere_Significant_GeneHandles.txt', 'Output_Files/ModerateVsSevere_GeneHandles.txt']

		Nodes_Preprocess = Nodes.Nodes_Process(infile_list[0], outfile_genehandles_list[0], infile_list[1], infile_list[2], infile_list[3], outfile_entrez_list[0], outfile_genehandles_list[1], outfile_entrez_list[1], outfile_entrez_list[2], infile_list[4], outfile_entrez_list[3], outfile_genehandles_list[2])

	def DEG_Cluster(infile, outfile):
		file1 = open(infile,'r')
		file4 = open(outfile,'w')

		for i in range(2):
			line1 = file1.readline()
		while line1:
			line1 = line1.rstrip()
			split_line1 = line1.split('\t')
			genes1 = split_line1[0]
			fold_change1 = split_line1[2]
			if ((float(fold_change1) >= 1) and (float(split_line1[6]) <= 0.05)):
				file4.write(genes1+"\t")
			if ((float(fold_change1) <= -1) and (float(split_line1[6]) <= 0.05)):
				file4.write(genes1+"\t")
			line1 = file1.readline()
		file4.write("\n")
		file1.close(); file4.close();

	def Module_Clusters(infile, outfile1, outfile2):
		file1 = open(infile,'r')
		line1 = file1.readline()
		line1 = file1.readline()
		genes, modules = [],[];
		while line1:
			line1 = (line1.rstrip()).replace('"','')
			split_line1 = line1.split()
			genes.append(split_line1[0])
			modules.append(int(split_line1[1]))
			line1 = file1.readline()
		print(len(genes),len(modules))

		mods = []
		for module in modules:
			if module not in mods:
				mods.append(int(module))
		mods.sort()
		x = 0; gene=[];
		file2 = open(outfile1,'w')
		file3 = open(outfile2,'w')
		for mod in mods:
			print(mod); file2.write(str(mod)+"\t")
			for x in range(len(modules)):
				if (modules[x] == mod):
					print(mod, modules[x], genes[x])
					gene.append(genes[x])
					file2.write(genes[x]+"\t")
					file3.write(genes[x]+"\t")
					#file2.write(str(mod)+"\t"+str(modules[x])+"\t"+genes[x]+"\n")
			file2.write("\n")
			file3.write("\n")

		file1.close(); file2.close(); file3.close();

	def DEGClusters_2_WGCNAModules(infile1, infile2, outfile1, infile3, outfile2):
		# All DEGs
		file1 = open(infile1,'r')
		line1 = file1.readline()
		DEG_transcripts = []; indices = [];
		while line1:
			line1 = line1.rstrip()
			split_line1 = line1.split('\t')
			for tids in split_line1:
				DEG_transcripts.append(tids)
			DEG_transcripts.append('\n')
			line1 = file1.readline()

		for items in range(len(DEG_transcripts)):
			if DEG_transcripts[items] == '\n':
				indices.append(items); print(items)

		# Modules using WGCNA
		file2 = open(infile2,'r')
		for i in range(2):
			line2 = file2.readline()
		modules, gene_id = [],[];
		while line2:
			line2 = line2.rstrip()
			split_line2 = line2.split('\t')
			gene_id.append(split_line2[0].replace('"',''))
			modules.append(split_line2[1])
			line2 = file2.readline()
		print(DEG_transcripts,"\n",gene_id)

		file3 = open(outfile1,'w')
		k = 0; DEG1_Module = [];
		for genes in DEG_transcripts:
			if (k < int(indices[0])):	# Cluster of Pathogenesis DEGs
				if (genes in gene_id):
					indexes = gene_id.index(genes); print(indexes, genes, gene_id[indexes]);
					DEG1_Module.append(modules[indexes])
					file3.write(gene_id[indexes]+"\t"+modules[indexes]+"\n")
			k += 1

		#print(DEG1_Module,"\n\n", DEG2_Module,"\n\n", DEG3_Module)
		file3.close(); file1.close(); file2.close();

		file4 = open(infile3,'r')
		line4 = file4.readline()
		geneid, mod = [],[];
		while line4:
			line4 = line4.rstrip()
			split_line4 = line4.split('\t')
			geneid.append(split_line4[0])
			mod.append(int(split_line4[1]))
			line4 = file4.readline()

		modset = sorted(set(mod), reverse=False); print(modset);

		file5 = open(outfile2,'w')
		for i in modset:
			file5.write(str(i)+'\t')
			indices = [index for index, element in enumerate(mod) if element == i]; #print(i, indices);
			m = 0
			for j in indices:
				m += 1; print(geneid[j], mod[j])
				if m < len(indices):
					file5.write(geneid[j]+',')
				elif m == len(indices):
					file5.write(geneid[j]+'\n')

		file4.close(); file5.close();


