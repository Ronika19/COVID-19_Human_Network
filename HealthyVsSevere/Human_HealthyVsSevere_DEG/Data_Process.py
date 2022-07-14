import re
from .Nodes_Preprocess import Nodes

class Data_Processer:
	def DataPreprocess():
		infile_list = ['Input_Files/COVID19_Human_Nodes_All.txt', 'Input_Files/Homo_sapiens.GRCh38.101.gtf', 'Input_Files/Homo_sapiens.GRCh38.103.entrez.tsv', 'Input_Files/Healthy_Vs_Severe_Significant.txt', 'Input_Files/Healthy_Vs_Severe.txt']
		outfile_entrez_list = ['Output_Files/HealthyVsSevere_Significant_EntrezIDs.txt', 'Output_Files/HealthyVsSevere_Up_Significant_EntrezIDs.txt', 'Output_Files/HealthyVsSevere_Down_Significant_EntrezIDs.txt', 'Output_Files/HealthyVsSevere_EntrezIDs.txt']
		outfile_genehandles_list = ['Output_Files/COVID19_Human_NodesAll.txt', 'Output_Files/HealthyVsSevere_Significant_GeneHandles.txt', 'Output_Files/HealthyVsSevere_GeneHandles.txt']

		Nodes_Preprocess = Nodes().Nodes_Process(infile_list[0], outfile_genehandles_list[0], infile_list[1], infile_list[2], infile_list[3], outfile_entrez_list[0], outfile_genehandles_list[1], outfile_entrez_list[1], outfile_entrez_list[2], infile_list[4], outfile_entrez_list[3], outfile_genehandles_list[2])

	def DEG_Cluster(infile, outfile):
		f = open(outfile,'w')
		dict1 = Nodes().data_extract(infile,1)
		gene_ids, fold_change, pval = dict1['arr_0'], dict1['arr_2'], dict1['arr_6']
		for i in range(len(fold_change)):
			if ((float(fold_change[i]) >= 1) and (float(pval[i]) <= 0.05)):
				f.write(gene_ids[i]+"\t")
			if ((float(fold_change[i]) <= -1) and (float(pval[i]) <= 0.05)):
				f.write(gene_ids[i]+"\t")
		f.write("\n"); f.close();

	def Module_Clusters(infile, outfile1, outfile2):
		dict2 = Nodes().data_extract(infile,1)
		genes_ids, modules = dict2['arr_0'], dict2['arr_1']
		genes = [ids.replace('"','') for i, ids in enumerate(genes_ids)]

		mods = []
		for module in modules:
			if int(module) not in mods:
				mods.append(int(module))
		mods.sort()
		x = 0; gene=[];
		file2 = open(outfile1,'w'); file3 = open(outfile2,'w');
		for mod in mods:
			file2.write(str(mod)+"\t")
			for x in range(len(modules)):
				if (int(modules[x]) == int(mod)):
					gene.append(genes[x])
					file2.write(genes[x]+"\t"); file3.write(genes[x]+"\t");
			file2.write("\n"); file3.write("\n");
		file2.close(); file3.close();

	def DEGClusters_2_WGCNAModules(infile1, infile2, outfile1, infile3, outfile2):
		# All DEGs
		file1 = open(infile1,'r'); file3 = open(outfile1,'w'); file5 = open(outfile2,'w');
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
		dict3 = Nodes().data_extract(infile2,1) 
		gene_ids, modules = dict3['arr_0'], dict3['arr_1']
		gene_id = [val.replace('"','') for i, val in enumerate(gene_ids)]
		
		k = 0; DEG1_Module = [];
		for genes in DEG_transcripts:
			if (k < int(indices[0])):	# Cluster of DEGs
				if (genes in gene_id):
					indexes = gene_id.index(genes); print(indexes, genes, gene_id[indexes]);
					DEG1_Module.append(modules[indexes])
					file3.write(gene_id[indexes]+"\t"+modules[indexes]+"\n")
			k += 1

		file3.close(); file1.close();

		dict4 = Nodes().data_extract(infile3,0) 
		geneid, mod = dict4['arr_0'], dict4['arr_1']

		modset = sorted(set(mod), reverse=False); print(modset);
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

		file5.close();


