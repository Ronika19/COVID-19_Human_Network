import re
import sys
from .DEG_Module_Enrichment import DEG_Enrichment

class NonStatistical_DEG:
	def UpDown_NonStatistical_DEGs(self, infile1, infile2, infile3, infile4, infile5, outfile1, outfile2):
		DEG_Enrichment().Enrichment('Input_Files/COVID19_Human_Nodes_All.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Output_Files/Human_DEG_Modules.txt', 'Output_Files/Human_DEGMod_Enrichment.txt')

		file1 = open(infile1,'r')
		for i in range(2):
			line1 = file1.readline()
		module, observed, expected, foldchange, fishertest = [],[],[],[],[]
		stats_module, stats_foldchange, stats_fishertest = [],[],[]
		while line1:
			split_line1 = (line1.rstrip()).split('\t')
			if (float(split_line1[3]) > 1): # DEG Enriched Modules with Fold Change is greater than 1
				module.append(split_line1[0])
				observed.append(split_line1[1])
				expected.append(split_line1[2])
				foldchange.append(split_line1[3])
				fishertest.append(split_line1[4])
				if (float(split_line1[4]) <= 0.05): # DEG Enriched Modules with p-value is less than equal to 0.05
					stats_module.append(split_line1[0])
					stats_foldchange.append(split_line1[3])
					stats_fishertest.append(split_line1[4])
			line1 = file1.readline()

		file5 = open(infile2,'r')
		m=0; mod, ltag = [],[];
		for line5 in file5:
			split_line5 = (line5.rstrip()).split('\t')
			mod.append(str(m))
			ltag.append(split_line5)
			m += 1

		file3 = open(infile3,'r')
		for i in range(2):
			line3 = file3.readline()
		locusid, foldchanges, pvalue = [],[],[];
		while line3:
			split_line3 = (line3.rstrip()).split('\t')
			locusid.append(split_line3[0])
			foldchanges.append(split_line3[2])
			pvalue.append(split_line3[6])
			line3 = file3.readline()

		file4 = open(infile4,'r')
		for i in range(2):
			line4 = file4.readline()
		locus_id, fold_change, p_value = [],[],[];
		while line4:
			split_line4 = (line4.rstrip()).split('\t')
			locus_id.append(split_line4[0])
			#gene_name.append(split_line4[1])
			fold_change.append(split_line4[2])
			p_value.append(split_line4[6])
			line4 = file4.readline()

		file6 = open(infile5,'r')
		for i in range(2):
			line6 = file6.readline()
		locus_tag, entrez_tag = [],[];
		while line6:
			split_line6 = (line6.rstrip()).split('\t')
			locus_tag.append(split_line6[1])
			entrez_tag.append(split_line6[4]); #Entrez_Ids
			line6 = file6.readline()

		outfile1 = open(outfile1,'w')
		outfile1.write('Module'+'\t'+'Transcript_Id'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
		outfile2 = open(outfile2,'w')
		outfile2.write('Module'+'\t'+'Transcript_Id'+'\t'+'Locus_Tag'+'\t'+'Expression_Fold_Change'+'\t'+'Expression_Pvalue'+'\t'+'FoldChange_Module'+'\t'+'FisherTest_Module'+'\n')
		modid, geneids = [],[];
		for x in range(len(module)): # DEG Enriched Modules
			if str(module[x]) in mod:
				index1 = int(mod.index(str(module[x])))
				mod_ltag = ltag[index1]; 
				for t in range(len(mod_ltag)): # All transcripts in DEG Enriched Modules
					outfile1.write(str(module[x])+'\t'+mod_ltag[t]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')
					modid.append(str(module[x]))
					geneids.append(mod_ltag[t])
					if (mod_ltag[t] in locus_tag):
						index2 = int(locus_tag.index(mod_ltag[t]))
						entrez_gene = entrez_tag[index2]
						if (mod_ltag[t] in locusid):
							if (mod_ltag[t] in locus_id): # Up & Down Regulated DEGs in paper
								index3 = int(locus_id.index(mod_ltag[t]))
								#outfile2.write(str(module[x])+'\t'+mod_ltag[t]+'\t'+locus_gene+'\t'+fold_change[index3]+'\t'+p_value[index3]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')
							elif (mod_ltag[t] not in locus_id): # Up & Down Regulated DEGs not in paper
								index4 = int(locusid.index(mod_ltag[t]))
								outfile2.write(str(module[x])+'\t'+mod_ltag[t]+'\t'+entrez_gene+'\t'+foldchanges[index4]+'\t'+pvalue[index4]+'\t'+foldchange[x]+'\t'+fishertest[x]+'\n')

if __name__ == "__main__":
	NonStatistical_DEG().UpDown_NonStatistical_DEGs('Output_Files/Human_DEGMod_Enrichment.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Input_Files/Healthy_Vs_Moderate.txt', 'Input_Files/Healthy_Vs_Moderate_Significant.txt', 'Input_Files/Homo_sapiens.GRCh38.103.entrez.tsv', 'Output_Files/All_Transcripts_In_DEGenrichedModules.txt', 'Output_Files/UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt')


