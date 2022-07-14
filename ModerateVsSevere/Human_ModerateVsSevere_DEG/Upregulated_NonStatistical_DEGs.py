import re
import sys
from .Nodes_Preprocess import Nodes
from .DEG_Module_Enrichment import DEG_Enrichment

class NonStatistical_DEG:
	def UpDown_NonStatistical_DEGs(self, infile1, infile2, infile3, infile4, infile5, outfile1, outfile2):
		DEG_Enrichment().Enrichment('Input_Files/COVID19_Human_Nodes_All.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Output_Files/Human_DEG_Modules.txt', 'Output_Files/Human_DEGMod_Enrichment.txt')

		dicts1 = Nodes().data_extract(infile1, 1)
		module, observed, expected, foldchange, fishertest = [],[],[],[],[]
		stats_module, stats_foldchange, stats_fishertest = [],[],[]
		modules, observed_degs, expected_degs, foldschange, fisher_test = dicts1['arr_0'], dicts1['arr_1'], dicts1['arr_2'], dicts1['arr_3'], dicts1['arr_4']
		for i in range(len(foldschange)):		
			if (float(foldschange[i]) > 1): # DEG Enriched Modules with Fold Change is greater than 1
				module.append(modules[i])
				observed.append(observed_degs[i])
				expected.append(expected_degs[i])
				foldchange.append(foldschange[i])
				fishertest.append(fisher_test[i])
				if (float(fisher_test[i]) <= 0.05): # DEG Enriched Modules with p-value is less than equal to 0.05
					stats_module.append(modules[i])
					stats_foldchange.append(foldschange[i])
					stats_fishertest.append(fisher_test[i])

		f = open(infile2,'r')
		m=0; mod, ltag = [],[];
		for line in f:
			split_line = (line.rstrip()).split('\t')
			mod.append(str(m))
			ltag.append(split_line)
			m += 1

		dicts2 = Nodes().data_extract(infile3, 1)
		locusid, foldchanges, pvalue = dicts2['arr_0'], dicts2['arr_2'], dicts2['arr_6']
		
		dicts3 = Nodes().data_extract(infile4, 1)
		locus_id, fold_change, p_value = dicts3['arr_0'], dicts3['arr_2'], dicts3['arr_6']
		
		dicts4 = Nodes().data_extract(infile5, 1)
		locus_tag, entrez_tag = dicts4['arr_1'], dicts4['arr_4']
		
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
	NonStatistical_DEG().UpDown_NonStatistical_DEGs('Output_Files/Human_DEGMod_Enrichment.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Input_Files/Moderate_Vs_Severe.txt', 'Input_Files/Moderate_Vs_Severe_Significant.txt', 'Input_Files/Homo_sapiens.GRCh38.103.entrez.tsv', 'Output_Files/All_Transcripts_In_DEGenrichedModules.txt', 'Output_Files/UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt')


