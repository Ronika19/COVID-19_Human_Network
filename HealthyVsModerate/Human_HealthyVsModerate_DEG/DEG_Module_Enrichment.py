import re
import scipy.stats as stats
from .Data_Process import Data_Processer

class DEG_Enrichment:
	def DataPreprocess(self):
		fl_list = ['Input_Files/Healthy_Vs_Moderate_Significant.txt', 'Output_Files/Human_Cluster_DEG_BiologicalReplicates.txt']
		Cluster = Data_Processer.DEG_Cluster(fl_list[0], fl_list[1])

		modclust_list = ['Input_Files/Module_Assignment_COVID19_Human_WGCNA.tsv', 'Output_Files/COVID19_Human_Modules_Clusters.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt']
		mod_cluster = Data_Processer.Module_Clusters(modclust_list[0], modclust_list[1], modclust_list[2])

		clust_modlist = ['Output_Files/Human_Cluster_DEG_BiologicalReplicates.txt', 'Input_Files/Module_Assignment_COVID19_Human_WGCNA.tsv', 'Output_Files/Human_DEG_Modules.txt', 'Output_Files/Human_DEG_Modules.txt', 'Output_Files/Human_DEG_Genes_Modules.txt']
		Cluster_Modules = Data_Processer.DEGClusters_2_WGCNAModules(clust_modlist[0], clust_modlist[1], clust_modlist[2], clust_modlist[3], clust_modlist[4])

	def Enrichment(self, infile1, infile2, infile3, infile4, outfile): 
		self.DataPreprocess()

		nodes_file = open(infile1,'r')
		line = nodes_file.readlines()
		len_nodes = len(line[1:]); #print(len_nodes);	# Total no of Genes in WGCNA Network

		file1_1 = open(infile2,'r')
		Total_Modules = len(file1_1.readlines()); print(Total_Modules)
		file1_1.close()
		file1 = open(infile3,'r')
		line1 = file1.readline()
		Mod_Size, Mod_Size_Percent = [],[]
		while line1:
			line1 = line1.rstrip()
			split_line1 = line1.split('\t'); #print(len(split_line1), (len(split_line1)/len_nodes)*100);
			Mod_Size.append(len(split_line1))	# No of Genes or Gene Count in each WGCNA Modules = Mod_Size[1:]
			Mod_Size_Percent.append((len(split_line1)/len_nodes)*100)	# Percentage of Genes in each WGCNA Modules = Mod_Size_Percent[1:]
			line1 = file1.readline()

		file2 = open(infile4,'r')
		line2 = file2.readline()
		PAO1_dict = {}; len_PAO1_degs = 0
		while line2:
			line2 = line2.rstrip()
			split_line2 = line2.split('\t')
			PAO1_deg_mods = int(split_line2[1]); #print(PAO1_deg_mods);
			for i in range(Total_Modules):
				if (PAO1_deg_mods == i):
					PAO1_dict[split_line2[0]] = i
			len_PAO1_degs += 1	# No of PAO1 DEGs
			line2 = file2.readline()
		print(PAO1_dict); print(len_PAO1_degs); print(Mod_Size_Percent);

		outfile2 = open(outfile,'w')
		outfile2.write('Module'+'\t'+'Observed_DEGs'+'\t'+'Expected_DEGs'+'\t'+'Fold_Change'+'\t'+'Fisher_Test'+'\n')
		PAO1_degs_count = {}; PAO1_degs_counter = []
		for j in range(Total_Modules):
			PAO1_count = sum(x == j for x in PAO1_dict.values()); #print(PAO1_count);
			PAO1_degs_count[j] = PAO1_count	# dictionary where key = modules & values = DEG count for Pathogenesis
			PAO1_degs_counter.append(PAO1_count)	# List of DEG counts in each module for Pathogenesis
			PAO1_deg_mod = float(PAO1_count)	# DEG in each Module
			PAO1_nondeg_mod = float(Mod_Size[j])-float(PAO1_count)	# Non-DEG in each Module
			PAO1_deg_nonmod = float(len_PAO1_degs)-float(PAO1_count)	# DEG not in each Module
			PAO1_nondeg_nonmod = float(len_nodes)-float(PAO1_nondeg_mod)	# Non-DEG not in each Module
			oddsratio, pvalue = stats.fisher_exact([[PAO1_deg_mod, PAO1_nondeg_mod], [PAO1_deg_nonmod, PAO1_nondeg_nonmod]]); # p-values of DEGs in each Module
			PAO1_observed_degs_mod = PAO1_count;	# DEGs observed in each Module
			PAO1_expected_degs_mod = (len_PAO1_degs*Mod_Size_Percent[j])/100	# DEGs in each Module, expected by random chance
			PAO1_fold_change = float(PAO1_observed_degs_mod)/float(PAO1_expected_degs_mod)	# Fold Change of DEGs in each Module
			print("P-value = ", pvalue);
			outfile2.write(str(j)+'\t'+str(PAO1_observed_degs_mod)+'\t'+str(PAO1_expected_degs_mod)+'\t'+str(PAO1_fold_change)+'\t'+str(pvalue)+'\n')

if __name__ == "__main__":
	DEG_Enrichment().Enrichment('Input_Files/COVID19_Human_Nodes_All.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Output_Files/WGCNA_COVID19_Human_Modules_Clusters.txt', 'Output_Files/Human_DEG_Modules.txt', 'Output_Files/Human_DEGMod_Enrichment.txt')




