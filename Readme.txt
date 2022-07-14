################################## Readme File ##############################

Requirements:
1. Ubuntu OS
2. Python version >= 3.6 


To run the python scripts and generate output you need to run the shell script (Network_Analysis.sh) by the following command:
sh Network_Analysis.sh

Input Files for Network analysis of COVID-19 infected patients:
1. 'KEGG_Human_AllPathways.txt' file contains the pathway annotations for Human
2. 'Module_Assignment_COVID19_Human_WGCNA.tsv' file contains the module assignment of Human genes involved in the network.
3. 'Healthy_Vs_ICU_Significant.txt' file contain the upregulated and downregulated DEGs respectively for ICU COVID-19 patients vs Healthy humans.
   'Healthy_Vs_Moderate_Significant.txt' file contain the upregulated and downregulated DEGs respectively for Moderate COVID-19 patients vs Healthy humans.
   'Healthy_Vs_Severe_Significant.txt' file contain the upregulated and downregulated DEGs respectively for Severe COVID-19 patients vs Healthy humans.
   'Moderate_Vs_Severe_Significant.txt' file contain the upregulated and downregulated DEGs respectively for Severe vs Moderate COVID-19 patients.
   'Severe_Vs_ICU_Significant.txt' file contain the upregulated and downregulated DEGs respectively for ICU vs Severe COVID-19 patients.
4. 'gene_result.txt' file contain the Entrez ids for Human genes and their corresponding gene names.
5. 'COVID19_Human_Nodes_All.txt' file contain the nodes (gene ids) of Human involved in the network.
6. 'Healthy_Vs_ICU.txt' file contains the fold changes and p-values of all Human genes for ICU COVID-19 patients vs Healthy humans.
   'Healthy_Vs_Moderate.txt' file contains the fold changes and p-values of all Human genes for Moderate COVID-19 patients vs Healthy humans.
   'Healthy_Vs_Severe.txt' file contains the fold changes and p-values of all Human genes for Severe COVID-19 patients vs Healthy humans.
   'Moderate_Vs_Severe.txt' file contains the fold changes and p-values of all Human genes for Severe vs Moderate COVID-19 patients.
   'Severe_Vs_ICU.txt' file contains the fold changes and p-values of all Human genes for ICU vs Severe COVID-19 patients.
7. 'Homo_sapiens.GRCh38.101.gtf' is the Human gtf file.
8. 'Homo_sapiens.GRCh38.103.entrez.tsv' and 'Human_GeneInfo.txt' files contain the Entrez ids and Locustags for Human genes.
9. 'COVID19_Human_Network.tsv' is the file that is taken as input by WGCNA.

Output Files:
The 'DEG_NonDEG_Module_Pathway.txt' file is the final output file containing the DEG count, LDEG count, and the percentages of DEGs and DEGs+LDEGs in each of the Modules.
'UpDownRegulatedNonDEGs_In_DEGenrichedModules.txt': This files contain the elevated and suppressed LDEGs identified in Modules enriched in DEGs.
'Human_DEGMod_Enrichment.txt': This file contains the Fold Change and P-value (Fisher test) for identification of Modules enriched in DEGs.
