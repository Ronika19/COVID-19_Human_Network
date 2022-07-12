# Display the current working directory
getwd()
# "." means current working directory
workingDir = "."
setwd(workingDir)
# Load the WGCNA package
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Read the COVID19_Human dataset. t() - transpose matrix/dataframe. For .csv file use read.csv("filename.csv")
datExpr0 = t(as.data.frame(read.table(file = 'COVID19_Human_Network.tsv', sep = '\t', header = TRUE, row.names=1)))
# Take a look at the data. Rows - Samples. Columns - Genes
# head(datExpr0, n=10); cat("\n\n"); tail(datExpr0, n=10); cat("\n\n");
dim(datExpr0); cat("\n\n");  row.names(datExpr0); cat("\n\n"); colnames(datExpr0); cat("\n\n");

# Checking data for excessive missing values and outlier RNASeq samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
# datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
# dim(datExpr0); cat("\n\n");
# If it returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples
if (gsg$allOK != TRUE) {
	if (gsg$goodGenes == FALSE){
		printFlush(paste("BAD Genes : "), paste(colnames(datExpr0)))
	}
	if (gsg$goodSamples == FALSE) {
		printFlush(paste("BAD Samples : "), paste(row.names(datExpr0)))
	}
	# Remove offending genes and sample (if any)
	datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
	# Check data again
	dim(datExpr0); cat("\n\n");
}

# Next we cluster the samples (not genes) to check for any obvious outliers
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9)	# size of graphical window
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main="Sample Clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(h=3e+05, col="red")	# Plot a line to show the cut

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from=12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose=5)
printFlush(paste("Power Estimate : "), paste(sft$powerEstimate)); cat("\n\n");
# Plot results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Use established soft-thresh power to build network with adjusted settings
net = blockwiseModules(datExpr0, power = 10,TOMType = "signed", minModuleSize = 20,reassignThreshold = 0, mergeCutHeight = 0.15,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "COVID19_Human_WGCNA",verbose = 3,maxBlockSize=40000)
print(net$blockGenes); cat("\n\n"); print(net$blocks); cat("\n\n"); print(net$blockOrder); cat("\n\n"); print(net$colors); cat("\n\n");
table(net$colors); cat("\n\n");
# Write modules to file
write.table(file="Module_Assignment_COVID19_Human_WGCNA.tsv",net$colors,sep="\t"); 
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
 
#i = 1:6; j = 1:6; k = 1:6;
i = 1; j = 1; k = 1;
rdata_files = as.vector(sprintf("%s%i%s", "COVID19_Human_WGCNA-block.",i,".RData")); print(rdata_files);
edges_files = as.vector(sprintf("%s%i%s", "COVID19_Human_WGCNA_Reduced_Edges",j,".txt"))
nodes_files = as.vector(sprintf("%s%i%s", "COVID19_Human_WGCNA_Reduced_Nodes",k,".txt"))
for (rdata in 1:length(rdata_files)) {
	# Reload the TOM
	load(rdata_files[rdata])
	block_nodes = character();
	gene_nos = colnames(datExpr0)
	typeof(gene_nos)
	is.data.frame(gene_nos)
	is.vector(gene_nos)
	is.list(gene_nos)
	for (g in net$blockGenes[rdata]) {
		print(gene_nos[g])
		block_nodes = c(block_nodes,gene_nos[g])
	}
	block_nodes

	# Export edges with threshold atleast 0.06 from network for further study
	cyt = exportNetworkToCytoscape(TOM, edgeFile=edges_files[rdata], nodeFile=nodes_files[rdata], weighted=TRUE, threshold=0.05, nodeNames=block_nodes)
}



