# STEP-BY-STEP NETWORK CONSTRUCTION

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

#Experiment name
experiment = "sp"

#set working dir
workingDir = "/home/hugo/Dropbox/Esalq/RNAseq48hai/RNAseq/WGCNA/";
setwd(workingDir); 

# Load the WGCNA package
library(WGCNA);

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read input table
RNAseqDATA = read.csv("wgcna_sp_input.csv", row.names=1);

#Transpose input table
datExpr = as.data.frame(t(RNAseqDATA))

#WRITE OUTPUT TO VERIFY FORMAT. NOT NECESSARY.
#write.table(datExpr, file = paste(experiment,"datExpr.table.txt"), sep = ",")

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()


#=====================================================================================
#
#  Code chunk 2 - Step-by-step network construction and module detection
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3 - Co-expression similarity and adjacency
#
#=====================================================================================

#We choose the power, which is the lowest power for which the scale-free topology fit
#index reaches 0.9 (or the closest value to 1)
softPower = 20;

# Construction of adjacency matrix based on selected softPower
adjacency = adjacency(datExpr, power = softPower, type = "signed");


#=====================================================================================
#
#  Code chunk 4 - Topological Overlap Matrix (TOM)
#
#=====================================================================================


# Turn adjacency into topological overlap
# Longest amount of time running
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

save(dissTOM, file = paste(experiment,"-TOM.RData"))

rm(adjacency)
rm(TOM) #release memory from TOM object

#=====================================================================================
#
#  Code chunk 5 - Clustering using TOM
#
#=====================================================================================

myTOM = exists("dissTOM")

if (myTOM == FALSE) {
	lnames = load(file = paste(experiment,"-TOM.RData"));
}

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk 8 - Merging of modules whose expression profiles are very similar
#
#=====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


MEDissThres = 0.25
# Plot the cut line into the dendrogram
# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge

abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


sizeGrWindow(12, 9)
#pdf(file = paste(experiment,"geneDendro-3.pdf"), wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = paste(experiment,"-networkConstruction-stepByStep.RData"))


#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

# EXPORT NETWORK

#=====================================================================================
#
#  Code chunk 4 - Exporting to Cytoscape
#
#=====================================================================================

rm(RNAseqDATA)

# Recalculate topological overlap if needed
#TOM = TOMsimilarityFromExpr(datExpr, power = 6);

# Select modules
modules = unique(moduleColors);

for (i in modules){
        print (i)
        # Select module probes
        probes = names(datExpr)
        inModule = is.finite(match(moduleColors, i));
        modProbes = probes[inModule];

        #modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
        # Select the corresponding Topological Overlap
        modTOM = dissTOM[inModule, inModule];
        dimnames(modTOM) = list(modProbes, modProbes)

        # Export the network into edge and node list files Cytoscape can read

        cyt = exportNetworkToCytoscape(modTOM,
                edgeFile = paste(experiment, i, "CytoscapeInput-edges-", collapse="-", ".txt", sep=""),
                nodeFile = paste(experiment, i, "CytoscapeInput-nodes-", collapse="-", ".txt", sep=""),
                weighted = TRUE,
                threshold = 0.02,
                nodeNames = modProbes,
                #altNodeNames = modGenes,
                nodeAttr = moduleColors[inModule]);
        rm(probes)
        rm(inModule)
        rm(modProbes)
        rm(modTOM)
}
