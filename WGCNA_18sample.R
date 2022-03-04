### ----------------------------- STEP I --------------------------------------------
##----------------------- Data Input and Preprocessing ------------------------------

setwd("C:/Users/thula/OneDrive/Documents/Thulani/University/Winter_semester/Research/Co-expression_network")
library(WGCNA)
options(stringsAsFactors = FALSE)

#Reading data
exp = read.csv('Bni.csv',header = T,row.names = 1) 

#Preprocessing
data = t(exp)
data = as.data.frame(data)
rm(exp)

#QC - checking for missing values
gsg = goodSamplesGenes(data, verbose = 2);
gsg$allOK

if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(data)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data)[!gsg$goodSamples], collapse = ", ")))
  data = data[gsg$goodSamples, gsg$goodGenes]
}
#Cluster samples
sampleTree = hclust(dist(data), method = "average")

# Plot the sample tree as a dendrogram
library(grDevices)
pdf(file = "Plots/1-sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
#No outliers, therefore no need of removing samples

#saving tree
library(ape)
nj_tree = nj(dist(data))
write.tree(nj_tree,"Samples-tree.nwk")
#plotClusterTreeSamples(data)

#Plot the samples as a heatmap
label = c("Bni Pol Sub1-a", "Bni Pol Sub1-b","Bni Pol Sub1-c","Bni Pol Sub2-a", "Bni Pol Sub2-b","Bni Pol Sub2-c",
          "Bni Pol Sub3-a","Bni Pol Sub3-b","Bni Pol Sub3-c","Bni Stig Sub1-a", "Bni Stig Sub1-b","Bni Stig Sub1-c","Bni Stig Sub2-a", "Bni Stig Sub2-b","Bni Stig Sub2-c",
          "Bni Stig Sub3-a","Bni Stig Sub3-b","Bni Stig Sub3-c") 
main = c("Bni Pol Sub1", "Bni Pol Sub1", "Bni Pol Sub1", "Bni Pol Sub2", "Bni Pol Sub2", "Bni Pol Sub2","Bni Pol Sub3", "Bni Pol Sub3","Bni Pol Sub3","Bni Stig Sub1", "Bni Stig Sub1", "Bni Stig Sub1", "Bni Stig Sub2", "Bni Stig Sub2", "Bni Stig Sub2","Bni Stig Sub3", "Bni Stig Sub3","Bni Stig Sub3")
dt = as.matrix(dist(data, method = 'euclidean'))
row.names(dt) = label
colnames(dt) = label

library(RColorBrewer)
pdf(file = "Plots/1-sampleClustering_heatmap.pdf", width = 20, height = 15)
par(cex = 1.6)
par(mar=c(1,4,4,1)) 

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
my_group <- as.numeric(factor(main))
colSide <- brewer.pal(4, "Pastel1")[my_group]
colMain <- colorRampPalette(brewer.pal(8, "BrBG"))(30)
heatmap(dt, col = colMain,  RowSideColors=colSide,  margins = c(8,2),cexRow = 1.7, cexCol = 1.7)
legend(x = "topleft", legend=c("Minimum", "Average", "Maximum"), 
       fill=colorRampPalette(brewer.pal(8, "BrBG"))(3))
dev.off()
collectGarbage()

#saving data
save(data, file = "dataInput.RData")


### --------------------------------- STEP II ---------------------------------------
##-------------------- Network Construction & Module Detection -----------------------

#rm(list = ls(all.names = TRUE))

setwd("C:/Users/thula/OneDrive/Documents/Thulani/University/Winter_semester/Research/Co-expression_network")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(grDevices)

# Load the data saved in the first part
load(file = "dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=3))
# Call the network topology analysis function
sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)

# Plot the results:
pdf(file = "Plots/2-thresholding.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

##We choose the power 15, which is the lowest power for which the scale-free topology 
#fit index curve flattens out upon reaching a high value (in this case, roughly 0.90)

#convert data to numeric
data[] <- lapply(data, as.numeric)
## One-step network construction and module detection
net = blockwiseModules(data, power = 15, corType = "pearson", networkType = "signed",
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)

#SIGNED NETWORK is selected as it is more accurate. 
#https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/

#number of modules and module sizes
table(net$colors)

#plotting modules
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
write.table(table(mergedColors), 'Modules.txt', quote = F, col.names = F)

# Plot the dendrogram and the module colors underneath
pdf(file = "Plots/3-Modules.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#saving the environment
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(data, moduleColors)$eigengenes
MEs = removeGreyME(MEs0,  greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction.RData")

### --------------------------------- STEP III ---------------------------------------
##----------------------- Network visualization -------------------------------------

#rm(list = ls(all.names = TRUE))

setwd("C:/Users/thula/OneDrive/Documents/Thulani/University/Winter_semester/Research/Co-expression_network")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(grDevices)

#load data
load(file = "dataInput.RData")

# Load network data saved in the second part.
load(file = "networkConstruction.RData")

#Number of genes and samples
nGenes = ncol(data)
nSamples = nrow(data)

#Visualization using Heatmaps

# Calculate topological overlap 
data[] <- lapply(data, as.numeric)
TOM = TOMsimilarityFromExpr(data, power = 15)
save(TOM,file = 'TOM_signed.RData')
load('TOM_signed.RData')
dissTOM = 1-TOM

# For reproducibility, we set the random seed
set.seed(10)
nSize = 1000
select = sample(nGenes, size = nSize)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

plotDiss = selectTOM^7           #makes the plot more informative
diag(plotDiss) = NA             #improves the clarity of the plot
pdf(file = "Plots/4-Heatmap.pdf", width = 12, height = 9)
TOMplot(plotDiss, selectTree, selectColors, 
        main = "Network heatmap plot, selected genes")
dev.off()


### ------------------------------- STEP V -------------------------------------------
##------------------------ Exporting the Network -------------------------------------

#rm(list = ls(all.names = TRUE))

setwd("C:/Users/thula/OneDrive/Documents/Thulani/University/Winter_semester/Research/Co-expression_network")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#load data
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")
load('TOM_signed.RData')

#Hub genes
hub = chooseTopHubInEachModule(data, moduleColors, omitColors = "grey", 
                               power = 15, type = "signed")
write.table(hub,'Hub genes.txt',quote = F, col.names = F, row.names = F)

#EXPORTING TO CYTOSCAPE
setwd("C:/Users/thula/OneDrive/Documents/Thulani/University/Winter_semester/Research/Co-expression_network")
clr = unique(moduleColors)
clr = sort(clr)
clr = clr[clr != "grey"]
probes = names(data)

# Select modules
for (modules in clr) {
  # Select module probes
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  modTOM = as.table(modTOM)
  write.csv(modTOM,file = paste("TOM",modules,".csv"))
  
  # Export the network into edge and node list files that Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", 
                                                  paste(modules, collapse="-"), 
                                                  ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.01,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule]);
}

#genes in the network after thresholding
count = numeric()
mod = character()

for (i in 1: length(clr)) {
  modules = clr[i]
  path = paste("CytoscapeInput-edges-",modules,'.txt',sep = "")
  x = read.delim(path)
  dt = unique(c(as.character(x$fromNode,x$toNode)))
  count[i] = length(dt)
  mod[i] = modules
}

summary = cbind(mod,count)
write.table(summary,file = 'Gene count.txt',row.names = F,quote = F)
head(summary)


### ------------------------------- STEP VI ------------------------------------------
##-------------------------- Validating the Network ----------------------------------

#rm(list = ls(all.names = TRUE))

setwd("C:/Users/thula/OneDrive/Documents/Thulani/University/Winter_semester/Research/Co-expression_network")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#load data
load(file = "dataInput.RData")
load(file = "networkConstruction.RData")

softPower = 15
adjacency = adjacency(data, power = softPower)

setLabels = c("Network", "Test")
multiExpr = list(Network = list(data = adjacency), Test = list(data = data))
multiColor = list(Network = moduleColors, Test = moduleColors);
nSets = 2

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = c(1:2),
                          nPermutations = 200,
                          randomSeed = 1,
                          verbose = 3)
} );
# Save the results
save(mp, file = "ModulePreservation.RData")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
mp_sum = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
write.csv(mp_sum,'Module Preservation.csv',quote = F)

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
library(grDevices)
pdf(file = "Plots/7-modulePreservation.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off()
