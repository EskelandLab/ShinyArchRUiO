#####################################################################
#                 Setting up                                #        
#####################################################################

# if (!requireNamespace("Shiny", quietly = TRUE)) install.packages("shiny")
# ##Install ArchR,and its dependencies
# #check for devtools
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# #Install ArchR
# if (!requireNamespace("ArchR",quietly = TRUE)) devtools::install_github("GreenleafLab/ArchR",ref="master",repos=BiocManager::repositories())
# #Install ArchR dependencies
# library(ArchR)
# ArchR::installExtraPackages()
# #Install presto
# if (!requireNamespace("presto",quietly = TRUE)) devtools::install_github("immunogenomics/presto")

####
packages <- c("Seurat","shinycssloaders","hexbin","magick",
             "gridExtra", "grid","patchwork","shinybusy","ArchR","ggseqlogo")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

options(repos = BiocManager::repositories())

#####################################################################
#                        Setting up ARCHR 
#####################################################################
#specify desired number of threads , If required
addArchRThreads(threads = 4) 
#specify genome version. Default hg19 set
addArchRGenome("hg19")
set.seed(1)


#####################################################################
#                         Load Data
#####################################################################
##Load the Saved projects folders from ArchR analysis as saved in ArchR full manual using  saveArchRProject() function
#Load Saved-project folders path e.g 'Save-ArchRProject2'<- loadArchRProject("path/to/your/Save-ArchRProject2"). Save project also after trajectory analysis e.g as Save-ArchRProject5
#Please see ArchR  full manual for saveArchRProject() function or use the ArchR.RMD for your analysis provided with the source code which follows the steps illustrated in ArchR full manual. Save-ArchRProject5
# 
savedArchRProject1 <- loadArchRProject("~/Save-ProjHeme2/")
savedArchRProject2 <- loadArchRProject("~/Save-ProjHeme3/")
savedArchRProject3 <- loadArchRProject("~/Save-ProjHeme5/")


########################################################################
#                         Add Metadata of Trajectory
########################################################################
#specify trajectory name as given in ARCHR analysis
#for more details please visit section 16: https://www.archrproject.com/
trajectory_name<-"LymphoidU"
########################################################################
#                         UMAP Visualization
########################################################################

cluster_umap <- plotEmbedding(
  ArchRProj = savedArchRProject1,
  baseSize=12,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP",
  rastr = FALSE,
  size=0.5,

  )+ggtitle("Colored by scATAC-seq clusters")+theme(text=element_text(size=12), legend.title = element_text(size = 12),legend.text = element_text(size = 6))

sample_umap <- plotEmbedding(
  ArchRProj = savedArchRProject1,
  baseSize=12,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP",
  rastr = FALSE,
  size=0.5
)+ ggtitle("Colored by original identity")+theme(text=element_text(size=12), legend.title = element_text( size = 12),legend.text = element_text(size = 6))

unconstrained_umap <- plotEmbedding(
  ArchRProj = savedArchRProject1,
  embedding = "UMAP",
  reducedDims = "IterativeLSI",
  colorBy = "cellColData",
  name = "predictedGroup_Un",
  baseSize=12,

  rastr = FALSE,
  size=0.5
  )+ggtitle("UMAP: unconstrained integration")+theme(text=element_text(size=12), legend.title = element_text( size = 12),legend.text = element_text(size = 6))

constrained_umap <- plotEmbedding(
  ArchRProj = savedArchRProject1,
  reducedDims = "IterativeLSI",
  colorBy = "cellColData",
  name = "predictedGroup_Co",
  rastr = FALSE,
  baseSize=12,
  size=0.5
)+ggtitle("UMAP: constrained integration")+theme(text=element_text(size=12), legend.title = element_text( size = 12),legend.text = element_text(size = 6))


constrained_remapped_umap <- plotEmbedding(
  savedArchRProject2,
  colorBy = "cellColData",
  name = "Clusters2",
  rastr = FALSE,
  )+ggtitle("UMAP: Constrained remapped clusters")+theme(text=element_text(size=12), legend.title = element_text( size = 12),legend.text = element_text(size = 6))

umaps<-list("Clusters"= cluster_umap,"Sample"= sample_umap,"Unconstrained"=unconstrained_umap,"Constrained"=constrained_umap,"Constrained remap"=constrained_remapped_umap)
########################################################################
#                         MarkerGenes
########################################################################
#please check FDR Threshold and Log2FC  values and use it accordingly.
#for more details please visit: https://www.archrproject.com/

#Find Marker Genes
markersGS <- getMarkerFeatures(
  ArchRProj = savedArchRProject3,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList_p1 <- getMarkers(markersGS, cutOff = "FDR <= 0.1 & Log2FC >= 0.1")

# Find the all the genes available in the data
gene_names <- markerList_p1[[names(markerList_p1)[1]]]$name
for (cn in names(markerList_p1)[-1]){
  gene_names <- union(gene_names,markerList_p1[[cn]]$name)
}

########################################################################
#                    motifs for feature comparison panel
########################################################################
motifMatrix_forShiny=getMatrixFromProject(
  ArchRProj = savedArchRProject3,
  useMatrix = "MotifMatrix",
  useSeqnames = NULL,
  verbose = FALSE,
  binarize = FALSE,
  threads = getArchRThreads()
)

#motifMatrix_dropdown=sapply(strsplit(motifMatrix_forShiny@NAMES, "_"), "[", 1)
 motifMatrix_dropdown=motifMatrix_forShiny@NAMES

#get PWM of motifs and convert them to probability matrix for seqlogo:Utilized function from utils.R of https://github.com/GreenleafLab/ChrAccR
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  (2^as(x, "matrix"))*TFBSTools::bg(x)/sum(TFBSTools::bg(x))
}

########################################################################
#                         Peak2Genelinks
########################################################################
# Plot of Heatmap of Peak To Gene Links
#for more details please visits section 15.3.2: https://www.archrproject.com/
savedArchRProject3 <- addPeak2GeneLinks(
ArchRProj <- savedArchRProject3,
reducedDims <- "IterativeLSI"
   )

########################################################################
#                         motif footprinting
########################################################################
motifPositions=getPositions(savedArchRProject3)
########################################################################
#                       Heatmap:Trajectory and peak2genelink
########################################################################
p_heatmap_peak_to_gene <- plotPeak2GeneHeatmap(ArchRProj = savedArchRProject3, groupBy = "Clusters2")

trajGSM <- getTrajectory(ArchRProj = savedArchRProject3, name = trajectory_name, useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p_GeneScore_traj<- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))

trajGIM  <- getTrajectory(ArchRProj = savedArchRProject3, name = trajectory_name, useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p_GIM_traj <- plotTrajectoryHeatmap(trajGIM, pal = paletteContinuous(set = "blueYellow"))

trajMM  <- getTrajectory(ArchRProj = savedArchRProject3, name = trajectory_name, useMatrix = "MotifMatrix", log2Norm = FALSE)
p_motifMatrix_traj <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

trajPM  <- getTrajectory(ArchRProj = savedArchRProject3, name = trajectory_name, useMatrix = "PeakMatrix", log2Norm = FALSE)
p_peakMatrix_traj <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
########################################################################
#                               End
########################################################################


