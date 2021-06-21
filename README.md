# ShinyArchrUiO
  ShinyArchR.UiO (ShinyArchR User interface Open) is an user-friendly, integrative, and open-source shiny-based web app using R programming for   
  visualization of massive single-cell chromatin accessibility data (scATAC-seq) based on [ArchR](https://archrproject.com), ([Granja et al, 
  2021](https://www.nature.com/articles/s41588-021-00790-6)).
 
***
## Salient Features
  ShinyArchR.UiO is written in R Programming using Shiny package, enabling its use locally as well as making it available to broader 
  audiences by hosting on Shiny Server.
  The web interface has a scalable low memory footprint due to the use of the Arrow file format used by ArchR given massive single-cell ATAC- 
  seq data. 
  Users can export manuscript-ready figures in PDF.  

 
## Table of Contents and Additional Tutorials
### This readme is divided into the following sections:
* Installation
* Quick Start Guide to rapidly deploy a ShinyArchRUiO
* Frequently Asked Question
* Links and Citation info 

## Installation:
 ```r 
   Download ShinyArchR.UiO from github.com/EskelandLab/ShinyArchR.UiO
   or 
   git clone https://github.com/EskelandLab/ShinyArchR.UiO.git
```
## Quick Start guide
  The analysis performed as shown in the [ArchR full manual](https://www.archrproject.com/bookdown/index.html) on a test dataset of 
  hematopoietic cells can be applied to users' datasets.  

## Installation of mandatory packages
Open R environment or R GUI of your choice and run the following code  
```r 
    #check for devtools
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

    #Install ArchR
    devtools::install_github("GreenleafLab/ArchR",ref="master",repos=BiocManager::repositories())

    #Install ArchR dependencies
    library(ArchR)
    ArchR::installExtraPackages()

    #Installing Required packages  
    install.packages(c("shiny","Seurat", "magick","hexbin","shinybusy","gridExtra", "grid","shinycssloaders") 
``` 
On command Line 

```r 
    #Installs devtools
    Rscript -e 'install.packages("devtools",repos="http://cran.r-project.org")'
    #Installs BiocManages
    Rscript -e 'install.packages("BiocManager",repos="http://cran.r-project.org")'
    #Installs ArchR
    Rscript -e  'devtools::install_github("GreenleafLab/ArchR",ref="master",repos=BiocManager::repositories())'
    #Installs ArchR dependencies
    Rscript -e 'ArchR::installExtraPackages()'

    #Installs packages 
    Rscript -e 'install.packages(c("shiny","magick","hexbin","Seurat","shinybusy","gridExtra", "grid","shinycssloaders")' 
``` 

### Adding folder paths and parameters 

## Setting up parameters
Open **global.R** file in a file editor and specify the following parameters:  

### ArchRThreads
```r
#Set ArchRThreads as per available computational resources in setting up the ArchR section of the file. 
#Default set to 1.  
ArchRThreads = 1
```


Provide the path to the saved folders in global.R**

    savedArchRProject1 <- loadArchRProject("path to projHeme2/")
    savedArchRProject2 <- loadArchRProject("path to ProjHeme3/")
    savedArchRProject3 <- loadArchRProject("path to projHeme5/")

Use trajectory name used in [getTrajectory](https://www.archrproject.com//bookdown/myeloid-trajectory-monocyte-differentiation.html) function instead of "LymphoidU" in global.R. Save the 
    
     trajectory_name<-"LymphoidU"

Save the global.R file.

### Running ShinyArchRUiO on Commandline

Navigate to the folder containing ShinyArchRUiO.  
       
    R -e "shiny::runApp('ShinyArchR.UiO',launch.browser =TRUE)" 

### Running ShinyArchRUiO on R GUI 

    shiny::runApp('ShinyArchRUiO')

***
![ShinyArchRUiO tab Information](https://github.uio.no/ankushs/ShinyArchrUiO/blob/main/example_data/tab_info_wiki.png)
***
## ShinyArchR.UiO Visualization
#### scClusters
   Users can select scATAC-seq clusters, unconstrained, constrained, and Remapped clusters and visualize and compare multi-dimensional 
   reduction UMAP plots side-by-side;
#### scATAC-seq peak browser  
   Visualize chromatin accessibility peak browser tracks for original samples or clusters on scATAC-seq modality;
#### Peak2GeneLinks. 
   Visualize Peaks2Genelinks peak browsers tracks for the selected feature of interest and a bottom panel showing peak co-accessibility 
   information for a shown feature of interest;
#### Feature of Interest UMAPS. 
   Feature Comparision tab shows Multi-dimensional reduction Umaps allowing users to compare features of interest for **GeneScoreMatrix** and **GeneIntegrationmatrix**. 
#### Pseudotime trajectories. 
  Visualize Pseudotime trajectory for  

**GeneScorematrix**,  
**GeneIntegrationmatrix**,  
**Motifmatrix**,  
**Peakmatrix**  
 
#### Peak2GeneLinks heatmaps. 
Visualize heatmaps for Peak2Genelinks top markers in scATAC-seq and scRNA-seq.  

**For a more detailed description, Please see supplementary information (biorXiv supplementary information link).**

## Frequently Asked Questions.  

Q: Which version of R programming is required?
* R version 4.0.0  and over is recommended.  

Q: What are the required packages for running ShinyArchRUiO?
* Please ensure shiny, ArchR, Seurat, Magick, hexbin, shinybusy, and other dependencies required are installed and loaded properly in the R environment if running within the R environment. More details on the procedure to install packages and their dependencies are in the detailed tutorial.

Q: How much memory/storage space does ShinyArchR.UiO and the Shiny app consume?
* A: The Shiny app itself is less memory intensive and is meant to be a heavy-duty app where multiple users can access the app at the same time. The memory required is dependent on the saved project files from ArchR. Simultaneously, ArchR employs Arrow files, an HDF5 file format, to store massive single-cell chromatin accessibility data on disk/ user’s server. Initial setup of ShinyArchr.UiO is computation-intensive.  This includes steps for computing marker genes for peak2genelinks analysis and other plots. A typical laptop with 8GB RAM can handle datasets from estimated 10k cells while 16GB RAM machines can handle around 20k-50k cells. Initialization of the app takes approximately 5-10 minutes for a downsampled example tutorial dataset of hematopoietic cells utilized in the ArchR manual. 

***
## Links and Citation info
` ShinyArchR.UiO example website is available at https://cancell.medisin.uio.no/ShinyArchR.UiO:`. 

` Please cite:: BiorXiV.com`
***
