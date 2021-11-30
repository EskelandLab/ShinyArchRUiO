# ShinyArchrUiO
  [ShinyArchR.UiO](https://cancell.medisin.uio.no/ShinyArchR.UiO)(ShinyArchR User interface Open)is an user-friendly, integrative, and open-source shiny-based web app using R programming for visualization of  massive single-cell chromatin accessibility data (scATAC-seq) based on [ArchR](https://archrproject.com), ([Granja et al, 2021](https://www.nature.com/articles/s41588-021-00790-6)).
  Example of web interface on tutorial dataset is available at [ShinyArchR.UiO](https://cancell.medisin.uio.no/ShinyArchR.UiO) website ([Sharma et al, Bioinformatics,2021](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab680/6377776)). Learn more from our videotutorial: https://youtu.be/gIUGgJWlWCw.
 
 ## Downsampled tutorial data
 We utilized the tutorial data downloaded using the **getTutorialData()** function for Shiny instance of ShinyArchR.UiO. The downsampled tutorial data of hematopoietic cells approximately 0.5 GB in size is used for the analysis with the steps described in full manual of ArchR toolkit (https://www.archrproject.com/bookdown/index.html). 
## Input Data for ArchR processing
ArchR can read a wide range of input formats, often in fragment files or BAM files, but it is also capable of reading scATAC-seq data. scATAC-seq fragment files contain the corresponding cell ID for each scATAC-seq fragment, sorted into tabix files. BAM files contain information about each scATAC-seq fragment, raw sequence, cellular barcode ID, and other information in tabularized format. The preprocessing pipeline defines what input format is used. The 10x Genomics Cell Ranger software, for example, returns fragment files, while sci-ATAC-seq applications use BAM files. To read fragment files, ArchR uses "scanTabix" and to read BAM files, it uses "scanBam". To support the input process, input data chunks are converted to a compressed table-based representation of fragments, which includes the fragment chromosome, offset-adjusted start and end positions, as well as the cellular barcode ID. To preserve memory consumption while maintaining quick access to chunks, chunks are stored in a temporary HDF5-formatted file. The final step involves reading, organizing, and rewriting all portions of each chromosome to an Arrow file within a single HDF5 group referred to as a "fragment". In this way, ArchR is able to handle extremely large input files efficiently and with very low memory usage, allowing it to fully utilize parallel processing.

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

### Installation of mandatory packages
Open R environment or R GUI of your choice and run the following code:
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
On command Line: 

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

#### Setting up parameters
* Open **global.R** file in a file editor and specify the following parameters:  
* ArchRThreads
```r
#Set ArchRThreads as per available computational resources in setting up the ArchR section of the file. 
#Default set to 1.  
ArchRThreads = 1
```

* Provide the path to the saved folders in global.R**
```r
savedArchRProject1 <- loadArchRProject("path to projHeme2/")
savedArchRProject2 <- loadArchRProject("path to ProjHeme3/")
savedArchRProject3 <- loadArchRProject("path to projHeme5/")
```
* Use trajectory name used in [getTrajectory](https://www.archrproject.com//bookdown/myeloid-trajectory-monocyte-differentiation.html) function instead of "LymphoidU" in global.R. Save the 
 ```r   
trajectory_name<-"LymphoidU"
```
* Save the global.R file.

### Running ShinyArchRUiO on Commandline

Navigate to the folder containing ShinyArchRUiO.  
       
    R -e "shiny::runApp('ShinyArchR.UiO',launch.browser =TRUE)" 

### Running ShinyArchRUiO on R GUI 

    shiny::runApp('ShinyArchRUiO')

***
![ShinyArchRUiO tab Information](https://github.com/EskelandLab/ShinyArchrUiO/blob/main/example_data/tab_info_wiki.png)
***
## ShinyArchR.UiO Visualization
#### 1. scClusters
* Users can select scATAC-seq clusters, unconstrained, constrained, and Remapped clusters and visualize and compare multi-dimensional 
   reduction UMAP plots side-by-side.
#### 2. scATAC-seq peak browser  
* Visualize chromatin accessibility peak browser tracks for original samples or clusters on scATAC-seq modality;
#### 3. Peak2GeneLinks 
* Visualize Peaks2Genelnks peak browsers tracks for the selected feature of interest and a bottom panel showing peak co-accessibility 
   information for a shown feature of interest;
#### 4. Feature of Interest UMAPS 
* Feature Comparision tab shows Multi-dimensional reduction Umaps allowing users to compare features of interest for **GeneScoreMatrix** and **GeneIntegrationmatrix**. 
#### 5. Pseudotime trajectories
* Visualize Pseudotime trajectory for:  
  * **GeneScorematrix**,  
  * **GeneIntegrationmatrix**,  
  * **Motifmatrix**,  
  * **Peakmatrix**  
 
#### 6. Peak2GeneLinks heatmaps. 
* Visualize heatmaps for Peak2Genelinks top markers in scATAC-seq and scRNA-seq.  

**For a more detailed description, Please see [supplementary information](https://www.biorxiv.org/content/10.1101/2021.06.21.449316v2.supplementary-material).**
## Frequently Asked Questions.  

Q: Which version of R programming is required?
* R version 4.0.0  and over is recommended.  

Q: What are the required packages for running ShinyArchRUiO?
* Please ensure shiny, ArchR, Seurat, Magick, hexbin, shinybusy, and other dependencies required are installed and loaded properly in the R environment if running within the R environment. More details about required packages and their dependencies along  can be found in session information file in Github repository .

Q: How to set reference genome other than hg19 ?
* ShinyArchR.UiO supports visualization of additional genome annotations and custom annotations. A genome is set as the basis for gene and genome annotations.  In our demo version, we utilized data aligned using hg19 genome version.  However, User can analyze data for any species by custom genome and gene annotations using the ```createGeneAnnotation()``` and ```createGenomeAnnotation()``` functions or ArchR. ArchR natively supports ```hg19```, ```hg38```, ```mm9```, and ```mm10``` and using ```addArchRGenome("hg38")``` will use ```hg38``` instead of ```hg19```. 

Q: How much memory/storage space does ShinyArchR.UiO and the Shiny app consume?
* The Shiny app itself is less memory intensive and is meant to be a heavy-duty app where multiple users can access the app at the same time. The memory required is dependent on the saved project files from ArchR. Simultaneously, ArchR employs Arrow files, an HDF5 file format, to store massive single-cell chromatin accessibility data on disk/ user’s server. Initial setup of ShinyArchr.UiO is computation-intensive.  This includes steps for computing marker genes for peak2genelinks analysis and other plots. A typical laptop with 8GB RAM can handle datasets from estimated 10k cells while 16GB RAM machines can handle around 20k-50k cells. Initialization of the app takes approximately 5-10 minutes for a downsampled example tutorial dataset of hematopoietic cells utilized in the ArchR manual.

Q: If you are getting *502 Bad Gateway* error on demo [ShinyArchR.UiO web interface](https://cancell.medisin.uio.no/ShinyArchR.UiO)? 
* Perform a hard refresh in your browser. Clear your browser cache and delete cookies. Your browser may be holding on to certain files that were saved once you visited the website with a 502 error. Please wait for 5-10 minutes,  This could be due to higher load on our server side , please let us know if the error still persists.  

Q: Specification of ShinyArchR.UiO server ?

```r 
RHEL system and we use SElinux, nginx and SSL
Model name:	Intel(R) Xeon(R) Platinum 8168 CPU @ 2.70GHz
Architecture:	 x86_64
CPU op-mode(s):	32-bit, 64-bit
CPU(s):	4
CPU family:	6
RAM:	31Gi
Icon name:	computer-vm
Virtualization:	vmware
Operating System:	Red Hat Enterprise Linux 8.4 (Ootpa)
CPE OS Name:	cpe:/o:redhat:enterprise_linux:8.4:GA
Kernel:	Linux 4.18.0-305.el8.x86_64
```

***
## Links
For a general introduction of the tool and how to setting up ShinyArchR.UiO locally. 

Please watch ShinyArchR.UiO's **[Introduction video tutorial](https://youtu.be/tSLj9CoNUrs)**.  


Please watch  ShinyArchR.UiO's **[Setup video tutorial](https://youtu.be/q1soZ4Tcyjg)**.  

## Citations information
Please cite **ShinyArchR.UiO** article published in **[ OUP Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab680/6377776)**
and **[BiorXiV preprint](https://www.biorxiv.org/content/10.1101/2021.06.21.449316v2).**

Ankush Sharma, Akshay Akshay, Marie Rogne, Ragnhild Eskeland, ShinyArchR.UiO: user-friendly, integrative and open-source tool for visualization of single-cell ATAC-seq data using ArchR, Bioinformatics, 2021;, btab680, https://doi.org/10.1093/bioinformatics/btab680
***
