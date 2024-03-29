###########################################################
#           This file contain UI widgets.                ##
###########################################################

###########################################################
#                  Umap plotting                         ##
###########################################################
umap_panel <- tabPanel(id="umap_panel",
  
    titlePanel(h5("scClusters")),
    sidebarPanel(
      titlePanel(h5('UMAP Name', align = 'center')),
      width = 3,
      h4(''),
      hr(style = "border-color: grey"),
      
      selectizeInput(
        'UMAP1_forComparison',
        label = 'UMAP 1',
        choices = c("Sample","Clusters","Unconstrained","Constrained","Constrained remap"),
        selected = "Clusters"
      ),
      
      selectizeInput(
        'UMAP2_forComparison',
        label = 'UMAP 2',
        choices = c("Clusters","Sample","Unconstrained","Constrained","Constrained remap"),
        selected = "Sample"
      ),
      
      hr(style = "border-color: grey"),
      
      splitLayout(cellWidths = c("30%","30%","40%"),
                  numericInput("UMAP1_plot_width", "Width", min = 0, max = 250, value = 8),
                  numericInput("UMAP1_plot_height", "Height", min = 0, max = 250, value = 12),
                  selectizeInput(
                    'plot_choice_download_UMAP1',
                    label = "Format",
                    choices = c(".pdf",".png",".tiff"),
                    selected = ".pdf"),    
                  tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
                  ),
      
      downloadButton(outputId = "download_UMAP1", label = "Download UMAP 1"),
      
      hr(style = "border-color: grey"),
      
      splitLayout(cellWidths = c("30%","30%","40%"),
                  numericInput("UMAP2_plot_width", "Width", min = 0, max = 250, value = 8),
                  numericInput("UMAP2_plot_height", "Height", min = 0, max = 250, value = 12),
                  selectizeInput(
                    'plot_choice_download_UMAP2',
                    label = "Format",
                    choices = c(".pdf",".png",".tiff"),
                    selected = ".pdf"),    
                  tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
      ),
      downloadButton(outputId = "download_UMAP2", label = "Download UMAP 2"),
      
    ),

mainPanel(
  
  
                fluidRow(h5("Dimension Reduction scClusters UMAPs"
                            )),
                fluidRow(helpText("Users can view and compare side-by-side UMAPs' representing identified scATAC-seq clusters,
                            origin of sample, unconstrained and constrained integration with scRNA-seq datasets, and integrated remapped clusters.", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
                ),
                fluidRow(
                column(6,plotOutput("UMAP_plot_1")),  ##%>% withSpinner(color="#0dc5c1")
                column(6,plotOutput("UMAP_plot_2"))
              )
        )
)
###########################################################
#                  Plot Browser:scATAC Clusters          ##
###########################################################

scATACbrowser_panel <- tabPanel(
    
    titlePanel(h5("scATAC-seq peak browser")),
    
    sidebarPanel(
      titlePanel(h5('Gene Name', align = 'center')),
      width = 3,
      h4(''),
      hr(style = "border-color: grey"),
      actionButton(inputId = "restartButton", label = "Plot Track", icon = icon("play-circle")),
      
      selectizeInput(
        'gene_name',
        label = 'Gene Name',
        choices = sort(gene_names),
        selected = sort(sort(gene_names))[1]
      ),
      selectizeInput(
        'group_by',
        label = 'Group By',
        choices = c('Clusters','Sample'),
        selected = 'Clusters'
      ),
      sliderInput("range", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("range_min", "Distance (-kb):", min = -250, max = 250, value = -50),
          numericInput("range_max", "Distance (+kb):", min = -250, max = 250, value = 50)
        ),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("tile_size", "TileSize:", min = 10, max = 5000, value = 250),
          numericInput("ymax", "Y-Max (0,1):", min = 0, max = 1, value = 0.99)
        ),
      
      hr(style = "border-color: grey"),

        splitLayout(cellWidths = c("30%","30%","40%"),
          numericInput("plot_width", "Width", min = 0, max = 250, value = 8),
          numericInput("plot_height", "Height", min = 0, max = 250, value = 12),
          selectizeInput(
            'plot_choice_download_peakBrowser',
            label = "Format",
            choices = c(".pdf",".png",".tiff"),
            selected = ".pdf"),    
          tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
        ),
      downloadButton(outputId = "down", label = "Download"),
      
    ),
      
    mainPanel(fluidRow(h5("Peak browser of scATAC-seq clusters"
                           )),
    fluidRow(helpText("Users can view and compare the single-cell chromatin accessibility data in scalable peak browser view along with co-accessibility of peaks on scATAC-seq modality.", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
    ),
            plotOutput("browser_atacClusters")
        )
)
###########################################################
#                  Peak2GeneLinks:BrowserView            ##
###########################################################



peak2gl_panel <- tabPanel(
    
    titlePanel(h5("Peak2GeneLinks")),
    sidebarPanel(
      titlePanel(h5('Gene Name', align = 'center')),
      width = 3,
      h4(''),
      hr(style = "border-color: grey"),
      actionButton(inputId = "restartButton_1", label = "Plot Track", icon = icon("play-circle")),
      
      selectizeInput(
        'gene_name_1',
        label = 'Gene Name',
        choices = sort(gene_names),
        selected = sort(gene_names)[1]
      ),

      sliderInput("range_1", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("range_min_1", "Distance (-kb):", min = -250, max = 250, value = -50),
          numericInput("range_max_1", "Distance (+kb):", min = -250, max = 250, value = 50)
        ),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("tile_size_1", "TileSize:", min = 10, max = 5000, value = 250),
          numericInput("ymax_1", "Y-Max (0,1):", min = 0, max = 1, value = 0.99)
        ),
      
      hr(style = "border-color: grey"),

        splitLayout(cellWidths = c("30%","30%","40%"),
          numericInput("plot_width_1", "Width", min = 0, max = 250, value = 8),
          numericInput("plot_height_1", "Height", min = 0, max = 250, value = 12),
          selectizeInput(
            'plot_choice_download_peak2GeneLink',
            label = "Format",
            choices = c(".pdf",".png",".tiff"),
            selected = ".pdf"),    
          tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
        ),
      
      downloadButton(outputId = "down_1", label = "Download"),
      
    ),
      
    mainPanel
              (fluidRow(h5("Browser view of Peak2GeneLinks"
      )),
      fluidRow(helpText("User can visualize genome accessibility tracks of marker genes with peak co-accessibility", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
      ),
            plotOutput("co_access_peaks")            
        )
)

###########################################################
#                  Feature comparison:UMAP plot          ##
###########################################################
feature_comparison_panel <- tabPanel(
  
  titlePanel(h5("Feature of interest UMAPs")),
  sidebarPanel(
    titlePanel(h5('Feature Name', align = 'center')),
    width = 3,
    h4(''),
    hr(style = "border-color: grey"),
    actionButton(inputId = "gene_to_gene_restartButton", label = "Plot", icon = icon("play-circle")),
    
    selectizeInput(
      'matrix_forComparison',
      label = 'Matrix Type',
      choices = c("GeneScoreMatrix","GeneIntegrationMatrix","MotifMatrix"),
      selected = "GeneScoreMatrix"
    ),
    
    selectizeInput(
      'gene_forComparison_1',
      label = 'Feature Name 1',
      choices = "",
      selected = NULL
    ),
    selectizeInput(
      'gene_forComparison_2',
      label = 'Feature Name 2',
      choices = "",
      selected = NULL
    ),
    
    hr(style = "border-color: grey"),
    splitLayout(cellWidths = c("30%","30%","40%"),
                numericInput("gene_Comparison_plot_width", "Width", min = 0, max = 250, value = 8),
                numericInput("gene_Comparison_plot_height", "Height", min = 0, max = 250, value = 12),
                selectizeInput(
                  'plot_choice_download_feature_comparison',
                  label = "Format",
                  choices = c(".pdf",".png",".tiff"),
                  selected = ".pdf"),    
                tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
    ),
    downloadButton(outputId = "download_feature_comparison", label = "Download"),
    
    hr(style = "border-color: grey"),
    selectizeInput(
      'motif_for_motifPos',
      label = 'Download Motif Position',
      choices = sort(motifMatrix_dropdown),
      selected = sort(motifMatrix_dropdown)[1]
    ),
    downloadButton(outputId = "download_motifPos", label = "Download"),
    
  ),
  
  mainPanel
    (
      fluidRow(h5("Feature of interest : Dimensionality Reduction UMAPs"
    )),
    fluidRow(helpText("Users can view and compare side-by-side UMAPs representing features of interest in GeneScoreMatrix, GeneIntegrationMatrix or MotifMatrix with a representative sequence logo. Download list of Motif Positions. ", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
    ),
    
    plotOutput("feature_comparison"))
)


##########################################################
              #   Motif footprinting  plot              ##
##########################################################
motif_Footprinting_panel<- tabPanel(
  titlePanel(h5("TF Footprinting")),
  sidebarPanel(
    titlePanel(h5('Gene ', align = 'center')),
    width = 3,
    h4(''),
    hr(style = "border-color: grey"),
    actionButton(inputId = "plotFootprints_button", label = "Plot Footprints", icon = icon("play-circle")),

    selectizeInput(
      'motifName_input',
      label = 'Motif Name',
      choices = sort(names(motifPositions)),
      selected = NULL,
      options = list(maxItems = 1)
      ),

    selectizeInput(
      'normMethod_Input',
      label = 'Normalization Method',
      choices = c('None','Divide',"Subtract"),
      selected = 'Subtract'
     ),

    splitLayout(cellWidths = c("30%","30%","40"),
                numericInput("motif_plot_width_1", "Width", min = 0, max = 250, value = 8),
                numericInput("motif_plot_height_1", "Height", min = 0, max = 250, value = 12),
                selectizeInput(
                  'plot_choice_download_motif_down_1',
                  label = "Format",
                  choices = c(".pdf",".png",".tiff"),
                  selected = ".pdf"),    
                tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
    ),
    hr(style = "border-color: grey"),
    downloadButton(outputId = "motif_down_1", label = "Download"),
  ),

  mainPanel(
    fluidRow(h5("TF Footprinting"
    )),
    fluidRow(helpText("Users can view Tn5 bias-adjusted TF footprints for identified transcription factor motifs of interest.", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
    ),
    plotOutput("motifPlot")
  )
)


###########################################################
#                  Trajectory Heatmap plot               ##
###########################################################
traj_heatmap_panel <- tabPanel(
  titlePanel(h5("Pseudotime trajectories")),
  sidebarPanel(
    titlePanel(h5('Pseudotime heatmaps ', align = 'center')),
    width = 3,
    h4(''),
    hr(style = "border-color: grey"),
    actionButton(inputId = "plotheat_traj_button", label = "Plot Heatmaps", icon = icon("play-circle")),
    
    selectizeInput(
      'matrix_forTraj1',
      label = 'Matrix Type',
      choices = c("GeneScoreMatrix","GeneIntegrationMatrix","MotifMatrix","PeakMatrix"),
      selected = "GeneScoreMatrix"
    ),
    
    splitLayout(cellWidths = c("30%","30%","40%"),
                numericInput("matrix_forTraj_width", "Width", min = 0, max = 250, value = 8),
                numericInput("matrix_forTraj_height", "Height", min = 0, max = 250, value = 12),
                selectizeInput(
                  'plot_choice_down_heatmap_traj1',
                  label = "Format",
                  choices = c(".pdf",".png",".tiff"),
                  selected = ".pdf"),    
                tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
    ),
    
    downloadButton(outputId = "down_heatmap_traj1", label = "Download"),
    
    hr(style = "border-color: grey"),
    
    selectizeInput(
      'matrix_forTraj2',
      label = 'Matrix Type',
      choices = c("GeneScoreMatrix","GeneIntegrationMatrix","MotifMatrix","PeakMatrix"),
      selected = "GeneIntegrationMatrix"
    ),
    
    splitLayout(cellWidths = c("30%","30%","40%"),
                numericInput("matrix_forTraj_width2", "Width", min = 0, max = 250, value = 8),
                numericInput("matrix_forTraj_height2", "Height", min = 0, max = 250, value = 12),
                selectizeInput(
                  'plot_choice_down_heatmap_traj2',
                  label = "Format",
                  choices = c(".pdf",".png",".tiff"),
                  selected = ".pdf"),    
                tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
    ),
    downloadButton(outputId = "down_heatmap_traj2", label = "Download"),
    
    #downloadButton(outputId = "download_UMAP11", label = "Download UMAP 1"),
    
  ),
  
  mainPanel(
    fluidRow(h5("ArchR defined trajectory heatmaps"
    )),
    fluidRow(helpText("Users can visualize heatmaps from ArchR defined trajectory analysis on four different computed matrices. ", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
    ),
    fluidRow(
      column(5,
             plotOutput("traj_heatmap1")),
      column(2,offset=0),
      column(5,
             plotOutput("traj_heatmap2")),
      
    )
  )
)

###########################################################
#                  peak2genelink Heatmap plot            ##
###########################################################
peak2GLheatmap_panel <- tabPanel(
  titlePanel(h5("Peak2GeneLink heatmaps")),
  sidebarPanel(
    titlePanel(h5('Peak2GeneLink heatmaps', align = 'center')),
    width = 3,
    h4(''),
    hr(style = "border-color: grey"),
  actionButton(inputId = "heatmap_peak2gl", label = "Plot Heatmap", icon = icon("play-circle")),
  
  hr(style = "border-color: grey"),
  splitLayout(cellWidths = c("30%","30%","40%"),
              
              numericInput("p2g_plot_width_1", "Width", min = 0, max = 250, value = 8),
              numericInput("p2g_plot_height_1", "Height", min = 0, max = 250, value = 12),
              selectizeInput(
                'plot_choice_down_heatmap_p2g',
                label = "Format",
                choices = c(".pdf",".png",".tiff"),
                selected = ".pdf"),    
              tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;}")))
  ),
  downloadButton(outputId = "heatmap_p2g", label = "Download"),
  
  ),
  mainPanel(
   
    fluidRow(h5("Peak2GeneLinks heatmaps"
    )),
    fluidRow(helpText("Users can view Peak2GeneLinks identified across the dataset with ArchR. ", style = "font-family: 'Helvetica Now Display Bold'; font-si20pt"),
    ),
    fluidRow(
      column(1,offset=2),
      column(10,
             plotOutput("heatmap_peak_to_gene")),
    )
  )
)
#Users can view and compare putative cis-regulatory elements based on correlated peak 
#accessibility and gene expression heatmaps
###########################################################
#                        about                          ##
###########################################################

about_panel <- tabPanel(

  titlePanel(h5("About")),
  # Various tabs.
  tabsetPanel(
    # General info.
    tabPanel(
      "Overview",
      tags$h3("Scope"),
      tags$p(HTML("ShinyArchR.UiO is a user-friendly, integrative open-source Shiny-based web app using R programming for visualization of massive single-cell chromatin accessibility data (scATAC-seq) based on <a href=\"https://www.archrproject.com\" target=\"_blank\">ArchR</a> (Corces et al., 2021).")),
      tags$h3("Approach"),
      
      tags$p(HTML(" The ArchR objects saved in folders along with HDF5 formatted Arrow files are used for input in ShinyArchR.UiO.")),
      tags$h5("Data Visualization of ShinyArchR.UiO:"),
      tags$ul(
        tags$li(HTML("scATAC-seq clusters, unconstrained and constrained clusters on integrated reduced dimensions UMAP from ArchR objects")),
        tags$li(HTML("Peaks using plot browser tracks on clusters on scATAC-seq modality")),
        tags$li(HTML("Peaks2Genelinks tracks on single-cell RNA sequencing (scRNA-seq) integrated data with scATAC-seq using plot browser tracks. The co-accessibility among the genes can be visualized in the bottom panel")),
        tags$li(HTML("Heatmaps of pseudo time trajectory")),
        tags$li(HTML("Heatmaps of top 50 markers Peak2Genelinks in scATAC-seq and scRNA-seq")),

        )
    ),
    

    # About us .
    tabPanel(
      "About us ",
      tags$h3("Contributions and Citation info"),
      tags$p(HTML("ShinyArchR.UiO software is developed at <a href=\"https://www.med.uio.no/cancell/english/groups/chromatin-biology\" target=\"_blank\">Chromatin Biology</a> Lab at <a href=\"https://www.uio.no\" target=\"_blank\">University of Oslo</a>, as an open-source project mainly under the GPL license version 3 (see source code for details).")),
      tags$p(HTML("If ShinyArchR.UiO in any way help you in visualizing and sharing your research work such that it warrants a citation, please cite the ShinyArchR.UiO preprint in BioRXiv or the final publication.")),
    ),
    
    # Contact
    tabPanel(
      "Contact",
      br(),br(),
      sidebarPanel(tags$h3("Prof. Dr. Ragnhild Eskeland"),
                  
                   fluidRow(
                     column(1,img(src='re.jpeg', align = "left")),
                     column(3), 
                     column(8,
                            tags$h4(HTML("Associate Professor")),
                            tags$h5(HTML("<a href=\"https://www.med.uio.no/imb/personer/vit/ragnhesk/index.html\" target=\"_blank\">Chromatin Biology Group</a>")),
                            tags$h5(HTML("<a href=\"https://www.uio.no\" target=\"_blank\">University of Oslo</a>")),
                            tags$h5(HTML("Email : <a href=\"mailto:ragnhild.eskeland@medisin.uio.no\" target=\"_blank\">ragnhild.eskeland@medisin.uio.no</a>")),
                     )
                   ),
                   width = 6),
      
      sidebarPanel(tags$h3("Dr. Ankush Sharma"),
                   
                   fluidRow(
                     column(1,img(src='as.jpg', align = "left",
                                  height="200px",style="object-fit:contain")),
                     column(3), 
                     column(7,
                            tags$h4(HTML("Researcher (Bioinformatics)")),
                            tags$h5(HTML("Chromatin Biology Group")),
                            tags$h5(HTML("<a href=\"https://www.uio.no\" target=\"_blank\">University of Oslo</a>")),
                            tags$h5(HTML("Email : <a href=\"mailto:ankush.sharma@medisin.uio.no\" target=\"_blank\">ankush.sharma@medisin.uio.no</a>")),
                     )
                   ),
                   width = 6),
      
      sidebarPanel(tags$h3("Akshay Suhag"),
                   
                   fluidRow(
                     column(1,img(src='aa.png', align = "left",
                                  width="150px",style="object-fit:contain")),
                     column(3), 
                     column(7,
                            tags$h4(HTML("Doctoral Researcher (Bioinformatics)")),
                            tags$h5(HTML("<a href=\"https://www.urofun.ch\" target=\"_blank\">Functional Urology Group</a>")),
                            tags$h5(HTML("<a href=\"https://www.dbmr.unibe.ch/about_us/staff/personenpool_index/akshay_akshay/index_eng.html\" target=\"_blank\">University of Bern</a>")),
                            tags$h5(HTML("Email : <a href=\"mailto:akshay.akshay@unibe.ch\" target=\"_blank\">akshay.akshay@unibe.ch</a>")),
                     )
                   ),
                     width = 6))

  )

)
ui <- shinyUI(fluidPage(
    # Use this function somewhere in UI
    #add_busy_spinner(spin = "cube-grid", color = "#CCCCCC", onstart = TRUE, height = "65px", width = "65px"),
  add_busy_spinner(spin = "radar", color = "#CCCCCC", onstart = TRUE, height = "55px", width = "55px"),
  
navbarPage( 
  umap_panel,
    scATACbrowser_panel,
    peak2gl_panel,
    feature_comparison_panel,
    #motif_Footprinting_panel,
    traj_heatmap_panel,
    peak2GLheatmap_panel,
    about_panel,
  # Application title.
   title ="ShinyArchR.UiO",
  tags$head(tags$style(".shiny-output-error{color: grey;}"))###showing error in grey color
),

tags$footer(HTML("<p><i>Webpage generated with</i> <a href='https://github.com/EskelandLab/ShinyArchRUiO' target=\"_blank\">ShinyArchr.UiO</a>.</p>"),
            align = "left", style = "
              position:relative;
              bottom:0;
              color: black;
              padding: 10px;
              z-index: 1000;")
)
)

