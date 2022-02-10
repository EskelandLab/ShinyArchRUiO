###########################################################
#This file contain server file functions for computation.##
###########################################################
shinyServer <- function(input,output, session){
 #observe_helpers() 

  
###########################################################
#                   UMAPS                                ##
###########################################################
#Observe the inputs for UMAPS

#Output Handler: Downloads UMAPS
  output$download_UMAP1<-downloadHandler(
    filename <- function(){
      paste0("UMAP-",input$UMAP1_forComparison,input$plot_choice_download_UMAP1)
    },
    content = function(file){
      if(input$plot_choice_download_UMAP1==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$UMAP1_plot_width, height = input$UMAP1_plot_height)}
      
      else if(input$plot_choice_download_UMAP1==".png")
      {png(file = file, width = input$UMAP1_plot_width, height = input$UMAP1_plot_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$UMAP1_plot_width, height = input$UMAP1_plot_height,units="in",res=1000)}
      
      # grid.arrange(umaps[[input$UMAP1_forComparison]]) # <-Original
      grid.arrange(umaps[[input$UMAP1_forComparison]])  # <- Modification
      dev.off()
    }
  )

  output$download_UMAP2<-downloadHandler(
    filename <- function(){
      paste0("UMAP-",input$UMAP2_forComparison,input$plot_choice_download_UMAP2)
    },
    content = function(file){
      
      if(input$plot_choice_download_UMAP2==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$UMAP2_plot_width, height = input$UMAP2_plot_height)}
      
      else if(input$plot_choice_download_UMAP2==".png")
      {png(file = file, width = input$UMAP2_plot_width, height = input$UMAP2_plot_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$UMAP2_plot_width, height = input$UMAP2_plot_height,units="in",res=1000)}
      
      grid.arrange(umaps[[input$UMAP2_forComparison]])  #  <- Modification
      dev.off()
    }
  )

  output$UMAP_plot_1 <- DT::renderDT(NULL)
  output$UMAP_plot_2 <- DT::renderDT(NULL)
  

    if (isolate(input$UMAP1_forComparison) == "" || 
        isolate(input$UMAP2_forComparison) == ""){
      
      output$UMAP_plot_1 <- renderPlot({
        p <- ggplot() +
          xlim(c(-5,5)) + ylim(c(-5,5)) +
          geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
        print(p)
      })
    }else{
      # Plots UMAPS
      output$UMAP_plot_1<- renderPlot({
        umaps[input$UMAP1_forComparison]},height = 450,width=450)
      output$UMAP_plot_2<- renderPlot({
        umaps[input$UMAP2_forComparison]},height = 450,width=450)
    }
   
###########################################################
#               scATACseq Cluster:Browserview            ##
###########################################################  
  # Observe the inputs for ATAC-Seq Explorer
  observeEvent(input$range_min, {
    updateSliderInput(session, "range",
                      value = c(input$range_min,max(input$range)))
  })

  observeEvent(input$range_max, {
    updateSliderInput(session, "range",
                      value = c(input$range_min,input$range_max))
  })

  observeEvent(input$range , {

    updateNumericInput(session, "range_min", value = min(input$range))
    updateNumericInput(session, "range_max", value = max(input$range))

  }, priority = 200)
  
  # Output Handler:downloads file 
  output$down<-downloadHandler(
    filename <- function(){
      paste0("ArchRBrowser-",input$gene_name,input$plot_choice_download_peakBrowser)
    },
    content = function(file){
      
      if(input$plot_choice_download_peakBrowser==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$plot_width, height = input$plot_height)}
      
      else if(input$plot_choice_download_peakBrowser==".png")
      {png(file = file, width = input$plot_width, height = input$plot_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$plot_width, height = input$plot_height,units="in",res=1000)}
      

      p_browser_atacClusters<- plotBrowserTrack(
        ArchRProj = savedArchRProject3,
        groupBy = isolate(input$group_by),
        baseSize = 11,
        facetbaseSize = 11,
        geneSymbol = isolate(input$gene_name),
        upstream = -min(isolate(input$range))*1000,
        downstream = max(isolate(input$range))*1000,
        tileSize = isolate(input$tile_size),
        ylim =  c(0, isolate(input$ymax)),
        loops = getCoAccessibility(savedArchRProject3)
        
        )[[input$gene_name]]
      
      # grid::grid.draw(p_browser_atacClusters)  # <- Orignal
      grid.arrange(p_browser_atacClusters)  # <- modification

      dev.off()
    }
  )
  output$browser_atacClusters <- DT::renderDT(NULL)
 
 ##handles error 
  restartFN <- observeEvent(input$restartButton, {
    if (isolate(input$gene_name) == ""){

      output$browser_atacClusters <- renderPlot({
        p <- ggplot() +
          xlim(c(-5,5)) + ylim(c(-5,5)) +
          geom_text(size=20, aes(x = 0, y = 0, label = "Please supply\na valid gene name!")) + theme_void()
        print(p)
      })
    }else{
         
  # Plots scATACSeq clusters 
      output$browser_atacClusters<- renderPlot({
        grid::grid.newpage()
        p_browser_atacClusters<- plotBrowserTrack(
          ArchRProj = savedArchRProject3,
          baseSize = 11,
          facetbaseSize = 11,
          groupBy = isolate(input$group_by),
          geneSymbol = isolate(input$gene_name),
          upstream = -min(isolate(input$range))*1000,
          downstream = max(isolate(input$range))*1000,
          tileSize = isolate(input$tile_size),
          ylim =  c(0, isolate(input$ymax)),
          loops = getCoAccessibility(savedArchRProject3)
          
          )[[input$gene_name]]
        
        grid::grid.draw(p_browser_atacClusters)
        
      },height = 900)

    }
  })
  
  
###########################################################
#scATACseq coaccessibility and peak2 genelinks:Browserview##
########################################################### 
  # Observe the inputs for ATAC-Seq  Co-accessibility
  observeEvent(input$range_min_1, {
    updateSliderInput(session, "range_1",
                      value = c(input$range_min_1,max(input$range_1)))
  })

  observeEvent(input$range_max_1, {
    updateSliderInput(session, "range_1",
                      value = c(input$range_min_1,input$range_max_1))
  })

  observeEvent(input$range_1 , {

    updateNumericInput(session, "range_min_1", value = min(input$range_1))
    updateNumericInput(session, "range_max_1", value = max(input$range_1))

  }, priority = 200)
  
  # Output Handler:downloads file
  output$down_1<-downloadHandler(
    filename <- function(){
      paste0("ArchRBrowser_rds-",input$gene_name_1,input$plot_choice_download_peak2GeneLink)
    },
    content = function(file){
    
      if(input$plot_choice_download_peak2GeneLink==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$plot_width_1, height = input$plot_height_1)}
      
      else if(input$plot_choice_download_peak2GeneLink==".png")
      {png(file = file, width = input$plot_width_1, height = input$plot_height_1,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$plot_width_1, height = input$plot_height_1,units="in",res=1000)}
      
      
      p_co_access_peaks <- plotBrowserTrack(
        ArchRProj = savedArchRProject3,
        groupBy = "Clusters2",
        baseSize = 11,
        facetbaseSize = 11,
        geneSymbol = isolate(input$gene_name_1),
        upstream =-min(isolate(input$range_1))*1000 ,
        downstream = max(isolate(input$range_1))*1000,
        tileSize = isolate(input$tile_size_1),
        ylim =  c(0, isolate(input$ymax_1)),
        loops = getPeak2GeneLinks(savedArchRProject3)
      )[[input$gene_name_1]]
      
      grid.arrange(p_co_access_peaks)  # <- Orignal
      dev.off()
    }

  )
  output$co_access_peaks <- DT::renderDT(NULL)
  
  restartFN_2 <- observeEvent(input$restartButton_1,{
    if (isolate(input$gene_name_1) == ""){

      output$co_access_peaks <- renderPlot({
        p <- ggplot() +
          xlim(c(-5,5)) + ylim(c(-5,5)) +
          geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
        print(p)
      })
    }else{
      
   
  # Plot from peak2Genelinks and coaccessibility plots
      output$co_access_peaks <- renderPlot({
        grid::grid.newpage()
        p_co_access_peaks <- plotBrowserTrack(
          ArchRProj = savedArchRProject3,
          groupBy = "Clusters2",
          baseSize = 11,
          facetbaseSize = 11,
          geneSymbol = isolate(input$gene_name_1),
          upstream =-min(isolate(input$range_1))*1000 ,
          downstream = max(isolate(input$range_1))*1000,
          tileSize = isolate(input$tile_size_1),
          ylim =  c(0, isolate(input$ymax_1)),
          loops = getPeak2GeneLinks(savedArchRProject3)
         )[[input$gene_name_1]]

        grid.arrange(p_co_access_peaks)
      },height = 1200)

    }
  })
  
###########################################################
#Feature comparison : plot UMAP2                         ##
########################################################### 
  #Obtaining pwm matrix
  motif_PWMatrix=savedArchRProject3@peakAnnotation@listData[["Motif"]][["motifs"]]
  motif_ProbMatrices <- lapply(motif_PWMatrix, PWMatrixToProbMatrix)
  motifSummary=readRDS(savedArchRProject3@peakAnnotation@listData[["Motif"]][["Positions"]])
  
  
  
  #Provide gene/motif names for drop down
  observe({
    if(isolate(input$matrix_forComparison)=="MotifMatrix")
    {
      updateSelectizeInput(session, 'gene_forComparison_1', label = 'Feature Name 1',
                           choices = sort(motifMatrix_dropdown), 
                           server = TRUE,selected =sort(motifMatrix_dropdown)[1])
      
      updateSelectizeInput(session, 'gene_forComparison_2', label = 'Feature Name 2',
                           choices = sort(motifMatrix_dropdown), 
                           server = TRUE,selected =sort(motifMatrix_dropdown)[2])
    }
    else{
      updateSelectizeInput(session, 'gene_forComparison_1', label = 'Feature Name 1',
                           choices = sort(gene_names), 
                           server = TRUE,sort(gene_names)[1])
      updateSelectizeInput(session, 'gene_forComparison_2', label = 'Feature Name 2',
                           choices = sort(gene_names), 
                           server = TRUE,sort(gene_names)[2])
    }
  })
  
  #change it with the dropdown option for matrix
  observeEvent(input$matrix_forComparison,{
    
    if(isolate(input$matrix_forComparison)=="MotifMatrix")
    {
      updateSelectizeInput(session, 'gene_forComparison_1',label = 'Feature Name 1',
                           choices = sort(motifMatrix_dropdown), 
                           server = TRUE,selected =sort(motifMatrix_dropdown)[1])
      updateSelectizeInput(session, 'gene_forComparison_2',label = 'Feature Name 2',
                           choices = sort(motifMatrix_dropdown), 
                           server = TRUE,selected =sort(motifMatrix_dropdown)[2])
    }
    else{
      updateSelectizeInput(session, 'gene_forComparison_1', label = 'Feature Name 1',
                           choices = sort(gene_names), 
                           server = TRUE,selected =sort(gene_names)[1])
      updateSelectizeInput(session, 'gene_forComparison_2', label = 'Feature Name 2',
                           choices = sort(gene_names), 
                           server = TRUE,selected =sort(gene_names)[2])
    }
    
  })
  
  # Output Handler : download plots
  output$download_feature_comparison<-downloadHandler(
    filename <- function(){
      paste0("ArchRBrowser-",input$gene_forComparison_1,"-VS-",input$gene_forComparison_2,input$plot_choice_download_feature_comparison)
    },
    content = function(file){
      
      if(input$plot_choice_download_feature_comparison==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$gene_Comparison_plot_width, height = input$gene_Comparison_plot_height)}
      
      else if(input$plot_choice_download_feature_comparison==".png")
      {png(file = file, width = input$gene_Comparison_plot_width, height = input$gene_Comparison_plot_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$gene_Comparison_plot_width, height = input$gene_Comparison_plot_height,units="in",res=1000)}
      
      if(isolate(input$matrix_forComparison)=="GeneScoreMatrix")
        
      {
        
        gene1_plot=plotEmbedding(
          ArchRProj = savedArchRProject3,
          colorBy = "GeneScoreMatrix",
          #continuousSet = "yellowBlue",
          name = isolate(input$gene_forComparison_1),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject3),
          
        )
        
        gene2_plot=plotEmbedding(
          ArchRProj = savedArchRProject3,
          colorBy = "GeneScoreMatrix",
          #continuousSet = "yellowBlue",
          name = isolate(input$gene_forComparison_2),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject3),
          
        )}
      
      else if(isolate(input$matrix_forComparison)=="GeneIntegrationMatrix")
      {
        gene1_plot=plotEmbedding(
          ArchRProj = savedArchRProject3,
          colorBy = "GeneIntegrationMatrix",
          name = isolate(input$gene_forComparison_1),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject3),
          
        )
        
        gene2_plot=plotEmbedding(
          ArchRProj = savedArchRProject3,
          colorBy = "GeneIntegrationMatrix",
          name = isolate(input$gene_forComparison_2),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject3),
          
        )}
      
      else
      {
        gene1_plot=plotEmbedding(
          ArchRProj = savedArchRProject3,
          colorBy = "MotifMatrix",
          name = getFeatures(savedArchRProject3, 
                             select = paste(isolate(input$gene_forComparison_1), collapse="|"), 
                             useMatrix = "MotifMatrix"),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject3))[[2]]  #get z-score and deviation plot
        
        gene2_plot=plotEmbedding(
          ArchRProj = savedArchRProject3,
          colorBy = "MotifMatrix",
          name = getFeatures(savedArchRProject3, 
                             select = paste(isolate(input$gene_forComparison_2), collapse="|"), 
                             useMatrix = "MotifMatrix"),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject3))[[2]] #get z-score and deviation plot
        
        #get seq logo
        motif1=unlist(strsplit(getFeatures(savedArchRProject3, 
                                           select = paste(isolate(input$gene_forComparison_1), collapse="|"), 
                                           useMatrix = "MotifMatrix")[1],":"))[2]
        motif1_seqlogo=ggseqlogo(motif_ProbMatrices[[motif1]],method = 'prob')+ggtitle(motif1)
        
        motif2=unlist(strsplit(getFeatures(savedArchRProject3, 
                                           select = paste(isolate(input$gene_forComparison_2), collapse="|"), 
                                           useMatrix = "MotifMatrix")[1],":"))[2]
        motif2_seqlogo=ggseqlogo(motif_ProbMatrices[[motif2]],method = 'prob')+ggtitle(motif2)
        
        
        gene1_plot=grid.arrange(gene1_plot,motif1_seqlogo, ncol=1,nrow=2)
        gene2_plot=grid.arrange(gene2_plot,motif2_seqlogo, ncol=1,nrow=2)
      }
      
      grid.arrange(gene1_plot,gene2_plot, ncol=2) # <- Orignal
      
      
      dev.off()
    }
  )
  
  # Output Handler : download motif position
  output$download_motifPos<-downloadHandler(
    filename <- function(){
      paste0("MotifPosition-",input$motif_for_motifPos,".csv")
    },
    content = function(file){
      motif_name_Temp=unlist(strsplit(getFeatures(savedArchRProject3, 
                                                  select = paste(isolate(input$motif_for_motifPos), collapse="|"), 
                                                  useMatrix = "MotifMatrix")[1],":"))[2]
      
      temp=as.data.frame(motifSummary[motif_name_Temp])
      write.csv(temp, file)
      
    }
  )
  
  #get height and width based on motif matrix
  getHeight_featComparison<-function()
  {
    if(isolate(input$matrix_forComparison)=="MotifMatrix")
    {return(800)}
    else{"auto"}
  }
  
  
  output$feature_comparison <- DT::renderDT(NULL)
  restartFN <- observeEvent(input$gene_to_gene_restartButton, {
    if (isolate(input$gene_forComparison_1) == "" || 
        isolate(input$gene_forComparison_2) == ""){
      
      output$feature_comparison <- renderPlot({
        p <- ggplot() +
          xlim(c(-5,5)) + ylim(c(-5,5)) +
          geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
        print(p)
      })
    }else{
      # Plot feature comparison
      
      output$feature_comparison<- renderPlot({
        grid::grid.newpage()
        
        if(isolate(input$matrix_forComparison)=="GeneScoreMatrix")
          
        {
          
          gene1_plot=plotEmbedding(
            ArchRProj = savedArchRProject3,
            colorBy = "GeneScoreMatrix",
            #continuousSet = "SolarExtra",
            name = isolate(input$gene_forComparison_1),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject3),
            
          )
          
          gene2_plot=plotEmbedding(
            ArchRProj = savedArchRProject3,
            colorBy = "GeneScoreMatrix",
            #continuousSet = "SolarExtra",
            name = isolate(input$gene_forComparison_2),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject3),
            
          )}
        
        else if(isolate(input$matrix_forComparison)=="GeneIntegrationMatrix")
        {
          gene1_plot=plotEmbedding(
            ArchRProj = savedArchRProject3,
            colorBy = "GeneIntegrationMatrix",
            name = isolate(input$gene_forComparison_1),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject3),
            
          )
          
          gene2_plot=plotEmbedding(
            ArchRProj = savedArchRProject3,
            colorBy = "GeneIntegrationMatrix",
            name = isolate(input$gene_forComparison_2),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject3),
            
          )}
        
        else
        {
          gene1_plot=plotEmbedding(
            ArchRProj = savedArchRProject3,
            colorBy = "MotifMatrix",
            name = getFeatures(savedArchRProject3, 
                               select = paste(isolate(input$gene_forComparison_1), collapse="|"), 
                               useMatrix = "MotifMatrix"),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject3))[[2]]  #get z-score and deviation plot
          
          gene2_plot=plotEmbedding(
            ArchRProj = savedArchRProject3,
            colorBy = "MotifMatrix",
            name = getFeatures(savedArchRProject3, 
                               select = paste(isolate(input$gene_forComparison_2), collapse="|"), 
                               useMatrix = "MotifMatrix"),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject3))[[2]] #get z-score and deviation plot
          
          #get seq logo
          motif1=unlist(strsplit(getFeatures(savedArchRProject3, 
                                      select = paste(isolate(input$gene_forComparison_1), collapse="|"), 
                                      useMatrix = "MotifMatrix")[1],":"))[2]
          motif1_seqlogo=ggseqlogo(motif_ProbMatrices[[motif1]],method = 'prob')+ggtitle(motif1)
          
          motif2=unlist(strsplit(getFeatures(savedArchRProject3, 
                                             select = paste(isolate(input$gene_forComparison_2), collapse="|"), 
                                             useMatrix = "MotifMatrix")[1],":"))[2]
          motif2_seqlogo=ggseqlogo(motif_ProbMatrices[[motif2]],method = 'prob')+ggtitle(motif2)
          
          
          gene1_plot=grid.arrange(gene1_plot,motif1_seqlogo, ncol=1,nrow=2)
          gene2_plot=grid.arrange(gene2_plot,motif2_seqlogo, ncol=1,nrow=2)
          
          
          
        }
        
        grid.arrange(gene1_plot,gene2_plot, ncol=2)
        
      },height=getHeight_featComparison())
      
    }
  })

###########################################################
#scATACseq Motif footprinting                            ##
###########################################################
##Motif footprinting can only plotted on local deployment of shinyapp  if and if the analysis is performed on the same system
##It is due to the fact that Group coverage step of ArchR uses absolute paths instead of relative paths:please see https://github.com/GreenleafLab/ArchR/issues/529 for more details. Please also uncomment motif_Footprinting_panel in section Ui fluid page of file ui.R
# Output Handler
    output$motif_down_1 <- downloadHandler(
    filename <- function(){
      paste0("motif-footPrint",input$motifName_input,input$plot_choice_download_motif_down_1)
    },
    content = function(file){

      seFoot <- getFootprints(
        ArchRProj = savedArchRProject3,
        positions = motifPositions[input$motifName_input],
        groupBy = "Clusters2"
      )
      
      if(input$plot_choice_download_motif_down_1==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$motif_plot_width_1, height = input$motif_plot_height_1)}
      
      else if(input$plot_choice_download_motif_down_1==".png")
      {png(file = file, width = input$motif_plot_width_1, height = input$motif_plot_height_1,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$motif_plot_width_1, height = input$motif_plot_height_1,units="in",res=1000)}
      
      motif_FP <- plotFootprints(
        seFoot = seFoot,
        ArchRProj = savedArchRProject3,
        normMethod = input$normMethod_Input,
        plotName = "Footprints-No-Normalization",
        addDOC = TRUE,
        plot=FALSE,
        smoothWindow = 5,
        baseSize = 8
        )

      # grid.draw(motif_FP[[names(motif_FP)]]) # <- Orignal
      grid.arrange(motif_FP[[names(motif_FP)]]) # <- Modified
      dev.off() # <-Modified
    }

  )
  output$motifPlot <- DT::renderDT(NULL)

  restartFN_motifPlot <- observeEvent(input$plotFootprints_button,{
    if (isolate(input$motifName_input) == ""){

      output$motifPlot <- renderPlot({
        p <- ggplot() +
          xlim(c(-5,5)) + ylim(c(-5,5)) +
          geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
        print(p)
      })
    }else{

  # Plot motif foot printing
      output$motifPlot<- renderPlot({

        seFoot <- getFootprints(
          ArchRProj = savedArchRProject3,
          positions = motifPositions[input$motifName_input],
          groupBy = "Clusters2"
        )

        grid::grid.newpage()
        motif_FP <- plotFootprints(
          seFoot = seFoot,
          ArchRProj = savedArchRProject3,
          normMethod = input$normMethod_Input,
          plotName = "Footprints-No-Normalization",
          addDOC = TRUE,
          plot=FALSE,
          smoothWindow = 5,
          baseSize = 8
        )

        grid.draw(motif_FP[[names(motif_FP)]])
      })
    }
  }
)

###########################################################
#               HEATMAP trajectories                     ##
########################################################### 

  output$down_heatmap_traj1 <- downloadHandler(
    filename <- function(){
      paste("HEATMAP-trajectories-",input$matrix_forTraj1,input$plot_choice_down_heatmap_traj1,sep="")
    },
    content <- function(file){
      
      
      if(input$plot_choice_down_heatmap_traj1==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$matrix_forTraj_width, height = input$matrix_forTraj_height)}
      
      else if(input$plot_choice_down_heatmap_traj1==".png")
      {png(file = file, width = input$matrix_forTraj_width, height = input$matrix_forTraj_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$matrix_forTraj_width, height = input$matrix_forTraj_height,units="in",res=1000)}
      
      #print(p_peakMatrix_traj@ht_list$PeakMatrix+p_peakMatrix_traj@ht_list$heatmap_annotation_18+p_peakMatrix_traj@layout)
      
      if (input$matrix_forTraj1 == "GeneScoreMatrix")
      {print(p_GeneScore_traj@ht_list[[1]]+p_GeneScore_traj@ht_list[[2]]+p_GeneScore_traj@layout)}
      
      else if (input$matrix_forTraj1 == "GeneIntegrationMatrix")
      {print(p_GIM_traj@ht_list[[1]]+p_GIM_traj@ht_list[[2]]+p_GIM_traj@layout)}
      
      else if (input$matrix_forTraj1 == "MotifMatrix")
      {print(p_motifMatrix_traj@ht_list[[1]]+p_motifMatrix_traj@ht_list[[2]]+p_motifMatrix_traj@layout)}
      
      else
      {print(p_peakMatrix_traj@ht_list[[1]]+p_peakMatrix_traj@ht_list[[2]]+p_peakMatrix_traj@layout)}
      
      dev.off()
      
    })
  
  output$down_heatmap_traj2 <- downloadHandler(
    filename <- function(){
      paste("HEATMAP-trajectories-",input$matrix_forTraj2,input$plot_choice_down_heatmap_traj2,sep="")
    },
    content <- function(file){
      
      
      if(input$plot_choice_down_heatmap_traj2==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$matrix_forTraj_width2, height = input$matrix_forTraj_height2)}
      
      else if(input$plot_choice_down_heatmap_traj2==".png")
      {png(file = file, width = input$matrix_forTraj_width2, height = input$matrix_forTraj_height2,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$matrix_forTraj_width2, height = input$matrix_forTraj_height2,units="in",res=1000)}
      
      #print(p_peakMatrix_traj@ht_list$PeakMatrix+p_peakMatrix_traj@ht_list$heatmap_annotation_18+p_peakMatrix_traj@layout)
      
      if (input$matrix_forTraj2 == "GeneScoreMatrix")
      {print(p_GeneScore_traj@ht_list[[1]]+p_GeneScore_traj@ht_list[[2]]+p_GeneScore_traj@layout)}
      
      else if (input$matrix_forTraj2 == "GeneIntegrationMatrix")
      {print(p_GIM_traj@ht_list[[1]]+p_GIM_traj@ht_list[[2]]+p_GIM_traj@layout)}
      
      else if (input$matrix_forTraj2 == "MotifMatrix")
      {print(p_motifMatrix_traj@ht_list[[1]]+p_motifMatrix_traj@ht_list[[2]]+p_motifMatrix_traj@layout)}
      
      else
      {print(p_peakMatrix_traj@ht_list[[1]]+p_peakMatrix_traj@ht_list[[2]]+p_peakMatrix_traj@layout)}
      
      dev.off()
      
    })
  
  
  output$traj_heatmap1 <- DT::renderDT(NULL)
  output$traj_heatmap2 <- DT::renderDT(NULL)
  restartFN <- observeEvent(input$plotheat_traj_button, {
    
  output$traj_heatmap1 <- renderPlot({
   
    if (input$matrix_forTraj1 == "GeneScoreMatrix")
    {p_GeneScore_traj}
    
    else if (input$matrix_forTraj1 == "GeneIntegrationMatrix")
    {p_GIM_traj}
    
    else if (input$matrix_forTraj1 == "MotifMatrix")
    {p_motifMatrix_traj}
    
    else
    {p_peakMatrix_traj}
  
    },height = 600,width = 500)
  
  output$traj_heatmap2 <- renderPlot({
    
    if (input$matrix_forTraj2 == "GeneScoreMatrix")
    {p_GeneScore_traj}
    
    else if (input$matrix_forTraj2 == "GeneIntegrationMatrix")
    {p_GIM_traj}
    
    else if (input$matrix_forTraj2 == "MotifMatrix")
    {p_motifMatrix_traj}
    
    else
    {p_peakMatrix_traj}
    },height = 600,width=500)

  })
  ###########################################################
  #         HEATMAP _INTEGRATED_peakto gene link plot      ##
  ###########################################################   
  
  # Plot of Heatmap of Peak To Gene Links in section 15.3.2
  
  output$heatmap_p2g <- downloadHandler(
    filename <- function(){
      paste0("HEATMAP-peak2genelinks",input$plot_choice_down_heatmap_p2g)
    },
    content=function(file) {
      
      if(input$plot_choice_down_heatmap_p2g==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$p2g_plot_width_1, height = input$p2g_plot_height_1)}
      
      else if(input$plot_choice_down_heatmap_p2g==".png")
      {png(file = file, width = input$p2g_plot_width_1, height = input$p2g_plot_height_1,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$p2g_plot_width_1, height = input$p2g_plot_height_1,units="in",res=1000)}
      
      print(p_heatmap_peak_to_gene@ht_list[[1]]+p_heatmap_peak_to_gene@ht_list[[2]]+p_heatmap_peak_to_gene@layout)
      dev.off()
    }
  )
  output$heatmap_peak_to_gene <- DT::renderDT(NULL)
  restartFN <- observeEvent(input$heatmap_peak2gl, {
  output$heatmap_peak_to_gene <- renderPlot({p_heatmap_peak_to_gene},height = 800,width=1000)
  })
 }
###########################################################
#             END OF FILE                                ##
########################################################### 
