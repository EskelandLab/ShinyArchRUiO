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
      paste0("UMAP-",input$UMAP1_forComparison,".pdf")
    },
    content = function(file){
      pdf(file = file,onefile=FALSE, width = input$UMAP1_plot_width, height = input$UMAP1_plot_height)
      # grid.arrange(umaps[[input$UMAP1_forComparison]]) # <-Original
      grid.arrange(umaps[[input$UMAP1_forComparison]])  # <- Modification
      dev.off()
    }
  )

  output$download_UMAP2<-downloadHandler(
    filename <- function(){
      paste0("UMAP-",input$UMAP2_forComparison,".pdf")
    },
    content = function(file){
      pdf(file = file,onefile=FALSE, width = input$UMAP2_plot_width, height = input$UMAP2_plot_height)
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
      paste0("ArchRBrowser-",input$gene_name,".pdf")
    },
    content = function(file){
      pdf(file = file,onefile=FALSE, width = input$plot_width, height = input$plot_height)
      p_browser_atacClusters<- plotBrowserTrack(
        ArchRProj = savedArchRProject2,
        groupBy = isolate(input$group_by),
        baseSize = 11,
        facetbaseSize = 11,
        geneSymbol = isolate(input$gene_name),
        upstream = -min(isolate(input$range))*1000,
        downstream = max(isolate(input$range))*1000,
        tileSize = isolate(input$tile_size),
        ylim =  c(0, isolate(input$ymax))
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
          ArchRProj = savedArchRProject2,
          baseSize = 11,
          facetbaseSize = 11,
          groupBy = isolate(input$group_by),
          geneSymbol = isolate(input$gene_name),
          upstream = -min(isolate(input$range))*1000,
          downstream = max(isolate(input$range))*1000,
          tileSize = isolate(input$tile_size),
          ylim =  c(0, isolate(input$ymax)),
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
      paste0("ArchRBrowser_rds-",input$gene_name_1,".pdf")
    },
    content = function(file){
    
      pdf(file = file, onefile=FALSE,width = input$plot_width_1, height = input$plot_height_1)
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
        loops = getCoAccessibility(savedArchRProject3)
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
          loops = getCoAccessibility(savedArchRProject3)
         )[[input$gene_name_1]]

        grid.arrange(p_co_access_peaks)
      },height = 1200)

    }
  })
  
###########################################################
#Feature comparison : plot UMAP2                         ##
########################################################### 
  
  # Output Handler : download plots
  output$download_feature_comparison<-downloadHandler(
    filename <- function(){
      paste0("ArchRBrowser-",input$gene_forComparison_1,"-VS-",input$gene_forComparison_2,".pdf")
    },
    content = function(file){
      
      
      pdf(file = file,onefile=FALSE, width = input$gene_Comparison_plot_width, height = input$gene_Comparison_plot_height)
      if(isolate(input$matrix_forComparison)=="GeneScoreMatrix")
        
      {
        
        gene1_plot=plotEmbedding(
          ArchRProj = savedArchRProject1,
          colorBy = "GeneScoreMatrix",
          name = isolate(input$gene_forComparison_1),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject1),
          
        )
        
        gene2_plot=plotEmbedding(
          ArchRProj = savedArchRProject1,
          colorBy = "GeneScoreMatrix",
          name = isolate(input$gene_forComparison_2),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject1),
          
        )}
      else
      {
        gene1_plot=plotEmbedding(
          ArchRProj = savedArchRProject2,
          colorBy = "GeneIntegrationMatrix",
          name = isolate(input$gene_forComparison_1),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject1),
          
        )
        
        gene2_plot=plotEmbedding(
          ArchRProj = savedArchRProject2,
          colorBy = "GeneIntegrationMatrix",
          name = isolate(input$gene_forComparison_2),
          embedding = "UMAP",
          imputeWeights = getImputeWeights(savedArchRProject1),
          
        )
      }
      
      grid.arrange(gene1_plot,gene2_plot, ncol=2) # <- Orignal
      
      
      dev.off()
    }
  )
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
            ArchRProj = savedArchRProject1,
            colorBy = "GeneScoreMatrix",
            name = isolate(input$gene_forComparison_1),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject1),
            
          )
          
          gene2_plot=plotEmbedding(
            ArchRProj = savedArchRProject1,
            colorBy = "GeneScoreMatrix",
            name = isolate(input$gene_forComparison_2),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject1),
            
          )}
        else
        {
          gene1_plot=plotEmbedding(
            ArchRProj = savedArchRProject2,
            colorBy = "GeneIntegrationMatrix",
            name = isolate(input$gene_forComparison_1),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject1),
            
          )
          
          gene2_plot=plotEmbedding(
            ArchRProj = savedArchRProject2,
            colorBy = "GeneIntegrationMatrix",
            name = isolate(input$gene_forComparison_2),
            embedding = "UMAP",
            imputeWeights = getImputeWeights(savedArchRProject1),
            
          )
        }
        
        grid.arrange(gene1_plot,gene2_plot, ncol=2)
      })
      
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
      paste0("motif-footPrint",input$motifName_input,".pdf")
    },
    content = function(file){

      seFoot <- getFootprints(
        ArchRProj = savedArchRProject3,
        positions = motifPositions[input$motifName_input],
        groupBy = "Clusters2"
      )

      pdf(file = file,onefile=FALSE, width = input$motif_plot_width_1, height = input$motif_plot_height_1)
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
      paste("HEATMAP-trajectories-",input$matrix_forTraj1,".pdf",sep="")
    },
    content <- function(file){
      
      pdf(file = file,onefile=TRUE, width = input$matrix_forTraj_width, height = input$matrix_forTraj_height)
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
      paste("HEATMAP-trajectories-",input$matrix_forTraj2,".pdf",sep="")
    },
    content <- function(file){
      
      pdf(file = file,onefile=TRUE, width = input$matrix_forTraj_width2, height = input$matrix_forTraj_height2)
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
      paste0("HEATMAP-peak2genelinks",".pdf")
    },
    content=function(file) {
      pdf(file = file,onefile=TRUE, width = input$p2g_plot_width_1, height = input$p2g_plot_height_1)
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
