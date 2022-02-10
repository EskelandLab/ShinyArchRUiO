# Load libraries so they are available
# Run the app through this file.
source("ui.R")
source("server.R")
shinyApp(ui:ui, server:shinyServer)
