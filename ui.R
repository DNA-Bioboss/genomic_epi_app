########################## USER INTERFACE ###################################
ui <- 
  # secure_app(head_auth = tags$script(inactivity),
                 fluidPage(  theme = shinytheme("cerulean"),
                navbarPage("Bioinformatic Tools",
                           tabPanel("Gantt Chart",
                                    sidebarPanel(
                                      #Date Range UI
                                      dateRangeInput('dateRange',
                                                     label = 'Date range input: mm-dd-yyyy',
                                                     start = min(df$date, na.rm = TRUE),
                                                     end = max(df$date, na.rm = TRUE),
                                                     min = min(df$date, na.rm = TRUE),
                                                     max = max(df$date, na.rm = TRUE),
                                                     format = "mm-dd-yyyy"),
                                      # Text Input for Multiple person_ids
                                      textInput("person_ids", "Enter person IDs (comma separated):",
                                                value = ""),
                                      actionButton("plot_button", "Plot", class = "btn btn-primary")
                                      
                                    ),#sidebarpanel
                                    # ), #tabpanel
                                    
                                    
                                    mainPanel(
                                      plotOutput(outputId = "gnattchart",height = "600px")
                                    )
                           ), #end of tabpanel #1 (if turning on tab 2 add comma)
                           ############ UI for tab 2 ######################################################
                           tabPanel("Network",
                                    #     sidebarPanel(),   
                                    mainPanel(
                                      visNetworkOutput(outputId = "network",height = "700px", width = "1300px")
                                    ) # end of mainpanel
                           ) # end of tab 2
                )# navbarpage
) # fluidPage
# )