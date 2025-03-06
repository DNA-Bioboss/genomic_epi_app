server <- function(input, output, session){
  
  # #security and password
  # result_auth <- secure_server(check_credentials = check_credentials(credentials))
  # output$res_auth <- renderPrint({
  #   reactiveValuesToList(result_auth)
  # })
  # Create a reactive value to track button clicks
  plot_ready <- reactiveValues(trigger = FALSE)
  
  
  # Create a reactive value to store entered person IDs and colors
  reactive_ids <- reactiveVal(NULL)
  reactive_colors <-reactiveVal(NULL)
  
  # Observe the action button to change the state
  observeEvent(input$plot_button, {
    plot_ready$trigger <- TRUE
    #split entered values
    if(length(input$person_ids)>0){
      entered_ids <- strsplit(input$person_ids, ",\\s*")[[1]]
      # Update person_id reactive value
      reactive_ids(entered_ids)
    }
  })
  
  output$gnattchart <-renderPlot({
    if (!plot_ready$trigger) {
      isolate({
        
        # create plot
        plot_x <-ggplot(df, aes(x=date,
                                y=person_id,
                                fill=person_facility_name))+
          geom_tile()+
          # scale_fill_manual(values = reduced_facility_colors)+
          scale_fill_manual(values =location_colors)+
          scale_x_date(date_breaks = "1 months",
                       labels=label_date_short()
          )+
          geom_point(data = all_DODs, size=1, aes(x = date, y = person_id), color = "Red")+
          labs(title="Initial Plot: Timeline of Diagnoses by Facility and Person ID",
               subtitle=paste0("Created on ",Sys.Date()), 
               x="Date",y="")+
          theme_Publication()
        
        
        # calculate the dimensions of the phylogenetic tree to align with chart
        bottom_value_tree <-tree_reduced %>% as_tibble() %>% na.omit()
        bottom_value <-((NROW(unique(bottom_value_tree$label))/NROW(unique(df$person_id)))*100)
        
        # observe({
        # Pause execution here and open the browser to inspect
        # 8621346
        #   browser()
        # 
        #   print(bottom_value)  # Inspect value
        # })
        
        # layout dimensions
        layout <- c(
          area(t = 0, l = 0, b = bottom_value, r = 33),
          area(t = 0, l = 12, b = 100, r = 100)
        )
        # render tree and chart
        updated_tree_viz + plot_x+
          plot_layout(design = layout)
        
      })# END of isolate
    } else {
      
      #if additional person_ids are added
      if(input$person_ids !=""){
        # Retrieve the stored entered_ids
        entered_ids <- reactive_ids()
        # find associated facilitates (aka links) for each sample
        v <-as.vector(entered_ids)
        additonal_links <-determine_links(v)
        
        # observe({
        #   # Pause execution here and open the browser to inspect
        #   #8621346
        #   browser()
        # 
        #   print(entered_ids)  # Inspect value
        # })
        
        # add in-between dates
        additonal_links_updated <-adding_dates(additonal_links)
        
        # capture additional ids
        additional_ids <-additonal_links_updated %>% 
          select(person_id) %>% distinct() %>% pull()
        
        # order person_ids
        viz_order_updated <-c(additional_ids,tree_viz_label_order)
        
        # bind and factor person_ids
        df <-rbind(additonal_links_updated,df) %>% 
          mutate(person_id=factor(person_id,levels=viz_order_updated))
        
        #Adding retrieving DODs for additional samples
        additonal_DODs <-determine_DOD(v)
        prepped_addi_DODs <-prep_DODs(additonal_DODs)
        all_DODs <-rbind(prepped_addi_DODs,all_DODs)
        
        #update color palette
        location_colors <-set_color_palette(df)
        ################## tab 2 ###############
        
      }# end of IF statement
      
      #create plot
      plot_x <-ggplot(df, aes(x=date,
                              y=person_id,
                              fill=person_facility_name))+
        geom_tile()+
        scale_fill_manual(values =location_colors)+
        # scale_fill_manual(values = reduced_facility_colors)+
        scale_x_date(date_breaks = "1 months",
                     labels=label_date_short(),
                     limits = as.Date(c(input$dateRange[1],input$dateRange[2]))
        )+
        geom_point(data = all_DODs, size=1, aes(x = date, y = person_id), color = "Red")+
        labs(title="Updated Plot: Timeline of Diagnoses by Facility and Person ID",
             subtitle=paste0("From ",input$dateRange[1]," to ",input$dateRange[2]), 
             x="Date",y="")+
        theme_Publication()
      
      # calculate the dimensions of the phylogenetic tree to align with chart
      bottom_value_tree <-tree_reduced %>% as_tibble() %>% na.omit()
      if(input$person_ids !=""){
        bottom_value <-((NROW(unique(bottom_value_tree$label))/NROW(unique(df$person_id)))*100)
        bottom_value <-ceiling(bottom_value)+2
      } else{
        bottom_value <-((NROW(unique(bottom_value_tree$label))/NROW(unique(df$person_id)))*100)
        
      }
      # layout dimensions
      layout <- c(
        area(t = 0, l = 0, b = bottom_value, r = 33),
        area(t = 0, l = 11, b = 100, r = 100)
      )
      # render tree and chart
      updated_tree_viz + plot_x+
        plot_layout(design = layout)
    } # END of else
  })# end of renderplot
  ############################## START OF TAB 2 #################################
  output$network <-renderVisNetwork({
  # output$network <-renderUI({
    if (!plot_ready$trigger) {
      isolate({
        network_ren <-vis_epicontactsMT(
          epic_object,
          thin = F,
          node_color = "type",
          # label = "ID",
          # annot = TRUE,
          node_shape = "node_type",
          shapes = c(case = "circle", facility = "plus"),
          # edge_label = NULL,
          edge_color = "location",
          legend = TRUE,
          # col_pal = cases_pal,
          # NA_col = "lightgrey",
          # edge_col_pal = edges_pal,
          # width = "90%",
          # width = "100%",
          width = "1300px",
          # height = "100%",
          height = "700px",
          # height = "800px",
          selector = TRUE,
          editor = FALSE,
          edge_width = 3,
          physics=T
        ) %>% 
          visEdges(physics=F, smooth = list(enabled=F,type="straightCross"))
        
        network_ren
      }) # end of isolate
    }else{
      if(input$person_ids !=""){
        # Retrieve the stored entered_ids
        entered_ids <- reactive_ids()
        # find associated facilitates (aka links) for each sample
        v <-as.vector(entered_ids)
        additonal_links <-determine_links(v)
        
        #Adding retrieving DODs for additional samples
        additonal_DODs <-determine_DOD(v)
        prepped_addi_DODs <-prep_DODs(additonal_DODs)
        all_DODs <-rbind(prepped_addi_DODs,all_DODs)
        
        
        addi_link_types <-determine_link_type(additonal_links,prepped_addi_DODs)
        updated_link_status <-rbind(addi_link_types,link_status)
        addi_AMR <-determine_AMR_genes(epitrax_data,v)
        updated_pcr_list <-rbind(addi_AMR,reduced_pcr_combined)
        updated_only_links <-rbind(additonal_links,links_reduced)
        updated_only_DODs <-rbind(additonal_DODs, DODs_reduced)
        updated_nodelist <-create_nodelist(updated_only_DODs,updated_pcr_list,updated_only_links)
        
        updated_epic_object <-generate_epicontacts_object(updated_link_status,updated_nodelist)
        epic_object <-updated_epic_object
      } # end of if statement
      #create network
      network_ren <-vis_epicontactsMT(
        epic_object,
        thin = F,
        node_color = "type",
        # label = "ID",
        # annot = TRUE,
        node_shape = "node_type",
        shapes = c(case = "circle", facility = "plus"),
        # edge_label = NULL,
        edge_color = "location",
        legend = TRUE,
        # col_pal = cases_pal,
        # NA_col = "lightgrey",
        # edge_col_pal = edges_pal,
        # width = "90%",
        # width = "100%",
        # width = "900px",
        width = "1300px",
        # width = "1000px",
        # height = "100%",
        height = "700px",
        # height = "700px",
        selector = TRUE,
        editor = FALSE,
        edge_width = 3,
        physics=T
      ) %>%
        visEdges(physics=F, smooth = list(enabled=F,type="straightCross"))
      
      #render network
      network_ren
      
    }# end of else statement
  }) # end of render plot
  
} #end of server
