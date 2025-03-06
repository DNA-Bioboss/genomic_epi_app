# remotes::install_github("YuLab-SMU/ggtree")

library(shiny)
library(shinythemes)
library(readxl)
library(grDevices)
library(ggplot2)
library(scales)
library(tidyverse)
library(ggtree)
library(ape)
library(patchwork)
library(visNetwork)
library(epicontacts)
library(bslib)
library(shinymanager)

library(treeio)
library(BiocManager)


# Load data
# contains HSNs and NCBI accession numbers
SAMN_data<-read_excel("fake_SAMN_data.xlsx")

# NCBI SNP cluster tree: tree file retrieved from pathogen detect
nw_tree<-read.tree("export1.newick") 

# NCBI isolate data
sample_data <-read_csv("isolates.csv")

# epitrax information
epitrax_data<-read_excel("fake_epitrax_data.xlsx")
########################### PARAMETERS ########################################
target <-"SAMN38334673"


########################## security and password #############################
# credentials <- data.frame(
#   user = c("Deloitte"),
#   password = c("Solutions"),
#   stringsAsFactors = FALSE
# )

# inactivity <- "function idleTimer() {
# var t = setTimeout(logout, 120000);
# window.onmousemove = resetTimer; // catches mouse movements
# window.onmousedown = resetTimer; // catches mouse movements
# window.onclick = resetTimer;     // catches mouse clicks
# window.onscroll = resetTimer;    // catches scrolling
# window.onkeypress = resetTimer;  //catches keyboard actions
# function logout() {
# window.close();  //close the window
# }
# function resetTimer() {
# clearTimeout(t);
# t = setTimeout(logout, 120000);  // time is in milliseconds (1000 is 1 second)
# }
# }
# idleTimer();"


################## Phylogenetic Cleaning ######################################
# The tree is large and must be trimmed
# For this script's purpose, we will assume that our sample of interest is
# SAMN38334673
# We will have to identify what tips to drop, by looking at the internal nodes

# Convert the tree file into a tibble df
tibble_nw_tree <-nw_tree %>% as_tibble() 

# Write a function that identifies all internal nodes related to a tip
find_related_nodes <- function(data) {
  library(dplyr)
  
  # Initialize a list to store results
  results <- list()
  
  # Iterate over each row in the data
  for (i in seq_len(nrow(data))) {
    parent <- data$parent[i]
    node <- data$node[i]
    related_nodes <- c() # Store related nodes
    
    current_node <- parent
    while (TRUE) {
      # Find the parent of the current node
      next_row <- data %>% filter(node == current_node)
      
      if (nrow(next_row) == 0) break # No more parents to find
      
      current_node <- next_row$parent
      related_nodes <- c(related_nodes, current_node)
      
      # Break the loop if current node loops back to itself
      if (current_node == next_row$node) break
    }
    
    # Store results for this row
    results[[i]] <- c(parent, node, related_nodes)
  }
  
  # Convert results to a data frame
  max_length <- max(sapply(results, length))
  results <- do.call(rbind, lapply(results, function(x) {
    c(x, rep(NA, max_length - length(x)))
  }))
  
  colnames(results) <- c("parent", "node", paste0("G", seq_len(ncol(results) - 2)))
  return(as.data.frame(results))
}

#call the function
result <- find_related_nodes(tibble_nw_tree) 

# Find the internal nodes related to target sample (SAMN38334673)
target_internal_nodes <-result %>% 
  right_join(tibble_nw_tree) %>% 
  select(-branch.length) %>%
  filter(label==target) %>%
  mutate(parent=as.numeric(parent)) %>% 
  pivot_longer(cols=c(1,3:(NCOL(result)-2)),names_to="node_gen",values_to = "internal_node") %>% 
  filter(!is.na(internal_node)) %>% 
  select(internal_node) %>% 
  distinct() %>% 
  pull()

# Identify how many tips are in each internal node
tip_totals <-result %>% 
  pivot_longer(col=c(1,3:NCOL(result)),names_to="node_gen",values_to = "internal_node") %>% 
  filter(!is.na(internal_node)) %>% 
  mutate(is_target=internal_node %in% target_internal_nodes) %>% 
  filter(is_target==T) %>% 
  select(node,internal_node) %>%
  distinct() %>% 
  group_by(internal_node) %>% 
  summarise(n=n())

# Select the node that has 20-30 tip labels
index_node <-tip_totals %>% 
  filter(n > 20 & n < 30) %>% 
  select(internal_node) %>% 
  pull()

# if there are multiple index nodes the internal node
# with the most tip labels between 20-30 will be selected
if(NROW(index_node)>1){
  index_nodes <-tip_totals %>% 
    filter(n > 20 & n < 30)  
  max_row <- index_nodes[which.max(index_nodes$n), ] 
  index_node <-max_row$internal_node
} 

# obtain all tips associated with the index node
associated_tips <-result %>% 
  pivot_longer(col=c(1,3:NCOL(result)),names_to="node_gen",values_to = "internal_node") %>% 
  filter(!is.na(internal_node)) %>% 
  filter(internal_node==index_node) %>% 
  select(node) %>% 
  pull()

#obtain all tips not associated with the index node
tips_to_drop <-tibble_nw_tree %>% 
  mutate(associated_tip=node %in% associated_tips) %>%
  filter(associated_tip==F) %>%
  filter(!is.na(label)) %>% 
  select(label) %>% 
  pull()

# drop tips from phylo tree file
tree_reduced <-drop.tip(nw_tree,tips_to_drop)

############ CLEANING DATA FROM EPITRAX ######################################
# pull labels from reduced tree to find which of them are in Epitrax
reduced_SAMN_labels <-tree_reduced %>% 
  as_tibble() %>% 
  filter(!is.na(label)) %>% 
  select(label) %>%
  pull()

# Obtain epitrax data associated with labels from "tree_reduced"
# This next step would be done using epitrax data b/c the SAMN ID is yet to be stored in epitrax
# So this step will be done using the excelsheet that contains lab_accession_no (aka HSN) and SAMN
# this step will capture the lab_accession_no (HSN), ideally this step would be skipped once SAMN ID
# is entered in Epitrax
# 
SAMNs_intree <-SAMN_data %>% 
  mutate(in_reduced_tree=SAMN %in% reduced_SAMN_labels) %>% 
  filter(in_reduced_tree==T) %>% 
  select(HSN) %>% 
  pull()


# Identify person_ids associated with lab_accession_no (HSN). Ideally we would want to use SAMN ID to find person_ids
person_ids_reduced <-epitrax_data %>% 
  mutate(in_tree=lab_accession_no %in% SAMNs_intree) %>% 
  filter(in_tree==T) %>% 
  select(person_id) %>% 
  distinct() %>% 
  pull()

# additional samples
if(exists("x")){
  person_ids_reduced <-c(person_ids_reduced,x)
}
# Identify DOD (date of diagnosis) by selecting the earliest patient_date_diagnosed date associated for each person_id
determine_DOD <-function(vector){
  petal <-vector
  DODs_reduced <-epitrax_data %>% 
    mutate(person_id=as.character(person_id)) %>%
    mutate(in_tree=person_id %in% petal) %>% 
    filter(in_tree==T) %>% 
    select(person_id,patient_date_diagnosed) %>% 
    mutate(patient_date_diagnosed=ifelse(patient_date_diagnosed=="",NA_real_,patient_date_diagnosed)) %>% 
    distinct() %>% 
    arrange(patient_date_diagnosed) %>% 
    group_by(person_id) %>% 
    filter(patient_date_diagnosed==patient_date_diagnosed[1]) %>% 
    ungroup() 
  # mutate(person_id=as.numeric(person_id))
  
  # if a person_id has NA as DOD, then the NA will be substituted with current date  
  if(any(is.na(DODs_reduced$patient_date_diagnosed)==T)){
    DODs_reduced <-DODs_reduced(is.na(patient_date_diagnosed),Sys.Date(),patient_date_diagnosed)
  }
  return(DODs_reduced)
}
#call "determine_DOD" fxn
DODs_reduced <-determine_DOD(person_ids_reduced)

#Using the person IDs find all associated facilities
determine_links <-function(vector){
  flower <-vector
  links_reduced <-epitrax_data %>% 
    mutate(in_tree=person_id %in% flower) %>% 
    filter(in_tree==T) %>% 
    select(person_id,person_facility_name,patient_visit_start_date,patient_visit_end_date) %>% 
    distinct() %>% 
    filter(person_facility_name!="") %>% 
    separate(patient_visit_start_date,c("start_date", "start_time"), " ") %>% 
    select(-start_time) %>% 
    mutate(start_date=format(as.Date(start_date),"%Y-%m-%d")) %>% 
    separate(patient_visit_end_date,c("end_date", "end_time"), " ") %>% 
    select(-end_time) %>% 
    mutate(end_date=format(as.Date(end_date),"%Y-%m-%d")) %>% 
    mutate(end_date=ifelse(is.na(end_date),start_date,end_date)) %>% 
    filter(!is.na(start_date)) %>% 
    filter(!is.na(end_date)) %>% 
    mutate(start_date=ymd(start_date),
           end_date=ymd(end_date))
  
  return(links_reduced)
}

# call determine_links function
links_reduced <-determine_links(person_ids_reduced)

# Identify AMR genes associated with person IDs
determine_AMR_genes<-function(info_df,river_vector){
  reduced_AMR <-info_df %>% 
    mutate(in_tree=person_id %in% river_vector) %>% 
    filter(in_tree==T) %>% 
    mutate(PCR_testing=str_detect(lab_test_type,"PCR for")) %>% 
    filter(patient_event_type=="Carbapenemase-Producing Acinetobacter baumannii (CP-CRAB)" & PCR_testing==T & lab_test_result.1=="Positive / Reactive") %>%
    unique() %>% 
    select(person_id,lab_test_type) %>% 
    distinct() %>% 
    mutate(AMR_gene=str_extract(lab_test_type,"(?<=for\\s).+(?=\\sDNA)")) %>% 
    select(-lab_test_type) %>% 
    arrange(person_id,AMR_gene)
  
  # aggregate identified AMR genes with person IDs
  pcr_combined <-aggregate(AMR_gene~ person_id,reduced_AMR,paste, collapse = " & ") %>% 
    rename(c("type"="AMR_gene"))
  return(pcr_combined)
}
#call function
reduced_pcr_combined <-determine_AMR_genes(epitrax_data,person_ids_reduced)

# Create SAMN,person_ID chart
SAMN_pID_tbl <-epitrax_data %>% 
  mutate(in_tree=lab_accession_no %in% SAMNs_intree) %>% 
  filter(in_tree==T) %>% 
  select(person_id,lab_accession_no) %>% 
  distinct() %>% 
  left_join(SAMN_data,by=c("lab_accession_no"="HSN")) %>% 
  select(person_id,SAMN)

###################### CREATE PHYLO FILE ######################################
# update phylogenetic tree file to reflect which tips have associated person_ids
# reflect that by pasting "patient" in front of person_id
updated_tree <-tree_reduced %>% 
  as_tibble() %>% 
  left_join(SAMN_pID_tbl,by=c("label"="SAMN")) %>% 
  mutate(new_label=ifelse(is.na(person_id),label,paste0("Patient ",person_id))) %>% 
  select(parent,node,branch.length,new_label) %>% 
  rename(c("label"="new_label")) %>% 
  as.phylo(length="branch.length", label="label")



################## Phylogenetic and Gnatt Visualization #######################
# Create a phylogenetic and gantt chart visualization containing epidemiological data.
# This script is to be ran after the "cleaning file" script.
# The double display would be due to the use of patchwork and the elements align would 
# align. Addition of elements that do not contain genetic data but have investigation data
# regarding location would appear towards the end of the gantt chart.

#theme
theme_Publication <- function(base_size=11, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(#face = "bold",
      size = rel(1.1), hjust = -0.05),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1.0)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(size = rel(1.0)), 
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      #legend.margin = unit(0, "cm"),
      #legend.title = element_text(face="italic"),
      plot.margin=unit(c(2,2,2,2),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold"),
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, vjust = 0, size = rel(1))
      
    ))
  ################### ADDING DATES FOR START AND END DATES #######################  
}

adding_dates <-function(dataframe){
  dfx <-dataframe
  # Loop to create dates for each row
  all_dates_links <-data.frame()
  
  for(i in 1:NROW(dfx)){
    dfi1<-dfx[i,] %>% 
      pivot_longer(cols=c("start_date","end_date"),names_to = "date_type",values_to = "date")
    
    dfi_btwn_dates <-data.frame("person_id"=dfi1$person_id[1],
                                "person_facility_name"=dfi1$person_facility_name[1],
                                # "DOD"=dfi1$DOD[1],
                                "date_type"="btwn_date",
                                "date" = seq(min(dfi1$date), max(dfi1$date), by = "day")) 
    all_dates_links <-rbind(all_dates_links,dfi_btwn_dates)
    
  }
  # update person_ids with "patient" in front
  all_dates_links_updated <-all_dates_links %>% 
    mutate(person_id=paste0("Patient ",person_id)) 
  return(all_dates_links_updated)
}
#call adding_dates fxn
all_dates_links_updated <-adding_dates(links_reduced)

################## ADDING TREE SAMPLES TO CHART ################
#gaterh person Ids in data "all_dates_links_updated"
patient_ids <-all_dates_links_updated %>% select(person_id) %>% pull()


# Identify which samples are not in the chart
# Prep to match chart data by adding cols: ORIGINAL
labels_not_inchart_df <-updated_tree %>% 
  as_tibble() %>% 
  mutate(in_chart=label %in% patient_ids) %>% 
  filter(in_chart==F, !is.na(label)) %>% 
  select(label) %>% 
  rename(c("person_id"="label")) %>% 
  mutate(person_facility_name="DOD",
         date_type="DOD",
         date=NA) %>% 
  mutate(date=as.Date(date,format="%Y-%m-%d"))



# bind chart data with tree samples
all_labels <-rbind(all_dates_links_updated,labels_not_inchart_df)


##### order the chart as tree visualization

# visualize phylo tree
updated_tree_viz<-ggtree(updated_tree)+
  # geom_tiplab()+
  scale_x_continuous(limits = c(0, 20))

# capture order of tips of tree
tree_viz_label_order <-rev(get_taxa_name(updated_tree_viz))

# order the samples on the chart as they appear on the tree by factoring
main_df <-all_labels %>% mutate(person_id=factor(person_id,levels=tree_viz_label_order))
df <-main_df
# ############################### COLOR PALETTE ##################################
##set color palette
# gather facilities
set_color_palette <-function(primary_df){
  pollen <-primary_df
  reduced_facilities_list <-unique(pollen$person_facility_name)
  
  # gather colors and exclude gray/grey
  reduce_color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  # randomly select color for facilities (-1 b/c DOD is a "facility" but will be assigned a clear color "#fffff00")
  reduced_updated_colors <-sample(reduce_color,NROW(reduced_facilities_list)-1)
  
  # set colors
  reduced_facility_colors <-setNames(c("#ffffff00",reduced_updated_colors),
                                     c("DOD", setdiff(reduced_facilities_list, "DOD"))
  )
  return(reduced_facility_colors)
}# end of color palette fxn

#call set_color_palettefxn
location_colors <-set_color_palette(main_df)
# ##############################geom point data####################################
# gather all DODs in a df

# prep DODs_reduced
prep_DODs<-function(DOD_df){
  stem <-DOD_df
  sepal <-stem %>%
    mutate(person_id=paste0("Patient ",person_id)) %>%
    rename(c("date"="patient_date_diagnosed")) %>%
    mutate(date=as.Date(date,format="%Y-%m-%d")) %>%
    mutate(person_facility_name="DOD",
           date_type="DOD") %>%
    select(person_id,person_facility_name,date_type,date)
  return(sepal)
}

# call prep fxn
DODs_reduced_updated<-prep_DODs(DODs_reduced)

#prep tree samples
tree_samples_DODs <-labels_not_inchart_df %>%
  left_join(subset(sample_data,select=c(`#BioSample`,`Collection date`)), by=c("person_id"="#BioSample")) %>%
  select(-date) %>%
  rename(c("date"="Collection date"))


# combine all DOD dfs
all_DODs <-rbind(DODs_reduced_updated,tree_samples_DODs)


############################## NETWORK DATA ###################################

#classify links
determine_link_type <-function(links_df,primary_DOD_df){
  ocean <-links_df
  pond <-primary_DOD_df
  link_type <-ocean %>%
    mutate(person_id=paste0("Patient ",person_id)) %>%
    left_join(subset(pond,select=c(person_id,date))) %>%
    mutate(link_type=case_when(date >=start_date & date <=end_date~"During",
                               date >=start_date~"Before",
                               date <=end_date~"After",
                               .default="Unknown"))
}
# call link_type function
link_status <-determine_link_type(links_reduced,all_DODs)

## create nodelist
# retrieve case info
create_nodelist <-function(DOD_df,amr_df,links_df){
  # retrieve case info
  sample_nodes <-left_join(DOD_df,amr_df) %>% 
    mutate(node_type="case",
           person_id=paste0("Patient ",person_id))%>% 
    rename(c("DOD"="patient_date_diagnosed"))%>% 
    select(person_id,node_type,type,DOD)
  
  # retrieve facility
  facility_nodes<-links_df %>% 
    select(person_facility_name) %>% 
    distinct() %>% 
    mutate(node_type="facility",
           DOD=NA_real_,
           type=case_when(str_detect(person_facility_name,"Ctr")~"Hospital",
                          str_detect(person_facility_name,"Resort")~"LTCF",
                          str_detect(person_facility_name,"Court")~"LTCF",
                          str_detect(person_facility_name,"Surgery")~"Clinic",
                          str_detect(person_facility_name,"Hospital")~"Hospital",
                          str_detect(person_facility_name,"Wound")~"Clinic",
                          str_detect(person_facility_name,"Rehab")~"LTCF",
                          str_detect(person_facility_name,"Hos")~"Hospital",
                          str_detect(person_facility_name,"Med")~"Hospital",
                          str_detect(person_facility_name,"Uni")~"Hospital",
                          str_detect(person_facility_name,"Patho")~"Clinic",
                          str_detect(person_facility_name,"Consult")~"Clinic",
                          str_detect(person_facility_name,"Home")~"LTCF",
                          str_detect(person_facility_name,"Health")~"Hospital",
                          str_detect(person_facility_name,"Corp")~"Clinic",
                          str_detect(person_facility_name,"CORP")~"Clinic",
                          str_detect(person_facility_name,"GROUP")~"Clinic",
                          .default="LTCF")) %>% 
    rename(c("person_id"="person_facility_name")) %>% 
    select(person_id,node_type,type,DOD)
  
  line_list <-rbind(sample_nodes,facility_nodes)
  
  return(line_list)
}
#call "create_nodelist" function
node_list <-create_nodelist(DODs_reduced,reduced_pcr_combined,links_reduced)

# create epicontacts object for viz network
generate_epicontacts_object<-function(links,list){
  contacts_network <-links %>% 
    transmute(
      infector = person_id,
      case_id = person_facility_name,
      location = link_type
    )  %>%
    drop_na(infector) 
  
  
  
  ## generate epicontacts object
  
  EC_object <-make_epicontacts(
    linelist = list,
    contacts = contacts_network,
    id = "person_id",
    from = "infector",
    to = "case_id",
  )
  
  return(EC_object)
}

#call "generate_epicontacts_object" fxn
epic_object <-generate_epicontacts_object(link_status,node_list)

#updated fxn
vis_epicontactsMT <-function (x, thin = TRUE, node_color = "id", label = "id", annot = TRUE, 
                              node_shape = NULL, shapes = NULL, edge_label = NULL, edge_color = NULL, 
                              legend = TRUE, legend_max = 10, x_axis = NULL, col_pal = cases_pal, 
                              NA_col = "lightgrey", edge_col_pal = edges_pal, width = "90%", 
                              height = "700px", selector = TRUE, editor = FALSE, edge_width = 3, 
                              ...) 
{
  if (thin) {
    x <- thin(x)
  }
  if (!is.null(x_axis)) 
    stop(paste("x_axis feature only available in development 'timeline' branch, which can be installed", 
               "via remotes::install_github('reconhub/epicontacts@timeline')."))
  node_color <- epicontacts:::assert_node_color(x, node_color)
  node_shape <- epicontacts:::assert_node_shape(x, node_shape)
  annot <- epicontacts:::assert_annot(x, annot)
  edge_label <- epicontacts:::assert_edge_label(x, edge_label)
  edge_color <- epicontacts:::assert_edge_color(x, edge_color)
  all_nodes <- get_id(x, which = "all")
  nodes <- data.frame(id = all_nodes, stringsAsFactors = FALSE)
  nodes <- merge(nodes, x$linelist, by = "id", all = TRUE)
  if (!is.null(label)) {
    labels <- apply(nodes[, label, drop = FALSE], 1, paste, 
                    collapse = "\n")
    nodes$label <- labels
  }
  if (!is.null(annot)) {
    temp <- nodes[, annot, drop = FALSE]
    temp <- vapply(names(temp), function(e) paste(e, temp[, 
                                                          e], sep = ": "), character(nrow(nodes)))
    nodes$title <- paste("<p>", apply(temp, 1, paste0, collapse = "<br>"), 
                         "</p>")
  }
  if (!is.null(node_color)) {
    node_col_info <- fac2col(factor(nodes[, node_color]), 
                             col_pal, NA_col, legend = TRUE)
    K <- length(node_col_info$leg_lab)
    nodes$color <- node_col_info$color
  }
  if (!is.null(node_shape)) {
    if (is.null(shapes)) {
      msg <- paste("'shapes' needed if 'node_shape' provided;", 
                   "to see codes, node_shape: codeawesome")
      stop(msg)
    }
    vec_node_shapes <- as.character(unlist(nodes[node_shape]))
    shapes["NA"] <- "question-circle"
    unknown_codes <- !shapes %in% names(codeawesome)
    if (any(unknown_codes)) {
      culprits <- paste(shapes[unknown_codes], collapse = ", ")
      msg <- sprintf("unknown icon codes: %s \nto see 'codeawesome'", 
                     culprits)
      stop(msg)
    }
    vec_node_shapes <- paste(vec_node_shapes)
    node_code <- codeawesome[shapes[vec_node_shapes]]
    nodes$shape <- "icon"
    nodes$icon.code <- node_code
    nodes$icon.color <- nodes$color
  }
  else {
    nodes$borderWidth <- 2
  }
  edges <- x$contacts
  edges$width <- edge_width
  if (x$directed) {
    edges$arrows <- "to"
  }
  if (!is.null(edge_label)) {
    edges$label <- edges[, edge_label]
  }
  if (!is.null(edge_color)) {
    edge_col_info <- fac2col(factor(edges[, edge_color]), 
                             edge_col_pal, NA_col, legend = TRUE)
    L <- length(edge_col_info$leg_lab)
    edges$color <- edge_col_info$color
  }
  out <- visNetwork::visNetwork(nodes, edges, width = width, 
                                height = height, ...)
  if (legend) {
    if (!is.null(node_color) && (K < legend_max)) {
      leg_nodes <- data.frame(label = node_col_info$leg_lab, 
                              color = node_col_info$leg_col, shape = "box", 
                              shadow = TRUE, font.size = 20)
    }
    else {
      leg_nodes <- NULL
    }
    if (!is.null(edge_color) && (L < legend_max)) {
      leg_edges <- data.frame(label = edge_col_info$leg_lab, 
                              color = edge_col_info$leg_col, font.size = 15)
    }
    else {
      leg_edges <- NULL
    }
    out <- visNetwork::visLegend(out, addNodes = leg_nodes, 
                                 addEdges = leg_edges, useGroups = FALSE)
  }
  enabled <- list(enabled = TRUE)
  arg_selec <- if (selector) 
    node_color
  else NULL
  out <- visNetwork::visOptions(out, highlightNearest = TRUE)
  out <- visNetwork::visOptions(out, selectedBy = arg_selec, 
                                manipulation = editor, highlightNearest = enabled)
  out <- visNetwork::visPhysics(out, stabilization = FALSE,solver = "barnesHut", barnesHut = list(springConstant = 0,damping=1,springLength=15, gravitationalConstant=-10000))
  out <- visNetwork::addFontAwesome(out)
  return(out)
}


