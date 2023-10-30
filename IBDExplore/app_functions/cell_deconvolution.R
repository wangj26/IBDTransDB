#==============================================================================================================
# Function used to perform cell deconvolution (using Scaden)
#==============================================================================================================
get_cell_fractions <- function(dataset_acc, sample_group_map, db) {
  
  # Format the dataset ID so that it can be used for the query
  dataset_acc <- paste0("'", dataset_acc, "'")
  
  # Pull the comparisons
  cell_fractions <- dbGetQuery(db,
                               stringr::str_interp(
                                 paste("SELECT sample_id, cell_type, fraction",
                                       "FROM cell_deconvolution",
                                       "WHERE dataset_acc = ${dataset_acc}")))
  
  # Convert the cell fractions data to wide format
  cell_fractions <- cell_fractions %>%
    pivot_wider(names_from = cell_type, values_from = fraction) %>% 
    as.data.frame()
  rownames(cell_fractions) <- cell_fractions$sample_id

  # Add a column specifying the group to the predictions
  sample_group_map$group <- gsub(";$", "", sample_group_map$group)  # MAKE SURE FORMATTING IS FIXED IN DB
  sample_group_map <- unique(sample_group_map)
  sample_group_map <- unique(separate_rows(sample_group_map, samples, sep = ";"))
  cell_fractions <- cell_fractions %>%
    left_join(sample_group_map, by = c("sample_id" = "samples")) %>%
    relocate(group)
  rownames(cell_fractions) <- cell_fractions$sample_id
  cell_fractions <- select(cell_fractions, -sample_id)
  
  return(cell_fractions)
}
#==============================================================================================================


#==============================================================================================================
# Function used to create a DT table from a dataframe of cell fractions
# Returns a DT table
#==============================================================================================================
generate_cell_fractions_table <- function(cell_fractions) {
  
  # Convert the cell fractions to scientific notation
  group_labels <- cell_fractions$group
  cell_fractions <- cell_fractions %>%
    select(-group) %>%
    format(signif(4), scientific = TRUE) %>%
    mutate(group = group_labels) %>%
    relocate(group)
  
  # Generate the DT table
  cell_fractions <- cell_fractions %>%
    arrange(group) %>%
    datatable(extensions = c("Scroller", "FixedColumns", "RowGroup"),
              selection = list(mode = "single", target = "column"),
              options = list(dom = "ti",
                             fixedColumns = list(leftColumns = 2),
                             rowGroup = list(dataSrc = 1),
                             autoWidth = TRUE,
                             scrollX = TRUE,
                             scrollY = 600,
                             scroller = TRUE)) %>%
    formatStyle(columns = 1:ncol(cell_fractions), fontSize = "90%")
  
  return(cell_fractions)
}
#==============================================================================================================


#==============================================================================================================
# Function used to run a test between cell fractions for each cell type given a selected comparison
# Returns a list of p-values
#==============================================================================================================
run_cell_fraction_tests <- function(cell_fractions) {
  
  # Perform a Wilcoxon rank sum test between the two groups for each cell type
  p_values <- sapply(select(cell_fractions, -group), function(x) wilcox.test(x ~ cell_fractions$group)$p.value)

  return(p_values)
}
#==============================================================================================================


#==============================================================================================================
# Function used to compute the mean cell fraction for each cell type
# Returns a list of mean cell fractions
#==============================================================================================================
compute_mean_cell_fraction <- function(cell_fractions) {
  
  # Compute the mean cell fraction for each cell type
  mean_fractions <- sapply(select(cell_fractions, -group), function(x) mean(x))
  
  return(mean_fractions)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a cell fraction box plot for each cell type that passes the specified criteria (p-value and min fraction)
# Returns a list of ggplot objects
#==============================================================================================================
generate_cell_fraction_plots <- function(selected_comparison, cell_fractions, max_p_val, min_fraction) {
  
  # Use the selected comparison to filter the cell fractions table and group the table
  selected_groups <- unlist(strsplit(selected_comparison, split = " vs "))
  cell_fractions <- cell_fractions %>%
    filter(group %in% selected_groups)
  
  # Compute p-values and mean cell fractions
  p_values <- run_cell_fraction_tests(cell_fractions)
  mean_fractions <- compute_mean_cell_fraction(cell_fractions)
  
  # Filter the cell fractions table based on the thresholds for the p-value and cell fraction
  p_values_to_keep <- p_values <= max_p_val
  columns_to_keep <- ifelse(is.na(p_values), FALSE, p_values_to_keep) & (mean_fractions >= min_fraction)
  group_labels <- pull(cell_fractions, group)
  cell_fractions <- cell_fractions %>%
    select(-group) %>%
    select_if(columns_to_keep) %>%
    signif(4)
  # Filter the p-values
  p_values <- p_values[columns_to_keep]
  
  if (sum(columns_to_keep) == 0) {
    return(NULL)  
  }
  
  # Iterate over each cell type and produce an expression plot
  plot_list <- lapply(1:ncol(cell_fractions), function(col_index) {
    
    # Pull the current cell type and associated p-value
    cell_type <- colnames(cell_fractions)[col_index]
    p_val <- p_values[col_index]
    p_val = formatC(p_val, format = "e", digits = 2)
    
    # Create a dataframe used to store the group labels and cell fractions
    data <- data.frame("group" = group_labels,
                       "proportion" = unname(cell_fractions[cell_type]))
    # Format the group labels for plotting
    data <- data %>%
      mutate(group = factor(group, levels = selected_groups),
             group = gsub(";", " ", group),
             group = stringr::str_wrap(group, 10))
    
    # # Generate the y coordinate for the p-value, increase the y position of the p-value, add the p-value, and format group columns
    # p_val_position <- get_y_position(data,
    #                                  proportion ~ group,
    #                                  comparisons = list(unique(group_labels)),
    #                                  step.increase = 0.1) %>%
    #   select(-c("groups")) %>%
    #   mutate(y.position = y.position + (max(data$proportion) - min(data$proportion))/20,
    #          p_value = formatC(p_val, format = "e", digits = 2)) %>%
    #   mutate(group1 = gsub(";", " ", group1),
    #          group2 = gsub(";", " ", group2)) %>%
    #   mutate(group1 = stringr::str_wrap(group1, 10), 
    #          group2 = stringr::str_wrap(group2, 10))
    
    # Generate the y coordinates for the p-value and p-value line
    p_val_line_position <- max(data$proportion) + (max(data$proportion) - min(data$proportion))*0.1
    p_val_position <- max(data$proportion) + (max(data$proportion) - min(data$proportion))*0.15
    
    # Define the color palette
    plot_colors <- brewer.pal(8, "Dark2")[1:2]
    
    # Boxplot of expression for the given comparison and gene
    cell_fraction_plot <- ggplot(data, aes(x = group, y = proportion, fill = group)) + 
      geom_boxplot() +
      annotate("text", x = 1.5, y = p_val_position, label = p_val, size = 4) +
      geom_segment(aes(x = 1, y = p_val_line_position, xend = 2, yend = p_val_line_position)) +
      xlab(NULL) +
      ylab("Proportion") +
      ggtitle(str_trunc(cell_type, width = 30, ellipsis = "...")) +
      scale_fill_manual(values = plot_colors) +
      theme(panel.background = element_rect(fill = "white"), 
            plot.title = element_text(size = 14),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
            axis.text.y = element_text(size = 12),
            legend.position = "None",
            panel.border = element_rect(colour = "black", fill = NA))
    
    return(cell_fraction_plot)
  })
  
  # Define the names for the plot list (necessary to create plotOutput objects)
  names(plot_list) <- paste0("cell_fraction_", colnames(cell_fractions), "_plot")
  
  return(plot_list)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate a cell fraction boxplot for each selected comparison (p-value is displayed)
# Returns a list of ggplot objects
#==============================================================================================================
# generate_cell_fraction_plots <- function(selected_comparisons, selected_column, cell_fractions) {
#   
#   # Pull the cell type associated with the selected column
#   selected_cell_type <- colnames(cell_fractions)[selected_column]
#   
#   # Keep the group identifier and the cell fractions corresponding to the selected cell type
#   cell_type_col <- which(colnames(cell_fractions) == selected_cell_type)
#   cell_fractions <- cell_fractions[,c(1, cell_type_col)]
#   colnames(cell_fractions)[2] <- "proportion"
#   cell_fractions$proportion <- signif(cell_fractions$proportion, 4)
#     
#   # Iterate over each selected comparison and produce an expression plot
#   plot_list <- lapply(selected_comparisons, function(x) {
#     
#     # Use the selected comparison to filter the cell fractions table
#     selected_groups <- unlist(strsplit(x, split = " vs "))
#     plot_data <- cell_fractions %>%
#       filter(group %in% selected_groups) %>%
#       mutate(group = gsub(";", "; ", group)) %>%
#       mutate(group = stringr::str_wrap(group, 10))
#     
#     # Define the color palette
#     plot_colors <- brewer.pal(8, "Dark2")[1:2]
#     
#     # Boxplot of expression for the given comparison and gene
#     cell_fraction_plot <- ggplot(plot_data, aes(x = group, y = proportion, fill = group)) + 
#       geom_boxplot() +
#       xlab(NULL) +
#       ylab("Proportion") +
#       ggtitle(str_trunc(selected_cell_type, width = 30, ellipsis = "...")) +
#       scale_fill_manual(values = plot_colors) +
#       theme(panel.background = element_rect(fill = "white"), 
#             plot.title = element_text(size = 14),
#             axis.title.y = element_text(size = 10),
#             axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
#             axis.text.y = element_text(size = 12),
#             legend.position = "None",
#             panel.border = element_rect(colour = "black", fill = NA))
#     
#     return(cell_fraction_plot)
#   })
#   
#   # Define the names for the plot list (necessary to create plotOutput objects)
#   names(plot_list) <- paste0("cell_fraction_", gsub(" ", "_", selected_comparisons), "_plot")
#   
#   return(plot_list)
# }
#==============================================================================================================


#==============================================================================================================
# Function used to determine the plot layout for the cell fraction box plots
#==============================================================================================================
determine_cell_fraction_layout <- function(num_plots) {
  
  # Maximum number of plots to display is 9
  # if (num_plots > 12) {
  #   num_plots <- 12
  # }
  
  # Determine the plot layout based on the number of selected rows in the comparison metrics table
  num_row <- num_plots%/%4
  num_remaining <- num_plots%%4
  if (num_remaining == 0) {
    col_1_indices <- 1:num_row
    col_2_indices <- (num_row + 1):(2*num_row)
    col_3_indices <- (2*num_row + 1):(3*num_row)
    col_4_indices <- (3*num_row + 1):num_plots
  } else if (num_remaining == 1) {
    col_1_indices <- 1:(num_row + 1)
    col_2_indices <- (num_row + 2):(2*num_row + 1)
    col_3_indices <- (2*num_row + 2):(3*num_row + 1)
    col_4_indices <- (3*num_row + 2):num_plots
  } else if (num_remaining == 2) {
    col_1_indices <- 1:(num_row + 1)
    col_2_indices <- (num_row + 2):(2*num_row + 2)
    col_3_indices <- (2*num_row + 3):(3*num_row + 2)
    col_4_indices <- (3*num_row + 3):(4*num_row + 2)
  } else {
    col_1_indices <- 1:(num_row + 1)
    col_2_indices <- (num_row + 2):(2*num_row + 2)
    col_3_indices <- (2*num_row + 3):(3*num_row + 3)
    col_4_indices <- (3*num_row + 4):(4*num_row + 3)
  }
  
  return(list(col_1_indices, col_2_indices, col_3_indices, col_4_indices))
}
#==============================================================================================================
