#==============================================================================================================
# Function used to perform PCA on the expression data
#==============================================================================================================
get_pca <- function(dataset_acc, sample_ann, db) {
  
  # Format the dataset ID so that it can be used for the query
  dataset_acc <- paste0("'", dataset_acc, "'")
  
  # Pull the PCA coordinates
  pca_df <- dbGetQuery(db,
                       stringr::str_interp(
                         paste("SELECT sample_id, pc1, pc2, pc3, pc4, pc5",
                               "FROM pca",
                               "WHERE dataset_acc = ${dataset_acc}")))
  colnames(pca_df)[2:6] <- toupper(colnames(pca_df)[2:6])
  
  ### Add sample level variables to the PCA dataframe ###
  sample_ann <- sample_ann %>%
    mutate(sample_ann_type = str_split(sample_ann_type, ";"),
           sample_ann_value = str_split(sample_ann_value, ";"))
  annotation_types <- unique(unlist(sample_ann$sample_ann_type))
  # Add a column to the PCA dataframe for the given annotation type if it applies to all samples
  for (type in annotation_types) {
    
    # If the given annotation type is not found for every sample, continue to the next iteration
    type_freq <- sum(unlist(sample_ann$sample_ann_type) == type)
    if (type_freq != nrow(sample_ann)) {
      next
    }
    
    type_indices <- sapply(sample_ann$sample_ann_type, function(x) {
      type_index <- which(x == type)
      return(type_index)
    })
    
    # Add a column to the PCA dataframe for the given annotation type (e.g. Disease)
    pca_df[type] <- mapply(function(x,y) y[x], 
                           type_indices, 
                           sample_ann$sample_ann_value)
  }
  ##########
  
  # Pull the PC variance
  pc_variance <- dbGetQuery(db,
                            stringr::str_interp(
                              paste("SELECT pc, variance",
                                    "FROM pca_variance",
                                    "WHERE dataset_acc = ${dataset_acc}")))
  colnames(pc_variance) <- c("PC", "percent_variance")
  
  return(list("pca_df" = pca_df, "pc_variance" = pc_variance))
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate scree plots 
#==============================================================================================================
generate_scree_plot <- function(pc_variance) {
  
  # Create scree plot
  scree_plot <- ggplot(data = pc_variance, aes(x = PC, y = percent_variance)) +
    geom_bar(stat = "identity") +
    geom_point(color = "red") + geom_line(aes(group = 1), color = "red") +
    ggtitle("Scree Plot") +
    xlab("PC") + ylab("Percent Variance Explained") +
    theme(legend.position = "none", panel.background = element_rect(fill = "white"), 
          plot.title = element_text(size = 18), axis.title = element_text(size = 12), 
          axis.text = element_text(size = 14), panel.border = element_blank(),
          axis.line.x.bottom = element_line(color = 'black'),
          axis.line.y.left   = element_line(color = 'black'))

  return(scree_plot)
}
#==============================================================================================================


#==============================================================================================================
# Function used to generate PCA plots (colored by the user specified variable)
#==============================================================================================================
generate_pca_plot <- function(pca_df, pc_variance, selected_var, selected_filter_feature, selected_filters, pc_x_axis, pc_y_axis, plot_ellipse, hide_legend) {
  
  if (is.null(selected_var)) {
    selected_var <- "---"
    selected_filter_feature <- "---"
    pc_x_axis <- 1
    pc_y_axis <- 2
    plot_ellipse <- FALSE
    hide_legend <- FALSE
  }
  
  pc_x_axis <- as.numeric(pc_x_axis)
  pc_y_axis <- as.numeric(pc_y_axis)
  
  # Define the column names for the chosen PCs
  pc_x_label <- paste0("PC", pc_x_axis)
  pc_y_label <- paste0("PC", pc_y_axis)
  
  # Apply a filter to mask datapoints if specified
  if (selected_filter_feature != "---" & length(selected_filters) > 0) {
    pca_df <- pca_df %>%
      filter(pca_df[[selected_filter_feature]] %in% selected_filters)
  }
  
  if (selected_var == "---") {
    # Create a PCA plot (no group specified)
    pca_plot <- ggplot(data = pca_df, aes_string(x = pc_x_label, y = pc_y_label)) +
      geom_point(size = 3) + 
      ggtitle("PCA") +
      xlab(paste0(pc_x_label, " - ", round(pc_variance$percent_variance[pc_x_axis], 2), "% variance explained")) +
      ylab(paste0(pc_y_label, " - ", round(pc_variance$percent_variance[pc_y_axis], 2), "% variance explained")) +
      theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 18), 
            axis.title = element_text(size = 12), axis.text = element_text(size = 14),
            panel.border = element_rect(color = "black", fill = NA))
  } else {
    # Color palette for plotting
    num_colors <- length(unique(pca_df[[selected_var]]))
    if (num_colors <= 8) {
      plot_colors <- brewer.pal(num_colors, "Dark2")
    } else {
      plot_colors <- c(brewer.pal(8, "Dark2"), brewer.pal(num_colors - 8, "Set1"))
    }
    
    # Create a PCA plot colored by the specified group
    pca_plot <- ggplot(data = pca_df, aes_string(x = pc_x_label, y = pc_y_label, color = selected_var)) +
      geom_point(size = 3) + 
      ggtitle("PCA") +
      xlab(paste0(pc_x_label, " - ", round(pc_variance$percent_variance[pc_x_axis], 2), "% variance explained")) +
      ylab(paste0(pc_y_label, " - ", round(pc_variance$percent_variance[pc_y_axis], 2), "% variance explained")) +
      scale_colour_manual(name = selected_var, values = plot_colors) + 
      theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 18), 
            axis.title = element_text(size = 12), axis.text = element_text(size = 14),
            panel.border = element_rect(color = "black", fill = NA), 
            legend.title = element_text(size = 10), legend.text = element_text(size = 10),
            legend.background = element_rect(size = 0.5, linetype = "solid", color = "black"),
            legend.key.size = unit(0.5, "cm"), legend.key = element_blank())
    
    # Plot ellipse around groups if specified
    if (plot_ellipse) {
      pca_plot <- pca_plot + stat_ellipse()
    } 
    
    # Hide the legend if specified
    if (hide_legend) {
      pca_plot <- pca_plot + theme(legend.position = "none")
    } 
  }
  
  return(pca_plot)
}
#==============================================================================================================


#==============================================================================================================
# Function used to determine the association between the selected PC and the selected feature
# Runs a Wilcoxon rank sum test
#==============================================================================================================
run_pc_test <- function(pca_df, selected_pc, selected_var, selected_filter_feature, selected_filters, group1, group2) {
  
  # Format the selected PC
  selected_pc <- paste0("PC", selected_pc)
  
  # If a filter was specified, filter the PCA data
  if (selected_filter_feature != "---" & length(selected_filters) > 0) {
    pca_df <- pca_df %>%
      filter(pca_df[[selected_filter_feature]] %in% selected_filters)
  }
  
  # Pull the PC coordinates for group 1 given the selected feature and PC
  group1_values <- pca_df %>%
    filter(pca_df[selected_var] == group1) %>%
    pull(selected_pc)
  # Pull the PC coordinates for group 2 given the selected feature and PC
  group2_values <- pca_df %>%
    filter(pca_df[selected_var] == group2) %>%
    pull(selected_pc)
  
  # If there are not sufficient values to run the test, return NULL
  if (length(group1_values) == 0 | length(group2_values) == 0) {
    return(NULL)
  }
  
  # Run a Wilcoxon rank sum test
  test_results <- wilcox.test(group1_values, group2_values, alternative = "two.sided")
  p_val <- signif(test_results$p.value, 4)
  
  return(p_val)
}
#==============================================================================================================


#==============================================================================================================
# Function used to create a box plot of PC coordinates for the selected PC and the selected feature
# Returns a ggplot object
#==============================================================================================================
generate_pc_boxplot <- function(pca_df, selected_pc, selected_var, selected_filter_feature, selected_filters, group1, group2, p_value) {
  
  # Format the selected PC
  selected_pc <- paste0("PC", selected_pc)
  
  # If a filter was specified, filter the PCA data
  if (selected_filter_feature != "---" & length(selected_filters) > 0) {
    pca_df <- pca_df %>%
      filter(pca_df[[selected_filter_feature]] %in% selected_filters)
  }
  
  # Pull the PC coordinates for the selected PC and the two selected groups
  pca_df <- pca_df %>%
    select(c(selected_pc, selected_var)) %>%
    filter(unlist(pca_df[selected_var]) %in% c(group1, group2))
  # Rename the columns
  colnames(pca_df) <- c("pc", "group")
  
  # # Generate the y coordinate for the p-value, increase the y position of the p-value, and add the p-value
  # p_val_position <- get_y_position(pca_df,
  #                                  pc ~ group,
  #                                  comparisons = list(c(group1, group2)),
  #                                  step.increase = 0.1) %>%
  #   select(-groups) %>%
  #   mutate(y.position = y.position + (max(pca_df$pc) - min(pca_df$pc))/20,
  #          p_val = p_value)
  
  # Generate the y coordinates for the p-value and p-value line
  pc_values <- pull(pca_df, pc)
  p_val_line_position <- max(pc_values) + (max(pc_values) - min(pc_values))*0.1
  p_val_position <- max(pc_values) + (max(pc_values) - min(pc_values))*0.18
  
  # Define the color palette
  plot_colors <- brewer.pal(8, "Dark2")[1:2]
  
  # Boxplot of PC coordinates
  pc_boxplot <- ggplot(pca_df, aes(x = group, y = pc, fill = group)) + 
    geom_boxplot() +
    annotate("text", x = 1.5, y = p_val_position, label = p_value, size = 4) +
    geom_segment(aes(x = 1, y = p_val_line_position, xend = 2, yend = p_val_line_position)) +
    xlab(NULL) +
    ylab(selected_pc) +
    ggtitle(NULL) +
    scale_fill_manual(values = plot_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme(panel.background = element_rect(fill = "white"), 
          plot.title = element_text(size = 14),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
          axis.text.y = element_text(size = 12),
          legend.position = "None",
          panel.border = element_rect(colour = "black", fill=NA),
          plot.margin = margin(1, 1, 1, 40),
          aspect.ratio = 0.5)
  
  return(pc_boxplot)
}
#==============================================================================================================
