#==============================================================================================================
# Function used to generate volcano plots
# Returns a ggplot object
#==============================================================================================================
volcano_plot <- EnhancedVolcano(voclano_data,
                                lab = voclano_data$gene,
                                x = "logFC",
                                y = "adj.P.Val",
                                selectLab = c('OSM', 'JAK1', "IL1B", "CCL7", "IFIT3"),
                                xlab = bquote(~Log[2]~ '(fold change)'),
                                pCutoff = 0.05,
                                FCcutoff = 1.5,
                                pointSize = 3,
                                labSize = 5,
                                col = c('black', 'black', 'black', 'red3'),
                                legendPosition = NULL,
                                drawConnectors = TRUE,
                                widthConnectors = 2)
#==============================================================================================================
