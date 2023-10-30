#==============================================================================================================
# Function used to pull genesets based on the user specified queries
# Returns one kable table and one dataframe:
# a kable table containing geneset descriptions
# a dataframe containing associated genes for each geneset
#==============================================================================================================
generate_genesets <- function(genes, diseases, organisms, db) {
  # Pull the list of data ID's
  genesetDesc <- get_genesetDesc(diseases, organisms, db)
  
  # If the disease was not found in the DB, return NULL
  if (nrow(genesetDesc) == 0) {
    return(NULL)
  }
  
  # Pull the geneset ID's from the geneset description
  unique_geneset_id <- genesetDesc$geneset_id
  # Pull the unique organisms
  unique_organisms <- unique(genesetDesc$organism)
  
  # Pull the geneset data
  genesetData <- get_genesetData(unlist(genesetDesc$geneset_id), db)
  
  # Create a column in the geneset description that states the number of user specified genes found in the given geneset
  genes_in_geneset <- sapply(unique_geneset_id, function(x) {
    geneset_genes <- genesetData[genesetData$geneset_id == x, 'gene']
    num_genes <- sum(genes %in% geneset_genes)
    return(paste0(num_genes, "/", length(genes)))
  })
  # genesetDesc['# of user genes in geneset'] <- genes_in_geneset
  
  # Create a column per gene in the geneset description that states whether or not the gene is in the given geneset
  for (gene in genes) {
    genesetDesc[gene] <- lapply(unique_geneset_id, function(x) {
      geneset_genes <- genesetData[genesetData$geneset_id == x, 'gene']
      if (gene %in% geneset_genes) {
        return(as.character(span(style = "color: green;", icon("check", class = "check-icon"))))
      } else {
        return(as.character(span(style = "color: darkred", icon("times", class = "x-icon"))))
      }
    })
  }
  
  # genes_in_geneset <- sapply(unique_geneset_id, function(x) {
  #   geneset_genes <- genesetData[genesetData$geneset_id == x, 'gene']
  #   genes_included <- genes[genes %in% geneset_genes]
  #   genes_included <- paste(genes_included, collapse = ", ")
  #   return(genes_included)
  # })
  # genesetDesc['user genes in geneset'] <- genes_in_geneset
  
  
  # Create a column that contains download buttons for downloading member genes
  # download_genes_buttons <- lapply(genesetDesc$geneset_id, function(x) {
  #   tags$button(id = paste0(x, "_download_genes"),
  #               type = "button",
  #               class = "btn action-button btn-primary",
  #               HTML("download genes"))
  # })
  genesetDesc["download genes"] <- paste0('<input type="button" id="', 
                                          genesetDesc$geneset_id,
                                          '_download_genes" value="download genes">')
  
  
  ### Pretty table that will be used as the output in the geneset tab ###
  
  # A disease factor used to group the table
  row_headers <- factor(genesetDesc$disease, unique(genesetDesc$disease))
  
  # Popover for the geneset names that displays the associated summary
  geneset_id_popover <- lapply(genesetDesc$geneset_id, function(x) {
    title <- genesetDesc[genesetDesc$geneset_id == x, "title"]
    summary <- genesetDesc[genesetDesc$geneset_id == x, "summary"]
    popover_summary <- spec_popover(title = title, 
                                    content = summary,
                                    trigger = "hover",
                                    position = "bottom")
    return(cell_spec(x, popover = popover_summary))
  })
  
  # Remove the title and summary columns from the geneset description and rename some columns
  output_table <- genesetDesc %>%
    arrange(disease) %>%
    mutate(geneset_id = geneset_id_popover) %>%
    select(-c("title", "summary", "disease")) %>%
    rename("geneset ID" = "geneset_id", "# of genes" = "gene_number") %>%
    relocate("download genes", .after = "geneset ID") %>%
    select(-c("# of genes"))
  
  # Create the table
  output_table <- output_table %>% kable('html', escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#26de5a", extra_css = "text-align: left") %>%
    row_spec(1:nrow(output_table), background = 'white') %>%
    column_spec(5, extra_css = "text-align: left") %>%
    pack_rows(index = table(row_headers), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black;")
  ###
  
  
  # Merge the geneset description and geneset data
  genesetDesc <- left_join(genesetDesc, genesetData, by = c('geneset_id' = 'geneset_id'))
  # Select the relevant information to be displayed in the geneset section next to the geneset table
  genesetDesc <- select(genesetDesc, c("geneset_id", "title", "summary", "gene"))
  
  return(list("table_1" = output_table, "table_2" = genesetDesc))
}
#==============================================================================================================
