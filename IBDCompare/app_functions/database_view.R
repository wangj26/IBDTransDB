#==============================================================================================================
# Function used to create a table of database links for each gene
# Returns one dataframe:
# a kable table containing database links per gene
#==============================================================================================================
generate_DB_links <- function(genes, db) {
  
  # Pull the database links from the DB
  db_links <- get_db_info(db)
  
  # Format the database name column to contain the full name and short name
  db_links$database_name <- mapply(function(x,y) if (x != y) paste0(x, " (", y, ")") else x, 
                                   db_links$database_name, 
                                   db_links$database_full_name)
  
  # Popover for the database names that displays the associated summary
  database_popover <- lapply(db_links$database_name, function(x) {
    title <- db_links[db_links$database_name == x, "database_full_name"]
    summary <- db_links[db_links$database_name == x, "summary"]
    popover_summary <- spec_popover(title = title, 
                                    content = summary,
                                    trigger = "hover",
                                    position = "right")
    return(cell_spec(x, popover = popover_summary))
  })
  
  # Create a dataframe used to store DB links for each gene
  db_df <- db_links %>%
    select(c(database_name, type)) %>%
    mutate(database_name = database_popover)
  for (gene in genes) {
    links <- paste0(db_links$link, gene)
    db_df[gene] <- sapply(links, function(x) {
      cell_spec("link", link = x, new_tab = TRUE)
    })
  }
  
  # Create a type factor used to group the table and remove the type column
  type_info <- factor(db_df$type, unique(db_df$type))
  db_df <- select(db_df, -c("type"))
  colnames(db_df)[1] <- ""
  
  ### Create a pretty table that will be used as the output in the Database Links tab ###
  db_df <- db_df %>%
    kable('html', escape = FALSE, row.names = FALSE) %>%
    kable_styling(font_size = 12) %>%
    row_spec(0, font_size = 14, color = "white", background = "#26de5a") %>%
    row_spec(1:nrow(db_df), background = 'white') %>%
    pack_rows(index = table(type_info), background = '#e1e2e3', 
              label_row_css = "border-top: 2px solid black;")
  ###
  
  return(db_df)
}
#==============================================================================================================


#==============================================================================================================
# Function used to create a table of database links for each gene
# Returns one dataframe:
# a kable table containing database links per gene
#==============================================================================================================
generate_DB_summary <- function(db) {
  
  # Pull the database links from the DB
  db_links <- get_db_info(db)
  
  # Create a dataframe used to store DB links for each gene
  db_df <- db_links %>%
    select(c(database_name, summary))
  
  return(db_df)
}
#==============================================================================================================
