#'Check Database structure
#'
#'This function will examine a SQLite3 database to make sure it's consistent
#'with the structure required for my Shiny applications.
#'
#'@param db.file a character object that gives the file location of a SQLite3 database.
#'@param verbose a logical value indicating if errors should be displayed. If set to FALSE, all messages will be suppressed.
#'
#'@return A logical value indicating if the database passes all required checks.
#'
#'@examples
#'database_file <- "//AIBSData2/mct-t200/ShinyApps/linnarsson/linnarsson2.db"
#'
#'check_db_structure(database_file)
#'
#'check_db_structure(database_file,verbose=F)
check_db_structure <- function(db.file,verbose=T) {
  library(DBI)
  library(RSQLite)

  # 5 checks need to pass for the database to work with shiny
  passed_checks <- 0
  
  con <- dbConnect(RSQLite::SQLite(),db.file)
  
  # Check to see if the 3 tables are available.
  get <- "SELECT name FROM sqlite_master WHERE type='table';"
  res <- dbSendQuery(con,get)
  db_tables <- dbFetch(res,n=-1)
  dbClearResult(res)
  
  expected_tables <- c("desc","anno","data")
  
  for(table_name in expected_tables) {
    if(table_name %in% db_tables$name) {
      if(verbose) { cat(paste0(table_name," table exists :)\n")) }
      passed_checks <- passed_checks + 1
    } else { 
      if(verbose) { cat(paste0("no ",table_name," table! :'(\n!!This table is required for the shiny app!!\n\n")) }
    } 
  }
  
  # If desc exists, read it to get the annotations we expect in anno.
  if(sum(c("desc","anno") %in% db_tables$name) == 2) {
    get <- paste("SELECT * FROM desc",sep="")
    res <- dbSendQuery(con,get)
    desc <- dbFetch(res,n=-1)
    dbClearResult(res)
    
    if(verbose) {
      cat("Description table:\n")
      print(desc)
      cat("\n")
    }
    
    expected_anno_columns <- c("sample_id",
                               paste0(desc$base,"_id"),
                               paste0(desc$base,"_label"),
                               paste0(desc$base,"_color"))
    
    get <- "PRAGMA table_info('anno')"
    res <- dbSendQuery(con,get)
    anno_columns <- dbFetch(res,n=-1)
    
    if(sum(expected_anno_columns %in% anno_columns$name) == length(expected_anno_columns)) {
      if(verbose) { cat("all described columns were found! :)\n") }
      passed_checks <- passed_checks + 1
    } else if (sum(expected_anno_columns %in% anno_columns$name) < length(expected_anno_columns)) {
      if(verbose) { cat("some described columns were not found! :'(\n")
      for(col_name in expected_anno_columns[!expected_anno_columns %in% anno_columns$name]) {
        cat(paste0(col_name,"\n\n"))
      }
      cat("!!this will cause errors in the shiny app!!\n\n")
      }
    }
    
    if(nrow(anno_columns) > length(expected_anno_columns)) {
      if(verbose) {
        cat("some annotations were not found in the desc table:\n")
        for(col_name in anno_columns$name[!anno_columns$name %in% expected_anno_columns]) {
          cat(paste0(col_name,"\n"))
        }
        cat("this will work with the shiny app, but these columns will be hidden.\n\n")
      }
    }
    
  }
  
  # Check to see if all filenames in anno are found in rpkm (and vice-versa)
  if(sum(c("anno","data") %in% db_tables$name) == 2) {
    
    get <- paste("SELECT sample_id FROM anno",sep="")
    res <- dbSendQuery(con,get)
    anno_filenames <- dbFetch(res,n=-1)
    anno_filenames <- anno_filenames[,1]
    dbClearResult(res)
    
    if(verbose) { cat(paste0("found annotations for ",length(anno_filenames)," cells.\n")) }
    
    expected_rpkm_columns <- c("gene",anno_filenames)
    
    get <- "PRAGMA table_info('data')"
    res <- dbSendQuery(con,get)
    rpkm_columns <- dbFetch(res,n=-1)
    dbClearResult(res)
    
    if(verbose) { cat(paste0("found data for ",nrow(rpkm_columns) - 1," cells.\n")) }
    
    if(sum(expected_rpkm_columns %in% rpkm_columns$name) == length(expected_rpkm_columns)) {
      if(verbose) { cat("all annotations have corresponding data! :)\n") }
      passed_checks <- passed_checks + 1
    } else if (sum(expected_rpkm_columns %in% rpkm_columns$name) < length(expected_rpkm_columns)) {
      if(verbose) {
        cat("some annotated cells don't have data! :'(\n")
        for(col_name in expected_rpkm_columns[!expected_rpkm_columns %in% rpkm_columns$name]) {
          cat(paste0(col_name,"\n"))
        }
        cat("!!this will cause errors in the shiny app!!\n\n")
      }
    }
    
    if(nrow(rpkm_columns) > length(expected_rpkm_columns)) {
      if(verbose) {
        cat("some cells have data, but are not annotated:\n")
        for(col_name in rpkm_columns$name[!rpkm_columns$name %in% expected_rpkm_columns]) {
          cat(paste0(col_name,"\n"))
        }
        cat("this will work with the shiny app, but the data for these cells will be hidden.\n\n")
        
      }
    }
    
  }
  
  dbDisconnect(con)

  if(passed_checks == 5) {
    if(verbose) { cat("The database has passed all essential checks! HOORAY!!\n") 
    } 
    return(TRUE)
  } else {
    if(verbose) { cat("The database has some problems that will prevent it from working with shiny.\n") 
    } 
    return(FALSE)
    
  }
  
}

write_database <- function(desc=NULL,anno=NULL,data=NULL,file = stop("'file' must be specified"),overwrite=TRUE) {
  library(DBI)
  library(RSQLite)
  
  # open connection to the database
  con <- dbConnect(RSQLite::SQLite(),file)
  
  if(!is.null(desc)) {
    dbWriteTable(con,"desc",desc,row.names=F,overwrite=overwrite)
  }
  
  if(!is.null(anno)) {
    dbWriteTable(con,"anno",anno,row.names=F,overwrite=overwrite)
  }
  
  if(!is.null(data)) {
    dbWriteTable(con,"data",data,row.names=F,overwrite=overwrite)
    idx <- "CREATE INDEX idx ON data (gene)"
    res <- dbSendQuery(con,idx)
    dbClearResult(res)
  }
  
  dbDisconnect(con)
  
}

auto_sample_names <- function(n) {
  
  sample_names <- character()
  
  for(i in 1:n) {
    sample_name <- "S"
    
    if(nchar(i) < nchar(n)) {
      sample_name <- paste0(sample_name,paste0(rep("0",nchar(n)-nchar(i)),collapse=""))
    }
    
    sample_name <- paste0(sample_name,i)
    
    sample_names <- c(sample_names,sample_name)
  }
  
  return(sample_names)
  
}

auto_format <- function(df,sample_id=NULL,anno_cols=NULL,data_cols=NULL) {
  library(dplyr)
  
  # Build the desc table based on the annotation columns
  if(!is.null(anno_cols)) {
    if(is.numeric(anno_cols)) {
      anno_names <- names(df)[anno_cols]
    } else if (is.character(anno_cols)) {
      anno_names <- anno_cols
    }
    desc <- data.frame(base=anno_names,name=anno_names)
  } else {
    desc <- data.frame(base=character(),name=character())
  }
  
  # Get the sample_id data set.
  if(!is.null(sample_id)) {
    sample_ids <- sample_id
  } else {
    sample_ids <- auto_sample_names(nrow(df))
  }

  anno <- data.frame(sample_id=sample_ids)
  
  if(!is.null(anno_cols)) {
    label_names <- paste0(anno_cols,"_label")
    id_names <- paste0(anno_cols,"_id")
    color_names <- paste0(anno_cols,"_color")
    
    anno <- cbind(anno,df[,anno_cols])
    names(anno)[2:length(anno)] <- label_names
    
    for(i in 1:length(anno_cols)) {
      anno_set <- anno %>% select_(label_names[i]) %>%
        arrange_(label_names[i]) %>%
        unique() %>%
        mutate(id=1:n(),
               color=rainbow(n()))
      names(anno_set)[2:3] <- c(id_names[i],color_names[i])
      
      anno <- anno %>% left_join(anno_set,by=label_names[i])
      df <- df %>% select_(paste0("-",anno_cols[i]))
    }
  }
  
  if(!is.null(data_cols)) {
  } else {
    data <- df[,sapply(df,is.numeric)]
    data <- cbind(data,df[,sapply(df,is.integer)])
    data <- data %>% t() %>% as.data.frame()
    names(data) <- sample_ids
    gene <- rownames(data)
    data <- cbind(gene,data)
  }
  
  return(list(desc=desc,anno=anno,data=data))
  
}