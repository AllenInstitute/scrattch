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
  db_tables <- dbListTables(con)

  expected_tables <- c("desc","anno","data")
  
  for(table_name in expected_tables) {
    if(table_name %in% db_tables) {
      if(verbose) { cat(paste0(table_name," table exists :)\n")) }
      passed_checks <- passed_checks + 1
    } else { 
      if(verbose) { cat(paste0("no ",table_name," table! :'(\n!!This table is required for the shiny app!!\n\n")) }
    } 
  }
  
  # If desc exists, read it to get the annotations we expect in anno.
  if(sum(c("desc","anno") %in% db_tables) == 2) {
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
  if(sum(c("anno","data") %in% db_tables) == 2) {
    
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

#'Read a SQLite database into a list
#'
#'This function will read a scrattch SQLite database into memory as a list.
#'
#'Each table (desc, anno, and data) will be a data frame within the list.
#'
#'@param db.file - a character object specifying the location of the database
#'@param genes - a character object with a list of genes to retrieve. NULL retrieves all genes. default = NULL
#'@param group_by - a character object with the base of the annotation column to filter on if only specific groups should be selected. NULL retrieves all groups. default = NULL
#'@param groups - a numeric object with the group indexes to be selected. NULL retrieves all groups. default = NULL
#'
#'@return A list object containing data.frames for each table in the database, filtered according to genes, group_by, and groups if provided.
#'
db_to_list <- function(db.file,get_tables=NULL,genes=NULL,group_by=NULL,groups=NULL) {
  library(DBI)
  library(RSQLite)
  library(dplyr)
  
  con <- dbConnect(RSQLite::SQLite(),db.file)
  
  out_list <- list()
  
  if(is.null(get_tables) | "desc" %in% get_tables) {
  
    desc <- dbReadTable(con,"desc")
    
    out_list[["desc"]] <- desc
    
  }
  
  anno <- dbReadTable(con,"anno")
  
  if(is.null(group_by) | is.null(groups)) {
    get_samples <- "*"
  } else {
    group_id <- paste0(group_by,"_id")
    anno <- anno[anno[,group_id] %in% groups,]
    get_samples <- paste(c("gene",anno$sample_id),collapse=",")
    
  }
  
  if(is.null(get_tables) | "anno" %in% get_tables) {
    
    out_list[["anno"]] <- anno 
    
  }
  
  if(is.null(get_tables) | "data" %in% get_tables) {
  
  if(is.null(genes)) {
    get <- paste("SELECT ",get_samples," FROM data",sep="")
    data <- dbGetQuery(con,get)
  } else {
    get_genes <- chr_to_sql(genes)
    get <- paste0("SELECT ",get_samples," FROM data WHERE gene IN ",get_genes)
    data <- dbGetQuery(con,get)
  }
  
  data <- data %>%
    select(one_of("gene",anno$sample_id))
    
    out_list[["data"]] <- data
  }
  
  dbDisconnect(con)
  
  return(out_list)
}



#'Check List structure
#'
#'This function will examine a list of tables to make sure it's consistent
#'with the structure required for scrattch applications.
#'
#'@param list_data - a list object that contains data, anno, and desc tables.
#'@param verbose a logical value indicating if errors should be displayed. If set to FALSE, all messages will be suppressed.
#'
#'@return A logical value indicating if the list data passes all required checks.
#'
#'@examples
#'dataset_list <- list(data = v1_data, anno = v1_anno, desc = v1_desc)
#'
#'check_list_structure(dataset_list)
#'
#'check_db_structure(dataset_list,verbose=F)
check_list_structure <- function(list_data,verbose=T) {
  # 5 checks need to pass for the dataset to work with shiny
  passed_checks <- 0
  
  # Check to see if the 3 tables are available.
  db_tables <- list_data
  
  expected_tables <- c("desc","anno","data")
  
  for(table_name in expected_tables) {
    if(table_name %in% names(db_tables)) {
      if(verbose) { cat(paste0(table_name," table exists :)\n")) }
      passed_checks <- passed_checks + 1
    } else { 
      if(verbose) { cat(paste0("no ",table_name," table! :'(\n!!This table is required for scrattch!!\n\n")) }
    } 
  }
  
  # If desc exists, read it to get the annotations we expect in anno.
  if(sum(c("desc","anno") %in% names(db_tables)) == 2) {
    desc <- db_tables$desc
    
    if(verbose) {
      cat("Description table:\n")
      print(desc)
      cat("\n")
    }
    
    expected_anno_columns <- c("sample_id",
                               paste0(desc$base,"_id"),
                               paste0(desc$base,"_label"),
                               paste0(desc$base,"_color"))
    
    anno_columns <- names(db_tables$anno)
    
    if(sum(expected_anno_columns %in% anno_columns) == length(expected_anno_columns)) {
      if(verbose) { cat("all described columns were found! :)\n") }
      passed_checks <- passed_checks + 1
    } else if (sum(expected_anno_columns %in% anno_columns) < length(expected_anno_columns)) {
      if(verbose) { cat("some described columns were not found! :'(\n")
        for(col_name in expected_anno_columns[!expected_anno_columns %in% anno_columns]) {
          cat(paste0(col_name,"\n\n"))
        }
        cat("!!this will cause errors in scrattch!!\n\n")
      }
    }
    
    if(length(anno_columns) > length(expected_anno_columns)) {
      if(verbose) {
        cat("some annotations were not found in the desc table:\n")
        for(col_name in anno_columns[!anno_columns %in% expected_anno_columns]) {
          cat(paste0(col_name,"\n"))
        }
        cat("this will work with scrattch, but these columns may not be fully functional.\n\n")
      }
    }
    
  }
  
  # Check to see if all filenames in anno are found in rpkm (and vice-versa)
  if(sum(c("anno","data") %in% names(db_tables)) == 2) {
    
    anno_filenames <- db_tables$anno$sample_id

    if(verbose) { cat(paste0("found annotations for ",length(anno_filenames)," cells.\n")) }
    
    expected_rpkm_columns <- c("gene",anno_filenames)
    
    rpkm_columns <- names(db_tables$data)

    if(verbose) { cat(paste0("found data for ",nrow(rpkm_columns) - 1," cells.\n")) }
    
    if(sum(expected_rpkm_columns %in% rpkm_columns) == length(expected_rpkm_columns)) {
      if(verbose) { cat("all annotations have corresponding data! :)\n") }
      passed_checks <- passed_checks + 1
    } else if (sum(expected_rpkm_columns %in% rpkm_columns) < length(expected_rpkm_columns)) {
      if(verbose) {
        cat("some annotated cells don't have data! :'(\n")
        for(col_name in expected_rpkm_columns[!expected_rpkm_columns %in% rpkm_columns]) {
          cat(paste0(col_name,"\n"))
        }
        cat("!!this will cause errors in scrattch!!\n\n")
      }
    }
    
    if(length(rpkm_columns) > length(expected_rpkm_columns)) {
      if(verbose) {
        cat("some cells have data, but are not annotated:\n")
        for(col_name in rpkm_columns[!rpkm_columns %in% expected_rpkm_columns]) {
          cat(paste0(col_name,"\n"))
        }
        cat("this will work with the shiny app, but the data for these cells will be hidden.\n\n")
        
      }
    }
    
  }
  
  if(passed_checks == 5) {
    if(verbose) { cat("The dataset has passed all essential checks! HOORAY!!\n") 
    } 
    return(TRUE)
  } else {
    if(verbose) { cat("The dataset has some problems that will prevent it from working with scrattch.\n") 
    } 
    return(FALSE)
    
  }
  
}


#'Write SQLite3 Database
#'
#'This function will write anno, data, and desc tables to a SQLite3 database.
#'
#'It can also be used to overwrite any one table in the database for updates.
#'
#'@param desc a data.frame with two columns: base, which correspondes to the beginning of
#'annotation columns; and label, which will be the label used to display related objects. Default: NULL
#'@param anno a data.frame with annotation values. This table should always have a sample_id column
#'that corresponds to the column names of the data table. Each annotation should have an _id, _label, and _color column.
#'@param data a data.frame with data values. This first column should be called "gene", and each subsequent column
#'should correspond to a single sample_id in the anno table. Default = NULL
#'@param file a string specifying the location of the output database file. Required.
#'@param overwrite a logical value indicating if existing tables in the database should be overwritten. Default = FALSE.
#'
#'@return Writes anno, desc, and data tables to the specified file.
#'
#'@examples
#'rpkm <- v1_data
#'anno <- v1_anno
#'desc <- v1_desc
#'
#'write_database(desc = desc, anno = anno, data = rpkm, file = "v1.db", overwrite=F)
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