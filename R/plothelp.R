#' Convert font sizes in pt to mm
#' @param pt A numeric font size in pt.
#' @return A numeric font size in mm.
#' 
#' @examples
#' pt2mm(12)
#' 
#' ggplot(mtcars) +
#'   geom_text(aes(x = mpg, y = wt, label = rownames(mtcars)),
#'             size = pt2mm(7))
pt2mm <- function(pt) {
  mm <- pt / 2.834645669
  return(mm)
}

#' Convert integers to scientific notation labels
#' 
#' @param in_num a numeric vector
#' @param sig_figs a number indicating how many significant figures should be displayed.
#' @return a character vector with numeric values reformatted in 1.2E3 format
#' 
#' @examples
#' my_numbers <- c(100,15.359,32687,.000468)
#' 
#' sci_label(my_numbers)
#' 
#' sci_label(my_numbers,sig_figs=3)
sci_label <- function(in_num,sig_figs=2) {
  labels <- character()
  for(i in 1:length(in_num)) {
    x <- in_num[i]
    if(x == 0) {
      first <- paste0("0",".",paste0(rep("0",sig_figs-1),collapse=""))
    } else if(log10(x) %% 1 == 0) {
      first <- substr(x,1,1)
      if(sig_figs > 1) {
        first <- paste0(first,".",paste0(rep("0",sig_figs-1),collapse=""))
      }
    } else {
      first <- round(x/(10^floor(log10(x))),sig_figs-1)
    }
    if(x == 0) {
      label <- paste0(first,"E0")
    } else {
      label <- paste0(first,"E",floor(log10(x)))
    }
    labels <- c(labels,label)
  }
  return(labels)
}

#' Remove the X-axis (and most other margins)
#' 
#' Makes plots more suitable for use with Illustrator by removing most margins
#' and the X-axis (which is usually replaced by something else in my plots).
#' 
#' To further remove the space below the x-axis, use labs(x = NULL)
#' 
#' Based on theme_classic() from ggplot2.
#' 
#' @examples
#' ggplot(mtcars) +
#'  geom_point(aes(x = mpg, y = wt)) +
#'  theme_no_x_margin() +
#'  labs(x = NULL)
theme_no_x <- function(base_size = 12, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(plot.margin = unit(c(rep(0,4)),"line"),
          axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.margin = unit(0,"cm"),
          axis.ticks.x = element_blank())
}

#' Split a character string by commas, spaces, and tabs
#' 
#' @param in_string a chr string containing commas, spaces, and/or tabs
#' @return a character vector with each object separated by any combination of commas, spaces, and/or tabs
#' 
#' @examples
#' test <- "Hspa8, Scnn1a,Rbp4    Ptgs2"
#' split_chr(test)
split_cst <- function(in_string) {
  out_chr <- strsplit(in_string,"[, \t]+")[[1]]
  return(out_chr)
}

#' Convert the case of objects in a character vector to Title Case
#' 
#' @param in_chr a character vector
#' @return a character vector with Each Object In Title Case
#' 
#' @examples
#' test <- c("hspa8","scnn1a","fhqwghads")
#' title_case(test)
title_case <- function(in_chr) {
    lower <- tolower(in_chr)
    s <- strsplit(lower, " ")
    result <- paste(toupper(substring(s, 1,1)), substring(s, 2), sep="")
    return(result)
}

#' Convert the case of Riken genes no matter input case
#' 
#' @param in_chr a character vector of Riken gene symbols
#' @return a character vector with correck Riken capitalization
#' 
#' @examples
#' test <- c("6330527o06RiK","A930038C07RIK","a330070k13rik")
#' riken_case(test)
riken_case <- function(in_chr) {
  upper <- toupper(in_chr)
  result <- sub("RIK","Rik",upper)
  return(result)
}

#' Correct case of mouse gene symbols
#' 
#' @param genes a character vector containing gene symbols
#' @return a character vector containing genes with proper case
#' 
#' @examples
#' test <- c("a330070k13rik","HDAC1","sncg","SeRPinB11")
#' fix_mouse_genes(test)
fix_mouse_genes <- function(genes) {
  genes[grepl("[Rr]ik$",genes)] <- riken_case(genes[grepl("[Rr]ik$",genes)])
  genes[!grepl("[Rr]ik$",genes)] <- title_case(genes[!grepl("[Rr]ik$",genes)])
  
  return(genes)
}

#' Convert a character object to a format compatible with SQL SELECT
#' 
#' @param in_chr a character vector
#' @return a character vector consisting of a single string with each object separated by 
#' escaped quotes, bracketted by parens.
#' 
#' @examples
#' test <- c("Hspa8","Cux2","Snap25","Rspo2")
#' chr_to_sql(test)
#' 
#' sql_query <- paste0("SELECT * FROM my_table WHERE gene IN ",chr_to_sql(test))
#' print(sql_query)
chr_to_sql <- function(in_chr) {
  result <- paste0("(",paste0("\"",in_chr,"\"",collapse=","),")")
  return(result)
}

#' Evaluate a character string to a numeric vector
chr_to_num <- function(in_chr) {
  result <- eval(parse(text=paste0("round(c(",in_chr,"),0)")))
  return(result)
}

get_internal_data <- function(genes,grouping,clusters) {
  data <- scrattch::v1_data
  all.anno <- scrattch::v1_anno
  
  cluster_order <- data.frame(clusters=clusters) %>%
    mutate(cluster_x=1:n())
  
  data <- data %>%
    filter(gene %in% genes)
  
  # Reformat the retrieved data
  row.names(data) <- data[,1]
  data <- data %>% 
    select(-1) %>% 
    t() %>% 
    as.data.frame()
  
  data <- data %>%
    mutate(sample_id=row.names(data)) %>%
    select(one_of(c("sample_id",genes)))
  
  genes[genes == "9630013A20Rik"] <- "X9630013A20Rik"
  names(data)[names(data) == "9630013A20Rik"] <- "X9630013A20Rik"
  
  # Filter and order the rows
  data <- left_join(data,all.anno,by="sample_id") %>%
    rename_("plot_id" = paste0(grouping,"_id"),
            "plot_label" = paste0(grouping,"_label"),
            "plot_color" = paste0(grouping,"_color")) %>%
    filter(plot_id %in% clusters) %>%
    left_join(cluster_order,by=c("plot_id"="clusters")) %>%
    arrange(cluster_x) %>%
    mutate(xpos=1:n()) %>%
    select(-plot_id) %>%
    rename(plot_id=cluster_x)
  
  return(data)
}

get_db_data <- function(data_source,genes,grouping,clusters) {
  library(DBI)
  library(RSQLite)
  
  cluster_order <- data.frame(clusters=clusters) %>%
    mutate(cluster_x=1:n())
  
  con <- dbConnect(RSQLite::SQLite(),data_source)
  get <- "SELECT * FROM anno;"
  res <- dbSendQuery(con,get)
  all.anno <- dbFetch(res,n=-1)
  dbClearResult(res)
  getgenes <- paste("\"",genes,"\"",collapse=",",sep="")
  get <- paste("SELECT * from data WHERE gene IN (",getgenes,")",sep="")
  res <- dbSendQuery(con,get)
  data <- dbFetch(res,n=-1)
  dbClearResult(res)
  row.names(data) <- data[,1]
  data <- data %>% 
    select(-1) %>% 
    t() %>% 
    as.data.frame()
  
  data <- data %>%
    mutate(sample_id=row.names(data)) %>%
    select(one_of(c("sample_id",genes)))
  
  # Filter and order the rows
  data <- left_join(data,all.anno,by="sample_id") %>%
    rename_("plot_id" = paste0(grouping,"_id"),
            "plot_label" = paste0(grouping,"_label"),
            "plot_color" = paste0(grouping,"_color")) %>%
    filter(plot_id %in% clusters) %>%
    left_join(cluster_order,by=c("plot_id"="clusters")) %>%
    arrange(cluster_x) %>%
    mutate(xpos=1:n()) %>%
    select(-plot_id) %>%
    rename_("plot_id" = "cluster_x")
  
  return(data)
}