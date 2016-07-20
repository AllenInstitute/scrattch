## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(dplyr)
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  
  names(data)[names(data) == measurevar] <- "measurevar"
  
  datac <- data %>%
    select(one_of(groupvars,"measurevar")) %>%
    filter(ifelse(na.rm == T, !is.na(measurevar), T)) %>%
    group_by_(c(groupvars)) %>%
    summarise(N = n(),
              median = median(measurevar),
              mean = mean(measurevar),
              max = max(measurevar),
              sd = sd(measurevar),
              q25 = as.numeric(quantile(measurevar, 0.25)),
              q75 = as.numeric(quantile(measurevar, 0.75))) %>%
    mutate(se = sd/sqrt(N)) %>%
    mutate(ci = se * qt(conf.interval/2 + 0.5, N-1))
  
  
  datac <- datac %>%
    mutate(xpos = 1:n())

  return(datac)
}

sci.label <- function(vec) {
  labels <- character()
  for(i in 1:length(vec)) {
    x <- round(vec[i],0)
    first <- round(x/(10^(nchar(x)-1)),1)
    if(first %% 1 == 0) {
      first <- paste0(first,".0")
    }
    label <- first
#    label <- paste0(first,"E",nchar(x)-1)
    labels <- c(labels,label)
  }
  return(labels)
}
