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