#' Rboretum Tandem Plotter
#' Simple cowplot wrapper to plot two or more ggtree/ggplot objects side-by-side, or vertically. Plots can be passed individually, or as a list
#' @param vertical OPTIONAL: If TRUE, draw plots as vertical stack. [Default: FALSE, print plots side-by-side]
#' @return Side-by-side ggplot object
#' @export

tandemPlotter <- function(...,vertical){
  
  if(missing(vertical)){
    vertical=FALSE
  } else if(!is.logical(vertical)){
    stop("'vertical' must be TRUE (stack plots vertically) or FALSE (add plots as horizontal side-by-side")
  }
  
  if(has_error(silent=TRUE,expr=is.list(...))){
    plotList <- list(...)
  } else{
    if(!is.list(...)){
      stop("Argument(s) must either be multiple plots, or a single list of plots.")
    }
    else{
      plotList <- list(...)[[1]]
    }
  }
  
  if(length(plotList)<2){
    stop("Two or more plots needed for tandemPlotter.")
  }
  
  plotCheck <- purrr::map(plotList,Rboretum::isPlot) %>% unlist()
  
  if(!all(plotCheck)){
    stop("One argument passed is not of class ggplot or ggtree")
  }
  
  plot_count <- length(plotList)
  plot_step <- 1/plot_count
  
  return_plot <- ggdraw()
  
  start_coord <- 0

  if(vertical){
    for(i in plot_count:1){
      return_plot <- return_plot + draw_plot(plotList[[i]],x=0,y=start_coord,height = plot_step)
      start_coord <- start_coord + plot_step
    }
  } else{
    for(i in 1:plot_count){
      return_plot <- return_plot + draw_plot(plotList[[i]],x=start_coord,y=0,width = plot_step)
      start_coord <- start_coord + plot_step
    }
  }
  
  return(return_plot)
}