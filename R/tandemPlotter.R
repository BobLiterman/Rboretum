#' Rboretum Tandem Plotter
#' Simple cowplot wrapper to plot two or more ggtree/ggplot objects side by side. Plots can be passed individually, or as a list
#' @return Side-by-side ggplot object
#' @export

tandemPlotter <- function(...){
  
  if(has_error(is.list(...),silent = TRUE)){
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
    stop("Two or more plots needed for tandem.treePlot.")
  }
  
  plotCheck <- purrr::map(plotList,Rboretum::isPlot) %>% unlist()
  
  if(!all(plotCheck)){
    stop("One argument passed is not of class ggplot or ggtree")
    }
  
  plot_count <- length(plotList)
  plot_step <- 1/plot_count
  
  return_plot <- ggdraw()
  
  start_x <- 0
  
  for(i in 1:plot_count){
    return_plot <- return_plot + draw_plot(plotList[[i]],x=start_x,y=0,width = plot_step)
    start_x <- start_x + plot_step
  }
  
  return(return_plot)
}