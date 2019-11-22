#' Rboretum Simple Tandem Plotter
#' Simple cowplot wrapper to plot two or more ggtree/ggplot objects side by side
#' @return Side-by-side ggplot object
#' @export
#' @examples
#' tandem.treePlot(ggtree1,ggtree2)

tandem.treePlot <- function(...){
  
  plotList <- list(...)
  
  if(length(plotList)<2){
    stop("Two or more plots needed for tandem.treePlot.")
  }
  
  if(any(map(plotList,Rboretum::is.plot))){
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