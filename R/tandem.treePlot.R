#' Rboretum Simple Tandem Plotter
#' Simple cowplot wrapper to plot two ggtree/ggplot objects side by side
#' @param plot1 ggtree/ggplot object to be plotted on the left
#' @param plot2  ggtree/ggplot object to be plotted on the right
#' @return Side-by-side ggplot object
#' @export
#' @examples
#' tandem.treePlot(ggtree1,ggtree2)

tandem.treePlot <- function(plot1,plot2){
  
  if(has_error(unlist(attributes(plot1)$class)) | has_error(unlist(attributes(plot2)$class))){ 
    stop("'plot1' and 'plot2' should be plottable objects (e.g. ggtree, ggplot, etc.")
  }
  
  plot1_class <- unlist(attributes(plot1)$class)
  plot2_class <- unlist(attributes(plot2)$class)
  
  if(!any(c('ggtree','ggplot') %in% plot1_class)){
    stop("'plot1' not a plottable object (ggtree or ggplot)")
  }
  
  if(!any(c('ggtree','ggplot') %in% plot2_class)){
    stop("'plot2' not a plottable object (ggtree or ggplot)")
  }
  
  return_plot <- ggdraw() + draw_plot(plot1, x = 0, y = 0, width = .5) + draw_plot(plot2, x = .5, y = 0, width = .5)
  
  return(return_plot)
}