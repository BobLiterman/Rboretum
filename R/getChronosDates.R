#' Get Chronos Node Dates
#'
#' This function takes a tree and root node age calibrations, and converts branch lengths into divergence times using Chronos in ape
#' @param tree Phylo object
#' @param iter_count Number of iterations for chronos
#' @param min_root  Minimum divergence time estimate for root node
#' @param max_root  Maximum divergence time estimate for root node
#' @return Ultrametric chronogram
#' @export
#' @examples
#'
#' # If root node is forced at 50MYA, with 2000 iterations
#' getChronosDates(tree,2000,50,50)
#'
#' # If root node is estimated between 30 and 50 MYA, with 1000 iterations
#' getChronosDates(tree,1000,30,50)
#'

getChronosDates <- function(tree,iter_count,min_root,max_root){

  if(missing(min_root) | missing(max_root)){
    stop("No minimum or maximum age given for root. Use 1 for each argument if relative dating. ie. getChronosDates(tree,iterations,1,1)")
  }

  edge_list <- ape::compute.brlen(tree,1)$edge.length
  branch_list <- ape::branching.times(compute.brlen(tree,1))

  chronos_cal <- ape::makeChronosCalib(tree,age.min = min_root,age.max = max_root)

  for(i in 1:iter_count){
    chronos_iter <- ape::chronos(tree,calibration = chronos_cal)
    edge_list <- rbind(edge_list, chronos_iter$edge.length)
    branch_list <- rbind(branch_list, ape::branching.times(chronos_iter))
  }

  edge_list <- edge_list[-1,]
  branch_list <- branch_list[-1,]

  branches <- colnames(branch_list)

  median_dates <- c()
  for(i in 1:ncol(branch_list)){
    median_dates[colnames(branch_list)[i]] <- median(branch_list[,i])
  }

  timetree <- ape::compute.brlen(tree,1)
  median_edge <- c()

  for(i in 1:ncol(edge_list)){
    median_edge<-c(median_edge,median(edge_list[,i]))
  }

  timetree$edge.length <- median_edge
  if(!ape::is.ultrametric(timetree)){
    timetree <- phytools::force.ultrametric(timetree,"extend")
  }
  return(timetree)
}
