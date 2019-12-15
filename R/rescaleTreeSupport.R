#' Rboretum Tree Support Rescaler
#' 
#' This function takes the output of getTreeSupport, and returns a numeric vector of total, summed support that has been rescaled
#' @param tree_support Output of getTreeSupport
#' @param scale Scaling factor. Options include:
#' \itemize{
#'   \item Single numeric value [return a repeated vector of this value (e.g. to set all geom_nodepoint to this size)]
#'   \item Two-item numeric vector [return summation values that have been rescaled to fit between c(min,max) (e.g. rescale millions of sites to range between c(1,20))]
#'   \item "log" [return log-transformed totals]
#' }
#' @param return_total OPTIONAL: If TRUE, return raw sum total support values for each clade
#' @return Numeric vector the length of the number of clades, and corresponding to rescaled total support values
#' @export

rescaleTreeSupport <- function(tree_support,scale,return_total){
  
  if(missing(tree_support)){
    stop("'tree_support' is required.")
  } 
  
  if(missing(return_total)){
    return_total <- FALSE
  } else if(!is.logical(return_total)){
    return_total <- FALSE
  }
  
  if(missing(scale)){
    if(!return_total){
      stop("'scale' is required if not returning raw totals.")
    }
  } else if(!is.data.frame(tree_support)){
    stop("'tree_support' is not a dataframe")
  } else if(names(tree_support)[1] != 'Clade'){
    stop("The first column of getTreeSupport output is called 'Clade")
  } else if(ncol(tree_support)==1){
    stop("No data columns.")
  } else{
    
    data_columns <- 2:ncol(tree_support)
    
    if(has_error(silent=TRUE,expr = tree_support[data_columns] %>% rowSums() )){
      stop("Cannot perform rowSums on second through the last columns. Are all columns integer/numeric?")
    } else{
      clade_totals <- tree_support[data_columns] %>% rowSums() # Get row support (total support summed across datasets)
      
      if(return_total){
        return(as.numeric(clade_totals))
      }
    }
  }
  
  new_df <- data.frame(Sort = as.integer(1:length(tree_support$Clade)),Clade=as.character(tree_support$Clade),Raw_Total_Support=as.numeric(clade_totals)) # Create dataframe
  
  zero_support <- new_df %>% filter(Raw_Total_Support == 0) # Filter out clades with no support 
  has_support <- new_df %>% filter(Raw_Total_Support > 0)
  
  if(nrow(has_support)==0){
    stop("No data found")
  }
  
  # Handle 0 support nodes
  if(nrow(zero_support)>0){
    zero_support <- zero_support %>%
      mutate(Scaled_Support = 0)
    add_zero <- TRUE
  } else{
    add_zero <- FALSE
  }
  
  if(is.character(scale)){ # If 'scale' is a chaaracter, it should be 'log'
    if(length(scale)!=1){
      stop("'scale' should be 'log' if a character is given")
    } else if(scale != 'log'){
      stop("'scale' should be 'log' if a character is given")
    } else{
      has_support <- has_support %>%
        mutate(Scaled_Support = log(Raw_Total_Support))
    }
  } else if(is.numeric(scale)){ # If 'scale' is a number...
    
    if(length(scale)==1){
      has_support <- has_support %>%
        mutate(Scaled_Support = as.numeric(scale))
      
    } else if(length(scale)==2){
      
      min <- scale[1]
      max <- scale[2]
      
      if(!(min < max)){
        stop("'scale' values should be given in terms of c(min,max)")
      }
      
      raw_values <- has_support %>% pull(Raw_Total_Support) %>% as.numeric()
      new_values <- scales::rescale(raw_values,to=scale)
      has_support$Scaled_Support <- new_values
      
    } else{
      stop("If 'scale' is numeric, it should be one or two items long.")
    }
  }
  
  if(add_zero){
    
    scaled_values <- rbind(zero_support,has_support) %>% arrange(Sort) %>% pull(Scaled_Support) %>% as.numeric()
    
  } else{
    
    scaled_values <- has_support %>% pull(Scaled_Support) %>% as.numeric()
    
  }
  
  return(scaled_values)
  
}