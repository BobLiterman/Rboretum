#' Rboretum Clade Label Converter
#'
#' This function takes a phylo, multiPhylo, a lits, or a character vector and changes names based on a supplied dataframe
#' @param to_convert R object containing taxon labels to be replaced. Options include: 
#' \itemize{
#'   \item A phylo object; or,
#'   \item A multiPhylo object; or,
#'   \item A character vector of taxon IDs [e.g. c('Species1','Species2',etc.)]; or,
#'   \item A list of character vectors containing taxon IDs [e.g. list(A=c('Species1','Species2'),B=c('Species_3','Species_4'))]; or,
#'   \item A character vector of semicolon-separated clades [e.g. c('Species_1;Species_2','Species_3;Species_4')]
#' }
#' @param name_df Dataframe with name equivalencies
#' @param from Column name with all current taxon tip IDs (if missing, first column of name_df used as default)
#' @param to Column name with desired IDs (if missing, second column of name_df used as default)
#' @return Matching object with new names in place
#' @export
#' 
convertLabels <- function(to_convert,name_df,from,to){
  
  # Ensure input is valid
  if(!Rboretum::isPhylo(to_convert) & !Rboretum::isMultiPhylo(to_convert) & !is.list(to_convert) & !is.character(to_convert)){
    stop("'to_convert' must be a phylo, multiPhylo, list, or character vector")
  }
  
  # Process name dataframe
  if(ncol(name_df)<2){
    stop("'name_df' must at least have two columns [current IDs, new IDs]")
  }
  
  # Get 'from' column [current IDs]
  if(missing(from)){
    current_ids <- pull(name_df,1) %>% as.character()
  } else{ 
    if(!as.character(from) %in% names(name_df)){
      stop("'from' column not found in name_df")
    } else{ current_ids <- pull(name_df,from) %>% as.character() }
  }
  
  # Ensure 'from' column IDs are unique
  if(any(duplicated(current_ids))){
    stop("Specified 'from' column contains duplicate IDs")
  }
  
  # Get 'to' column [desired IDs]
  if(missing(to)){
    new_ids <- pull(name_df,2) %>% as.character()
  } else{ 
    if(!as.character(to) %in% names(name_df)){
      stop("'to' column not found in name_df")
    } else{ new_ids <- pull(name_df,to) %>% as.character() }
  }
  
  # Ensure 'to' column IDs are unique
  if(any(duplicated(new_ids))){
    warning("Specified 'to' column contains duplicate IDs")
  }
  
  # Ensure all tips are in 'from' dataset
  if(Rboretum::isPhylo(to_convert)){
    
    convert_type <- 'phylo'
    old_taxa <- to_convert$tip.label # Get taxa from tree
    
  } else if(Rboretum::isMultiPhylo(to_convert)){
    
    tree_count <- length(to_convert)
    convert_type <- 'multiPhylo'
    old_taxa <- purrr::map(.x=to_convert,.f=function(x){x$tip.label}) %>% unlist() %>% as.character() %>% unique()
    
  } else if(is.list(to_convert)){
    
    convert_type <- 'list'
    old_taxa <- to_convert %>% unlist() %>% as.character() %>% unique()
    
  } else if(is.character(to_convert)){
    
    if(any(str_detect(to_convert,";"))){ # If character, check for semicolon separated clades
      
      if(!all(str_detect(to_convert,";"))){ stop("'to_convert' appears to have a mix of taxon IDs and semicolon separated clades") } 
      
      else{
        
        convert_type <- 'clade'
        old_taxa <- purrr::map(.x=to_convert,.f= function(x){ semiVector(x) }) %>% unlist() %>% as.character() %>% unique() # Break up clades into taxon IDs
      }
    } else{
      
      convert_type <- 'taxa'
      old_taxa <- to_convert %>% unique()
    }
  }
  
  if(!all(old_taxa %in% current_ids)){
    stop("Specified 'from' column does not contain all IDs from object")
  } 
  
  if(convert_type == 'phylo'){ # If to_convert is a phylo, replace tips...
    
    new_id_list <- c()
    
    for(old_id in to_convert$tip.label){
      new_id_list <- c(new_id_list,new_ids[match(old_id,current_ids)]) # Get matching ID
    }
    
    to_convert$tip.label <- new_id_list # Replace tip labels and return
    return(to_convert)
  }
  
  if(convert_type == 'multiPhylo'){ # If to_convert is a multiPhylo, repeat process above for each tree in multiPhylo...
    
    for(i in 1:tree_count){ 
      
      new_id_list <- c()
      
      for(old_id in to_convert[[i]]$tip.label){
        new_id_list <- c(new_id_list,new_ids[match(old_id,current_ids)])
      }
      
      to_convert[[i]]$tip.label <- new_id_list
    }
    
    return(to_convert)
    
  } else if(convert_type == 'list'){ # If converting a list of raw taxon IDs
    
    new_list <- list()
    
    for(i in 1:length(to_convert)){
      
      new_id_list <- c()
      
      for(old_id in to_convert[[i]]){
        new_id_list <- c(new_id_list,new_ids[match(old_id,current_ids)])
      }
      
      new_list[[i]] <- new_id_list
    }
    
    if(!is.null(names(to_convert))){
      names(new_list) <- names(to_convert)
    }
    
    return(new_list)
    
  } else if(convert_type == 'clade'){ # If converting semicolon-separated clades, return sorted, converted clade
    
    clade_list <- lapply(to_convert,function(x) Rboretum::semiVector(x))
    
    for(i in 1:length(clade_list)){
      
      new_id_list <- c()
      
      for(old_id in clade_list[[i]]){
        new_id_list <- c(new_id_list,new_ids[match(old_id,current_ids)])
      }
      
      
      new_id_list <- new_id_list %>% naturalsort() %>% paste(collapse = ";")
      
      clade_list[[i]] <- new_id_list
    }
    
    return(unlist(clade_list))
    
  } else if(convert_type == 'taxa'){ # If converting a vector of raw taxon IDs
    
    new_id_list <- c()
    
    for(old_id in to_convert){
      new_id_list <- c(new_id_list,new_ids[match(old_id,current_ids)])
    }
    
    return(new_id_list)
  }
}