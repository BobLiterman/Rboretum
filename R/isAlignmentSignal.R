#' Rboretum Alignment Signal Checker
#'
#' This function returns TRUE if the passed object is the output of getAlignmentSignal; otherwise, FALSE
#' @param test_object R object to check
#' @param species_info OPTIONAL: Check that passed signal contains information about a specific set of taxa. Set of 3+ species that can be passed as either:
#' \itemize{
#'   \item phylo object from which species will be extracted; or
#'   \item Character vector of desired taxa
#' }
#' @return TRUE if output of getAlignmentSignal(); otherwise, FALSE
#' @export

isAlignmentSignal <- function(test_object,species_info){
  
  if(!is.data.frame(test_object)){
    return(FALSE)
  } else if(nrow(test_object) == 0){
    return(FALSE)
  } else if(has_error(all(names(test_object) == c('Alignment_Name','Alignment_Position','Site_Pattern','Gap','Singleton','Singleton_Taxa','Non_Base_Taxa','Non_Base_Count','Split_1','Split_2','Split_3','Split_4','Split_5')))){
    return(FALSE)
  } else if(!all(names(test_object) == c('Alignment_Name','Alignment_Position','Site_Pattern','Gap','Singleton','Singleton_Taxa','Non_Base_Taxa','Non_Base_Count','Split_1','Split_2','Split_3','Split_4','Split_5'))){
    return(FALSE)
  }
  
  if(missing(species_info)){
    return(TRUE)
  } else{
    
    signal_taxa <- test_object %>%
      filter(!is.na(Split_1)) %>% 
      head(1) %>% 
      select(Singleton_Taxa,Non_Base_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5) %>% 
      select_if(~ !any(is.na(.))) %>% 
      unite(col = "Taxa",sep = ";") %>% 
      semiVector() %>% 
      sort()
    
    if(Rboretum::isPhylo(species_info)){
      spp_list = sort(species_info$tip.label)
    } else if(typeof(species_info)=="character" && length(species_info) > 3){
      spp_list = sort(species_info)
    } else { stop("'species_info' is not a phylo object or character vector 3+ species IDs") }  
    
    if(all(signal_taxa == species_info)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}
