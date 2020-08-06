#' Rboretum Alignment Signal Checker
#'
#' This function returns TRUE if the passed object is the output of getAlignmentSignal; otherwise, FALSE
#' @param test_object R object to check
#' @param species_info OPTIONAL: Check that passed signal contains information about a specific set of taxa. Set of 3+ species that can be passed as either:
#' \itemize{
#'   \item phylo object from which species will be extracted; or
#'   \item multiPhylo object from which common species will be extracted; or
#'   \item Character vector of desired taxa
#' }
#' @param return_taxa OPTIONAL: If TRUE, if object appears to be signal, return signal taxa rather than TRUE or FALSE [Default: FALSE, return logical]
#' @return TRUE if output of getAlignmentSignal(); otherwise, FALSE
#' @export

isAlignmentSignal <- function(test_object,species_info,return_taxa){
  
  # Ensure dataframe columns match expected and that data exists
  if(!is.data.frame(test_object)){
    return(FALSE)
  } else if(nrow(test_object) == 0){
    return(FALSE)
  } else if(has_error(silent=TRUE,expr=all(names(test_object) == c("Alignment_Position","All_Site_Bases","Site_Pattern","Non_Base_Taxa","Non_Base_Count","Singleton_Taxa","Singleton_Count","Gap_Taxa","Split_1","Split_2","Split_3","Split_4","Split_5","Alignment_Name")))){
    return(FALSE)
  } else if(!all(names(test_object) == c("Alignment_Position","All_Site_Bases","Site_Pattern","Non_Base_Taxa","Non_Base_Count","Singleton_Taxa","Singleton_Count","Gap_Taxa","Split_1","Split_2","Split_3","Split_4","Split_5","Alignment_Name"))){
    return(FALSE)
  }
  
  # Return logical or taxa list?
  if(missing(return_taxa)){
    return_taxa <- FALSE
  } else if(!is.logical(return_taxa)){
    return_taxa <- FALSE
  }
  
  # Get species info for each alignment
  
  alignments <- unique(test_object$Alignment_Name)
  signal_taxa <- c()
  
  for(i in alignments){
    
    alignment_taxa <- test_object %>%
      filter(Alignment_Name == i) %>%
      filter(!is.na(Split_1)) %>% 
      head(1) %>% 
      select(Singleton_Taxa,Non_Base_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5) %>% 
      select_if(~ !any(is.na(.))) %>% 
      unite(col = "Taxa",sep = ";") %>% 
      pull(Taxa) %>%
      semiSorter()
    
    signal_taxa <- c(signal_taxa,alignment_taxa)
  }
  
  if(length(unique(signal_taxa))!=1){ # If entries for different alignments return different species lists, return FALSE
    return(FALSE)
  }
  
  signal_taxa <-semiVector(alignment_taxa) %>% unlist() %>% naturalsort()
  
  if(return_taxa){
    return(signal_taxa)
  }
  
  if(missing(species_info)){ # If no species tests are requested,  return TRUE
      return(TRUE)
  } else{ # Ensure species from 'test_object' match species from 'species_info'
    
    if(Rboretum::isPhylo(species_info)){
      spp_list = species_info$tip.label %>% naturalsort()
    } else if(Rboretum::isMultiPhylo(species_info,check_three_taxa = TRUE)){
      spp_list = Rboretum::getSharedTaxa(species_info) %>% naturalsort()
    } else if(is.character(species_info) & length(species_info) > 3 & !any(duplicated(species_info))){
      spp_list = species_info %>% naturalsort()
    } else { stop("'species_info' is not a valid phylo object, multiPhylo object where all trees share three taxa, or a character vector with 3+ unique species IDs") }
    
    if(all(spp_list %in% signal_taxa) & all(signal_taxa %in% spp_list)){
        return(TRUE)
      } else{
        return(FALSE)
      }
  }
}
