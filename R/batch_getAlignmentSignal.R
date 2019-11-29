#' Rboretum Batch Alignment Signal Fetcher
#'
#' Given the paths to multiple alignments and a common list of taxa, this script returns site patterns for each site, in each alignment (after pruning if needed)
#' @param alignment_paths Vector of paths to alignment files (absolute or relative)
#' @param species_info Set of 3+ species that can be passed as either:
#' \itemize{
#'   \item phylo object from which species will be extracted; or
#'   \item Character vector of desired taxa
#' }
#' @param use_gaps
#' \itemize{
#'   \item FALSE (Default: Treat all gaps (-) in all alignments as missing data)
#'   \item TRUE (Treat all gaps in all alignments as potentially informative signal)
#'   \item Logical vecor indicating TRUE or FALSE for each alignment
#' }
#' @param alignment_names Chacter vector of names for each alignment. If missing or incomplete, the base filename is used
#' @return Dataframe containing split pattern for each site, in each alignment, relative to the given set of taxa
#' @export

batch_getAlignmentSignal <- function(alignment_paths,species_info,use_gaps,alignment_names){
  
  if(length(alignment_paths)==1){
    stop("Single alignment path passed. Use getAlignmentSignal() for single alignments")
  }
  
  alignment_count <- length(alignment_paths)
  
  # Check paths
  for(i in 1:alignment_count){
    if(has_error(file_path_as_absolute(alignment_paths[i]))){
      print(alignment_paths[i])
      stop("The above path does not point to a valid file.")
    }
  }
  
  # Get gap data
  if(missing(use_gaps)){
    gap_list <- rep(FALSE,alignment_count)
  } else if(is.logical(use_gaps)){
   if(length(use_gaps)==1){
     if(use_gaps){
       gap_list <- rep(TRUE,alignment_count)
     }else{
       gap_list <- rep(FALSE,alignment_count)
     }
   } else{
     if(length(use_gaps)!=alignment_count){
       stop("Arguments to 'alignment_paths' and 'use_gaps' differ in count.")
     } else{
       gap_list <- use_gaps
     }
   } 
  } else{
    stop("'use_gaps' must be logical (single value or vector)")
  }
  
  if(missing(alignment_names)){
    alignment_names <- c()
    for(i in 1:alignment_count){
      alignment_names <- c(alignment_names,file_path_sans_ext(basename(file_path_as_absolute(alignment_paths[i]))))
    }
  } else if(length(alignment_names) != alignment_count){
    print("Not enough alignment names provided. Generating names from filenames...")
    alignment_names <- c()
    for(i in 1:alignment_count){
      alignment_names <- c(alignment_names,file_path_sans_ext(basename(file_path_as_absolute(alignment_paths[i]))))
    }
  }
  
  if(any(duplicated(alignment_names))){
    stop("'alignment_names' contains duplicated values. All names must be unique.")
  }
  
  return_table <- tibble(Alignment_Name=character(),Alignment_Position=integer(),Site_Pattern=character(),Gap=logical(),Singleton=logical(),Singleton_Taxa=character(),Non_Base_Taxa=character(),Non_Base_Count=integer(),Split_1=character(),Split_2=character(),Split_3=character(),Split_4=character(),Split_5=character())
  
  for(i in 1:alignment_count){
    new_df <- Rboretum::getAlignmentSignal(alignment_paths[i],species_info,gap_list[i],alignment_names[i])
    if(nrow(new_df)==0){
      print(alignment_paths[i])
      print("The above alignment returned an error.")
    } else{
      return_table <- rbind(return_table,new_df)
    }
  }
  
  if(nrow(return_table)==0){
    stop("All alignments returned errors in getAlignmentSignal.")
  } else{
    return(return_table)
  }
}