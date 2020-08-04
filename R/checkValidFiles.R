#' Rboretum File Checker
#'
#' This function takes one or more file paths, and returns TRUE if all paths point to valid files, and FALSE otherwise. If return_full_path = TRUE and all paths are valid, function returns absolute paths. If return_invalid = TRUE and there are invalid paths, function returns invalid paths.
#' @param alignment_path A character vector of one or more alignment file paths (relative or absolute)
#' @param return_full_path OPTIONAL: If TRUE and all paths are valid, function returns absolute paths [Default: FALSE]
#' @param return_invalid OPTIONAL: If TRUE and any paths are invalid, function returns invalid paths [Default: FALSE]
#' @return TRUE if all paths point to valid files, otherwise FALSE. If return_full_path = TRUE and all paths are valid, function returns absolute paths. If return_invalid = TRUE and there are invalid paths, function returns invalid paths.
#' @export

checkValidFiles <- function(alignment_path,return_full_path,return_invalid){
  
  # Ensure alignment_path is provided
  if(missing(alignment_path)){
    stop("'checkValidFiles' requires an 'alignment_path' argument...")
  } else if(!is.character(alignment_path)){
    stop("'alignment_path' should be a character vector...")
  }
  
  # Return logical or full paths?
  if(missing(return_full_path)){
    return_full_path <- FALSE
  } else if(!is.logical(return_full_path)){
    return_full_path <- FALSE
  } else if(length(return_full_path)!=1){
    return_full_path <- FALSE    
  }
  
  # Return invalid paths?
  if(missing(return_invalid)){
    return_invalid <- FALSE
  } else if(!is.logical(return_invalid)){
    return_invalid <- FALSE
  } else if(length(return_invalid)!=1){
    return_invalid <- FALSE    
  }
  
  # If one path is given...
  if(length(alignment_path)==1){
    
    if(file.exists(alignment_path) & !dir.exists(alignment_path)){ 
      
      if(return_full_path){
        return(file_path_as_absolute(alignment_path))
      } else{
        return(TRUE)
      }
      
    } else{
      
      if(return_invalid){
        return(alignment_path)
      } else{
        return(FALSE)
      }
    }
  } else{ # 'alignment_path' is a list of file paths
    
    file_check <- purrr::map(.x = alignment_path,.f=function(x){file.exists(x) & !dir.exists(x)}) %>% unlist() # Get T/F for each file
    
    # If all files are present, return TRUE or paths
    if(all(file_check)){
      if(return_full_path){
        file_paths <- purrr::map(.x = alignment_path,.f=function(x){file_path_as_absolute(x)}) %>% unlist() # Get full path for each file
        return(file_paths)
      } else{
        return(TRUE)
      }
    } else{ # If any paths point to invalid files...
      if(return_invalid){
        # Get invalid paths
        invalid_paths <- alignment_path[!file_check]
        return(invalid_paths)      
      } else{
        return(FALSE)
      }
    }
  }
}