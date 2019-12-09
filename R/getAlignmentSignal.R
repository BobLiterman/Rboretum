#' Rboretum Alignment Signal Fetcher
#'
#' Given the path(s) to multiple alignments and a common list of taxa, this script returns site patterns for each site, in each alignment (after pruning if needed)
#' @param alignment_path Where to find alignment files. Options include:
#' \itemize{
#'   \item A character vector of one or more alignment file paths  (relative or absolute)
#'   \item A path to a single directory containing all alignment files (relative or absolute)
#' }
#' @param species_info List of taxa to analyze. Can be provided as:
#' \itemize{
#'   \item phylo object from which species will be extracted; or
#'   \item Character vector with 3+ taxon IDs
#' }
#' @param use_gaps OPTIONAL: Options include:
#' \itemize{
#'   \item FALSE (Default: Treat all gaps (-) in all alignments as missing data)
#'   \item TRUE (Treat all gaps in all alignments as potentially informative signal)
#'   \item Logical vecor indicating TRUE or FALSE for each alignment
#' }
#' @param alignment_name Chacter vector of names for each alignment. If missing or incomplete, the base filename is used
#' @param prefix OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file prefixes (e.g. all alignment files start with "Mafft_")
#' @param suffix OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file suffixes (e.g. all alignment files end with ".phy")
#' @return Dataframe containing split pattern for each site, in each alignment, relative to the given set of taxa
#' @export

getAlignmentSignal <- function(alignment_path,species_info,use_gaps,alignment_name,prefix,suffix){
  
  # Ensure that an alignment path and species info are provided
  if(missing(alignment_path)){
    stop("No alignment file or directories indicated with 'alignment_path'")
  } else if(!is.character(alignment_path)){
    stop("'alignment_path' should be a character vector of file paths or the path to an alignment directory.")
  } else if(missing(species_info)){
    stop("No 'species_info' provided")
  } else if(!is.character(species_info) & !Rboretum::isPhylo(species_info)){
    stop("'species_info' must be either a character vector of taxa, or a phylo object.")
  }
  
  # Create regex search pattern in case a directory is given

  if(missing(prefix)){
    prefix <- c()
  } else if(!is.character(prefix)){
    stop("'prefix' must be a character vector")
  } else{
    prefix <- unlist(purrr::map(.x=prefix,.f=function(x){paste(c("^",x),collapse = '')}))
    prefix <- paste(c("(",paste(prefix,collapse = "|"),")"),collapse = '')
  }
  
  if(missing(suffix)){
    suffix <- c()
  } else if(!is.character(suffix)){
    stop("'suffix' must be a character vector")
  } else{
    suffix <- unlist(purrr::map(.x=suffix,.f=function(x){ifelse(substr(x,start = 1,stop = 1)==".",paste(c("\\",x,"$"),collapse = ''),paste(c(x,"$"),collapse = ''))}))
    suffix <- paste(c("(",paste(suffix,collapse = "|"),")"),collapse = '')
  }
  
  if(length(prefix)==0 & length(suffix)==0){
    align_regex <- ''
  } else if(length(prefix)>0 & length(suffix)==0){
    align_regex <- prefix
  } else if(length(prefix)==0 & length(suffix)>0){
    align_regex <- suffix
  } else if(length(prefix)>0 & length(suffix)>0){
    align_regex <- paste(paste(c(prefix,"(.*)",suffix),collapse = ""))
  }
  
  # If 'alignment_path' is a single item...
  if(length(alignment_path)==1){
    isFile <- file.exists(alignment_path) & !dir.exists(alignment_path)
    isDir <- dir.exists(alignment_path) & !isFile
    
    if(!isFile & !isDir){
      stop("'alignment_path' doesn't point to a valid file or directory path.")
    } else if(isFile){ # 'alignment_path' points to a single alignment file
      alignment_path <- file_path_as_absolute(alignment_path)
      alignment_count <- 1
    } else if(isDir){ # 'alignment_path' points to a directory
      
      if(has_error(silent=TRUE,list.files(path=to_root,pattern=align_regex,full.names = TRUE,include.dirs = FALSE))){
        stop("Can't process file fetch. Check path or regex?")
      } else{
        
        alignment_path <- list.files(path=alignment_path,pattern=align_regex,full.names = TRUE,include.dirs = FALSE)
        
        if(length(alignment_path)==0){
          stop("No files found in 'alignment_path' using provided path and regex")
        } else{
          alignment_count <- length(alignment_path)
        }
      }
    }
  } else{ # If 'alignment_path' is a vector of multiple paths
    
    # Check if all files are valid
    path_check <- purrr::map(.x=alignment_path,.f=function(x){file.exists(x) & !dir.exists(x)}) %>% unlist() %>% all()
    
    if(path_check){
      alignment_count <- length(alignment_path)
    } else{
      stop("At lease one path in 'alignment_path' points to a non-existent file. Check working directory?")
    }
  }
  
  
  
  
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#   
#   # Get gap data
#   if(missing(use_gaps)){
#     gap_list <- rep(FALSE,alignment_count)
#   } else if(is.logical(use_gaps)){
#    if(length(use_gaps)==1){
#      if(use_gaps){
#        gap_list <- rep(TRUE,alignment_count)
#      }else{
#        gap_list <- rep(FALSE,alignment_count)
#      }
#    } else{
#      if(length(use_gaps)!=alignment_count){
#        stop("Arguments to 'alignment_paths' and 'use_gaps' differ in count.")
#      } else{
#        gap_list <- use_gaps
#      }
#    } 
#   } else{
#     stop("'use_gaps' must be logical (single value or vector)")
#   }
#   
#   if(missing(alignment_names)){
#     alignment_names <- c()
#     for(i in 1:alignment_count){
#       alignment_names <- c(alignment_names,file_path_sans_ext(basename(file_path_as_absolute(alignment_paths[i]))))
#     }
#   } else if(length(alignment_names) != alignment_count){
#     print("Not enough alignment names provided. Generating names from filenames...")
#     alignment_names <- c()
#     for(i in 1:alignment_count){
#       alignment_names <- c(alignment_names,file_path_sans_ext(basename(file_path_as_absolute(alignment_paths[i]))))
#     }
#   }
#   
#   if(any(duplicated(alignment_names))){
#     stop("'alignment_names' contains duplicated values. All names must be unique.")
#   }
#   
#   return_table <- tibble(Alignment_Name=character(),Alignment_Position=integer(),Site_Pattern=character(),Gap=logical(),Singleton=logical(),Singleton_Taxa=character(),Non_Base_Taxa=character(),Non_Base_Count=integer(),Split_1=character(),Split_2=character(),Split_3=character(),Split_4=character(),Split_5=character())
#   
#   for(i in 1:alignment_count){
#     new_df <- Rboretum::getAlignmentSignal(alignment_paths[i],species_info,gap_list[i],alignment_names[i])
#     if(nrow(new_df)==0){
#       print(alignment_paths[i])
#       print("The above alignment returned an error.")
#     } else{
#       return_table <- rbind(return_table,new_df)
#     }
#   }
#   
#   if(nrow(return_table)==0){
#     stop("All alignments returned errors in getAlignmentSignal.")
#   } else{
#     return(return_table)
#   }
# }
# 
