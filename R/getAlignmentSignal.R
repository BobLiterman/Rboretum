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
#'   \item multiPhylo object from which common species will be extracted; or
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
#' @param existing_signal OPTIONAL: Append these results to the output from getAlignmentSignal() run with the same species_info and a different alignment
#' @return Dataframe containing split pattern for each site, in each alignment, relative to the given set of taxa
#' @export

getAlignmentSignal <- function(alignment_path,species_info,use_gaps,alignment_name,prefix,suffix,exisiting_signal){
  
  # Ensure that a path and root taxa are provided as character vectors
  if(missing(alignment_path)){
    stop("No alignment file or directories indicated with 'alignment_path'")
  } else if(!is.character(alignment_path)){
    stop("'alignment_path' should be a character vector of file paths or the path to an alignment directory.")
  } else if(missing(species_info)){
    stop("No 'species_info' provided")
  } else if(!is.character(species_info) & !Rboretum::isPhylo(species_info) & !Rboretum::isMultiPhylo(species_info,check_three_taxa = TRUE)){
    stop("'species_info' should be a character vector of tip labels, a phylo object, or a multiPhylo where all trees share three taxa.")
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
  
  # Figure out how many files are being read in
  if(length(alignment_path)==1){
    
    isFile <- file.exists(alignment_path) & !dir.exists(alignment_path)
    isDir <- dir.exists(alignment_path) & !isFile
    
    if(isFile){ # 'alignment_path' points to a single valid file
      
      alignment_count <- 1
      alignment_path <- file_path_as_absolute(alignment_path)
      default_name <- basename(alignment_path)
      
    } else if(isDir){ # 'alignment_path' points to a valid directory
      
      if(has_error(silent=TRUE,list.files(path=alignment_path,pattern=align_regex,full.names = TRUE,include.dirs = FALSE))){
        stop("Can't process file fetch. Check path or regex?")
      } else{
        
        alignment_path <- list.files(path=alignment_path,pattern=align_regex,full.names = TRUE,include.dirs = FALSE)
        
        if(length(alignment_path)==0){
          stop("Directory found, but no files identified in 'alignment_path'. Check regex?")
        } else if(length(alignment_path)==1){
          alignment_count <- 1
          default_name <-  basename(alignment_path)
        } else{
          alignment_count <- length(alignment_path)
          default_name <- purrr::map(alignment_path,.f = function(x){basename(x)}) %>% unlist() %>% as.character()
        }
      }
    } else{ stop("'alignment_path' points to neither a valid file or directory.") }
    
  } else{ # 'alignment_path' is a list of file paths
    
    file_check <- purrr::map(.x = alignment_path,.f=function(x){ file.exists(x) & !dir.exists(x)}) %>% unlist() %>% all() # Check that all paths in 'alignment_path' point to valid files
    
    if(!file_check){
      stop("At least one file from 'alignment_path' points to an invalid path.")
    } else{
      alignment_path <- purrr::map(.x=alignment_path,.f=function(x){file_path_as_absolute(x)}) %>% unlist()
      alignment_count <- length(alignment_path)
      default_name <- purrr::map(alignment_path,.f = function(x){basename(x)}) %>% unlist() %>% as.character()
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
        print(paste(c("'use_gaps' (",length(use_gaps),") and number of alignments (",alignment_count,") do not match...using default setting of ignoring gaps..."),collapse = ''))
        gap_list <- rep(FALSE,alignment_count)
      } else{
        gap_list <- use_gaps
      }
    } 
  } else{
    stop("'use_gaps' must be logical (single value or vector)")
  }
  
  if(missing(alignment_name)){
    alignment_name <- default_name
  } else if(length(alignment_name) != alignment_count){
    print(paste(c("'alignment_names' (",length(alignment_names),") and number of alignments (",alignment_count,") do not match...using default names..."),collapse = ''))
    alignment_name <- default_name
  }
  
  if(any(duplicated(alignment_name))){
    stop("'alignment_name' cannot contain duplicate values")
  }
  
  if(alignment_count == 1){
    signal_df <- Rboretum::getAlignmentSignal_Worker(alignment_path,species_info,gap_list,alignment_name)
  } else if(alignment_count > 1){
    
    signal_df <- tibble(Alignment_Name=character(),Alignment_Position=integer(),Site_Pattern=character(),Gap=logical(),Singleton=logical(),Singleton_Taxa=character(),Non_Base_Taxa=character(),Non_Base_Count=integer(),Split_1=character(),Split_2=character(),Split_3=character(),Split_4=character(),Split_5=character())
    for(i in 1:alignment_count){
      temp_df <- Rboretum::getAlignmentSignal_Worker(alignment_path[i],species_info,gap_list[i],alignment_name[i])
      if(nrow(temp_df)>0){
        signal_df <- rbind(signal_df,temp_df)
      }
    }
  }
  
  if(nrow(signal_df)==0){
    stop("No data returned. Check alignments and taxa?")
  } else if(missing(exisiting_signal)){ # If  not appending results, return results
    return(signal_df)
  } else if(!Rboretum::isAlignmentSignal(exisiting_signal,species_info)){ # Ensure exisiting_signal contains information about the same taxa
      print("Object passed as 'exisiting_signal' is failing the test isAlignmentSignal(exisiting_signal,species_info)...Returning results without appending...")
      return(signal_df)
  } else if(any(alignment_name %in% existing_signal$Alignment_Name)){ # Ensure that alignment names are all unique
    print("Object passed as 'exisiting_signal' contains data on an alignment named identically one processed here...Returning results without appending...")
    return(signal_df)
  } else{
    return(rbind(exisiting_signal,signal_df))
  }
}