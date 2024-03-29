#' Rboretum Alignment Pattern Fetcher
#'
#' Given the path(s) to one or more alignments (and an optional list of taxa), this script returns the site variation pattern for each site (e.g. invariant, biallelic, etc.)
#' @param alignment_path Where to find alignment files. Options include:
#' \itemize{
#'   \item A character vector of one or more alignment file paths  (relative or absolute)
#'   \item A path to a single directory containing all alignment files (relative or absolute)
#' }
#' @param species_info OPTIONAL: List of taxa to analyze [Default: Process alignment using all taxa if solo, all shared in multiple]. Can be provided as:
#' \itemize{
#'   \item phylo object from which species will be extracted; or
#'   \item multiPhylo object from which common species will be extracted; or
#'   \item Character vector with 3+ taxon IDs
#'   \item Semicolon-separated list of taxon IDs
#' }
#' @param use_gaps OPTIONAL: If FALSE, treat gaps as missing data (like N) [Default: TRUE, treat gaps as indel data]
#' @param alignment_name OPTIONAL: Chacter vector of names for each alignment. If missing or incomplete, the base filename of the alignment is used
#' @param prefix OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file prefixes (e.g. all alignment files start with "Mafft_")
#' @param suffix OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file suffixes (e.g. all alignment files end with ".phy")
#' @return Dataframe containing variation pattern each site in the alignment
#' @export

getAlignmentPatterns <- function(alignment_path,species_info,use_gaps,alignment_name,prefix,suffix){
  
  # Ensure that a path and root taxa are provided as character vectors
  if(missing(alignment_path)){
    stop("No alignment file or directories indicated with 'alignment_path'")
  } else if(!is.character(alignment_path)){
    stop("'alignment_path' should be a character vector of file paths or the path to an alignment directory.")
  } 
  
  # Create regex search pattern in case a directory is given
  if(missing(prefix)){
    prefix <- c()
  } else if(is.null(prefix)){
    prefix <- c()    
  } else if(!is.character(prefix)){
    stop("'prefix' must be a character vector")
  } else{
    prefix <- unlist(purrr::map(.x=prefix,.f=function(x){paste(c("^",x),collapse = '')}))
    prefix <- paste(c("(",paste(prefix,collapse = "|"),")"),collapse = '')
  }
  
  if(missing(suffix)){
    suffix <- c()
  } else if(is.null(suffix)){
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
  
  isFile <- file.exists(alignment_path) & !dir.exists(alignment_path)
  isDir <- dir.exists(alignment_path) & !isFile
  
  # Ensure files all exist
  if(length(alignment_path)==1){
    
    if(Rboretum::checkValidFiles(alignment_path)){ # 'alignment_path' points to a single valid file
      
      alignment_count <- 1
      alignment_path <- Rboretum::checkValidFiles(alignment_path,return_full_path = TRUE)
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
    
    file_check <- Rboretum::checkValidFiles(alignment_path) # Check that all paths in 'alignment_path' point to valid files
    
    if(!file_check){
      invalid_paths <- Rboretum::checkValidFiles(alignment_path,return_invalid = TRUE)
      print(invalid_paths)
      stop("The above paths from 'alignment_path' do not point to a valid file...")
    } else{
      alignment_path <- Rboretum::checkValidFiles(alignment_path,return_full_path = TRUE)
      alignment_count <- length(alignment_path)
      default_name <- purrr::map(alignment_path,.f = function(x){basename(x)}) %>% unlist() %>% as.character()
    }
  }
  
  # Set alignment names if requested, otherwise use default basenames
  if(missing(alignment_name)){
    alignment_name <- default_name
  } else if(length(alignment_name) != alignment_count){
    print(paste(c("'alignment_names' (",length(alignment_name),") and number of alignments (",alignment_count,") do not match...using default names..."),collapse = ''))
    alignment_name <- default_name
  }
  
  if(any(duplicated(alignment_name))){
    stop("'alignment_name' cannot contain duplicate values")
  }
  
  # Establish species of interest
  
  # If no species_info is provided, process alignments with taxa shared among all alignments
  if(missing(species_info)){
    if(alignment_count == 1){
      species_info = Rboretum::getAlignmentSpecies(alignment_path)
    } else{
      all_species <- purrr::map(.x=alignment_path,.f=function(x){Rboretum::semiVector(Rboretum::getAlignmentSpecies(x))})
      species_info <- Reduce(intersect, all_species) %>% Rboretum::vectorSemi() %>% rep(alignment_count)
    }
  } else if(Rboretum::isPhylo(species_info)){ # Get species from phylo
    species_info <- rep(Rboretum::semiSorter(species_info$tip.label),alignment_count)
  } else if(Rboretum::isMultiPhylo(species_info,check_three_taxa=TRUE)){ # Get shared species from multiPhylo
    species_info <- rep(Rboretum::semiSorter(Rboretum::getSharedTaxa(species_info)),alignment_count)
  } else if(is.character(species_info)){ # Get species from character vectors
    if(!Rboretum::semiChecker(species_info)){
      if(length(species_info)>=3){
        species_info <- rep(semiSorter(species_info),alignment_count)
      } else{ 
        stop("'species_info' contains fewer than 3 taxa...")
        }
    } else{
      if(length(species_info)==1){
        species_info <- semiSorter(species_info)
      } else if(length(species_info)>1){
        species_info <- purrr::map(.x=species_info,.f=function(x){semiVector(x)}) %>% unlist() %>% as.character() %>% unique() %>% naturalsort() %>% semiSorter()
      }
      if(length(semiVector(species_info))<3){
        stop("'species_info' contains fewer than 3 taxa...")
      } else{
        species_info <- rep(species_info,alignment_count)
      }
    }
  } else{
    stop("'species_info' should be a phylo, multiPhylo where trees share 3+ taxa, or a character vector with 3+ taxa")
  }
  
  # Get gap data
  if(missing(use_gaps)){
    gap_list <- rep("1",alignment_count) # Default, include gap data
  } else if(is.logical(use_gaps)){
    if(length(use_gaps)==1){
      if(use_gaps){
        gap_list <- rep("1",alignment_count)
      }else{
        gap_list <- rep("0",alignment_count)
      }
    } else{
      if(length(use_gaps)!=alignment_count){
        print(paste(c("'use_gaps' (",length(use_gaps),") and number of alignments (",alignment_count,") do not match...using default setting of including gaps..."),collapse = ''))
        gap_list <- rep("1",alignment_count)
      } else{
        gap_list <- purrr::map(.x=use_gaps,.f=function(x){ifelse(x,"1","0")}) %>% unlist()
      }
    } 
  } else{
    stop("'use_gaps' must be logical (single value or vector)")
  }
  
  if(alignment_count == 1){
    pattern_df <- rb_run_align_patterns(user_python,rb_lib_dir,alignment_path,species_info,gap_list,alignment_name) %>%
      mutate_if(is.list,as.character) %>%
      select(-index)
  } else{
    pattern_df = purrr::map(.x=1:alignment_count,.f=function(x){rb_run_align_patterns(user_python,rb_lib_dir,alignment_path[x],species_info[x],gap_list[x],alignment_name[x])}) %>% do.call(rbind, .) %>%
      mutate_if(is.list,as.character) %>%
      select(-index)
  }
  
  # Python NaN
  pattern_df <- pattern_df %>%
    mutate_all(~na_if(as.character(.), 'NaN'))
  
  pattern_df[is.na(pattern_df)] <- NA
  
  if(nrow(pattern_df)==0){
    stop("No data returned. Check alignments and taxa?")
  } else{
  return(pattern_df)
  }
}
