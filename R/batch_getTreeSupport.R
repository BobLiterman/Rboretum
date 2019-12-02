#' Rboretum Batch Alignment Signal Support Mapper
#'
#' This function computes the support for a specified tree or a set of trees from an alignment or set of alignments. To use this function, you must supply either multiple trees or multiple alignment signals. 
#' @param signal Output table from getAlignmentSignal() or batch_getAlignmentSignal()
#' @param tree Tree(s) on which to map alignment support. Options include:
#' \itemize{
#'   \item A single phylo object where taxa EXACTLY match those from signal [Returns results as a single dataframe with multiple support columns]; or, 
#'   \item A named multiPhylo object, where all trees have (1) taxa that EXACTLY match those from signal and (2) a unique topology [Returns results as a named list of dataframes]
#' }
#' @param alignment_name OPTIONAL: Character vector of alignment names for each 'signal' passed [Default: Alignment name from signal dataframe]
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column [Default: 0]
#' @param include_gap OPTIONAL: TRUE or FALSE; Count sites with gap positions ('-') as part of total support [Default: TRUE]
#' @param include_singleton OPTIONAL: TRUE or FALSE; Count sites with singletons as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: TRUE or FALSE; Count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: TRUE or FALSE; Count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: TRUE or FALSE; Count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: TRUE or FALSE; Count sites with pentallelic variation as part of total support [Default: TRUE]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') as part of total support [Default: FALSE]
#' @return The same split table from getTreeSplits(tree), but with a support column for the specfied alignment/missing combination
#' @export

batch_getTreeSupport <- function(tree,signal,max_missing,alignment_name,include_gap,include_singleton,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic,only_gap,existing_splits){

  if(!Rboretum::isMultiPhylo(tree,check_rooted = TRUE) & !Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'tree' does not appear to be a valid + rooted phylo or multiPhylo object")
  }
  
  if(Rboretum::isPhylo(tree)){
    tree_taxa <- sort(tree$tip.label)
    tree_count <- 1
  }
  
  if(Rboretum::isMultiPhylo(tree)){
    if(is.null(names(tree))){
      stop("'tree' multiPhlyo must have names assigned via names(tree) <- c('Name1','Name2',etc.)")
    } else if(!checkSameTaxa(tree)){
        stop("All trees in 'tree' must contain identical taxa.")
      } else{
        compare_vector <- Rboretum::compareTrees(tree) %>% pull(Comparable)
        if(!all(compare_vector)){
          Rboretum::compareTrees(tree) %>% select(Tree_1,Tree_2,Comparable)
          stop("'Some trees aren't comparable, meaning they have different taxa or the same topology.")
        } else{
          tree_taxa <- getSharedTaxa(tree)
          tree_count <- length(tree)
          tree_names <- names(tree)
          if(any(duplicated(tree_names))){
            stop("'tree' multiPhlyo contains trees with identical names.")
          }
        }
      }
    }
  
  if(!Rboretum::isAlignmentSignal(signal,tree_taxa)){
    stop("'signal' does not appear to be the ouput from getAlignmentSignal() or batch_getAlignmentSignal, and/or 'signal' taxa don't match those from 'tree'")
  } else{
      signal_taxa <- signal %>%
        filter(!is.na(Split_1)) %>%
        head(1) %>%
        select(Singleton_Taxa,Non_Base_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5) %>%
        select_if(~ !any(is.na(.))) %>%
        unite(col = "Taxa",sep = ";") %>%
        semiVector() %>%
        sort()
      
      alignment_count <- length(unique(signal$Alignment_Name))
      raw_names <- unique(signal$Alignment_Name)
      
      if(alignment_count == 1 & tree_count == 1){
        stop("Only a single tree and signal provided. Use getTreeSupport().")
      } else{
          one_count <- length(which(signal$Alignment_Position == 1))
          if(alignment_count != one_count){
            print(paste(c('Alignment Names:',alignment_count),collapse = ' '))
            print(paste(c('Position 1 Counts:',one_count),collapse = ' '))
            stop("'signal' contains a mismatched number of names and first positions.")
          }
      }
  }
  
  if(missing(max_missing)){
    max_missing <- as.integer(0)
  } else if(has_error(silent=TRUE,as.integer(max_missing))){
    print("Invalid 'max_missing', allowing 0 missing taxa...")
    max_missing <- as.integer(0)
  } else if(length(signal_taxa) - as.integer(max_missing) < 3){
    print("'max_missing' value too high, fewer than 3 taxa required. Using largest possible value...")
    max_missing <- as.integer(length(signal_taxa) - 3)
  } else{
    max_missing <- as.integer(max_missing)
  }
  
  if(missing(alignment_name)){
    if(alignment_count == 1){
      alignment_name <-  paste(c(as.character(signal$Alignment_Name[1]),'_m',as.character(max_missing)),collapse = '')
    } else{
      alignment_name <- c()
      for(i in 1:alignment_count){
        alignment_name <- c(alignment_name,paste(c(raw_names[i],'_m',as.character(max_missing)),collapse = ''))
      }
    }
  } else if(length(alignment_name)!=alignment_count){
    print("'alignment_name' vector and number of alignments detected do not match. Generating default names...")
    if(alignment_count == 1){
      alignment_name <-  paste(c(as.character(signal$Alignment_Name[1]),'_m',as.character(max_missing)),collapse = '')
    } else{
      alignment_name <- c()
      for(i in 1:alignment_count){
        alignment_name <- c(alignment_name,paste(c(raw_names[i],'_m',as.character(max_missing)),collapse = ''))
      }
    }
  } else{
    alignment_name <- as.character(alignment_name)
    print('Alignment Name Information:')
    for(i in 1:alignment_count){
      print(paste(c(raw_names[i],":",alignment_name[i]),collapse = ' '))
    }
  }
  
  if(any(duplicated(alignment_name))){
    stop("'alignment_name' includes duplicate entries. Fix input data such that each alignment file has a unique ID.")
  }
  
  if(missing(include_gap)){
    include_gap <- TRUE
  } else if (!is.logical(include_gap)){
    include_gap <- TRUE
  }
  
  if(missing(include_singleton)){
    include_singleton <- TRUE
  } else if (!is.logical(include_singleton)){
    include_singleton <- TRUE
  }

  if(missing(include_biallelic)){
    include_biallelic <- TRUE
  } else if (!is.logical(include_biallelic)){
    include_biallelic <- TRUE
  }

  if(missing(include_triallelic)){
    include_triallelic <- TRUE
  } else if (!is.logical(include_triallelic)){
    include_triallelic <- TRUE
  }

  if(missing(include_quadallelic)){
    include_quadallelic <- TRUE
  } else if (!is.logical(include_quadallelic)){
    include_quadallelic <- TRUE
  }

  if(missing(include_pentallelic)){
    include_pentallelic <- TRUE
  } else if (!is.logical(include_pentallelic)){
    include_pentallelic <- TRUE
  }

  if(missing(only_gap)){
    only_gap <- FALSE
  } else if (!is.logical(only_gap)){
    only_gap <- FALSE
  }
 
  informative_patterns <- c('biallelic','triallelic','quadallelic','pentallelic')

  signal <- signal %>%
    filter(Non_Base_Count <= max_missing) %>%
    filter(Site_Pattern %in% informative_patterns)

  if(!include_gap){
    if(only_gap){
      stop("Cannot only use (only_gap) and exclude (include_gap) gap positions.") }
    signal <- signal %>%
      filter(Gap==FALSE)
  }

  if(!include_singleton){
    signal <- signal %>%
      filter(Singleton==FALSE)
  }

  if(!include_biallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'biallelic'))
  }

  if(!include_triallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'triallelic'))
  }

  if(!include_quadallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'quadallelic'))
  }

  if(!include_pentallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'pentallelic'))
  }

  if(only_gap){
    if(!include_gap){
      stop("Cannot only use (only_gap) and exclude (include_gap) gap positions.")
    } else if(signal %>% filter(Gap==TRUE) %>% nrow() == 0){
        stop("Data contains no gap positions, but 'only_gap' was specified.")
      } else{
    signal <- signal %>%
      filter(Gap==TRUE)
      }
  }
  
  if(tree_count == 1){
    support_df <- Rboretum::getTreeSplits(tree) %>%
      filter(!is.na(Split_Node))%>%
      mutate(Clade = as.character(Clade),Mirror_Clade = as.character(Mirror_Clade))

    clades <- support_df %>% pull(Clade) %>% as.character() %>% sort()

    for(i in 1:alignment_count){
      
      temp_name <- alignment_name[i]
      
      all_signal_splits <- signal %>% filter(Alignment_Name == raw_names[i]) %>% select(starts_with('Split_')) %>% unlist() %>% table()
  
      clade_support <- c()
      for(clade in clades){
        clade_support <- c(clade_support,tableCount(all_signal_splits,clade))
      }
  
      temp_df <- data.frame(Clade = clades,Support = clade_support) %>%
        mutate(Clade = as.character(Clade),Support = as.integer(Support)) %>%
        rename(!!temp_name := Support)
      
      support_df <- support_df %>%
        left_join(temp_df,by='Clade')
    }
    
    return(support_df)
  } else{
    
    support_list <- list()
    
    for(i in 1:tree_count){
    
      support_df <- Rboretum::getTreeSplits(tree[[i]]) %>%
        filter(!is.na(Split_Node))%>%
        mutate(Clade = as.character(Clade),Mirror_Clade = as.character(Mirror_Clade))
      
      clades <- support_df %>% pull(Clade) %>% as.character() %>% sort()

      for(j in 1:alignment_count){
        
        temp_name <- alignment_name[j]
        
        all_signal_splits <- signal %>% filter(Alignment_Name == raw_names[j]) %>% select(starts_with('Split_')) %>% unlist() %>% table()
        
        clade_support <- c()
        for(clade in clades){
          clade_support <- c(clade_support,tableCount(all_signal_splits,clade))
        }
        
        temp_df <- data.frame(Clade = clades,Support = clade_support) %>%
          mutate(Clade = as.character(Clade),Support = as.integer(Support)) %>%
          rename(!!temp_name := Support)
        
        support_df <- support_df %>%
          left_join(temp_df,by='Clade')
        
        support_list[[i]] <- support_df
      }
    }
    names(support_list) <- tree_names
    return(support_list)
  }
}