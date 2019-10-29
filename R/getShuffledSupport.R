#' Plot Alignment Support for Specific Topology
#'
#' DESCRIPTION
#' @param signal Output table from getAlignmentSignal run with same taxa as in tree
#' @param existing_support OPTIONAL: Should be output of getShuffledSupport from the same comparison, run with different signal
#' @return Tree with nodepoint geoms scaled to support, and colored by number of trees containing that split
#' @export
#' @examples
#' getShuffledSupport(comparison,signal,existing_support)
#'

getShuffledSupport <- function(comparison,signal,existing_support,max_missing){

  if(missing(existing_support)){
    existing_support <- FALSE
  } else if(is.data.frame(existing_support)){
    old_support <- existing_support
    existing_support <- TRUE
  } else{
    print("existing_support argument not a dataframe, processing results separately.")
    existing_support <- FALSE
  }
  
  if(missing(max_missing)){
    max_missing <- 0
  }
  
  signal <- signal %>% 
    filter(Non_Base_Count <= max_missing)

  # Get taxa from signal analysis
  signal_taxa <- signal %>%
    filter(!grepl("non_base",Site_Pattern)) %>%
    filter(!is.na(Split_1)) %>%
    head(1) %>%
    select(starts_with('Split_')) %>%
    select_if(~ !any(is.na(.))) %>%
    unite(col = "Taxa",sep = ";") %>%
    pull() %>% as.character() %>% str_split(pattern = ";") %>% unlist() %>% sort()

  # Get taxa from multiphylo comparison
  comparison_taxa <- comparison %>%
    pull(Clade) %>%
    paste(collapse = ';') %>%
    semiVector() %>%
    unique() %>%
    sort()

  # Ensure taxa sets match
  if(!all(sapply(list(comparison_taxa), FUN = identical, signal_taxa))){
    print("Signal Taxa:")
    print(signal_taxa)
    print("Comparison Taxa:")
    print(comparison_taxa)
    stop("Taxa from arguments don't match")
  }

  informative_patterns <- c('non_base_biallelic','non_base_gap_biallelic',
                            'non_base_triallelic','non_base_gap_triallelic',
                            'non_base_quadallelic','non_base_gap_quadallelic',
                            'non_base_pentallelic','non_base_gap_pentallelic',
                            'biallelic','gap_biallelic',
                            'triallelic','gap_triallelic',
                            'quadallelic','gap_quadallelic',
                            'pentallelic','gap_pentallelic')

  signal_counts <- signal %>%
    filter(Site_Pattern %in% informative_patterns) %>% 
    select(starts_with('Split_')) %>%
    unlist() %>% table()

  alignment_name <- paste(c(signal$Alignment_Name[1],"_m",as.character(max_missing)),collapse = '')

  support_list <- c()
  for(clade in comparison$Clade){
    support_list <- c(support_list,tableCount(signal_counts,clade))
  }

  comparison$Support <- as.integer(support_list)

  comparison <- comparison %>%
    rename(!!alignment_name := Support)

  if(!existing_support){
    return(comparison)
    } else{
      if('Clade' %in% names(old_support)){
        new_clades <- pull(comparison,Clade) %>% as.character() %>% sort()
        old_clades <- pull(old_support,Clade) %>% as.character() %>% sort()

        if(length(new_clades) == length(old_clades)){

          if(all(new_clades == old_clades)){
            comparison <- left_join(old_support,comparison)
            print(paste(c("Added tree split data from",alignment_name),collapse = " "))
            return(comparison)
          } else{
            print("Clades don't match between pre-existing splits and tree provided. Returning results for tree and alignment provided.")
            return(comparison)
          }
        } else{
          print("Clades don't match between pre-existing splits and tree provided. Returning results for tree and alignment provided.")
          return(comparison)
        }
      } else{
        print("Unknown dataframe provided as exisitng splits. Returning results for tree and alignment provided.")
        return(comparison)
      }
    }
}
