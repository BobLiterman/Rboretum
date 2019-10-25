#' Get Alignment Support for Specific Topology
#'
#' Given a (1) rooted tree, and
#'         (2) the signal output from getAlignmentSignal from an alignment and matched set of taxa, this function calculates alignment support for each split
#' @param tree Rooted phylo object
#' @param signal Output table from getAlignmentSignal run with same taxa as in tree
#' @param existing_splits OPTIONAL: Should be output of getTreeSupport from the same tree, run with other alignments
#' @return The same split table from getAllSplits(tree), but with a support column listing total site counts supporting each split
#' @export
#' @examples
#' getTreeSupport(tree,signal)
#'

getTreeSupport <- function(tree,signal,existing_splits){

  if(missing(existing_splits)){
    existing_splits <- FALSE
  } else if(is.data.frame(existing_splits)){
    old_splits <- existing_splits
    existing_splits <- TRUE
  } else{
    print("existing_splits argument not a dataframe, processing results separately.")
    existing_splits <- FALSE
  }

  signal_taxa <- signal %>%
    filter(!is.na(Split_1)) %>%
    head(1) %>%
    select(starts_with('Split_')) %>%
    select_if(~ !any(is.na(.))) %>%
    unite(col = "Taxa",sep = ";") %>%
    pull() %>% as.character() %>% str_split(pattern = ";") %>% unlist() %>% sort()

  alignment_name <-  as.character(signal$Alignment_Name[1])

  tree_taxa <- sort(tree$tip.label)

  if(!(all(tree_taxa %in% signal_taxa) & all(signal_taxa %in% tree_taxa))){
    print("Tree Taxa:")
    print(tree_taxa)
    print("Signal Taxa:")
    print(signal_taxa)
    stop("Taxa from signal analysis doesn't match that from tree.")
    }

  splits <- Rboretum::getAllSplits(rooted_tree = tree) %>% mutate(Clade = as.character(Clade),Mirror_Clade = as.character(Mirror_Clade))

  biallelic <- c('biallelic','gap_biallelic')
  biallelic_splits <- signal %>% filter(Site_Pattern %in% biallelic) %>% select(starts_with('Split_')) %>% unlist() %>% table()

  triplus <- c('triallelic','gap_triallelic','quadallelic','gap_quadallelic','pentallelic','gap_pentallelic')
  triplus_splits <- signal %>% filter(Site_Pattern %in% triplus) %>% select(starts_with('Split_')) %>% unlist() %>% table()

  root_split <- splits %>% filter(is.na(Split_Node))
  non_root_splits <- splits %>% filter(!is.na(Split_Node))

  root_clade <- root_split$Clade %>% as.character()
  root_clades <- sort(c(root_clade,root_split$Mirror_Clade %>% as.character()))

  non_root_clades <- non_root_splits %>% pull(Clade) %>% as.character()

  root_support <- tableCount(biallelic_splits,root_clades[1])+tableCount(triplus_splits,root_clades[1])+tableCount(triplus_splits,root_clades[2])

  non_root_support <- c()

  for(clade in non_root_clades){
    non_root_support <- c(non_root_support,(sum(tableCount(biallelic_splits,clade),tableCount(triplus_splits,clade))))
  }

  support_df <- data.frame(Clade = c(root_clade,non_root_clades),Support = c(root_support,non_root_support)) %>%
    mutate(Clade = as.character(Clade),Support=as.integer(Support)) %>%
    left_join(splits) %>%
    select(Clade,Mirror_Clade,Split_Node,Support) %>%
    rename(!!alignment_name := Support)

  if(!existing_splits){
    return(support_df)
  } else{
    if(alignment_name %in% names(old_splits)){
      print("Data from an alignment with that same name already exists in the supplied dataframe. Returning results for tree and alignment provided.")
      return(support_df)
    } else{
      if('Clade' %in% names(old_splits)){
        new_clades <- pull(support_df,Clade) %>% as.character() %>% sort()
        old_clades <- pull(old_splits,Clade) %>% as.character() %>% sort()

        if(length(new_clades) == length(old_clades)){

          if(all(new_clades == old_clades)){
            support_df <- left_join(old_splits,support_df)
            print(paste(c("Added tree split data from",alignment_name),collapse = " "))
            return(support_df)
          } else{
            print("Clades don't match between pre-existing splits and tree provided. Returning results for tree and alignment provided.")
            return(support_df)
          }
        } else{
          print("Clades don't match between pre-existing splits and tree provided. Returning results for tree and alignment provided.")
          return(support_df)
        }
      } else{
          print("Unknown dataframe provided as exisitng splits. Returning results for tree and alignment provided.")
          return(support_df)
        }
      }
  }
}
