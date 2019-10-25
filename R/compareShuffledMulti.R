#' Check Split Occurences Among Trees Containing Identical Species
#'
#' This function takes a multiphylo object where each tree contains identical taxa, and returns a summary of split occurrences
#' @param trees Multiphylo object
#' @param tree_names Vector of tree names
#' @param return_shared_only OPTIONAL: If TRUE, returns only splits shared in all trees; Anything else, including missing, returns all splits and their counts (Default: FALSE)
#' @return Dataframe with information about splits among trees
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#'
#' # Return information on all possible splits
#' compareShuffledMulti(trees)
#'
#' # Return information on only those splits shared in all trees
#' compareShuffledMulti(trees,return_shared_only=TRUE)
#'

compareShuffledMulti <- function(trees,tree_names,return_shared_only){

  if(missing(tree_names)){
    tree_names <- unlist(lapply(X = 1:length(trees),function(x) paste(c("tree",x),collapse = "_")))
  }

  if(missing(return_shared_only)){
    return_shared_only <- FALSE
  }

  if(length(trees)!=length(tree_names)){
    print("Multiphylo object and name vector are not the same length. Each tree must be named. Using numbers based on multiphlyo order")
    tree_names <- unlist(lapply(X = 1:length(trees),function(x) paste(c("tree",x),collapse = "_")))
  }

  # Default: Return all splits
  if(return_shared_only != TRUE){
    return_shared_only <- FALSE
  }

  # Ensure multiphylo has >= 2 trees
  tree_count <- length(trees)

  if(!tree_count>=2){
    stop("At least two trees are required for comparison.")
  }

  # Check that all trees contain all species
  species_check <- c()
  for(i in 1:tree_count){
    temp_tree <- trees[[i]]
    species_check <- c(species_check,temp_tree$tip.label)
  }

  species_check_df <- as.data.frame(table(species_check))

  if(any(species_check_df$Freq < tree_count)){
    stop("Trees do not contain the same species.")
  }

  # Tally splits by tree and summarize
  species_list <- sort(trees[[1]]$tip.label)
  all_splits <- c()
  for(i in 1:tree_count){
    temp_splits <- Rboretum::getAllSplits(trees[[i]])
    temp_nonroot <- temp_splits %>% filter(!is.na(Split_Node)) %>% pull(Clade) %>% unlist() %>% as.character()
    temp_root_1 <- temp_splits %>% filter(is.na(Split_Node)) %>% pull(Clade) %>% unlist() %>% as.character()
    temp_root_2 <- temp_splits %>% filter(is.na(Split_Node)) %>% pull(Mirror_Clade) %>% unlist() %>% as.character()
    all_splits <- c(all_splits,temp_nonroot,temp_root_1,temp_root_2)
  }

  tallied_splits <- as.data.frame(table(all_splits))

  split_df <- tallied_splits %>%
    rename(Clade = 'all_splits',Tree_Count = 'Freq') %>%
    mutate(Clade_Size = (str_count(Clade,';')+1)) %>%
    mutate(Clade = as.character(Clade),Tree_Count = as.integer(Tree_Count),Clade_Size = as.integer(Clade_Size))

  # Ensure there are splits unique to some trees (i.e. All trees can't be equal)
  unshared_count <- split_df %>% filter(Tree_Count<tree_count) %>% nrow()

   if(unshared_count==0){
    stop("All trees have the same topology. All splits shared.")
  }

  split_df <- split_df %>% mutate(Tree_Percent = as.numeric(Tree_Count)/as.numeric(tree_count))

  # If return_shared_only = TRUE, just return clades shared by all trees. Otherwise, return all clades
  if(return_shared_only){
    split_df <- split_df %>% filter(Tree_Percent == 1)
    return(split_df)
  } else{
    tree_column <- c()

    for(clade in split_df$Clade){
      tree_list <- c()
      clade_vector <- Rboretum::semiVector(as.character(clade))

      for(i in 1:length(trees)){
        tree_check <- trees[[i]]
        tree_name <- tree_names[i]

        if(checkMonophyletic(tree_check,clade_vector)){
          tree_list <- c(tree_list,tree_name)
        }
      }
      tree_column <- c(tree_column,paste(tree_list,collapse = ";"))
    }
    split_df$Trees_with_Split <- tree_column
    return(split_df)
  }
}
