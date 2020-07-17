context("getTreeClades")

# Read in single phylos
tree_1 <- readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))
tree_1_subtree <- subtrees(tree_1)
tree_1_subtree_clades <- purrr::map(.x=tree_1_subtree[-1],.f=function(x){semiSorter(x$tip.label)}) %>% unlist() %>% naturalsort()
clades_1_noroot <- getTreeClades(tree_1)
clades_1_root <- getTreeClades(tree_1,include_root = TRUE)

tree_2 <- readRooted("inst/extdata/unrootedTrees/Gene_4.nwk",root_taxa = c('Species_C','Species_H'))
tree_2_subtree <- subtrees(tree_2)
tree_2_subtree_clades <- purrr::map(.x=tree_2_subtree[-1],.f=function(x){semiSorter(x$tip.label)}) %>% unlist() %>% naturalsort()
clades_2_noroot <- getTreeClades(tree_2)
clades_2_root <- getTreeClades(tree_2,include_root = TRUE)

tree_3 <- readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = 'Species_B')
tree_3_subtree <- subtrees(tree_3)
tree_3_subtree_clades <- purrr::map(.x=tree_3_subtree[-1],.f=function(x){semiSorter(x$tip.label)}) %>% unlist() %>% naturalsort()
clades_3_noroot <- getTreeClades(tree_3)
clades_3_root <- getTreeClades(tree_3,include_root = TRUE)

test_that('getTreeClades returning character vectors for a phylo...', {
  expect_type(clades_1_noroot, 'character')
  expect_type(clades_1_root, 'character')
  expect_type(clades_2_noroot, 'character')
  expect_type(clades_2_root, 'character')
  expect_type(clades_3_noroot, 'character')
  expect_type(clades_3_root, 'character')
})

test_that('getTreeClades returning the correct clades for a phylo...', {
  expect_true(all(clades_1_noroot %in% tree_1_clades))
  expect_true(all(clades_1_root %in% tree_1_clades))
  expect_true(length(clades_1_root) - length(clades_1_noroot)==2)
  expect_true(all(clades_2_noroot %in% tree_2_clades))
  expect_true(all(clades_2_root %in% tree_2_clades))
  expect_true(length(clades_2_root) - length(clades_2_noroot)==2)
  expect_true(all(clades_3_noroot %in% tree_3_clades))
  expect_true(all(clades_3_root %in% tree_3_clades))
  expect_true(length(clades_3_root) - length(clades_3_noroot)==0)
})

# Read in multiPhylo
file_names <- c('inst/extdata/unrootedTrees/Gene_1.nwk','inst/extdata/unrootedTrees/Gene_4.nwk')
trees <- readRooted(file_names,c('Species_C','Species_H'))
tree_names <- c(names(trees),'Single_Root')
trees <- c(trees,tree_3) %>% `names<-`(tree_names)

clades_with_root <- getTreeClades(trees,include_root = TRUE)
clades_without_root <- getTreeClades(trees)
clade_counts <- getTreeClades(trees,return_counts = TRUE)

all_no_root_clades <- c(clades_1_noroot,clades_2_noroot,clades_3_noroot) %>% unique() %>% naturalsort()
all_root_clades <- c(clades_1_root,clades_2_root,clades_3_root) %>% unique() %>% naturalsort()

test_that('getTreeClades returning character vectors and dataframes for a multiPhylo...', {
  expect_type(clades_with_root, 'character')
  expect_type(clades_without_root, 'character')
  expect_s3_class(clade_counts, 'data.frame')
  expect_type(clade_counts, 'list')
  expect_true(all(names(clade_counts) == c('Clade','Count','Trees')))
  expect_type(clade_counts$Clade,'character')
  expect_type(clade_counts$Count,'integer')
  expect_type(clade_counts$Trees,'character')
})

test_that('getTreeClades returning the correct clades for a multiPhylo...', {
  expect_true(all(clades_with_root == all_root_clades))
  expect_true(all(clades_without_root == all_no_root_clades))
})