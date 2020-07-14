context("treeTrimmer")

# Read in phylo
test_tree <- readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))
test_tree_noG <- treeTrimmer(test_tree,'Species_G',remove=TRUE) # Remove 'Species_G'
test_tree_clade <- treeTrimmer(test_tree,c("Species_A","Species_E","Species_F","Species_G","Species_I","Species_J","Species_K","Species_L","Species_M","Species_N")) # Reduce to this group
clade_lost <- c("Species_B","Species_C","Species_D","Species_H","Species_O")

test_tree_tips <- test_tree$tip.label
test_tree_noG_tips <- test_tree_noG$tip.label
test_tree_clade_tips <- test_tree_clade$tip.label

# Check that treeTrimmer returns a phylo
test_that("treeTrimmer returns phylo...", {
  expect_s3_class(test_tree_noG, "phylo")
  expect_s3_class(test_tree_clade, "phylo")
})

# Check tips
test_that("treeTrimmer removes the correct tips from phylo...", {
  expect_true("Species_G" %in% test_tree_tips)
  expect_false("Species_G" %in% test_tree_noG_tips)
  expect_true("Species_G" %in% test_tree_clade_tips)
  expect_true(all(clade_lost %in% test_tree_tips) & all(clade_lost %in% test_tree_noG_tips))
  expect_false(any(clade_lost %in% test_tree_clade_tips))
})

# Read in multiPhylo
file_names <- c('inst/extdata/unrootedTrees/Gene_1.nwk','inst/extdata/unrootedTrees/Gene_2.nwk','inst/extdata/unrootedTrees/Gene_3.nwk')
test_trees <- readRooted(file_names,root_taxa = c('Species_C','Species_H'))
test_trees_noG <- treeTrimmer(test_trees,'Species_G',remove=TRUE)
test_trees_clade <- treeTrimmer(test_trees,c("Species_A","Species_E","Species_F","Species_G","Species_I","Species_J","Species_K","Species_L","Species_M","Species_N"))

# Create a multiPhylo where trees don't share all taxa
mixed_trees <- test_trees
mixed_trees[[1]] <- test_trees_noG
mixed_trimmed <- treeTrimmer(mixed_trees) # With no taxa list, treeTrimmer reduces multiPhylo to common taxa

# Get all tips from all multiPhylo
test_trees_tips <- purrr::map(.x=test_trees,.f=function(x){x$tip.label})
test_trees_noG_tips <- purrr::map(.x=test_trees_noG,.f=function(x){x$tip.label})
test_trees_clade_tips <- purrr::map(.x=test_trees_clade,.f=function(x){x$tip.label})
mixed_trimmed_tips <- purrr::map(.x=mixed_trimmed,.f=function(x){x$tip.label})

# Check that treeTrimmer returns a multiPhylo
test_that("treeTrimmer returns multiPhylo...", {
  expect_s3_class(test_trees_noG, "multiPhylo")
  expect_s3_class(test_trees_clade, "multiPhylo")
  expect_s3_class(mixed_trimmed, "multiPhylo")
})

# Check tips
test_that("treeTrimmer removes the correct tips from multiPhylo...", {
  expect_true(purrr::map(.x=test_trees_tips,.f=function(x){"Species_G" %in% x}) %>% unlist() %>% all()) # Species_G should be in all trees
  expect_false(purrr::map(.x=test_trees_noG_tips,.f=function(x){"Species_G" %in% x}) %>% unlist() %>% any()) # Species_G should be in no trees
  expect_false(purrr::map(.x=mixed_trimmed_tips,.f=function(x){"Species_G" %in% x}) %>% unlist() %>% any()) # Species_G should be in no trees
  
  expect_true(purrr::map(.x=test_trees_tips,.f=function(x){all(clade_lost) %in% x}) %>% unlist() %>% all())
  expect_true(purrr::map(.x=test_trees_noG_tips,.f=function(x){all(clade_lost) %in% x}) %>% unlist() %>% all())
  expect_true(purrr::map(.x=mixed_trimmed_tips,.f=function(x){all(clade_lost) %in% x}) %>% unlist() %>% all())
  expect_false(purrr::map(.x=test_trees_clade_tips,.f=function(x){any(clade_lost) %in% x}) %>% unlist() %>% any())
})