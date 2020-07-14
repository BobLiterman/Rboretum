context("readRooted")

# Read in phylo rooted with 1 or 2 taxa
one_root <- readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = 'Species_H')
two_root <- readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))

one_mirror <- one_root$tip.label[one_root$tip.label != 'Species_H']
two_mirror <- two_root$tip.label[!two_root$tip.label %in% c('Species_C','Species_H')]

# Test that readRooted returns a phylo
test_that("Trees read as phylo...", {
  expect_s3_class(one_root, "phylo")
  expect_s3_class(two_root, "phylo")
})

# Test that readRooted returns a rooted phylo
test_that("Trees read as rooted phylo...", {
  expect_true(ape::is.rooted.phylo(one_root))
  expect_true(ape::is.rooted.phylo(two_root))
})

# Test that readRooted returns a phylo rooted at the desired taxa
test_that("readRooted placed roots at the correct nodes for phylo...", {
  expect_true(ape::is.monophyletic(one_root,'Species_H'))
  expect_true(ape::is.monophyletic(one_root,one_mirror))
  
  expect_true(ape::is.monophyletic(two_root,c('Species_C','Species_H')))
  expect_true(ape::is.monophyletic(two_root,two_mirror))
})

# Read in multiPhylo
file_names <- c('inst/extdata/unrootedTrees/Gene_1.nwk','inst/extdata/unrootedTrees/Gene_2.nwk','inst/extdata/unrootedTrees/Gene_3.nwk')

one_root_trees <- readRooted(to_root = file_names,root_taxa = 'Species_H')
two_root_trees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))

# Test that readRooted returns a multiPhylo
test_that("Trees read as multiPhylo...", {
  expect_s3_class(one_root_trees, "multiPhylo")
  expect_s3_class(two_root_trees, "multiPhylo")
})

# Test that readRooted returns a rooted multiPhylo
test_that("Trees read as rooted multiPhylo...", {
  expect_true(purrr::map(.x=one_root_trees,.f=function(x){ape::is.rooted.phylo(x)}) %>% unlist() %>% all())
  expect_true(purrr::map(.x=two_root_trees,.f=function(x){ape::is.rooted.phylo(x)}) %>% unlist() %>% all())
})

# Test that readRooted returns a phylo rooted at the desired taxa
test_that("readRooted placed roots at the correct nodes for multiPhylo...", {
  expect_true(purrr::map(.x=one_root_trees,.f=function(x){ape::is.monophyletic(x,'Species_H')}) %>% unlist() %>% all())
  expect_true(purrr::map(.x=one_root_trees,.f=function(x){ape::is.monophyletic(x,one_mirror)}) %>% unlist() %>% all())
  
  expect_true(purrr::map(.x=two_root_trees,.f=function(x){ape::is.monophyletic(x,c('Species_C','Species_H'))}) %>% unlist() %>% all())
  expect_true(purrr::map(.x=two_root_trees,.f=function(x){ape::is.monophyletic(x,two_mirror)}) %>% unlist() %>% all())
})
