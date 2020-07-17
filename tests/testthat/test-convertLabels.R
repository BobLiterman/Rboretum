context("convertLabels")

# Read in test tree and name conversion table
test_tree <- Rboretum::readRooted("inst/extdata/unrootedTrees/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))
name_table <- read_tsv('inst/extdata/Name_Conversion_Table.tsv')

# Convert labels
renamed_tree <- Rboretum::convertLabels(to_convert = test_tree,name_df = name_table,from = 'Alignment_IDs',to='Common_Names')

test_that("Renamed tree is a phylo...", {
  expect_type(renamed_tree, "phylo")
})

test_that("Phylo names substituted correctly...", {
  expect_true(checkTips(test_tree,'Species_A'))
  expect_false(checkTips(renamed_tree,'Species_A'))
  expect_true(checkTips(renamed_tree,'Fuzzy Wollop'))
  expect_false(checkTips(test_tree,'Fuzzy Wollop'))
})

# Read in multiPhylo
file_names <- c('inst/extdata/unrootedTrees/Gene_1.nwk','inst/extdata/unrootedTrees/Gene_2.nwk','inst/extdata/unrootedTrees/Gene_3.nwk')
three_trees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))
renamed_trees <- Rboretum::convertLabels(to_convert = three_trees,name_df = name_table,from = 'Alignment_IDs',to='Common_Names')

test_that("Renamed tree is a multiPhylo...", {
  expect_type(renamed_tree, "multiPhylo")
})

test_that("multiPhylo names substituted correctly...", {
  expect_true(checkTips(three_trees,'Species_A'))
  expect_false(checkTips(renamed_trees,'Species_A'))
  expect_true(checkTips(renamed_trees,'Fuzzy Wollop'))
  expect_false(checkTips(three_trees,'Fuzzy Wollop'))  
})