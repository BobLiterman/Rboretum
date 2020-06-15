context("treeTrimmer")

#read tree
tree2 <- readRooted("inst/extdata/Gene_1.nwk",root_taxa = c('Species_C','Species_H'))

test_that("trimming returns correct subset", {
  
  trimmedTree <- treeTrimmer(tree2,taxa = paste0("Species_", LETTERS[2:14])) #base
  expect_equal(sort(trimmedTree$tip.label), paste0("Species_", LETTERS[2:14]))
  
  trimmedTree <- treeTrimmer(tree2,taxa = c('Species_A','Species_O'), remove = TRUE) #using remove
  expect_equal(sort(trimmedTree$tip.label), paste0("Species_", LETTERS[2:14]))
  
  trimmedTree <- treeTrimmer(tree2,taxa = 'Species_C', remove = TRUE) #remove a root sp
  expect_equal(sort(trimmedTree$tip.label), paste0("Species_", c(LETTERS[1:2],LETTERS[4:15])))
                             
})

#read multiphylo
file_names <- c('inst/extdata/Gene_1.nwk','inst/extdata/Gene_2.nwk','inst/extdata/Gene_3.nwk')
allTrees <- readRooted(to_root = file_names,root_taxa = c('Species_C','Species_H'))

test_that("trimming return correct subset - multiphylo", {
  
  trimmedTrees <- treeTrimmer(allTrees,taxa = paste0("Species_", LETTERS[2:14])) #base
  map(trimmedTrees, "tip.label") %>% map(sort) %>% 
    map(~expect_equal(.x,paste0("Species_", LETTERS[2:14]))) #each tree has correct species
  
  trimmedTrees <- treeTrimmer(allTrees,taxa = c('Species_A','Species_O'), remove = TRUE) #using remove
  map(trimmedTrees, "tip.label") %>% map(sort) %>% 
    map(~expect_equal(.x,paste0("Species_", LETTERS[2:14]))) #each tree has correct species
  
  trimmedTrees <- treeTrimmer(allTrees,taxa = 'Species_C', remove = TRUE) #remove a root sp
  map(trimmedTrees, "tip.label") %>% map(sort) %>% 
    map(~expect_equal(.x,paste0("Species_", c(LETTERS[1:2],LETTERS[4:15])))) #each tree has correct species
  
})