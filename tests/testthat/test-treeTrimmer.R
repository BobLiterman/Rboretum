context("treeTrimmer")

#read tree

test_that("trimming return correct subset", {
  
  #subset for something
  expect_equal(tree$tip.label, paste0("t", 4:6)) #include expected tips
  
})

test_that("trimming return correct subset - multiphylo", {
  
  #subset for something
  expect_equal(tree[1]$tip.label, paste0("t", 4:6)) #include expected tips
  
})

test_that("trimming return correct subset remove=TRUE", {
  
  #subset for something
  expect_equal(tree$tip.label, paste0("t", 4:6)) #include expected tips
  
})

test_that("trimming return correct subset - multiphylo remove=TRUE", {
  
  #subset for something
  expect_equal(tree[1]$tip.label, paste0("t", 4:6)) #include expected tips
  
})