context("readRooted")

#read tree

test_that("test multiphylo", {
  
  expect_s3_class(tree, "multiPhylo")
  
})

test_that("test single root", {
  
  
  
})

test_that("multi root", {
  
  
  
})

test_that("warning when there is no data", {
  expect_warning(, "  ")
})