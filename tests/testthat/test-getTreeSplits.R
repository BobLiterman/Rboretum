context("treeTrimmer")

#read tree

#generate data frame

#test ideas from https://katherinemwood.github.io/post/testthat/

test_that('data types correct', {
  expect_is(testing_data,'data.frame')
  expect_is(testing_data$numbers, 'integer')
})

test_that('data dimensions correct', {
  expect_equal(ncol(testing_data), 2) #fix numbers
  expect_equal(nrow(testing_data), 4) #fix numbers
})

test_that('no missing values', {
  expect_identical(testing_data, na.omit(testing_data))
})

#repeat for multiphylo