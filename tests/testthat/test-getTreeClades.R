context("getTreeClades")

#read tree

test_that('data types correct', {
  expect_is(testing_data,'vector')
  expect_is(testing_data, 'character')
})

#repeat for include_root

test_that('include counts returns table', {
  expect_is(testing_data,'table')
  expect_is(testing_data$, 'character')
  expect_is(testing_data$, 'integer')
})

#repeat for multiphylo