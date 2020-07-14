context("extractNodeAges")

# Read in single chronograms

chronogram_1 <- readRooted('inst/extdata/timeTrees/Chronogram_1.nwk',root_taxa = c('Species_C','Species_H'))
chronogram_2 <- readRooted('inst/extdata/timeTrees/Chronogram_2.nwk',root_taxa = c('Species_C','Species_H'))
chronogram_3 <- readRooted('inst/extdata/timeTrees/Chronogram_3.nwk',root_taxa = c('Species_C','Species_H'))

dates_1 <- extractNodeAges(chronogram_1)
dates_2 <- extractNodeAges(chronogram_2)
dates_3 <- extractNodeAges(chronogram_3)

#run

test_that('data types correct', {
  expect_is(testing_data,'data.frame')
  expect_is(testing_data$numbers, 'double') #fix
})



#repeat for multiphylo

test_that('node summary summarizes', {
  expect_is(testing_data,'data.frame')
  expect_is(testing_data$numbers, 'double') #fix
})