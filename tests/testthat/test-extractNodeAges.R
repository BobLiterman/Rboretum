context("extractNodeAges")

# Read in single chronograms

chronogram_1 <- readRooted('inst/extdata/timeTrees/Chronogram_1.nwk',root_taxa = c('Species_C','Species_H'))
chronogram_2 <- readRooted('inst/extdata/timeTrees/Chronogram_2.nwk',root_taxa = c('Species_C','Species_H'))
chronogram_3 <- readRooted('inst/extdata/timeTrees/Chronogram_3.nwk',root_taxa = c('Species_C','Species_H'))

chronogram_1_times <- branching.times(chronogram_1)
chronogram_2_times <- branching.times(chronogram_2)
chronogram_3_times <- branching.times(chronogram_3)

dates_1 <- extractNodeAges(chronogram_1)
dates_2 <- extractNodeAges(chronogram_2)
dates_3 <- extractNodeAges(chronogram_3)

test_that('extractNodeAges is returning expected dataframes for a phylo...', {
  expect_type(dates_1,'data.frame')
  expect_type(dates_2,'data.frame')
  expect_type(dates_3,'data.frame')
  
  expect_true(all(names(dates_1)==c('Clade','Node_Age')))
  expect_true(all(names(dates_2)==c('Clade','Node_Age')))
  expect_true(all(names(dates_3)==c('Clade','Node_Age')))
  
  expect_type(dates_1$Clade, 'character')
  expect_type(dates_2$Clade, 'character')
  expect_type(dates_3$Clade, 'character')
  
  expect_type(dates_1$Node_Age, 'double')
  expect_type(dates_2$Node_Age, 'double')
  expect_type(dates_3$Node_Age, 'double')
})

test_that('extractNodeAges is returning the correct node ages for a phylo...', {
  expect_true(all(chronogram_1_times == dates_1$Node_Age))
  expect_true(all(chronogram_2_times == dates_2$Node_Age))
  expect_true(all(chronogram_3_times == dates_3$Node_Age))
})

# Read in multiple chronograms
chronograms <- readRooted(c('inst/extdata/timeTrees/Chronogram_1.nwk','inst/extdata/timeTrees/Chronogram_2.nwk','inst/extdata/timeTrees/Chronogram_3.nwk'),c('Species_C','Species_H'))
chronogram_names <- names(chronograms)

# Strip node labels for comparison
chronograms <- Rboretum::stripNodeLabels(chronograms)
multi_branching_times <- purrr::map(.x=chronograms,.f=function(x){branching.times(x)})

clade_df <- tibble(Clade=character(),Node_Age=double(),Tree_Name=character())

# Get Node/Clade/Age data manually
for(i in 1:3){
  temp_subtree <- subtrees(chronograms[[i]])
  temp_name <- chronogram_names[[i]]
  node_labels <- purrr::map(.x=temp_subtree,.f=function(x){as.character(x$node.label[1])}) %>% unlist()
  clades <- purrr::map(.x=temp_subtree,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()
  clade_df <- rbind(clade_df,tibble(Clade=clades,Node=node_labels,Tree_Name=temp_name) %>% rowwise() %>% mutate(Node_Age = multi_branching_times[[i]][Node]) %>% select(-Node))
}

clade_df <- select(clade_df,Clade,Node_Age,Tree_Name) %>% arrange(Clade,Node_Age)

mean_df <- clade_df %>% group_by(Clade) %>% summarize(Mean_Node_Age = mean(Node_Age))
manual_mean_values <- setNames(mean_df$Mean_Node_Age, mean_df$Clade)

median_df <- clade_df %>% group_by(Clade) %>% summarize(Median_Node_Age = median(Node_Age))
manual_median_values <- setNames(median_df$Median_Node_Age, median_df$Clade)

# Get date information from function
raw_dates <- extractNodeAges(chronograms) %>% arrange(Clade,Node_Age)

mean_dates <- extractNodeAges(chronograms,return_summary = 'mean')
function_mean_values <- setNames(mean_dates$Node_Age, mean_dates$Clade)

median_dates <- extractNodeAges(chronograms,return_summary = 'median')
function_median_values <- setNames(median_dates$Node_Age, median_dates$Clade)

both_dates <- extractNodeAges(chronograms,return_summary = 'both')
both_mean_values <- setNames(both_dates$Mean_Node_Age, both_dates$Clade)
both_median_values <- setNames(both_dates$Median_Node_Age, both_dates$Clade)

test_that('extractNodeAges is returning expected dataframes for a multiPhylo...', {
  
  expect_type(raw_dates,'data.frame')
  expect_type(mean_dates,'data.frame')
  expect_type(median_dates,'data.frame')
  expect_type(both_dates,'data.frame')
  
  expect_true(all(names(raw_dates)==c('Clade','Node_Age','Tree_Name')))
  expect_true(all(names(mean_dates)==c('Clade','Node_Age')))
  expect_true(all(names(median_dates)==c('Clade','Node_Age')))
  expect_true(all(names(both_dates)==c("Clade","Mean_Node_Age","Median_Node_Age","StdDev_Node_Age","MAD_Node_Age")))
  
  expect_type(raw_dates$Clade, 'character')
  expect_type(mean_dates$Clade, 'character')
  expect_type(median_dates$Clade, 'character')
  expect_type(both_dates$Clade, 'character')
  
  expect_type(raw_dates$Node_Age, 'double')
  expect_type(mean_dates$Node_Age, 'double')
  expect_type(median_dates$Node_Age, 'double')

  expect_type(raw_dates$Tree_Name, 'character')
  
  expect_type(both_dates$Mean_Node_Age, 'double')
  expect_type(both_dates$Median_Node_Age, 'double')
  expect_type(both_dates$StdDev_Node_Age, 'double')
  expect_type(both_dates$MAD_Node_Age, 'double')
})

test_that('extractNodeAges is returning the correct node ages for a multiPhylo...', {
  expect_true(all(clade_df == raw_dates))
  expect_true(all(manual_mean_values == function_mean_values))
  expect_true(all(manual_median_values == function_median_values))
})