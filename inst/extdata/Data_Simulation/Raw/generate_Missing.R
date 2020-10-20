# Simulate data for 10 'genes' using Rborteum trees
missing_val <- sample(50:500,15,replace=F)
missing_df <- tibble(Species = character(),Percent_N = numeric())

for(x in 1:15){
  missing_df <- missing_df %>% add_row(Species = paste(c('Species_',LETTERS[x]),collapse = ''),Percent_N = missing_val[x]/10000)  
}

write_tsv(missing_df,'inst/extdata/Data_Simulation/Raw/Missing_Table.tsv')