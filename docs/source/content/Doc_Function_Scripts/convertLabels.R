library(Rboretum)
sourceRboretum()

# Read in table of name equivalencies
name_information <- read_tsv(rb_name_file)
head(name_information)

# Read in tree with 'Alignment_IDs' as tip labels
alignmentTree <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
alignmentTree$tip.label

# Convert 'Alignment_IDs' [First column in name table] to 'Number_IDs' [Second column in name table]
numberTree <- convertLabels(to_convert = alignmentTree,name_df = name_information)
numberTree$tip.label

# Convert 'Number_IDs' to 'Common_Names'
commonTree <- convertLabels(to_convert = numberTree,name_df = name_information,from='Number_IDs',to='Common_Names')
commonTree$tip.label

# Convert a set of clades
alignment_clades <- getTreeClades(alignmentTree)
alignment_clades

number_clades <- convertLabels(to_convert = alignment_clades,name_df = name_information,from='Alignment_IDs',to='Number_IDs')
number_clades

# Convert columns in a dataframe
alignment_splits <- getTreeSplits(alignmentTree)
alignment_splits

number_splits <- alignment_splits %>%
  rowwise() %>%
  mutate(Clade=convertLabels(Clade,name_information),
         Mirror_Clade=convertLabels(Mirror_Clade,name_information)) %>% ungroup()

number_splits

# Convert a list of IDs
root_taxa <- alignment_splits %>% filter(Root) %>% pull(Clade) %>% semiVector()

nonroot_taxa <- alignment_splits %>% filter(!Root) %>% pull(Clade) %>% semiVector() %>% unlist() %>% unique()

taxa_list <- list('Root'=root_taxa,'Non_Root'=nonroot_taxa)
taxa_list

number_list <- convertLabels(to_convert = taxa_list ,name_df = name_information,from = 'Alignment_IDs',to='Number_IDs')
number_list