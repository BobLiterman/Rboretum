library(Rboretum)
sourceRboretum()

rb_name_file

name_information <- read_tsv(rb_name_file)

myTree <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
