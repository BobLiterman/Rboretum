library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
mySpecies <- getAlignmentSpecies(myAlignmentFile)
myAlignmentDir <- rb_alignment_dir

# Get alignment Patterns information for a single alignment
myPatterns <- getAlignmentPatterns(alignment_path = myAlignmentFile,species_info = mySpecies)

# Get alignment Patterns information from all .phylip files in a directory, providing new names, consider gaps as missing data
myPatternss <- getAlignmentPatterns(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment Patterns from dummy alignment
print(getAlignmentPatterns(alignment_path = rb_dummy_align_path,use_gaps = TRUE)[c(8,9,14,19,24,27,29),])

