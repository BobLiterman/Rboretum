library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
mySpecies <- getAlignmentSpecies(myAlignmentFile)
mySubspecies <- semiVector(mySpecies)[1:5]
myAlignmentDir <- rb_alignment_dir

# Get alignment Patterns information for a single alignment
getAlignmentPatterns(alignment_path = myAlignmentFile)

# Get alignment Patterns information from all .phy files in a directory, providing new names, consider gaps as missing data
getAlignmentPatterns(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment Patterns information from all .phy files in a directory, providing new names, consider gaps as missing data, consider only Species A - E
getAlignmentPatterns(alignment_path = myAlignmentDir,species_info = mySubspecies,use_gaps = FALSE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment Patterns from dummy alignment, with and without gap support
getAlignmentPatterns(alignment_path = rb_dummy_align_path)

getAlignmentPatterns(alignment_path = rb_dummy_align_path,use_gaps = FALSE)