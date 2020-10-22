library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
mySpecies <- getAlignmentSpecies(myAlignmentFile)
mySubspecies <- semiVector(mySpecies)[1:5]
myAlignmentDir <- rb_alignment_dir

# Get alignment signal information for a single alignment
getAlignmentSignal(alignment_path = myAlignmentFile)

# Get alignment signal information from all .phylip files in a directory, providing new names, consider gaps as missing data
getAlignmentSignal(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment signal information from all .phylip files in a directory, providing new names, consider gaps as missing data
getAlignmentSignal(alignment_path = myAlignmentDir,species_info = mySubspecies,use_gaps = FALSE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment signal from dummy alignment, with and without gap support
getAlignmentSignal(alignment_path = rb_dummy_align_path)

getAlignmentSignal(alignment_path = rb_dummy_align_path,use_gaps = FALSE)

