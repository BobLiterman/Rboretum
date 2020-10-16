library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
mySpecies <- getAlignmentSpecies(myAlignmentFile)
myAlignmentDir <- rb_alignment_dir

# Get alignment Breakdown information for a single alignment
getAlignmentBreakdown(alignment_path = myAlignmentFile)

# Get alignment Breakdown information from all .phylip files in a directory, providing new names, consider gaps as missing data
getAlignmentBreakdown(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment Breakdown from dummy alignment, with and without gap support
getAlignmentBreakdown(alignment_path = rb_dummy_align_path)

getAlignmentBreakdown(alignment_path = rb_dummy_align_path,use_gaps = FALSE)