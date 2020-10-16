library(Rboretum)
sourceRboretum()

# Set path to alignment data
myAlignmentFile <- rb_align1_path
mySpecies <- getAlignmentSpecies(myAlignmentFile)
myAlignmentDir <- rb_alignment_dir

# Get alignment signal information for a single alignment
mySignal <- getAlignmentSignal(alignment_path = myAlignmentFile,species_info = mySpecies)

# Get alignment signal information from all .phylip files in a directory, providing new names, consider gaps as missing data
mySignals <- getAlignmentSignal(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

# Get alignment signal from dummy alignment
print(getAlignmentSignal(alignment_path = rb_dummy_align_path,use_gaps = TRUE)[c(8,9,14,19,24,27,29),])

