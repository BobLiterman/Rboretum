#' Rboretum Python Script Sourcer
#'
#' This function sources all Rborteum Python scripts for 'reticulate', and globally assigns the directory containing the test data to 'test_data_dir'
#' @export
#' 
sourceRboretum <- function(){
  
  # Source python scripts
  source_python(system.file("", "Alignment_Species.py", package = "Rboretum"),envir=globalenv())
  source_python(system.file("", "Splits_Processor.py", package = "Rboretum"),envir=globalenv())
  source_python(system.file("", "Alignment_Patterns.py", package = "Rboretum"),envir=globalenv())
  source_python(system.file("", "Alignment_Composition.py", package = "Rboretum"),envir=globalenv())
  source_python(system.file("", "Species_Composition.py", package = "Rboretum"),envir=globalenv())
  
  test_data_dir <<- system.file("extdata",package = 'Rboretum')
}