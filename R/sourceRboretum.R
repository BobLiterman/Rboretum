#' Rboretum Python Script Sourcer
#'
#' This function sources all Rborteum Python scripts for 'reticulate'
#' @export
#' 
sourceRboretum <- function(){
  source_python(system.file("", "Alignment_Species.py", package = "Rboretum"))
  source_python(system.file("", "Splits_Processor.py", package = "Rboretum"))
  source_python(system.file("", "Alignment_Patterns.py", package = "Rboretum"))
  source_python(system.file("", "Alignment_Composition.py", package = "Rboretum"))
  source_python(system.file("", "Species_Composition.py", package = "Rboretum"))
}