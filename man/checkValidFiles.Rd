% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkValidFiles.R
\name{checkValidFiles}
\alias{checkValidFiles}
\title{Rboretum File Checker}
\usage{
checkValidFiles(alignment_path, return_full_path, return_invalid)
}
\arguments{
\item{alignment_path}{A character vector of one or more alignment file paths (relative or absolute)}

\item{return_full_path}{OPTIONAL: If TRUE and all paths are valid, function returns absolute paths [Default: FALSE]}

\item{return_invalid}{OPTIONAL: If TRUE and any paths are invalid, function returns invalid paths [Default: FALSE]}
}
\value{
TRUE if all paths point to valid files, otherwise FALSE. If return_full_path = TRUE and all paths are valid, function returns absolute paths. If return_invalid = TRUE and there are invalid paths, function returns invalid paths.
}
\description{
This function takes one or more file paths, and returns TRUE if all paths point to valid files, and FALSE otherwise. If return_full_path = TRUE and all paths are valid, function returns absolute paths. If return_invalid = TRUE and there are invalid paths, function returns invalid paths.
}
