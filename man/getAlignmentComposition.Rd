% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAlignmentComposition.R
\name{getAlignmentComposition}
\alias{getAlignmentComposition}
\title{Rboretum Alignment Composition Fetcher}
\usage{
getAlignmentComposition(
  alignment_path,
  species_info,
  alignment_name,
  prefix,
  suffix
)
}
\arguments{
\item{alignment_path}{Where to find alignment files. Options include:
\itemize{
  \item A character vector of one or more alignment file paths  (relative or absolute)
  \item A path to a single directory containing all alignment files (relative or absolute)
}}

\item{species_info}{OPTIONAL: List of taxa to analyze [Default: Process entire alignment]. Can be provided as:
\itemize{
  \item phylo object from which species will be extracted; or
  \item multiPhylo object from which common species will be extracted; or
  \item Character vector with 3+ taxon IDs
  \item Semicolon-separated list of taxon IDs
}}

\item{alignment_name}{OPTIONAL: Chacter vector of names for each alignment. If missing or incomplete, the base filename of the alignment is used}

\item{prefix}{OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file prefixes (e.g. all alignment files start with "Mafft_")}

\item{suffix}{OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file suffixes (e.g. all alignment files end with ".phy")}
}
\value{
Dataframe containing alignment length, %GC, %N, and %Gap for each alignment
}
\description{
Given the path(s) to one or more alignments (and an optional list of taxa), this script returns alignment lengths, %GC, %N, and %Gap for the alignment
}
