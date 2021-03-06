% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineBootstraps.R
\name{combineBootstraps}
\alias{combineBootstraps}
\title{Rboretum Bootstrap Combiner}
\usage{
combineBootstraps(trees)
}
\arguments{
\item{trees}{A multiPhylo object where all trees share all taxa and a common topology}
}
\value{
A phylo object with concatenated bootstrap values/node labels
}
\description{
This function takes a multiPhylo of trees that all share a common topolgy, and returns a single Phylo object with combined node labels (e.g. bootstrap values)
}
