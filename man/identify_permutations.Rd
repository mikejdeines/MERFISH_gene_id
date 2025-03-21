% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneid.R
\name{identify_permutations}
\alias{identify_permutations}
\title{Finds all permutations of 2 clusters in a Seurat object.
The cluster identities should be in the Idents slot of the Seurat object.}
\usage{
identify_permutations(seuratobject)
}
\arguments{
\item{seuratobject}{a Seurat object with cluster identities in the Idents slot}
}
\value{
an n-by-2 integer of permutations between all pairs of clusters
}
\description{
Finds all permutations of 2 clusters in a Seurat object.
The cluster identities should be in the Idents slot of the Seurat object.
}
