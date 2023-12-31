% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scripts.R
\name{network_function}
\alias{network_function}
\title{Create networks from gene-sets}
\usage{
network_function(
  gsea_result,
  direction = "none",
  edges = "intersection",
  layout,
  cutoff_intersection = 0,
  cutoff_jaccard = 0.7,
  low = "#1465AC",
  mid = "white",
  high = "#B31B21"
)
}
\arguments{
\item{gsea_result}{A Data Frame of a gene-set enrichment analysis result. Required columns: 
Description, NES, core_enrichment.}

\item{direction}{"none" (default) if we want gene-sets with both positive and negative normalized enrichment scores (NESs), "up" if we want gene-sets with only positive NESs, and "down" if we want gene-sets with only negative NESs}

\item{edges}{defaults to "intersection" otherwise Jaccard similarity}

\item{layout}{Layout types available in the ggraph package or created with create_layout()}

\item{low}{Color of the gene-sets with negative NES. Default color: "#1465AC"}

\item{mid}{Color of the gene-sets with zero NES. Default color: "white"}

\item{high}{Color of the gene-sets with positive NES. Default color: "#B31B21"}
}
\description{
Expects a gene-set enrichment analysis result and creates a network graph from it.
The edges of the graph can be either based on the intersection of genes in each set
or Jaccard similarity.
Holds together all the functions from the package.
}
