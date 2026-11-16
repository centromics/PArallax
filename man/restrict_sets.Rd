\name{restrict_sets}
\alias{restrict_sets}
\alias{convert_sets}
\title{Gene sets restricted to a gene universe}
\description{\code{restrict_sets} keeps the genes in \code{univ} (same identifiers) and the sets of \code{min_size} to
\code{max_size} genes. \code{convert_sets} first converts between gene symbols and Entrez IDs (needs AnnotationDbi and org.Hs.eg.db).}
\usage{
restrict_sets(sets, univ, min_size = 5L, max_size = 500L)

convert_sets(sets, univ, min_size = 5L, max_size = 500L)
}
\arguments{
\item{sets}{named list of gene sets.}
\item{univ}{character vector of gene identifiers.}
\item{min_size, max_size}{set size limits after restriction.}
}
\value{Named list of gene sets.}
\examples{
toy <- parallax_toy()
s <- restrict_sets(toy$pathways$KEGG, rownames(toy$library)[1:300], min_size = 5)
head(lengths(s))
}
