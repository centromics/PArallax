\name{make_signature}
\alias{make_signature}
\title{Query signature from a two-group matrix}
\description{Top up and down features of case vs control samples.}
\usage{
make_signature(x, group, n, method = c("rankprod", "meandiff"), nperm = 100L, seed = 1L)
}
\arguments{
\item{x}{numeric matrix, features (genes or pathways) x samples.}
\item{group}{two levels, one per column; the first level is the control.}
\item{n}{features per direction.}
\item{method}{\code{"rankprod"}: single-study rank product (the MetaDE 1.0.5 algorithm used in the paper; built in, no MetaDE needed); \code{"meandiff"}: case mean minus control mean.}
\item{nperm, seed}{rank-product permutations and seed.}
}
\value{A list: \code{up}, \code{dn}, \code{table}, \code{method}.}
\examples{
toy <- parallax_toy()
s <- make_signature(toy$expr, toy$group, n = 20, method = "rankprod", nperm = 20)
head(s$up); head(s$dn)
}
