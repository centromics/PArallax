\name{coh_cor}
\alias{coh_cor}
\title{Exploratory correlation inside the connecting pathways}
\description{Spearman correlation of the query change and the compound change on genes \code{cg}, against the same number of random genes.}
\usage{
coh_cor(cg, z, dchg, n_draw, n_min, seed, tail = c("lower", "upper"))
}
\arguments{
\item{cg}{genes of the connecting pathways.}
\item{z}{compound profile, named by gene.}
\item{dchg}{query change (case mean minus control mean), named by gene.}
\item{n_draw}{random draws.}
\item{n_min}{minimum number of genes.}
\item{seed}{random seed.}
\item{tail}{\code{"lower"}: P of a correlation this low or lower (reversal); \code{"upper"}: this high or higher (similarity).}
}
\value{Named numeric: \code{n}, \code{rho}, \code{bg} (mean of the random correlations), \code{P}.}
