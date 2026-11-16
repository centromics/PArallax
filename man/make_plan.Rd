\name{make_plan}
\alias{make_plan}
\alias{eval_pair}
\title{Same-compound benchmark across cell lines}
\description{\code{make_plan}: for each cell-line pair (A, B), the queries are compounds measured in both (up to \code{qmax}, fixed seed)
and the library is the queries plus up to \code{extra} other compounds of B. \code{eval_pair}: G, P, tiers, confirmation counts,
explanation share and A/B overlap for every query of one pair, in every pathway space.}
\usage{
make_plan(IDX, pairs, full_pairs, qmax, extra, seed)

eval_pair(pl, gm, pw, spaces, prm)
}
\arguments{
\item{IDX}{list per cell line of named vectors; names = compound IDs in their original order.}
\item{pairs, full_pairs}{two-column matrices of cell-line pairs; the seed index follows \code{full_pairs}.}
\item{qmax, extra, seed}{plan size and seed.}
\item{pl}{one element of a plan.}
\item{gm}{function(cell, ids) returning genes x compounds.}
\item{pw}{function(space, cell, ids, weight_cell) returning pathways x compounds.}
\item{spaces}{list of pathway spaces, each with \code{name}, \code{ks}, \code{sets}, \code{tiers}.}
\item{prm}{list of parameters as in \code{\link{rank_list}}, plus \code{ex_top} and \code{small_max}.}
}
\value{\code{make_plan}: a list of pairs. \code{eval_pair}: a data.frame, one row per query and pathway space.}
\examples{
IDX <- list(A = stats::setNames(1:5, paste0("c", 1:5)), B = stats::setNames(1:6, paste0("c", 2:7)))
pr <- rbind(c("A", "B"))
make_plan(IDX, pr, pr, qmax = 3, extra = 2, seed = 1)
}
