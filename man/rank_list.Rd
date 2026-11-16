\name{rank_list}
\alias{rank_list}
\title{Rank a library for one query: G, P, tiers and connecting pathways}
\description{Scores every compound with WTCS on genes and on pathways and assigns tiers.}
\usage{
rank_list(gq, pq, LG, LP, sets, prm, direction = "similar", split = TRUE)
}
\arguments{
\item{gq, pq}{lists with elements \code{up} and \code{dn}: query genes and query pathways, ordered.}
\item{LG}{library, genes x compounds.}
\item{LP}{library pathway scores, pathways x compounds.}
\item{sets}{gene sets named as \code{rownames(LP)}, on the genes of \code{LG}.}
\item{prm}{list with \code{n_gene}, \code{cut}, \code{n_top}, \code{spec_max}, \code{coh_min}, \code{coh_ratio}.}
\item{direction}{\code{"similar"} or \code{"reverse"}.}
\item{split}{split tier 1 into 1a (coherent) and 1b.}
}
\value{A data.frame, one row per compound, with attributes \code{conn} (connecting pathways) and \code{freq} (library frequency of each pathway in the top \code{n_top}).}
\examples{
toy <- parallax_toy()
k <- get_pathways(toy$pathways, "KEGG", rownames(toy$library))
lp <- map_pathways(toy$library, k)
qg <- make_signature(toy$expr, toy$group, 30, "meandiff")
qp <- make_signature(map_pathways(toy$expr, k), toy$group, 6, "meandiff")
prm <- list(n_gene = 30, cut = 0.1, n_top = 5, spec_max = 0.1, coh_min = 3, coh_ratio = 1.5)
r <- rank_list(qg[c("up", "dn")], qp[c("up", "dn")], toy$library, lp, k, prm, direction = "reverse")
head(r[order(r$G_rank), c("drug", "G_rank", "P_rank", "tier")])
}
