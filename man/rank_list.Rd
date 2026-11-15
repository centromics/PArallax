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
