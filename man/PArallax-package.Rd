\name{PArallax-package}
\alias{PArallax-package}
\alias{PArallax}
\docType{package}
\title{PArallax: connectivity scoring in pathway space on both sides}
\description{
Ranks library compounds against a query signature with the weighted connectivity score (WTCS)
computed twice by the same statistic: on genes (G) and on pathways (P), after both the query and the
library are mapped to pathway scores (sgses_max). G stays the base ranking; P adds evidence. Compounds
are placed in tiers (1a/1b or 1: G and P both near the top; 3: G only), each with the pathways that
connect it to the query. Main entry point: \code{\link{parallax_query}}.
}
\keyword{package}
