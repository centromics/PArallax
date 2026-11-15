\name{map_pathways}
\alias{map_pathways}
\title{Map profiles to pathway scores}
\description{Single-sample, signed pathway scores (positive = up). \code{"sgses_max"}, the PArallax default, is the signed
maximum deviation of the weighted Kolmogorov-Smirnov running sum of the pathway's genes.}
\usage{
map_pathways(M, sets, method = c("sgses_max", "avgz", "zscore", "ssgsea", "gsva", "singscore", "ora_signed"))
}
\arguments{
\item{M}{numeric matrix, genes x profiles.}
\item{sets}{named list of gene sets with the identifiers of \code{rownames(M)}.}
\item{method}{mapping method; \code{"gsva"} needs GSVA.}
}
\value{Numeric matrix, pathways x profiles, with attribute \code{mapping} (checked by \code{\link{parallax_query}} when given as \code{library_pw}).}
