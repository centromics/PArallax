\name{prep_lib}
\alias{prep_lib}
\alias{wtcs_query}
\title{Weighted connectivity score of a query against a library}
\description{\code{prep_lib} ranks every library profile once; \code{wtcs_query} scores a query (up and down sets) against all profiles.}
\usage{
prep_lib(L)

wtcs_query(lib, up, dn)
}
\arguments{
\item{L}{numeric matrix, features x profiles (genes or pathways).}
\item{lib}{output of \code{prep_lib}.}
\item{up, dn}{character vectors: the query's up and down features.}
}
\value{\code{prep_lib}: a list. \code{wtcs_query}: named numeric vector of WTCS, one per profile.}
