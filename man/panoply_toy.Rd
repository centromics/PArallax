\name{panoply_toy}
\alias{panoply_toy}
\title{A small simulated data set}
\description{
600 genes, 8 samples (4 control, 4 case), a library of 150 compound profiles and two pathway collections
(\code{KEGG}: 40 sets, \code{REACTOME}: 60 sets; invented gene sets with these names, for trying \code{db}).
The case raises the genes of KEGG pathways 1-3 and lowers those of pathways 4-6. Seven compounds are built in
(\code{truth}): four reversers of different breadth and strength, one that reverses a few scattered genes strongly,
one that shifts two pathways weakly, and one that mimics the case. The other compounds are noise; a few move
unrelated pathways. A longer walk-through: \code{system.file("examples", "toy_example.R", package = "PAnoply")}.
}
\usage{
panoply_toy(seed = 1L)
}
\arguments{
\item{seed}{random seed.}
}
\value{List: \code{expr} (genes x samples), \code{group} (factor, control first), \code{library} (genes x compounds),
\code{pathways} (list of two collections), \code{annotation} (\code{id}, \code{name}, \code{moa}), \code{truth}.}
\examples{
toy <- panoply_toy()
toy$truth
}
