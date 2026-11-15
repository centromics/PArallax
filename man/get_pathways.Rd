\name{get_pathways}
\alias{get_pathways}
\title{Pathway gene sets of one database, on the library's genes}
\description{
Picks one pathway database and returns its gene sets restricted to the gene universe of the library
(sets of \code{min_size}-\code{max_size} genes). Gene symbols and Entrez IDs are converted when the sets
and the library differ (needs AnnotationDbi and org.Hs.eg.db). Databases are used one at a time and never
merged: they are curated independently and overlap heavily. The default (\code{pathways = NULL}, \code{db = NULL}) is the latest
KEGG without disease pathways (hsa05...), downloaded from the KEGG REST API (free for academic use) and named as in the paper,
e.g. \code{"GLYCOLYSIS / GLUCONEOGENESIS [KEGG: hsa00010]"}; it is downloaded once per session. KEGG data are not shipped with
the package.
}
\usage{
get_pathways(pathways = NULL, db = NULL, univ, min_size = 5L, max_size = 500L)
}
\arguments{
\item{pathways}{one of: \code{NULL} (download \code{db}: KEGG from the KEGG REST API, the others from MSigDB with msigdbr);
a named list of gene sets (one collection; with \code{db}, the name pattern selects that database, e.g. \code{"REACTOME_"},
\code{"[KEGG: "}); a named list of collections, each a named list of gene sets (\code{db} names the collection);
a data.frame such as msigdbr output; or a \code{GeneSetCollection}.}
\item{db}{name of one database: a collection name in \code{pathways}; \code{"KEGG"} (default: latest KEGG, human,
disease pathways excluded, as in the paper) or \code{"KEGG_ALL"} (with them); or an MSigDB C2 canonical-pathway subcollection:
\code{"KEGG_MEDICUS"}, \code{"KEGG_LEGACY"}, \code{"REACTOME"}, \code{"BIOCARTA"}, \code{"PID"}, \code{"WIKIPATHWAYS"}.}
\item{univ}{gene identifiers of the library (\code{rownames(library)}).}
\item{min_size, max_size}{set size limits after restriction.}
}
\value{Named list of gene sets with attribute \code{db}.}
\examples{
toy <- parallax_toy()
k <- get_pathways(toy$pathways, "KEGG", rownames(toy$library))
length(k); attr(k, "db")
}
