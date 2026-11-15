\name{plot_query}
\alias{plot_query}
\alias{panoply_figures}
\alias{print.panoply_figure}
\title{Figures of a panoply_result}
\description{
\code{plot_query} draws one figure: \code{"figure"} (all case panels on one page: a query signature in pathway space,
b every compound by G rank and P rank, c coherence of tier 1, d gene-level change of the 1a compounds inside their
connecting pathways, e connecting pathways of the 1a compounds; the format of Figure 7 of the paper),
\code{"list"} (the ranked list a researcher receives; the format of Figures 3c and 3d), or one panel
(\code{"signature"}, \code{"landscape"}, \code{"coherence"}, \code{"cases"}, \code{"pathways"}).
\code{panoply_figures} writes the list and the figure as PDF and PNG (600 dpi) at the printed width of 183 mm.
}
\usage{
plot_query(x, type = c("figure", "list", "landscape", "signature", "coherence", "cases", "pathways"),
           n = NULL, highlight = NULL, query_label = "query", n_case = 6L)

panoply_figures(x, prefix = "panoply", query_label = "query", highlight = NULL, n_list = 25L, n_case = 6L)

\method{print}{panoply_figure}(x, ...)
}
\arguments{
\item{x}{a \code{panoply_result} (\code{panoply_figures} also takes a \code{panoply_multi}: one pair of files per database);
for \code{print}, a \code{panoply_figure}.}
\item{type}{figure type.}
\item{n}{rows of the list (default 25) or pathways per direction of the signature (default \code{n_top}).}
\item{highlight}{compound IDs or names to highlight in the list.}
\item{query_label}{word for the query in titles and axes (e.g. \code{"senescence"}).}
\item{n_case}{compounds in the case panels (d, e).}
\item{prefix}{file name prefix (may include a folder).}
\item{n_list}{rows of the list.}
\item{...}{unused.}
}
\value{\code{plot_query}: a ggplot (a grob of class \code{panoply_figure} for \code{"figure"}; print it to draw).
\code{panoply_figures}: the file names, invisibly.}
\examples{
toy <- panoply_toy()
res <- panoply_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = "KEGG",
                     annotation = toy$annotation, signature = "meandiff",
                     n_gene = 30, k_pw = 6, n_top = 5, n_draw = 200, explore_min = 5)
plot_query(res, "list", highlight = "revertin")
plot_query(res, "figure", query_label = "toy case")
}
