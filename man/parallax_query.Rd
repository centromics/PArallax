\name{parallax_query}
\alias{parallax_query}
\alias{print.parallax_result}
\alias{print.parallax_multi}
\title{Rank compounds for a two-group expression matrix}
\description{
Builds the query signature (genes and pathways) from case vs control samples, maps the query and the library to
pathway scores with the chosen \code{mapping} on the gene sets of one pathway database (\code{db}), scores every
compound with the weighted connectivity score (WTCS) on genes (G) and on pathways (P), and returns the tiers with the
connecting pathways of each compound. Tier 1 = G and P both in the top \code{cut} of the library; 1a / 1b split it by
coherence; 3 = G only; 2 / 2x = P only.
}
\usage{
parallax_query(expr, group, library, pathways = NULL, db = NULL,
              mapping = c("sgses_max", "avgz", "zscore", "ssgsea", "gsva",
                          "singscore", "ora_signed"),
              annotation = NULL, library_pw = NULL,
              direction = c("reverse", "similar"), signature = c("rankprod", "meandiff"),
              split = c("explore", "rule", "none"),
              n_gene = 100L, k_pw = 80L, cut = 0.10, n_top = 10L, spec_max = 0.10,
              coh_min = 3L, coh_ratio = 1.5, coh_p = 0.05, min_size = 5L, max_size = 500L,
              nperm = 100L, seed = 1L, explore = TRUE, n_draw = 1000L, explore_min = 10L)

\method{print}{parallax_result}(x, n = 20, ...)

\method{print}{parallax_multi}(x, n = 10, ...)
}
\arguments{
\item{expr}{numeric matrix, genes x samples, with the same gene identifiers as \code{library}.}
\item{group}{two levels, one per column of \code{expr}; first level = control (use a factor to set the order).}
\item{library}{numeric matrix, genes x compounds (e.g. LINCS signatures).}
\item{pathways}{gene sets: a named list of gene sets, a named list of collections (each a named list of gene sets),
a data.frame such as msigdbr output, or \code{NULL} to download \code{db}. See \code{\link{get_pathways}}.}
\item{db}{pathway database. Default (\code{pathways} and \code{db} both \code{NULL}): the latest KEGG, human, without the
disease pathways (hsa05...), downloaded from the KEGG REST API, as in the paper. Otherwise a collection name in \code{pathways},
\code{"KEGG"}, \code{"KEGG_ALL"} (with disease pathways), or an MSigDB C2 canonical-pathway subcollection fetched with msigdbr
(\code{"KEGG_MEDICUS"}, \code{"KEGG_LEGACY"}, \code{"REACTOME"}, \code{"BIOCARTA"}, \code{"PID"}, \code{"WIKIPATHWAYS"}).
Several names: each database is scored on its own (never merged) and a \code{parallax_multi} list is returned.}
\item{mapping}{pathway mapping of the query samples and the library: \code{"sgses_max"} (default; signed maximum deviation of the
weighted Kolmogorov-Smirnov running sum), \code{"avgz"}, \code{"zscore"}, \code{"ssgsea"}, \code{"gsva"} (needs GSVA),
\code{"singscore"}, \code{"ora_signed"}. The score is WTCS whatever the mapping.}
\item{annotation}{optional data.frame of compound information: your own table, LINCS pert_info, or
\code{signatureSearch::lincs_pert_info2} (install signatureSearch yourself). The ID column (\code{pert_id}, \code{id}, \code{drug},
\code{sig} or the first column) is matched to the library column names, their part before "__", or the BRD compound ID inside them;
name (\code{name}, \code{pert_iname}, \code{cmap_name}), \code{moa} and \code{target} columns are used when present.}
\item{library_pw}{optional precomputed library pathway scores (pathways x compounds, from \code{\link{map_pathways}} with the same mapping);
for several \code{db}, a list named by database.}
\item{direction}{\code{"reverse"} (compounds that reverse the case signature) or \code{"similar"}.}
\item{signature}{\code{"rankprod"} (single-study rank product, the MetaDE 1.0.5 algorithm used in the paper; built in) or \code{"meandiff"}.}
\item{split}{how tier 1 is split: \code{"explore"} (1a = change correlation inside the connecting pathways beyond random genes,
empirical P < \code{coh_p}; as in the senescence case), \code{"rule"} (1a = the genes behind G lie inside the connecting pathways,
\code{coh_min} genes and \code{coh_ratio} times the expected share; as in the same-compound benchmark), or \code{"none"}.
\code{TRUE} / \code{FALSE} mean \code{"rule"} / \code{"none"}.}
\item{n_gene, k_pw}{query genes and pathways per direction.}
\item{cut}{top fraction of the library for G and P.}
\item{n_top}{pathways per profile for connecting pathways.}
\item{spec_max}{a connecting pathway is uncommon if in the top \code{n_top} of fewer than this share of the library.}
\item{coh_min, coh_ratio}{thresholds of \code{split = "rule"}.}
\item{coh_p}{threshold of \code{split = "explore"}.}
\item{min_size, max_size}{pathway size limits.}
\item{nperm, seed}{rank-product permutations and seed.}
\item{explore}{compute the exploratory correlation for tier 1 (always on with \code{split = "explore"}).}
\item{n_draw, explore_min}{random draws and minimum genes of the exploratory correlation.}
\item{x}{a \code{parallax_result} or \code{parallax_multi}.}
\item{n}{rows to print.}
\item{...}{unused.}
}
\value{A \code{parallax_result}: \code{list} (tiers 1a, 1b, 1, then 3; G order inside each tier), \code{all}, \code{signature},
\code{library_pw}, \code{sets}, \code{change}, \code{change_lib}, \code{conn}, \code{freq}, \code{settings}.
For several \code{db}, a \code{parallax_multi} list of them, named by database.}
\details{The tiers rank candidates; they are not a measure of efficacy.}
\examples{
toy <- parallax_toy()
res <- parallax_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = "KEGG",
                     annotation = toy$annotation, signature = "meandiff",
                     n_gene = 30, k_pw = 6, n_top = 5, n_draw = 200, explore_min = 5)
res
}
