#' PAnoply: connectivity scoring in pathway space on both sides
#'
#' Ranks library compounds against a query signature with two scores computed by the
#' same statistic, the weighted connectivity score (WTCS): G on genes and P on pathways,
#' after both the query and the library are mapped to pathway scores (seven mapping methods,
#' sgses_max by default; one pathway database at a time, never merged).
#' G stays the base ranking; P adds evidence. Compounds are placed in tiers
#' (1a/1b or 1: G and P both near the top; 3: G only), with the pathways that connect each
#' compound to the query.
#'
#' Main entry point: \code{\link{panoply_query}}.
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("G_rank", "P_rank", "grp", "name", "v", "lab", "dir", "plab", "score", "spec", "connecting",
                         "d", "z", "pdir", "case", "coh_bg", "coh_rho", "priority", "G_pos", "is_X", "y", "x", "hj", "tl",
                         "sym", "hit"))
