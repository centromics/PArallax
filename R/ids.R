id_type <- function(v) { v <- utils::head(stats::na.omit(as.character(v)), 1000); if (mean(grepl("^[0-9]+$", v)) > 0.8) "entrez" else "symbol" }
.VK <- new.env()
conv <- function(v, from, to) {
  v <- unique(stats::na.omit(as.character(v))); if (from == to || !length(v)) return(v)
  kt <- if (from == "entrez") "ENTREZID" else "SYMBOL"; cl <- if (to == "entrez") "ENTREZID" else "SYMBOL"
  if (is.null(.VK[[kt]])) .VK[[kt]] <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = kt)
  v <- v[v %in% .VK[[kt]]]; if (!length(v)) return(character(0))
  unique(stats::na.omit(unname(suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = v, column = cl,
                                                                      keytype = kt, multiVals = "first")))))
}
conv_map <- function(v, from, to) {
  if (from == to) return(stats::setNames(as.character(v), v))
  kt <- if (from == "entrez") "ENTREZID" else "SYMBOL"; cl <- if (to == "entrez") "ENTREZID" else "SYMBOL"
  suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(v), column = cl, keytype = kt, multiVals = "first"))
}
sets_on <- function(sets, univ, min_size = 5L, max_size = 500L) {
  from <- id_type(unlist(utils::head(sets, 20))); to <- id_type(univ)
  if (from != to) { gall <- unique(stats::na.omit(as.character(unlist(sets)))); mp <- stats::setNames(rep(NA_character_, length(gall)), gall)
    kt <- if (from == "entrez") "ENTREZID" else "SYMBOL"
    if (is.null(.VK[[kt]])) .VK[[kt]] <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = kt)
    ok <- gall[gall %in% .VK[[kt]]]; if (length(ok)) { m <- conv_map(ok, from, to); mp[ok] <- unname(m[ok]) }   # one lookup for all genes
    k <- lapply(sets, function(g) intersect(unique(stats::na.omit(unname(mp[as.character(g)]))), univ))
  } else k <- lapply(sets, function(g) intersect(unique(as.character(g)), univ))
  k[lengths(k) >= min_size & lengths(k) <= max_size]
}

#' Gene sets restricted to a gene universe
#'
#' \code{restrict_sets} keeps the genes present in \code{univ} (same identifiers) and the
#' sets of \code{min_size}-\code{max_size} genes. \code{convert_sets} first converts
#' between gene symbols and Entrez IDs (needs AnnotationDbi and org.Hs.eg.db).
#' @param sets named list of gene sets.
#' @param univ character vector of gene identifiers.
#' @param min_size,max_size set size limits after restriction.
#' @return named list of gene sets.
#' @export
restrict_sets <- function(sets, univ, min_size = 5L, max_size = 500L) {
  k <- lapply(sets, function(g) intersect(unique(as.character(g)), univ))
  k[lengths(k) >= min_size & lengths(k) <= max_size]
}

#' @rdname restrict_sets
#' @export
convert_sets <- function(sets, univ, min_size = 5L, max_size = 500L) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("convert_sets needs AnnotationDbi and org.Hs.eg.db")
  sets_on(sets, univ, min_size, max_size)
}
norm_sets <- function(r) {
  if (inherits(r, "GeneSetCollection")) return(GSEABase::geneIds(r))
  if (is.data.frame(r)) { g <- intersect(c("gene_symbol", "human_gene_symbol", "gene", "symbol"), names(r))[1]
    n <- intersect(c("gs_name", "geneset", "term", "pathway"), names(r))[1]; return(split(as.character(r[[g]]), r[[n]])) }
  if (is.list(r) && all(vapply(r, is.character, logical(1)))) return(r)
  unlist(unname(r), recursive = FALSE)
}
short <- function(x) { y <- trimws(sub("^(REACTOME|BIOCARTA|PID|WP)_", "", sub(" - Homo sapiens.*$", "", sub(" \\[KEGG.*\\]", "", x))))
  ifelse(nchar(y) == 0, x, y) }
