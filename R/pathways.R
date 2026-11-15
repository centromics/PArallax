## ---- pathway databases ---------------------------------------------------------------------------------
## Default: the latest KEGG pathways (human), downloaded from the KEGG REST API at run time and named as in the paper
## ("GLYCOLYSIS / GLUCONEOGENESIS [KEGG: hsa00010]"), disease pathways (hsa05...) excluded. KEGG data are not shipped
## with the package (KEGG copyright). Other databases: MSigDB C2 canonical-pathway subcollections fetched with msigdbr,
## or the user's own gene sets. Databases are always used one at a time, never merged.
.MSIG <- data.frame(db     = c("KEGG_MEDICUS", "KEGG_LEGACY", "REACTOME", "BIOCARTA", "PID", "WIKIPATHWAYS"),
                    sub    = c("CP:KEGG_MEDICUS", "CP:KEGG_LEGACY", "CP:REACTOME", "CP:BIOCARTA", "CP:PID", "CP:WIKIPATHWAYS"),
                    prefix = c("^KEGG_MEDICUS_", "^KEGG_(?!MEDICUS_)", "^REACTOME_", "^BIOCARTA_", "^PID_", "^WP_"),
                    stringsAsFactors = FALSE)
.KEGG_DB <- c("KEGG", "KEGG_ALL")                       # latest KEGG without / with the disease pathways
.db_key <- function(db) {
  d <- toupper(gsub("[ -]", "_", db))
  if (d %in% c("WP", "WIKIPATHWAY")) d <- "WIKIPATHWAYS"
  if (!d %in% c(.KEGG_DB, .MSIG$db)) stop("unknown db '", db, "'; one of: ", paste(c(.KEGG_DB, .MSIG$db), collapse = ", "),
                                          " (or a name of your own collections)")
  d
}
## latest KEGG pathways from the KEGG REST API (free for academic use; https://www.kegg.jp/kegg/legal.html)
## the same source and naming as the paper's pathway list (clusterProfiler prepare_KEGG + rest.kegg.jp/list/pathway)
.kegg_latest <- function(exclude_disease = TRUE, org = "hsa") {
  key <- paste0("kegg_", org)
  if (is.null(.VK[[key]])) {
    message("downloading the latest KEGG pathways (", org, ") from rest.kegg.jp - free for academic use, see https://www.kegg.jp/kegg/legal.html")
    rd <- function(u) tryCatch(utils::read.delim(u, header = FALSE, colClasses = "character", quote = "", comment.char = ""),
                               error = function(e) stop("KEGG download failed (", conditionMessage(e), "); give `pathways` yourself"))
    lk <- rd(paste0("https://rest.kegg.jp/link/", org, "/pathway"))
    ls <- rd(paste0("https://rest.kegg.jp/list/pathway/", org))
    pid <- sub("^path:", "", lk[[1]]); gid <- sub(paste0("^", org, ":"), "", lk[[2]])
    desc <- stats::setNames(ls[[2]], sub("^path:", "", ls[[1]]))
    sets <- lapply(split(gid, pid), unique)
    names(sets) <- paste0(toupper(sub(" - [^-]*\\([^)]*\\)$", "", desc[names(sets)])), " [KEGG: ", names(sets), "]")
    attr(sets, "downloaded") <- format(Sys.Date())
    .VK[[key]] <- sets
  }
  sets <- .VK[[key]]
  if (exclude_disease) sets <- sets[!grepl(paste0("\\[KEGG: ", org, "05"), names(sets))]
  sets
}
.df_sets <- function(df, id) {
  gcol <- if (id == "entrez") intersect(c("ncbi_gene", "entrez_gene", "gene"), names(df)) else
    intersect(c("gene_symbol", "human_gene_symbol", "symbol", "gene"), names(df))
  ncol_ <- intersect(c("gs_name", "geneset", "term", "pathway"), names(df))
  if (!length(gcol) || !length(ncol_)) stop("a data.frame of gene sets needs a set-name column (gs_name) and a gene column")
  split(as.character(df[[gcol[1]]]), as.character(df[[ncol_[1]]]))
}
.msig_sets <- function(db, id) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("db = '", db, "' without `pathways` needs the msigdbr package")
  sub <- .MSIG$sub[.MSIG$db == db]
  df <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = sub),          # msigdbr >= 10
                 error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = sub))   # older msigdbr
  if (!nrow(df)) stop("msigdbr returned no gene set for ", sub, " (older msigdbr versions name the subcollections differently)")
  .df_sets(as.data.frame(df), id)
}
.as_sets <- function(x, id) {
  if (inherits(x, "GeneSetCollection")) return(GSEABase::geneIds(x))
  if (is.data.frame(x)) return(.df_sets(x, id))
  x
}

#' Pathway gene sets of one database, on the library's genes
#'
#' Picks one pathway database and returns its gene sets restricted to the gene universe of the
#' library (sets of \code{min_size}-\code{max_size} genes). Gene symbols and Entrez IDs are converted
#' when the sets and the library differ (needs AnnotationDbi and org.Hs.eg.db).
#' Databases are used one at a time and never merged: they are curated independently and overlap heavily.
#' The default (\code{pathways = NULL}, \code{db = NULL}) is the latest KEGG without disease pathways, as in the paper.
#' @param pathways one of: \code{NULL} (download \code{db}: KEGG from the KEGG REST API, the others from MSigDB with msigdbr);
#'   a named list of gene sets (one collection; with \code{db}, the name pattern selects that database, e.g. \code{"REACTOME_"},
#'   \code{"[KEGG: "}); a named list of collections, each a named list of gene sets (\code{db} names the collection);
#'   a data.frame such as msigdbr output; or a \code{GeneSetCollection}.
#' @param db name of one database: a collection name in \code{pathways}; \code{"KEGG"} (default: latest KEGG, human,
#'   disease pathways hsa05... excluded, as in the paper) or \code{"KEGG_ALL"} (with them); or an MSigDB C2 canonical-pathway
#'   subcollection: \code{"KEGG_MEDICUS"}, \code{"KEGG_LEGACY"}, \code{"REACTOME"}, \code{"BIOCARTA"}, \code{"PID"}, \code{"WIKIPATHWAYS"}.
#' @param univ gene identifiers of the library (\code{rownames(library)}).
#' @param min_size,max_size set size limits after restriction.
#' @return named list of gene sets with attribute \code{db}.
#' @export
get_pathways <- function(pathways = NULL, db = NULL, univ, min_size = 5L, max_size = 500L) {
  if (length(db) > 1L) stop("get_pathways takes one db at a time (databases are never merged)")
  id <- id_type(univ)
  if (is.null(pathways)) {
    label <- if (is.null(db)) "KEGG" else .db_key(db)
    sets <- if (label %in% .KEGG_DB) .kegg_latest(exclude_disease = label == "KEGG") else .msig_sets(label, id)
  } else {
    pathways <- .as_sets(pathways, id)
    if (!is.list(pathways) || is.null(names(pathways))) stop("`pathways` must be a named list of gene sets (or of collections)")
    nested <- all(vapply(pathways, function(e) is.list(e) || is.data.frame(e) || inherits(e, "GeneSetCollection"), logical(1)))
    if (nested) {
      if (is.null(db)) {
        if (length(pathways) > 1L) stop("`pathways` holds several collections (", paste(names(pathways), collapse = ", "), "); choose one with db")
        db <- names(pathways)
      }
      i <- match(toupper(db), toupper(names(pathways)))
      if (is.na(i)) stop("db '", db, "' is not a collection of `pathways` (", paste(names(pathways), collapse = ", "), ")")
      label <- names(pathways)[i]; sets <- .as_sets(pathways[[i]], id)
    } else if (!is.null(db)) {
      label <- .db_key(db)
      pat <- if (label %in% .KEGG_DB) "\\[KEGG: |^KEGG_(?!MEDICUS_)" else .MSIG$prefix[.MSIG$db == label]
      sets <- pathways[grepl(pat, names(pathways), perl = TRUE)]
      if (label == "KEGG") sets <- sets[!grepl("\\[KEGG: [a-z]+05", names(sets))]
      if (!length(sets)) stop("no gene set of `pathways` has a ", label, " name")
    } else { label <- "custom"; sets <- pathways }
  }
  sets <- lapply(sets, function(g) unique(as.character(g[!is.na(g)])))
  sets <- sets[lengths(sets) > 0]
  from <- id_type(unlist(utils::head(sets, 20)))
  out <- if (from == id) restrict_sets(sets, univ, min_size, max_size) else convert_sets(sets, univ, min_size, max_size)
  if (!length(out)) stop("no pathway of ", label, " has ", min_size, "-", max_size, " genes in the library")
  attr(out, "db") <- label
  out
}
