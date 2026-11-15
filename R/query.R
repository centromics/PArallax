#' Rank compounds for a two-group expression matrix
#'
#' Builds the query signature (genes and pathways) from case vs control samples, maps the query and the library
#' to pathway scores with the chosen \code{mapping} on the gene sets of one pathway database (\code{db}), scores
#' every compound with the weighted connectivity score (WTCS) on genes (G) and on pathways (P), and returns the
#' tiers with the connecting pathways of each compound. Tier 1 = G and P both in the top \code{cut} of the library;
#' 1a / 1b split it by coherence; 3 = G only; 2 / 2x = P only.
#' @param expr numeric matrix, genes x samples, same gene identifiers as \code{library}.
#' @param group two levels, one per column of \code{expr}; first level = control (use a factor to set the order).
#' @param library numeric matrix, genes x compounds (e.g. LINCS signatures).
#' @param pathways gene sets: a named list of gene sets, a named list of collections (each a named list of gene sets),
#'   a data.frame such as msigdbr output, or \code{NULL} to download \code{db}. See \code{\link{get_pathways}}.
#' @param db pathway database. Default (\code{pathways} and \code{db} both \code{NULL}): the latest KEGG, human, without the
#'   disease pathways (hsa05...), downloaded from the KEGG REST API, as in the paper. Otherwise a collection name in \code{pathways},
#'   \code{"KEGG"}, \code{"KEGG_ALL"} (with disease pathways), or an MSigDB C2 canonical-pathway subcollection fetched with msigdbr
#'   (\code{"KEGG_MEDICUS"}, \code{"KEGG_LEGACY"}, \code{"REACTOME"}, \code{"BIOCARTA"}, \code{"PID"}, \code{"WIKIPATHWAYS"}).
#'   Several names: each database is scored on its own (never merged) and a \code{panoply_multi} list is returned.
#' @param mapping pathway mapping of the query samples and the library: \code{"sgses_max"} (default; signed maximum deviation of the
#'   weighted Kolmogorov-Smirnov running sum), \code{"avgz"}, \code{"zscore"}, \code{"ssgsea"}, \code{"gsva"} (needs GSVA),
#'   \code{"singscore"}, \code{"ora_signed"}. The score is WTCS whatever the mapping.
#' @param annotation optional data.frame of compound information: your own table, LINCS pert_info, or
#'   \code{signatureSearch::lincs_pert_info2} (install signatureSearch yourself). The ID column (\code{pert_id}, \code{id}, \code{drug},
#'   \code{sig} or the first column) is matched to the library column names, their part before "__", or the BRD compound ID inside them;
#'   name (\code{name}, \code{pert_iname}, \code{cmap_name}), \code{moa} and \code{target} columns are used when present.
#' @param library_pw optional precomputed library pathway scores (pathways x compounds, from \code{\link{map_pathways}} with the same mapping);
#'   for several \code{db}, a list named by database.
#' @param direction \code{"reverse"} (compounds that reverse the case signature) or \code{"similar"}.
#' @param signature \code{"rankprod"} (single-study rank product, the MetaDE 1.0.5 algorithm used in the paper; built in) or \code{"meandiff"}.
#' @param split how tier 1 is split: \code{"explore"} (1a = change correlation inside the connecting pathways beyond random genes,
#'   empirical P < \code{coh_p}; as in the senescence case), \code{"rule"} (1a = the genes behind G lie inside the connecting pathways,
#'   \code{coh_min} genes and \code{coh_ratio} times the expected share; as in the same-compound benchmark), or \code{"none"}.
#'   \code{TRUE} / \code{FALSE} mean \code{"rule"} / \code{"none"}.
#' @param n_gene,k_pw query genes and pathways per direction.
#' @param cut top fraction of the library for G and P.
#' @param n_top pathways per profile for connecting pathways.
#' @param spec_max a connecting pathway is uncommon if in the top \code{n_top} of fewer than this share of the library.
#' @param coh_min,coh_ratio thresholds of \code{split = "rule"}.
#' @param coh_p threshold of \code{split = "explore"}.
#' @param min_size,max_size pathway size limits.
#' @param nperm,seed rank-product permutations and seed.
#' @param explore compute the exploratory correlation for tier 1 (always on with \code{split = "explore"}).
#' @param n_draw,explore_min random draws and minimum genes of the exploratory correlation.
#' @param x a \code{panoply_result} or \code{panoply_multi}.
#' @param n rows to print.
#' @param ... unused.
#' @return a \code{panoply_result}: \code{list} (tiers 1a, 1b, 1, then 3; G order inside each tier), \code{all}, \code{signature},
#'   \code{library_pw}, \code{sets}, \code{change}, \code{change_lib}, \code{conn}, \code{freq}, \code{settings};
#'   for several \code{db}, a \code{panoply_multi} list of them, named by database.
#' @export
panoply_query <- function(expr, group, library, pathways = NULL, db = NULL,
                          mapping = c("sgses_max", "avgz", "zscore", "ssgsea", "gsva", "singscore", "ora_signed"),
                          annotation = NULL, library_pw = NULL,
                          direction = c("reverse", "similar"), signature = c("rankprod", "meandiff"),
                          split = c("explore", "rule", "none"),
                          n_gene = 100L, k_pw = 80L, cut = 0.10, n_top = 10L, spec_max = 0.10,
                          coh_min = 3L, coh_ratio = 1.5, coh_p = 0.05, min_size = 5L, max_size = 500L,
                          nperm = 100L, seed = 1L, explore = TRUE, n_draw = 1000L, explore_min = 10L) {
  ## several databases: each one on its own, never merged
  if (length(db) > 1L) {
    if (!is.null(library_pw) && (is.matrix(library_pw) || is.data.frame(library_pw)))
      warning("one library_pw matrix given for several databases: ignored (give a list named by database)")
    cl <- match.call(); out <- list()
    for (d in db) {
      cl$db <- d
      cl$library_pw <- if (is.list(library_pw) && !is.data.frame(library_pw) && !is.matrix(library_pw)) library_pw[[d]] else NULL
      out[[d]] <- eval(cl, parent.frame())
    }
    return(structure(out, class = "panoply_multi"))
  }
  mapping <- match.arg(mapping); direction <- match.arg(direction); signature <- match.arg(signature)
  if (is.logical(split)) split <- if (isTRUE(split)) "rule" else "none"
  split <- match.arg(split, c("explore", "rule", "none"))
  if (split == "explore") explore <- TRUE
  if (mapping == "gsva" && !requireNamespace("GSVA", quietly = TRUE)) stop("mapping = 'gsva' needs the GSVA package")
  expr <- as.matrix(expr); library <- as.matrix(library)
  y <- as.integer(as.factor(group)) - 1L
  if (length(y) != ncol(expr) || any(!y %in% 0:1)) stop("group needs two levels, one value per column of expr (first level = control)")
  if (sum(y == 0) < 2 || sum(y == 1) < 2) stop("need at least 2 control and 2 case samples")
  if (!length(intersect(rownames(expr), rownames(library)))) stop("expr and library share no gene identifiers")
  library[!is.finite(library)] <- 0; x0 <- expr; x0[!is.finite(x0)] <- 0

  ## pathways of one database, mapping of the library and of the query samples
  sets <- get_pathways(pathways, db, rownames(library), min_size, max_size); db_lab <- attr(sets, "db")
  ks <- restrict_sets(sets, rownames(x0), min_size, Inf)
  if (!length(ks)) stop("no pathway of ", db_lab, " has ", min_size, " or more genes in expr")
  if (is.null(library_pw)) LP <- TF[[mapping]](library, sets) else {
    m0 <- attr(library_pw, "mapping")
    if (!is.null(m0) && m0 != mapping) stop("library_pw was mapped with '", m0, "', not '", mapping, "'")
    LP <- as.matrix(library_pw)
    if (!all(colnames(library) %in% colnames(LP))) stop("library_pw does not cover every column of library")
  }
  PS <- TF[[mapping]](x0, ks)
  LP[!is.finite(LP)] <- 0; PS[!is.finite(PS)] <- 0
  pw <- intersect(intersect(rownames(LP), names(sets)), rownames(PS))
  if (length(pw) < 2L * n_top) warning("only ", length(pw), " pathways of ", db_lab, " are usable; n_top = ", n_top, " per direction")
  LP <- LP[pw, colnames(library), drop = FALSE]; sets <- sets[pw]

  ## query signature in both spaces
  sg <- make_signature(x0, y, n_gene, signature, nperm, seed)
  sp <- make_signature(PS[pw, , drop = FALSE], y, k_pw, signature, nperm, seed)
  md <- rowMeans(expr[, y == 1, drop = FALSE], na.rm = TRUE) - rowMeans(expr[, y == 0, drop = FALSE], na.rm = TRUE)
  agree <- mean(c(md[sg$up] > 0, md[sg$dn] < 0), na.rm = TRUE)
  if (!is.finite(agree) || agree < 0.8) warning("only ", round(agree, 2), " of the signature genes point the case-minus-control way; check group coding")

  ## G, P, tiers
  prm <- list(n_gene = n_gene, cut = cut, n_top = n_top, spec_max = spec_max, coh_min = coh_min, coh_ratio = coh_ratio)
  R <- rank_list(sg[c("up", "dn")], sp[c("up", "dn")], library, LP, sets, prm, direction, split == "rule")
  conn <- attr(R, "conn"); freq <- attr(R, "freq")
  R <- annotate_compounds(R, annotation)

  ## exploratory correlation inside the connecting pathways (tier 1)
  dchg <- md[is.finite(md)]; tail <- if (direction == "reverse") "lower" else "upper"
  R$coh_genes <- NA_real_; R$coh_rho <- NA_real_; R$coh_bg <- NA_real_; R$coh_P <- NA_real_
  if (explore) for (i in which(R$tier %in% c("1", "1a", "1b"))) { j <- R$drug[i]
    cg <- intersect(unique(unlist(sets[c(conn[[j]]$up, conn[[j]]$dn)])), rownames(library))
    v <- coh_cor(cg, library[, j], dchg, n_draw, explore_min, seed + match(j, colnames(library)), tail)
    R$coh_genes[i] <- v[["n"]]; R$coh_rho[i] <- v[["rho"]]; R$coh_bg[i] <- v[["bg"]]; R$coh_P[i] <- v[["P"]] }
  if (split == "explore") {
    t1 <- R$tier == "1"
    beyond <- if (direction == "reverse") R$coh_rho < R$coh_bg else R$coh_rho > R$coh_bg
    pass <- t1 & !is.na(beyond) & beyond & !is.na(R$coh_P) & R$coh_P < coh_p
    R$tier[t1] <- ifelse(pass[t1], "1a", "1b")
  }

  ord <- c("1a", "1b", "1", "3")
  L <- R[R$tier %in% ord, , drop = FALSE]; L <- L[order(match(L$tier, ord), L$G_rank), , drop = FALSE]; L$priority <- seq_len(nrow(L))
  t1 <- L$drug[L$tier %in% c("1a", "1b", "1")]
  structure(list(list = L, all = R, signature = list(genes = sg, pathways = sp, agreement = agree), library_pw = LP, sets = sets,
                 change = dchg, change_lib = lapply(stats::setNames(t1, t1), function(j) library[, j]), conn = conn, freq = freq,
                 settings = c(prm, list(direction = direction, signature = signature, split = split, coh_p = coh_p, k_pw = k_pw,
                                        nperm = nperm, seed = seed, explore = explore, n_draw = n_draw, explore_min = explore_min,
                                        score = "WTCS", mapping = mapping, db = db_lab, n_pathways = length(pw)))),
            class = "panoply_result")
}

## compound names, MOA and targets from an annotation table (own table, LINCS pert_info, or signatureSearch::lincs_pert_info2)
## ID: first matching column of pert_id / id / drug / sig (else the first column); matched on the full library column name,
## its part before "__", or the BRD compound ID inside it (BRD-X00000000)
.acol <- function(a, cand) { k <- match(tolower(cand), tolower(names(a))); k <- k[!is.na(k)]; if (length(k)) names(a)[k[1]] else NA_character_ }
.brd <- function(x) { m <- regmatches(x, regexpr("BRD-[A-Z][0-9]{8}", x)); out <- rep(NA_character_, length(x)); out[grepl("BRD-[A-Z][0-9]{8}", x)] <- m; out }
annotate_compounds <- function(R, annotation) {
  R$name <- R$drug; R$moa <- NA_character_; R$target <- NA_character_
  if (is.null(annotation)) return(R)
  a <- as.data.frame(annotation, stringsAsFactors = FALSE)
  idc <- .acol(a, c("pert_id", "id", "drug", "sig", "sig_id")); if (is.na(idc)) idc <- names(a)[1]
  ids <- as.character(a[[idc]])
  i <- match(R$drug, ids)
  j <- match(sub("__.*$", "", R$drug), ids); i[is.na(i)] <- j[is.na(i)]
  j <- match(.brd(R$drug), .brd(ids)); i[is.na(i)] <- j[is.na(i)]
  if (all(is.na(i))) { warning("annotation: no compound ID of the library found in column '", idc, "'"); return(R) }
  pick <- function(cand) { cc <- .acol(a, cand); if (is.na(cc)) return(rep(NA_character_, nrow(R)))
    v <- a[[cc]]; if (is.list(v)) v <- vapply(v, function(e) paste(unique(stats::na.omit(as.character(e))), collapse = ", "), "")
    as.character(v)[i] }
  nm <- pick(c("name", "pert_iname", "cmap_name", "pert_name", "compound"))
  R$name <- ifelse(is.na(nm) | nm == "" | nm == "-666", R$drug, nm)
  R$moa <- pick(c("moa", "mechanism")); R$target <- pick(c("target", "targets"))
  R
}

#' @rdname panoply_query
#' @export
print.panoply_result <- function(x, n = 20, ...) {
  s <- x$settings; tb <- table(factor(x$all$tier, levels = c("1a", "1b", "1", "3", "2", "2x", "-"))); tb <- tb[tb > 0]
  cat("PAnoply result:", nrow(x$all), "library compounds |", s$direction, "| signature:", s$signature, "| score WTCS\n")
  cat("pathways:", s$db, "(", s$n_pathways, "sets ) | mapping:", s$mapping, "| tier-1 split:", s$split, "\n")
  cat("query:", length(x$signature$genes$up), "+", length(x$signature$genes$dn), "genes,", length(x$signature$pathways$up), "+",
      length(x$signature$pathways$dn), "pathways | direction agreement", round(x$signature$agreement, 2), "\n")
  cat("tiers:", paste(names(tb), tb, sep = " ", collapse = " | "), "\n\n")
  print(utils::head(x$list[, c("priority", "tier", "name", "G_rank", "P_rank", "specific_pathways")], n), row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @rdname panoply_query
#' @export
print.panoply_multi <- function(x, n = 10, ...) {
  cat("PAnoply results for", length(x), "pathway databases, each scored on its own (never merged):", paste(names(x), collapse = ", "), "\n")
  for (d in names(x)) { cat("\n=====", d, "=====\n"); print(x[[d]], n = n) }
  invisible(x)
}
