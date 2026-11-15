#' A small simulated data set
#'
#' 600 genes, 8 samples (4 control, 4 case), a library of 150 compound profiles and two pathway collections
#' (\code{KEGG}: 40 sets, \code{REACTOME}: 60 sets; invented gene sets with these names, for trying \code{db}).
#' The case raises the genes of KEGG pathways 1-3 and lowers those of pathways 4-6. Seven compounds are built in
#' (\code{truth}): four reversers of different breadth and strength, one that reverses a few scattered genes strongly,
#' one that shifts two pathways weakly, and one that mimics the case. The other compounds are noise; a few move
#' unrelated pathways.
#' @param seed random seed.
#' @return list: \code{expr} (genes x samples), \code{group} (factor, control first), \code{library} (genes x compounds),
#'   \code{pathways} (list of two collections), \code{annotation} (\code{id}, \code{name}, \code{moa}), \code{truth}.
#' @examples
#' toy <- parallax_toy()
#' res <- parallax_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = "KEGG",
#'                      annotation = toy$annotation, signature = "meandiff",
#'                      n_gene = 30, k_pw = 6, n_top = 5, n_draw = 200, explore_min = 5)
#' res
#' @export
parallax_toy <- function(seed = 1L) {
  set.seed(seed)
  genes <- sprintf("G%03d", 1:600)
  mk <- function(prefix, n) stats::setNames(lapply(seq_len(n), function(i) sample(genes, sample(15:40, 1))), sprintf("%s_PATHWAY_%02d", prefix, seq_len(n)))
  pw <- list(KEGG = mk("KEGG", 40), REACTOME = mk("REACTOME", 60))
  K <- pw$KEGG; up <- unique(unlist(K[1:3])); dn <- setdiff(unique(unlist(K[4:6])), up)
  eff <- stats::setNames(rep(0, length(genes)), genes); eff[up] <- 1.5; eff[dn] <- -1.5
  noise <- function(n) matrix(stats::rnorm(length(genes) * n), length(genes), dimnames = list(genes, NULL))
  expr <- cbind(noise(4), noise(4) + eff); colnames(expr) <- c(paste0("control", 1:4), paste0("case", 1:4))
  group <- factor(rep(c("control", "case"), each = 4), levels = c("control", "case"))
  lib <- noise(150); colnames(lib) <- sprintf("TOY%03d", 1:150)
  lib[, "TOY001"] <- lib[, "TOY001"] - 1.2 * eff                               # reverses the whole case signature
  lib[, "TOY002"] <- lib[, "TOY002"] - 1.5 * pmax(eff, 0)                      # reverses the raised pathways only
  lib[, "TOY003"] <- lib[, "TOY003"] - 1.5 * pmin(eff, 0)                      # reverses the lowered pathways only
  lib[, "TOY004"] <- lib[, "TOY004"] - 0.6 * eff                               # weak full reversal
  few <- sample(names(eff)[eff != 0], 15)
  lib[few, "TOY005"] <- lib[few, "TOY005"] - 3 * sign(eff[few])                # a few scattered genes, strongly
  sh <- unique(unlist(K[c(2, 5)]))
  lib[sh, "TOY006"] <- lib[sh, "TOY006"] - 0.5 * sign(eff[sh])                 # two pathways, weakly
  lib[, "TOY007"] <- lib[, "TOY007"] + 1.2 * eff                               # mimics the case
  for (j in sprintf("TOY%03d", 8:15)) { g <- unlist(K[sample(10:40, 2)]); lib[g, j] <- lib[g, j] + 1.2 }   # unrelated pathways
  role <- c("reverser, full", "reverser, raised pathways", "reverser, lowered pathways", "reverser, weak",
            "few genes, strong", "two pathways, weak", "mimic")
  nm <- c("revertin", "upsetrol", "downstatin", "mildocin", "spotinib", "pathomab", "mimicane")
  moa <- sample(c("kinase inhibitor", "HDAC inhibitor", "GPCR agonist", "ion channel blocker", "tubulin binder", "HSP90 inhibitor"), 150, TRUE)
  ann <- data.frame(id = colnames(lib), name = c(nm, sprintf("compound %03d", 8:150)), moa = moa, stringsAsFactors = FALSE)
  list(expr = expr, group = group, library = lib, pathways = pw, annotation = ann,
       truth = data.frame(id = sprintf("TOY%03d", 1:7), name = nm, role = role, stringsAsFactors = FALSE))
}
