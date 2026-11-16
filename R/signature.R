## single-study rank product, the algorithm of MetaDE 1.0.5 "rankProd" (Hong et al. 2006, Bioinformatics 22:2825),
## written out so that the package needs no MetaDE (archived on CRAN in 2018):
## every control-case sample pair gives a fold change; genes are ranked within each pair; RP = geometric mean of the ranks;
## null = the same after shuffling the genes within each sample, nperm times; P = share of the pooled null RP <= observed;
## FDR = Benjamini-Hochberg. AveFC = control mean - case mean (MetaDE sign: AveFC < 0 = higher in the case group).
## Same algorithm as the paper's MetaDE runs; the permutation draws differ, so FDR values can differ slightly.
rp_single <- function(x, y, nperm, seed) {
  x <- as.matrix(x); x[!is.finite(x)] <- NA
  i0 <- which(y == 0); i1 <- which(y == 1)
  pr <- expand.grid(a = i0, b = i1)
  rp <- function(X) {
    FC <- X[, pr$a, drop = FALSE] - X[, pr$b, drop = FALSE]                    # control - case, one column per pair
    ru <- apply(FC, 2, rank, ties.method = "average", na.last = "keep")        # rank 1 = most negative = most up in case
    rd <- apply(-FC, 2, rank, ties.method = "average", na.last = "keep")       # rank 1 = most down in case
    list(up = exp(rowMeans(log(ru), na.rm = TRUE)), dn = exp(rowMeans(log(rd), na.rm = TRUE)))
  }
  obs <- rp(x)
  NU <- ND <- matrix(NA_real_, nrow(x), nperm)
  withr::with_seed(seed, for (b in seq_len(nperm)) { r <- rp(apply(x, 2, sample)); NU[, b] <- r$up; ND[, b] <- r$dn })
  pnull <- function(o, N) { v <- sort(as.vector(N[is.finite(N)])); p <- findInterval(o, v) / length(v); p[!is.finite(o)] <- NA; p }
  fc <- rowMeans(x[, i0, drop = FALSE], na.rm = TRUE) - rowMeans(x[, i1, drop = FALSE], na.rm = TRUE)
  data.frame(feature = rownames(x), AveFC = unname(fc),
             FDR.up = stats::p.adjust(pnull(obs$up, NU), "BH"), FDR.down = stats::p.adjust(pnull(obs$dn, ND), "BH"),
             stringsAsFactors = FALSE)
}
## AveFC < 0 = up in the case group (MetaDE rank product)
select_topk <- function(tab, k) {                     # AveFC < 0 = up in senescence
  up <- tab[tab$AveFC < 0, ]; up <- up[order(up$FDR.up, up$AveFC), ]
  dn <- tab[tab$AveFC > 0, ]; dn <- dn[order(dn$FDR.down, -dn$AveFC), ]
  list(up = utils::head(up$feature, k), dn = utils::head(dn$feature, k))
}

#' Query signature from a two-group expression matrix
#'
#' @param x numeric matrix, features x samples.
#' @param group two levels, first = control.
#' @param n features per direction.
#' @param method \code{"rankprod"} (single-study rank product, the MetaDE 1.0.5 algorithm used in the paper) or \code{"meandiff"}.
#' @param nperm,seed permutations and seed of the rank product.
#' @return list: \code{up}, \code{dn}, \code{table}, \code{method}.
#' @export
make_signature <- function(x, group, n, method = c("rankprod", "meandiff"), nperm = 100L, seed = 1L) {
  method <- match.arg(method); x <- as.matrix(x); y <- as.integer(as.factor(group)) - 1L
  if (method == "rankprod") {
    tab <- rp_single(x, y, nperm, seed); s <- select_topk(tab, n)
  } else {
    d <- rowMeans(x[, y == 1, drop = FALSE], na.rm = TRUE) - rowMeans(x[, y == 0, drop = FALSE], na.rm = TRUE); d <- d[is.finite(d)]
    o <- order(d, decreasing = TRUE); tab <- data.frame(feature = names(d), diff = unname(d), stringsAsFactors = FALSE)
    s <- list(up = utils::head(names(d)[o][d[o] > 0], n), dn = utils::head(names(d)[rev(o)][d[rev(o)] < 0], n))
  }
  list(up = s$up, dn = s$dn, table = tab, method = method)
}

#' Exploratory check: query change vs compound change on the connecting-pathway genes
#'
#' Spearman correlation on genes \code{cg} against the same number of random genes.
#' @param cg genes.
#' @param z compound profile (named).
#' @param dchg query change (named).
#' @param n_draw random draws.
#' @param n_min minimum number of genes.
#' @param seed seed.
#' @param tail \code{"lower"}: P of a correlation this low or lower (reversal); \code{"upper"}: this high or higher (similarity).
#' @return named numeric: \code{n}, \code{rho}, \code{bg} (mean random), \code{P}.
#' @export
coh_cor <- function(cg, z, dchg, n_draw, n_min, seed, tail = c("lower", "upper")) {
  tail <- match.arg(tail)
  com <- intersect(names(z), names(dchg)); cg <- intersect(cg, com); n <- length(cg)
  if (n < n_min) return(c(n = n, rho = NA_real_, bg = NA_real_, P = NA_real_))
  rho <- stats::cor(dchg[cg], z[cg], method = "spearman"); dv <- dchg[com]; zv <- z[com]
  rr <- withr::with_seed(seed, replicate(n_draw, { k <- sample.int(length(com), n); stats::cor(dv[k], zv[k], method = "spearman") }))
  c(n = n, rho = rho, bg = mean(rr), P = (1 + sum(if (tail == "lower") rr <= rho else rr >= rho)) / (n_draw + 1))
}
