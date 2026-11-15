## WTCS from the up- and down-set enrichment scores
comb <- function(a, b) ifelse(a != 0 & b != 0 & sign(a) != sign(b), (a - b) / 2, 0)   # WTCS from ES_up, ES_down
pos_mat <- function(M) { P <- apply(-M, 2, rank, ties.method = "first"); storage.mode(P) <- "integer"; P }
## signed maximum deviation of the weighted Kolmogorov-Smirnov running sum of set idx, for every column of VAL
ks_cols <- function(POS, VAL, idx, s = 1) {
  N <- nrow(POS); M <- ncol(POS); idx <- idx[idx > 0L]; n <- length(idx)
  if (n == 0L || n >= N) return(rep(NA_real_, M))
  P <- POS[idx, , drop = FALSE]; W <- abs(VAL[idx, , drop = FALSE])^s
  o <- order(rep(seq_len(M), each = n), P); P <- matrix(P[o], n); W <- matrix(W[o], n)
  tot <- colSums(W); H <- if (n == 1L) W else apply(W, 2, cumsum); H <- H / rep(tot, each = n)
  miss <- (P - seq_len(n)) / (N - n)
  mx <- apply(H - miss, 2, max); mn <- apply(rbind(0, H[-n, , drop = FALSE]) - miss, 2, min)
  out <- ifelse(abs(mx) >= abs(mn), mx, mn); out[!is.finite(tot) | tot == 0] <- NA_real_; out
}
#' Weighted connectivity score (WTCS) of a query against a library
#'
#' \code{prep_lib} ranks every library profile once; \code{wtcs_query} scores a query
#' (up and down feature sets) against all profiles.
#' @param L numeric matrix, features x profiles (genes or pathways).
#' @param lib output of \code{prep_lib}.
#' @param up,dn character vectors: the query's up and down features.
#' @return \code{prep_lib}: a list; \code{wtcs_query}: named numeric vector of WTCS, one per profile.
#' @export
prep_lib <- function(L) list(L = L, POS = pos_mat(L), feat = rownames(L), ids = colnames(L))
#' @rdname prep_lib
#' @export
wtcs_query <- function(lib, up, dn) {
  a <- ks_cols(lib$POS, lib$L, match(up, lib$feat, 0L)); b <- ks_cols(lib$POS, lib$L, match(dn, lib$feat, 0L))
  stats::setNames(comb(as.numeric(a), as.numeric(b)), lib$ids)
}
## reference implementation (score_p(m, 0)), used in the tests
ref_wtcs <- function(L, up, dn) {
  one <- function(v, S) { v <- sort(v, decreasing = TRUE); N <- length(v); hit <- names(v) %in% S; n <- sum(hit)
    if (n == 0L || n == N) return(NA_real_)
    p <- numeric(N); p[hit] <- abs(v[hit]) / sum(abs(v[hit]))
    dd <- cumsum(p) - cumsum(!hit) / (N - n); dd[which.max(abs(dd))] }
  sapply(colnames(L), function(j) { v <- stats::setNames(L[, j], rownames(L)); comb(one(v, up), one(v, dn)) })
}
