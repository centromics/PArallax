member_mat <- function(M, sets, scale) {
  A <- matrix(0, length(sets), nrow(M), dimnames = list(names(sets), rownames(M)))
  for (k in seq_along(sets)) { j <- match(sets[[k]], rownames(M)); j <- j[!is.na(j)]; if (length(j)) A[k, j] <- scale(length(j)) }
  A
}
tf_sgses_max <- function(M, sets) {
  POS <- pos_mat(M); idx <- lapply(sets, match, rownames(M))
  out <- t(vapply(idx, function(i) as.numeric(ks_cols(POS, M, i)), numeric(ncol(M))))
  dimnames(out) <- list(names(sets), colnames(M)); out
}
tf_avgz   <- function(M, sets) { M[!is.finite(M)] <- 0; member_mat(M, sets, function(n) 1 / n) %*% M }
tf_zscore <- function(M, sets) { M[!is.finite(M)] <- 0; member_mat(M, sets, function(n) 1 / sqrt(n)) %*% M }
tf_ssgsea <- function(M, sets, tau = 0.25) {
  N <- nrow(M); POS <- pos_mat(M); W <- (N - POS + 1)^tau
  out <- t(vapply(sets, function(g) { i <- match(g, rownames(M)); i <- i[!is.na(i)]; n <- length(i)
    if (n == 0L || n >= N) return(rep(NA_real_, ncol(M)))
    P <- POS[i, , drop = FALSE]; w <- W[i, , drop = FALSE]
    colSums(w * (N - P + 1)) / colSums(w) - (N * (N + 1) / 2 - colSums(N - P + 1)) / (N - n) }, numeric(ncol(M))))
  dimnames(out) <- list(names(sets), colnames(M)); out
}
tf_gsva <- function(M, sets) {
  M[!is.finite(M)] <- 0; sets <- lapply(sets, intersect, rownames(M))
  out <- if (exists("gsvaParam", asNamespace("GSVA")))
    GSVA::gsva(GSVA::gsvaParam(M, sets, kcdf = "Gaussian", minSize = 1, maxSize = Inf), verbose = FALSE) else
    GSVA::gsva(M, sets, method = "gsva", kcdf = "Gaussian", min.sz = 1, max.sz = Inf, verbose = FALSE)
  out <- as.matrix(out); out[names(sets)[names(sets) %in% rownames(out)], , drop = FALSE]
}
tf_singscore <- function(M, sets) {
  N <- nrow(M); R <- (N - pos_mat(M) + 1) / N
  out <- t(vapply(sets, function(g) { i <- match(g, rownames(M)); i <- i[!is.na(i)]
    if (!length(i)) return(rep(NA_real_, ncol(M))); colMeans(R[i, , drop = FALSE]) - (N + 1) / (2 * N) }, numeric(ncol(M))))
  dimnames(out) <- list(names(sets), colnames(M)); out
}
tf_ora_signed <- function(M, sets, tail = 0.05) {
  N <- nrow(M); POS <- pos_mat(M); t_n <- max(1L, round(tail * N)); TOP <- POS <= t_n; BOT <- POS > N - t_n
  out <- t(vapply(sets, function(g) { i <- match(g, rownames(M)); i <- i[!is.na(i)]; n <- length(i)
    if (!n) return(rep(NA_real_, ncol(M)))
    pu <- stats::phyper(colSums(TOP[i, , drop = FALSE]) - 1, n, N - n, t_n, lower.tail = FALSE, log.p = TRUE)
    pd <- stats::phyper(colSums(BOT[i, , drop = FALSE]) - 1, n, N - n, t_n, lower.tail = FALSE, log.p = TRUE)
    (pd - pu) / log(10) }, numeric(ncol(M))))
  dimnames(out) <- list(names(sets), colnames(M)); out
}
TF <- list(sgses_max = tf_sgses_max, avgz = tf_avgz, zscore = tf_zscore, ssgsea = tf_ssgsea, gsva = tf_gsva,
           singscore = tf_singscore, ora_signed = tf_ora_signed)

#' Map profiles to pathway scores
#'
#' Single-sample, signed pathway scores (positive = up). \code{sgses_max} is the signed
#' maximum deviation of the weighted Kolmogorov-Smirnov running sum (the PArallax default).
#' @param M numeric matrix, genes x profiles.
#' @param sets named list of gene sets with the identifiers of \code{rownames(M)}.
#' @param method mapping method.
#' @return numeric matrix, pathways x profiles, with attribute \code{mapping}.
#' @export
map_pathways <- function(M, sets, method = c("sgses_max", "avgz", "zscore", "ssgsea", "gsva", "singscore", "ora_signed")) {
  method <- match.arg(method)
  if (method == "gsva" && !requireNamespace("GSVA", quietly = TRUE)) stop("method = 'gsva' needs the GSVA package")
  out <- TF[[method]](as.matrix(M), sets); attr(out, "mapping") <- method; out
}
## ensemble of methods, weighted by their correlation on the library
scale_prof <- function(P) { P[!is.finite(P)] <- 0; s <- apply(P, 2, stats::sd); s[!is.finite(s) | s == 0] <- 1; sweep(P, 2, s, "/") }
make_weights <- function(SC, cor_cut) {
  R <- stats::cor(SC); pc <- stats::prcomp(SC, center = TRUE, scale. = TRUE)
  grp <- stats::cutree(stats::hclust(stats::as.dist(1 - R), method = "average"), h = 1 - cor_cut)
  w_cl <- 1 / (length(unique(grp)) * as.numeric(table(grp)[as.character(grp)])); names(w_cl) <- names(grp)
  l1 <- pc$rotation[, 1]; if (sum(l1) < 0) l1 <- -l1
  list(R = R, var = summary(pc)$importance[2, ], load = pc$rotation[, 1:min(3, ncol(SC))], group = grp,
       ENS_CLUSTER = w_cl, ENS_PC1 = l1 / sum(abs(l1)))
}
ens_mat <- function(Plist, w) Reduce(`+`, lapply(names(w), function(tf) w[[tf]] * scale_prof(Plist[[tf]])))
