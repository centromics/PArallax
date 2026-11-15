q_top <- function(v, n) { v <- v[is.finite(v)]; o <- order(v, decreasing = TRUE)
  list(up = utils::head(names(v)[o], n), dn = utils::head(names(v)[rev(o)], n)) }
q_pw <- function(v, k) { v[!is.finite(v)] <- 0; o <- order(v, decreasing = TRUE)
  list(up = utils::head(names(v)[o][v[o] > 0], k), dn = utils::head(names(v)[rev(o)][v[rev(o)] < 0], k)) }
nrk <- function(sc) (rank(-sc, ties.method = "first") - 1) / (length(sc) - 1)
top_ind <- function(M, n, up = TRUE) {
  out <- matrix(FALSE, nrow(M), ncol(M), dimnames = dimnames(M))
  for (j in seq_len(ncol(M))) { v <- M[, j]; o <- if (up) order(v, decreasing = TRUE) else order(v)
    o <- o[if (up) v[o] > 0 else v[o] < 0]; o <- o[!is.na(o)]; out[utils::head(o, n), j] <- TRUE }
  out
}
jac <- function(a, b) if (!length(union(a, b))) NA_real_ else length(intersect(a, b)) / length(union(a, b))
## the genes that make the G connection lie inside the connecting pathways
coherent <- function(drv, cg, qgenes, prm) {
  if (!length(drv) || !length(cg)) return(FALSE)
  nin <- sum(drv %in% cg); base <- mean(qgenes %in% cg)
  nin >= prm$coh_min && base > 0 && (nin / length(drv)) / base >= prm$coh_ratio
}

#' Rank a library for one query: G, P, tiers and connecting pathways
#'
#' @param gq,pq lists with elements \code{up} and \code{dn}: query genes and query pathways (ordered).
#' @param LG library, genes x compounds.
#' @param LP library pathway scores, pathways x compounds.
#' @param sets gene sets named as \code{rownames(LP)}, on the genes of \code{LG}.
#' @param prm list: \code{n_gene}, \code{cut}, \code{n_top}, \code{spec_max}, \code{coh_min}, \code{coh_ratio}.
#' @param direction \code{"similar"} (mimic the query) or \code{"reverse"} (reverse it).
#' @param split split tier 1 into 1a (coherent) and 1b.
#' @return data.frame, one row per compound; attributes \code{conn} (connecting pathways) and \code{freq} (library frequency).
#' @export
rank_list <- function(gq, pq, LG, LP, sets, prm, direction = "similar", split = TRUE) {
  LG[!is.finite(LG)] <- 0; sgn <- if (direction == "similar") 1 else -1
  sG <- sgn * wtcs_query(prep_lib(LG), gq$up, gq$dn); sP <- sgn * wtcs_query(prep_lib(LP), pq$up, pq$dn)
  ok <- colnames(LG)[is.finite(sG) & is.finite(sP)]; rG <- nrk(sG[ok]); rP <- nrk(sP[ok]); pwn <- rownames(LP)
  U <- top_ind(LP, prm$n_top, TRUE); D <- top_ind(LP, prm$n_top, FALSE); fu <- rowMeans(U); fd <- rowMeans(D)
  if (direction == "similar") { A <- U; B <- D; fA <- fu; fB <- fd } else { A <- D; B <- U; fA <- fd; fB <- fu }
  qtU <- pwn %in% utils::head(pq$up, prm$n_top); qtD <- pwn %in% utils::head(pq$dn, prm$n_top)
  conn <- lapply(stats::setNames(ok, ok), function(j) list(up = pwn[qtU & A[, j]], dn = pwn[qtD & B[, j]]))
  lab <- vapply(conn, function(cj) { f <- c(fA[cj$up], fB[cj$dn])
    l <- c(if (length(cj$up)) paste0("query up: ", short(cj$up)), if (length(cj$dn)) paste0("query down: ", short(cj$dn)))
    o <- order(f); sp <- f[o] < prm$spec_max
    paste(sprintf("%s (%.0f%%)", l[o][sp], 100 * f[o][sp]), collapse = "; ") }, "")
  n_spec <- vapply(conn, function(cj) sum(fA[cj$up] < prm$spec_max) + sum(fB[cj$dn] < prm$spec_max), numeric(1))
  hiG <- rG <= prm$cut; hiP <- rP <= prm$cut
  tier <- ifelse(hiG & hiP, "1", ifelse(hiP, ifelse(n_spec > 0, "2", "2x"), ifelse(hiG, "3", "-")))
  if (split) {
    KIDX <- lapply(sets[pwn], function(g) { i <- match(g, rownames(LG)); i[!is.na(i)] })
    TGU <- apply(LG[, ok, drop = FALSE], 2, function(v) order(v, decreasing = TRUE)[seq_len(prm$n_gene)])
    TGD <- apply(LG[, ok, drop = FALSE], 2, function(v) order(v)[seq_len(prm$n_gene)])
    qUi <- match(gq$up, rownames(LG)); qDi <- match(gq$dn, rownames(LG)); qUi <- qUi[!is.na(qUi)]; qDi <- qDi[!is.na(qDi)]
    for (i in which(tier == "1")) { j <- ok[i]
      drv <- if (direction == "similar") c(intersect(qUi, TGU[, j]), intersect(qDi, TGD[, j])) else
                                         c(intersect(qUi, TGD[, j]), intersect(qDi, TGU[, j]))
      cg <- unique(unlist(KIDX[c(conn[[j]]$up, conn[[j]]$dn)]))
      tier[i] <- if (coherent(drv, cg, c(qUi, qDi), prm)) "1a" else "1b" }
  }
  out <- data.frame(drug = ok, G_rank = rG, P_rank = rP, WTCS_G = sgn * sG[ok], WTCS_P = sgn * sP[ok], tier = tier,
                    n_spec = n_spec, specific_pathways = lab, stringsAsFactors = FALSE)
  attr(out, "conn") <- conn; attr(out, "freq") <- list(A = fA, B = fB)
  out
}
