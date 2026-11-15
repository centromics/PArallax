#' Same-compound benchmark across cell lines
#'
#' \code{make_plan}: for each cell-line pair (A, B), queries = compounds measured in both
#' (up to \code{qmax}, fixed seed), library = the queries plus up to \code{extra} other compounds of B.
#' \code{eval_pair}: G, P, tiers, confirmation counts, explanation share and A/B overlap for every
#' query of one pair, in every pathway space.
#' @param IDX list per cell line of named vectors (names = compound IDs, original order).
#' @param pairs,full_pairs two-column matrices of cell-line pairs (the seed index follows \code{full_pairs}).
#' @param qmax,extra,seed plan size and seed.
#' @param pl one element of a plan.
#' @param gm function(cell, ids) returning genes x compounds.
#' @param pw function(space, cell, ids, weight_cell) returning pathways x compounds.
#' @param spaces list of pathway spaces (\code{name}, \code{ks}, \code{sets}, \code{tiers}).
#' @param prm list of parameters (see \code{rank_list}; also \code{ex_top}, \code{small_max}).
#' @return \code{make_plan}: list of pairs; \code{eval_pair}: data.frame, one row per query and space.
#' @export
make_plan <- function(IDX, pairs, full_pairs, qmax, extra, seed) lapply(seq_len(nrow(pairs)), function(i) {
  A <- pairs[i, 1]; B <- pairs[i, 2]; ii <- which(full_pairs[, 1] == A & full_pairs[, 2] == B)
  common <- sort(intersect(names(IDX[[A]]), names(IDX[[B]])))
  set.seed(seed + ii)
  q <- if (length(common) > qmax) sort(sample(common, qmax)) else common
  oth <- setdiff(names(IDX[[B]]), q); oth <- oth[sample.int(length(oth), min(extra, length(oth)))]
  list(A = A, B = B, query = q, lib = sort(c(q, oth)))
})
TIERS <- c("1a", "1b", "2", "2x", "3", "-")
#' @rdname make_plan
#' @export
eval_pair <- function(pl, gm, pw, spaces, prm) {
  LG <- gm(pl$B, pl$lib); LG[!is.finite(LG)] <- 0
  libG <- prep_lib(LG); ids <- colnames(LG)
  TGU <- apply(LG, 2, function(v) order(v, decreasing = TRUE)[seq_len(prm$n_gene)])
  TGD <- apply(LG, 2, function(v) order(v)[seq_len(prm$n_gene)])
  GA <- gm(pl$A, pl$query); GBq <- gm(pl$B, pl$query)
  CX <- lapply(spaces, function(s) {
    PB <- pw(s, pl$B, pl$lib, pl$B); PA <- pw(s, pl$A, pl$query, pl$B)[rownames(PB), , drop = FALSE]
    U <- top_ind(PB, prm$n_top, TRUE); D <- top_ind(PB, prm$n_top, FALSE)
    list(PA = PA, PB = PB, lib = prep_lib(PB), U = U, D = D, fu = rowMeans(U), fd = rowMeans(D),
         small = lengths(s$sets[rownames(PB)]) <= prm$small_max,
         KIDX = lapply(s$sets[rownames(PB)], function(g) { i <- match(g, rownames(LG)); i[!is.na(i)] }))
  })
  blank <- data.frame(space = NA_character_, k = NA_integer_, A = pl$A, B = pl$B, drug = NA_character_, lib = NA_integer_,
                      rank_G = NA_real_, rank_P = NA_real_, gene_ovl = NA_real_, pw_ovl = NA_real_,
                      n_list10 = NA_integer_, n_sup10 = NA_integer_, X_in10 = NA, X_sup10 = NA,
                      n_list05 = NA_integer_, n_sup05 = NA_integer_, X_in05 = NA, X_sup05 = NA,
                      spec_share = NA_real_, spec_n = NA_real_, small_n = NA_real_, X_tier = NA_character_,
                      n_1a = NA_real_, n_1b = NA_real_, n_2 = NA_real_, n_2x = NA_real_, n_3 = NA_real_, n_rest = NA_real_,
                      stringsAsFactors = FALSE)
  do.call(rbind, lapply(pl$query, function(d) {
    gq <- q_top(GA[, d], prm$n_gene); gb <- q_top(GBq[, d], prm$n_gene)
    sG <- wtcs_query(libG, gq$up, gq$dn); govl <- mean(c(jac(gq$up, gb$up), jac(gq$dn, gb$dn)))
    qUi <- match(gq$up, rownames(LG)); qDi <- match(gq$dn, rownames(LG)); qUi <- qUi[!is.na(qUi)]; qDi <- qDi[!is.na(qDi)]
    do.call(rbind, lapply(seq_along(spaces), function(si) { s <- spaces[[si]]; cc <- CX[[si]]; pwn <- rownames(cc$PB)
      do.call(rbind, lapply(seq_along(s$ks), function(ki) { k <- s$ks[ki]
        pq <- q_pw(cc$PA[, d], k); if (!length(pq$up) || !length(pq$dn)) return(NULL)
        sP <- wtcs_query(cc$lib, pq$up, pq$dn)
        ok <- ids[is.finite(sG) & is.finite(sP)]; if (!d %in% ok) return(NULL)
        rG <- nrk(sG[ok]); rP <- nrk(sP[ok]); x <- ok == d; hiP <- rP <= prm$cut
        ravg <- function(sc) (rank(-sc, ties.method = "average")[x] - 1) / (length(sc) - 1)
        y <- q_pw(cc$PB[, d], k)
        r <- blank; r$space <- s$name; r$k <- k; r$drug <- d; r$lib <- length(ok)
        r$rank_G <- ravg(sG[ok]); r$rank_P <- ravg(sP[ok]); r$gene_ovl <- govl
        r$pw_ovl <- mean(c(jac(pq$up, y$up), jac(pq$dn, y$dn)))
        for (cut in c(0.10, 0.05)) { tg <- sprintf("%02d", round(100 * cut)); L <- rG <= cut
          r[[paste0("n_list", tg)]] <- sum(L); r[[paste0("n_sup", tg)]] <- sum(L & hiP)
          r[[paste0("X_in", tg)]] <- L[x]; r[[paste0("X_sup", tg)]] <- L[x] && hiP[x] }
        if (ki == 1L) {
          qtU <- pwn %in% utils::head(pq$up, prm$n_top); qtD <- pwn %in% utils::head(pq$dn, prm$n_top)
          sA <- qtU & cc$fu < prm$spec_max; sB <- qtD & cc$fd < prm$spec_max
          cs <- function(ind, sel) if (!any(sel)) rep(0, length(ok)) else colSums(ind[sel, ok, drop = FALSE])
          n_spec <- cs(cc$U, sA) + cs(cc$D, sB); n_small <- cs(cc$U, sA & cc$small) + cs(cc$D, sB & cc$small)
          top <- order(rG)[seq_len(min(prm$ex_top, length(ok)))]
          r$spec_share <- mean(n_spec[top] > 0); r$spec_n <- sum(n_spec[top]); r$small_n <- sum(n_small[top])
          if (isTRUE(s$tiers)) {
            hiG <- rG <= prm$cut
            tier <- ifelse(hiG & hiP, "1", ifelse(hiP, ifelse(n_spec > 0, "2", "2x"), ifelse(hiG, "3", "-")))
            for (i in which(tier == "1")) { j <- ok[i]
              drv <- c(intersect(qUi, TGU[, j]), intersect(qDi, TGD[, j]))
              cg <- unique(unlist(cc$KIDX[c(pwn[qtU & cc$U[, j]], pwn[qtD & cc$D[, j]])]))
              tier[i] <- if (coherent(drv, cg, c(qUi, qDi), prm)) "1a" else "1b" }
            r$X_tier <- tier[x]; r[paste0("n_", c("1a", "1b", "2", "2x", "3", "rest"))] <- as.list(as.numeric(table(factor(tier, levels = TIERS))))
          }
        }
        r
      }))
    }))
  }))
}
