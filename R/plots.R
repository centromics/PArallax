## ---- figure style: journal sizes (8 pt text, 7 pt axes) at the printed width (183 mm) -----------------------------
.pcol <- list(blue = "#2a78d6", lblue = "#9fc1ec", green = "#1baf7a", orange = "#eb6834", violet = "#6250d6", yellow = "#fdf3d6",
              ink = "#0b0b0b", ink2 = "#52514e", ink3 = "#898781", grid = "#e1e0d9", paper = "#f4f3ee", rest = "#d4d3cc", seg = "#c3c2b7")
.ptheme <- function() ggplot2::theme_minimal(base_size = 8) + ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 9), plot.subtitle = ggplot2::element_text(colour = .pcol$ink2, size = 7),
  plot.caption = ggplot2::element_text(colour = .pcol$ink3, hjust = 0, size = 6.5), panel.grid.minor = ggplot2::element_blank(),
  panel.grid.major = ggplot2::element_line(colour = .pcol$grid, linewidth = 0.3), legend.position = "bottom",
  legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 7), legend.key.size = grid::unit(3.5, "mm"),
  axis.text = ggplot2::element_text(colour = .pcol$ink2, size = 7), axis.title = ggplot2::element_text(size = 8),
  strip.text = ggplot2::element_text(face = "bold", hjust = 0, size = 7.5),
  plot.margin = ggplot2::margin(4, 6, 4, 4), legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "mm"))
.prepel <- function(data, mapping, size = 2.3, ...) ggrepel::geom_text_repel(data = data, mapping = mapping, size = size, seed = 1,
                                                                             min.segment.length = 0, box.padding = 0.3, max.overlaps = 50, ...)
.TL <- c(`1a` = "1a  G + P + coherent", `1b` = "1b  G + P, not coherent", `1` = "1  G + P", `3` = "3  G only", P = "P only (2, 2x)", rest = "rest")
.TC <- stats::setNames(c(.pcol$blue, .pcol$orange, .pcol$blue, .pcol$green, .pcol$violet, .pcol$rest), .TL)
.labs_dir <- function(x, ql) list(up = paste("up in", ql), dn = paste("down in", ql), rev = x$settings$direction == "reverse")
.symbols <- function(g) {
  g <- as.character(g)
  if (!length(g) || id_type(g) != "entrez" || !requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(g)
  s <- conv_map(g, "entrez", "symbol"); s <- unname(s[g]); ifelse(is.na(s), g, s)
}
## compounds shown in the case panels: 1a (in coherence order when available), else tier 1 (G order)
.case_rows <- function(x, n_case) {
  L <- x$list; A <- L[L$tier == "1a", , drop = FALSE]
  if (nrow(A)) { if (any(!is.na(A$coh_P))) A <- A[order(A$coh_P, A$G_rank), , drop = FALSE] } else A <- L[L$tier %in% c("1", "1b"), , drop = FALSE]
  utils::head(A, n_case)
}

## ---- panels -------------------------------------------------------------------------------------------------------------
.p_signature <- function(x, n, ql) {
  ld <- .labs_dir(x, ql); sp <- x$signature$pathways; up <- utils::head(sp$up, n); dn <- utils::head(sp$dn, n)
  D <- sp$table[match(c(up, dn), sp$table$feature), , drop = FALSE]; D$dir <- rep(c(ld$up, ld$dn), c(length(up), length(dn)))
  if ("FDR.up" %in% names(D)) {
    D$v <- -log10(pmax(ifelse(D$dir == ld$up, D$FDR.up, D$FDR.down), 1e-4)) * ifelse(D$dir == ld$up, 1, -1)
    xl <- "-log10 FDR (dashed: FDR 0.1)"; ref <- c(-1, 1)
  } else { D$v <- D$diff; xl <- paste0("mean pathway score difference (", ql, " - control)"); ref <- numeric(0) }
  D$lab <- factor(short(D$feature), levels = unique(short(D$feature)[order(D$v)]))
  ggplot2::ggplot(D, ggplot2::aes(v, lab, fill = dir)) + ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_vline(xintercept = ref, colour = .pcol$ink3, linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = stats::setNames(c(.pcol$orange, .pcol$blue), c(ld$up, ld$dn))) +
    ggplot2::scale_x_continuous(labels = function(v) abs(v)) +
    ggplot2::labs(title = paste("The", ql, "signature in pathway space"),
                  subtitle = sprintf("Top %d %s pathways up and down (%s mapping, %s)", n, x$settings$db, x$settings$mapping,
                                     if (x$settings$signature == "rankprod") "rank product" else "mean difference"),
                  x = xl, y = NULL) + .ptheme()
}
.p_landscape <- function(x, ql) {
  R <- x$all; cut <- x$settings$cut
  key <- ifelse(R$tier %in% c("2", "2x"), "P", ifelse(R$tier == "-", "rest", R$tier))
  lv <- intersect(c("1a", "1b", "1", "3", "P", "rest"), unique(key))
  o <- order(match(key, rev(lv))); R <- R[o, , drop = FALSE]; key <- key[o]
  R$grp <- factor(unname(.TL[key]), levels = unname(.TL[lv]))
  A <- if (any(key == "1a")) R[key == "1a", , drop = FALSE] else R[key %in% c("1", "1b"), , drop = FALSE]
  A <- utils::head(A[order(A$G_rank), , drop = FALSE], 15)
  lim <- min(1, max(0.35, 3.5 * cut))
  sz <- stats::setNames(c(2.2, 1.7, 2.0, 1.4, 1.4, 0.8), .TL)
  ggplot2::ggplot(R, ggplot2::aes(G_rank, P_rank)) +
    ggplot2::annotate("rect", xmin = 0, xmax = cut, ymin = 0, ymax = cut, fill = .pcol$blue, alpha = 0.06) +
    ggplot2::geom_hline(yintercept = cut, colour = .pcol$ink3, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = cut, colour = .pcol$ink3, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = grp, size = grp), alpha = 0.85, stroke = 0) +
    ggplot2::geom_point(data = A, shape = 21, size = 2.8, colour = .pcol$ink, fill = NA, stroke = 0.4) +
    .prepel(A, ggplot2::aes(label = name), size = 2.2, fontface = "bold", colour = .pcol$blue) +
    ggplot2::scale_colour_manual(values = .TC[levels(R$grp)], drop = FALSE) +
    ggplot2::scale_size_manual(values = sz[levels(R$grp)], guide = "none") +
    ggplot2::scale_x_continuous(labels = scales::percent) + ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::coord_cartesian(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, override.aes = list(size = 2.5))) +
    ggplot2::labs(title = "Every compound placed by two connectivity scores",
                  subtitle = sprintf("%d library compounds, %s the %s signature; shaded: tier 1 (G and P both top %.0f%%)",
                                     nrow(R), if (x$settings$direction == "reverse") "reversing" else "mimicking", ql, 100 * cut),
                  x = "G rank (gene-level WTCS)", y = sprintf("P rank (pathway-space WTCS, %s)", x$settings$db)) + .ptheme()
}
.p_coherence <- function(x, ql) {
  D <- x$list[x$list$tier %in% c("1a", "1b", "1") & !is.na(x$list$coh_rho), , drop = FALSE]
  if (!nrow(D)) return(NULL)
  D$name <- factor(D$name, levels = unique(D$name[order(-D$coh_P)]))
  D$grp <- factor(unname(.TL[D$tier]), levels = unname(.TL[intersect(c("1a", "1b", "1"), D$tier)]))
  D$plab <- ifelse(D$coh_P < 0.001, "P<0.001", sprintf("P=%.3f", D$coh_P))
  ggplot2::ggplot(D, ggplot2::aes(y = name)) + ggplot2::geom_vline(xintercept = 0, colour = .pcol$ink3, linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = coh_bg, xend = coh_rho, yend = name), colour = .pcol$seg, linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(x = coh_bg), colour = .pcol$ink3, size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(x = coh_rho, colour = grp), size = 2.6) +
    ggplot2::geom_text(ggplot2::aes(x = pmax(coh_rho, coh_bg), label = plab), hjust = -0.25, size = 2.2, colour = .pcol$ink2) +
    ggplot2::scale_colour_manual(values = .TC[levels(D$grp)], drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.25))) +
    ggplot2::labs(title = "Coherence of the tier-1 candidates",
                  subtitle = sprintf("Spearman correlation of %s change vs compound change: connecting-pathway genes\n(colour) vs random genes (grey, %s draws)",
                                     ql, format(x$settings$n_draw, big.mark = ",")),
                  x = sprintf("Spearman correlation (%s)", if (x$settings$direction == "reverse") "negative = reversal" else "positive = similarity"),
                  y = NULL) + .ptheme()
}
.p_cases <- function(x, n_case, ql) {
  ld <- .labs_dir(x, ql); A <- .case_rows(x, n_case); if (!nrow(A)) return(NULL)
  pu <- paste("pathway", ld$up); pd <- paste("pathway", ld$dn)
  D <- do.call(rbind, lapply(seq_len(nrow(A)), function(i) { j <- A$drug[i]; z <- x$change_lib[[j]]; if (is.null(z)) return(NULL)
    gu <- unlist(x$sets[x$conn[[j]]$up]); gd <- unlist(x$sets[x$conn[[j]]$dn])
    gs <- intersect(intersect(unique(c(gu, gd)), names(x$change)), names(z)); if (!length(gs)) return(NULL)
    lab <- if (!is.na(A$coh_rho[i])) sprintf("%s  (rho %.2f vs %.2f)", A$name[i], A$coh_rho[i], A$coh_bg[i]) else
      sprintf("%s  (G %.1f%%, P %.1f%%)", A$name[i], 100 * A$G_rank[i], 100 * A$P_rank[i])
    data.frame(case = lab, gene = gs, d = unname(x$change[gs]), z = unname(z[gs]),
               pdir = ifelse(gs %in% gu & gs %in% gd, "both", ifelse(gs %in% gu, pu, pd)), stringsAsFactors = FALSE) }))
  if (is.null(D)) return(NULL)
  D$case <- factor(D$case, levels = unique(D$case)); D$hit <- if (ld$rev) -D$d * D$z else D$d * D$z
  LAB <- do.call(rbind, lapply(split(D, D$case), function(e) { e <- e[e$hit > 0, , drop = FALSE]; utils::head(e[order(-e$hit), , drop = FALSE], 4) }))
  q <- if (ld$rev) list(c(0, Inf, -Inf, 0), c(-Inf, 0, 0, Inf)) else list(c(0, Inf, 0, Inf), c(-Inf, 0, -Inf, 0))
  p <- ggplot2::ggplot(D, ggplot2::aes(d, z)) +
    ggplot2::annotate("rect", xmin = q[[1]][1], xmax = q[[1]][2], ymin = q[[1]][3], ymax = q[[1]][4], fill = .pcol$blue, alpha = 0.05) +
    ggplot2::annotate("rect", xmin = q[[2]][1], xmax = q[[2]][2], ymin = q[[2]][3], ymax = q[[2]][4], fill = .pcol$blue, alpha = 0.05) +
    ggplot2::geom_hline(yintercept = 0, colour = .pcol$ink3, linewidth = 0.3) + ggplot2::geom_vline(xintercept = 0, colour = .pcol$ink3, linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = pdir), size = 1.2, alpha = 0.75, stroke = 0) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = .pcol$ink, linewidth = 0.4) +
    ggplot2::facet_wrap(~ case, scales = "free", ncol = 3) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(.pcol$orange, .pcol$blue, .pcol$ink3), c(pu, pd, "both"))) +
    ggplot2::labs(title = sprintf("%s cases: the %s at gene level, inside the connecting pathways",
                                  if (any(A$tier == "1a")) "1a" else "Tier-1", if (ld$rev) "reversal" else "similarity"),
                  subtitle = sprintf("Each point is a gene of the connecting pathways; shaded quadrants = %s; labelled: the four strongest genes",
                                     if (ld$rev) "reversed" else "same direction"),
                  x = sprintf("change in %s (%s - control)", ql, ql), y = "change with the compound") + .ptheme()
  if (!is.null(LAB) && nrow(LAB)) { LAB$sym <- .symbols(LAB$gene); p <- p + .prepel(LAB, ggplot2::aes(label = sym), size = 2.1, colour = .pcol$ink) }
  p
}
.p_pathways <- function(x, n_case, ql) {
  A <- .case_rows(x, n_case); if (!nrow(A)) return(NULL)
  cn <- function(j) c(x$conn[[j]]$up, x$conn[[j]]$dn); CP <- unique(unlist(lapply(A$drug, cn))); if (!length(CP)) return(NULL)
  qup <- utils::head(x$signature$pathways$up, x$settings$n_top); pw <- rownames(x$library_pw)
  D <- expand.grid(pathway = CP, drug = A$drug, stringsAsFactors = FALSE)
  D$score <- x$library_pw[cbind(match(D$pathway, pw), match(D$drug, colnames(x$library_pw)))]
  D$connecting <- mapply(function(a, j) a %in% cn(j), D$pathway, D$drug)
  nconn <- tapply(D$connecting, D$pathway, sum); keep <- names(sort(-nconn))[seq_len(min(25, length(nconn)))]
  D <- D[D$pathway %in% keep, , drop = FALSE]
  D$spec <- ifelse(D$pathway %in% qup, x$freq$A[D$pathway], x$freq$B[D$pathway])
  D$name <- factor(A$name[match(D$drug, A$drug)], levels = unique(A$name))
  D$plab <- paste0(ifelse(D$pathway %in% qup, "up  ", "down  "), short(D$pathway))
  D$plab <- factor(D$plab, levels = rev(unique(D$plab[order(!D$pathway %in% qup, -stats::ave(as.numeric(D$connecting), D$pathway, FUN = sum))])))
  lim <- max(abs(D$score), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
  ggplot2::ggplot(D, ggplot2::aes(name, plab)) + ggplot2::geom_tile(ggplot2::aes(fill = score), colour = "white", linewidth = 0.6) +
    ggplot2::geom_tile(data = D[D$connecting, , drop = FALSE], fill = NA, colour = .pcol$ink, linewidth = 0.5) +
    ggplot2::geom_text(data = D[D$connecting, , drop = FALSE], ggplot2::aes(label = sprintf("%.0f%%", 100 * spec)), size = 2.1, colour = .pcol$ink) +
    ggplot2::scale_fill_gradient2(low = .pcol$blue, mid = .pcol$paper, high = .pcol$orange, midpoint = 0, limits = c(-lim, lim),
                                  name = sprintf("compound's pathway score (%s)", x$settings$mapping)) +
    ggplot2::labs(title = sprintf("%s pathways each %s compound %s", paste0(toupper(substr(ql, 1, 1)), substring(ql, 2)), if (any(A$tier == "1a")) "1a" else "tier-1",
                                  if (x$settings$direction == "reverse") "reverses" else "shares"),
                  subtitle = sprintf("Outlined: connecting (%s direction, both top %d); number: share of the library\nwith this pathway in its top %d (uncommon if < %.0f%%)",
                                     if (x$settings$direction == "reverse") "opposite" else "same", x$settings$n_top, x$settings$n_top, 100 * x$settings$spec_max),
                  x = NULL, y = NULL) + .ptheme() +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(size = 7, colour = .pcol$ink2),
                   legend.key.width = grid::unit(8, "mm"), legend.key.height = grid::unit(2.5, "mm"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, face = "bold"), panel.grid = ggplot2::element_blank())
}
## the list a researcher receives (format of Figures 3c and 3d)
.p_list <- function(x, n, highlight) {
  L <- x$list; L$G_pos <- rank(L$G_rank, ties.method = "first"); L <- utils::head(L[order(L$priority), , drop = FALSE], n)
  if (!nrow(L)) stop("the list is empty (no compound in the top ", 100 * x$settings$cut, "% by G)")
  L$is_X <- if (is.null(highlight)) rep(FALSE, nrow(L)) else L$drug %in% highlight | L$name %in% highlight
  L$y <- rev(seq_len(nrow(L))); lv <- intersect(c("1a", "1b", "1", "3"), L$tier); L$tl <- factor(unname(.TL[L$tier]), levels = unname(.TL[lv]))
  sp <- ifelse(is.na(L$specific_pathways), "", L$specific_pathways)
  first <- sub("; .*$", "", sp); more <- ifelse(sp == "", 0, lengths(strsplit(sp, "; ")) - 1)
  L$spec <- paste(vapply(first, function(s) paste(strwrap(s, 52), collapse = "\n"), ""), ifelse(more > 0, paste0("+", more), ""))
  X0 <- c(prio = 0.4, tier = 1.0, name = 1.5, gpos = 4.7, grank = 5.5, prank = 6.3, spec = 7.0)
  H <- data.frame(x = X0, lab = c("#", "tier", "compound", "G position", "G rank", "P rank", "rarest uncommon connecting pathway (% of library)"),
                  hj = c(0.5, 0.5, 0, 0.5, 0.5, 0.5, 0))
  fillc <- stats::setNames(c(.pcol$blue, .pcol$lblue, .pcol$lblue, .pcol$green), .TL[c("1a", "1b", "1", "3")])
  hl <- if (any(L$is_X)) sprintf(" Highlighted: priority %s.", paste(L$priority[L$is_X], collapse = ", ")) else ""
  ggplot2::ggplot(L) +
    ggplot2::geom_rect(data = L[L$is_X, , drop = FALSE], ggplot2::aes(ymin = y - 0.5, ymax = y + 0.5), xmin = 0, xmax = 13, fill = .pcol$yellow) +
    ggplot2::geom_tile(ggplot2::aes(x = X0[["tier"]], y = y, fill = tl), width = 0.42, height = 0.72) +
    ggplot2::geom_text(ggplot2::aes(x = X0[["prio"]], y = y, label = priority), size = 2.4, colour = .pcol$ink2) +
    ggplot2::geom_text(ggplot2::aes(x = X0[["tier"]], y = y, label = tier), size = 2.3, fontface = "bold", colour = ifelse(L$tier == "1a", "white", .pcol$ink)) +
    ggplot2::geom_text(ggplot2::aes(x = X0[["name"]], y = y, label = name, fontface = ifelse(is_X, "bold", "plain")), hjust = 0, size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(x = X0[["gpos"]], y = y, label = G_pos, colour = G_pos != priority), size = 2.4) +
    ggplot2::geom_text(ggplot2::aes(x = X0[["grank"]], y = y, label = sprintf("%.1f%%", 100 * G_rank)), size = 2.4, colour = .pcol$ink2) +
    ggplot2::geom_text(ggplot2::aes(x = X0[["prank"]], y = y, label = sprintf("%.1f%%", 100 * P_rank)), size = 2.4, colour = .pcol$ink2) +
    ggplot2::geom_text(ggplot2::aes(x = X0[["spec"]], y = y, label = spec), hjust = 0, size = 2.15, lineheight = 0.85) +
    ggplot2::geom_text(data = H, ggplot2::aes(x = x, y = max(L$y) + 1, label = lab, hjust = hj), size = 2.3, fontface = "bold", colour = .pcol$ink2) +
    ggplot2::scale_fill_manual(values = fillc[levels(L$tl)], drop = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#c2410c", `FALSE` = .pcol$ink2), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 13), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.4, max(L$y) + 1.6), expand = c(0, 0)) +
    ggplot2::labs(title = sprintf("Ranked compounds: %s pathways, %s mapping", x$settings$db, x$settings$mapping),
                  subtitle = sprintf("First %d rows (1a, then 1b, then 3; G order inside each tier).%s", nrow(L), hl),
                  caption = "G position: place in the plain gene-level ranking (orange when the tier order moved it).", x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 9), plot.subtitle = ggplot2::element_text(colour = .pcol$ink2, size = 7),
                   plot.caption = ggplot2::element_text(colour = .pcol$ink3, hjust = 0, size = 6.5), legend.position = "bottom",
                   legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 7), legend.key.size = grid::unit(3.5, "mm"),
                   plot.margin = ggplot2::margin(6, 8, 6, 8), plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}
## all case panels on one page (format of Figure 7): a signature, b landscape, c coherence, d cases, e pathways
.p_figure <- function(x, n, n_case, ql) {
  ttl <- function(p, letter) if (is.null(p)) grid::nullGrob() else p + ggplot2::labs(title = paste0(letter, "  ", p$labels$title))
  P <- list(ttl(.p_signature(x, n, ql), "a"), ttl(.p_landscape(x, ql), "b"), ttl(.p_coherence(x, ql), "c"),
            ttl(.p_cases(x, n_case, ql), "d"), ttl(.p_pathways(x, n_case, ql), "e"))
  P <- lapply(P, function(p) if (inherits(p, "ggplot")) ggplot2::ggplotGrob(p) else p)
  g <- gridExtra::arrangeGrob(grobs = P, layout_matrix = rbind(c(1, 2), c(3, 5), c(4, 4)), heights = c(1, 1.05, 1.15), widths = c(1, 1.05))
  class(g) <- c("parallax_figure", class(g)); g
}

#' Figures of a parallax_result
#'
#' \code{plot_query} draws one figure: \code{"figure"} (all case panels on one page: a query signature in pathway space,
#' b every compound by G rank and P rank, c coherence of tier 1, d gene-level change of the 1a compounds inside their
#' connecting pathways, e connecting pathways of the 1a compounds; the format of Figure 7 of the paper),
#' \code{"list"} (the ranked list a researcher receives; the format of Figures 3c and 3d), or one panel
#' (\code{"signature"}, \code{"landscape"}, \code{"coherence"}, \code{"cases"}, \code{"pathways"}).
#' \code{parallax_figures} writes the list and the figure as PDF and PNG (600 dpi) at the printed width of 183 mm.
#' @param x a \code{parallax_result} (\code{parallax_figures} also takes a \code{parallax_multi}: one pair of files per database).
#' @param type figure type.
#' @param n rows of the list (default 25) or pathways per direction of the signature (default \code{n_top}).
#' @param highlight compound IDs or names to highlight in the list.
#' @param query_label word for the query in titles and axes (e.g. \code{"senescence"}).
#' @param n_case compounds in the case panels (d, e).
#' @param prefix file name prefix (may include a folder).
#' @param n_list rows of the list.
#' @return \code{plot_query}: a ggplot (a grob of class \code{parallax_figure} for \code{"figure"}; print it to draw);
#'   \code{parallax_figures}: the file names, invisibly.
#' @export
plot_query <- function(x, type = c("figure", "list", "landscape", "signature", "coherence", "cases", "pathways"),
                       n = NULL, highlight = NULL, query_label = "query", n_case = 6L) {
  type <- match.arg(type)
  if (inherits(x, "parallax_multi")) stop("x holds several databases: plot one of them, e.g. plot_query(x[[\"", names(x)[1], "\"]])")
  stopifnot(inherits(x, "parallax_result"))
  if (is.null(n)) n <- if (type == "list") 25L else x$settings$n_top
  need1 <- function(p, what) if (is.null(p)) stop("no ", what, " to draw") else p
  switch(type,
         figure    = .p_figure(x, n, n_case, query_label),
         list      = .p_list(x, n, highlight),
         landscape = .p_landscape(x, query_label),
         signature = .p_signature(x, n, query_label),
         coherence = need1(.p_coherence(x, query_label), "exploratory coherence (explore = FALSE, or too few genes)"),
         cases     = need1(.p_cases(x, n_case, query_label), "tier-1 compound with connecting-pathway genes"),
         pathways  = need1(.p_pathways(x, n_case, query_label), "connecting pathway of a tier-1 compound"))
}

#' @rdname plot_query
#' @export
parallax_figures <- function(x, prefix = "parallax", query_label = "query", highlight = NULL, n_list = 25L, n_case = 6L) {
  if (inherits(x, "parallax_multi"))
    return(invisible(unlist(lapply(names(x), function(d) parallax_figures(x[[d]], paste0(prefix, "_", gsub("[^A-Za-z0-9]", "", d)),
                                                                         query_label, highlight, n_list, n_case)))))
  stopifnot(inherits(x, "parallax_result"))
  mm <- function(v) v / 25.4
  pdf_dev <- function(f, w, h) if (capabilities("cairo")) grDevices::cairo_pdf(f, width = w, height = h) else grDevices::pdf(f, width = w, height = h)
  out <- character(0)
  draw <- function(g, name, w, h) {
    f <- paste0(prefix, "_", name, c(".pdf", ".png"))
    pdf_dev(f[1], w, h); grid::grid.newpage(); grid::grid.draw(g); grDevices::dev.off()
    grDevices::png(f[2], width = w, height = h, units = "in", res = 600, bg = "white"); grid::grid.newpage(); grid::grid.draw(g); grDevices::dev.off()
    out <<- c(out, f)
  }
  pl <- .p_list(x, n_list, highlight); nr <- min(n_list, nrow(x$list))
  draw(ggplot2::ggplotGrob(pl), "list", mm(183), mm(30 + 4.6 * nr))
  draw(.p_figure(x, x$settings$n_top, n_case, query_label), "figure", mm(183), mm(240))
  cat("written:", paste(out, collapse = ", "), "\n")
  invisible(out)
}

#' @export
print.parallax_figure <- function(x, ...) { grid::grid.newpage(); grid::grid.draw(x); invisible(x) }
