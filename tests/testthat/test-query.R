sim <- function() {
  set.seed(1)
  genes <- paste0("g", 1:400)
  sets <- stats::setNames(lapply(1:40, function(i) sample(genes, 20)), paste0("PW", 1:40))
  eff <- stats::setNames(rep(0, 400), genes); eff[sets$PW1] <- 2; eff[sets$PW2] <- -2
  expr <- cbind(matrix(rnorm(1200), 400), matrix(rnorm(1200), 400) + eff); rownames(expr) <- genes
  lib <- matrix(rnorm(400 * 60), 400, dimnames = list(genes, paste0("cpd", 1:60)))
  lib[, "cpd1"] <- -eff + rnorm(400, sd = 0.5)
  list(expr = expr, lib = lib, sets = sets)
}

test_that("panoply_query puts the reversing compound in tier 1", {
  d <- sim()
  res <- panoply_query(d$expr, rep(0:1, each = 3), d$lib, d$sets, signature = "meandiff", split = "none",
                       n_gene = 30, k_pw = 5, explore = FALSE)
  expect_s3_class(res, "panoply_result")
  r <- res$all[res$all$drug == "cpd1", ]
  expect_lte(r$G_rank, 0.1); expect_lte(r$P_rank, 0.1); expect_equal(r$tier, "1")
  expect_equal(res$list$drug[1], "cpd1")
  expect_true(all(c("priority", "tier", "G_rank", "P_rank", "specific_pathways") %in% names(res$list)))
  expect_equal(res$settings$mapping, "sgses_max"); expect_equal(res$settings$db, "custom")
})

test_that("every mapping method runs and the score stays WTCS", {
  d <- sim()
  for (m in c("sgses_max", "avgz", "zscore", "ssgsea", "singscore", "ora_signed")) {
    res <- panoply_query(d$expr, rep(0:1, each = 3), d$lib, d$sets, mapping = m, signature = "meandiff", split = "none",
                         n_gene = 30, k_pw = 5, explore = FALSE)
    expect_equal(res$settings$mapping, m); expect_equal(res$settings$score, "WTCS")
  }
})

test_that("split modes, explore and every figure run", {
  d <- sim()
  for (sp in c("explore", "rule")) {
    res <- panoply_query(d$expr, rep(0:1, each = 3), d$lib, d$sets, signature = "meandiff", n_gene = 30, k_pw = 5,
                         split = sp, n_draw = 20, explore_min = 5)
    expect_true(all(res$list$tier %in% c("1a", "1b", "3")))
  }
  for (t in c("list", "landscape", "signature", "cases", "pathways")) expect_s3_class(plot_query(res, t), "ggplot")
  expect_s3_class(plot_query(res, "figure"), "panoply_figure")
})

test_that("db picks one collection; several db are scored separately", {
  toy <- panoply_toy()
  k <- get_pathways(toy$pathways, "kegg", rownames(toy$library))
  expect_equal(attr(k, "db"), "KEGG"); expect_true(all(grepl("^KEGG_", names(k))))
  expect_error(get_pathways(toy$pathways, NULL, rownames(toy$library)))
  m <- panoply_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = c("KEGG", "REACTOME"),
                     signature = "meandiff", n_gene = 30, k_pw = 6, n_top = 5, n_draw = 50, explore_min = 5)
  expect_s3_class(m, "panoply_multi"); expect_equal(names(m), c("KEGG", "REACTOME"))
  expect_true(all(grepl("^REACTOME_", rownames(m$REACTOME$library_pw))))
})

test_that("the toy reverser is found and a mismatched library_pw is refused", {
  toy <- panoply_toy()
  res <- panoply_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = "KEGG", annotation = toy$annotation,
                       signature = "meandiff", n_gene = 30, k_pw = 6, n_top = 5, n_draw = 50, explore_min = 5)
  expect_true(res$all$tier[res$all$drug == "TOY001"] %in% c("1a", "1b"))
  lp <- map_pathways(toy$library, get_pathways(toy$pathways, "KEGG", rownames(toy$library)), "avgz")
  expect_error(panoply_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = "KEGG", library_pw = lp,
                             signature = "meandiff", n_gene = 30, k_pw = 6, n_top = 5, explore = FALSE, split = "none"))
})

test_that("bad group coding is refused", {
  d <- sim()
  expect_error(panoply_query(d$expr, rep(0, 6), d$lib, d$sets, signature = "meandiff"))
})

test_that("annotation columns are found by name and BRD IDs are matched", {
  R <- data.frame(drug = c("CPC006_A549_6H:BRD-K12345678-001-01-1:10", "TOY002__x", "cpd3"), stringsAsFactors = FALSE)
  ann <- data.frame(pert_id = c("BRD-K12345678", "TOY002"), pert_iname = c("drugA", "drugB"), moa = c("m1", "m2"),
                    target = c("T1", "T2"), stringsAsFactors = FALSE)
  out <- PAnoply:::annotate_compounds(R, ann)
  expect_equal(out$name, c("drugA", "drugB", "cpd3")); expect_equal(out$target, c("T1", "T2", NA))
})

test_that("the built-in rank product finds the raised and lowered genes", {
  d <- sim()
  s <- make_signature(d$expr, rep(0:1, each = 3), 20, "rankprod", nperm = 20, seed = 1)
  expect_gt(mean(s$up %in% d$sets$PW1), 0.5); expect_gt(mean(s$dn %in% d$sets$PW2), 0.5)
  expect_true(all(c("feature", "AveFC", "FDR.up", "FDR.down") %in% names(s$table)))
})
