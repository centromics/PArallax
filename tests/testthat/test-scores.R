test_that("fast WTCS reproduces the reference WTCS", {
  set.seed(1)
  L <- matrix(rnorm(2000), 200, 10, dimnames = list(paste0("g", 1:200), paste0("c", 1:10)))
  up <- paste0("g", 1:20); dn <- paste0("g", 21:40)
  expect_equal(unname(wtcs_query(prep_lib(L), up, dn)), unname(PAnoply:::ref_wtcs(L, up, dn)), tolerance = 1e-10)
})

test_that("sgses_max is the signed maximum deviation of the weighted KS running sum", {
  set.seed(2)
  L <- matrix(rnorm(600), 200, 3, dimnames = list(paste0("g", 1:200), paste0("c", 1:3)))
  S <- paste0("g", seq(1, 200, by = 7))
  v <- sort(stats::setNames(L[, 2], rownames(L)), decreasing = TRUE); hit <- names(v) %in% S
  p <- numeric(length(v)); p[hit] <- abs(v[hit]) / sum(abs(v[hit]))
  dd <- cumsum(p) - cumsum(!hit) / (length(v) - sum(hit))
  expect_equal(unname(map_pathways(L, list(S = S))["S", 2]), unname(dd[which.max(abs(dd))]), tolerance = 1e-10)
})
