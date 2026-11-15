## PArallax toy example: run with  source(system.file("examples", "toy_example.R", package = "PArallax"))
## Takes a few seconds. Writes into the working folder: toy_KEGG_list / toy_KEGG_figure (.pdf, .png), the same for the
## two-database run (toy_singscore_*), and toy_REACTOME_landscape.pdf. Open the PDFs to see the figures.
library(PArallax)
toy <- parallax_toy()
str(toy, max.level = 1)
print(toy$truth)                                    # the seven built-in compounds

## 1. one database, default mapping (sgses_max), tier 1 split by the exploratory coherence
res <- parallax_query(toy$expr, toy$group, toy$library,
                     pathways = toy$pathways, db = "KEGG",            # db: which collection of `pathways`
                     mapping = "sgses_max",                           # pathway mapping of query and library
                     annotation = toy$annotation, direction = "reverse",
                     signature = "meandiff",                          # "rankprod" (default) is slower
                     split = "explore",
                     n_gene = 30, k_pw = 6, n_top = 5, n_draw = 200, explore_min = 5)
print(res)                                          # tiers and the ranked list
## the list (format of Figures 3c / 3d) and the one-page figure (format of Figure 7), written at their printed size
## (183 mm wide); drawing them in a small plot window can fail with "Viewport has zero dimension(s)"
parallax_figures(res, "toy_KEGG", query_label = "toy case", highlight = "revertin")

## 2. another mapping and two databases: each database is scored on its own, never merged
res2 <- parallax_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = c("KEGG", "REACTOME"),
                      mapping = "singscore", annotation = toy$annotation, signature = "meandiff",
                      n_gene = 30, k_pw = 6, n_top = 5, n_draw = 200, explore_min = 5)
print(res2)
ggplot2::ggsave("toy_REACTOME_landscape.pdf", plot_query(res2[["REACTOME"]], "landscape", query_label = "toy case"),
                width = 120, height = 110, units = "mm")
parallax_figures(res2, "toy_singscore", query_label = "toy case")

## 3. real data (not run here): Entrez-ID query and library, e.g. LINCS signatures
## res3 <- parallax_query(expr, group, lincs)                          # default: latest KEGG from rest.kegg.jp, no disease pathways
## res4 <- parallax_query(expr, group, lincs, db = "REACTOME")         # MSigDB Reactome, fetched with msigdbr
## ann  <- signatureSearch::lincs_pert_info2                          # compound names / MOA (install signatureSearch yourself)
## kegg <- get_pathways(univ = rownames(lincs))                       # latest KEGG once, then map a large library once
## lp   <- map_pathways(lincs, kegg, "sgses_max")
## res5 <- parallax_query(expr, group, lincs, pathways = kegg, library_pw = lp, annotation = ann)
