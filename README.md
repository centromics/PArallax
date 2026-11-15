# PArallax 1.1.0

**PArallax: connectivity mapping with depth from gene and pathway perspectives.**

Connectivity scoring in pathway space on both sides. The same statistic (WTCS) scores every library compound twice:
on genes (G) and on pathways (P), after the query and the library are both mapped to pathway scores.
G stays the base ranking; P adds evidence. Each compound gets a tier and the pathways that connect it to the query.

| tier | meaning |
|---|---|
| 1a | G and P both in the top `cut` (default 10%) of the library, and coherent |
| 1b | G and P both in the top, not coherent |
| 1 | G and P both in the top (`split = "none"`) |
| 3 | G in the top, P not |
| 2 / 2x | P in the top, G not (with / without an uncommon connecting pathway) |

Coherence (`split`): `"explore"` (default) = the change correlation inside the connecting pathways goes beyond random
genes (empirical P < 0.05), as in the paper's senescence case; `"rule"` = the genes behind G lie inside the connecting
pathways, as in the paper's same-compound benchmark; `"none"`.

## Install

```r
# from CRAN (once accepted)
install.packages("PArallax")

# from GitHub (development version)
install.packages("remotes")
remotes::install_github("OWNER/PArallax")
# while the repository is private: remotes::install_github("OWNER/PArallax", auth_token = "<your GitHub token>")

# optional packages
install.packages(c("msigdbr"))                                        # MSigDB databases
BiocManager::install(c("GSVA", "AnnotationDbi", "org.Hs.eg.db", "signatureSearch"))   # gsva mapping, ID conversion, LINCS compound information
```

## Use

```r
library(PArallax)
# expr: genes x samples (control and case); lincs: genes x compounds; same gene IDs (Entrez)
grp <- factor(c("ctrl", "ctrl", "ctrl", "case", "case", "case"), levels = c("ctrl", "case"))
ann <- signatureSearch::lincs_pert_info2             # compound name, MOA, target (or your own table)
res <- parallax_query(expr, grp, lincs, annotation = ann)   # default: latest KEGG, sgses_max, reverse, rank product
res                                                  # tiers and the ranked list
plot_query(res, "list")                              # the ranked list (format of Figures 3c/3d)
plot_query(res, "figure", query_label = "senescence")   # panels a-e on one page (format of Figure 7)
parallax_figures(res, "my_query", query_label = "senescence")   # both, as PDF + PNG at 183 mm
```

Main arguments of `parallax_query()`:

* **`db`, pathway database** (one at a time; several names are scored one by one, never merged, and returned as a list)
  * default: the latest KEGG, human, disease pathways (hsa05...) excluded, downloaded from the KEGG REST API once per session,
    named as in the paper (`"GLYCOLYSIS / GLUCONEOGENESIS [KEGG: hsa00010]"`); `"KEGG_ALL"` keeps the disease pathways
  * MSigDB C2 canonical pathways, fetched with msigdbr: `"REACTOME"`, `"WIKIPATHWAYS"`, `"PID"`, `"BIOCARTA"`, `"KEGG_MEDICUS"`, `"KEGG_LEGACY"`
  * your own gene sets: `pathways = list(A = setsA, B = setsB), db = "A"`, or a flat list without `db`
* **`mapping`**: `"sgses_max"` (default), `"avgz"`, `"zscore"`, `"ssgsea"`, `"gsva"`, `"singscore"`, `"ora_signed"`.
  The connectivity score is WTCS whatever the mapping.
* **`annotation`**: your own table, LINCS pert_info, or `signatureSearch::lincs_pert_info2`; the ID, name, MOA and target
  columns are found by name, and library columns are matched by full ID, the part before `"__"`, or the BRD compound ID.
* **Large library**: download and map once, then reuse:
  `kegg <- get_pathways(univ = rownames(lincs)); lp <- map_pathways(lincs, kegg, "sgses_max")`,
  then `parallax_query(..., pathways = kegg, library_pw = lp)`.

## Data

No pathway or compound data are shipped with the package.

| data | source | terms |
|---|---|---|
| KEGG (default) | KEGG REST API (rest.kegg.jp), at run time | free for academic use; see https://www.kegg.jp/kegg/legal.html |
| Reactome, WikiPathways, PID | MSigDB via msigdbr, at run time | MSigDB: CC BY 4.0 |
| KEGG_MEDICUS / KEGG_LEGACY / BioCarta | MSigDB via msigdbr, at run time | additional restrictions; see the MSigDB license terms |
| compound information | `signatureSearch::lincs_pert_info2` or your own table | LINCS / CMap terms |
| toy data | `parallax_toy()` generates it | — |

## Toy example

```r
toy <- parallax_toy()                                 # 600 genes, 8 samples, 150 compounds, two invented collections
res <- parallax_query(toy$expr, toy$group, toy$library, pathways = toy$pathways, db = "KEGG",
                     annotation = toy$annotation, signature = "meandiff",
                     n_gene = 30, k_pw = 6, n_top = 5, n_draw = 200, explore_min = 5)
res
plot_query(res, "figure", query_label = "toy case")
source(system.file("examples", "toy_example.R", package = "PArallax"))   # the full walk-through
```

## Notes

* Gene IDs of `expr` and `library` must match; pathway gene sets are converted between symbols and Entrez IDs when needed
  (AnnotationDbi + org.Hs.eg.db).
* Pathway scores use the standard weighted score (WTCS) only.
* The tiers rank candidates; they are not a measure of efficacy.
* `make_plan()` and `eval_pair()` run the same-compound benchmark across cell lines used in the paper.

## Author

Sung Young Kim (palelamp@kku.ac.kr), School of Medicine, Konkuk University, South Korea.
