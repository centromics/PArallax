## Submission

First submission of PAnoply 1.1.0.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for CRAN

* No data are shipped. KEGG pathways are downloaded from the KEGG REST API only when the user asks for them
  (default database of `panoply_query()`); examples and tests use simulated data (`panoply_toy()`) and need no
  internet access. A failed download stops with an informative message.
* Suggested Bioconductor packages (GSVA, AnnotationDbi, org.Hs.eg.db, GSEABase, signatureSearch) are used
  conditionally via `requireNamespace()`.
