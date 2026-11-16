## Submission

First submission of PArallax 1.1.0.

## Test environments

* GitHub Actions: ubuntu-latest (R release and R devel), macos-latest (R release), windows-latest (R release),
  and ubuntu-latest (R release) without the suggested packages.
* win-builder (R devel).

## R CMD check results

<!-- fill in after the checks above have been run -->

## Notes for CRAN

* No data are shipped. KEGG pathways are downloaded from the KEGG REST API only when the user asks for them
  (default database of `parallax_query()`); examples and tests use simulated data (`parallax_toy()`) and need no
  internet access. A failed download stops with an informative message.
* Suggested packages (GSVA, AnnotationDbi, org.Hs.eg.db, GSEABase, msigdbr, signatureSearch) are used
  conditionally via `requireNamespace()`; examples and tests do not need them.
* Functions with a `seed` argument restore the user's random-number state on exit.
* Examples that write files write to `tempdir()` and remove the files.
