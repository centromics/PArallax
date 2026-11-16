## Submission

First submission of PArallax 1.1.0.

## Test environments

* GitHub Actions: ubuntu-latest (R release and R devel), windows-latest (R release),
  and ubuntu-latest (R release) without the suggested packages: 0 errors | 0 warnings | 0 notes.
* win-builder (R devel, 2026-11-14 r90890 ucrt).

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  Possibly misspelled words in DESCRIPTION: WTCS (the weighted connectivity score, spelled out where it
  first appears), Subramanian, et al (author names in the method references).

## Notes for CRAN

* No data are shipped. KEGG pathways are downloaded from the KEGG REST API only when the user asks for them
  (default database of `parallax_query()`); examples and tests use simulated data (`parallax_toy()`) and need no
  internet access. A failed download stops with an informative message.
* Suggested packages (GSVA, AnnotationDbi, org.Hs.eg.db, GSEABase, msigdbr, signatureSearch) are used
  conditionally via `requireNamespace()`; examples and tests do not need them.
* Functions with a `seed` argument restore the user's random-number state on exit.
* Examples that write files write to `tempdir()` and remove the files.
