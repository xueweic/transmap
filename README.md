# transmap

Development R package for TransMap, a framework for mapping distal genetic
effects on individual genes and co-regulatory programs across cell types.

This repository is an installable package scaffold. It does not yet contain the
validated TransMap inference functions or claim to reproduce manuscript results.
Those implementations, examples, and tests will be added as the analysis code
is reviewed and moved into the package.

## Local development

From the repository root:

```sh
R CMD INSTALL .
R CMD build .
R CMD check --no-manual transmap_*.tar.gz
```

After the repository has code and has been pushed, a private GitHub installation
can use `remotes::install_github("xueweic/transmap")` with appropriate access.

## Repository layout

- `R/`: package source and roxygen2 documentation.
- `man/`: generated help files; update with `roxygen2::roxygenise()`.
- `tests/testthat/`: package tests.
- `vignettes/`: longer-form usage and methods documentation.
- `inst/`: future small, distributable example files only.
- `.github/workflows/`: package checks, site builds, and an opt-in Pages workflow.

The package website is configured but not published. The Pages workflow requires
an explicit manual opt-in and GitHub Pages configuration. Review all content
before publishing: a Pages site can be public even when its repository is
private.

The related figure and analysis repository is
[`transmap_manuscript`](https://github.com/xueweic/transmap_manuscript).
