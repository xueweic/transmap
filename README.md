# transmap

R package for TransMap, a framework for mapping distal genetic effects on
individual genes and co-regulatory programs across cell types. transmap builds
program phenotypes for trans-QTL mapping at two scales. All functions take
in-memory matrices.

| Scale | Step 1: find programs | Step 2: program phenotypes |
|-------|-----------------------|----------------------------|
| Local | `r_hotspot(Y)`: gene modules (R implementation of Hotspot) | `transmap_local(Y, module, method = "WPC")`: WPC (default), PC1, aPC or AUCell |
| Global | `r_transcca(data_list)`: sparse programs from a transBoost-style `data_list` (Z + N, Beta, Beta + N + P/se, log2FC, individual-level G + E, or Sigma_GE; optional whitening and cis masking), built into W by `w_construction()` | `transmap_global(Y, program)`: loading-weighted expression |

```r
library(transmap)
data(transmap_example)

# Local programs
hs  <- r_hotspot(transmap_example$Y, min_gene = 5)
loc <- transmap_local(transmap_example$Y, hs, null_z = transmap_example$null_z)

# Global programs
prog <- r_transcca(list(Z = transmap_example$z, N = transmap_example$n), K = 3)
glob <- transmap_global(transmap_example$Y, prog)

# Null distribution and empirical p-values of the global programs
null <- transcca_null(n_null = 20, n = transmap_example$n,
                      p = nrow(transmap_example$z), maf = runif(1000, 0.05, 0.5),
                      sigma_ee = cor(transmap_example$Y), K = 3,
                      sumabs = prog$sparsity)
transcca_pvalue(prog, null)
transcca_plot(prog, null)
```

See the vignettes `transmap_local` and `transmap_global` for details.

## Installation

```r
remotes::install_github("xueweic/transmap")
```

AUCell (Bioconductor) is needed only for `method = "AUCell"`.

## Local development

From the repository root:

```sh
R CMD INSTALL .
R CMD build .
R CMD check --no-manual transmap_*.tar.gz
```

Regenerate documentation with `roxygen2::roxygenise()` and the example data
with `Rscript data-raw/transmap_example.R`.

## Repository layout

- `R/`: package source and roxygen2 documentation.
- `man/`: generated help files.
- `data/`, `data-raw/`: simulated example data and the script that makes it.
- `tests/testthat/`: package tests.
- `vignettes/`: local and global program workflows.
- `.github/workflows/`: package checks, site builds, and an opt-in Pages workflow.

The package website is configured but not published. The Pages workflow requires
an explicit manual opt-in and GitHub Pages configuration. Review all content
before publishing: a Pages site can be public even when its repository is
private.

The related figure and analysis repository is
[`transmap_manuscript`](https://github.com/xueweic/transmap_manuscript).
