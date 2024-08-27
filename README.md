<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- First time: run usethis::use_readme_rmd() to create a pre-commit hook that 
prevents from committing if the README.Rmd has changed, but has not been 
re-knitted to generate an updated README.md -->

## bayesSYNC - A Bayesian functional factor model for high-dimensional molecular curves

<!-- Run for the R CMD checks, run usethis::use_github_actions() to set up the pipeline, possibly modify the .yaml file and then: -->
<!-- [![R build status](https://github.com/hruffieux/bayesSYNC/workflows/R-CMD-check/badge.svg)](https://github.com/hruffieux/bayesSYNC/actions) -->
<!-- [![](https://travis-ci.org/hruffieux/bayesSYNC.svg?branch=master)](https://travis-ci.org/hruffieux/bayesSYNC) -->

[![License: GPL
v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://img.shields.io/badge/devel%20version-0.1.0-blue.svg)](https://github.com/hruffieux/bayesSYNC)
<!-- [![](https://img.shields.io/github/languages/code-size/hruffieux/bayesSYNC.svg)](https://github.com/hruffieux/bayesSYNC) -->

Welcome to **bayesSYNC**, an R package that implements a Bayesian
functional factor model for analysing high-dimensional curves. This
package is particularly suited for handling longitudinal measurements of
variables, offering a tool for understanding complex temporal dynamics
in high-dimensional data.

### Overview

The increasing availability of longitudinal data allows for the study of
temporal patterns across numerous variables. However, modelling
coordinated temporal variations in high-dimensional settings remains
challenging. **bayesSYNC** addresses this challenge by introducing a
Bayesian framework that combines latent factor modelling with functional
principal component analysis (FPCA).

This approach allows for: - capturing correlations across both variables
and time; - representing high-dimensional curves using a small number of
FPCA expansions, which reflect underlying latent factors; - modelling
individual variability through functional principal components, each
characterised by smoothly varying temporal functions.

For example, in a biological context, the variables might represent gene
expression levels, and the latent factors could correspond to biological
pathways. **bayesSYNC** can help elucidate how groups of genes work
together over time, providing insights into the biological processes
underlying diseases.

**Key features:** - *Variable-specific loadings*: estimation of loadings
for each variable that contributes to the FPCA expansions, providing
insights into how subsets of variables influence latent factors. -
*Subject-specific component scores*: estimation of subject-specific
scores for each functional principal component, allowing for
personalised analyses of temporal dynamics. - *Variational inference*:
variational inference algorithm with analytical updates, ensuring both
computational efficiency and robust uncertainty quantification. -
*Scalability*: the method is capable of handling realistic data sizes
(e.g., longitudinal measurements for each of ~20 000 genes in hundreds
of individuals), making it practical for large-scale studies.

### Installation

Then, to install the package in R, run the following command:

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("hruffieux/bayesSYNC")
```

### Authors and license

Selima Jaoua (University of Zürich, Switzerland), Daniel Temko & Hélène
Ruffieux (University of Cambridge, UK).

This software uses the GPL v3 license. Authors and copyright are also
provided in [DESCRIPTION](DESCRIPTION).

### Issues

To report an issue, please use the [bayesSYNC issue
tracker](https://github.com/hruffieux/bayesSYNC/issues) at github.com.
