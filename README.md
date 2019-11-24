# ACTOR: A latent dirichlet model to Compare expressed isoform proportions TO a Reference panel

Overview
=====
`actor` is an R package which takes transcript expression as input and uses a latent Dirichlet model to compare to GTEx tissues. `actor` provides gene level similarities to the GTEx tissues and identifies collections of genes, or gene classes, which align similarly to GTEx tissues. 

Installing ACTOR
================
ACTOR is available as an R package which can be installed using `devtools` as follows.
```{r}
library(devtools)
devtools::install_github("mccabes292/actor")
```



Vignettes
=========
A vignette for running your own analysis using `actor` is available [here](https://htmlpreview.github.io/https://github.com/mccabes292/actor/tree/master/vignettes/actor.Rmd)
