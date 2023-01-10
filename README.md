# StanProt
Bayesian Modelling for Label Free MS Based Proteomics

## What is StanProt ?

`StanProt` is a R package that provides Bayesian models for Label 
Free MS based proteomics data. It allows a.o. to model the data missingness 
process simultaneously to the MS intensities.

On the one hand, `StanProt` proposes some S4 classes to standardise the access, 
manipulation and interpretation of Bayesian models generated with the Stan
MCMC Engine (exact inference). See [here](https://mc-stan.org/).

On the other hand, `StanProt` provides also support for approximate 
inference, i.e. MFVB approximation (using built-in code), and Laplace 
approximation (using also Stan as en engine).

## License

The `StanProt` code is provided under [GPL license version 3.0 or 
higher](https://opensource.org/licenses/GPL-3.0).  
The documentation, 
including the manual pages and the vignettes, are distributed under a [CC BY-SA 
4.0 license](https://creativecommons.org/licenses/by-sa/4.0/).

## Citing StanProt

If you use `StanProt` in your research, please cite the following paper:  

Hauchamps, Philippe. A Bayesian Approach combining Peptide Intensity and Missingness Modelling to analyse Label Free Mass Spectrometry based Proteomics Data. Faculté des sciences, Université catholique de Louvain, 2021. Prom. : Lambert, Philippe ; Gatto, Laurent. You can find a link to the manuscript [here](https://dial.uclouvain.be/memoire/ucl/object/thesis:29911).


## Installation Instructions 

It is recommended to install first the `rstan` package and to make sure it works. 
See the [rstan quick start guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

As soon as this is done, and the installation is verified by running a test model, 
the next step is to install the `StanProt` package from the github repo:

```
devtools::install_github("https://github.com/lgatto/StanProt")
```