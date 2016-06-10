# umMDR
A unified model based MDR framework for detecting gene-gene interaction


## Description
 A two-step unified model based MDR approach (UM-MDR), in which, the significance of a multi-locus model, even a high-order model, can be easily obtained through a regression framework, and a semiparametric correction procedure is designed for control type I error rates


Key idea: assume the null distributioin is non-central chisquare for testing association between phenotype and S;
the non-cetral parameter is estimated through a small number of permulation (say 5 or 10 times)

## Installation
* Download the package umMAD_0.1.tar.gz
* To install, in R type: 
  ```
   install.packages("path/umMDR_0.1.tar.gz", repos = NULL, type = "source")
   
  ```
## Usage 
  ```
  umMDR(snp.all, phe, K = 2, cova = NULL, classm = "mean",
  adj.main = "FALSE", nperm = 5)
  ```
 
## Example 
  ```
   library(umMDR)
    ## an example
    snps <- matrix(rbinom(100 * 5, 1, 0.2), nrow = 100)  ## generate 5 snps
    phe <- rnorm(100)  ## generate phenotype
    umMDR(snps, phe, 2)
  ```

## More information of usage
 
```
library(umMDR)
help(umMDR)
```



## Citation
Yu, Wenbao and Park, Taesung "A unified model based multifactor dimensionality reduction framework for detecting gene-gene interactions", submitted, 2016

## Other Sources for the manuscript

gene_simu_nLD.R --- generate simulation data under 70 penetrace function and a 3-way interaction model

umMDR_simu.R --- run simulation Case1-case4

qqplot_simu.R --- generate qqplot in the manuscript

M_QMDR.R --- run MDR, QMDR or Multi-QMDR as a comparison

