# umMDR
A unified model based MDR framework for detecting gene-gene interaction


## Description:
 A two-step unified model based MDR approach (UM-MDR), in which, the significance of a multi-locus model, even a high-order model, can be easily obtained through a regression framework, and a semiparametric correction procedure is designed for control type I error rates


Key idea: assume the null distributioin is non-central chisquare for testing association between phenotype and S;
the non-cetral parameter is estimated through a small number of permulation (say 5 or 10 times)

## Installation and Usage
* Download the package MaxmzpAUC_0.1.tar.gz
* In R, 
  * To install: 
  ```
   install.packages("path/umMDR_0.1.tar.gz", repos = NULL, type = "source")
   
  ```
  
  * To use: 
  ```
   library(umMDR)
    help(umMDR)
    
  ```
  
## Citation
Yu, Wenbao and Park, Taesung "A unified model based multifactor dimensionality reduction framework for detecting gene-gene interactions", submitted, 2016

## Notice
gene_simu_nLD.R --- generate simulation data under 70 penetrance models and no LD

gene_simu_3way.R --- generate 3 way interaction data under the model defined in the manuscript

umMDR_simu.R --- give some examples of executing UMDR for calculating type I error and power, the 
                        powers of QMDR or MDR are also provided

M-QMDR-simp.R --- run MDR, QMDR or M-QMDR as a comparison
