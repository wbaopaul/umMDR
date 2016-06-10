# generate many snps -phenotypes under null for drawing qqplot

library(compiler)
library(MASS)
library(glmnet)
source('umMDR_source.R')
source('gene_simu_nLD.R')

## simu.type should be simu0

N = 400
repts = 2000
cvts = NULL
K = 2           # K-way interaction

classm = 'mean'
trait.type = 'quantitative'
maf = 0.2
cal_pv <- function(maf){
  output = NULL
  set.seed(1)
  model.no = 1
  pv.raw = pv = rep(0, repts)
  for(run in 1:repts){
    #generate data

    dat0 <- simu_popu(N, M1=2, M0=0, p.f1=rep(maf, 2), p.f0=maf,
                      pen1 = pen[, model.no], mu1 = 1, trait.type)

    dat = dat0$data

    snp.mats = dat[, -1]
    p = ncol(snp.mats)
    snp.combs <- combn(p, K)
    ns = ncol(snp.combs)
    if(trait.type == 'binary'){
      phe = as.matrix(rbinom(N, 1, 0.5))  ## generate null data
    }else{
      phe = as.matrix(rnorm(N))  ## generate null data
    }


    res = umMDR(snp.mats, phe, K, cvts, classm, adj.main)

    pv[run] = res$pvs[1]

    pv.raw[run] =res$raw.pvs[1]
  }


  res = list('pv' = pv, 'rpv' = pv.raw)
  #save(pvs, file = 'pvs-null.Rdata')
  return(pvs)
}
cal_pv = cmpfun(cal_pv)

res <- cal_pv(maf)



postscript('qq-plot.eps', width = 900, height = 350)
par(mfrow = c(1, 2))
qqplot(y = res$rpv, qunif(ppoints(1000), 0, 1), xlab = 'Expected P-value',
       ylab = 'Empirical P-value', col=2)
qqline(y = res$rpv, distribution = function(p) qunif(p, 0, 1), lwd=2)
#abline(0, 1, lwd = 2)

qqplot(y = res$pv, qunif(ppoints(1000), 0, 1), xlab = 'Expected P-value',
       ylab = 'Empirical P-value', col=3)
qqline(y = res$pv, distribution = function(p) qunif(p, 0, 1), prob = c(0.25, 0.75),
        lwd = 2)
#abline(0, 1, lwd = 2)
dev.off()

