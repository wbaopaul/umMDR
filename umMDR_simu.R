## simulation study in the manuscript

library(umMDR) ## source codes for umdr
source('gene_simu_nLD.R')    ## generate data without LD, Case1-3
source('M-QMDR-simp.R')      ## for running MDR/QMDR
source('gene_simu_3way.R')   ## generate data for case4


## calculate type I error--no causal SNPs
## input: --- N: sample size
##----------- nsnp: number of snps
##------------trait.type: type of trait, either "quantitative" or "binary"
##------------cova: covariates
##------------adj.main: adjust main effect of not
##------------maf : 0.2 or 0.4
##------------repts: total number of run for each penatrace model
##------------alpha: significant level
##------------nperm: small number of permutation times to estimate non-central parameter
##----Output: the type I error and the raw type I error (without non-central correction)
CalType1 <- function(N, nsnp = 10, trait.type = 'quantitative',
            cova = NULL, maf = 0.2, repts = 100, alpha = 0.05, nperm = 10){
    K = 2   # 2 way interaction, first 2 are causal
    output = NULL

    set.seed(1)
    model.no = which(MAF == maf)[1]
    error.raw = error = rep(0, repts)
    for(run in 1:repts){
      #generate data under the null
      dat <- simu_popu(N, M1 = K, M0 = nsnp-K, p.f1 = rep(maf, K), p.f0 = maf,
                        pen1 = pen[, model.no], mu1 = 1, trait.type)$data
      snp.mats = dat[, -1]
      p = ncol(snp.mats)

      snp.combs <- combn(p, K)   ## all SNP pairs
      ns = ncol(snp.combs)
      if(trait.type == 'binary'){
        phe = as.matrix(rbinom(N, 1, 0.5))  ## generate null data
        classm = 'obs-ratio'
      }else{
        phe = as.matrix(rnorm(N))  ## generate null data
        classm = 'mean'
      }

      res = umMDR(snp.mats, K, phe, cova = cova, classm, adj.main = FALSE, nperm)


      error[run] = ifelse(res$pv[1] < alpha/ns, 1, 0)         ## error under bonferoni correction

      error.raw[run] = ifelse(res$rpv[1] < alpha/ns, 1, 0)    ## error (without non-central correction) under bonferoni correction

    }

    merror = mean(error)
    merror.raw = mean(error.raw)


    output = rbind(output, c(maf, merror,  merror.raw))
    colnames(output) = c("maf", "error", "error.raw")
  return(output)
}
CalType1 = cmpfun(CalType1)



## calculate power under 70 penetrance models
## input: ---- N: sample size
##------------nsnp: number of snps, first 2 are causal
##------------trait.type: type of trait, either "quantitative" or "binary"
##------------classm: classification of H/L, either "mean" or "obs-ratio"
##------------cova: covariates
##------------adj.main: adjust main effect of not
##------------repts: total number of run for each penatrace model
##------------simu.type: 'Case1'---SNP1,SNP2 have causal effect, no marginal effect, binary
##---------------------: 'Case2'---SNP1,SNP2 have causal effect, no marginal effect, binary
##---------------------: 'Case3'---SNP1,SNP2 are causal and SNP3 has marginal effect (marg.efsize))
##------------mu: effect size of snp1 and snp2;
##------------alpha: significant level
##------------marg.efsize: marginal effect size
##------------nperm: small number of permutation to estimate non-central parameter
##------------cand.pmodels: candidate penetrance model, subset of 1:70
##----Output: power and the raw power (without non-central correction), both corrected
##----------  by bonferroni correction
##------------also power.rank (the causal model been ranked as no.1), power of qmdr(power.qmdr)
##------------cvc of qmdr (cvc.qmdr)
##------------note that power.raw is meaningless since it does not control type I error
CalPower <- function(N, nsnp = 10,  cova = NULL, adj.main = FALSE, repts = 100, simu.type = 'Case1', mu = 1,
                     marg.efsize = 0.5, alpha = 0.05, nperm = 10, cand.pmodels = 1:70){

  K = 2   # 2 way interaction, first 2 are causal
  output = NULL
  for(model.no in cand.pmodels){
    set.seed(1)
    maf = MAF[model.no]
    pow.raw = pow = pow.rank = pow.qmdr = cvcs = rep(0, repts)
    for(run in 1:repts){
      #generate data under the null

      trait.type = 'quantitative'
      if(simu.type == 'Case1') trait.type = 'binary'

      dat0 <- simu_popu(N, M1 = K, M0 = nsnp-K, p.f1 = rep(maf, K), p.f0 = maf,
                        pen1 = pen[, model.no], mu1 = mu, trait.type)

      dat = dat0$data

      phe = as.matrix(dat[, 1])
      snp.mats = dat[, -1]
      p = ncol(snp.mats)

      snp.combs <- combn(p, K)   ## all SNP pairs
      ns = ncol(snp.combs)


      if(simu.type == 'Case3') phe = as.matrix(phe + rnorm(N, marg.efsize * dat[, 4]))  ## s3 marginal

      classm = 'mean'
      if(simu.type == 'Case1') classm = 'obs-ratio'

      res = umMDR(snp.mats, K, phe, cova = NULL, classm, adj.main, nperm)


      pow[run] = ifelse(res$pv[1] < alpha/ns, 1, 0)         ## power under bonferoni correction
      pow.rank[run] = ifelse(which.min(res$pv) == 1, 1, 0) ## only care about pv rank

      pow.raw[run] = ifelse(res$rpv[1] < alpha/ns, 1, 0)    ## power (without non-central correction) under bonferoni correction

      # for qmdr

      if(trait.type == 'binary'){
        res = kway_MDR('mdr', phe, cova, K, snp.mats, sele.type = 'cvc')
      }else{
        res = kway_QMDR(phe, K, snp.mats, sele.type = 'cvc')
      }

      pow.qmdr[run] = ifelse(res[1] == 1, 1, 0)
      cvcs[run] = res[4]


    }

    mpow = mean(pow)
    mpow.rank = mean(pow.rank)
    mpow.raw = mean(pow.raw)
    mpow.qmdr = mean(pow.qmdr)
    mcvc = mean(cvcs)

    output = rbind(output, c(model.no, mpow, mpow.rank,  mpow.raw, mpow.qmdr, mcvc))
    colnames(output) = c("model.no", "power", "power.rank", "power.raw", "power.qmdr", "cvc.qmdr")

  }
  return(output)
}
CalPower = cmpfun(CalPower)



## examples to calculate type I error
CalType1(N = 100, nsnp = 2, maf = 0.2, trait.type = 'quantitative',
         repts = 100, nperm = 5, alpha = 0.05)

CalType1(N = 100, nsnp = 2, maf = 0.4, trait.type = 'quantitative',
         repts = 100, nperm = 5, alpha = 0.05)

CalType1(N = 100, nsnp = 2, maf = 0.2, trait.type = 'binary',
         repts = 100, nperm = 5, alpha = 0.05)
CalType1(N = 100, nsnp = 2, maf = 0.4, trait.type = 'binary',
         repts = 100, nperm = 5, alpha = 0.05)



## examples to calculate power on some of the 70 models with 100 repeats for each one
## this may take a while if you try all models at a time, so try a few,
## for example: cand.pmodels = 1:2
## note that power.raw is meaningless since it does not control type I error
CalPower(N = 400, nsnp = 10, simu.type = 'Case1',  adj.main = FALSE,
         repts = 100, nperm = 5, alpha = 0.05, cand.pmodels = 1:2)
CalPower(N = 400, nsnp = 10, simu.type = 'Case2',  adj.main = FALSE,
         repts = 100, nperm = 5, alpha = 0.05, cand.pmodels = 1:2)

## adjusting marginal effects
CalPower(N = 400, nsnp = 10, simu.type = 'Case3', adj.main = TRUE,
         repts = 100, nperm = 5, alpha = 0.05, cand.pmodels = 1:2)

