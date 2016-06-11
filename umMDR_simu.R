library(umMDR)

source('gene_simu_nLD.R')    ## generate simulation data without LD
source('M_QMDR.R')           ## for running MDR/QMDR



## calculate type I error--no causal SNPs
## input: --- N: sample size
##----------- nsnp: number of snps
##------------trait.type: type of trait, either "quantitative" or "binary"
##------------cova: covariates
##------------maf : 0.2 or 0.4
##------------repts: total number of run for each penatrace model
##------------alpha: significant level
##------------nperm: small number of permutation times to estimate non-central parameter
##----Output: the type I error and the raw type I error (without non-central correction)
CalType1 <- function(N, nsnp = 2, trait.type = 'quantitative',
                     cova = NULL, maf = 0.2, repts = 100, alpha = 0.05, nperm = 5){
  K = 2   # 2 way interaction, first 2 are causal
  output = NULL

  set.seed(1)
  error.raw = error = rep(0, repts)
  snp.dat = list()
  for(run in 1:repts) snp.dat[[run]] =  simu_popu(N, M1 = K, M0 = nsnp-K, p.f1 = rep(maf, K), p.f0 = maf,
                                                  pen1 = pen[, 1], mu1 = 1, trait.type)$data[, -1]
  for(run in 1:repts){
    #generate data under the null
    snp.mats = snp.dat[[run]]
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

    res = umMDR(snp.mats, K = K, phe = phe, cova = cova, classm, adj.main = FALSE, nperm)


    error[run] = ifelse(res$pvs[1] < alpha, 1, 0)         ## error

    error.raw[run] = ifelse(res$raw.pvs[1] < alpha, 1, 0)    ## error (without non-central correction) under bonferoni correction

  }

  merror = mean(error)
  merror.raw = mean(error.raw)


  output = rbind(output, c(maf, merror,  merror.raw))
  colnames(output) = c("maf", "error", "error.raw")
  return(output)
}
CalType1 = cmpfun(CalType1)



## calculate power under 70 penetrance models for Case1-3
## input: ---- N: sample size
##------------nsnp: number of snps, first 2 are causal
##------------classm: classification of H/L, either "mean" or "obs-ratio"
##------------cova: covariates
##------------adj.main: adjust main effect of not
##------------repts: total number of run for each penatrace model
##------------simu.type: 'Case1'---SNP1,SNP2 have causal effect, no marginal effect, binary
##---------------------: 'Case2'---SNP1,SNP2 have causal effect, no marginal effect, binary
##---------------------: 'Case3'---SNP1,SNP2 are causal and SNP3 has marginal effect (marg.efsize))
##------------alpha: significant level
##------------marg.efsize: marginal effect size
##------------nperm: small number of permutation to estimate non-central parameter
##------------cand.pmodels: candidate penetrance model, subset of 1:70
##----Output: power and the raw power (without non-central correction), both corrected
##----------  by bonferroni correction
##------------also power.PRank (the causal model been ranked as no.1), power of mdr (or qmdr)(power.mdr)
##------------cvc of mdr (or qmdr) (cvc.mdr)
##------------note that power.raw is meaningless since it does not control type I error
CalPower <- function(N, nsnp = 10,  cova = NULL, adj.main = FALSE, repts = 100, simu.type = 'Case1', mu = 1,
                     marg.efsize = 1, alpha = 0.05, nperm = 5, cand.pmodels = 1){

  K = 2   # 2 way interaction, first 2 are causal
  output = NULL
  for(model.no in cand.pmodels){
    set.seed(1)
    trait.type = 'quantitative'
    if(simu.type == 'Case1') trait.type = 'binary'
    maf = MAF[model.no]
    pow.raw = pow = pow.rank = pow.qmdr = cvcs = rep(0, repts)
    dats = list()
    for(run in 1:repts) {
      dats[[run]] = simu_popu(N, M1 = K, M0 = nsnp-K, p.f1 = rep(maf, K), p.f0 = maf,
                              pen1 = pen[, model.no], mu1 = mu, trait.type)$data
    }

    for(run in 1:repts){
      #generate data under the null

      dat = dats[[run]]

      phe = as.matrix(dat[, 1])
      snp.mats = dat[, -1]
      p = ncol(snp.mats)

      snp.combs <- combn(p, K)   ## all SNP pairs
      ns = ncol(snp.combs)


      if(simu.type == 'Case3') phe = as.matrix(phe + rnorm(N, marg.efsize * dat[, 4]))  ## s3 marginal

      classm = 'mean'
      if(simu.type == 'Case1') classm = 'obs-ratio'

      res = umMDR(snp.mats, K = K, phe = phe, cova = NULL, classm, adj.main, nperm)


      pow[run] = ifelse(res$pvs[1] < alpha/ns, 1, 0)         ## power under bonferoni correction
      pow.rank[run] = ifelse(which.min(res$pvs) == 1, 1, 0) ## only care about pv rank

      pow.raw[run] = ifelse(res$raw.pvs[1] < alpha/ns, 1, 0)    ## power (without non-central correction) under bonferoni correction

      # for qmdr

      if(trait.type == 'binary'){
        res = kway_MDR(phe,  K, snp.mats, sele.type = 'cvc')
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
    colnames(output) = c("model.id", "power.PBonf", "power.PRank", "power.raw",
                         "power.mdr", "cvc.mdr")

  }
  return(output)
}
CalPower = cmpfun(CalPower)


## calculate power for Case4
CalPower_3way <- function(N, nsnp = 10,  cova = NULL, adj.main = TRUE, repts = 100,
                          mu1 = 0.2, nperm = 5, maf = 0.1, marg.eff = 1, trait.type = 'quantitative'){

  K = 3   # 3 way interaction, first 3 are causal
  output = NULL
  set.seed(1)
  classm = ifelse(trait.type == 'quantitative', 'mean', 'obs-ratio')
  pow.raw = pow = pow.rank = pow.qmdr = cvcs = rep(0, repts)
  dats = list()

  for(run in 1:repts) {
    dats[[run]] <- simu_popu3(N, M1=K, M0=nsnp-K, p.f1=rep(maf, K), p.f0=maf,
                              mu1 = mu1, trait.type)$data
  }

  for(run in 1:repts){
    #generate data under the null

    dat = dats[[run]]

    phe = as.matrix(dat[, 1])
    snp.mats = dat[, -1]
    p = ncol(snp.mats)

    snp.combs <- combn(p, K)   ## all SNP pairs
    ns = ncol(snp.combs)

    if(trait.type == 'quantitative') phe = as.matrix(phe + rnorm(N, marg.eff * snp.mats[, 4]))  ## snp4 marginal
    res = umMDR(snp.mats, K = K, phe = phe, cova = cova, classm = classm,
                adj.main = adj.main, nperm = nperm)


    pow[run] = ifelse(res$pvs[1] < 0.05/ns, 1, 0)         ## power under bonferoni correction
    pow.rank[run] = ifelse(which.min(res$pvs) == 1, 1, 0) ## only care about pv rank

    pow.raw[run] = ifelse(res$raw.pvs[1] < 0.05/ns, 1, 0)    ## power (without non-central correction) under bonferoni correction

    # for mdr or qmdr
    if(classm == 'obs-ratio'){
      res = kway_MDR(phe, K, snp.mats, sele.type = 'cvc')
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

  output = rbind(output, c(maf, mpow, mpow.rank,  mpow.raw, mpow.qmdr, mcvc))
  colnames(output) = c("MAF", "power.PBonf", "power.PRank", "power.raw",
                       "power.mdr", "cvc.mdr")

  return(output)
}
CalPower_3way = cmpfun(CalPower_3way)




## examples to calculate type I error
CalType1(N = 400, nsnp = 2, maf = 0.1, trait.type = 'quantitative',
         repts = 500, nperm = 5, alpha = 0.05)

CalType1(N = 400, nsnp = 2, maf = 0.2, trait.type = 'quantitative',
         repts = 500, nperm = 5, alpha = 0.05)

CalType1(N = 400, nsnp = 2, maf = 0.1, trait.type = 'binary',
         repts = 500, nperm = 5, alpha = 0.05)
CalType1(N = 400, nsnp = 2, maf = 0.2, trait.type = 'binary',
         repts = 500, nperm = 5, alpha = 0.05)



## examples to calculate power on some of the 70 models with 100 repeats for each one
## this may take a while if you try all models at a time, so try a few,
## for example: cand.pmodels = 11:12
## note that power.raw is meaningless since it does not control type I error


CalPower(N = 400, nsnp = 10, simu.type = 'Case1',  adj.main = FALSE,
         repts = 100, nperm = 5, alpha = 0.05, cand.pmodels = 11:12)

CalPower(N = 400, nsnp = 10, simu.type = 'Case2',  adj.main = FALSE,
         repts = 100, nperm = 5, alpha = 0.05, cand.pmodels = 11:12)

## adjusting marginal effects
CalPower(N = 400, nsnp = 10, simu.type = 'Case3', adj.main = TRUE,
         repts = 100, nperm = 5, alpha = 0.05, cand.pmodels = 11:12)

## Case4
CalPower_3way(N = 400, nsnp = 10,  adj.main = TRUE, repts = 100,
                     nperm = 5, maf = 0.4, mu1 = 0.5, marg.eff = 1)

