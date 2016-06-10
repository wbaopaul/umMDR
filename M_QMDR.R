##=======================================================================##
## **** this version is the simple version of M-QMDR-reuse.R
## **** only consider MDR and QMDR, not calculate the pv by permutation
## **** note that the best model of MDR by cvc and BA
##=======================================================================##


library(compiler)
library(MASS)

##------------------------------------------------------------------------##
## list functions
## calculating t score
cal_tstat <- function(ids, high.all, S){
  high.ids = intersect(ids, high.all)
  low.ids = setdiff(ids, high.ids)

  s1 = S[high.ids]

  s2 = S[low.ids]

  if(length(high.ids) == 0 || length(low.ids) == 0) return(0)
  stat = t.test(s1, y = s2, var.equal = TRUE)$statistic

  return(abs(stat))
}
cal_tstat = cmpfun(cal_tstat)

## calculating ht2 score
cal_ht2 <- function(ids, high.all, phes){

  high.ids = intersect(ids, high.all)
  low.ids = setdiff(ids, high.ids)
  d = ncol(phes)

  s1 = phes[high.ids, ]

  s2 = phes[low.ids, ]

  if(length(high.ids) == 0 || length(low.ids) == 0) return(0)

  s1 = as.matrix(s1)
  s2 = as.matrix(s2)

  if(ncol(s1) == 1) s1 = t(s1)
  if(ncol(s2) == 1) s2 = t(s2)


  stat = dire_ht2(s1, s2, phes)$fstat                   ## another version
  ## stat is scaled that it follows a F distribution with degree d, n-1-d under the null
  return(stat)
}


## calculate training score for a given train samples id and
## snp.dat, n by k, -- k column correspond to the specific k-snps
## snp.dat could have missing values and test score was calculated for each model
qmdr <- function(train.ids, test.ids, snp.dat, S, phes, test.type = 'ht2'){
  snp.dat = as.matrix(snp.dat)

  ## remove missing values
  fids = which(complete.cases(snp.dat))
  snp.dat = as.matrix(snp.dat[fids, ])

  test.ids = intersect(test.ids, fids)
  train.ids = intersect(train.ids, fids)

  #sbar = mean(S[train.ids])
  sbar = mean(S)
  k = ncol(snp.dat)

  ## split data into cells
  tlist = vector('list', k)
  for(i in 1:k) tlist[[i]] = snp.dat[, i]
  cells = split(data.frame(cbind(fids, snp.dat)), tlist)

  ## delete NULL cells
  obs.cell = sapply(cells, function(x) nrow(x))
  cell.null = which(obs.cell == 0)
  if(length(cell.null) > 0 ) cells = cells[-cell.null]

  ## get trainid in each cell
  cells.trainid = lapply(cells, function(x) return(intersect(x[, 1], train.ids)))

  cells.num = length(cells)
  #cells.label = rep(0, cells.num)
  high.all  = NULL
  for(i in 1:cells.num){
    temp.ids = cells.trainid[[i]]
    if(length(temp.ids) == 0) next
    if (mean(S[temp.ids]) >= sbar){
      high.all = c(high.all, cells[[i]][, 1])
    }
  }

  if(test.type == 't'){
    train.stat = cal_tstat(train.ids, high.all, S)
    test.stat = cal_tstat(test.ids, high.all, S)
  }else{
    train.stat = cal_ht2(train.ids, high.all, phes)
    test.stat = cal_ht2(test.ids, high.all, phes)
  }
  return(list('train.stat' = train.stat, 'test.stat' = test.stat ))
}
qmdr = cmpfun(qmdr)


mdr <- function(train.ids, test.ids, snp.dat, HL){
  snp.dat = as.matrix(snp.dat)

  ## remove missing values
  fids = which(complete.cases(snp.dat))
  snp.dat = as.matrix(snp.dat[fids, ])

  test.ids = intersect(test.ids, fids)
  train.ids = intersect(train.ids, fids)

  gratio = mean(HL)
  k = ncol(snp.dat)
  n = length(fids)

  ## split data into cells
  tlist = vector('list', k)
  for(i in 1:k) tlist[[i]] = snp.dat[, i]
  cells = split(data.frame(cbind(fids, snp.dat)), tlist)

  ## delete NULL cells
  obs.cell = sapply(cells, function(x) nrow(x))
  cell.null = which(obs.cell == 0)
  if(length(cell.null) > 0 ) cells = cells[-cell.null]

  ## get trainid in each cell
  cells.trainid = lapply(cells, function(x) return(intersect(x[, 1], train.ids)))

  cells.num = length(cells)
  #cells.label = rep(0, cells.num)
  high.all  = NULL
  for(i in 1:cells.num){
    temp.ids = cells.trainid[[i]]
    if(length(temp.ids) == 0) next
    if (mean(HL[temp.ids]) >= gratio){
      high.all = c(high.all, cells[[i]][, 1])
    }
  }

  pHL = rep(0, n)
  pHL[high.all] = 1

  table1 = table(pHL[train.ids], HL[train.ids])
  table2 = table(pHL[test.ids], HL[test.ids])

  if(any(dim(table1) == 1)) {
    train.ba = length(which(pHL[train.ids]==HL[train.ids]))/length(train.ids)
  }else{
    train.ba = 0.5*(table1[1, 1]/ sum(table1[, 1]) +
                      table1[2, 2]/ sum(table1[, 2]))
  }

  if(any(dim(table2) == 1)) {
    test.ba = length(which(pHL[test.ids]==HL[test.ids]))/length(test.ids)
  }else{
    test.ba = 0.5*(table2[1, 1]/ sum(table2[, 1]) +
                     table2[2, 2]/ sum(table2[, 2]))
  }

  return(list('train.ba' = train.ba, 'test.ba' = test.ba ))
}
mdr = cmpfun(mdr)


## choose the best k-way interaction by cross validation
## return the snp combination and its test score
## snp.combs includes a possible k-way snps combinations in each column
## using cv consistency or testing score to decide best model
## S -- summary score (for multi-MDR)
## phes -- all phenotypes
CV_QMDR <- function(folds = 10, snp.all, S, phes,
                             test.type = 'ht2', sele.type = 'cvc', snp.combs){
  ns = ncol(snp.combs)
  n = length(S)

  ## split the whole data into folds
  cvlen = floor(n/folds)
  cc = 1:n
  test.stats = train.stats = rep(0, ns)
  temptest.stats = matrix(0, ns, folds)
  best.comb = rep(0, folds)

  ## select best model(i.e. snp combination)
  for(i in 1:folds){
    testid = ((i - 1) * cvlen + 1) : (i * cvlen)
    trainid = cc[-testid]

    for(j in 1:ns){
      temp.result = qmdr(trainid, testid, snp.all[, snp.combs[, j]], S, phes, test.type)
      train.stats[j] = temp.result$train.stat
      temptest.stats[j, i] = temp.result$test.stat
    }
    # which snp pair has best training stat for each trainind set
    best.comb[i] = j0 = which.max(train.stats)
  }

  test.stats = rowMeans(temptest.stats, na.rm = TRUE)  ## average testing stats for all snp pairs

  if(sele.type == 'cvc'){
    ta = table(best.comb)
    cvc = max(ta)          ## the largest cvc
    best.pair = as.numeric(names(which(ta == cvc)))[1]  ## the pair gets largest cvc
  }
  if(sele.type == 'score'){
    best.pair = which.max(test.stats)  ## the pair gives larges test score
    cvc = length(which(best.comb == best.pair))
  }

  sele.score = test.stats[best.pair]

  ## sele.score -- corresponding to the test score of the final selected model
  ## test.stats -- record test.scores for all possible k-way model (snp interactions)

  ## Save the test score of the best model in each cv
  scores.cv = temptest.stats[best.pair, ]

  return(list('cvc' = cvc, 'score' = sele.score, 'scores.cv' = scores.cv,
              'best.pair' = best.pair, 'test.stats' = test.stats))

}
CV_QMDR = cmpfun(CV_QMDR)


CV_MDR <- function(folds = 10, snp.all, snp.combs, sele.type = 'cvc', HL){
  ns = ncol(snp.combs)
  n = length(HL)

  ## split the whole data into folds
  cvlen = floor(n/folds)
  cc = 1:n
  test.ba = train.ba = rep(0, ns)
  temptest.ba = matrix(0, ns, folds)
  best.comb = rep(0, folds)

  ## select best model(i.e. snp combination)
  for(i in 1:folds){
    testid = ((i - 1) * cvlen + 1) : (i * cvlen)
    trainid = cc[-testid]

    for(j in 1:ns){
      temp.result = mdr(trainid, testid, snp.all[, snp.combs[, j]], HL)
      train.ba[j] = temp.result$train.ba
      temptest.ba[j, i] = temp.result$test.ba
    }
    # which snp pair has best training stat for each trainind set
    best.comb[i] = j0 = which.max(train.ba)
  }

  test.ba = rowMeans(temptest.ba, na.rm = TRUE)  ## average testing stats for all snp pairs

  if(sele.type == 'cvc'){
    ta = table(best.comb)
    cvc = max(ta)          ## the largest cvc
    best.pair = as.numeric(names(which(ta == cvc)))[1]  ## the pair gets largest cvc
  }
  if(sele.type == 'score'){
    best.pair = which.max(test.ba)  ## the pair gives larges test score
    cvc = length(which(best.comb == best.pair))
  }

  sele.score = test.ba[best.pair]

  ## sele.score -- corresponding to the average test score of the final selected model
  ## test.ba -- record average test.scores for all possible k-way model (snp interactions)

  ## Save the test score of the best model in each cv
  scores.cv = temptest.ba[best.pair, ]

  return(list('cvc' = cvc, 'score' = sele.score, 'scores.cv' = scores.cv,
              'best.pair' = best.pair, 'test.ba' = test.ba))

}
CV_MDR = cmpfun(CV_MDR)


kway_QMDR <- function(phe, K, snp.mats, test.type = 't', sele.type = 'cvc'){

  #set.seed(1)
  n = length(phe)
  p = ncol(snp.mats)
  snp.combs <- combn(p, K)  ## all possible combinatory pairs
  ns  = ncol(snp.combs)
  test.stats = rep(0, ns)

  aa = 1:n
  #aa = sample(1:n, n)  ## shuffle samples

  result = CV_QMDR(10, snp.mats[aa, ], phe[aa], phe[aa], test.type,
                            sele.type, snp.combs)

  model.cons = result$cvc
  model.sele = result$best.pair
  model.score = result$score
  test.stats = result$test.stats
  scores.cv = result$scores.cv

  best.ksnps = snp.combs[, model.sele]


  return(c(model.sele, best.ksnps, model.cons))

}
kway_QMDR = cmpfun(kway_QMDR)


## for mdr
kway_MDR <- function(SS, K, snp.mats, sele.type = 'cvc'){

  gm = mean(SS)
  HL = ifelse(SS > gm, 1, 0)   ## change phenotypes into two groups (for gmdr)

  if(length(unique(SS)) == 2) HL = SS
  #set.seed(1)
  n = length(SS)
  p = ncol(snp.mats)
  snp.combs <- combn(p, K)  ## all possible combinatory pairs
  ns  = ncol(snp.combs)
  test.ba = rep(0, ns)  ## average test score for each pair

  aa = 1:n
  #aa = sample(1:n, n)  ## shuffle samples

  result = CV_MDR(10, snp.mats[aa, ], snp.combs, sele.type, HL[aa])

  model.cons = result$cvc
  model.sele = result$best.pair
  model.score = result$score
  test.ba = result$test.ba
  scores.cv = result$scores.cv

  best.ksnps = snp.combs[, model.sele]

  return(c(model.sele, best.ksnps, model.cons))

}
kway_MDR = cmpfun(kway_MDR)




