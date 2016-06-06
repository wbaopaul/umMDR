## ---------------------------------------------------------------##
## this program is written by Wenbao yu
## to generate simulation data for MDR work
## It is modified based on Min-Seok Kwon's code for Multi-QMDR
##----------------------------------------------------------------##

require(MASS)
library(compiler)

## Generate genotype of a single SNP with fixed MAF
simu_sSNP_genotype <- function(n, p){
  ## n: sample size, p: MAF
  ## Additive mode (0,1,2), under HWE, no LD
  x <- rmultinom(n, 1, c(p^2, 2 * p * (1-p), (1-p)^2))
  y <- x[2, ] + 2 * x[1, ]
  c(y, p)
}
simu_sSNP_genotype = cmpfun(simu_sSNP_genotype)


## Generate phenotype data via penetrance function of M-locus model
simu_trait <- function(dG, pen1, z, mu1,  trait.type = 'quantitative'){
  ## dG: genotypes of causal SNPs  (col=SNP, row=individual)
  ## pen: penetrance function
  n <- nrow(dG)
  M <- ncol(dG)

  geno <- rep(1,n)

  for(j in 1:M) {
    geno <- geno + dG[, j] * (3^(j-1))  # geno with be range from 1 to 9
  } ## geno[j] = 3 means the j-th observation locates in the 3-rd cell


  y0 <- NULL
  pens = pen1[geno]

  if(trait.type == 'binary'){
    for(i in 1:n) y0 <- rbind(y0, rbinom(1, 1, pens[i]))
  }else{
    #sigma <- matrix(c(1, rho, rho, 1), 2, 2)

    #for (i in 1:n) y0 <- rbind(y0, mvrnorm(1, mu = c(mu1 * pens[i], mu2 * pens[i]), sigma))
    for (i in 1:n) y0 <- rbind(y0, rnorm(1, mean = mu1 * pens[i]))

  }

  return(y0)
}
simu_trait = cmpfun(simu_trait)


## Generate a population data
simu_popu <- function(N, M1 = 2, M0 = 0, p.f1, p.f0, pen1, mu1, trait.type = 'quantitative'){
  ## N: population size
  ## M1: number of causal SNPs
  ## M0: number of non-causal SNPs
  ## p.f1: fixed MAF of causal SNPs
  ## p.f0: fixed MAF of non-causal SNPs
  ## pen1: penetrance function being used
  ## rslt$data: (col1:ID, col2:Trait, other cols:SNPs)

  #z <- rnorm(N, cov_mean, sqrt(cov_var))  ## covariates
  z <- 0;

  ## Causal SNP genotype generation
  dGen1 <- matrix(rep(0, N * M1), nrow = N, ncol = M1)
  dP1 <- rep(0, M1)
  for(i in 1:M1){
    res <- simu_sSNP_genotype(n = N, p = p.f1[i])
    dGen1[, i] <- res[1:N]
    dP1[i] <- res[(N+1)]  ## MAF

  }

  ## Non-causal SNP genotype generation
  dGen0 = NULL
  if(M0 > 0){
    dGen0 <- matrix(rep(0, N * M0), nrow = N, ncol = M0)
    dP0 <- rep(0, M0)
    for(i in 1:M0){
      res <- simu_sSNP_genotype(n = N, p = p.f0)
      dGen0[, i] <- res[1:N]
      dP0[i] <- res[(N+1)]
    }
  }

  snpnames <- c()
  for (i in 1:(M1 + M0)){
    snpnames <- c(snpnames, paste("S", i, sep=""))
  }


  ## Trait generation
  dTrt <- simu_trait(dG = dGen1, pen = pen1, z, mu1, trait.type)

  rslt <- list()
  if(all(z == 0)){
    # no covariate
    rslt$data <- cbind(dTrt, dGen1, dGen0)
    #colnames(rslt$data) = c("Y1","Y2", snpnames)
    colnames(rslt$data) = c("Y", snpnames)
  }else{
    rslt$data <- cbind(dTrt, z, dGen1, dGen0) # col1:ID, col2:Trait, other cols:SNPs
    #colnames(rslt$data) = c("Y1", "Y2", "Cov", snpnames)
    colnames(rslt$data) = c("Y", "Cov", snpnames)
  }

  a <- cbind(dGen1, dGen0)
  rslt$MAF <- (apply((a==1), 2, sum) + 2 * apply((a==2), 2, sum))/2/nrow(a)

  return (rslt)
}
simu_popu = cmpfun(simu_popu)

## penetrace function for 70 models ####
pen <- c ()
pen <- c (pen, 0.486, 0.960, 0.538, 0.947, 0.004, 0.811, 0.640, 0.606, 0.909 )
pen <- c (pen, 0.469, 0.956, 0.697, 0.945, 0.019, 0.585, 0.786, 0.407, 0.013 )
pen <- c (pen, 0.498, 0.954, 0.786, 0.978, 0.038, 0.428, 0.590, 0.821, 0.380 )
pen <- c (pen, 0.505, 0.988, 0.624, 0.945, 0.085, 0.807, 0.969, 0.116, 0.159 )
pen <- c (pen, 0.486, 0.963, 0.512, 0.941, 0.006, 0.899, 0.691, 0.541, 0.614 )
pen <- c (pen, 0.428, 0.757, 0.812, 0.788, 0.132, 0.044, 0.559, 0.548, 0.373 )
pen <- c (pen, 0.507, 0.842, 0.605, 0.845, 0.162, 0.629, 0.581, 0.678, 0.729 )
pen <- c (pen, 0.577, 0.247, 0.428, 0.227, 0.928, 0.578, 0.586, 0.262, 0.158 )
pen <- c (pen, 0.340, 0.637, 0.654, 0.689, 0.017, 0.041, 0.242, 0.866, 0.403 )
pen <- c (pen, 0.387, 0.726, 0.734, 0.749, 0.090, 0.034, 0.551, 0.401, 0.724 )
pen <- c (pen, 0.463, 0.703, 0.431, 0.653, 0.277, 0.806, 0.830, 0.008, 0.129 )
pen <- c (pen, 0.319, 0.507, 0.569, 0.553, 0.105, 0.045, 0.203, 0.777, 0.280 )
pen <- c (pen, 0.627, 0.393, 0.335, 0.396, 0.779, 0.953, 0.314, 0.997, 0.530 )
pen <- c (pen, 0.297, 0.540, 0.441, 0.541, 0.072, 0.278, 0.434, 0.293, 0.228 )
pen <- c (pen, 0.332, 0.562, 0.573, 0.583, 0.112, 0.147, 0.399, 0.496, 0.033 )
pen <- c (pen, 0.492, 0.664, 0.481, 0.642, 0.330, 0.746, 0.656, 0.396, 0.000 )
pen <- c (pen, 0.499, 0.639, 0.765, 0.666, 0.389, 0.083, 0.543, 0.527, 0.953 )
pen <- c (pen, 0.212, 0.350, 0.116, 0.336, 0.054, 0.495, 0.227, 0.273, 0.495 )
pen <- c (pen, 0.805, 0.683, 0.638, 0.657, 0.936, 0.989, 0.850, 0.564, 0.866 )
pen <- c (pen, 0.638, 0.488, 0.383, 0.464, 0.765, 0.957, 0.580, 0.562, 0.719 )
pen <- c (pen, 0.500, 0.926, 0.615, 0.895, 0.131, 0.647, 0.858, 0.160, 0.999 )
pen <- c (pen, 0.413, 0.851, 0.535, 0.831, 0.008, 0.580, 0.692, 0.268, 0.736 )
pen <- c (pen, 0.455, 0.848, 0.897, 0.890, 0.088, 0.016, 0.562, 0.686, 0.467 )
pen <- c (pen, 0.609, 0.980, 0.980, 0.993, 0.300, 0.275, 0.876, 0.483, 0.683 )
pen <- c (pen, 0.446, 0.844, 0.774, 0.879, 0.044, 0.233, 0.492, 0.796, 0.410 )
pen <- c (pen, 0.077, 0.656, 0.880, 0.892, 0.235, 0.312, 0.174, 0.842, 0.106 )
pen <- c (pen, 0.895, 0.323, 0.161, 0.068, 0.728, 0.806, 0.925, 0.233, 0.362 )
pen <- c (pen, 0.805, 0.251, 0.085, 0.002, 0.668, 0.638, 0.830, 0.079, 0.542 )
pen <- c (pen, 0.307, 0.682, 0.958, 0.997, 0.390, 0.281, 0.012, 0.990, 0.698 )
pen <- c (pen, 0.083, 0.891, 0.037, 0.619, 0.271, 0.691, 0.853, 0.079, 0.742 )
pen <- c (pen, 0.356, 0.891, 0.809, 0.955, 0.508, 0.611, 0.617, 0.755, 0.630 )
pen <- c (pen, 0.086, 0.536, 0.641, 0.677, 0.275, 0.096, 0.219, 0.413, 0.712 )
pen <- c (pen, 0.855, 0.339, 0.772, 0.513, 0.651, 0.607, 0.250, 0.999, 0.154 )
pen <- c (pen, 0.506, 0.838, 0.024, 0.603, 0.454, 0.957, 0.729, 0.427, 0.753 )
pen <- c (pen, 0.393, 0.764, 0.664, 0.850, 0.398, 0.733, 0.406, 0.927, 0.147 )
pen <- c (pen, 0.137, 0.484, 0.187, 0.482, 0.166, 0.365, 0.193, 0.361, 0.430 )
pen <- c (pen, 0.469, 0.198, 0.754, 0.337, 0.502, 0.141, 0.339, 0.453, 0.285 )
pen <- c (pen, 0.478, 0.311, 0.864, 0.387, 0.579, 0.263, 0.634, 0.436, 0.138 )
pen <- c (pen, 0.068, 0.299, 0.017, 0.289, 0.044, 0.285, 0.048, 0.262, 0.174 )
pen <- c (pen, 0.539, 0.120, 0.258, 0.165, 0.378, 0.325, 0.123, 0.426, 0.276 )
pen <- c (pen, 0.002, 0.155, 0.214, 0.199, 0.071, 0.022, 0.081, 0.122, 0.135 )
pen <- c (pen, 0.188, 0.020, 0.171, 0.032, 0.174, 0.059, 0.134, 0.087, 0.092 )
pen <- c (pen, 0.005, 0.179, 0.251, 0.211, 0.100, 0.026, 0.156, 0.098, 0.156 )
pen <- c (pen, 0.174, 0.321, 0.154, 0.223, 0.254, 0.245, 0.448, 0.025, 0.424 )
pen <- c (pen, 0.098, 0.219, 0.302, 0.302, 0.126, 0.121, 0.053, 0.308, 0.136 )
pen <- c (pen, 0.891, 0.362, 0.480, 0.213, 0.829, 0.601, 0.925, 0.267, 0.685 )
pen <- c (pen, 0.077, 0.689, 0.417, 0.763, 0.150, 0.491, 0.196, 0.657, 0.247 )
pen <- c (pen, 0.132, 0.793, 0.274, 0.799, 0.213, 0.514, 0.255, 0.528, 0.793 )
pen <- c (pen, 0.611, 0.104, 0.759, 0.180, 0.674, 0.019, 0.532, 0.189, 0.681 )
pen <- c (pen, 0.091, 0.827, 0.863, 0.869, 0.393, 0.415, 0.738, 0.508, 0.363 )
pen <- c (pen, 0.495, 0.415, 0.657, 0.429, 0.616, 0.121, 0.552, 0.331, 0.419 )
pen <- c (pen, 0.592, 0.691, 0.743, 0.712, 0.493, 0.419, 0.580, 0.746, 0.504 )
pen <- c (pen, 0.108, 0.194, 0.186, 0.196, 0.037, 0.045, 0.172, 0.073, 0.130 )
pen <- c (pen, 0.112, 0.186, 0.128, 0.193, 0.024, 0.138, 0.079, 0.236, 0.251 )
pen <- c (pen, 0.272, 0.192, 0.185, 0.172, 0.367, 0.390, 0.345, 0.069, 0.005 )
pen <- c (pen, 0.247, 0.301, 0.205, 0.300, 0.173, 0.378, 0.215, 0.357, 0.268 )
pen <- c (pen, 0.222, 0.276, 0.141, 0.259, 0.169, 0.401, 0.278, 0.128, 0.420 )
pen <- c (pen, 0.260, 0.221, 0.201, 0.204, 0.315, 0.348, 0.339, 0.074, 0.128 )
pen <- c (pen, 0.139, 0.188, 0.221, 0.190, 0.111, 0.020, 0.206, 0.051, 0.253 )
pen <- c (pen, 0.558, 0.616, 0.674, 0.632, 0.500, 0.418, 0.546, 0.674, 0.395 )
pen <- c (pen, 0.166, 0.165, 0.128, 0.114, 0.199, 0.143, 0.281, 0.028, 0.281 )
pen <- c (pen, 0.108, 0.006, 0.080, 0.026, 0.079, 0.046, 0.021, 0.090, 0.025 )
pen <- c (pen, 0.006, 0.094, 0.008, 0.079, 0.016, 0.076, 0.052, 0.043, 0.057 )
pen <- c (pen, 0.199, 0.072, 0.168, 0.086, 0.187, 0.076, 0.125, 0.108, 0.226 )
pen <- c (pen, 0.165, 0.096, 0.262, 0.166, 0.151, 0.091, 0.050, 0.250, 0.056 )
pen <- c (pen, 0.103, 0.063, 0.124, 0.098, 0.086, 0.069, 0.021, 0.147, 0.059 )
pen <- c (pen, 0.185, 0.291, 0.234, 0.286, 0.201, 0.277, 0.249, 0.266, 0.166 )
pen <- c (pen, 0.073, 0.042, 0.015, 0.024, 0.064, 0.059, 0.068, 0.019, 0.095 )
pen <- c (pen, 0.046, 0.127, 0.069, 0.115, 0.067, 0.097, 0.107, 0.069, 0.108 )
pen <- c (pen, 0.095, 0.122, 0.127, 0.097, 0.129, 0.100, 0.201, 0.044, 0.122 )

MAF <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
         0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
         0.4, 0.4, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4,
         0.4, 0.4, 0.4, 0.4, 0.4, 0.4)

dim(pen) = c(9, 70)


# ================================================================================= #


#  An: Example for first 2-locus interaction causal model    ####
# ================================================================================= #


#R CMD BATCH '--args arg_1 arg_2 arg_3 arg_4 arg_5 arg_6' ex_script.R result.Rout
#args = (commandArgs(TRUE))
#model.no = as.numeric(args[1])
#mu1 = as.numeric(args[2])
#rho = as.numeric(args[3])
#k = as.numeric(args[4])
#seed = as.numeric(args[5]) + model.no * 1000
#outfile = args[6]


#M1 = 2           # Number of causal SNPs with fixed MAFs, fisrt M1 as causeal
#M0 = 18          # Number of non-causal SNPs

#cov_mean = 0
#cov_var = 0

#pen1 <- pen[, model.no]   ## input model.no

#maf = MAF[model.no]
#MAF1 <- rep(maf, M1)
#MAF0 <- maf












