# looped demo 
library(ape) # for phylogenetic data
library(mvtnorm) # to compute multivariate probabilities
#library(cmdstanr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


nom = "demo_rho_"
model_name = "jsdm_stan_demo_noloop" #"poisson_binomial_pobs" # 
#model_name = "poisson_binomial_pobs"

mod <- stan_model(file = paste(model_name,  '.stan', sep=""), 
                  verbose = FALSE, 
                  auto_write = rstan_options("auto_write" = TRUE))


tmp = matrix(NA, 1, 8)
colnames(tmp) = c("par", "replicate", "mean", "2.5%", "97.5%", "n_eff", "Rhat", "true")
write.table(tmp, file = paste("res_", nom,"_", model_name, ".csv", sep = ""), 
            sep = ",", append = TRUE, quote = FALSE)

repeats <- 50

n_sp  = 20  # number of species
n_env = 2   # environmental covariates 
n_t   = 2   # number of species specific traits 

sd_z  = 0.5 # sd for the Z matrix (params for traits)
sd_s  = 0.3 # sd for the Sigma matrix (covariance among parameters)

n_sites = 100
n_pars= n_env+1 #parameters (environmental cov + intercept)


for (k in 1:repeats){
  
  rho   = runif(1, 0, 1) # strength of phylogenetic signal
  
  
  # defining species specific  traits 
  dgrass = rbeta(n_sp, 1, 1) # fraction of grass in diet 
  log_bm = rnorm(n_sp, 0, 1) # log of body mass (centered)
  
  #arrange species specific traits in a matrix + intercept
  TT = as.matrix(cbind(rep(1, n_sp), scale(dgrass), scale(log_bm)))
  
  # simulate phylogeny 
  tree = rtree(n_sp)
  CC = vcv(tree, corr=T) 
  
  # sort species and re-arrange phylogenetic correlation matrix
  tmp = dimnames(CC)
  ids = as.numeric(as.factor(tmp[[1]]))
  
  C = matrix(NA, ncol(CC), ncol(CC))
  for(i in 1:ncol(CC)){
    for(j in 1:ncol(CC)){
      C[ids[i],ids[j]] = CC[i,j]
    }
  }
  
  #Z = matrix(rnorm((n_t + 1) * n_pars, 0 , sd_z),  (n_t + 1), n_pars)
  Z = matrix(runif((n_t + 1) * n_pars, -0.6 , 0.6),  (n_t + 1), n_pars)
  M = TT %*% Z
  
  Sigma = diag(n_pars) * sd_s
  
  betas = rmvnorm(1, mean=as.vector(M), kronecker(Sigma, rho*C + (1-rho)*diag(n_sp)))
  Beta = matrix(betas[1,],n_sp,n_pars)
  
  X = cbind(rep(1, n_sites), matrix(rnorm(n_sites*n_env), n_sites, n_env))
  
  
  
  # simulate real abundances
  N = matrix(NA, n_sites, n_sp)
  for(i in 1: n_sp){
    N[,i] = rpois(n_sites, lambda = exp(X %*% Beta[i, ]))
  }
  
  p = rbeta(n_sp, 5, 3) # detection probability
  
  
  data = NULL
  for (a in 1:n_sites) {
    for(b in 1:n_sp){
      if(N[a,b] > 0){
        #d = runif(N[i,j], 0, B)
        #gs = rpois(N[i,j], lambda.group[j]) + 1 # group size
        #sigma = exp(Beta[j,1] + gs * Beta[j,2] + TT[j,2] * Beta[j,3])
        #p = exp(-d * d/(2 * (sigma^2)))
        y = rbinom(N[a,b], 1, p[b])
        #d = d[y == 1]
        #gs = gs[y == 1]
        y = y[y == 1]
        if (sum(y) > 0){
          data = rbind(data, cbind(rep(a, sum(y)), rep(b, sum(y)), y))
        }
      }
    }
  }
  
  colnames(data) = c("site","sp", "y")
  datos = as.data.frame(data)
  
  
  
  stan_dat <- list(
    n_obs = nrow(datos),
    area = rep(1.0, n_sites),
    n_sites = n_sites,
    site = as.integer(datos$site), # (c(s, rep(1:n_sites, each = nzs ))),
    K = ncol(X), 
    X = X,
    n_max = rep(100, n_sp),
    n_s = as.integer(n_sp),
    n_t = dim(TT)[2],
    TT = TT,
    C = C,
    ones = numeric(n_sp) + 1,
    sp = datos$sp,
    p_obs = p
  )
  
  
  fit <- sampling(mod,
                  data = stan_dat,
                  iter = 1000, 
                  thin = 1,
                  chains = 3)
  
  fit_summary <- summary(fit)$summary
  
  rs = as.matrix(fit_summary[grepl("rho", rownames(fit_summary)),c(1,4,8,9,10)])
  colnames(rs) = "rho"
  bs = as.matrix(fit_summary[grepl("betas", rownames(fit_summary)),c(1,4,8,9,10)])
  zs = as.matrix(fit_summary[grepl("z", rownames(fit_summary)),c(1,4,8,9,10)])
  
  tmp = rbind(t(rs), zs,bs)
  trues = rbind(rho, matrix(c(Z),length(Z),1), t(betas)) 
  res = cbind(matrix(rep(k, nrow(tmp)), nrow(tmp),1), tmp, trues)
  
  write.table(res, file = paste("res_", nom,"_", model_name, ".csv", sep = ""), 
              sep = ",", append = TRUE, quote = FALSE, col.names = FALSE) #, row.names = FALSE)
  
}
