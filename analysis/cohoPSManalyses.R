options(device=windows)
library(rstan)
library(loo)
library(shinystan)
library(Hmisc)
library(corrplot)

# Convenience functions to simplify extracting posterior means
# or single parameters from stanfit objects
stan_mean <- function(object, pars)
{
  mm <- get_posterior_mean(object, pars)
  return(mm[,ncol(mm)])
}

extract1 <- function(object, par)
{
  extract(object, par)[[1]]
}


#==================================================================
# DATA
#==================================================================

#----------------------------------------------
# PSM SAMPLING SITES FOR ANALYSIS
#----------------------------------------------

# read in spawner data
spawner.data <- read.csv("spawner.data.csv", header=T)

# read in landscape data
spatial.data <- read.csv("spatial.data.csv", header=T)

# LU/LC and roads data
lulc.roads.data <- spatial.data[,c("site","area","ccap.decid","ccap.open","ccap.evgrn","ccap.hdev","ccap.ldev",
                                   "ccap.mdev","ccap.mxforest","ccap.wetland","ccap.ag","ccap.other",
                                   "roads.1","roads.2","roads.3","roads.4","roads.5","traffic","nlcd.imperv",
                                   "restoration","pop.census","pop.lscan")]

# convert basin area from m2 to ha
lulc.roads.data$area <- lulc.roads.data$area/10e4

# convert to km of road/km2
lulc.roads.data[,substring(names(lulc.roads.data),1,5)=="roads"] <- lulc.roads.data[,substring(names(lulc.roads.data),1,5)=="roads"]/1000

# Pretty variable names for LU/LC and roads
lulc.roads.labels <- read.csv("lulc.roads.labels.csv", header=T)
lulc.roads.labels$plot.label <- ordered(as.character(lulc.roads.labels$plot.label),
                                        levels = c("Imperviousness","High developed","Medium developed",
                                                   "Low developed","Evergreen","Deciduous","Mixed forest",
                                                   "Open space","Wetland","Agriculture","Interstate",
                                                   "Principal arterial","Minor arterial","Collector arterial",
                                                   "Local roads","Traffic intensity","Restoration",
                                                   "U.S. Census","LandScan"))

# matrix of annual summer and fall precip with sites and years as rows
ppt.data <- spatial.data[,substring(names(spatial.data),1,8) %in% c("site","ppt.su.2","ppt.fa.2")]
ppt.su.cols <- substring(names(ppt.data),1,8)=="ppt.su.2"
ppt.fa.cols <- substring(names(ppt.data),1,8)=="ppt.fa.2"
n.yrs.ppt <- sum(ppt.su.cols)

ppt.data <- data.frame(site=rep(ppt.data$site, each=n.yrs.ppt), 
                       year=rep(substring(names(ppt.data)[ppt.su.cols], 8), nrow(ppt.data)),
                       ppt.su=as.vector(t(as.matrix(ppt.data[,ppt.su.cols]))),
                       ppt.fa=as.vector(t(as.matrix(ppt.data[,ppt.fa.cols]))))

# merge basin-specific variables and annual su and fa ppt into spawner data
psm <- merge(spawner.data, ppt.data, by=c("site","year"))
spatial.data.cols <- substring(names(spatial.data),1,3) != "ppt"
psm <- merge(psm, spatial.data[,spatial.data.cols], by=c("ID","site"))
psm <- data.frame(psm[,c("ID","site","watershed","year","n","n.psm","ppt.su","ppt.fa","area")], psm[,-(1:10)])
psm <- psm[order(psm$site,psm$year),]

rm(list=c("ppt.su.cols","ppt.fa.cols","n.yrs.ppt","spatial.data.cols"))

#----------------------------------------------
# ADDITIONAL SITES FOR PSM PREDICTIONS
#----------------------------------------------

# read in landscape data
spatial.data.pre <- read.csv("spatial.data.predict.csv", header=T)
names(spatial.data.pre)[names(spatial.data.pre)=="ID"] <- "site"

# convert basin area from m2 to ha
spatial.data.pre$area <- spatial.data.pre$area/10e4

# convert to km of road/km2
spatial.data.pre[,substring(names(spatial.data.pre),1,5)=="roads"] <- spatial.data.pre[,substring(names(spatial.data.pre),1,5)=="roads"]/1000

spatial.data.pre <- spatial.data.pre[spatial.data.pre$coho==1,] # only include basins w/ coho 
spatial.data.pre$site <- factor(spatial.data.pre$site)

lulc.roads.data.pre <- spatial.data.pre[,names(lulc.roads.data)]

# Combine PSM-sampling and unsampled-site data frames
lulc.roads.data.all <- rbind(data.frame(data="psm", lulc.roads.data), 
                             data.frame(data="pre", lulc.roads.data.pre))

# Add rows to psm for unobserved sites, plus one extra row for each observed site
# containing no PSM data and having ppt set to the mean. These extra rows are 
# used to generate predictions for monitored sites under baseline precip conditions.
psm.all <- as.data.frame(matrix(nrow=nrow(spatial.data.pre) + nrow(spatial.data), ncol=ncol(psm)))
names(psm.all) <- names(psm)
psm.all$site <- c(as.character(lulc.roads.data$site), as.character(lulc.roads.data.pre$site))
psm.all$n <- 0
psm.all$n.psm <- 0
psm.all$ppt.su <- mean(psm$ppt.su)
psm.all$ppt.fa <- mean(psm$ppt.fa)
psm.all <- rbind(data.frame(data="psm", psm), 
                 data.frame(data="pre", psm.all))


#==================================================================
# STRUCTURAL EQUATION MODELS FOR LANDSCAPE DATA AND PSM
#==================================================================

#------------------------------------------------------
# Fit full model to sites with PSM observations
#------------------------------------------------------

# site-level covariates

# composition data (multiple categories):
# bound away from 0 and 1 and re-standardize to sum to 1,
# then transform to log ratios
X1 <- as.matrix(lulc.roads.data[,c("ccap.hdev","ccap.ldev","ccap.mdev",
                                   "ccap.decid","ccap.evgrn","ccap.mxforest",
                                   "ccap.open","ccap.wetland","ccap.ag", "ccap.other")])
X1[X1==0] <- 1e-4
X1[X1==1] <- 1 - 1e-4
X1 <- sweep(X1, 1, rowSums(X1), "/")
X1 <- sweep(log(X1[,-ncol(X1)]), 1, log(X1[,ncol(X1)]), "-") 
X1 <- sweep(sweep(X1, 2, colMeans(X1), "-"), 2, apply(X1, 2, sd), "/")

# proportion data:
# bound away from 0 and 1 and logit-transform
X2 <- as.matrix(lulc.roads.data[,"nlcd.imperv",drop=F])
X2[X2==0] <- 1e-4
X2[X2==1] <- 1 - 1e-4
X2 <- qlogis(X2)
X2 <- sweep(sweep(X2, 2, colMeans(X2), "-"), 2, apply(X2, 2, sd), "/")

# nonnegative continuous data:
# scale to SD = 1, bound away from 0
X3 <- as.matrix(lulc.roads.data[,c("roads.1","roads.2","roads.3","roads.4","roads.5","traffic",
                                   "restoration","pop.census","pop.lscan")])
X3 <- sweep(X3, 2, apply(X3, 2, sd), "/")
X3[X3==0] <- 1e-4
X3 <- sweep(X3, 2, apply(X3, 2, sd), "/")

# reassemble all variables into matrix
X <- cbind(X2,X1,X3)
rm(list = c("X1","X2","X3"))

# indices of covariates assumed to be normally (gamma) distributed
normal.indx <- which(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd"))
gamma.indx <- which(!(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd")))

# Data for Stan
stan.dat <- list(S = nrow(X), 
                 D_normal = length(normal.indx), D_gamma = length(gamma.indx),
                 X = X, 
                 L = 1,  # user-specified!
                 N = nrow(psm), 
                 site = as.numeric(psm$site),
                 ppt_su = as.vector(scale(psm$ppt.su/10, scale=F)),
                 ppt_fa = as.vector(scale(psm$ppt.fa/10, scale=F)),
                 I0_Z = 1,
                 I_su = 1,
                 I_su_Z = 1,
                 I_fa = 1,
                 I_fa_Z = 1,
                 n = psm$n,
                 n_psm = psm$n.psm,
                 I_fit = rep(1, nrow(psm)),
                 I_lpd = rep(1, nrow(psm)))


# Function to generate initial values for chains
stan.init <- function() 
{
  if(exists("stan.dat"))
    for(i in names(stan.dat)) assign(i, stan.dat[[i]])
  # X <- t(X)
  D <- D_normal + D_gamma
  
  list(a0 = rnorm(D, c(colMeans(X[,1:D_normal]), colMeans(log(X[,-(1:D_normal)]))), 1),
       A_nid_vec = array(rnorm(D*L - L*(L-1)/2, 0, 1), dim = D*L - L*(L-1)/2),
       Z_nid = matrix(rnorm(S*L, 0, 1), nrow = S, ncol = L),
       phi = runif(D, 0.5, 1),
       # mu_b0 = rnorm(1, qlogis(mean(n_psm/n, na.rm=T)), 1),
       mu_b0 = rnorm(1, -3, 1),
       b0_Z_nid = array(rnorm(L, 1, 0.5), dim = L),
       sigma_b0 = runif(1, 1, 2),
       b0_std = array(rnorm(S, 0, 0.1), dim = S),
       mu_b_su = rnorm(1, 0.3, 0.1),
       b_su_Z_nid = array(rnorm(L, -0.2, 0.1), dim = L),
       sigma_b_su = runif(1, 0.1, 1),
       b_su_std = array(rnorm(S, 0, 0.1), dim = S),
       mu_b_fa = rnorm(1, 0.03, 0.02),
       b_fa_Z_nid = array(rnorm(L, -0.03, 0.02), dim = L),
       sigma_b_fa = runif(1, 0.02, 0.1),
       b_fa_std = array(rnorm(S, 0, 0.1), dim = S),
       sigma_psm = runif(1, 0.5, 2),
       logit_p_psm_std = array(rnorm(N, 0, 0.1), dim = N))
}

# Fit it!
stan.psm <- stan(file = "cohoPSM_SEM_Stan.stan",
                 data = stan.dat, 
                 init = stan.init,
                 pars = c("a0","A","Z","phi","g_mu_X",
                          "mu_b0","b0_Z","sigma_b0","b0",
                          "mu_b_su","b_su_Z","sigma_b_su","b_su",
                          "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                          "sigma_psm","p_psm","ll_psm"), 
                 chains = 3, iter = 50000, warmup = 10000, thin = 40, cores = 3,
                 control = list(adapt_delta = 0.8))

# Inspect summary and a few traceplots and bivariate posteriors
print(stan.psm, pars = c("g_mu_X","ll_psm"), include = F)
traceplot(stan.psm, pars = c("A","lp__"), inc_warmup = F)
traceplot(stan.psm, pars = c("Z"), inc_warmup = F)
pairs(stan.psm, pars = c("a0[1]", "a0[2]", "A[1,1]", "A[2,1]", "Z[1,1]"), condition = "accept_stat__")
pairs(stan.psm, pars = paste("A[", normal.indx, ",1]", sep=""), labels = colnames(X)[normal.indx])
pairs(stan.psm, pars = paste("A[", gamma.indx, ",1]", sep=""), labels = colnames(X)[gamma.indx])


#------------------------------------------------------------------
# In-sample model selection:
# Fit all candidate models to sites with PSM observations,
# compare expected out-of-sample performance via WAIC and IS-LOO
#------------------------------------------------------------------

# Create list to store results
stan.psm.list <- vector("list", 18)
k <- 1
for(b0_mod in c(1,"Z"))
  for(b_su_mod in c(0,1,"Z"))
    for(b_fa_mod in c(0,1,"Z"))
    {
      names(stan.psm.list)[k] <- paste(b0_mod, b_su_mod, b_fa_mod, sep="")
      k <- k + 1
    }

# Loop over candidate models, fit to data, and store stan objects
# (Note that this assumes stan.dat and stan.init have been assigned in
# the previous code block)
for(i in 1:length(stan.psm.list))
{
  # Assign binary in/out indicators
  stan.dat$I0_Z <- as.numeric(substring(names(stan.psm.list)[i],1,1) == "Z")
  stan.dat$I_su <- as.numeric(substring(names(stan.psm.list)[i],2,2) != "0")
  stan.dat$I_su_Z <- as.numeric(substring(names(stan.psm.list)[i],2,2) == "Z")
  stan.dat$I_fa <- as.numeric(substring(names(stan.psm.list)[i],3,3) != "0")
  stan.dat$I_fa_Z <- as.numeric(substring(names(stan.psm.list)[i],3,3) == "Z")
  
  # Fit model
  cat("Working on model", i, "(see Viewer for progress) \n")
  fit <- stan(file = "cohoPSM_SEM_Stan.stan",
              data = stan.dat, 
              init = stan.init,
              pars = c("a0","A","Z","phi",
                       "mu_b0","b0_Z","sigma_b0","b0",
                       "mu_b_su","b_su_Z","sigma_b_su","b_su",
                       "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                       "sigma_psm","p_psm","ll_psm"), 
              chains = 3, iter = 50000, warmup = 10000, thin = 40, cores = 3)
  
  # Store fitted object
  stan.psm.list[[i]] <- list(fit = fit, 
                             WAIC = waic(extract1(fit,"ll_psm")), 
                             LOO = loo(extract1(fit,"ll_psm")))
}

# Summarize results in data frame
stan.psm.mods <- data.frame(model = names(stan.psm.list),
                            Dbar = sapply(stan.psm.list, function(x) -2*sum(stan_mean(x$fit,"ll_psm"))),
                            p_WAIC = sapply(stan.psm.list, function(x) x$WAIC$p_waic),
                            WAIC = sapply(stan.psm.list, function(x) x$WAIC$waic),
                            dWAIC = NA,
                            se_dWAIC = NA,
                            p_LOO = sapply(stan.psm.list, function(x) x$LOO$p_loo),
                            LOO = sapply(stan.psm.list, function(x) x$LOO$looic),
                            dLOO = NA,
                            se_dLOO = NA)

for(i in 1:length(stan.psm.list))
{
  compare_waic <- compare(stan.psm.list[[i]]$WAIC, 
                          stan.psm.list[[which.min(stan.psm.mods$WAIC)]]$WAIC)
  stan.psm.mods$dWAIC[i] <- 2*compare_waic["elpd_diff"]
  stan.psm.mods$se_dWAIC[i] <- 2*compare_waic["se"]
  compare_loo <- compare(stan.psm.list[[i]]$LOO, 
                         stan.psm.list[[which.min(stan.psm.mods$LOO)]]$LOO)
  stan.psm.mods$dLOO[i] <- 2*compare_loo["elpd_diff"]
  stan.psm.mods$se_dLOO[i] <- 2*compare_loo["se"]
  rm(compare_waic); rm(compare_loo)
}

stan.psm.mods <- stan.psm.mods[order(stan.psm.mods$dLOO),]                            
stan.psm.mods[order(stan.psm.mods$dWAIC),]


#---------------------------------------------------------------------
# K-fold cross-validation over YEARS:
# Leave out one year of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations
#---------------------------------------------------------------------

# Shortlist of candidate models chosen by consensus between in-sample WAIC and LOO
stan.psm.cv.year.list <- vector("list", 4)
names(stan.psm.cv.year.list) <- c("100","1Z0","ZZ1","ZZZ")

# Fit "full" model to complete data once to generate good starting values
# for leave-one-out runs
fit1 <- stan(file = "cohoPSM_SEM_Stan.stan",
             data = stan.dat, 
             init = stan.init,
             pars = names(stan.init()), 
             chains = 3, iter = 10000, warmup = 5000, thin = 1, cores = 3,
             control = list(adapt_delta = 0.8))

stan.init.cv <- function()
{
  samples <- extract(fit1)
  M <- nrow(samples$a0)
  i <- sample(M,1)
  
  list(a0 = samples$a0[i,],
       A_nid_vec = array(samples$A_nid_vec[i,], dim = ncol(samples$A_nid_vec)),
       Z_nid = matrix(samples$Z[i,,], nrow = dim(samples$Z)[2], ncol = dim(samples$Z)[3]),
       phi = samples$phi[i,],
       mu_b0 = samples$mu_b0[i],
       b0_Z_nid = array(samples$b0_Z_nid[i,], dim = ncol(samples$b0_Z_nid)),
       sigma_b0 = samples$sigma_b0[i],
       b0_std = array(samples$b0_std[i,], dim = ncol(samples$b0_std)),
       mu_b_su = samples$mu_b_su[i],
       b_su_Z_nid = array(samples$b_su_Z_nid[i,], dim = ncol(samples$b_su_Z_nid)),
       sigma_b_su = samples$sigma_b_su[i],
       b_su_std = array(samples$b_su_std[i,], dim = ncol(samples$b_su_std)),
       mu_b_fa = samples$mu_b_fa[i],
       b_fa_Z_nid = array(samples$b_fa_Z_nid[i,], dim = ncol(samples$b_fa_Z_nid)),
       sigma_b_fa = samples$sigma_b_fa[i],
       b_fa_std = array(samples$b_fa_std[i,], dim = ncol(samples$b_fa_std)),
       sigma_psm = samples$sigma_psm[i],
       logit_p_psm_std = array(samples$logit_p_psm_std[i,], dim = ncol(samples$logit_p_psm_std)))
}

# Loop over candidate models and then over years, fit to data, predict hold-out data,
# and store stan objects
# (Note that this assumes stan.dat was assigned in the previous code block)
for(i in 1:length(stan.psm.cv.year.list))
{
  stan.psm.cv.year.list[[i]] <- list(fit = vector("list",length(unique(psm$year))), 
                                     ll_psm = matrix(NA, 3000, nrow(psm)),
                                     elpd = NULL)
  names(stan.psm.cv.year.list[[i]]$fit) <- sort(unique(psm$year))
  stan.dat.cv.year <- stan.dat
  
  # Assign binary in/out indicators
  stan.dat.cv.year$I0_Z <- as.numeric(substring(names(stan.psm.cv.year.list)[i],1,1) == "Z")
  stan.dat.cv.year$I_su <- as.numeric(substring(names(stan.psm.cv.year.list)[i],2,2) != "0")
  stan.dat.cv.year$I_su_Z <- as.numeric(substring(names(stan.psm.cv.year.list)[i],2,2) == "Z")
  stan.dat.cv.year$I_fa <- as.numeric(substring(names(stan.psm.cv.year.list)[i],3,3) != "0")
  stan.dat.cv.year$I_fa_Z <- as.numeric(substring(names(stan.psm.cv.year.list)[i],3,3) == "Z")
  
  for(j in sort(unique(psm$year)))
  {
    stan.dat.cv.year$I_fit <- as.numeric(psm$year != j)  # training data
    stan.dat.cv.year$I_lpd <- as.numeric(psm$year == j)  # hold-out data
    
    # Fit model
    cat("Working on model", i, "and hold-out year", j, "(see Viewer for progress) \n")
    fit <- stan(file = "cohoPSM_SEM_Stan.stan",
                data = stan.dat.cv.year, 
                init = lapply(1:3,function(x) stan.init.cv()),
                pars = c("a0","A","Z","phi",
                         "mu_b0","b0_Z","sigma_b0","b0",
                         "mu_b_su","b_su_Z","sigma_b_su","b_su",
                         "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                         "sigma_psm","p_psm","ll_psm"), 
                chains = 3, iter = 50000, warmup = 10000, thin = 40, cores = 3)
    
    # Store fitted object and log-likelihood matrix
    stan.psm.cv.year.list[[i]]$fit[[as.character(j)]] <- fit
    stan.psm.cv.year.list[[i]]$ll_psm[,psm$year == j] <- matrix(extract1(fit,"ll_psm"), nrow = 3000)
  }
  
  # Average the posterior predictive density (not log density!) across posterior samples,
  # take the log, then sum across observations to get expected log predictive density
  stan.psm.cv.year.list[[i]]$elpd <- -2*sum(log(colMeans(exp(stan.psm.cv.year.list[[i]]$ll_psm))))
}

# Summarize model comparison results
stan.psm.cv.year.mods <- data.frame(model = names(stan.psm.cv.year.list),
                                    elpd = sapply(stan.psm.cv.year.list,
                                                  function(x) x$elpd),
                                    se_elpd = sapply(stan.psm.cv.year.list, function(x) 
                                      2*sqrt(ncol(x$ll_psm))*sd(log(colMeans(exp(x$ll_psm))))),
                                    d_elpd = NA,
                                    se_d_elpd = NA)

for(i in 1:length(stan.psm.cv.year.list))
{
  stan.psm.cv.year.mods$d_elpd[i] <- stan.psm.cv.year.mods$elpd[i] - min(stan.psm.cv.year.mods$elpd)
  m <- which.min(stan.psm.cv.year.mods$elpd)
  devs <- -2*(log(colMeans(exp(stan.psm.cv.year.list[[i]]$ll_psm))) - 
                log(colMeans(exp(stan.psm.cv.year.list[[m]]$ll_psm))))
  stan.psm.cv.year.mods$se_d_elpd[i] <- sqrt(length(devs))*sd(devs)
  rm(m);rm(devs)
}

stan.psm.cv.year.mods


#---------------------------------------------------------------------
# K-fold cross-validation over SITES:
# Leave out one or more sites of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations.
# Sites are randomly partitioned into K = 10 groups that are
# roughly similar in size (i.e., number of observations).
#---------------------------------------------------------------------

# Randomly partition sites into groups
KfoldCV.partition <- function(psm.dat, K, N.random = 1000)
{
  grps <- matrix(NA, length(levels(psm$site)), N.random)
  pos <- 1:nrow(grps)
  target <- nrow(psm)/K
  for(j in 1:ncol(grps))
  {
    pos.permuted <- sample(nrow(grps), nrow(grps), replace = FALSE)
    N.site <- table(psm$site)[pos.permuted]
    grp.permuted <- c(1, rep(0, nrow(grps) - 1))
    for(i in 2:nrow(grps))
    {
      if(abs(sum(N.site[grp.permuted==grp.permuted[i-1]]) + N.site[i] - target) <
         abs(sum(N.site[grp.permuted==grp.permuted[i-1]]) - target))
      {
        grp.permuted[i] <- grp.permuted[i-1]
      } else {
        if(grp.permuted[i] < 10) grp.permuted[i] <- grp.permuted[i-1] + 1
      }
    }
    grps[,j] <- grp.permuted[order(pos.permuted)]
  }
  N.site <- table(psm$site)
  range.N.grp <- apply(grps, 2, function(grp) diff(range(tapply(N.site, grp, sum))))
  grp <- grps[,which.min(range.N.grp)]
  return(list(N.group = tapply(N.site, grp, sum), group = grp[match(psm$site, levels(psm$site))]))
}

partitions <- KfoldCV.partition(psm.dat = psm, K = 10)
partitions  # check that the procedure found a "good" partition (roughly equal group sizes)
site.group <- partitions$group

# Shortlist of candidate models chosen by consensus between in-sample WAIC and LOO
stan.psm.cv.site.list <- vector("list", 4)
names(stan.psm.cv.site.list) <- c("100","1Z0","ZZ1","ZZZ")

# Fit "full" model to complete data once to generate good starting values
# for leave-one-out runs
fit1 <- stan(file = "cohoPSM_SEM_Stan.stan",
             data = stan.dat, 
             init = stan.init,
             pars = names(stan.init()), 
             chains = 3, iter = 10000, warmup = 5000, thin = 1, cores = 3,
             control = list(adapt_delta = 0.8))

stan.init.cv <- function()
{
  samples <- extract(fit1)
  M <- nrow(samples$a0)
  i <- sample(M,1)
  
  list(a0 = samples$a0[i,],
       A_nid_vec = array(samples$A_nid_vec[i,], dim = ncol(samples$A_nid_vec)),
       Z_nid = matrix(samples$Z[i,,], nrow = dim(samples$Z)[2], ncol = dim(samples$Z)[3]),
       phi = samples$phi[i,],
       mu_b0 = samples$mu_b0[i],
       b0_Z_nid = array(samples$b0_Z_nid[i,], dim = ncol(samples$b0_Z_nid)),
       sigma_b0 = samples$sigma_b0[i],
       b0_std = array(samples$b0_std[i,], dim = ncol(samples$b0_std)),
       mu_b_su = samples$mu_b_su[i],
       b_su_Z_nid = array(samples$b_su_Z_nid[i,], dim = ncol(samples$b_su_Z_nid)),
       sigma_b_su = samples$sigma_b_su[i],
       b_su_std = array(samples$b_su_std[i,], dim = ncol(samples$b_su_std)),
       mu_b_fa = samples$mu_b_fa[i],
       b_fa_Z_nid = array(samples$b_fa_Z_nid[i,], dim = ncol(samples$b_fa_Z_nid)),
       sigma_b_fa = samples$sigma_b_fa[i],
       b_fa_std = array(samples$b_fa_std[i,], dim = ncol(samples$b_fa_std)),
       sigma_psm = samples$sigma_psm[i],
       logit_p_psm_std = array(samples$logit_p_psm_std[i,], dim = ncol(samples$logit_p_psm_std)))
}

# Loop over candidate models and then over groups of sites, fit to data, 
# predict hold-out data, and store stan objects
# (Note that this assumes stan.dat has been assigned in a previous code block)
for(i in 1:length(stan.psm.cv.site.list))
{
  stan.psm.cv.site.list[[i]] <- list(fit = vector("list",length(unique(site.group))),
                                     site.group = site.group,
                                     ll_psm = matrix(NA, 3000, nrow(psm)),
                                     elpd = NULL)
  names(stan.psm.cv.site.list[[i]]$fit) <- sort(unique(site.group))
  stan.dat.cv.site <- stan.dat
  
  # Assign binary in/out indicators
  stan.dat.cv.site$I0_Z <- as.numeric(substring(names(stan.psm.cv.site.list)[i],1,1) == "Z")
  stan.dat.cv.site$I_su <- as.numeric(substring(names(stan.psm.cv.site.list)[i],2,2) != "0")
  stan.dat.cv.site$I_su_Z <- as.numeric(substring(names(stan.psm.cv.site.list)[i],2,2) == "Z")
  stan.dat.cv.site$I_fa <- as.numeric(substring(names(stan.psm.cv.site.list)[i],3,3) != "0")
  stan.dat.cv.site$I_fa_Z <- as.numeric(substring(names(stan.psm.cv.site.list)[i],3,3) == "Z")
  
  for(j in sort(unique(site.group)))
  {
    stan.dat.cv.site$I_fit <- as.numeric(site.group != j)  # training data
    stan.dat.cv.site$I_lpd <- as.numeric(site.group == j)  # hold-out data
    
    # Fit model
    cat("Working on model", i, "and hold-out group", j, "(see Viewer for progress) \n")
    fit <- stan(file = "cohoPSM_SEM_Stan.stan",
                data = stan.dat.cv.site, 
                init = lapply(1:3,function(x) stan.init.cv()),
                pars = c("a0","A","Z","phi",
                         "mu_b0","b0_Z","sigma_b0","b0",
                         "mu_b_su","b_su_Z","sigma_b_su","b_su",
                         "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                         "sigma_psm","p_psm","ll_psm"), 
                chains = 3, iter = 50000, warmup = 10000, thin = 40, cores = 3)
    
    # Store fitted object and log-likelihood matrix
    stan.psm.cv.site.list[[i]]$fit[[as.character(j)]] <- fit
    stan.psm.cv.site.list[[i]]$ll_psm[,site.group == j] <- matrix(extract1(fit,"ll_psm"), nrow = 3000)
  }
  
  # Average the posterior predictive density (not log density!) across posterior samples,
  # take the log, then sum across observations to get expected log predictive density
  stan.psm.cv.site.list[[i]]$elpd <- -2*sum(log(colMeans(exp(stan.psm.cv.site.list[[i]]$ll_psm))))
}

# Summarize model comparison results
stan.psm.cv.site.mods <- data.frame(model = names(stan.psm.cv.site.list),
                                    elpd = sapply(stan.psm.cv.site.list,
                                                  function(x) x$elpd),
                                    se_elpd = sapply(stan.psm.cv.site.list, function(x) 
                                      2*sqrt(ncol(x$ll_psm))*sd(log(colMeans(exp(x$ll_psm))))),
                                    d_elpd = NA,
                                    se_d_elpd = NA)

for(i in 1:length(stan.psm.cv.site.list))
{
  stan.psm.cv.site.mods$d_elpd[i] <- stan.psm.cv.site.mods$elpd[i] - min(stan.psm.cv.site.mods$elpd)
  m <- which.min(stan.psm.cv.site.mods$elpd)
  devs <- -2*(log(colMeans(exp(stan.psm.cv.site.list[[i]]$ll_psm))) - 
                log(colMeans(exp(stan.psm.cv.site.list[[m]]$ll_psm))))
  stan.psm.cv.site.mods$se_d_elpd[i] <- sqrt(length(devs))*sd(devs)
  rm(m);rm(devs)
}

stan.psm.cv.site.mods


#---------------------------------------------------------
# Fit "best" model to sites with PSM observations plus
# unsampled sites to generate predictions for the latter
#---------------------------------------------------------

# site-level covariates

# composition data (multiple categories):
# bound away from 0 and 1 and re-standardize to sum to 1,
# then transform to log ratios
X1.all <- as.matrix(lulc.roads.data.all[,c("ccap.hdev","ccap.ldev","ccap.mdev",
                                           "ccap.decid","ccap.evgrn","ccap.mxforest",
                                           "ccap.open","ccap.wetland","ccap.ag", "ccap.other")])
X1.all[X1.all==0] <- 1e-4
X1.all[X1.all==1] <- 1 - 1e-4
X1.all <- sweep(X1.all, 1, rowSums(X1.all), "/")
X1.all <- sweep(log(X1.all[,-ncol(X1.all)]), 1, log(X1.all[,ncol(X1.all)]), "-") 
X1.all <- sweep(sweep(X1.all, 2, colMeans(X1.all), "-"), 2, apply(X1.all, 2, sd), "/")

# proportion data:
# bound away from 0 and 1 and logit-transform
X2.all <- as.matrix(lulc.roads.data.all[,"nlcd.imperv",drop=F])
X2.all[X2.all==0] <- 1e-4
X2.all[X2.all==1] <- 1 - 1e-4
X2.all <- qlogis(X2.all)
X2.all <- sweep(sweep(X2.all, 2, colMeans(X2.all), "-"), 2, apply(X2.all, 2, sd), "/")

# nonnegative continuous data:
# scale to SD = 1, bound away from 0
X3.all <- as.matrix(lulc.roads.data.all[,c("roads.1","roads.2","roads.3","roads.4","roads.5",
                                           "traffic", "restoration","pop.census","pop.lscan")])
X3.all <- sweep(X3.all, 2, apply(X3.all, 2, sd), "/")
X3.all[X3.all==0] <- 1e-4
X3.all <- sweep(X3.all, 2, apply(X3.all, 2, sd), "/")

# reassemble all variables into matrix
X.all <- cbind(X2.all,X1.all,X3.all)
rm(list = c("X1.all","X2.all","X3.all"))

# indices of covariates assumed to be normally (gamma) distributed
normal.indx <- which(substring(dimnames(X.all)[[2]],1,4) %in% c("ccap","nlcd"))
gamma.indx <- which(!(substring(dimnames(X.all)[[2]],1,4) %in% c("ccap","nlcd")))

# Data for Stan
stan.dat.all <- list(S = nrow(X.all), 
                     D_normal = length(normal.indx), D_gamma = length(gamma.indx),
                     # X = t(X.all), 
                     X = X.all, 
                     L = 1,  # user-specified!
                     N = nrow(psm.all), 
                     site = as.numeric(psm.all$site),
                     ppt_su = as.vector(scale(psm.all$ppt.su/10, scale=F)),
                     ppt_fa = as.vector(scale(psm.all$ppt.fa/10, scale=F)),
                     I0_Z = 1,
                     I_su = 1,
                     I_su_Z = 1,
                     I_fa = 1,
                     I_fa_Z = 1,
                     n = psm.all$n,
                     n_psm = psm.all$n.psm,
                     I_fit = as.numeric(psm.all$data=="psm"),
                     I_lpd = rep(0, nrow(psm.all)))


# Function to generate initial values for chains
stan.init.all <- function() 
{
  if(exists("stan.dat.all"))
    for(i in names(stan.dat.all)) assign(i, stan.dat.all[[i]])
  # X <- t(X)
  D <- D_normal + D_gamma
  
  list(a0 = rnorm(D, c(colMeans(X[,1:D_normal]), colMeans(log(X[,-(1:D_normal)]))), 1),
       A_nid_vec = array(rnorm(D*L - L*(L-1)/2, 0, 1), dim = D*L - L*(L-1)/2),
       Z_nid = matrix(rnorm(S*L, 0, 1), nrow = S, ncol = L),
       phi = runif(D, 0.5, 1),
       mu_b0 = rnorm(1, -3, 1),
       b0_Z_nid = array(rnorm(L, 1, 0.5), dim = L),
       sigma_b0 = runif(1, 1, 2),
       b0_std = array(rnorm(S, 0, 0.1), dim = S),
       mu_b_su = rnorm(1, 0.3, 0.1),
       b_su_Z_nid = array(rnorm(L, -0.2, 0.1), dim = L),
       sigma_b_su = runif(1, 0.1, 1),
       b_su_std = array(rnorm(S, 0, 0.1), dim = S),
       mu_b_fa = rnorm(1, 0.03, 0.02),
       b_fa_Z_nid = array(rnorm(L, -0.03, 0.02), dim = L),
       sigma_b_fa = runif(1, 0.02, 0.1),
       b_fa_std = array(rnorm(S, 0, 0.1), dim = S),
       sigma_psm = runif(1, 0.5, 2),
       logit_p_psm_std = array(rnorm(N, 0, 0.1), dim = N))
}

# Fit it!
stan.psm.all <- stan(file = "cohoPSM_SEM_Stan.stan",
                     data = stan.dat.all, 
                     init = stan.init.all,
                     pars = c("a0","A","Z","phi",
                              "mu_b0","b0_Z","sigma_b0",
                              "mu_b_su","b_su_Z","sigma_b_su",
                              "mu_b_fa","b_fa_Z","sigma_b_fa",
                              "sigma_psm","p_psm"), 
                     chains = 3, iter = 50000, warmup = 10000, thin = 40, cores = 3,
                     control = list(adapt_delta = 0.8))

print(stan.psm.all, pars = c("p_psm","Z"), include = F)
traceplot(stan.psm.all, pars = c("A","lp__"), inc_warmup = F)
pairs(stan.psm.all, pars = c("a0[1]", "a0[2]", "A[1,1]", "A[2,1]", "Z[1,1]"), condition = "accept_stat__")
pairs(stan.psm.all, pars = paste("A[", normal.indx, ",1]", sep=""), labels = colnames(X.all)[normal.indx])
pairs(stan.psm.all, pars = paste("A[", gamma.indx, ",1]", sep=""), labels = colnames(X.all)[gamma.indx])

# Store Z and predicted P(PSM) in matrix
Z.all <- extract1(stan.psm.all,"Z")
Z.all <- Z.all[,stan.dat.all$site[psm.all$data == "pre"],1]
psm.pre <- extract1(stan.psm.all,"p_psm")
psm.pre <- psm.pre[, psm.all$data == "pre"]
site.names <- psm.all$site[psm.all$data=="pre"]
psm.pre <- data.frame(site = site.names,
                      Z.mean = colMeans(Z.all),
                      Z.se = apply(Z.all, 2, sd),
                      psm.observed = site.names %in% psm$site,
                      p.psm.mean = colMeans(psm.pre),
                      p.psm.025 = apply(psm.pre, 2, quantile, 0.025),
                      p.psm.975 = apply(psm.pre, 2, quantile, 0.975),
                      logit.p.psm.mean = colMeans(qlogis(psm.pre)),
                      logit.p.psm.se = apply(qlogis(psm.pre), 2, sd))
write.table(psm.pre, "PSM_predictions.txt", sep="\t", row.names=FALSE)
rm(site.names)


#==================================================================
# FIGURES
#==================================================================

# #-----------------------------------------
# # Posteriors of loadings
# #-----------------------------------------
# 
# dev.new(width=10, height=10)
# # png(filename="loadings.png", width=7, height=7, units = "in", res=300, type="cairo-png")
# par(mar=c(5.1,9.5,2.1,2.1))
# 
# A <- extract1(stan.psm, "A")[,,1]  
# bxp.dat <- boxplot(A, plot=F)
# bxp.dat$stats[3,] <- colMeans(A)
# bxp.dat$stats[c(2,4),] <- apply(A,2,quantile,c(0.05,0.95))
# bxp.dat$stats[c(1,5),] <- apply(A,2,quantile,c(0.025,0.975))
# 
# bxp(bxp.dat, xlab="Factor loading", ylab="", ylim=range(bxp.dat$stats), yaxt="n", 
#     horizontal=T, boxwex=0.4, outpch="", whisklty=1, staplewex = 0, 
#     las=1, cex.axis=1.2, cex.lab=1.5)
# axis(2, at=1:ncol(A), las=1, cex.axis=1.2,
#      labels=lulc.roads.labels$plot.label[match(colnames(X), lulc.roads.labels$data.label)])
# abline(v=0, lwd=1, lty=2)
# rm(A);rm(bxp.dat)
# # dev.off()

#-----------------------------------------
# Posteriors of loadings
#-----------------------------------------

dev.new(width=10, height=10)
# png(filename="loadings.png", width=10, height=10, units = "in", res=300, type="cairo-png")
par(mfcol=c(2,1), mar=c(1,14,2.1,2.1), oma=c(3.5,0,0,0))

gamma.labels <- lulc.roads.labels$plot.label[match(colnames(X)[gamma.indx], lulc.roads.labels$data.label)]
A <- extract1(stan.psm, "A")[,gamma.indx[order(gamma.labels)],1]  
bxp.dat <- boxplot(A, plot=F)
bxp.dat$stats[3,] <- colMeans(A)
bxp.dat$stats[c(2,4),] <- apply(A,2,quantile,c(0.05,0.95))
bxp.dat$stats[c(1,5),] <- apply(A,2,quantile,c(0.025,0.975))

bxp(bxp.dat, xlab="", ylab="", ylim=range(bxp.dat$stats), yaxt="n", 
    horizontal=T, boxwex=0.4, outpch="", whisklty=1, staplewex = 0, 
    las=1, cex.axis=1.5, cex.lab=2)
axis(2, at=1:ncol(A), las=1, cex.axis=1.8, labels=sort(gamma.labels))
abline(v=0, lwd=1, lty=2)
legend("topleft", legend="", title="A", bty="n", cex=2)

normal.labels <- lulc.roads.labels$plot.label[match(colnames(X)[normal.indx], lulc.roads.labels$data.label)]
A <- extract1(stan.psm, "A")[,normal.indx[order(normal.labels)],1]  
bxp.dat <- boxplot(A, plot=F)
bxp.dat$stats[3,] <- colMeans(A)
bxp.dat$stats[c(2,4),] <- apply(A,2,quantile,c(0.05,0.95))
bxp.dat$stats[c(1,5),] <- apply(A,2,quantile,c(0.025,0.975))

bxp(bxp.dat, xlab="", ylab="", ylim=range(bxp.dat$stats), yaxt="n", 
    horizontal=T, boxwex=0.4, outpch="", whisklty=1, staplewex = 0, 
    las=1, cex.axis=1.5, cex.lab=2)
axis(2, at=1:ncol(A), las=1, cex.axis=1.8, labels=sort(normal.labels))
abline(v=0, lwd=1, lty=2)
mtext("Factor loading", side=1, line=3, cex=2)
legend("topleft", legend="", title="B", bty="n", cex=2)

rm(A);rm(bxp.dat);rm(gamma.labels);rm(normal.labels)
# dev.off()

#---------------------------------------------------------
# Observed vs. posterior mean P(PSM) under "full" model
#---------------------------------------------------------

dev.new(width=10, height=10)
# png(filename="PSM_obs_vs_fit.png", width=10, height=10, units="in", res=300, type="cairo-png")
par(mar=c(5.1,4.5,4.1,2.1))

cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], maxColorValue = 255, alpha = 0.7*255)
plot(stan_mean(stan.psm,"p_psm"), psm$n.psm/psm$n, las=1, pch="", cex.lab=1.8, cex.axis=1.5,
     xlab="Fitted mortality", ylab="Observed mortality", xlim=c(0,1), ylim=c(0,1))
abline(0,1,lwd=2)
points(stan_mean(stan.psm,"p_psm"), psm$n.psm/psm$n, pch = 16, #pch=ifelse(psm$n==1, 1, 16),
       cex=0.4*sqrt(psm$n + 10), col=cc)
points(stan_mean(stan.psm,"p_psm"), psm$n.psm/psm$n, pch = 1, #pch=ifelse(psm$n==1, 1, 16),
       cex=0.4*sqrt(psm$n + 10), col="darkgray")
segments(apply(extract1(stan.psm,"p_psm"), 2, quantile, 0.025), psm$n.psm/psm$n,
         apply(extract1(stan.psm,"p_psm"), 2, quantile, 0.975), psm$n.psm/psm$n, col=cc)
segments(stan_mean(stan.psm,"p_psm"), binconf(psm$n.psm, psm$n)[,"Lower"],
         stan_mean(stan.psm,"p_psm"), binconf(psm$n.psm, psm$n)[,"Upper"], col=cc)
lgd <- legend("bottomright", inset = c(0.03,0), legend=c("1","50","100","150"), bty = "n", cex = 1.2,
              y.intersp=1.8, x.intersp = 1.3, pch=16, pt.cex=0.4*sqrt(c(1,50,100,150) + 10), col = cc)
text(lgd$rect$left, lgd$rect$top, pos = 4, "N females", cex = 1.2*par("cex"))
rm(cc);rm(lgd)
# dev.off()

#----------------------------------------------------------------------------------------------
# Plots of site-level intercept and slope estimates against site-level LU/LC factor scores
#----------------------------------------------------------------------------------------------

dev.new(width=7,height=15)
# png(filename="psm_site-level_regressions.png", width=7*0.75, height=15*0.75, units="in", res=300, type="cairo-png")
par(mfcol=c(3,1), mar=c(1.1,7,3,1.1), oma=c(4,0,0,0))

mod <- stan.psm.list[["ZZZ"]]$fit
site <- as.numeric(psm$site)

plot(stan_mean(mod,"Z"), stan_mean(mod,"b0"), pch="", las=1, cex.lab=2.2, cex.axis=1.8,
     xlim=range(summary(mod, "Z", probs = c(0.025,0.975))$summary[,c("2.5%","97.5%")]),
     ylim=range(summary(mod, "b0", probs = c(0.025,0.975))$summary[,c("2.5%","97.5%")]),
     xlab="", ylab="")
fits <- sapply(1:nrow(extract1(mod,"lp__")), 
               function(i) 
                 extract1(mod,"mu_b0")[i] + extract1(mod,"b0_Z")[i]*seq(par()$usr[1], par()$usr[2], length=50))
cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], alpha=0.5*255, maxColorValue=255)
polygon(c(seq(par()$usr[1], par()$usr[2], length=50), seq(par()$usr[2], par()$usr[1], length=50)),
        c(apply(fits,1,quantile,0.025), rev(apply(fits,1,quantile,0.975))),
        col=cc, border=NA)
abline(stan_mean(mod,"mu_b0"), stan_mean(mod,"b0_Z"), lwd=2)
points(stan_mean(mod,"Z"), stan_mean(mod,"b0"), pch=18, 
       cex=1.5*log(table(site)/min(table(site)) + 1))
segments(stan_mean(mod,"Z"), summary(mod, "b0", probs = 0.025)$summary[,"2.5%"], 
         stan_mean(mod,"Z"), summary(mod, "b0", probs = 0.975)$summary[,"97.5%"])
segments(summary(mod, "Z", probs = 0.025)$summary[,"2.5%"], stan_mean(mod,"b0"), 
         summary(mod, "Z", probs = 0.975)$summary[,"97.5%"], stan_mean(mod,"b0"))
mtext("Baseline", side=2, line=4.5, cex=2.2*par()$cex)
legend("topleft", legend = "", title = "A", cex=2.2, bty="n")
legend("bottomright", legend=c("1","2","5","10"), title="N years", cex = 1.5,
       y.intersp=1, bty = "n", pch=18, pt.cex=1.5*log(c(1,2,5,10) + 1))

plot(stan_mean(mod,"Z"), stan_mean(mod,"b_su"), pch="", las=1, cex.lab=2.2, cex.axis=1.8,
     xlim=range(summary(mod, "Z", probs = c(0.025,0.975))$summary[,c("2.5%","97.5%")]),
     ylim=range(summary(mod, "b_su", probs = c(0.025,0.975))$summary[,c("2.5%","97.5%")]),
     xlab="", ylab="")
fits <- sapply(1:nrow(extract1(mod,"lp__")), 
               function(i) 
                 extract1(mod,"mu_b_su")[i] + extract1(mod,"b_su_Z")[i]*seq(par()$usr[1], par()$usr[2], length=50))
cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], alpha=0.5*255, maxColorValue=255)
polygon(c(seq(par()$usr[1], par()$usr[2], length=50), seq(par()$usr[2], par()$usr[1], length=50)),
        c(apply(fits,1,quantile,0.025), rev(apply(fits,1,quantile,0.975))),
        col=cc, border=NA)
abline(stan_mean(mod,"mu_b_su"), stan_mean(mod,"b_su_Z"), lwd=2)
points(stan_mean(mod,"Z"), stan_mean(mod,"b_su"), pch=18,
       cex=1.5*log(table(site)/min(table(site)) + 1))
segments(stan_mean(mod,"Z"), summary(mod, "b_su", probs = 0.025)$summary[,"2.5%"], 
         stan_mean(mod,"Z"), summary(mod, "b_su", probs = 0.975)$summary[,"97.5%"])
segments(summary(mod, "Z", probs = 0.025)$summary[,"2.5%"], stan_mean(mod,"b_su"), 
         summary(mod, "Z", probs = 0.975)$summary[,"97.5%"], stan_mean(mod,"b_su"))
mtext(expression("Summer rain effect (cm" * {}^-1 * ")"), side=2, line=4.5, cex=2.2*par()$cex)
legend("topleft", legend = "", title = "B", cex=2.2, bty="n")

plot(stan_mean(mod,"Z"), stan_mean(mod,"b_fa"), pch="", las=1, cex.lab=2.2, cex.axis=1.8,
     xlim=range(summary(mod, "Z", probs = c(0.025,0.975))$summary[,c("2.5%","97.5%")]),
     ylim=range(summary(mod, "b_fa", probs = c(0.025,0.975))$summary[,c("2.5%","97.5%")]),
     xlab="", ylab="")
fits <- sapply(1:nrow(extract1(mod,"lp__")), 
               function(i) 
                 extract1(mod,"mu_b_fa")[i] + extract1(mod,"b_fa_Z")[i]*seq(par()$usr[1], par()$usr[2], length=50))
cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], alpha=0.5*255, maxColorValue=255)
polygon(c(seq(par()$usr[1], par()$usr[2], length=50), seq(par()$usr[2], par()$usr[1], length=50)),
        c(apply(fits,1,quantile,0.025), rev(apply(fits,1,quantile,0.975))),
        col=cc, border=NA)
abline(stan_mean(mod,"mu_b_fa"), stan_mean(mod,"b_fa_Z"), lwd=2)
points(stan_mean(mod,"Z"), stan_mean(mod,"b_fa"), pch=18,
       cex=1.5*log(table(site)/min(table(site)) + 1))
segments(stan_mean(mod,"Z"), summary(mod, "b_fa", probs = 0.025)$summary[,"2.5%"], 
         stan_mean(mod,"Z"), summary(mod, "b_fa", probs = 0.975)$summary[,"97.5%"])
segments(summary(mod, "Z", probs = 0.025)$summary[,"2.5%"], stan_mean(mod,"b_fa"), 
         summary(mod, "Z", probs = 0.975)$summary[,"97.5%"], stan_mean(mod,"b_fa"))
mtext(expression("Fall rain effect (cm" * {}^-1 * ")"), side=2, line=4.5, cex=2.2*par()$cex)
mtext("Urbanization", side=1, line=4, cex=2.2*par()$cex)
legend("topleft", legend = "", title =  "C", cex=2.2, bty="n")
rm(list=c("mod","site","fits","cc"))
# dev.off()

#---------------------------------------------------------
# Posterior predictive checking:
# Simulate data (n.psm and X) from the "full" model,
# then plot the actual data against the predictions
# (with posterior predictive credible intervals)
#---------------------------------------------------------

# PSM data
mod <- stan.psm
pp.p.psm <- extract1(mod, "p_psm")
pp.n.psm <- sapply(1:ncol(pp.p.psm), 
                   function(j) rbinom(nrow(pp.p.psm), size = psm$n[j], prob = pp.p.psm[,j]))

dev.new()
plot(psm$n.psm, colMeans(pp.n.psm), xlab = "Observed", ylab = "Posterior predictive", las = 1,
     main = expression(n[PSM]), cex.lab = 1.5, cex.axis = 1.2, cex.main=1.5, cex = 1.2,
     ylim = range(apply(pp.n.psm, 2, quantile, c(0.025,0.975))))
segments(psm$n.psm, apply(pp.n.psm, 2, quantile, 0.025),
         psm$n.psm, apply(pp.n.psm, 2, quantile, 0.975))
abline(0,1,col="gray")


# LU/LC data (transformed, i.e. the X matrix)
pp.g.mu.X <- extract1(mod, "g_mu_X")
pp.phi <- extract1(mod, "phi")
pp.X <- array(NA, dim(pp.g.mu.X))
for(i in 1:dim(pp.X)[2])
  for(j in 1:dim(pp.X)[3])
    pp.X[,i,j] <- ifelse(rep(j, dim(pp.X)[1]) %in% normal.indx,
                         rnorm(dim(pp.X)[1], pp.g.mu.X[,i,j], pp.phi[,j]),
                         rgamma(dim(pp.X)[1], 
                                shape = pp.phi[,j], 
                                rate = pp.phi[,j]/exp(pp.g.mu.X[,i,j])))

dev.new(height = 10, width = 12)
par(mfrow = c(4,5))
for(j in 1:ncol(X))
{
  plot(X[,j], colMeans(pp.X[,,j]), xlab = "Observed", ylab = "Posterior predictive",
       main = colnames(X)[j], cex.lab = 1.5, cex.axis = 1.2, cex.main=1.5, cex = 1.2,
       ylim = range(apply(pp.X[,,j], 2, quantile, c(0.025,0.975))), 
       log = ifelse(j %in% gamma.indx, "xy", ""))
  segments(X[,j], apply(pp.X[,,j], 2, quantile, 0.025),
           X[,j], apply(pp.X[,,j], 2, quantile, 0.975))
  abline(0,1,col="gray")
}

rm(list = c("mod","pp.p.psm","pp.n.psm","pp.g.mu.X","pp.phi","pp.X"))

#---------------------------------------------------------
# Correlation plot of landscape variables
#---------------------------------------------------------

c1 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))(200)

X.plot <- X
X.plot[,gamma.indx] <- log(X.plot[,gamma.indx])
dimnames(X.plot)[[2]] <- lulc.roads.labels$plot.label[match(dimnames(X)[[2]], lulc.roads.labels$data.label)]
R <- cor(X.plot)

dev.new(width = 10, height = 10, mar = c(1,1,1,1))
png(filename="landscape.corrplot.ugly.png", width=10*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
# corrplot(R, diag = F, method = "ellipse", order = "original", 
#          col = c1, tl.col = "black", tl.cex = 1.2) 
corrplot(R, diag = T, method = "color", order = "original",
         col = c1, tl.col = "black", tl.cex = 1.2)
dev.off()
rm(c1);rm(X.plot);rm(R)




