setwd(file.path("~","cohoPSM","analysis"))
data_path <- file.path("~","cohoPSM","data")
options(device=windows)
library(rstan)
library(loo)
library(shinystan)
library(rstanarm)
library(Hmisc)
library(corrplot)
source("stan_mean.R")
source("extract1.R")

#==================================================================
# DATA
#==================================================================

#----------------------------------------------
# PSM SAMPLING SITES FOR ANALYSIS
#----------------------------------------------

# read in spawner data
spawner_data <- read.csv(file.path(data_path,"spawner_data.csv"), header=T)

# read in landscape data
spatial_data <- read.csv(file.path(data_path,"spatial_data.csv"), header=T)

# LU/LC and roads data
lulc_roads_data <- spatial_data[,c("site","watershed","area","ccap_decid","ccap_open","ccap_evgrn",
                                   "ccap_hdev","ccap_ldev","ccap_mdev","ccap_mxforest","ccap_wetland","ccap_ag",
                                   "ccap_other","roads1","roads2","roads3","roads4","roads5","traffic",
                                   "nlcd_imperv","restoration","pop_census","pop_lscan")]

# convert basin area from m2 to ha
lulc_roads_data$area <- lulc_roads_data$area/1e4

# convert to km of road/km2
lulc_roads_data[,grep("roads", names(lulc_roads_data))] <- lulc_roads_data[,grep("roads", names(lulc_roads_data))]/1000

# Pretty variable names for LU/LC and roads
lulc_roads_labels <- read.csv(file.path(data_path,"lulc_roads_labels.csv"), header=T)
lulc_roads_labels$plot_label <- ordered(as.character(lulc_roads_labels$plot_label),
                                        levels = c("Imperviousness","High developed","Medium developed",
                                                   "Low developed","Evergreen","Deciduous","Mixed forest",
                                                   "Open space","Wetland","Agriculture","Interstate",
                                                   "Principal arterial","Minor arterial","Collector arterial",
                                                   "Local roads","Traffic intensity","Restoration",
                                                   "U.S. Census","LandScan"))

# matrix of annual summer and fall precip with sites and years as rows
ppt_data <- spatial_data[,grep("site|ppt_su_2|ppt_fa_2", names(spatial_data))]
ppt_su_cols <- grepl("ppt_su_2", names(ppt_data))
ppt_fa_cols <- grepl("ppt_fa_2", names(ppt_data))
n_yrs_ppt <- sum(ppt_su_cols)

ppt_data <- data.frame(site = rep(ppt_data$site, each = n_yrs_ppt), 
                       year = rep(substring(names(ppt_data)[ppt_su_cols], 8), nrow(ppt_data)),
                       ppt_su = as.vector(t(as.matrix(ppt_data[,ppt_su_cols]))),
                       ppt_fa = as.vector(t(as.matrix(ppt_data[,ppt_fa_cols]))))

rm(list=c("ppt_su_cols","ppt_fa_cols","n_yrs_ppt"))

# merge basin-specific variables and annual su and fa ppt into spawner data
psm <- merge(spawner_data[,-1], ppt_data, by = c("site","year"))
psm <- merge(psm, lulc_roads_data, by = "site")
psm <- data.frame(psm[,c("site","watershed","year","n","n_psm","ppt_su","ppt_fa","area")], psm[,-(1:9)])
psm <- psm[order(psm$site,psm$year),]

#----------------------------------------------
# ADDITIONAL SITES FOR PSM PREDICTIONS
#----------------------------------------------

# read in landscape data
spatial_data_pre <- read.csv(file.path(data_path,"spatial_data_predict.csv"), header=T)
names(spatial_data_pre) <- gsub("ID", "site", names(spatial_data_pre))
spatial_data_pre$watershed <- NA

# convert basin area from m2 to ha
spatial_data_pre$area <- spatial_data_pre$area/1e4

# convert to km of road/km2
spatial_data_pre[,grep("roads", names(spatial_data_pre))] <- spatial_data_pre[,grep("roads", names(spatial_data_pre))]/1000

# only include basins w/ coho
spatial_data_pre <- spatial_data_pre[spatial_data_pre$coho==1,]  
spatial_data_pre$site <- factor(spatial_data_pre$site)

lulc_roads_data_pre <- spatial_data_pre[,names(lulc_roads_data)]

# Combine PSM-sampling and unsampled-site data frames
lulc_roads_data_all <- rbind(data.frame(data = "psm", lulc_roads_data), 
                             data.frame(data = "pre", lulc_roads_data_pre))

# Add rows to psm for unobserved sites, plus one extra row for each observed site
# containing no PSM data and having ppt set to the mean. These extra rows are 
# used to generate predictions for monitored sites under baseline precip conditions.
psm_all <- data.frame(lulc_roads_data_all, year = NA, n = 0, n_psm = 0,
                      ppt_su = mean(psm$ppt_su), ppt_fa = mean(psm$ppt_fa))
psm_all$data <- "pre"
psm_all <- psm_all[,c("data", names(psm))]
psm_all <- rbind(data.frame(data = "psm", psm), psm_all)
psm_all <- data.frame(ID = 1:nrow(psm_all), psm_all)


#==================================================================
# HIERARCHICAL REGRESSION MODELS FOR LANDSCAPE DATA AND PSM
#==================================================================

# Modify dataset
# precip, traffic and log(traffic) centered and scaled to SD = 1
psm_all_reg <- psm_all
psm_all_reg <- transform(psm_all_reg, ppt_su = scale(ppt_su), ppt_fa = scale(ppt_fa),
                         traffic = scale(traffic), log_traffic = scale(log(pmax(traffic, 0.1))))
##
psm_all_reg$Z <- stan_mean(stan.psm.all,"Z")[as.numeric(psm_all_reg$site)]
##

# Fit models
###
glmm_psm_Z <- stan_glmer(cbind(n_psm, n - n_psm) ~ (ppt_su + ppt_fa) * Z + 
                           (ppt_su + ppt_fa || site) + (1 | ID),
                       data = psm_all_reg, subset = data == "psm",
                       family = binomial("logit"),
                       prior_intercept = normal(0,3),
                       prior = normal(0,3),
                       prior_covariance = decov(),
                       chains = 3, cores = 3, iter = 2000, warmup = 1000,
                       control = list(adapt_delta = 0.9))

glmm_psm_Z
summary(glmm_psm_Z, prob = c(0.025, 0.5, 0.975), pars = "beta", include = FALSE)
###

glmm_psm <- stan_glmer(cbind(n_psm, n - n_psm) ~ (ppt_su + ppt_fa) * log_traffic + 
                         (ppt_su + ppt_fa || site) + (1 | ID),
                       data = psm_all_reg, subset = data == "psm",
                       family = binomial("logit"),
                       prior_intercept = normal(0,3),
                       prior = normal(0,3),
                       prior_covariance = decov(),
                       chains = 3, cores = 3, iter = 2000, warmup = 1000,
                       control = list(adapt_delta = 0.9))

glmm_psm
summary(glmm_psm, prob = c(0.025, 0.5, 0.975), pars = "beta", include = FALSE)


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
X1 <- as.matrix(lulc_roads_data[,c("ccap_hdev","ccap_ldev","ccap_mdev",
                                   "ccap_decid","ccap_evgrn","ccap_mxforest",
                                   "ccap_open","ccap_wetland","ccap_ag", "ccap_other")])
X1[X1==0] <- 1e-4
X1[X1==1] <- 1 - 1e-4
X1 <- sweep(X1, 1, rowSums(X1), "/")
X1 <- sweep(log(X1[,-ncol(X1)]), 1, log(X1[,ncol(X1)]), "-") 
X1 <- sweep(sweep(X1, 2, colMeans(X1), "-"), 2, apply(X1, 2, sd), "/")

# proportion data:
# bound away from 0 and 1 and logit-transform
X2 <- as.matrix(lulc_roads_data[,"nlcd_imperv",drop=F])
X2[X2==0] <- 1e-4
X2[X2==1] <- 1 - 1e-4
X2 <- qlogis(X2)
X2 <- sweep(sweep(X2, 2, colMeans(X2), "-"), 2, apply(X2, 2, sd), "/")

# nonnegative continuous data:
# scale to SD = 1, bound away from 0
X3 <- as.matrix(lulc_roads_data[,c("roads1","roads2","roads3","roads4","roads5","traffic",
                                   "restoration","pop_census","pop_lscan")])
X3 <- sweep(X3, 2, apply(X3, 2, sd), "/")
X3[X3==0] <- 1e-4
X3 <- sweep(X3, 2, apply(X3, 2, sd), "/")

# reassemble all variables into matrix
X <- cbind(X2,X1,X3)
rm(list = c("X1","X2","X3"))

# indices of covariates assumed to be normally (gamma) distributed
normal_indx <- grep("ccap|nlcd", colnames(X))
gamma_indx <- which(!grepl("ccap|nlcd", colnames(X)))

# Data for Stan
stan_dat <- list(S = nrow(X), 
                 D_normal = length(normal_indx), D_gamma = length(gamma_indx),
                 X = X, 
                 L = 1,  # user-specified!
                 N = nrow(psm), 
                 site = as.numeric(psm$site),
                 ppt_su = array(as.vector(scale(psm$ppt_su/10, scale = FALSE)), dim = nrow(psm)),
                 ppt_fa = array(as.vector(scale(psm$ppt_fa/10, scale = FALSE)), dim = nrow(psm)),
                 I0_Z = 1,
                 I_su = 1,
                 I_su_Z = 1,
                 I_fa = 1,
                 I_fa_Z = 1,
                 n = psm$n,
                 n_psm = psm$n_psm,
                 I_fit = rep(1, nrow(psm)),
                 I_lpd = rep(1, nrow(psm)))


# Function to generate initial values for chains
stan_init <- function(stan_dat) 
{
  with(stan_dat, {
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
  })
}

# Fit it!
stan_psm <- stan(file = "cohoPSM_SEM_Stan.stan",
                 data = stan_dat, 
                 init = lapply(1:3, function(i) stan_init(stan_dat)),
                 pars = c("a0","A","Z","phi","g_mu_X",
                          "mu_b0","b0_Z","sigma_b0","b0",
                          "mu_b_su","b_su_Z","sigma_b_su","b_su",
                          "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                          "sigma_psm","p_psm","ll_psm"), 
                 chains = 3, iter = 12000, warmup = 2000, thin = 5, cores = 3)

# Inspect summary use shinystan to explore samples
print(stan_psm, pars = c("g_mu_X","ll_psm","Z"), include = F)
launch_shinystan(stan_psm)

# Save stanfit
save(stan_psm, file = "stan_psm.RData")

#------------------------------------------------------------------
# In-sample model selection:
# Fit all candidate models to sites with PSM observations,
# compare expected out-of-sample performance via WAIC and PSIS-LOO
#------------------------------------------------------------------

# Create list to store results
stan_psm_list <- vector("list", 18)
k <- 1
for(b0_mod in c(1,"Z"))
  for(b_su_mod in c(0,1,"Z"))
    for(b_fa_mod in c(0,1,"Z"))
    {
      names(stan_psm_list)[k] <- paste(b0_mod, b_su_mod, b_fa_mod, sep="")
      k <- k + 1
    }

# Loop over candidate models, fit to data, and store stan objects
# (Note that this assumes stan_dat and stan_init have been assigned in
# the previous code block)
for(i in 1:length(stan_psm_list))
{
  # Assign binary in/out indicators
  stan_dat$I0_Z <- as.numeric(substring(names(stan_psm_list)[i],1,1) == "Z")
  stan_dat$I_su <- as.numeric(substring(names(stan_psm_list)[i],2,2) != "0")
  stan_dat$I_su_Z <- as.numeric(substring(names(stan_psm_list)[i],2,2) == "Z")
  stan_dat$I_fa <- as.numeric(substring(names(stan_psm_list)[i],3,3) != "0")
  stan_dat$I_fa_Z <- as.numeric(substring(names(stan_psm_list)[i],3,3) == "Z")
  
  # Fit model
  cat("Working on model", i, "(see Viewer for progress) \n")
  fit <- stan(file = "cohoPSM_SEM_Stan.stan",
              data = stan_dat, 
              init = stan_init,
              pars = c("a0","A","Z","phi",
                       "mu_b0","b0_Z","sigma_b0","b0",
                       "mu_b_su","b_su_Z","sigma_b_su","b_su",
                       "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                       "sigma_psm","p_psm","ll_psm"), 
              chains = 3, iter = 12000, warmup = 2000, thin = 5, cores = 3)
  
  # Store fitted object
  stan_psm_list[[i]] <- list(fit = fit, 
                             WAIC = waic(extract1(fit,"ll_psm")), 
                             LOO = loo(extract1(fit,"ll_psm")))
}

# Summarize results in data frame
stan_psm_mods <- data.frame(model = names(stan_psm_list),
                            Dbar = sapply(stan_psm_list, function(x) -2*sum(stan_mean(x$fit,"ll_psm"))),
                            p_WAIC = sapply(stan_psm_list, function(x) x$WAIC$p_waic),
                            WAIC = sapply(stan_psm_list, function(x) x$WAIC$waic),
                            dWAIC = NA,
                            se_dWAIC = NA,
                            p_LOO = sapply(stan_psm_list, function(x) x$LOO$p_loo),
                            LOO = sapply(stan_psm_list, function(x) x$LOO$looic),
                            dLOO = NA,
                            se_dLOO = NA)

for(i in 1:length(stan_psm_list))
{
  compare_waic <- compare(stan_psm_list[[i]]$WAIC, 
                          stan_psm_list[[which.min(stan_psm_mods$WAIC)]]$WAIC)
  stan_psm_mods$dWAIC[i] <- 2*compare_waic["elpd_diff"]
  stan_psm_mods$se_dWAIC[i] <- 2*compare_waic["se"]
  compare_loo <- compare(stan_psm_list[[i]]$LOO, 
                         stan_psm_list[[which.min(stan_psm_mods$LOO)]]$LOO)
  stan_psm_mods$dLOO[i] <- 2*compare_loo["elpd_diff"]
  stan_psm_mods$se_dLOO[i] <- 2*compare_loo["se"]
  rm(compare_waic); rm(compare_loo)
}

stan_psm_mods <- stan_psm_mods[order(stan_psm_mods$dLOO),]                            
stan_psm_mods[order(stan_psm_mods$dWAIC),]

# Save objects
save(stan_psm_list, stan_psm_mods, file = "stan_psm_WAIC_LOO.RData")


#---------------------------------------------------------------------
# K-fold cross-validation over YEARS:
# Leave out one year of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations
#---------------------------------------------------------------------

# Shortlist of candidate models chosen by consensus between in-sample WAIC and LOO
stan_psm_cv_year_list <- vector("list", 4)
names(stan_psm_cv_year_list) <- c("100","1Z0","ZZ1","ZZZ")

# Fit "full" model to complete data once to generate good starting values
# for leave-one-out runs
fit1 <- stan(file = "cohoPSM_SEM_Stan.stan",
             data = stan_dat, 
             init = stan_init,
             pars = names(stan_init()), 
             chains = 3, iter = 12000, warmup = 2000, thin = 5, cores = 3)

stan_init_cv <- function(fit)
{
  samples <- extract(fit)
  M <- nrow(samples$a0)
  i <- sample(M,1)
  
  with(samples, 
       list(a0 = a0[i,],
            A_nid_vec = array(A_nid_vec[i,], dim = ncol(A_nid_vec)),
            Z_nid = matrix(Z[i,,], nrow = dim(Z)[2], ncol = dim(Z)[3]),
            phi = phi[i,],
            mu_b0 = mu_b0[i],
            b0_Z_nid = array(b0_Z_nid[i,], dim = ncol(b0_Z_nid)),
            sigma_b0 = sigma_b0[i],
            b0_std = array(b0_std[i,], dim = ncol(b0_std)),
            mu_b_su = mu_b_su[i],
            b_su_Z_nid = array(b_su_Z_nid[i,], dim = ncol(b_su_Z_nid)),
            sigma_b_su = sigma_b_su[i],
            b_su_std = array(b_su_std[i,], dim = ncol(b_su_std)),
            mu_b_fa = mu_b_fa[i],
            b_fa_Z_nid = array(b_fa_Z_nid[i,], dim = ncol(b_fa_Z_nid)),
            sigma_b_fa = sigma_b_fa[i],
            b_fa_std = array(b_fa_std[i,], dim = ncol(b_fa_std)),
            sigma_psm = sigma_psm[i],
            logit_p_psm_std = array(logit_p_psm_std[i,], dim = ncol(logit_p_psm_std))))
}

# Loop over candidate models and then over years, fit to data, predict hold-out data,
# and store stan objects
# (Note that this assumes stan_dat was assigned in the previous code block)
for(i in 1:length(stan_psm_cv_year_list))
{
  stan_psm_cv_year_list[[i]] <- list(fit = vector("list",length(unique(psm$year))), 
                                     ll_psm = matrix(NA, 3000, nrow(psm)),
                                     elpd = NULL)
  names(stan_psm_cv_year_list[[i]]$fit) <- sort(unique(psm$year))
  stan_dat_cv_year <- stan_dat
  
  # Assign binary in/out indicators
  stan_dat_cv_year$I0_Z <- as.numeric(substring(names(stan_psm_cv_year_list)[i],1,1) == "Z")
  stan_dat_cv_year$I_su <- as.numeric(substring(names(stan_psm_cv_year_list)[i],2,2) != "0")
  stan_dat_cv_year$I_su_Z <- as.numeric(substring(names(stan_psm_cv_year_list)[i],2,2) == "Z")
  stan_dat_cv_year$I_fa <- as.numeric(substring(names(stan_psm_cv_year_list)[i],3,3) != "0")
  stan_dat_cv_year$I_fa_Z <- as.numeric(substring(names(stan_psm_cv_year_list)[i],3,3) == "Z")
  
  for(j in sort(unique(psm$year)))
  {
    stan_dat_cv_year$I_fit <- as.numeric(psm$year != j)  # training data
    stan_dat_cv_year$I_lpd <- as.numeric(psm$year == j)  # hold-out data
    
    # Fit model
    cat("Working on model", i, "and hold-out year", j, "(see Viewer for progress) \n")
    fit <- stan(file = "cohoPSM_SEM_Stan.stan",
                data = stan_dat_cv_year, 
                init = lapply(1:3,function(i) stan_init_cv(fit1)),
                pars = c("a0","A","Z","phi",
                         "mu_b0","b0_Z","sigma_b0","b0",
                         "mu_b_su","b_su_Z","sigma_b_su","b_su",
                         "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                         "sigma_psm","p_psm","ll_psm"), 
                chains = 3, iter = 12000, warmup = 2000, thin = 5, cores = 3)
    
    # Store fitted object and log-likelihood matrix
    stan_psm_cv_year_list[[i]]$fit[[as.character(j)]] <- fit
    stan_psm_cv_year_list[[i]]$ll_psm[,psm$year == j] <- matrix(extract1(fit,"ll_psm"), nrow = 3000)
  }
  
  # Average the posterior predictive density (not log density!) across posterior samples,
  # take the log, then sum across observations to get expected log predictive density
  stan_psm_cv_year_list[[i]]$elpd <- -2*sum(log(colMeans(exp(stan_psm_cv_year_list[[i]]$ll_psm))))
}

# Summarize model comparison results
stan_psm_cv_year_mods <- data.frame(model = names(stan_psm_cv_year_list),
                                    elpd = sapply(stan_psm_cv_year_list,
                                                  function(x) x$elpd),
                                    se_elpd = sapply(stan_psm_cv_year_list, function(x) 
                                      2*sqrt(ncol(x$ll_psm))*sd(log(colMeans(exp(x$ll_psm))))),
                                    d_elpd = NA,
                                    se_d_elpd = NA)

for(i in 1:length(stan_psm_cv_year_list))
{
  stan_psm_cv_year_mods$d_elpd[i] <- stan_psm_cv_year_mods$elpd[i] - min(stan_psm_cv_year_mods$elpd)
  m <- which.min(stan_psm_cv_year_mods$elpd)
  devs <- -2*(log(colMeans(exp(stan_psm_cv_year_list[[i]]$ll_psm))) - 
                log(colMeans(exp(stan_psm_cv_year_list[[m]]$ll_psm))))
  stan_psm_cv_year_mods$se_d_elpd[i] <- sqrt(length(devs))*sd(devs)
  rm(m);rm(devs)
}

stan_psm_cv_year_mods

# Save objects
save(stan_psm_cv_year_list, stan_psm_cv_year_mods, file = "stan_psm_cv_year.RData")


#---------------------------------------------------------------------
# K-fold cross-validation over SITES:
# Leave out one or more sites of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations.
# Sites are randomly partitioned into K = 10 groups that are
# roughly similar in size (i.e., number of observations).
#---------------------------------------------------------------------

# Randomly partition sites into groups
KfoldCV_partition <- function(psm_dat, K, N_random = 1000)
{
  grps <- matrix(NA, length(levels(psm$site)), N_random)
  pos <- 1:nrow(grps)
  target <- nrow(psm)/K
  for(j in 1:ncol(grps))
  {
    pos_permuted <- sample(nrow(grps), nrow(grps), replace = FALSE)
    N_site <- table(psm$site)[pos_permuted]
    grp_permuted <- c(1, rep(0, nrow(grps) - 1))
    for(i in 2:nrow(grps))
    {
      if(abs(sum(N_site[grp_permuted==grp_permuted[i-1]]) + N_site[i] - target) <
         abs(sum(N_site[grp_permuted==grp_permuted[i-1]]) - target))
      {
        grp_permuted[i] <- grp_permuted[i-1]
      } else {
        if(grp_permuted[i] < 10) grp_permuted[i] <- grp_permuted[i-1] + 1
      }
    }
    grps[,j] <- grp_permuted[order(pos_permuted)]
  }
  N_site <- table(psm$site)
  range_N_grp <- apply(grps, 2, function(grp) diff(range(tapply(N_site, grp, sum))))
  grp <- grps[,which.min(range_N_grp)]
  return(list(N_group = tapply(N_site, grp, sum), group = grp[match(psm$site, levels(psm$site))]))
}

partitions <- KfoldCV_partition(psm_dat = psm, K = 10)
partitions  # check that the procedure found a "good" partition (roughly equal group sizes)
site_group <- partitions$group

# Shortlist of candidate models chosen by consensus between in-sample WAIC and LOO
stan_psm_cv_site_list <- vector("list", 4)
names(stan_psm_cv_site_list) <- c("100","1Z0","ZZ1","ZZZ")

# Fit "full" model to complete data once to generate good starting values
# for leave-one-out runs
fit1 <- stan(file = "cohoPSM_SEM_Stan.stan",
             data = stan_dat, 
             init = stan_init,
             pars = names(stan_init()), 
             chains = 3, iter = 12000, warmup = 2000, thin = 5, cores = 3)

stan_init_cv <- function(fit)
{
  samples <- extract(fit)
  M <- nrow(samples$a0)
  i <- sample(M,1)
  
  with(samples,
       list(a0 = a0[i,],
            A_nid_vec = array(A_nid_vec[i,], dim = ncol(A_nid_vec)),
            Z_nid = matrix(Z[i,,], nrow = dim(Z)[2], ncol = dim(Z)[3]),
            phi = phi[i,],
            mu_b0 = mu_b0[i],
            b0_Z_nid = array(b0_Z_nid[i,], dim = ncol(b0_Z_nid)),
            sigma_b0 = sigma_b0[i],
            b0_std = array(b0_std[i,], dim = ncol(b0_std)),
            mu_b_su = mu_b_su[i],
            b_su_Z_nid = array(b_su_Z_nid[i,], dim = ncol(b_su_Z_nid)),
            sigma_b_su = sigma_b_su[i],
            b_su_std = array(b_su_std[i,], dim = ncol(b_su_std)),
            mu_b_fa = mu_b_fa[i],
            b_fa_Z_nid = array(b_fa_Z_nid[i,], dim = ncol(b_fa_Z_nid)),
            sigma_b_fa = sigma_b_fa[i],
            b_fa_std = array(b_fa_std[i,], dim = ncol(b_fa_std)),
            sigma_psm = sigma_psm[i],
            logit_p_psm_std = array(logit_p_psm_std[i,], dim = ncol(logit_p_psm_std))))
}

# Loop over candidate models and then over groups of sites, fit to data, 
# predict hold-out data, and store stan objects
# (Note that this assumes stan_dat has been assigned in a previous code block)
for(i in 1:length(stan_psm_cv_site_list))
{
  stan_psm_cv_site_list[[i]] <- list(fit = vector("list",length(unique(site_group))),
                                     site_group = site_group,
                                     ll_psm = matrix(NA, 3000, nrow(psm)),
                                     elpd = NULL)
  names(stan_psm_cv_site_list[[i]]$fit) <- sort(unique(site_group))
  stan_dat_cv_site <- stan_dat
  
  # Assign binary in/out indicators
  stan_dat_cv_site$I0_Z <- as.numeric(substring(names(stan_psm_cv_site_list)[i],1,1) == "Z")
  stan_dat_cv_site$I_su <- as.numeric(substring(names(stan_psm_cv_site_list)[i],2,2) != "0")
  stan_dat_cv_site$I_su_Z <- as.numeric(substring(names(stan_psm_cv_site_list)[i],2,2) == "Z")
  stan_dat_cv_site$I_fa <- as.numeric(substring(names(stan_psm_cv_site_list)[i],3,3) != "0")
  stan_dat_cv_site$I_fa_Z <- as.numeric(substring(names(stan_psm_cv_site_list)[i],3,3) == "Z")
  
  for(j in sort(unique(site_group)))
  {
    stan_dat_cv_site$I_fit <- as.numeric(site_group != j)  # training data
    stan_dat_cv_site$I_lpd <- as.numeric(site_group == j)  # hold-out data
    
    # Fit model
    cat("Working on model", i, "and hold-out group", j, "(see Viewer for progress) \n")
    fit <- stan(file = "cohoPSM_SEM_Stan.stan",
                data = stan_dat_cv_site, 
                init = lapply(1:3,function(i) stan_init_cv(fit1)),
                pars = c("a0","A","Z","phi",
                         "mu_b0","b0_Z","sigma_b0","b0",
                         "mu_b_su","b_su_Z","sigma_b_su","b_su",
                         "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                         "sigma_psm","p_psm","ll_psm"), 
                chains = 3, iter = 12000, warmup = 2000, thin = 5, cores = 3)
    
    # Store fitted object and log-likelihood matrix
    stan_psm_cv_site_list[[i]]$fit[[as.character(j)]] <- fit
    stan_psm_cv_site_list[[i]]$ll_psm[,site_group == j] <- matrix(extract1(fit,"ll_psm"), nrow = 3000)
  }
  
  # Average the posterior predictive density (not log density!) across posterior samples,
  # take the log, then sum across observations to get expected log predictive density
  stan_psm_cv_site_list[[i]]$elpd <- -2*sum(log(colMeans(exp(stan_psm_cv_site_list[[i]]$ll_psm))))
}

# Summarize model comparison results
stan_psm_cv_site_mods <- data.frame(model = names(stan_psm_cv_site_list),
                                    elpd = sapply(stan_psm_cv_site_list,
                                                  function(x) x$elpd),
                                    se_elpd = sapply(stan_psm_cv_site_list, function(x) 
                                      2*sqrt(ncol(x$ll_psm))*sd(log(colMeans(exp(x$ll_psm))))),
                                    d_elpd = NA,
                                    se_d_elpd = NA)

for(i in 1:length(stan_psm_cv_site_list))
{
  stan_psm_cv_site_mods$d_elpd[i] <- stan_psm_cv_site_mods$elpd[i] - min(stan_psm_cv_site_mods$elpd)
  m <- which.min(stan_psm_cv_site_mods$elpd)
  devs <- -2*(log(colMeans(exp(stan_psm_cv_site_list[[i]]$ll_psm))) - 
                log(colMeans(exp(stan_psm_cv_site_list[[m]]$ll_psm))))
  stan_psm_cv_site_mods$se_d_elpd[i] <- sqrt(length(devs))*sd(devs)
  rm(m);rm(devs)
}

stan_psm_cv_site_mods

# Save objects
save(stan_psm_cv_site_list, stan_psm_cv_site_mods, file = "stan_psm_cv_site.RData")


#---------------------------------------------------------
# Fit "best" model to sites with PSM observations plus
# unsampled sites to generate predictions for the latter
#---------------------------------------------------------

# site-level covariates

# composition data (multiple categories):
# bound away from 0 and 1 and re-standardize to sum to 1,
# then transform to log ratios
X1_all <- as.matrix(lulc_roads_data_all[,c("ccap_hdev","ccap_ldev","ccap_mdev",
                                           "ccap_decid","ccap_evgrn","ccap_mxforest",
                                           "ccap_open","ccap_wetland","ccap_ag", "ccap_other")])
X1_all[X1_all==0] <- 1e-4
X1_all[X1_all==1] <- 1 - 1e-4
X1_all <- sweep(X1_all, 1, rowSums(X1_all), "/")
X1_all <- sweep(log(X1_all[,-ncol(X1_all)]), 1, log(X1_all[,ncol(X1_all)]), "-") 
X1_all <- sweep(sweep(X1_all, 2, colMeans(X1_all), "-"), 2, apply(X1_all, 2, sd), "/")

# proportion data:
# bound away from 0 and 1 and logit-transform
X2_all <- as.matrix(lulc_roads_data_all[,"nlcd_imperv",drop=F])
X2_all[X2_all==0] <- 1e-4
X2_all[X2_all==1] <- 1 - 1e-4
X2_all <- qlogis(X2_all)
X2_all <- sweep(sweep(X2_all, 2, colMeans(X2_all), "-"), 2, apply(X2_all, 2, sd), "/")

# nonnegative continuous data:
# scale to SD = 1, bound away from 0
X3_all <- as.matrix(lulc_roads_data_all[,c("roads1","roads2","roads3","roads4","roads5",
                                           "traffic", "restoration","pop_census","pop_lscan")])
X3_all <- sweep(X3_all, 2, apply(X3_all, 2, sd), "/")
X3_all[X3_all==0] <- 1e-4
X3_all <- sweep(X3_all, 2, apply(X3_all, 2, sd), "/")

# reassemble all variables into matrix
X_all <- cbind(X2_all, X1_all, X3_all)
rm(list = c("X1_all","X2_all","X3_all"))

# indices of covariates assumed to be normally (gamma) distributed
normal_indx <- grep("ccap|nlcd", colnames(X_all))
gamma_indx <- which(!grepl("ccap|nlcd", colnames(X_all)))

# Data for Stan
stan_dat_all <- list(S = nrow(X_all), 
                     D_normal = length(normal_indx), D_gamma = length(gamma_indx),
                     # X = t(X_all), 
                     X = X_all, 
                     L = 1,  # user-specified!
                     N = nrow(psm_all), 
                     site = as.numeric(psm_all$site),
                     ppt_su = array(as.vector(scale(psm_all$ppt_su/10, scale = FALSE)), dim = nrow(psm_all)),
                     ppt_fa = array(as.vector(scale(psm_all$ppt_fa/10, scale = FALSE)), dim = nrow(psm_all)),
                     I0_Z = 1,
                     I_su = 1,
                     I_su_Z = 1,
                     I_fa = 1,
                     I_fa_Z = 1,
                     n = psm_all$n,
                     n_psm = psm_all$n_psm,
                     I_fit = as.numeric(psm_all$data=="psm"),
                     I_lpd = rep(0, nrow(psm_all)))


# Function to generate initial values for chains
stan_init_all <- function(stan_dat_all) 
{
  with(stan_dat_all, {
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
  })
}

# Fit it!
stan_psm_all <- stan(file = "cohoPSM_SEM_Stan.stan",
                     data = stan_dat_all, 
                     init = lapply(1:3, function(i) stan_init_all(stan_dat_all)),
                     pars = c("a0","A","Z","phi",
                              "mu_b0","b0_Z","sigma_b0",
                              "mu_b_su","b_su_Z","sigma_b_su",
                              "mu_b_fa","b_fa_Z","sigma_b_fa",
                              "sigma_psm","p_psm"), 
                     chains = 3, iter = 12000, warmup = 2000, thin = 5, cores = 3)


# Print and explore fit in shinystan
print(stan_psm_all, pars = c("p_psm","Z"), include = F)
launch_shinystan(stan_psm_all)

# Store Z and predicted P(PSM) in matrix
Z_all <- extract1(stan_psm_all,"Z")
Z_all <- Z_all[,stan_dat_all$site[psm_all$data == "pre"],1]
psm_pre <- extract1(stan_psm_all,"p_psm")
psm_pre <- psm_pre[, psm_all$data == "pre"]
site_names <- psm_all$site[psm_all$data=="pre"]
psm_pre <- data.frame(site = site_names,
                      Z_mean = colMeans(Z_all),
                      Z_se = apply(Z_all, 2, sd),
                      psm_observed = site_names %in% psm$site,
                      p_psm_mean = colMeans(psm_pre),
                      p_psm_lo = apply(psm_pre, 2, quantile, 0.025),
                      p_psm_up = apply(psm_pre, 2, quantile, 0.975),
                      logit_p_psm_mean = colMeans(qlogis(psm_pre)),
                      logit_p_psm_se = apply(qlogis(psm_pre), 2, sd))

# Save objects
save(stan_psm_all, psm_pre, file = "stan_psm_allbasins.RData")
write.table(psm_pre, "PSM_predictions.txt", sep="\t", row.names=FALSE)


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
# A <- extract1(stan_psm, "A")[,,1]  
# bxp_dat <- boxplot(A, plot=F)
# bxp_dat$stats[3,] <- colMeans(A)
# bxp_dat$stats[c(2,4),] <- apply(A,2,quantile,c(0.05,0.95))
# bxp_dat$stats[c(1,5),] <- apply(A,2,quantile,c(0.025,0.975))
# 
# bxp(bxp_dat, xlab="Factor loading", ylab="", ylim=range(bxp_dat$stats), yaxt="n", 
#     horizontal=T, boxwex=0.4, outpch="", whisklty=1, staplewex = 0, 
#     las=1, cex.axis=1.2, cex.lab=1.5)
# axis(2, at=1:ncol(A), las=1, cex.axis=1.2,
#      labels=lulc_roads_labels$plot_label[match(colnames(X), lulc_roads_labels$data.label)])
# abline(v=0, lwd=1, lty=2)
# rm(A);rm(bxp_dat)
# # dev.off()

#-----------------------------------------
# Posteriors of loadings
#-----------------------------------------

dev.new(width=10, height=10)
# png(filename="loadings.png", width=10, height=10, units = "in", res=300, type="cairo-png")
par(mfcol=c(2,1), mar=c(1,14,2.1,2.1), oma=c(3.5,0,0,0))

gamma_labels <- lulc_roads_labels$plot_label[match(colnames(X)[gamma_indx], lulc_roads_labels$data.label)]
A <- extract1(stan_psm, "A")[,gamma_indx[order(gamma_labels)],1]  
bxp_dat <- boxplot(A, plot=F)
bxp_dat$stats[3,] <- colMeans(A)
bxp_dat$stats[c(2,4),] <- apply(A,2,quantile,c(0.05,0.95))
bxp_dat$stats[c(1,5),] <- apply(A,2,quantile,c(0.025,0.975))

bxp(bxp_dat, xlab="", ylab="", ylim=range(bxp_dat$stats), yaxt="n", 
    horizontal=T, boxwex=0.4, outpch="", whisklty=1, staplewex = 0, 
    las=1, cex.axis=1.5, cex.lab=2)
axis(2, at=1:ncol(A), las=1, cex.axis=1.8, labels=sort(gamma_labels))
abline(v=0, lwd=1, lty=2)
legend("topleft", legend="", title="A", bty="n", cex=2)

normal_labels <- lulc_roads_labels$plot_label[match(colnames(X)[normal_indx], lulc_roads_labels$data.label)]
A <- extract1(stan_psm, "A")[,normal_indx[order(normal_labels)],1]  
bxp_dat <- boxplot(A, plot=F)
bxp_dat$stats[3,] <- colMeans(A)
bxp_dat$stats[c(2,4),] <- apply(A,2,quantile,c(0.05,0.95))
bxp_dat$stats[c(1,5),] <- apply(A,2,quantile,c(0.025,0.975))

bxp(bxp_dat, xlab="", ylab="", ylim=range(bxp_dat$stats), yaxt="n", 
    horizontal=T, boxwex=0.4, outpch="", whisklty=1, staplewex = 0, 
    las=1, cex.axis=1.5, cex.lab=2)
axis(2, at=1:ncol(A), las=1, cex.axis=1.8, labels=sort(normal_labels))
abline(v=0, lwd=1, lty=2)
mtext("Factor loading", side=1, line=3, cex=2)
legend("topleft", legend="", title="B", bty="n", cex=2)

rm(A);rm(bxp_dat);rm(gamma_labels);rm(normal_labels)
# dev.off()

#---------------------------------------------------------
# Observed vs. posterior mean P(PSM) under "full" model
#---------------------------------------------------------

dev.new(width=10, height=10)
# png(filename="PSM_obs_vs_fit.png", width=10, height=10, units="in", res=300, type="cairo-png")
par(mar=c(5.1,4.5,4.1,2.1))

cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], maxColorValue = 255, alpha = 0.7*255)
plot(stan_mean(stan_psm,"p_psm"), psm$n_psm/psm$n, las=1, pch="", cex.lab=1.8, cex.axis=1.5,
     xlab="Fitted mortality", ylab="Observed mortality", xlim=c(0,1), ylim=c(0,1))
abline(0,1,lwd=2)
points(stan_mean(stan_psm,"p_psm"), psm$n_psm/psm$n, pch = 16, #pch=ifelse(psm$n==1, 1, 16),
       cex=0.4*sqrt(psm$n + 10), col=cc)
points(stan_mean(stan_psm,"p_psm"), psm$n_psm/psm$n, pch = 1, #pch=ifelse(psm$n==1, 1, 16),
       cex=0.4*sqrt(psm$n + 10), col="darkgray")
segments(apply(extract1(stan_psm,"p_psm"), 2, quantile, 0.025), psm$n_psm/psm$n,
         apply(extract1(stan_psm,"p_psm"), 2, quantile, 0.975), psm$n_psm/psm$n, col=cc)
segments(stan_mean(stan_psm,"p_psm"), binconf(psm$n_psm, psm$n)[,"Lower"],
         stan_mean(stan_psm,"p_psm"), binconf(psm$n_psm, psm$n)[,"Upper"], col=cc)
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

mod <- stan_psm_list[["ZZZ"]]$fit
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
# Simulate data (n_psm and X) from the "full" model,
# then plot the actual data against the predictions
# (with posterior predictive credible intervals)
#---------------------------------------------------------

# PSM data
mod <- stan_psm
pp_p_psm <- extract1(mod, "p_psm")
pp_n_psm <- sapply(1:ncol(pp_p_psm), 
                   function(j) rbinom(nrow(pp_p_psm), size = psm$n[j], prob = pp_p_psm[,j]))

dev.new()
plot(psm$n_psm, colMeans(pp_n_psm), xlab = "Observed", ylab = "Posterior predictive", las = 1,
     main = expression(n[PSM]), cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5, cex = 1.2,
     ylim = range(apply(pp_n_psm, 2, quantile, c(0.025,0.975))))
segments(psm$n_psm, apply(pp_n_psm, 2, quantile, 0.025),
         psm$n_psm, apply(pp_n_psm, 2, quantile, 0.975))
abline(0,1,col="gray")


# LU/LC data (transformed, i.e. the X matrix)
pp_g_mu_X <- extract1(mod, "g_mu_X")
pp_phi <- extract1(mod, "phi")
pp_X <- array(NA, dim(pp_g_mu_X))
for(i in 1:dim(pp_X)[2])
  for(j in 1:dim(pp_X)[3])
    pp_X[,i,j] <- ifelse(rep(j, dim(pp_X)[1]) %in% normal_indx,
                         rnorm(dim(pp_X)[1], pp_g_mu_X[,i,j], pp_phi[,j]),
                         rgamma(dim(pp_X)[1], 
                                shape = pp_phi[,j], 
                                rate = pp_phi[,j]/exp(pp_g_mu_X[,i,j])))

dev.new(height = 10, width = 12)
par(mfrow = c(4,5))
for(j in 1:ncol(X))
{
  plot(X[,j], colMeans(pp_X[,,j]), xlab = "Observed", ylab = "Posterior predictive",
       main = colnames(X)[j], cex.lab = 1.5, cex.axis = 1.2, cex.main=1.5, cex = 1.2,
       ylim = range(apply(pp_X[,,j], 2, quantile, c(0.025,0.975))), 
       log = ifelse(j %in% gamma_indx, "xy", ""))
  segments(X[,j], apply(pp_X[,,j], 2, quantile, 0.025),
           X[,j], apply(pp_X[,,j], 2, quantile, 0.975))
  abline(0,1,col="gray")
}

rm(list = c("mod","pp_p_psm","pp_n_psm","pp_g_mu_X","pp_phi","pp_X"))

#---------------------------------------------------------
# Correlation plot of landscape variables
#---------------------------------------------------------

c1 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))(200)

X_plot <- X
X_plot[,gamma_indx] <- log(X_plot[,gamma_indx])
dimnames(X_plot)[[2]] <- lulc_roads_labels$plot_label[match(dimnames(X)[[2]], lulc_roads_labels$data_label)]
R <- cor(X_plot)

dev.new(width = 10, height = 10, mar = c(1,1,1,1))
png(filename="landscape_corrplot_ugly.png", width=10*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
# corrplot(R, diag = F, method = "ellipse", order = "original", 
#          col = c1, tl.col = "black", tl.cex = 1.2) 
corrplot(R, diag = T, method = "color", order = "original",
         col = c1, tl.col = "black", tl.cex = 1.2)
dev.off()
rm(c1);rm(X_plot);rm(R)




