#==================================================================
# SETUP
#==================================================================

options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
library(rstan)
library(loo)
library(shinystan)
library(rstanarm)
library(Hmisc)
library(corrplot)
library(here)

# load functions
source(here("analysis","functions","stan_mean.R"))
source(here("analysis","functions","extract1.R"))
source(here("analysis","functions","stan_init.R"))
source(here("analysis","functions","stan_init_cv.R"))
source(here("analysis","functions","kfold_partition.R"))

# read and wrangle data
source(here("analysis","cohoPSM1_data.R"))  

# load previously saved stanfit objects
if(file.exists(here("analysis","results","stan_psm.RData")))
  load(here("analysis","results","stan_psm.RData"))
if(file.exists(here("analysis","results","stan_psm_WAIC_LOO.RData")))
  load(here("analysis","results","stan_psm_WAIC_LOO.RData"))
if(file.exists(here("analysis","results","stan_psm_cv_year.RData")))
  load(here("analysis","results","stan_psm_cv_year.RData"))
if(file.exists(here("analysis","results","stan_psm_cv_site.RData")))
  load(here("analysis","results","stan_psm_cv_site.RData"))
if(file.exists(here("analysis","results","stan_psm_all.RData")))
  load(here("analysis","results","stan_psm_all.RData"))

#==================================================================
# STRUCTURAL EQUATION MODELS FOR LANDSCAPE DATA AND PSM
#==================================================================

#------------------------------------------------------
# Fit full model to sites with PSM observations
#------------------------------------------------------

## @knitr stan_psm_r2r
# Fit it!
stan_psm <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
                 data = stan_dat, 
                 init = lapply(1:3, function(i) stan_init(stan_dat)),
                 pars = c("a0","A","Z","phi","g_mu_X",
                          "mu_b0","b0_Z","sigma_b0","b0",
                          "mu_b_su","b_su_Z","sigma_b_su","b_su",
                          "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                          "sigma_psm","p_psm","ll_psm"), 
                 chains = 3, iter = 12000, warmup = 2000, thin = 5)

# Inspect and use shinystan to explore samples
print(stan_psm, prob = c(0.025, 0.5, 0.975), 
      pars = c("b0","b_su","b_ppt_su","b_fa","b_ppt_fa","g_mu_X","p_psm","ll_psm","Z"), include = F)
# launch_shinystan(stan_psm)

# Save stanfit
save(stan_psm, file = here("analysis","results","stan_psm.RData"))
## @knitr ignore

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
  cat("|", rep("*",i), rep(" ", length(stan_psm_list) - i), "| Working on model ", i, "/", 
      length(stan_psm_list), " (see Viewer for progress) \n\n", sep = "")
  
  fit <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
              data = stan_dat,
              init = lapply(1:3, function(i) stan_init(stan_dat)),
              pars = with(stan_dat, {
                c("a0","A","Z","phi",
                  "mu_b0", switch(I0_Z,NULL,"b0_Z"), "sigma_b0","b0",
                  switch(I_su,NULL,"mu_b_su"), switch(I_su_Z,NULL,"b_su_Z"),
                  switch(I_su,NULL,"sigma_b_su"), switch(I_su,NULL,"b_su"),
                  switch(I_fa,NULL,"mu_b_fa"), switch(I_fa_Z,NULL,"b_fa_Z"),
                  switch(I_fa,NULL,"sigma_b_fa"), switch(I_fa,NULL,"b_fa"),
                  "sigma_psm","p_psm","ll_psm")
              }),
              chains = 3, iter = 12000, warmup = 2000, thin = 5,
              control = list(stepsize = 0.05))
  
  # Store fitted object
  stan_psm_list[[i]] <- list(fit = fit, 
                             WAIC = waic(extract1(fit,"ll_psm")), 
                             LOO = loo(fit, pars = "ll_psm"))
}

# Summarize results in data frame
stan_psm_mods <- data.frame(model = names(stan_psm_list),
                            Dbar = sapply(stan_psm_list, function(x) -2*sum(stan_mean(x$fit,"ll_psm"))),
                            p_WAIC = sapply(stan_psm_list, function(x) x$WAIC$estimates["p_waic",1]),
                            WAIC = sapply(stan_psm_list, function(x) x$WAIC$estimates["waic",1]),
                            dWAIC = NA,
                            se_dWAIC = NA,
                            p_LOO = sapply(stan_psm_list, function(x) x$LOO$estimates["p_loo",1]),
                            LOO = sapply(stan_psm_list, function(x) x$LOO$estimates["looic",1]),
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
format(stan_psm_mods, nsmall = 2, digits = 2)

# Save objects
save(stan_psm_list, stan_psm_mods, file = here("analysis","results","stan_psm_WAIC_LOO.RData"))


#---------------------------------------------------------------------
# K-fold cross-validation over YEARS:
# Leave out one year of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations
#---------------------------------------------------------------------

# Shortlist of candidate models chosen by consensus between in-sample WAIC and LOO
stan_psm_cv_year_list <- vector("list", 4)
names(stan_psm_cv_year_list) <- c("100","1Z0","ZZ1","ZZZ")

# Loop over candidate models and then over years, fit to data, predict hold-out data,
# and store stan objects
for(i in 1:length(stan_psm_cv_year_list))
{
  stan_psm_cv_year_list[[i]] <- list(fit = vector("list",length(unique(psm$year))), 
                                     ll_psm = matrix(NA, 6000, nrow(psm)),
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
    cat("Working on model ", i, "/", length(stan_psm_cv_year_list), " and hold-out year ", 
        grep(j, sort(unique(psm$year))), "/", length(unique(psm$year)), 
        " (see Viewer for progress) \n\n", sep = "")
    fit <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
                data = stan_dat_cv_year, 
                init = lapply(1:3,function(i) stan_init_cv(stan_psm)),
                pars = with(stan_dat_cv_year, {
                  c("a0","A","Z","phi",
                    "mu_b0", switch(I0_Z,NULL,"b0_Z"), "sigma_b0","b0",
                    switch(I_su,NULL,"mu_b_su"), switch(I_su_Z,NULL,"b_su_Z"),
                    switch(I_su,NULL,"sigma_b_su"), switch(I_su,NULL,"b_su"),
                    switch(I_fa,NULL,"mu_b_fa"), switch(I_fa_Z,NULL,"b_fa_Z"),
                    switch(I_fa,NULL,"sigma_b_fa"), switch(I_fa,NULL,"b_fa"),
                    "sigma_psm","p_psm","ll_psm")
                }), 
                chains = 3, iter = 12000, warmup = 2000, thin = 5,
                control = list(stepsize = 0.05))
    
    # Store fitted object and log-likelihood matrix
    stan_psm_cv_year_list[[i]]$fit[[as.character(j)]] <- fit
    stan_psm_cv_year_list[[i]]$ll_psm[,psm$year == j] <- as.matrix(fit,"ll_psm")
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
save(stan_psm_cv_year_list, stan_psm_cv_year_mods, file = here("analysis","results","stan_psm_cv_year.RData"))


#---------------------------------------------------------------------
# K-fold cross-validation over SITES:
# Leave out one or more sites of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations.
# Sites are randomly partitioned into K = 10 groups that are
# roughly similar in size (i.e., number of observations).
#---------------------------------------------------------------------

partitions <- kfold_partition(psm_dat = psm, K = 10)
partitions  # check that the procedure found a "good" partition (roughly equal group sizes)
site_group <- partitions$group

# Shortlist of candidate models chosen by consensus between in-sample WAIC and LOO
stan_psm_cv_site_list <- vector("list", 4)
names(stan_psm_cv_site_list) <- c("100","1Z0","ZZ1","ZZZ")

# Loop over candidate models and then over groups of sites, fit to data, 
# predict hold-out data, and store stan objects
# (Note that this assumes stan_dat has been assigned in a previous code block)
for(i in 1:length(stan_psm_cv_site_list))
{
  stan_psm_cv_site_list[[i]] <- list(fit = vector("list",length(unique(site_group))),
                                     site_group = site_group,
                                     ll_psm = matrix(NA, 6000, nrow(psm)),
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
    cat("Working on model ", i, "/", length(stan_psm_cv_site_list), " and hold-out group ", 
        grep(j, sort(unique(site_group))), "/", length(unique(site_group)), 
        " (see Viewer for progress) \n\n", sep = "")
    fit <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
                data = stan_dat_cv_site, 
                init = lapply(1:3,function(i) stan_init_cv(stan_psm)),
                pars = with(stan_dat_cv_site, {
                  c("a0","A","Z","phi",
                    "mu_b0", switch(I0_Z,NULL,"b0_Z"), "sigma_b0","b0",
                    switch(I_su,NULL,"mu_b_su"), switch(I_su_Z,NULL,"b_su_Z"),
                    switch(I_su,NULL,"sigma_b_su"), switch(I_su,NULL,"b_su"),
                    switch(I_fa,NULL,"mu_b_fa"), switch(I_fa_Z,NULL,"b_fa_Z"),
                    switch(I_fa,NULL,"sigma_b_fa"), switch(I_fa,NULL,"b_fa"),
                    "sigma_psm","p_psm","ll_psm")
                }),
                chains = 3, iter = 12000, warmup = 2000, thin = 5,
                control = list(stepsize = 0.05))
    
    # Store fitted object and log-likelihood matrix
    stan_psm_cv_site_list[[i]]$fit[[as.character(j)]] <- fit
    stan_psm_cv_site_list[[i]]$ll_psm[,site_group == j] <- as.matrix(fit,"ll_psm")
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
save(stan_psm_cv_site_list, stan_psm_cv_site_mods, file = here("analysis","results","stan_psm_cv_site.RData"))


#---------------------------------------------------------
# Fit "best" (in fact, the full) model to sites with 
# PSM observations plus unsampled sites to generate 
# predictions for the latter
#---------------------------------------------------------

## @knitr stan_psm_all_r2r
# Fit it!
stan_psm_all <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
                     data = stan_dat_all, 
                     init = lapply(1:3, function(i) stan_init(stan_dat_all)),
                     pars = c("a0","A","Z","phi",
                              "mu_b0","b0_Z","sigma_b0","b0_std",
                              "mu_b_su","b_su_Z","sigma_b_su",
                              "mu_b_fa","b_fa_Z","sigma_b_fa",
                              "sigma_psm","p_psm"), 
                     chains = 3, iter = 12000, warmup = 2000, thin = 5,
                     control = list(stepsize = 0.05))

# Print and explore fit in shinystan
print(stan_psm_all, prob = c(0.025, 0.5, 0.975), pars = c("p_psm","Z","b0_std"), include = F)
# launch_shinystan(stan_psm_all)

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
save(stan_psm_all, psm_pre, file = here("analysis","results","stan_psm_all.RData"))
write.table(psm_pre, here("analysis","results","PSM_predictions.txt"), sep="\t", row.names=FALSE)
## @knitr ignore


#==================================================================
# FIGURES
#==================================================================

#-----------------------------------------
# Posteriors of loadings
#-----------------------------------------

dev.new(width=10, height=10)
# png(filename=here("analysis","results","figures","loadings.png"), 
#     width=10, height=10, units = "in", res=300, type="cairo-png")
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
# png(filename=here("analysis","results","figures","PSM_obs_vs_fit.png"), 
#     width=10, height=10, units="in", res=300, type="cairo-png")
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
# png(filename=here("analysis","results","figures","psm_site_level_regressions.png"), 
#     width=7*0.75, height=15*0.75, units="in", res=300, type="cairo-png")
par(mfcol=c(3,1), mar=c(1.1,7,3,1.1), oma=c(4,0,0,0))

mod <- stan_psm_list[["ZZZ"]]$fit
site <- as.numeric(psm$site)

# intercept                     
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

# summer precip effect
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

# fall precip effect
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

dev.new()
# png(filename=here("analysis","results","figures","landscape_corrplot.png"), 
#     width=10*0.9, height=10*0.9, units="in", res=200, type="cairo-png")
par(mar = c(1,1,1,1))
corrplot(R, diag = F, method = "ellipse", order = "original", 
         col = c1, tl.col = "black", tl.cex = 1.2)
# dev.off()
rm(c1);rm(X_plot);rm(R)




