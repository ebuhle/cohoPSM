if(!require(SEMPSM)) devtools::install_github("ebuhle/SEMPSM")
library(SEMPSM)
library(matrixStats)
library(dplyr)
library(here)

#==================================================================
# DATA
#==================================================================

#----------------------------------------------
# PSM SAMPLING SITES FOR ANALYSIS
#----------------------------------------------
## @knitr data_psm_sites

# read in spawner data
spawner_data <- read.csv(here("data","spawner_data.csv"), header=T)

# read in landscape data
spatial_data <- read.csv(here("data","spatial_data.csv"), header=T)

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
lulc_roads_labels <- read.csv(here("data","lulc_roads_labels.csv"), header=T)
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

## @ knitr ignore
#----------------------------------------------
# ADDITIONAL SITES FOR PSM PREDICTIONS
#----------------------------------------------
## @ knitr data_predict_sites

# read in landscape data
spatial_data_pre <- read.csv(here("data","spatial_data_predict.csv"), header=T)
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

## @ knitr ignore
#------------------------------------------------------
# Assemble data for sites with PSM observations
# in Stan-friendly format
#------------------------------------------------------
## @knitr stan_data_psm_sites

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
stan_dat <- stan_data(psm = psm, X = X, normal_indx = normal_indx, gamma_indx = gamma_indx,
                      L = 1, I0_Z = 1, I_su = 1, I_su_Z = 1, I_fa = 1, I_fa_Z = 1,
                      I_fit = rep(1, nrow(psm)), I_lpd = rep(1, nrow(psm)))

## @knitr ignore
#------------------------------------------------------
# Assemble data for all sites, including out-of sample
# subbasins (no PSM observations) to be predicted,
# in Stan-friendly format
#------------------------------------------------------
## @knitr stan_data_all_sites

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
stan_dat_all <- stan_data(psm = psm_all, X = X_all, 
                          normal_indx = normal_indx, gamma_indx = gamma_indx,
                          L = 1, I0_Z = 1, I_su = 1, I_su_Z = 1, I_fa = 1, I_fa_Z = 1,
                          I_fit = as.numeric(psm_all$data=="psm"), 
                          I_lpd = rep(0, nrow(psm_all)))

## @knitr ignore
#-----------------------------------------------------------------
# CONSERVATION ATTRIBUTE DATA FOR TNC PRIORITIZATION ANALYSIS
#-----------------------------------------------------------------
## @knitr tnc_prioritization_data

# Species occurrence in WADOE basins from SalmonScape
salmonscape <- read.csv(here("data","salmonscape","WA_integrated_Fish_coho_chinook_chum_with_WADOE.csv"), 
                        header = TRUE, skip = 1)
names(salmonscape)[1] <- "site"
salmonscape$site <- as.character(salmonscape$site)
names(salmonscape)[-1] <- tolower(names(salmonscape)[-1])
names(salmonscape)[-1] <- sapply(strsplit(names(salmonscape)[-1], ".", fixed = TRUE), 
                                 function(x) paste(c(rev(x), "km"), collapse = "_"))
salmonscape[is.na(salmonscape)] <- 0  # NA means not present
salmonscape[,-1] <- salmonscape[,-1]/1000 # convert m to km

salmonscape$N_spp <- as.numeric(salmonscape$coho_total_km > 0) +
  as.numeric(rowAnys(select_at(salmonscape, vars(contains("chin_total"))) > 0)) + 
  as.numeric(rowAnys(select_at(salmonscape, vars(contains("chum_total"))) > 0))

# 66 WADOE basins in the original analysis are missing from the SalmonScape data.
# Pad SalmonScape data with NAs for those basins 
salmonscape <- psm_all %>% select(c(data, site)) %>% mutate(site = as.character(site)) %>% 
  left_join(salmonscape)

# Pretty names for selected variables
salmonscape_labels <- read.csv(here("data","salmonscape_labels.csv"), header=T, stringsAsFactors = FALSE)








