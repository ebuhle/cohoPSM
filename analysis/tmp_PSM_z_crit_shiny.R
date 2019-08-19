#----------------------------------------------------------------------
# PSM Threshold vs. Critical Urbanization Decision Analysis
#
# Plot site-specific posterior predictive distributions of PSM risk 
# as a function of Z, with posterior medians of current conditions.
# Show psm_crit and corresponding site-specific z_crit values for a
# specified confidence level alpha.
# Highlight a selected site and show detail
#----------------------------------------------------------------------

psm_crit <- 0.3   # PSM threshold
alpha <- 0.9  # credibility level
prediction_level <- "site"  # "site" or "year"-within-site
WADOE_sites <- sort(unique(setdiff(psm_all$site, psm_all$site[stan_dat_all$I_fit==1])))
WADOE_nums <- sort(unique(setdiff(stan_dat_all$site, stan_dat_all$site[stan_dat_all$I_fit==1])))
show_site <- "400"
show_num <- which(WADOE_sites == show_site)
show_crv <- stan_dat_all$S  # arbitrary non-sampled site

Z_draws <- extract1(stan_psm_all, "Z")[,WADOE_nums,1]
Z <- colMedians(Z_draws)
cur <- sem_psm_predict(stan_psm_all, data = stan_dat_all, newsites = WADOE_nums, 
                       level = prediction_level, transform = TRUE)  # use estimated Z
PSM <- colMedians(cur$est)
z_out <- sem_z_crit(stan_psm_all, data = stan_dat_all, psm_crit = psm_crit, 
                    level = prediction_level, alpha = alpha)
newsites <- rep(show_crv, each = 50)
newZ <- matrix(rep(seq(min(Z, z_out$z_crit) - 0.1, max(Z, z_out$z_crit) + 0.1, length = 50), 
                   length(unique(newsites))), ncol = 1)
psm_pred <- sem_psm_predict(stan_psm_all, data = stan_dat_all, newsites = newsites, newZ = newZ,
                            level = prediction_level, transform = TRUE)$est
cur_show_site <- sem_psm_predict(stan_psm_all, data = stan_dat_all, newsites = show_num,
                                 newZ = z_out$z_crit[show_num], transform = TRUE)
psm_pred_show_site <- cur_show_site$est

c1 <- transparent("darkgray", 0.3)  # posterior predictive PSM curves, all sites
c2 <- "gray"  # highlighted site
c2t <- transparent(c2, 0.6)
dzcols <- color_values(z_out$delta_z, palette = t(col2rgb(cividis(256, direction = -1))))
dzcolst <- transparent(dzcols, 0.1)

dev.new(width = 7.5, height = 7)

par(mar = c(5.1, 4.2, 4.1, 4))

plot(Z, PSM, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlim = range(newZ), ylim = c(0,1), xaxs = "i",
     xlab = bquote("Urbanization (" * italic(Z) * ")"), ylab = "Predicted PSM")
# PSM vs. Z curve and current conditions
for(j in unique(newsites))
  lines(newZ[newsites==j], colMedians(psm_pred[,newsites==j]), col = c1)
points(Z, PSM, pch = 1, cex = 1.5, col = dzcolst)
# selected site: PSM vs. Z curve and gradient 
polygon(c(newZ[newsites==show_crv], rev(newZ[newsites==show_crv])),
        c(colQuantiles(psm_pred[,newsites==show_crv], probs = 0.05),
          rev(colQuantiles(psm_pred[,newsites==show_crv], probs = 0.95))),
        col = c2t, border = NA)
lines(newZ[newsites==show_crv], colMedians(psm_pred[,newsites==show_crv]), lwd = 3)
# selected site: posterior density of PSM at z_crit
vioplot2(psm_pred_show_site, at = z_out$z_crit[show_num], add = TRUE, 
         col = NULL, border = "black", lwd = 2, wex = 0.15, drawRect = FALSE, pchMed = "")
abline(v = z_out$z_crit[show_num], col = "red")
text(z_out$z_crit[show_num], par("usr")[3] - 0.03, bquote(italic(z)[crit]), adj = c(0.5,0), 
     xpd = TRUE, col = "red")
# selected site: delta_PSM (change in median PSM risk at z_crit relative to current)
segments(x0 = z_out$z_crit[show_num], x1 = Z[show_num], y0 = median(psm_pred_show_site), lty = 2)
arrows(x0 = Z[show_num], y0 = PSM[show_num], y1 = median(psm_pred_show_site), length = 0.1)
text(Z[show_num] + 0.02, 0.75*PSM[show_num] + 0.25*median(psm_pred_show_site), expression(Delta * "PSM"), adj = 0)
# selected site: posterior density of z; delta_z
vioplot2(Z_draws[,show_num], quantiles = alpha, horizontal = TRUE, at = PSM[show_num],
         add = TRUE, col = NULL, border = "black", lwd = 2, lwd.quantile = 2, wex = 0.05, 
         drawRect = FALSE, pchMed = "")
qz <- quantile(Z_draws[,show_num], alpha)
arrows(x0 = qz, x1 = z_out$z_crit[show_num], y0 = PSM[show_num], 
       col = dzcols[show_num], length = 0.1, lwd = 2)
text(0.25*qz + 0.75*z_out$z_crit[show_num], PSM[show_num] + 0.01, expression(Delta * italic(z)), 
     adj = c(0.5,0), col = dzcols[show_num])
# selected site: current conditions
points(Z[show_num], PSM[show_num], pch = ifelse(show_site %in% psm$site, 16, 1), 
       col = dzcols[show_num], cex = 1.5)
# PSM threshold and all z_crit values
abline(h = psm_crit, col = "red", lwd = 2)
text(par("usr")[1] - 0.02, psm_crit, bquote(PSM[crit]), adj = c(1,0.5), col = "red", xpd = TRUE)
rug(z_out$z_crit[WADOE_sites], col = transparent("red", 0.3))
shape::colorlegend(cividis(100, direction = -1, alpha = 0.9), 
                   zlim = round(range(z_out$delta_z)), dz = 1,
                   digit = 0, main = expression(Delta * italic(z)), main.cex = 1.5)

