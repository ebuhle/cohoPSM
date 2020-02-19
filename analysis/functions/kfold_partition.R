# Randomly partition sites into groups for cross-validation
# Sites are randomly partitioned into K groups that are
# roughly similar in size (i.e., number of observations).

kfold_partition <- function(psm_dat, K, N_random = 1000)
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
