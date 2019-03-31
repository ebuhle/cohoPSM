pairs(cbind(fitted(psm.glmer.list[[8]]), fitted(psm.glmer.list[[97]]), fitted(psm.glmer.list[[84]]),
fitted(psm.glmer.list[[102]]),fitted(psm.glmer.list[[117]]),fitted(psm.glmer.list[[69]]),
fitted(psm.glmer.list[[86]]),fitted(psm.glmer.list[[25]])))


plot(fitted(psm.glmer.list[[8]]), psm.glmer.list[[8]]@y)
abline(0,1)

summary(psm.glmer.list[[8]])
summary(update(psm.glmer.list[[8]], subset = site!="Fortson"))



blah <- cor(psm[,row.names(psm.var.wts)], method="pearson")
write.table(blah, "blah.txt", sep="\t")
blah2 <- cor(psm.glmer.tab[,row.names(psm.var.wts)], use="pairwise.complete.obs")
write.table(blah2, "blah2.txt", sep="\t")




blah <- glm(psm ~ Impervious + dense.urban + light.med.urb + grass.shrub.crops + 
  conifer.deciduous + Apt.condo + Commercial + Industrial + Parks.open.space +
  Residential + FC00.Road + Nonlocal.Road, 
  data = psm, weights = n.total, family = binomial)

blah <- glmer(psm ~ dense.urban + Parks.open.space + (1|site), 
  data = psm, weights = n.total, family = binomial)



blah <- lapply(psm.glmer.list, function(x) 
  {m <- vcov(x); 
   dimnames(m) <- list(names(fixef(x)), names(fixef(x))); 
   m <- sweep(m, 1, sqrt(diag(m)), "/");
   m <- sweep(m, 2, diag(m), "/")
   m})

blah2 <- matrix(0,12,12)
dimnames(blah2) <- list(row.names(psm.var.wts), row.names(psm.var.wts))

for(i in dimnames(blah2)[[1]])
  for(j in dimnames(blah2)[[2]])
    for(k in as.numeric(row.names(psm.glmer.tab)[1:8]))
    {
      if(any(row.names(blah[[k]])==i) & any(row.names(blah[[k]])==j))
         if(abs(blah[[k]][i,j]) > abs(blah2[i,j]))
            blah2[i,j] <- blah[[k]][i,j]
    }

write.table(blah2, "max.param.corr.txt",sep="\t")