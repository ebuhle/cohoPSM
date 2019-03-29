model { 
	 C <- 10000 
	 ##Likelihood 
	 for(i in 1:n) {
		 eta[i,1] <- inprod(all.params[1,2:(num.lv+1)],lvs[i,])
		 y[i,1] ~ dbeta(ilogit(all.params[1,1] + eta[i,1])*all.params[1,num.lv+2],(1-ilogit(all.params[1,1] + eta[i,1]))*all.params[1,num.lv+2])

		 eta[i,2] <- inprod(all.params[2,2:(num.lv+1)],lvs[i,])
		 y[i,2] ~ dbeta(ilogit(all.params[2,1] + eta[i,2])*all.params[2,num.lv+2],(1-ilogit(all.params[2,1] + eta[i,2]))*all.params[2,num.lv+2])

		 eta[i,3] <- inprod(all.params[3,2:(num.lv+1)],lvs[i,])
		 y[i,3] ~ dbeta(ilogit(all.params[3,1] + eta[i,3])*all.params[3,num.lv+2],(1-ilogit(all.params[3,1] + eta[i,3]))*all.params[3,num.lv+2])

		 eta[i,4] <- inprod(all.params[4,2:(num.lv+1)],lvs[i,])
		 y[i,4] ~ dbeta(ilogit(all.params[4,1] + eta[i,4])*all.params[4,num.lv+2],(1-ilogit(all.params[4,1] + eta[i,4]))*all.params[4,num.lv+2])

		 eta[i,5] <- inprod(all.params[5,2:(num.lv+1)],lvs[i,])
		 y[i,5] ~ dbeta(ilogit(all.params[5,1] + eta[i,5])*all.params[5,num.lv+2],(1-ilogit(all.params[5,1] + eta[i,5]))*all.params[5,num.lv+2])

		 eta[i,6] <- inprod(all.params[6,2:(num.lv+1)],lvs[i,])
		 y[i,6] ~ dbeta(ilogit(all.params[6,1] + eta[i,6])*all.params[6,num.lv+2],(1-ilogit(all.params[6,1] + eta[i,6]))*all.params[6,num.lv+2])

		 eta[i,7] <- inprod(all.params[7,2:(num.lv+1)],lvs[i,])
		 y[i,7] ~ dgamma(exp(all.params[7,1] + eta[i,7])*all.params[7,num.lv+2], all.params[7,num.lv+2])

		 eta[i,8] <- inprod(all.params[8,2:(num.lv+1)],lvs[i,])
		 y[i,8] ~ dgamma(exp(all.params[8,1] + eta[i,8])*all.params[8,num.lv+2], all.params[8,num.lv+2])

		 eta[i,9] <- inprod(all.params[9,2:(num.lv+1)],lvs[i,])
		 y[i,9] ~ dgamma(exp(all.params[9,1] + eta[i,9])*all.params[9,num.lv+2], all.params[9,num.lv+2])

		 eta[i,10] <- inprod(all.params[10,2:(num.lv+1)],lvs[i,])
		 y[i,10] ~ dgamma(exp(all.params[10,1] + eta[i,10])*all.params[10,num.lv+2], all.params[10,num.lv+2])

		 eta[i,11] <- inprod(all.params[11,2:(num.lv+1)],lvs[i,])
		 y[i,11] ~ dgamma(exp(all.params[11,1] + eta[i,11])*all.params[11,num.lv+2], all.params[11,num.lv+2])

		 eta[i,12] <- inprod(all.params[12,2:(num.lv+1)],lvs[i,])
		 y[i,12] ~ dgamma(exp(all.params[12,1] + eta[i,12])*all.params[12,num.lv+2], all.params[12,num.lv+2])

		 eta[i,13] <- inprod(all.params[13,2:(num.lv+1)],lvs[i,])
		 y[i,13] ~ dgamma(exp(all.params[13,1] + eta[i,13])*all.params[13,num.lv+2], all.params[13,num.lv+2])

		 eta[i,14] <- inprod(all.params[14,2:(num.lv+1)],lvs[i,])
		 y[i,14] ~ dgamma(exp(all.params[14,1] + eta[i,14])*all.params[14,num.lv+2], all.params[14,num.lv+2])



	 } for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } } 

	 ## Latent variables
	 for(i in 1:p) { 
		 all.params[i,1] ~ dnorm(0,0.01) } ## Species intercept 

	 for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { 
		 all.params[i,j] <- 0 } } ## Constraints to 0 on upper diagonal
	 for(i in 1:num.lv) { 
		 all.params[i,i+1] ~ dunif(0,100) } ## Sign constraints on diagonal elements
	 for(i in 2:num.lv) { for(j in 2:i) { 
		 all.params[i,j] ~ dnorm(0,0.01) } } ## Free lower diagonals
	 for(i in (num.lv+1):p) { for(j in 2:(num.lv+1)) { 
		 all.params[i,j] ~ dnorm(0,0.01) } } ## All other elements
	 for(j in 1:p) { all.params[j,num.lv+2] ~ dunif(0,1e+05) }

	 }
