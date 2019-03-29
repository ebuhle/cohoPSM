data {
  int<lower=1> K;  //number of indicators             
  int<lower=1> N; //number of observations
  int<lower=1> D; //number of dimensions
  int<lower=1> Q;  // number of off-diagonal elements              
  vector[K]  y[N]; //data matrix
 
}
parameters {    
  vector[D] F[N]; //factors
  vector[Q] offdiag;
  positive_ordered[D] diag;
  real<lower=0> sigma_L;
}

transformed parameters{
matrix[K,D]  L;// the loading matrix
{
int index;

for (j in 1:D) {
    L[j,j] <- diag[j]; //constrains the diagonal to positive ordered
    for (i in (j+1):K) {
      index <- index + 1;
      L[i,j] <- offdiag[index];
      
    }
}

for(i in 1:K){
     for(j in (i+1):D){
     L[i,j] <- 0; //constrains the upper off diagonal elements to zero
   }
}

}
}

model {

  offdiag ~normal(0,sigma_L); // priors of the loadings
  diag ~normal(0,sigma_L);          
  sigma_L ~ cauchy(0,5); 
     
  for (n in 1:N){
   F[n] ~ normal(0,1); //factor constraints 
   y[n] ~ normal(L*F[n],1);//the likelihood where the variance of residuals are constrained to 1.
 }

}