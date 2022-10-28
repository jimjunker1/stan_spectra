//custom likelihood from Edwards 2020.
functions{
  real paretocounts_lpdf(real x, real b_exp, real xmin, real xmax, real counts){
    return(counts*(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(x)));
  }
}
    
data {
	int<lower=0> N;
	real x[N];
	real xmin[N];
	real xmax[N];
	real counts[N];
	// real mat_s[N];
	int<lower = 1> n_sites;
	int<lower = 1, upper = n_sites> site[N];
}

parameters {
	real<lower = 0> sigma_site;
	// real beta;
	real intercept;
	real alpha_raw[n_sites];
}

transformed parameters{
  real b_exp[n_sites];
  for (n in 1:n_sites){
    b_exp[n] = 0 + sigma_site*alpha_raw[n]; // non-centered parameterization
  }
}

model {
	// likelihood
	for (i in 1:N){
	  int asite;
	  asite = site[i];
	  x[i] ~ paretocounts(intercept + b_exp[site[i]], xmin[i], xmax[i], counts[i]); 
	  }
	  
	//priors
  // target += normal_lpdf(intercept|-1.5, 0.5);
  alpha_raw ~ std_normal();
  
  //hyperpriors
  target += normal_lpdf(sigma_site|0, 0.1); //half-normal since lower is zero. (via McElreath 2020 page 407-408)
	// target += normal_lpdf(beta|0, 0.1);
	
}



