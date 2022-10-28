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
	real mat_s[N];
	int<lower = 1> n_years;
	int<lower = 1, upper = n_years> year[N];
}

parameters {
	real<lower = 0> sigma_year;
	real beta;
	real intercept;
	real alpha_raw[n_years];
}

transformed parameters{
  real b_exp[n_years];
  for (n in 1:n_years){
    b_exp[n] = 0 + sigma_year*alpha_raw[n] + beta*mat_s[n]; // non-centered parameterization
  }
}

model {
	// likelihood
	for (i in 1:N){
	  int ayear;
	  ayear = year[i];
	  x[i] ~ paretocounts(b_exp[ayear], xmin[i], xmax[i], counts[i]); 
	  }
	  
	//priors
  target += normal_lpdf(intercept|-1.5, 0.3);
  alpha_raw ~ std_normal();
  
  //hyperpriors
  target += normal_lpdf(sigma_year|0, 0.1); //half-normal since lower is zero. (via McElreath 2020 page 407-408)
	target += normal_lpdf(beta|0, 0.1);
	
}



