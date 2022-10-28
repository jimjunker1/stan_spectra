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
	int<lower = 1> n_groups;
	int<lower = 1, upper = n_groups> year[N];
}

parameters {
	// vector[n_groups] b_exp;
	real<lower = 0> sigma_year;
	// real mu_year;
	real beta;
	real intercept;
	real alpha_raw[n_groups];
	real sigma;
	real sigma_x;
}

transformed parameters{
  real b_exp[n_groups];
  for (n in 1:n_groups){
    b_exp[n] = 0 + sigma_year*alpha_raw[n]; // non-centered parameterization
  }
}

model {
	// likelihood
	for (i in 1:N){
	  int ayear;
	  ayear = year[i];
	  x[i] ~ paretocounts(intercept + b_exp[ayear] + beta*mat_s[i] + sigma*mat_s[i], xmin[i], xmax[i], counts[i]); 
	  }
	  
	//priors
  // target += normal_lpdf(b_exp|0, sigma_year);
  // target += normal_lpdf(b_exp|0, 0.5);
  target += normal_lpdf(intercept|-1.5, 0.5);
  alpha_raw ~ std_normal();
  sigma ~ normal(0, sigma_x);
  
  //hyperpriors
  // target += normal_lpdf(mu_year|1.5, 0.5);
  target += normal_lpdf(sigma_year|0, 0.1); //half-normal since lower is zero. (via McElreath 2020 page 407-408)
	target += normal_lpdf(beta|0, 0.2);
	target += exponential_lpdf(sigma_x|3);
	
}



