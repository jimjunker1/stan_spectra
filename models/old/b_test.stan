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
	// int<lower = 1> n_sites;
	// int<lower = 1, upper = n_sites> site[N]; //index of sites
	int<lower = 1> n_years;
	int<lower = 1, upper = n_years> year[N]; //index of years
	// real mat[N];
}

parameters {
	vector[n_years] b_exp;
	real mu;
	real<lower = 0> sigma;
	// real b_exp;
	// real beta;
}


model {
	for (i in 1:N){
	  x[i] ~ paretocounts(b_exp[year[i]], xmin[i], xmax[i], counts[i]);
	  }
	  b_exp ~ normal(mu, sigma);
	  mu ~ normal(-1.5, 0.5);
	  sigma ~ exponential(5);
	  // beta ~ normal(0, 1);
}
