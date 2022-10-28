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
	int<lower = 1> n_sites;
	int<lower = 1, upper = n_sites> site[N];
}

parameters {
	real beta;
	real a;
	vector[n_years] alpha_year;
	vector[n_sites] alpha_site;
  real<lower=0> sigma_year;
  real<lower=0> sigma_site;
}



model {
	// likelihood
	for (i in 1:N){
	  x[i] ~ paretocounts(a + alpha_year[year[i]] + alpha_site[site[i]] + beta*mat_s[i], xmin[i], xmax[i], counts[i]); 
	  }
	  
	//priors
  target += normal_lpdf(a|-1.5, 0.5);
  target += normal_lpdf(alpha_year|0, sigma_year);
  target += normal_lpdf(alpha_site|0, sigma_site);
  
  //hyperpriors
  target += normal_lpdf(sigma_year|0, 1); //half-normal since lower is zero. (via McElreath 2020 page 407-408)
	target += normal_lpdf(beta|0, 1);
	target += normal_lpdf(sigma_site|0, 1);
	
}
