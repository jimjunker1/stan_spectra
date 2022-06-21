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
	int<lower = 1> n_sites;
	int<lower = 1, upper = n_sites> site_no[N];
}

parameters {
	vector[n_groups] b_exp_year;
	real<lower = 0> sigma_year;
	real mu_year;
	real beta;
	vector[n_sites] b_exp_site;
	real<lower = 0> sigma_site;
	real mu_site;
}


model {
  //likelihood
	for (i in 1:N){
	  x[i] ~ paretocounts(b_exp_year[year[i]] + b_exp_site[site_no[i]] + beta*mat_s[i], xmin[i], xmax[i], counts[i]); 
	  }
	  
	//priors  
	target += normal_lpdf(b_exp_year|mu_year, sigma_year);
	target += normal_lpdf(b_exp_site|mu_site, sigma_site);
	
	//hyperpriors
  target += normal_lpdf(mu_year|0, 1);
  target += normal_lpdf(sigma_year|0, 0.1); //half-normal since lower is zero. (via McElreath 2020 page 407-408)
	target += normal_lpdf(beta|0, 0.1); 
  target += normal_lpdf(mu_site|0, 1);
  target += normal_lpdf(sigma_site|0, 0.1); //half-normal since lower is zero. (via McElreath 2020 page 407-408)
	target += normal_lpdf(beta|0, 0.1);   
}

generated quantities{
  real bexp_year;
  real bexp_site;
  bexp_year = normal_rng(mu_year, sigma_year);
  bexp_site = normal_rng(mu_site, sigma_site);
}

