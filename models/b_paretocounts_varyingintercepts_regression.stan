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
	//real mat[N];
	int<lower = 1> n_groups;
	int<lower = 1, upper = n_groups> year[N]; 
	// real mat[N];
}

parameters {
	vector[n_groups] b_exp;
	real<lower = 0> sigma_year;
	real mu_year;
	real beta;
}


model {
  target += normal_lpdf(b_exp|mu_year, sigma_year);
  target += normal_lpdf(mu_year|0, 1);
  target += normal_lpdf(sigma_year|0, 0.1); //half-normal since lower is zero. (via McElreath 2020 page 407-408)
	target += normal_lpdf(beta|0, 0.1);
	for (i in 1:N){
	  x[i] ~ paretocounts(b_exp[year[i]] + beta*mat_s[i], xmin[i], xmax[i], counts[i]); 
	  }
}

generated quantities{
  real bexp_overall;
  bexp_overall = normal_rng(mu_year, sigma_year);
}

