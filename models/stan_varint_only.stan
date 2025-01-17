functions{
  real paretocounts_lpdf(real x, real b_exp, real xmin, real xmax, real counts){
    if(b_exp != -1)
    return(counts*(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(x)));
    else
    return(counts*(log(log(xmin) - log(xmax)) + b_exp*log(x)));
  }
}
    
data {
	int<lower=0> N;
	real x[N];
	real xmin[N];
	real xmax[N];
	real counts[N];
	// real mat_s[N];
	// real gpp_s[N];
	int<lower = 1> n_group;
	int<lower = 1, upper = n_group> group[N];
	// int<lower = 1> n_sites;
	// int<lower = 1, upper = n_sites> site[N];
}

parameters {
	// real beta_mat;
	// real beta_gpp;
	// real beta_gpp_mat;
	real a;
	real alpha_raw_group[n_group];
	// real alpha_raw_site[n_sites];
  real<lower=0> sigma_group;
  // real<lower=0> sigma_site;
}


model {
	// likelihood
	for (i in 1:N){
	  x[i] ~ paretocounts(a + sigma_group*alpha_raw_group[group[i]], xmin[i], xmax[i], counts[i]); 
	  }
	  
	//priors
  target += normal_lpdf(a|-1.5, 0.5);
  alpha_raw_group ~ std_normal();
  // alpha_raw_site ~ std_normal();
  
  //hyperpriors
  target += exponential_lpdf(sigma_group|8); 
	// target += normal_lpdf(beta_mat|0, 0.1);
	// target += normal_lpdf(beta_gpp|0, 0.1);
	// target += normal_lpdf(beta_gpp_mat|0, 0.1);
	// target += exponential_lpdf(sigma_site|9);
	
}

