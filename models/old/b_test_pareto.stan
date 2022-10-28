// functions{
//   real paretocounts_lpdf(real x, real b_exp, real xmin, real xmax, real counts){
//     return(counts*(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(x)));
//   }
// }
//     
data {
	int<lower=0> N;
	real x[N];
	real xmin[N];
	// real xmax[N];
	// real counts[N];
	real mat_s[N];
	// int<lower = 1> n_sites;
	// int<lower = 1, upper = n_sites> site[N]; //index of sites
	int<lower = 1> n_years;
	int<lower = 1, upper = n_years> year[N]; //index of years
	// real mat[N];
}

parameters {
	vector<lower = 0>[n_years] b_exp;
	real mu;
	real<lower = 0> sigma;
	// real b_exp;
	// real beta;
}


model {
	for (i in 1:N){
	  x[i] ~ pareto(b_exp[year[i]], xmin[i]);
	  }
	  b_exp ~ lognormal(1,1);
	  // mu ~ lognormal(1, 1);
	  // sigma ~ exponential(1);
	  // beta ~ normal(0, 1);
}
