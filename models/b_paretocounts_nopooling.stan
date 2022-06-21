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
	//real mat[N];
	int<lower = 1> n_groups;
	int<lower = 1, upper = n_groups> site[N]; //index of sites
	// real mat[N];
}

parameters {
	vector[n_groups] b_exp;
}


model {
	for (i in 1:N){
	  x[i] ~ paretocounts(b_exp[site[i]], xmin[i], xmax[i], counts[i]); //non-centered parameterization of random effect
	  }
}

