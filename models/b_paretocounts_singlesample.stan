functions{
  real paretocounts_lpdf(real x, real b_exp, real xmin, real xmax, real counts){
    return(counts*(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(x)));
  }
}
    
data {
	int<lower=0> N;
	vector <lower = 0>[N] x;
	vector <lower = 0>[N] xmin;
	vector <lower = 0>[N] xmax;
	vector <lower = 0>[N] counts;
}

parameters {
	real b_exp;
}


model {
	b_exp ~ normal(-2, 1);
	for (i in 1:N){
	  x[i] ~ paretocounts(b_exp, xmin[i], xmax[i], counts[i]);
	  }
}
