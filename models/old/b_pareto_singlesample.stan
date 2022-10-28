functions{
  real paretocustom_lpdf(real x, real b_exp, real xmin, real xmax){
    if(b_exp != -1)
    return(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(x));
    else
    return(log(log(xmin) - log(xmax)) + b_exp*log(x));
  }
}
    
data {
	int<lower=0> N;
	vector <lower = 0>[N] x;
	vector <lower = 0>[N] xmin;
	vector <lower = 0>[N] xmax;
	// vector <lower = 0>[N] counts;
}

parameters {
	real b_exp;
}


model {
	b_exp ~ normal(-2, 1);
	for (i in 1:N){
	  x[i] ~ paretocustom(b_exp, xmin[i], xmax[i]);
	  }
}
