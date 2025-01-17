functions{
  real paretocounts_lpdf(real x, real lambda, real xmin, real xmax, real counts){
    if(lambda != -1)
    return(counts*(log((lambda+1) / ( xmax^(lambda+1) - xmin^(lambda+1))) + lambda*log(x)));
    else
    return(counts*(log(log(xmin) - log(xmax)) + lambda*log(x)));
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
	real lambda;
}

model {
	lambda ~ normal(-2, 1);
	for (i in 1:N){
	  x[i] ~ paretocounts(lambda, xmin[i], xmax[i], counts[i]);
	  }
}

generated quantities{
  vector[N] log_lik;
  vector[N] temp;{
        for(i in 1:N) {
        temp[i] = paretocounts_lpdf(x[i] | lambda, xmin[i], xmax[i], counts[i]);
        }
      }
    log_lik = to_vector(temp);
}


