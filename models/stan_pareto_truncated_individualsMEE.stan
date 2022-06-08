functions{
  real paretocounts_lpdf(vector x, real b, real xmin, real xmax){
    vector[num_elements(x)] temp; //array that stores densities
        real sumx;
        for(i in 1:num_elements(x)){  //loop over data and store each density
            temp[i] = log(x[i]); 
        }
    sumx = sum(temp);
    return(num_elements(x)*log((b + 1)/(xmax^(b+1) - xmin^(b+1))) + b*sumx);
    }
}
    
data {
	int<lower=0> N;
	real<lower=0> xmin;
	real<lower=0> xmax;
	vector[N] x;
}

parameters {
	real b;
}

model {
	b ~ normal(-2, 1);
	target += paretocounts_lpdf(x|b, xmin, xmax);
}



