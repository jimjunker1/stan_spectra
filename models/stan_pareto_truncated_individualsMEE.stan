functions{
  real paretocounts_lpdf(vector x, real b){
    vector[num_elements(x)] temp; //array that stores densities
        real sumx;
        for(i in 1:num_elements(x)){  //loop over data and store each density
            temp[i] = log(x[i]); 
        }
    sumx = sum(temp);
    real xmin;
    real xmax;
    xmin = min(x);
    xmax = max(x);
    return(num_elements(x)*log((b + 1)/(xmax^(b+1) - xmin^(b+1))) + b*sumx);
    }
}
    
data {
	int<lower=0> N;
	vector[N] x;
}

parameters {
	real b;
}

model {
	b ~ normal(-2, 1);
	target += paretocounts_lpdf(x|b);
}



