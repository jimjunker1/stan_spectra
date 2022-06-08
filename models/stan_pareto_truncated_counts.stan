functions{
  real paretocounts_lpdf(vector x, real b, real xmin, real xmax, vector c){
    vector[num_elements(x)] temp; //array that stores densities
        real sumcx;
        for(i in 1:num_elements(x)){  //loop over data and store each density
            temp[i] = c[i]*log(x[i]); 
        }
    sumcx = sum(temp);
    return(sum(c)*log((b + 1)/(xmax^(b+1) - xmin^(b+1))) + b*sumcx);
    }
}
    
data {
	int<lower=0> N;
	real<lower=0> xmin;
	real<lower=0> xmax;
	vector[N] x;
	vector[N] c;
}

parameters {
	real b;
}

model {
	b ~ normal(-2, 1);
	target += paretocounts_lpdf(x|b, xmin, xmax, c);
}



