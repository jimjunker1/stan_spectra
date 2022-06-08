functions{
  real paretocounts_lpdf(vector x, real b, vector c){
    vector[num_elements(x)] temp; //array that stores densities
        real sumcx;
        for(i in 1:num_elements(x)){  //loop over data and store each density
            temp[i] = c[i]*log(x[i]); 
        }
    sumcx = sum(temp);
    real xmin;
    real xmax;
    xmin = min(x);
    xmax = max(x);
    return(sum(c)*log((b + 1)/(xmax^(b+1) - xmin^(b+1))) + b*sumcx);
    }
}
    
data {
	int<lower=0> N;
	vector[N] x;
	vector[N] c;
}

parameters {
	real b;
}

model {
	b ~ normal(-2, 1);
	target += paretocounts_lpdf(x|b, c);
}



