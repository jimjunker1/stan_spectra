library(tidyr)
library(dplyr)


calc_exponent_maxlik = function(data, bmin = -2.5, bmax = -0.5, length.out = 100){
  negLL = data %>% 
    tidyr::expand_grid(b_exp = seq(-2.5, -0.5, length.out = 100)) %>% 
    dplyr::mutate(negLL = -no_m2*(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(dw))) %>% 
    dplyr::select(ID, negLL, b_exp) %>% 
    dplyr::group_by(ID, b_exp) %>% 
    dplyr::summarize(sumnegLL = sum(negLL)) %>% 
    dplyr::group_by(ID) %>% 
    dplyr::filter(sumnegLL == min(sumnegLL))
}

