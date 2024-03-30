data {                          
  int<lower=0> N;       
  real<lower=0> sigma;     
  real<lower=0> alpha_mu; 
  real<lower=0> beta_mu;                
  real<lower=0> theta_mu; 
  real<lower=0> gamma_mu;                

  vector<lower=0>[N] weights;
  int<lower=0,upper=1> y[N];  
  vector[N] xalpha;            
  vector[N] xbeta;                 
  vector[N] xtheta;                 
  vector[N] xgamma;                 

  
}
parameters {
  real alpha;                 
  real beta;                
  real theta;                    
  real gamma;                    
  
  
}

model {
  
  alpha ~ normal(alpha_mu,sigma);    
  beta ~ normal(beta_mu,sigma);        
  theta ~ normal(theta_mu,sigma);        


  for (i in 1:N) {
    target += weights[i] * (bernoulli_logit_lpmf(y[i] | alpha*xalpha[i] + beta*xbeta[i] + theta*xtheta[i] + gamma*xgamma[i]));
  }
  
  gamma ~ normal(gamma_mu,sigma);   
}

generated quantities {        
  vector[4] DTR;
  
  DTR[1] = exp(alpha + beta + theta + gamma) / (1 + exp(alpha + beta + theta + gamma)); // AAC00
  DTR[2] = exp(alpha + beta - theta - gamma) / (1 + exp(alpha + beta - theta - gamma)); // AAD00
  DTR[3] = exp(alpha - beta + theta - gamma) / (1 + exp(alpha - beta + theta - gamma)); // BBC00
  DTR[4] = exp(alpha - beta - theta + gamma) / (1 + exp(alpha - beta - theta + gamma)); // BBD00
  
}
