data {                          
  int<lower=0> N;       
  real<lower=0> sigma;     
  real<lower=0> alpha_mu00;                
  real<lower=0> alpha_mu01;                
  real<lower=0> alpha_mu10;                
  real<lower=0> alpha_mu11; 
  real<lower=0> beta0_mu;                
  real<lower=0> beta1_mu;                
  real<lower=0> theta0_mu;                
  real<lower=0> theta1_mu;                
  real<lower=0> gamma_mu;                
  

  vector<lower=0>[N] weights;
  int<lower=0,upper=1> y[N];  
  vector[N] xalpha00;          
  vector[N] xalpha01;            
  vector[N] xalpha10;                
  vector[N] xalpha11;               
  vector[N] xbeta0;                 
  vector[N] xbeta1;                
  vector[N] xtheta0;               
  vector[N] xtheta1;                 
  vector[N] xgamma;                 

  
}
parameters {
  real alpha00;                   
  real alpha01;                  
  real alpha10;                    
  real alpha11;                   
  real gamma; 
  real beta0;
  real beta1;
  real theta0;
  real theta1;
  
  
}

model {
  
  alpha00 ~ normal(alpha_mu00,sigma);        
  alpha01 ~ normal(alpha_mu01,sigma);        
  alpha10 ~ normal(alpha_mu10,sigma);        
  alpha11 ~ normal(alpha_mu11,sigma);   
  
  for (i in 1:N) {
    target += weights[i] * (bernoulli_logit_lpmf(y[i] | alpha00*xalpha00[i] + alpha01*xalpha01[i] + alpha10*xalpha10[i] + alpha11*xalpha11[i] + beta0*xbeta0[i] + beta1*xbeta1[i] + theta0*xtheta0[i] + theta1*xtheta1[i] + gamma*xgamma[i]));
  }
  
  beta0 ~ normal(beta0_mu,sigma);  
  beta1 ~ normal(beta1_mu,sigma);        
  theta0 ~ normal(theta0_mu,sigma); 
  theta1 ~ normal(theta1_mu,sigma); 
  gamma ~ normal(gamma_mu,sigma); 

}

generated quantities {         
  vector[16] DTR;
  
  DTR[1] = exp(alpha00 + beta0 + theta0 + gamma) / (1 + exp(alpha00 + beta0 + theta0 + gamma)); // AAC00
  DTR[2] = exp(alpha00 + beta0 - theta0 - gamma) / (1 + exp(alpha00 + beta0 - theta0 - gamma)); // AAD00
  DTR[3] = exp(alpha00 - beta0 + theta0 - gamma) / (1 + exp(alpha00 - beta0 + theta0 - gamma)); // BBC00
  DTR[4] = exp(alpha00 - beta0 - theta0 + gamma) / (1 + exp(alpha00 - beta0 - theta0 + gamma)); // BBD00
  DTR[5] = exp(alpha01 + beta0 + theta1 + gamma) / (1 + exp(alpha01 + beta0 + theta1 + gamma)); // AAC01
  DTR[6] = exp(alpha01 + beta0 - theta1 - gamma) / (1 + exp(alpha01 + beta0 - theta1 - gamma)); // AAD01
  DTR[7] = exp(alpha01 - beta0 + theta1 - gamma) / (1 + exp(alpha01 - beta0 + theta1 - gamma)); // BBC01
  DTR[8] = exp(alpha01 - beta0 - theta1 + gamma) / (1 + exp(alpha01 - beta0 - theta1 + gamma)); // BBD01
  DTR[9] = exp(alpha10 + beta1 + theta0 + gamma) / (1 + exp(alpha10 + beta1 + theta0 + gamma)); // AAC10
  DTR[10] = exp(alpha10 + beta1 - theta0 - gamma) / (1 + exp(alpha10 + beta1 - theta0 - gamma)); // AAD10
  DTR[11] = exp(alpha10 - beta1 + theta0 - gamma) / (1 + exp(alpha10 - beta1 + theta0 - gamma)); // BBC10
  DTR[12] = exp(alpha10 - beta1 - theta0 + gamma) / (1 + exp(alpha10 - beta1 - theta0 + gamma)); // BBD10
  DTR[13] = exp(alpha11 + beta1 + theta1 + gamma) / (1 + exp(alpha11 + beta1 + theta1 + gamma)); // AAC11
  DTR[14] = exp(alpha11 + beta1 - theta1 - gamma) / (1 + exp(alpha11 + beta1 - theta1 - gamma)); // AAD11
  DTR[15] = exp(alpha11 - beta1 + theta1 - gamma) / (1 + exp(alpha11 - beta1 + theta1 - gamma)); // BBC11
  DTR[16] = exp(alpha11 - beta1 - theta1 + gamma) / (1 + exp(alpha11 - beta1 - theta1 + gamma)); // BBD11
  
}
