// This is pseudo Stan code for a
// Bayesian weighted and replicated regression model (WRRM)
// to estimate DTRs from a PRPP-SMART.
// Note: this Stan code is cosistent with
// PRPP-SMART data of the form described in
// Supplementary Material 4: Example R Code
// in the paper titled "A Partially Randomized Patient Preference, Sequential,
// Multiple-Assignment, Randomized Trial Design Analyzed via Weighted and
// Replicated Frequentist and Bayesian Methods"


// Specify data names. Note, you may name the variables anyway you wish but they must correspond to your R code where you set the data for Stan
data {                          
  int<lower=0> N;             // replicated PRPP-SMART data sample size      
  real<lower=0> sigma;        // sqrt(variance) hyperparameter for all priors
  real<lower=0> alpha_mu00;   // mu hyperparameter for alpha_00 prior            
  real<lower=0> alpha_mu01;   // mu hyperparameter for alpha_01 prior             
  real<lower=0> alpha_mu10;   // mu hyperparameter for alpha_10 prior                
  real<lower=0> alpha_mu11;   // mu hyperparameter for alpha_11 prior   
  real<lower=0> beta0_mu;     // mu hyperparameter for beta_0 prior            
  real<lower=0> beta1_mu;     // mu hyperparameter for beta_1 prior            
  real<lower=0> theta0_mu;    // mu hyperparameter for theta_0 prior              
  real<lower=0> theta1_mu;    // mu hyperparameter for theta_1 prior              
  real<lower=0> gamma_mu;     // mu hyperparameter for gamma prior              
  

  vector<lower=0>[N] weights; // PRPP-SMART weights 
  int<lower=0,upper=1> y[N];  // outcome variable Y
  vector[N] xalpha00;         // data for alpha_00 intercept: (1 - P1) * (1 - P2)     
  vector[N] xalpha01;         // data for alpha_01 intercept: (1 - P1) * P2   
  vector[N] xalpha10;         // data for alpha_10 intercept: P1 * (1 - P2)       
  vector[N] xalpha11;         // data for alpha_11 intercept: P1 * P2      
  vector[N] xbeta0;           // data for beta_0: T1*(1 - P1),       
  vector[N] xbeta1;           // data for beta_1: T1 * P1      
  vector[N] xtheta0;          // data for theta_0: T2 * (1 - P2)
  vector[N] xtheta1;          // data for theta_0: T2 * P2       
  vector[N] xgamma;           // data for gamma: T1 * T2     

  
}

// Specify parameters to put priors on
parameters {
  real alpha00; // parameter alpha_00                
  real alpha01; // parameter alpha_01                 
  real alpha10; // parameter alpha_10                   
  real alpha11; // parameter alpha_11                 
  real gamma;   // parameter gamma
  real beta0;   // parameter beta_0
  real beta1;   // parameter beta_1
  real theta0;  // parameter theta_0
  real theta1;  // parameter theta_1
  
  
}

// Specify priors and model 
model {
  
  // Priors
  alpha00 ~ normal(alpha_mu00,sigma);        
  alpha01 ~ normal(alpha_mu01,sigma);        
  alpha10 ~ normal(alpha_mu10,sigma);        
  alpha11 ~ normal(alpha_mu11,sigma);   
  beta0 ~ normal(beta0_mu,sigma);  
  beta1 ~ normal(beta1_mu,sigma);        
  theta0 ~ normal(theta0_mu,sigma); 
  theta1 ~ normal(theta1_mu,sigma); 
  gamma ~ normal(gamma_mu,sigma); 

  // Model - parameter names must be exactly the same as those given in parameter section above
  for (i in 1:N) {
    target += weights[i] * (bernoulli_logit_lpmf(y[i] | alpha00*xalpha00[i] + alpha01*xalpha01[i] + alpha10*xalpha10[i] + alpha11*xalpha11[i] + beta0*xbeta0[i] + beta1*xbeta1[i] + theta0*xtheta0[i] + theta1*xtheta1[i] + gamma*xgamma[i]));
  }
  


}

// Compute DTRs so they can be monitored for posterior draws - again must use same parameter names as given above
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
