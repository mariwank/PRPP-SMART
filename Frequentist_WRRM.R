# install/load libraries 
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
p_load(geepack, geometry, dplyr, MASS, data.table, broom, msm, Rlab, tidyverse, utils) 


### Data generation settings
# Load in data generation function
source("DataGeneration.R")

## USER SETTINGS:
# specify what scenario you want to run (1,2,3)
scenario <- 1

if (scenario == 1){
  pNP1=0.50 #  desired proportion of individuals expressing No Preference in stage 1
  pNP2=0.50 # desired proportion of patients expressing No Preference in stage 2 (among non-responders)

} else if (scenario == 2){
  pNP1=1/3
  pNP2=1/3
}  else if (scenario == 3){
  pNP1=0.5
  pNP2=1/3
}

# Specify number of subjects in trial
N <- 500

# Specify theta targets
pTheta_A=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
pTheta_C=0.4 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)

# Data gen parameters 
beta1_link <- c(1,1) # Beta1A (responders to A) Beta1B (resopnders to B)
alphaP_link <- c(1.1, 1.05) # alpha_p1 (1st stage preference), alpha_p2 (2nd stage preference)


# First stage treatment response rates
pi_A <- 0.6 # stage 1 response rate to randomize A: Pr(R=1|T1=A,P1=0)
pi_B <- 0.45 # stage 1 response rate to randomize B: Pr(R=1|T1=B,P1=0)
pi_A1 <- pi_A*alphaP_link[1] # stage 1 response rate to prefer A: Pr(R=1|T1=A,P1=1)
pi_B1 <- pi_B*alphaP_link[1] # stage 1 response rate to prefer B: Pr(R=1|T1=B,P1=1)

# Second stage treatment response rates
pi_AC <- 0.5   # Second stage response rate of non-responders to randomized A who receive randomized C in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=C)        
pi_AD <- 0.4   # Second stage response rate of non-responders to randomized A who receive randomized D in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=D)
pi_BC <- 0.3   # Second stage response rate of non-responders to randomized B who receive randomized C in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=C)
pi_BD <- 0.2   # Second stage response rate of non-responders to randomized B who receive randomized D in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=D)
pA0A <- pi_A * beta1_link[1]                         # Second stage response rate of responders to randomized A
pB0B <- pi_B * beta1_link[2]                         # Second stage response rate of responders to randomized B
pA1A <- pi_A * alphaP_link[1] * beta1_link[1]        # Second stage response rate of responders to preferred A
pB1B <- pi_B * alphaP_link[1] * beta1_link[2]        # Second stage response rate of responders to preferred B
pA0C1 <- alphaP_link[2] * pi_AC                      # Second stage response rate of non-responders to randomized A who receive preferred C in the second stage
pA0D1 <- alphaP_link[2] * pi_AD                      # Second stage response rate of non-responders to randomized A who receive preferred D in the second stage
pA1C0 <- alphaP_link[1] * pi_AC                      # Second stage response rate of non-responders to preferred A who receive randomized C in the second stage
pA1D0 <- alphaP_link[1] * pi_AD                      # Second stage response rate of non-responders to preferred A who receive randomized D in the second stage
pA1C1 <- alphaP_link[1] * alphaP_link[2] * pi_AC     # Second stage response rate of non-responders to preferred A who receive preferred C in the second stage
pA1D1 <- alphaP_link[1] * alphaP_link[2] * pi_AD     # Second stage response rate of non-responders to preferred A who receive preferred D in the second stage
pB0C1 <- alphaP_link[2] * pi_BC                      # Second stage response rate of non-responders to randomized B who receive preferred C in the second stage
pB0D1 <- alphaP_link[2] * pi_BD                      # Second stage response rate of non-responders to randomized B who receive preferred D in the second stage
pB1C0 <- alphaP_link[1] * pi_BC                      # Second stage response rate of non-responders to preferred B who receive randomized C in the second stage
pB1D0 <- alphaP_link[1] * pi_BD                      # Second stage response rate of non-responders to preferred B who receive randomized D in the second stage
pB1C1 <- alphaP_link[1] * alphaP_link[2] * pi_BC     # Second stage response rate of non-responders to preferred B who receive preferred C in the second stage
pB1D1 <- alphaP_link[1] * alphaP_link[2] * pi_BD     # Second stage response rate of non-responders to preferred B who receive preferred D in the second stage


### TRUE DTRS ###

expected_pref <- c()  # expected DTR response rates from our simulated data
expected_pref[1] <- pi_A * pA0A + (1 - pi_A) * pi_AC  #AAC00
expected_pref[2] <- pi_A * pA0A + (1 - pi_A) * pi_AD  #AAD00
expected_pref[3] <- pi_B * pB0B + (1 - pi_B) * pi_BC  #BBC00
expected_pref[4] <- pi_B * pB0B + (1 - pi_B) * pi_BD  #BBD00
expected_pref[5] <- pi_A * pA0A + (1 - pi_A) * pA0C1 #AAC01
expected_pref[6] <- pi_A * pA0A + (1 - pi_A) * pA0D1 #AAD01
expected_pref[7] <- pi_B * pB0B + (1 - pi_B) * pB0C1 #BBC01
expected_pref[8] <- pi_B * pB0B + (1 - pi_B) * pB0D1 #BBD01
expected_pref[9] <- pi_A1 * pA1A + (1 - pi_A1) * pA1C0 #AAC10
expected_pref[10] <- pi_A1 * pA1A + (1 - pi_A1) * pA1D0 #AAD10
expected_pref[11] <- pi_B1 * pB1B + (1 - pi_B1) * pB1C0 #BBC10
expected_pref[12] <- pi_B1 * pB1B + (1 - pi_B1) * pB1D0 #BBD10
expected_pref[13] <- pi_A1 * pA1A + (1 - pi_A1) * pA1C1 #AAC11
expected_pref[14] <- pi_A1 * pA1A + (1 - pi_A1) * pA1D1 #AAD11
expected_pref[15] <- pi_B1 * pB1B + (1 - pi_B1) * pB1C1 #BBC11
expected_pref[16] <- pi_B1 * pB1B + (1 - pi_B1) * pB1D1 #BBD11

# Define contrast matrix for DTR calculation
contrast_dtr <- matrix(c(1,0,0,0,1,0,1,0,1, #AC00
                         1,0,0,0,1,0,-1,0,-1, #AD00
                         1,0,0,0,-1,0,1,0,-1, #BC00
                         1,0,0,0,-1,0,-1,0,1, #BD00
                         0,1,0,0,1,0,0,1,1, #AC01
                         0,1,0,0,1,0,0,-1,-1, #AD01
                         0,1,0,0,-1,0,0,1,-1, #BC01
                         0,1,0,0,-1,0,0,-1,1, #BD01
                         0,0,1,0,0,1,1,0,1, #AC10
                         0,0,1,0,0,1,-1,0,-1, #AD10
                         0,0,1,0,0,-1,1,0,-1, #BC10
                         0,0,1,0,0,-1,-1,0,1, #BD10
                         0,0,0,1,0,1,0,1,1, #AC11
                         0,0,0,1,0,1,0,-1,-1, #AD11
                         0,0,0,1,0,-1,0,1,-1, #BC11
                         0,0,0,1,0,-1,0,-1,1 #BD11
),nrow = 16, ncol=9, byrow = TRUE)

true_DTR_mat <- matrix(c(expected_pref), ncol = 1)
rownames(true_DTR_mat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")

## TRUE DTRs TRADITIONAL ##

# Define contrast matrix for DTR calculation
contrast_dtr_traditional <- matrix(c(1,1,1,1, #AC00
                                     1,1,-1,-1, #AD00
                                     1,-1,1,-1, #BC00
                                     1,-1,-1,1 #BD00
),nrow = 4, ncol=4, byrow = TRUE)

true_DTR_mat_traditional <- matrix(expected_pref[1:4], ncol = 1)
rownames(true_DTR_mat_traditional) <- c("AAC00", "AAD00", "BBC00", "BBD00")

### SIMULATION ###
# find number of sims needed to get 500 given the above settings 
source("nsim_toget_500.R")
iterations_needed <- count_iterations()
print(iterations_needed)
n.sim <- iterations_needed
num_skip <- rep(0, n.sim) # number simulations skipped

## Store Full WRRM PRPP-SMART results 
DTR_hat <- matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation 
parameter_hat <- matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter estimates per simulation 
variance_param_hat <- matrix(NA, nrow=9, ncol=n.sim) # matrix to store sandwich variance of parameter estimates per simulation
variance_dtr_hat <- matrix(NA, nrow=16, ncol=n.sim) # matrix to store delta method DTR variance estimates per simulation
ci_hat <- matrix(NA, nrow=16, ncol=n.sim) # matrix to store whether ci covers truth per simulation  
n.DTR <- matrix(NA,nrow = 16, ncol = n.sim) # matrix to store sample size for each DTR path per simulation

## Store Traditional WRRM analysis results 
DTR_hat_t <- matrix(NA,nrow = 4, ncol = n.sim) # matrix to store preference DTR estimate per simulation 
parameter_hat_t <- matrix(NA, nrow=4, ncol=n.sim) # matrix to store parameter estimates per simulation 
variance_param_hat_t <- matrix(NA, nrow=4, ncol=n.sim) # matrix to store sandwich variance of parameter estimates per simulation
variance_dtr_hat_t <- matrix(NA, nrow=4, ncol=n.sim) # matrix to store delta method DTR variance estimates per simulation
ci_hat_t <- matrix(NA, nrow=4, ncol=n.sim) # matrix to store whether ci covers truth per simulation  


## N DTR track
fullrep_n <- rep(NA, n.sim) # full analysis replicated data n
#00
n_A0A <- rep(NA, n.sim)  # treatment pathway n: randomized A responded
n_B0B <- rep(NA, n.sim)  # treatment pathway n: randomized B responded
n_A0C0 <- rep(NA, n.sim) # treatment pathway n: randomized A no response randomized C
n_A0D0 <- rep(NA, n.sim) # treatment pathway n: randomized A no response randomized D
n_B0C0 <- rep(NA, n.sim) # treatment pathway n: randomized B no response randomized C
n_B0D0 <- rep(NA, n.sim) # treatment pathway n: randomized B no response randomized D
#01
n_A0C1 <- rep(NA, n.sim) # treatment pathway n: randomized A no response preferred C
n_A0D1 <- rep(NA, n.sim) # treatment pathway n: randomized A no response preferred D
n_B0C1 <- rep(NA, n.sim) # treatment pathway n: randomized B no response preferred C
n_B0D1 <- rep(NA, n.sim) # treatment pathway n: randomized B no response preferred D
#10
n_A1A <- rep(NA, n.sim)  # treatment pathway n: preferred A responded
n_B1B <- rep(NA, n.sim)  # treatment pathway n: preferred B responded
n_A1C0 <- rep(NA, n.sim) # treatment pathway n: preferred A no response randomized C
n_A1D0 <- rep(NA, n.sim) # treatment pathway n: preferred A no response randomized D
n_B1C0 <- rep(NA, n.sim) # treatment pathway n: preferred B no response randomized C
n_B1D0 <- rep(NA, n.sim) # treatment pathway n: preferred B no response randomized D
#11
n_A1C1 <- rep(NA, n.sim) # treatment pathway n: preferred A no response preferred C
n_A1D1 <- rep(NA, n.sim) # treatment pathway n: preferred A no response preferred D
n_B1C1 <- rep(NA, n.sim) # treatment pathway n: preferred B no response preferred C
n_B1D1 <- rep(NA, n.sim) # treatment pathway n: preferred B no response preferred D

# N DTR traditional track
trad_n <- rep(NA, n.sim)       # traditional analysis n
tradrep_n <- rep(NA, n.sim)    # traditional analysis replicated data n
trad_n_A0A <- rep(NA, n.sim)   # treatment pathway n: traditional prpp-smart analysis randomized A responded
trad_n_A0C0 <- rep(NA, n.sim)  # treatment pathway n: traditional prpp-smart analysis randomized A non-response randomized C
trad_n_A0D0 <- rep(NA, n.sim)  # treatment pathway n: traditional prpp-smart analysis randomized A non-response randomized D
trad_n_B0B <- rep(NA, n.sim)   # treatment pathway n: traditional prpp-smart analysis randomized B responded
trad_n_B0C0 <- rep(NA, n.sim)  # treatment pathway n: traditional prpp-smart analysis randomized B non-response randomized C
trad_n_B0D0 <- rep(NA, n.sim)  # treatment pathway n: traditional prpp-smart analysis randomized B non-response randomized D

# simulation start
for (i in 1:n.sim){
  set.seed(i+1000)
  
  # Generate data
  data <- generate_data(N=N, pNP1=pNP1, pTheta_A=pTheta_A, pNP2=pNP2, pTheta_C=pTheta_C, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
  df <- data[[2]] # data for prpp_smart full analysis replicated and has calculated weights
  
  df <-df[order(df$id),] # sort data by id for gee statement

  ## create traditional analysis data ## 
  
  # pull out randomized only subjects
  df2 <- data[[1]] # relabeled raw data (not replicated but has calculated weights)
  rand.ind <- which(df2$Treatment_Path == "A0A" | df2$Treatment_Path == "A0C0" | df2$Treatment_Path == "A0D0"| df2$Treatment_Path == "B0B" | df2$Treatment_Path == "B0C0" | df2$Treatment_Path == "B0D0") # only want traditional randomized SMART sujbects
  df_rand <- df2[rand.ind,]
  
  # replicate data - only twice for traditional analysis 
  # 1st dataset of responders setting T2=1, 
  datareps1 <- df_rand[df_rand$R==1,] 
  datareps1$T2 <- 1
  
  # 2nd dataset of responders setting T2=1, 
  datareps2 <- df_rand[df_rand$R==1,] 
  datareps2$T2 <- -1
  
  # dataset for non-responders
  datanoresp <- df_rand[df_rand$R==0,]
  
  # replicated data
  replicated_dat_trad <- rbind(datareps1,datareps2,datanoresp)
  
  # create data used in traditional analysis 
  analysis_data_trad <- replicated_dat_trad %>% 
    dplyr::select(id, Y, w, T1, T2, Treatment_Path)
  
  analysis_data_trad <-analysis_data_trad[order(analysis_data_trad$id),] # sort data by id for gee statement
  
  # check to make sure at least three subjects per treatment path in PRPP-SMART data if not skip simulation
  trialpath_df <- data[[1]] %>% dplyr::group_by(Treatment_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
  if (nrow(trialpath_df) < 20) {num_skip[i]=1 ; next} # check to make sure data in each pathway

  # DTR N Track
  fullrep_n[i] <- nrow(df) # n (replicated) full prpp-smart data
  n.DTR[,i] = data[[3]]
  #00
  n_A0A[i] <- sum(data[[1]]$Treatment_Path == "A0A")
  n_B0B[i] <- sum(data[[1]]$Treatment_Path == "B0B")
  n_A0C0[i] <- sum(data[[1]]$Treatment_Path == "A0C0")
  n_A0D0[i] <- sum(data[[1]]$Treatment_Path == "A0D0")
  n_B0C0[i] <- sum(data[[1]]$Treatment_Path == "B0C0")
  n_B0D0[i] <- sum(data[[1]]$Treatment_Path == "B0D0")
  #01
  n_A0C1[i] <- sum(data[[1]]$Treatment_Path == "A0C1")
  n_A0D1[i] <- sum(data[[1]]$Treatment_Path == "A0D1")
  n_B0C1[i] <- sum(data[[1]]$Treatment_Path == "B0C1")
  n_B0D1[i] <- sum(data[[1]]$Treatment_Path == "B0D1")
  #10
  n_A1A[i] <- sum(data[[1]]$Treatment_Path == "A1A")
  n_B1B[i] <- sum(data[[1]]$Treatment_Path == "B1B")
  n_A1C0[i] <- sum(data[[1]]$Treatment_Path == "A1C0")
  n_A1D0[i] <- sum(data[[1]]$Treatment_Path == "A1D0")
  n_B1C0[i] <- sum(data[[1]]$Treatment_Path == "B1C0")
  n_B1D0[i] <- sum(data[[1]]$Treatment_Path == "B1D0")
  #11
  n_A1C1[i] <- sum(data[[1]]$Treatment_Path == "A1C1")
  n_A1D1[i] <- sum(data[[1]]$Treatment_Path == "A1D1")
  n_B1C1[i] <- sum(data[[1]]$Treatment_Path == "B1C1")
  n_B1D1[i] <- sum(data[[1]]$Treatment_Path == "B1D1")
  
  # DTR N Track Traditional
  trad_n[i] <- nrow(df_rand) # n (non-replicated) randomized randomized data
  tradrep_n[i] <- nrow(analysis_data_trad)
  trad_n_A0A[i] <- sum(df_rand$Treatment_Path == "A0A")
  trad_n_A0C0[i] <- sum(df_rand$Treatment_Path == "A0C0")
  trad_n_A0D0[i] <- sum(df_rand$Treatment_Path == "A0D0")
  trad_n_B0B[i] <- sum(df_rand$Treatment_Path == "B0B")
  trad_n_B0C0[i] <- sum(df_rand$Treatment_Path == "B0C0")
  trad_n_B0D0[i] <- sum(df_rand$Treatment_Path == "B0D0")
  
  ### FULL WRRM PRPP-SMART ANALYSIS ###
  
  # fit 
  gee.fit <- geeglm(Y ~ xalpha00 + xalpha01 + xalpha10 + xalpha11 + xbeta0 + xbeta1 + xtheta0 + xtheta1 + xgamma-1, id=id, weights=w, family = "binomial", corstr = "independence", data = df) # warning is normal since using non-integer weights
  
  ## parameter estimates
  parameter_hat[,i] <- gee.fit$coefficients
  rownames(parameter_hat) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1","theta0", "theta1","gamma")
  
  ## DTR estimates 
  for(j in 1:16){
    DTR_hat[j,i] = exp((contrast_dtr[j,]%*%parameter_hat[,i]))/(1 + exp((contrast_dtr[j,]%*%parameter_hat[,i])))
    
  }
  rownames(DTR_hat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
  
  
  ## Parameter Variances 
  
  # robust variances of each parameter 
  variance_param_hat[,i] <- diag(gee.fit$geese$vbeta)
  rownames(variance_param_hat) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1","theta0", "theta1","gamma")
  
  ## Delta method for DTR variances
  
  a00 <- parameter_hat[1,i]
  a01 <- parameter_hat[2,i]
  a10 <- parameter_hat[3,i]
  a11 <- parameter_hat[4,i]
  b0 <- parameter_hat[5,i]
  b1 <- parameter_hat[6,i]
  t0 <- parameter_hat[7,i]
  t1 <- parameter_hat[8,i]
  g <- parameter_hat[9,i]
  
  # Extract variance covariance matrix from geeglm
  varcov <- gee.fit$geese$vbeta
  
  # delta method for dtr variances
  for (d in 1:16){
    
    if (d == 1){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
    }
    
    if (d == 2){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
      
    }
    
    if (d == 3){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
      
    }
    
    if (d == 4){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
      
    }
    
    if (d == 5){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 6){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 7){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 8){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 9){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 10){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 11){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 12){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 13){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
    
    if (d == 14){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
    
    if (d == 15){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
    
    if (d == 16){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
  }
  
  rownames(variance_dtr_hat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
  
  
  ## calculate CI coverage of DTRs

  lwr_dtr <- c()
  upper_dtr <- c()
  
  for(d in 1:16){
    lwr_dtr[d] <- DTR_hat[d,i] - qnorm((1+0.95)/2)*sqrt(variance_dtr_hat[d,i])
    upper_dtr[d] <- DTR_hat[d,i] + qnorm((1+0.95)/2)*sqrt(variance_dtr_hat[d,i])
    
  }
  dtr_ci_mat <- cbind(lwr_dtr, upper_dtr)
  rownames(dtr_ci_mat) <-  c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")

  param_in_ci <- data.table::between(true_DTR_mat[,1], dtr_ci_mat[, 1], dtr_ci_mat[, 2])
  ci_hat[,i] <- param_in_ci
  rownames(ci_hat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
  
  
  ### Traditional PRPP-SMART ANALYSIS ###
  
  # fit 
  gee.fit.t <- geeglm(Y ~ T1 + T2 + T1:T2, id=id, weights=w, family = "binomial", corstr = "independence", data = analysis_data_trad) # warning is normal since using non-integer weights
  
  
  ## parameter estimates
  parameter_hat_t[,i] <- gee.fit.t$coefficients
  rownames(parameter_hat_t) <- c("alpha", "beta", "theta", "gamma")
  
  ## DTR estimates 
  for(j in 1:4){
    DTR_hat_t[j,i] = exp((contrast_dtr_traditional[j,]%*%parameter_hat_t[,i]))/(1 + exp((contrast_dtr_traditional[j,]%*%parameter_hat_t[,i])))
    
  }
  rownames(DTR_hat_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
  
  ## Parameter Variances 
  
  # robust variances of each parameter 
  variance_param_hat_t[,i] <- diag(gee.fit.t$geese$vbeta)
  rownames(variance_param_hat_t) <- c("alpha", "beta", "theta", "gamma")
  
  ## Delta method for DTR variances
  
  a <- parameter_hat_t[1,i]
  b <- parameter_hat_t[2,i]
  t <- parameter_hat_t[3,i]
  g <- parameter_hat_t[4,i]
  
  
  # Extract variance covariance matrix from geeglm
  varcov_t <- gee.fit.t$geese$vbeta 
  
  # delta method for dtr variances
  for (d in 1:4){
    
    if (d == 1){
      variance_dtr_hat_t[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a, b, t, g), varcov_t))^2
    }
    
    if (d == 2){
      variance_dtr_hat_t[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a, b, t, g), varcov_t))^2
      
    }
    
    if (d == 3){
      variance_dtr_hat_t[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a, b, t, g), varcov_t))^2
      
    }
    
    if (d == 4){
      variance_dtr_hat_t[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a, b, t, g), varcov_t))^2 # deltamethod returns standard error
      
    }
    
  }
  rownames(variance_dtr_hat_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
  
  ## calculate DTR CI coverage 
  
  lwr_dtr_t <- c()
  upper_dtr_t <- c()
  
  for(d in 1:4){
    lwr_dtr_t[d] <- DTR_hat_t[d,i] - qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_t[d,i])
    upper_dtr_t[d] <- DTR_hat_t[d,i] + qnorm((1+0.95)/2)*sqrt(variance_dtr_hat_t[d,i])
    
  }
  dtr_ci_mat_t <- cbind(lwr_dtr_t, upper_dtr_t)
  rownames(dtr_ci_mat_t) <-  c("AAC00", "AAD00", "BBC00", "BBD00")

  param_in_ci_t <- data.table::between(true_DTR_mat_traditional[,1], dtr_ci_mat_t[, 1], dtr_ci_mat_t[, 2])
  ci_hat_t[,i] <- param_in_ci_t
  rownames(ci_hat_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
  
}

######################### EVALUATION ############################

### Save Settings
# set date and time for file saving 
st<-format(Sys.time(), "%Y_%m_%d_%H_%M")

# Define the folder path where you want your results saved
folder_path <- "path to your folder where you want to store frequentist results"

# Create the folder if it doesn't exist
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created at", folder_path, "\n")
} else {
  cat("Folder already exists at", folder_path, "\n")
}

### DTR Full Model Results 
DTR = c("AAC00","AAD00","BBC00","BBD00","AAC01","AAD01","BBC01","BBD01",
        "AAC10","AAD10","BBC10","BBD10","AAC11","AAD11","BBC11","BBD11")
DTR_hat_avg = c()
DTR_sd_hat = c()
DTR_avg_sd = c()
DTR_avg_n = c()

for(i in 1:16){
  DTR_hat_avg[i] <- round(mean(DTR_hat[i,], na.rm = TRUE),4)
  DTR_sd_hat[i] <- round(sd(DTR_hat[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 DTR estimates)
  DTR_avg_sd[i] <- round(sqrt(mean(variance_dtr_hat[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(dtr hat))) equivalent to mean(sqrt(var(dtr hat)))
  DTR_avg_n[i] <- round(mean(n.DTR[i,], na.rm = TRUE),1)
  
}

# calculate bias
DTR_bias = round(DTR_hat_avg - expected_pref,4)

# calculate rMSE
rMSE_DTR <- sqrt(DTR_sd_hat^2 + DTR_bias^2)


DTR_results_tbl <- data.frame(DTR=DTR,True_DTR=expected_pref,
                              DTR_Hat_Avg=DTR_hat_avg, 
                              Bias=DTR_bias, 
                              SE = DTR_sd_hat, 
                              Avg_se = DTR_avg_sd,
                              rMSE = rMSE_DTR) 

# Define the file name
file_name <- paste0("DTR_Results_Full_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(DTR_results_tbl, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

## nominal coverage 
ci_coverage = c()

for(i in 1:nrow(ci_hat)){
  ci_coverage[i] <- mean(ci_hat[i,], na.rm = TRUE)
  
}

ci_results_tbl <- data.frame(DTR, 
                             CI_Coverage=ci_coverage)

# Define the file name
file_name <- paste0("DTR_CI_Results_Full_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(ci_results_tbl, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

## Pathway N
#00
n_A0A_avg <- mean(n_A0A, na.rm=T)
n_B0B_avg <- mean(n_B0B, na.rm=T)
n_A0C0_avg <- mean(n_A0C0, na.rm=T)
n_A0D0_avg <- mean(n_A0D0, na.rm=T)
n_B0C0_avg <- mean(n_B0C0, na.rm=T)
n_B0D0_avg <- mean(n_B0D0, na.rm=T)
#01
n_A0C1_avg <- mean(n_A0C1, na.rm=T)
n_A0D1_avg <- mean(n_A0D1, na.rm=T)
n_B0C1_avg <- mean(n_B0C1, na.rm=T)
n_B0D1_avg <- mean(n_B0D1, na.rm=T)
#10
n_A1A_avg <- mean(n_A1A, na.rm=T)
n_B1B_avg <- mean(n_B1B, na.rm=T)
n_A1C0_avg <- mean(n_A1C0, na.rm=T)
n_A1D0_avg <- mean(n_A1D0, na.rm=T)
n_B1C0_avg <- mean(n_B1C0, na.rm=T)
n_B1D0_avg <- mean(n_B1D0, na.rm=T)
#11
n_A1C1_avg <- mean(n_A1C1, na.rm=T)
n_A1D1_avg <- mean(n_A1D1, na.rm=T)
n_B1C1_avg <- mean(n_B1C1, na.rm=T)
n_B1D1_avg <- mean(n_B1D1, na.rm=T)

pathway_n_dtr_df <- data.frame(Pathway = c("A0A",
                                           "B0B",
                                           "A1A",
                                           "B1B",
                                           "A0C0", 
                                           "A0D0", 
                                           "B0C0",
                                           "B0D0",
                                           "A0C1", 
                                           "A0D1", 
                                           "B0C1",
                                           "B0D1",
                                           "A1C0", 
                                           "A1D0", 
                                           "B1C0",
                                           "B1D0",
                                           "A1C1", 
                                           "A1D1", 
                                           "B1C1",
                                           "B1D1",
                                           "Traditional_N"),
                               Avg_N = c(n_A0A_avg,
                                         n_B0B_avg,
                                         n_A1A_avg,
                                         n_B1B_avg,
                                         n_A0C0_avg,
                                         n_A0D0_avg,
                                         n_B0C0_avg,
                                         n_B0D0_avg,
                                         n_A0C1_avg,
                                         n_A0D1_avg,
                                         n_B0C1_avg,
                                         n_B0D1_avg,
                                         n_A1C0_avg,
                                         n_A1D0_avg,
                                         n_B1C0_avg,
                                         n_B1D0_avg,
                                         n_A1C1_avg,
                                         n_A1D1_avg,
                                         n_B1C1_avg,
                                         n_B1D1_avg,
                                         mean(trad_n, na.rm = TRUE)))

# Define the file name
file_name <- paste0("Treatment_Path_N_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(pathway_n_dtr_df, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

## DTR N

dtr_avgn_df <- data.frame(DTR, DTR_avg_n)

# Define the file name
file_name <- paste0("DTR_N_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(dtr_avgn_df, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")


## Parameter Results 

## Parameters
Parameter = c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1","theta0", "theta1","gamma")
param_hat_avg = c()
param_sd_hat = c()
param_avg_sd = c()

for(i in 1:length(Parameter)){
  param_hat_avg[i] <- round(mean(parameter_hat[i,], na.rm = TRUE),4)
  param_sd_hat[i] <- round(sd(parameter_hat[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 parameter estimates)
  param_avg_sd[i] <- round(sqrt(mean(variance_param_hat[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(parameter hat))) equivalent to mean(sqrt(var(parameter hat)))
  
}

param_results_tbl <- data.frame(Parameter=Parameter, 
                                Param_Hat_Avg=param_hat_avg, 
                                SE = param_sd_hat, # sd of 500 parameter estimates
                                Avg_se = param_avg_sd) # mean of 500 parameter variances then square root


# Define the file name
file_name <- paste0("Parameter_Results_Full_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(param_results_tbl, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

### DTR TRADITIONAL Results

DTR_t = c("AAC00","AAD00","BBC00","BBD00") # DTRs from traditional randomized randomized smart
DTR_hat_avg_t = c()
DTR_sd_hat_t = c()
DTR_avg_sd_t = c()

for(i in 1:4){
  DTR_hat_avg_t[i] <- round(mean(DTR_hat_t[i,], na.rm = TRUE),4)
  DTR_sd_hat_t[i] <- round(sd(DTR_hat_t[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 DTR estimates)
  DTR_avg_sd_t[i] <- round(sqrt(mean(variance_dtr_hat_t[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(dtr hat))) equivalent to mean(sqrt(var(dtr hat)))
}


# calculate bias
DTR_bias_t = round(DTR_hat_avg_t - expected_pref[1:4],4)

# calculate rMSE
rMSE_DTR_t <- sqrt(DTR_sd_hat_t^2 + DTR_bias_t^2)

DTR_results_tbl_t <- data.frame(DTR=DTR_t,
                                True_DTR=expected_pref[1:4],
                                DTR_Hat_Avg=DTR_hat_avg_t, 
                                Bias=DTR_bias_t, 
                                SE = DTR_sd_hat_t, 
                                Avg_se = DTR_avg_sd_t, 
                                rMSE = rMSE_DTR_t) 

# Define the file name
file_name <- paste0("DTR_Results_Traditional_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(DTR_results_tbl_t, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

## nominal coverage 
ci_coverage_t = c()

for(i in 1:nrow(ci_hat_t)){
  ci_coverage_t[i] <- mean(ci_hat_t[i,], na.rm = TRUE)
  
}

ci_results_tbl_t <- data.frame(DTR=DTR_t, 
                               CI_Coverage=ci_coverage_t)


# Define the file name
file_name <- paste0("DTR_CI_Results_Traditional_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(ci_results_tbl_t, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")



## Parameter Traditional 

## Parameters
Parameter_t = c("alpha", "beta", "theta", "gamma")
param_hat_avg_t = c()
param_sd_hat_t = c()
param_avg_sd_t = c()

for(i in 1:length(Parameter_t)){
  param_hat_avg_t[i] <- round(mean(parameter_hat_t[i,], na.rm = TRUE),4)
  param_sd_hat_t[i] <- round(sd(parameter_hat_t[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 parameter estimates)
  param_avg_sd_t[i] <- round(sqrt(mean(variance_param_hat_t[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(parameter hat))) equivalent to mean(sqrt(var(parameter hat)))
  
}

param_results_tbl_t <- data.frame(Parameter=Parameter_t, 
                                  Param_Hat_Avg=param_hat_avg_t, 
                                  SE = param_sd_hat_t, # sd of 500 parameter estimates
                                  Avg_se = param_avg_sd_t) # mean of 500 parameter variances then square root  


# Define the file name
file_name <- paste0("Parameter_Results_Traditional_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(param_results_tbl_t, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")


## number simulations skipped due to positivity assumption
num_skip_total <- sum(num_skip)

mean_skip <- mean(num_skip)
skip_df <- data.frame(num_sim = n.sim, 
                      num_skip_total = num_skip_total, 
                      avg_num_skip = mean_skip) 

# Define the file name
file_name <- paste0("Num_Sim_Skip_Summary_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(skip_df, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")


### EVALUATION - COMPARE FULL VS TRADITIONAL RESULTS ###

## DTR compare
DTR_results_tbl_c <- data.frame(DTR=DTR_t,
                                True_DTR=expected_pref[1:4],
                                DTR_Hat_Avg_all = DTR_hat_avg[1:4], # full prpp-smart analysis (all subjects)
                                DTR_Hat_Avg_trad = DTR_hat_avg_t, # traditional analysis (randomized randomized subjects only)
                                Bias_all = DTR_bias[1:4], # full prpp-smart analysis (all subjects)
                                Bias_trad = DTR_bias_t, # traditional analysis (randomized randomized subjects only)
                                Sd_all = DTR_sd_hat[1:4], 
                                Sd_trad = DTR_sd_hat_t, 
                                Avg_sd_all = DTR_avg_sd[1:4], 
                                Avg_sd_trad = DTR_avg_sd_t, 
                                rMSE_all = rMSE_DTR[1:4], 
                                rMSE_trad = rMSE_DTR_t) 


# Define the file name
file_name <- paste0("DTR_CompareResults_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(DTR_results_tbl_c, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

