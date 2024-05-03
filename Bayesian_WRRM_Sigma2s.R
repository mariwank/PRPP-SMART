# install/load libraries 
library(pacman)
p_load(survey, rstan, tidyverse, utils, matrixStats, Rlab) 

### STAN FUNCTIONS/CODE
# Cluster 
source("cs_samp_func.R") 
source("cs_samp_brms_func.R") 
source("grad_par_func.R") 
source("row_subset_func.R") 
source("DEadj_func.R") 

# load stan models
# Cluster
mod_dm <- stan_model('bayes_prpp_smart.stan', auto_write = TRUE) # prpp-smart stan
mod_dm_traditional <- stan_model('bayes_prpp_smart_traditional.stan', auto_write = TRUE) # traditional randomzied randomized only SMART from prpp-smart

### STAN SETTINGS

MCMC_SAMPLE <- 5500
BURN.IN <- 500
n_MCMC_chain <- 1
n_MCMC_thin <- 1

## Prior distributions hyperparameter specifications 
alpha_mu00 <- 0 
alpha_mu01 <- 0
alpha_mu10 <- 0
alpha_mu11 <- 0
beta0_mu <- 0
beta1_mu <- 0
theta0_mu <- 0
theta1_mu <- 0
gamma_mu <- 0 # used for traditional and full prpp-smart analysis 

alpha_mu <- 0 # prior on alpha for traditional analysis 
beta_mu <- 0 # prior on beta for traditional analysis 
theta_mu <- 0 # prior on theta for traditional analysis

### Data generation settings
# Load in data generation function
source("DataGeneration.R")

## USER SETTINGS:
# specify what scenario you want to run (1,2,3)
scenario = 1

if (scenario == 1){
  pNP1=0.50 #  desired proportion of individuals expressing No Preference in stage 1
  pNP2=0.50 # desired proportion of patients expressing No Preference in stage 2 (among non-responders)
  n.sim <- 527 # number of simulations to get 500
  
} else if (scenario == 2){
  pNP1=1/3
  pNP2=1/3
  n.sim <- 734 # number of simulations to get 500
  
}  else if (scenario == 3){
  pNP1=0.5
  pNP2=1/3
  n.sim <- 553 # number of simulations to get 500
  
} 

# Specify number of subjects in trial
N=500

# Specify theta targets
pTheta_A=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
pTheta_C=0.4 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)

# Data gen parameters 
beta1_link <- c(1,1) # Beta1A (responders to A) Beta1B (resopnders to B)
alphaP_link <- c(1.1, 1.05) # alpha_p1 (1st stage preference), alpha_p2 (2nd stage preference)


# First stage treatment response rates
pi_A=0.6 # stage 1 response rate to randomize A: Pr(R=1|T1=A,P1=0)
pi_B=0.45 # stage 1 response rate to randomize B: Pr(R=1|T1=B,P1=0)
pi_A1=pi_A*alphaP_link[1] # stage 1 response rate to prefer A: Pr(R=1|T1=A,P1=1)
pi_B1=pi_B*alphaP_link[1] # stage 1 response rate to prefer B: Pr(R=1|T1=B,P1=1)

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
sigma2_vec <- c(0.1, 0.5, 0.9) # different values of sigma2 to try
for (s in 1:length(sigma2_vec)) {
  
  sigma2 <- sigma2_vec[s]
  sigma <- sqrt(sigma2) 
  
  num_skip <- rep(0, n.sim) # number simulations skipped
  
  ## Store Full WRRM PRPP-SMART results 
  n.DTR <- matrix(NA,nrow = 16, ncol = n.sim) # matrix to store sample size for each DTR path per simulation
  DTR_hat_org = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation from original MCMC sample
  DTR_hat_adj = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation from adjusted MCMC sample
  DTR_var_org = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR variance estimate per simulation from original MCMC sample
  DTR_var_adj = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR variance estimate per simulation from adjusted MCMC sample
  parameter_hat_org = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter estimates per simulation from original MCMC sample
  parameter_hat_adj = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter estimates per simulation from adjusted MCMC sample
  parameter_var_org = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter variance estimates per simulation from original MCMC sample
  parameter_var_adj = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter variance estimates per simulation from adjusted MCMC sample
  ci_hat_org = matrix(NA, nrow=16, ncol=n.sim) # matrix to store whether ci covers truth per simulation from original MCMC sample 
  ci_hat_adj = matrix(NA, nrow=16, ncol=n.sim) # matrix to store whether ci covers truth per simulation from adjusted MCMC sample
  
  ## Store Traditional WRRM analysis results 
  DTR_hat_org_t = matrix(NA,nrow = 4, ncol = n.sim) # matrix to store preference DTR estimate per simulation from original MCMC sample
  DTR_hat_adj_t = matrix(NA,nrow = 4, ncol = n.sim) # matrix to store preference DTR estimate per simulation from adjusted MCMC sample
  DTR_var_org_t = matrix(NA,nrow = 4, ncol = n.sim) # matrix to store preference DTR variance estimate per simulation from original MCMC sample
  DTR_var_adj_t = matrix(NA,nrow = 4, ncol = n.sim) # matrix to store preference DTR variance estimate per simulation from adjusted MCMC sample
  parameter_hat_org_t = matrix(NA, nrow=4, ncol=n.sim) # matrix to store parameter estimates per simulation from original MCMC sample
  parameter_hat_adj_t = matrix(NA, nrow=4, ncol=n.sim) # matrix to store parameter estimates per simulation from adjusted MCMC sample
  parameter_var_org_t = matrix(NA, nrow=4, ncol=n.sim) # matrix to store parameter variance estimates per simulation from original MCMC sample
  parameter_var_adj_t = matrix(NA, nrow=4, ncol=n.sim) # matrix to store parameter variance estimates per simulation from adjusted MCMC sample
  ci_hat_org_t = matrix(NA, nrow=4, ncol=n.sim) # matrix to store whether ci covers truth per simulation from original MCMC sample 
  ci_hat_adj_t = matrix(NA, nrow=4, ncol=n.sim) # matrix to store whether ci covers truth per simulation from adjusted MCMC sample
  
  
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
  
  for (i in 1:n.sim){
    set.seed(i+100000)
    
    # Generate data
    data <- generate_data(N=N, pNP1=pNP1, pTheta_A=pTheta_A, pNP2=pNP2, pTheta_C=pTheta_C, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
    df <- data[[2]] # data for prpp_smart full analysis replicated and has calculated weights
    
    df$norm_w <- df$w/mean(df$w) # normalize weights to replicated sample size
    
    
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
    analysis_data_trad <- replicated_dat_trad %>% mutate(xalpha = 1,
                                                         xbeta = T1,
                                                         xtheta = T2,
                                                         xgamma = T1*T2,
                                                         norm_w = w/mean(w)) %>%
      dplyr::select(id, Y, xalpha, xbeta, xtheta, xgamma, norm_w, Treatment_Path)
    
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
    
    ### FULL PRPP-SMART ANALYSIS ###
    
    # create survey design object
    dclus1<-svydesign(id=~id, weights=~norm_w, data=df)
    # create replicated survey design object
    rclus1 <- survey::as.svrepdesign(design = dclus1, type = "JK1") # jacknife1 method 
    
    #Set the Data for Stan
    dfcs <- rclus1$variables
    weights_use <- rclus1$pweights
    data_stan <- with(dfcs, list(y=Y, xalpha00 = xalpha00, xalpha01 = xalpha01, xalpha10=xalpha10, xalpha11=xalpha11, xbeta0 = xbeta0, xbeta1 = xbeta1, xtheta0 = xtheta0, xtheta1 = xtheta1, xgamma=xgamma, N=nrow(dfcs), alpha_mu00=alpha_mu00, alpha_mu01=alpha_mu01, alpha_mu10=alpha_mu10, alpha_mu11=alpha_mu11, beta0_mu=beta0_mu, beta1_mu=beta1_mu, theta0_mu=theta0_mu, theta1_mu=theta1_mu, sigma=sigma, gamma_mu=gamma_mu, weights=weights_use))
    ctrl_stan<-list("chains"=n_MCMC_chain,"iter"=MCMC_SAMPLE,"warmup"=BURN.IN,"thin"=n_MCMC_thin)
    pars_stan <- c("alpha00","alpha01", "alpha10", "alpha11","beta0", "beta1", "theta0", "theta1", "gamma","DTR") # parameters to monitor 
    
    
    mod1 <- cs_sampling(svydes = rclus1, mod_stan = mod_dm,
                        data_stan = data_stan, par_stan = pars_stan, ctrl_stan = ctrl_stan, rep_design = TRUE) # if supplying rep survey object then set rep_design = TRUE if supply srvey design object set rep_design = FALSE
    
    
    # estimates from original sample
    DTR_hat_org[,i] = unname(colMeans(mod1$sampled_parms[,10:25]))
    rownames(DTR_hat_org) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
    
    parameter_hat_org[,i] = unname(colMeans(mod1$sampled_parms[,1:9]))
    rownames(parameter_hat_org) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
    
    DTR_var_org[,i] = unname(colVars(mod1$sampled_parms[,10:25]))
    rownames(DTR_var_org) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
    
    parameter_var_org[,i] = unname(colVars(mod1$sampled_parms[,1:9]))
    rownames(parameter_var_org) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
    
    # estimates from adjusted sample
    DTR_hat_adj[,i] = unname(colMeans(mod1$adjusted_parms[,10:25]))
    rownames(DTR_hat_adj) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
    
    parameter_hat_adj[,i] = unname(colMeans(mod1$adjusted_parms[,1:9]))
    rownames(parameter_hat_adj) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
    
    DTR_var_adj[,i] = unname(colVars(mod1$adjusted_parms[,10:25]))
    rownames(DTR_var_adj) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
    
    parameter_var_adj[,i] = unname(colVars(mod1$adjusted_parms[,1:9]))
    rownames(parameter_var_adj) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
    
    # check nominal coverage - original sample
    dtr1_quantiles <- c(quantile(mod1$sampled_parms[,10], .025), quantile(mod1$sampled_parms[,10], .975))
    dtr2_quantiles <- c(quantile(mod1$sampled_parms[,11], .025), quantile(mod1$sampled_parms[,11], .975))
    dtr3_quantiles <- c(quantile(mod1$sampled_parms[,12], .025), quantile(mod1$sampled_parms[,12], .975))
    dtr4_quantiles <- c(quantile(mod1$sampled_parms[,13], .025), quantile(mod1$sampled_parms[,13], .975))
    dtr5_quantiles <- c(quantile(mod1$sampled_parms[,14], .025), quantile(mod1$sampled_parms[,14], .975))
    dtr6_quantiles <- c(quantile(mod1$sampled_parms[,15], .025), quantile(mod1$sampled_parms[,15], .975))
    dtr7_quantiles <- c(quantile(mod1$sampled_parms[,16], .025), quantile(mod1$sampled_parms[,16], .975))
    dtr8_quantiles <- c(quantile(mod1$sampled_parms[,17], .025), quantile(mod1$sampled_parms[,17], .975))
    dtr9_quantiles <- c(quantile(mod1$sampled_parms[,18], .025), quantile(mod1$sampled_parms[,18], .975))
    dtr10_quantiles <- c(quantile(mod1$sampled_parms[,19], .025), quantile(mod1$sampled_parms[,19], .975))
    dtr11_quantiles <- c(quantile(mod1$sampled_parms[,20], .025), quantile(mod1$sampled_parms[,20], .975))
    dtr12_quantiles <- c(quantile(mod1$sampled_parms[,21], .025), quantile(mod1$sampled_parms[,21], .975))
    dtr13_quantiles <- c(quantile(mod1$sampled_parms[,22], .025), quantile(mod1$sampled_parms[,22], .975))
    dtr14_quantiles <- c(quantile(mod1$sampled_parms[,23], .025), quantile(mod1$sampled_parms[,23], .975))
    dtr15_quantiles <- c(quantile(mod1$sampled_parms[,24], .025), quantile(mod1$sampled_parms[,24], .975))
    dtr16_quantiles <- c(quantile(mod1$sampled_parms[,25], .025), quantile(mod1$sampled_parms[,25], .975))
    
    quantile_mat <- rbind(dtr1_quantiles, dtr2_quantiles, dtr3_quantiles, dtr4_quantiles, dtr5_quantiles, dtr6_quantiles, dtr7_quantiles, dtr8_quantiles, dtr9_quantiles, dtr10_quantiles, dtr11_quantiles, dtr12_quantiles, dtr13_quantiles, dtr14_quantiles, dtr15_quantiles, dtr16_quantiles)
    param_in_ci <- data.table::between(true_DTR_mat[,1], quantile_mat[, 1], quantile_mat[, 2])
    ci_hat_org[,i] <- param_in_ci
    rownames(ci_hat_org) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
    
    # check nominal coverage - adjusted sample
    dtr1_quantiles2 <- c(quantile(mod1$adjusted_parms[,10], .025), quantile(mod1$adjusted_parms[,10], .975))
    dtr2_quantiles2 <- c(quantile(mod1$adjusted_parms[,11], .025), quantile(mod1$adjusted_parms[,11], .975))
    dtr3_quantiles2 <- c(quantile(mod1$adjusted_parms[,12], .025), quantile(mod1$adjusted_parms[,12], .975))
    dtr4_quantiles2 <- c(quantile(mod1$adjusted_parms[,13], .025), quantile(mod1$adjusted_parms[,13], .975))
    dtr5_quantiles2 <- c(quantile(mod1$adjusted_parms[,14], .025), quantile(mod1$adjusted_parms[,14], .975))
    dtr6_quantiles2 <- c(quantile(mod1$adjusted_parms[,15], .025), quantile(mod1$adjusted_parms[,15], .975))
    dtr7_quantiles2 <- c(quantile(mod1$adjusted_parms[,16], .025), quantile(mod1$adjusted_parms[,16], .975))
    dtr8_quantiles2 <- c(quantile(mod1$adjusted_parms[,17], .025), quantile(mod1$adjusted_parms[,17], .975))
    dtr9_quantiles2 <- c(quantile(mod1$adjusted_parms[,18], .025), quantile(mod1$adjusted_parms[,18], .975))
    dtr10_quantiles2 <- c(quantile(mod1$adjusted_parms[,19], .025), quantile(mod1$adjusted_parms[,19], .975))
    dtr11_quantiles2 <- c(quantile(mod1$adjusted_parms[,20], .025), quantile(mod1$adjusted_parms[,20], .975))
    dtr12_quantiles2 <- c(quantile(mod1$adjusted_parms[,21], .025), quantile(mod1$adjusted_parms[,21], .975))
    dtr13_quantiles2 <- c(quantile(mod1$adjusted_parms[,22], .025), quantile(mod1$adjusted_parms[,22], .975))
    dtr14_quantiles2 <- c(quantile(mod1$adjusted_parms[,23], .025), quantile(mod1$adjusted_parms[,23], .975))
    dtr15_quantiles2 <- c(quantile(mod1$adjusted_parms[,24], .025), quantile(mod1$adjusted_parms[,24], .975))
    dtr16_quantiles2 <- c(quantile(mod1$adjusted_parms[,25], .025), quantile(mod1$adjusted_parms[,25], .975))
    
    quantile_mat2 <- rbind(dtr1_quantiles2, dtr2_quantiles2, dtr3_quantiles2, dtr4_quantiles2, dtr5_quantiles2, dtr6_quantiles2, dtr7_quantiles2, dtr8_quantiles2, dtr9_quantiles2, dtr10_quantiles2, dtr11_quantiles2, dtr12_quantiles2, dtr13_quantiles2, dtr14_quantiles2, dtr15_quantiles2, dtr16_quantiles2)
    param_in_ci2 <- data.table::between(true_DTR_mat[,1], quantile_mat2[, 1], quantile_mat2[, 2])
    ci_hat_adj[,i] <- param_in_ci2
    rownames(ci_hat_adj) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
    
    ### Traditional PRPP-SMART ANALYSIS ###
    
    # create survey design object
    dclust<-svydesign(id=~id, weights=~norm_w, data=analysis_data_trad)
    rclust <- survey::as.svrepdesign(design = dclust, type = "JK1") # jacknife1 method
    
    #Set the Data for Stan
    
    dfcst <- rclust$variables
    weights_use_t <- rclust$pweights
    data_stan_t <- with(dfcst, list(y=Y, xalpha = xalpha, xbeta = xbeta, xtheta = xtheta, xgamma=xgamma, N=nrow(dfcst), alpha_mu=alpha_mu, sigma=sigma, beta_mu = beta_mu, theta_mu = theta_mu, gamma_mu=gamma_mu, weights=weights_use_t))
    ctrl_stan<-list("chains"=n_MCMC_chain,"iter"=MCMC_SAMPLE,"warmup"=BURN.IN,"thin"=n_MCMC_thin)
    pars_stan_t <- c("alpha","beta", "theta", "gamma","DTR") # parameters to monitor 
    
    
    modt <- cs_sampling(svydes = rclust, mod_stan = mod_dm_traditional,
                        data_stan = data_stan_t, par_stan = pars_stan_t, ctrl_stan = ctrl_stan, rep_design = TRUE) # if supplying rep survey object then set rep_design = TRUE if supply srvey design object set rep_design = FALSE
    
    
    # estimates from original sample
    DTR_hat_org_t[,i] = unname(colMeans(modt$sampled_parms[,5:8]))
    rownames(DTR_hat_org_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
    
    parameter_hat_org_t[,i] = unname(colMeans(modt$sampled_parms[,1:4]))
    rownames(parameter_hat_org_t) <- c("alpha", "beta", "theta", "gamma")
    
    DTR_var_org_t[,i] = unname(colVars(modt$sampled_parms[,5:8]))
    rownames(DTR_var_org_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
    
    parameter_var_org_t[,i] = unname(colVars(modt$sampled_parms[,1:4]))
    rownames(parameter_var_org_t) <- c("alpha", "beta", "theta", "gamma")
    
    # estimates from adjusted sample
    DTR_hat_adj_t[,i] = unname(colMeans(modt$adjusted_parms[,5:8]))
    rownames(DTR_hat_adj_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
    
    parameter_hat_adj_t[,i] = unname(colMeans(modt$adjusted_parms[,1:4]))
    rownames(parameter_hat_adj_t) <- c("alpha", "beta", "theta", "gamma")
    
    DTR_var_adj_t[,i] = unname(colVars(modt$adjusted_parms[,5:8]))
    rownames(DTR_var_adj_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
    
    parameter_var_adj_t[,i] = unname(colVars(modt$adjusted_parms[,1:4]))
    rownames(parameter_var_adj_t) <- c("alpha", "beta", "theta", "gamma")
    
    # check nominal coverage - original sample
    dtr1_quantiles_t <- c(quantile(modt$sampled_parms[,5], .025), quantile(modt$sampled_parms[,5], .975))
    dtr2_quantiles_t <- c(quantile(modt$sampled_parms[,6], .025), quantile(modt$sampled_parms[,6], .975))
    dtr3_quantiles_t <- c(quantile(modt$sampled_parms[,7], .025), quantile(modt$sampled_parms[,7], .975))
    dtr4_quantiles_t <- c(quantile(modt$sampled_parms[,8], .025), quantile(modt$sampled_parms[,8], .975))
    
    quantile_mat_t <- rbind(dtr1_quantiles_t, dtr2_quantiles_t, dtr3_quantiles_t, dtr4_quantiles_t)
    param_in_ci_t <- data.table::between(true_DTR_mat_traditional[,1], quantile_mat_t[, 1], quantile_mat_t[, 2])
    ci_hat_org_t[,i] <- param_in_ci_t
    rownames(ci_hat_org_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
    
    # check nominal coverage - adjusted sample
    dtr1_quantiles2_t <- c(quantile(modt$adjusted_parms[,5], .025), quantile(modt$adjusted_parms[,5], .975))
    dtr2_quantiles2_t <- c(quantile(modt$adjusted_parms[,6], .025), quantile(modt$adjusted_parms[,6], .975))
    dtr3_quantiles2_t <- c(quantile(modt$adjusted_parms[,7], .025), quantile(modt$adjusted_parms[,7], .975))
    dtr4_quantiles2_t <- c(quantile(modt$adjusted_parms[,8], .025), quantile(modt$adjusted_parms[,8], .975))
    
    quantile_mat2_t <- rbind(dtr1_quantiles2_t, dtr2_quantiles2_t, dtr3_quantiles2_t, dtr4_quantiles2_t)
    param_in_ci2_t <- data.table::between(true_DTR_mat_traditional[,1], quantile_mat2_t[, 1], quantile_mat2_t[, 2])
    ci_hat_adj_t[,i] <- param_in_ci2_t
    rownames(ci_hat_adj_t) <- c("AAC00", "AAD00", "BBC00", "BBD00")
    
  }
  
  
  ######################### EVALUATION ############################
  
  ### Save Settings
  # set date and time for file saving 
  st<-format(Sys.time(), "%Y_%m_%d_%H_%M")
  
  # Define the folder path where you want your results saved
  folder_path <- "folder path here"
  
  # Example
  # folder_path <- paste0("/home/mariwank/PRPP_SMART_WRRM_Paper/SimResults/Bayesian/Sigma2_Results/Scenario", scenario, "/Sigma2", sigma2)
  
  # Create the folder if it doesn't exist
  if (!file.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Folder created at", folder_path, "\n")
  } else {
    cat("Folder already exists at", folder_path, "\n")
  }
  
  ### Parameter Full Model Results 
  ## Parameters
  Parameter = c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1","theta0", "theta1","gamma")
  param_hat_avg_org = c()
  param_hat_avg_adj = c()
  param_sd_hat_org = c()
  param_sd_hat_adj = c()
  param_avg_sd_org = c()
  param_avg_sd_adj = c()
  
  for(i in 1:length(Parameter)){
    param_hat_avg_org[i] <- round(mean(parameter_hat_org[i,], na.rm = TRUE),4)
    param_hat_avg_adj[i] <- round(mean(parameter_hat_adj[i,], na.rm = TRUE),4)
    param_sd_hat_org[i] <- round(sd(parameter_hat_org[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 posterior means)
    param_sd_hat_adj[i] <- round(sd(parameter_hat_adj[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 posterior means)
    param_avg_sd_org[i] <- round(sqrt(mean(parameter_var_org[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 posterior variances then square root to get sd) sqrt(mean(var(posterior draws))) equivalent to mean(sqrt(var(posterior draws)))
    param_avg_sd_adj[i] <- round(sqrt(mean(parameter_var_adj[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
  }
  
  param_results_tbl <- data.frame(Parameter=Parameter, 
                                  Param_Hat_Avg_org=param_hat_avg_org, 
                                  Param_Hat_Avg_adj=param_hat_avg_adj, 
                                  Sd_org = param_sd_hat_org, # sd of 500 posterior means
                                  Sd_adj = param_sd_hat_adj, # sd of 500 posterior means
                                  Avg_sd_org = param_avg_sd_org, # mean of 500 posterior variances then square root
                                  Avg_sd_adj = param_avg_sd_adj) # mean of 500 posterior variances then square root
  
  
  # Define the file name
  file_name <- paste0("Parameter_Results_Full_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(param_results_tbl, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ### DTR Full Model Results 
  DTR = c("AAC00","AAD00","BBC00","BBD00","AAC01","AAD01","BBC01","BBD01",
          "AAC10","AAD10","BBC10","BBD10","AAC11","AAD11","BBC11","BBD11")
  
  DTR_hat_avg_org = c()
  DTR_hat_avg_adj = c()
  DTR_sd_hat_org = c()
  DTR_sd_hat_adj = c()
  DTR_avg_sd_org = c()
  DTR_avg_sd_adj = c()
  DTR_avg_n = c()
  
  #DTR_avg_n = c()
  for(i in 1:16){
    DTR_hat_avg_org[i] <- round(mean(DTR_hat_org[i,], na.rm = TRUE),4)
    DTR_hat_avg_adj[i] <- round(mean(DTR_hat_adj[i,], na.rm = TRUE),4)
    DTR_sd_hat_org[i] <- round(sd(DTR_hat_org[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 posterior means)
    DTR_sd_hat_adj[i] <- round(sd(DTR_hat_adj[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 posterior means)
    DTR_avg_sd_org[i] <- round(sqrt(mean(DTR_var_org[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
    DTR_avg_sd_adj[i] <- round(sqrt(mean(DTR_var_adj[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
    
    DTR_avg_n[i] <- round(mean(n.DTR[i,], na.rm = TRUE),1)
  }
  
  # calculate bias
  DTR_bias_org = round(DTR_hat_avg_org - expected_pref,4)
  DTR_bias_adj = round(DTR_hat_avg_adj - expected_pref,4)
  
  # calculate rMSE
  rMSE_DTR_org <- sqrt(DTR_sd_hat_org^2 + DTR_bias_org^2)
  rMSE_DTR_adj <- sqrt(DTR_sd_hat_adj^2 + DTR_bias_adj^2)
  
  DTR_results_tbl <- data.frame(DTR=DTR,
                                True_DTR=expected_pref,
                                DTR_Hat_Avg_org=DTR_hat_avg_org, 
                                DTR_Hat_Avg_adj=DTR_hat_avg_adj,
                                Bias_org=DTR_bias_org, 
                                Bias_adj=DTR_bias_adj,
                                Sd_org = DTR_sd_hat_org, 
                                Sd_adj = DTR_sd_hat_adj, 
                                Avg_sd_org = DTR_avg_sd_org, 
                                Avg_sd_adj = DTR_avg_sd_adj, 
                                rMSE_org = rMSE_DTR_org,
                                rMSE_adj = rMSE_DTR_adj) 
  
  
  # Define the file name
  file_name <- paste0("DTR_Results_Full_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(DTR_results_tbl, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  
  ## nominal coverage 
  ci_coverage_org = c()
  ci_coverage_adj = c()
  
  for(i in 1:nrow(ci_hat_org)){
    ci_coverage_org[i] <- round(mean(ci_hat_org[i,], na.rm = TRUE),2)
    ci_coverage_adj[i] <- round(mean(ci_hat_adj[i,], na.rm = TRUE),2)
    
  }
  
  ci_results_tbl <- data.frame(DTR, 
                               CI_Coverage_org=ci_coverage_org, 
                               CI_Coverage_adj=ci_coverage_adj)
  
  
  
  # Define the file name
  file_name <- paste0("DTR_CI_Results_Full_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(ci_results_tbl, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  
  ### Parameter Traditional Results 
  ## Parameters
  Parameter_t = c("alpha", "beta", "theta", "gamma")
  param_hat_avg_org_t = c()
  param_hat_avg_adj_t = c()
  param_sd_hat_org_t = c()
  param_sd_hat_adj_t = c()
  param_avg_sd_org_t = c()
  param_avg_sd_adj_t = c()
  
  for(i in 1:length(Parameter_t)){
    param_hat_avg_org_t[i] <- round(mean(parameter_hat_org_t[i,], na.rm = TRUE),4)
    param_hat_avg_adj_t[i] <- round(mean(parameter_hat_adj_t[i,], na.rm = TRUE),4)
    param_sd_hat_org_t[i] <- round(sd(parameter_hat_org_t[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 posterior means)
    param_sd_hat_adj_t[i] <- round(sd(parameter_hat_adj_t[i,], na.rm = TRUE),4) # sd of mean parameter estimates over 500 sims (sd of 500 posterior means)
    param_avg_sd_org_t[i] <- round(sqrt(mean(parameter_var_org_t[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
    param_avg_sd_adj_t[i] <- round(sqrt(mean(parameter_var_adj_t[i,], na.rm = TRUE)),4) # mean of variance of the parameter estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
  }
  
  param_results_tbl_t <- data.frame(Parameter=Parameter_t, 
                                    Param_Hat_Avg_org=param_hat_avg_org_t, 
                                    Param_Hat_Avg_adj=param_hat_avg_adj_t, 
                                    Sd_org = param_sd_hat_org_t, # sd of 500 posterior means
                                    Sd_adj = param_sd_hat_adj_t, # sd of 500 posterior means
                                    Avg_sd_org = param_avg_sd_org_t, # mean of 500 posterior variances then square root
                                    Avg_sd_adj = param_avg_sd_adj_t) # mean of 500 posterior variances then square root
  
  # Define the file name
  file_name <- paste0("Parameter_Results_Traditional_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(param_results_tbl_t, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ### DTR TRADITIONAL Results
  
  DTR_t = c("AAC00","AAD00","BBC00","BBD00") # DTRs from traditional randomized randomized smart
  DTR_hat_avg_org_t = c()
  DTR_hat_avg_adj_t = c()
  DTR_sd_hat_org_t = c()
  DTR_sd_hat_adj_t = c()
  DTR_avg_sd_org_t = c()
  DTR_avg_sd_adj_t = c()
  
  #DTR_avg_n = c()
  for(i in 1:4){
    DTR_hat_avg_org_t[i] <- round(mean(DTR_hat_org_t[i,], na.rm = TRUE),4)
    DTR_hat_avg_adj_t[i] <- round(mean(DTR_hat_adj_t[i,], na.rm = TRUE),4)
    DTR_sd_hat_org_t[i] <- round(sd(DTR_hat_org_t[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 posterior means)
    DTR_sd_hat_adj_t[i] <- round(sd(DTR_hat_adj_t[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 posterior means)
    DTR_avg_sd_org_t[i] <- round(sqrt(mean(DTR_var_org_t[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
    DTR_avg_sd_adj_t[i] <- round(sqrt(mean(DTR_var_adj_t[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
    
  }
  
  # calculate bias
  DTR_bias_org_t = round(DTR_hat_avg_org_t - expected_pref[1:4],4)
  DTR_bias_adj_t = round(DTR_hat_avg_adj_t - expected_pref[1:4],4)
  
  # calculate rMSE
  rMSE_DTR_org_t <- sqrt(DTR_sd_hat_org_t^2 + DTR_bias_org_t^2)
  rMSE_DTR_adj_t <- sqrt(DTR_sd_hat_adj_t^2 + DTR_bias_adj_t^2)
  
  DTR_results_tbl_t <- data.frame(DTR=DTR_t,
                                  True_DTR=expected_pref[1:4],
                                  DTR_Hat_Avg_org=DTR_hat_avg_org_t, 
                                  Bias_org=DTR_bias_org_t, 
                                  DTR_Hat_Avg_adj=DTR_hat_avg_adj_t, 
                                  Bias_adj=DTR_bias_adj_t,
                                  Sd_org = DTR_sd_hat_org_t, 
                                  Sd_adj = DTR_sd_hat_adj_t, 
                                  Avg_sd_org = DTR_avg_sd_org_t, 
                                  Avg_sd_adj = DTR_avg_sd_adj_t, 
                                  rMSE_org = rMSE_DTR_org_t, 
                                  rMSE_adj = rMSE_DTR_adj_t) 
  
  # Define the file name
  file_name <- paste0("DTR_Results_Traditional_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(DTR_results_tbl_t, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  ## nominal coverage 
  ci_coverage_org_t = c()
  ci_coverage_adj_t = c()
  
  for(i in 1:nrow(ci_hat_org_t)){
    ci_coverage_org_t[i] <- round(mean(ci_hat_org_t[i,], na.rm = TRUE),2)
    ci_coverage_adj_t[i] <- round(mean(ci_hat_adj_t[i,], na.rm = TRUE),2)
    
  }
  ci_results_tbl_t <- data.frame(DTR_t, 
                                 CI_Coverage_org=ci_coverage_org_t, 
                                 CI_Coverage_adj=ci_coverage_adj_t)
  
  
  # Define the file name
  file_name <- paste0("DTR_CI_Results_Traditional_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(ci_results_tbl_t, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  
  ### EVALUATION - COMPARE FULL VS TRADITIONAL RESULTS ###
  
  ## DTR compare
  # based off of adjusted samples
  DTR_results_tbl_c <- data.frame(DTR=DTR_t,
                                  True_DTR=expected_pref[1:4],
                                  DTR_Hat_Avg_full = DTR_hat_avg_adj[1:4], # full prpp-smart analysis (all subjects)
                                  DTR_Hat_Avg_trad = DTR_hat_avg_adj_t, # traditional analysis (randomized randomized subjects only)
                                  Bias_full = DTR_bias_adj[1:4], # full prpp-smart analysis (all subjects)
                                  Bias_trad = DTR_bias_adj_t, # traditional analysis (randomized randomized subjects only)
                                  Sd_full = DTR_sd_hat_adj[1:4], 
                                  Sd_trad = DTR_sd_hat_adj_t, 
                                  Avg_sd_full = DTR_avg_sd_adj[1:4], 
                                  Avg_sd_trad = DTR_avg_sd_adj_t, 
                                  rMSE_full = rMSE_DTR_adj[1:4], 
                                  rMSE_trad = rMSE_DTR_adj_t) 
  
  
  # Define the file name
  file_name <- paste0("DTR_CompareResults_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(DTR_results_tbl_c, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  
  # CI DTR compare
  DTR_ci_tbl_c <- data.frame(DTR=DTR_t,
                             Coverage_full = ci_coverage_adj[1:4],
                             Coverage_trad = ci_coverage_adj_t)
  
  # Define the file name
  file_name <- paste0("DTR_CI_CompareResults_Sigma2", sigma2, "_", st, ".csv")
  
  # Create the file path
  file_path <- file.path(folder_path, file_name)
  
  # Write the data to a CSV file
  write.csv(DTR_ci_tbl_c, file = file_path, row.names = FALSE)
  cat("CSV file saved to", file_path, "\n")
  
  
}







