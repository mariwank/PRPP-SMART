###########################################################################################################################

# Data Simulation code for PRPP-SMART

# Author: Mari Wank
###########################################################################################################################
# Description: 
#    Code simulates a dataset from a two-stage PRPP-SMART design with a binary end of stage outcomes. 

###########################################################################################################################

#Function: gendata


#Purpose: This function generates one simulation dataset for a two stage PRPP-SMART trial with a binary outcome


# Required Parameters: 
#         N: Number of individuals in the generated dataset
#         pNP1: Desired proportion of individuals expressing No Preference in stage 1
#         pTheta_A: Desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
#         pNP2: Desired proportion of patients expressing No Preference in stage 2 (among non-responders)
#         pTheta_C: Desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)
#         pi_A:  Probability of responding in stage 1 given randomized treatment A
#         pi_B:  Probability of responding in stage 1 given randomized treatment B
#         pi_A1: Probability of responding in stage 1 given preferred treatment A
#         pi_B1: Probability of responding in stage 1 given preferred treatment B
#         pi_AC: Probability of Y=1 given treatment path A0C0
#         pi_AD: Probability of Y=1 given treatment path A0D0
#         pi_BC: Probability of Y=1 given treatment path B0C0
#         pi_BD: Probability of Y=1 given treatment path B0D0
#         pA0A: Probability of Y=1 given treatment path A0A
#         pB0B: Probability of Y=1 given treatment path B0B
#         pA1A: Probability of Y=1 given treatment path A1A
#         pB1B: Probability of Y=1 given treatment path B1B
#         pA1C1: Probability of Y=1 given treatment path A1C1
#         pA1C0: Probability of Y=1 given treatment path A1C0
#         pA1D1: Probability of Y=1 given treatment path A1D1
#         pA1D0: Probability of Y=1 given treatment path A1D0
#         pA0C1: Probability of Y=1 given treatment path A0C1
#         pA0D1: Probability of Y=1 given treatment path A0D1
#         pB1C1: Probability of Y=1 given treatment path B1C1
#         pB1C0: Probability of Y=1 given treatment path B1C0
#         pB1D1: Probability of Y=1 given treatment path B1D1
#         pB1D0: Probability of Y=1 given treatment path B1D0
#         pB0C1: Probability of Y=1 given treatment path B0C1
#         pB0D1: Probability of Y=1 given treatment path B0D1

#Output:   
#       (1) A dataset of N subjects before replication 
#             The dataset has the following variables:
#               ID: Numeric subject ID variable 
#               X1: A continuous baseline variable generated from N(0,1)
#               X2: A binary baseline variable generated from Bern(0.5)
#               S1_Preference: Categorical variable; Preference for first stage treatment; takes values A, B, NP; A-Prefer A, B-Prefer B, NP-No preference
#               P1: Binary Variable: Indicator for whether an individual had a preference in stage 1. 1-had a preference, 0-had no preference
#               T1: Binary variable; Individuals treatment assigned at Stage 1; takes Values A(1) or B(-1). 
#               R: Binary variable; Individual's stage 1 response. 1 denotes response, 0 denotes no response
#               S2_Preference: Categorical variable; Preference for second stage treatment; takes values C, D, NP, 999. C-Prefer C, D-Prefer D, NP-No preference, 999-responders (stage 1 responders are not asked stage 2 preference)
#               P2: Binary Variable: Indicator for whether an individual had a preference in stage 2. 1-had a preference, 0-had no preference, 999-responder in stage 1
#               T2: Categorical variable; Individuals treatment assigned at Stage 2; takes values C(1), D(-1), 999-responders (stage 1 resopnder continue on initial treatment) 
#               Y: Binary variable; Individual's binary endpoint. 1 denotes response, 0 denotes no response
#               Treatment_Path: Treatment path the participant follows

#       (2) A weighted and replicated dataset to be used in a WRRM model
#             The dataset has the following variables:
#               ID: Numeric subject ID variable 
#               xalpha00: intercept alpha_00 in WRRM
#               xalpha01: intercept alpha_01 in WRRM
#               xalpha10: intercept alpha_10 in WRRM
#               xalpha11: intercept alpha_11 in WRRM
#               xbeta0: stage 1 randomized treatment beta0 in WRRM
#               xbeta1: stage 1 preference treatment beta1 in WRRM
#               xtheta0: stage 2 randomized treatment theta0 in WRRM
#               xtheta1: stage 2 preference treatment theta1 in WRRM
#               xgamma: stage 1 and 2 treatment interaction gamma in WRRM
#               Y: Binary variable; Individual's binary endpoint. 1 denotes response, 0 denotes no response
#               w: PRPP-SMART weight
#               Treatment_Path: Treatment path the participant follows
#          
#       (3) A dataset of sample sizes for each embedded DTR in the PRPP-SMART. That is, the number of participants (responders and non-responders) that constitute each embedded DTR in the simulated data


###########################################################################################################################

generate_data <- function(N, pNP1, pTheta_A, pNP2, pTheta_C, pi_A, pi_B, pi_A1, pi_B1, pi_AC, pi_AD, pi_BC, pi_BD, pA0A, pB0B, pA1A, pB1B, pA1C1, pA1C0, pA1D1, pA1D0, pA0C1, pA0D1, pB1C1, pB1C0, pB1D1, pB1D0, pB0C1, pB0D1){
  
  # load libraries 
  library(Rlab)
  library(tidyverse)
  
  expit <- function(y) exp(y)/(1+exp(y))
  
  
  #Function: preference 
  #Purpose: These functions assign the preference of each individual based on the true propensity for exhibiting 
  #         a preference for treatment A, B, or having no preference in stage 1 and treatment C, D, or having no
  #         no preference in stage 2. Note, pref_probx is a Nx3 matrix with the
  #         true probabilities of preference for each individual
  
  
  preference_stage1 <- function(i){
    sample(c("A", "B", "NP"), size=1, replace = TRUE, prob=pref_prob1[i,])
  }
  
  preference_stage2 <- function(i){
    sample(c("C", "D", "NP"), size=1, replace = TRUE, prob=pref_prob2[i,])
  }

  
  # generate baseline covariates 
  X1<-rnorm(N,mean=0,sd=1)  # generate X1
  X2<-rbinom(N,1,0.5)  # generate X2: Home vs clinic/doctor’s office for treatment

  
  ### Preference Model: Stage 1 Propensity Scores for Preference ###
  #    Here we generate preference for first stage treatment for each patient. The true propensity for exhibiting a 
  #    preference for treatment A, B, or having no preference for each subject in the simulation population 
  #    is modeled using a logit model approach where first-stage propensity scores for preference are conditional on the 
  #    first-stage baseline covariates X1 and X2. In our preference model, 
  #    we search for an intercept via the uniroot function in r which allows us to control the proportion of subjects exhibiting 
  #    no preference/preference for A/B
  
  
  ## No Preference ##
  
  # Coefficients on X1 and X2 
  a1 <- 0.2
  a2 <- 0.133 
  
  # Find intercept to get probability of no preference close to target, on average
  search0_NP1 <- function(a0)
  {
    alpha <- rbind(a0, a1, a2)
    lnp <- cbind(1, X1, X2) %*% alpha # log scale
    pNP <- expit(lnp) # probability scale
    mean(pNP-pNP1) # want this to be close to 0
  }
  a0_star <- uniroot(search0_NP1, c(-4, 4))$root # intercept value that will produce a marginal no preference allocation close to target
  alpha <- c(a0_star, a1, a2)
  prob_NP <- expit(cbind(1, X1, X2) %*% alpha) # calculate no preference probability for each patient
  
  
  ## Model Theta: Prefer A among those with a preference ##
  pA_marginal <- pTheta_A*(1-pNP1) # marginal A probability in simulated data
  
  # Coefficients on X1 and X2 
  b1 <- 0.05
  b2 <- 0.1
  
  # Find intercept to get theta close to target, on average
  search0_ThetaA <- function(b0)
  {
    beta <- rbind(b0, b1, b2)
    lptheta<- cbind(1, X1, X2) %*% beta # log scale
    ptheta <- expit(lptheta) # probability scale
    mean(ptheta-pTheta_A) # want this to be close to 0
  }
  
  b0_star <- uniroot(search0_ThetaA, c(-4,4))$root 
  beta <- c(b0_star, b1, b2)
  prob_A <- expit(cbind(1, X1, X2) %*% beta) * (1-prob_NP)
  prob_B <- (1-expit(cbind(1, X1, X2) %*% beta)) * (1-prob_NP) 
  
  ## generate preference of first stage treatment of each individual ##
  
  # put probabilities into a matrix
  pref_prob1 <- cbind(prob_A,prob_B,prob_NP) 
  colnames(pref_prob1) <- c("Prob Prefer A", "Prob Prefer B", "Prob No Preference") # each row is a subject's probability to prefer A, prefer B, or have no preference
  # apply(pref_prob1,1,sum) # check
  
  # sample based on probabilities 
  S1_Preference <- sapply(1:N, preference_stage1) # stage 1 preference 
  
  ### CHECKS ###
  # sum(prop.table(table(S1_Preference))) # check to make sure prob_A, prob_B, prob_NP sum to 1
  # prop.table(table(S1_Preference)) # check if observed proportions close to targets
  # length(which(S1_Preference == "A")) / length(which(S1_Preference == "A" | S1_Preference == "B")) # should be close to pTheta_A
  
  
  ### Stage 1 treatment assigned ###
  #    For each subject in the simulated population we generate their stage 1 assigned treatment
  #    based on the true propensities for exhibiting a preference for treatment A, B, or having no preference
  #    Prefer A: get A
  #    Prefer B: get B
  #    No preference: equal probability of getting treatment A/B
  
  
  # generate actual stage 1 treatment that each individual received
  # 1-treatment A, 0-treatment B
  T1 <- 
    (S1_Preference == "A") * 1 + # prefer A get A
    (S1_Preference == "B") * 0 + # prefer B get B
    (S1_Preference == "NP") * sample(c(1, 0), N, replace = T, prob=c(0.5,0.5)) # no preference randomly assign A/B
  
  
  ### Generate Preference Indicator ###
  #    Here we create an indicator, P1, for whether an individual had a preference in stage 1 or not. Specifically, 1
  #    if an individual had a preference in stage 1, 0 if an individual had no preference in stage 1
  
  P1 <- ifelse(S1_Preference == "A" | S1_Preference == "B", 1, 0)
  
  ### Generate Stage 1 Response variable ###
  #    For each subject we generate a binary response where 1 indicates the subjects responds
  #    to the assigned Stage 1 treatment and 0 indicates the subject does not response to the assigned 
  #    stage 1 treatment. 
  
  # Generate stage 1 response variable 
  R <- rep(NA,N) # 0-no response to treatment, 1-responds to treatment
  
  # response paths
  R[P1 == 1 & T1==1] <- rbinom(sum(P1 == 1 & T1==1),1, pi_A1) # got preferred treatment A in Stage 1
  R[P1 == 0 & T1==1] <- rbinom(sum(P1 == 0 & T1==1),1, pi_A) # got randomized treatment A in Stage 1
  R[P1 == 1 & T1==0] <- rbinom(sum(P1 == 1 & T1==0),1, pi_B1) # got preferred treatment B in Stage 1
  R[P1 == 0 & T1==0] <- rbinom(sum(P1 == 0 & T1==0),1, pi_B) # got randomized treatment B in Stage 1
  
  ### Prefernce Model: Stage 2 Propensity Scores for Preference ###
  #    Here we generate preference of second stage treatment only for first-stage non-responders. Patients who respond
  #    continue on the stage 1 treatment. The true propensity for exhibiting a 
  #    preference for treatment C, D, or having no preference for non-responders in the simulation population 
  #    is modeled using a logit model approach where second-stage propensity scores for preference are conditional on whether
  #    the patient had a preference in stage 1. In our preference model, 
  #    we search for an intercept via the uniroot function in r which allows us to control the proportion of non-responder subjects 
  #    exhibiting no preference/preference for C/D
  
  # find non-responders 
  nr_index <- which(R == 0)
  
  ## No Preference ##
  
  # Coefficients on P1 
  c1 <- 0.3
  
  # Find intercept to get probability of no preference close to target, on average
  search0_NP2 <- function(c0){
    alpha <- rbind(c0, c1)
    lnp <- cbind(1, P1[nr_index]) %*% alpha # log scale
    pNP <- expit(lnp) # probability scale
    mean(pNP-pNP2) # want this to be close to 0
  }
  c0_star <- uniroot(search0_NP2, c(-4, 4))$root # intercept value that will produce a marginal no preference allocation close to target
  gamma <- c(c0_star, c1)
  prob_NP2 <- expit(cbind(1, P1[nr_index]) %*% gamma) # calculate no preference probability for each patient
  
  
  ## Model Theta: Prefer C among those with a preference ##
  pC_marginal <- pTheta_C*(1-pNP2) # marginal C probability in simulated data
  
  # Coefficients on S1_P
  d1 <- 0.1
  
  # Find intercept to get theta close to target, on average
  search0_ThetaC <- function(d0)
  {
    beta <- rbind(d0, d1)
    lptheta<- cbind(1, P1[nr_index]) %*% beta # log scale
    ptheta <- expit(lptheta) # probability scale
    mean(ptheta-pTheta_C) # want this to be close to 0
  }
  
  d0_star <- uniroot(search0_ThetaC, c(-4,4))$root 
  phi <- c(d0_star, d1)
  prob_C <- expit(cbind(1, P1[nr_index]) %*% phi) * (1-prob_NP2)
  prob_D <- (1-expit(cbind(1, P1[nr_index]) %*% phi)) * (1-prob_NP2) 
  
  ## generate preference of second stage treatment for non-responders ##
  
  # put probabilities into a matrix
  pref_prob2 <- cbind(prob_C,prob_D,prob_NP2) 
  colnames(pref_prob2) <- c("Prob Prefer C", "Prob Prefer D", "Prob No Preference") # each row is a subject's probability to prefer C, prefer D, or have no preference
  # apply(pref_prob2,1,sum) # check
  
  # sample based on probabilities 
  n2 <- length(nr_index)
  S2_Preference <- sapply(1:n2, preference_stage2) # stage 2 preference only generated for stage 1 non-responders
  
  ### CHECKS ###
  #sum(prop.table(table(S2_Preference))) # check to make sure prob_C, prob_D, prob_NP2 sum to 1
  # prop.table(table(S2_Preference)) # check if observed proportions close to targets
  # length(which(S2_Preference == "C)) / length(which(S2_Preference == "C" | S2_Preference == "D")) # should be close to pTheta_C
  
  
  # create final stage 2 preference variable 
  S2_Preference_final <- rep(999, N) # 999 for stage 1 responders since they continue on stage 1 treatment 
  
  S2_Preference_final[nr_index] <- S2_Preference
  
  ### Stage 2 treatment assigned ###
  #    For each subject in the simulated population we generate their stage 2 assigned treatment.
  #    For non-responders this is based off the true propensities for exhibiting a preference for treatment C, D, or having no preference
  #    For responders they continue on their stage 1 assigned treatment.
  #    Non-responders:
  #      Prefer C: get C
  #      Prefer D: get D
  #      No preference: equal probability of getting treatment C/D
  
  
  # generate actual stage 2 treatment that each individual received
  # 1-treatment C, 0-treatment D
  T2 <- c() 
  for(i in 1:N){
    
    if (R[i]==1){ # responds to Stage 1 treatment
      T2[i] <- T1[i] # continues on initial treatment received 
      
    }
    if (R[i] == 0 & S2_Preference_final[i] == "NP"){# no response to Stage 1 treatment & has no preference in stage 2
      T2[i] <- rbern(1,0.5) # equal chance for treatment C/D
    }
    
    if (R[i] == 0 &  S2_Preference_final[i] == "C"){# no response to Stage 1 treatment & prefer C
      T2[i] <- 1 # get preferred treatment 
    }
    
    if (R[i] == 0 &  S2_Preference_final[i] == "D"){# no response to Stage 1 treatment & prefer D
      T2[i] <- 0 # get preferred treatment 
    }
  }
  
  
  ### Generate Stage 2 Preference Indicator ###
  #    Here we create an indicator, S2_P, for whether an individual had a preference in stage 2 or not. Specifically, 1
  #    if an individual had a preference in stage 2, 0 if an individual had no preference in stage 2
  
  P2 <- ifelse(S2_Preference_final == "C" | S2_Preference_final == "D", 1, 0)
  
  P2[which(R==1)] <- 999
  
  # Relabel T1 and T2 with actual treatment 
  T1[T1==1] <- "A"
  T1[T1==0] <- "B"
  
  T2[which(R ==1)] <- 999 # assign stage 1 responders 999 
  T2[T2==1] <- "C"
  T2[T2==0] <- "D"
  
  
  ### Generate Trial Outcome Variable ###
  #    For each subject we generate the final binary trial outcome given their study path where 1 indicates the subjects responds
  #    0 indicates the subject does not respond.
  

  Y<-rep(NA,N) # 0-no response to treatment, 1-responds to treatment
  
  # response paths
  Y[P1 == 1 & T1=="A" & R==1]<-rbinom(sum(P1 == 1 & T1=="A" & R==1),1,pA1A)
  Y[P1 == 0 & T1=="A" & R==1]<-rbinom(sum(P1 == 0 & T1=="A" & R==1),1,pA0A)
  Y[P1 == 1 & T1=="B" & R==1]<-rbinom(sum(P1 == 1 & T1=="B" & R==1),1,pB1B)
  Y[P1 == 0 & T1=="B" & R==1]<-rbinom(sum(P1 == 0 & T1=="B" & R==1),1,pB0B)
  
  # paths that begin with prefer A in stage 1
  Y[P1 == 1 & T1=="A" & R==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="A" & R==0 & P2 == 1 & T2 == "C"),1,pA1C1)
  Y[P1 == 1 & T1=="A" & R==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="A" & R==0 & P2 == 0 & T2 == "C"),1,pA1C0)
  Y[P1 == 1 & T1=="A" & R==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="A" & R==0 & P2 == 1 & T2 == "D"),1,pA1D1)
  Y[P1 == 1 & T1=="A" & R==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="A" & R==0 & P2 == 0 & T2 == "D"),1,pA1D0)
  
  # paths that begin with randomize A in stage 1
  Y[P1 == 0 & T1=="A" & R==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="A" & R==0 & P2 == 1 & T2 == "C"),1,pA0C1)
  Y[P1 == 0 & T1=="A" & R==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="A" & R==0 & P2 == 0 & T2 == "C"),1,pi_AC)
  Y[P1 == 0 & T1=="A" & R==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="A" & R==0 & P2 == 1 & T2 == "D"),1,pA0D1)
  Y[P1 == 0 & T1=="A" & R==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="A" & R==0 & P2 == 0 & T2 == "D"),1,pi_AD)
  
  # paths that begin with prefer B in stage 1
  Y[P1 == 1 & T1=="B" & R==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="B" & R==0 & P2 == 1 & T2 == "C"),1,pB1C1)
  Y[P1 == 1 & T1=="B" & R==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 1 & T1=="B" & R==0 & P2 == 0 & T2 == "C"),1,pB1C0)
  Y[P1 == 1 & T1=="B" & R==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="B" & R==0 & P2 == 1 & T2 == "D"),1,pB1D1)
  Y[P1 == 1 & T1=="B" & R==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 1 & T1=="B" & R==0 & P2 == 0 & T2 == "D"),1,pB1D0)
  
  
  # paths that begin with randomize B in stage 1
  Y[P1 == 0 & T1=="B" & R==0 & P2 == 1 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="B" & R==0 & P2 == 1 & T2 == "C"),1,pB0C1)
  Y[P1 == 0 & T1=="B" & R==0 & P2 == 0 & T2 == "C"]<-rbinom(sum(P1 == 0 & T1=="B" & R==0 & P2 == 0 & T2 == "C"),1,pi_BC)
  Y[P1 == 0 & T1=="B" & R==0 & P2 == 1 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="B" & R==0 & P2 == 1 & T2 == "D"),1,pB0D1)
  Y[P1 == 0 & T1=="B" & R==0 & P2 == 0 & T2 == "D"]<-rbinom(sum(P1 == 0 & T1=="B" & R==0 & P2 == 0 & T2 == "D"),1,pi_BD)
  
  
  ### Create final dataframe ### 
  
  
  ### Create trial path variable ###
  treatment_path <- rep(NA, N)
  
  # responder paths
  treatment_path[which(P1==1 & T1=="A" & R==1)] <- "A1A"
  treatment_path[which(P1==0 & T1=="A" & R==1)] <- "A0A"
  treatment_path[which(P1==1 & T1=="B" & R==1)] <- "B1B"
  treatment_path[which(P1==0 & T1=="B" & R==1)] <- "B0B"
  
  # non-responder paths
  treatment_path[which(P1==1 & T1=="A" & R==0 & P2==1 & T2=="C")] <- "A1C1"
  treatment_path[which(P1==1 & T1=="A" & R==0 & P2==0 & T2=="C")] <- "A1C0"
  treatment_path[which(P1==1 & T1=="A" & R==0 & P2==1 & T2=="D")] <- "A1D1"
  treatment_path[which(P1==1 & T1=="A" & R==0 & P2==0 & T2=="D")] <- "A1D0"
  
  treatment_path[which(P1==1 & T1=="B" & R==0 & P2==1 & T2=="C")] <- "B1C1"
  treatment_path[which(P1==1 & T1=="B" & R==0 & P2==0 & T2=="C")] <- "B1C0"
  treatment_path[which(P1==1 & T1=="B" & R==0 & P2==1 & T2=="D")] <- "B1D1"
  treatment_path[which(P1==1 & T1=="B" & R==0 & P2==0 & T2=="D")] <- "B1D0"
  
  treatment_path[which(P1==0 & T1=="A" & R==0 & P2==1 & T2=="C")] <- "A0C1"
  treatment_path[which(P1==0 & T1=="A" & R==0 & P2==0 & T2=="C")] <- "A0C0"
  treatment_path[which(P1==0 & T1=="A" & R==0 & P2==1 & T2=="D")] <- "A0D1"
  treatment_path[which(P1==0 & T1=="A" & R==0 & P2==0 & T2=="D")] <- "A0D0"
  
  treatment_path[which(P1==0 & T1=="B" & R==0 & P2==1 & T2=="C")] <- "B0C1"
  treatment_path[which(P1==0 & T1=="B" & R==0 & P2==0 & T2=="C")] <- "B0C0"
  treatment_path[which(P1==0 & T1=="B" & R==0 & P2==1 & T2=="D")] <- "B0D1"
  treatment_path[which(P1==0 & T1=="B" & R==0 & P2==0 & T2=="D")] <- "B0D0"
  
  data_output <- 
    tibble(
      id = 1:N,
      X1 = X1,
      X2 = X2,
      S1_Preference = S1_Preference,
      P1 = P1,
      T1 = T1,
      R = R,
      S2_Preference = S2_Preference_final,
      P2 = P2,
      T2 = T2,
      Y = Y,
      Treatment_Path = treatment_path
    )
  
  # create list of sample sizes for each DTR
  DTR_n = c()
  DTR_n[1] = sum(data_output$Treatment_Path == "A0A") + 
    sum(data_output$Treatment_Path == "A0C0") #AC00
  DTR_n[2] = sum(data_output$Treatment_Path == "A0A") +
    sum(data_output$Treatment_Path == "A0D0") #AD00
  DTR_n[3] = sum(data_output$Treatment_Path == "B0B") + 
    sum(data_output$Treatment_Path == "B0C0") #BC00
  DTR_n[4] = sum(data_output$Treatment_Path == "B0B") +
    sum(data_output$Treatment_Path == "B0D0") #BD00
  DTR_n[5] = sum(data_output$Treatment_Path == "A0A") + 
    sum(data_output$Treatment_Path == "A0C1") #AC01
  DTR_n[6] = sum(data_output$Treatment_Path == "A0A") +
    sum(data_output$Treatment_Path == "A0D1") #AD01
  DTR_n[7] = sum(data_output$Treatment_Path == "B0B") + 
    sum(data_output$Treatment_Path == "B0C1") #BC01
  DTR_n[8] = sum(data_output$Treatment_Path == "B0B") +
    sum(data_output$Treatment_Path == "B0D1") #BD01
  DTR_n[9] = sum(data_output$Treatment_Path == "A1A") + 
    sum(data_output$Treatment_Path == "A1C0") #AC10
  DTR_n[10] = sum(data_output$Treatment_Path == "A1A") +
    sum(data_output$Treatment_Path == "A1D0") #AD10
  DTR_n[11] = sum(data_output$Treatment_Path == "B1B") + 
    sum(data_output$Treatment_Path == "B1C0") #BC10
  DTR_n[12] = sum(data_output$Treatment_Path == "B1B") +
    sum(data_output$Treatment_Path == "B1D0") #BD10
  DTR_n[13] = sum(data_output$Treatment_Path == "A1A") + 
    sum(data_output$Treatment_Path == "A1C1") #AC11
  DTR_n[14] = sum(data_output$Treatment_Path == "A1A") +
    sum(data_output$Treatment_Path == "A1D1") #AD11
  DTR_n[15] = sum(data_output$Treatment_Path == "B1B") + 
    sum(data_output$Treatment_Path == "B1C1") #BC11
  DTR_n[16] = sum(data_output$Treatment_Path == "B1B") +
    sum(data_output$Treatment_Path == "B1D1") #BD11
  
  
  # relabel data 
  data_output_relabel <- data_output %>% mutate(T1 = replace(T1, T1 == "A", 1),
                                                T1 = replace(T1, T1 == "B", -1),
                                                T1 = as.numeric(T1),
                                                T2 = replace(T2, T2 == "C", 1),
                                                T2 = replace(T2, T2 == "D", -1),
                                                T2 = as.numeric(T2),
                                                P2 = as.numeric(P2))
  
  
  ## Replicate and weight Data - need four datasets because responders consistent with four DTRs
  # get empirical proportions for weight calculation
  pNP1_hat <- sum(data_output_relabel$P1 == 0)/nrow(data_output_relabel)
  pNP2_hat <- sum(data_output_relabel$P2 == 0 & data_output_relabel$R == 0)/sum(data_output_relabel$R == 0)
  pTheta_A_hat <- sum(data_output_relabel$T1 == 1 & data_output_relabel$P1 == 1) / sum(data_output_relabel$P1 == 1)
  pTheta_C_hat <- sum(data_output_relabel$T2 == 1 & data_output_relabel$R == 0 & data_output_relabel$P2 == 1) / sum(data_output_relabel$P2 == 1)
  
  data_output_relabel$w=rep(0)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A0A")] <- 1/(pNP1_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A0C1")] <- 1/(pNP1_hat*0.5*(1-pNP2_hat)*pTheta_C_hat)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A0D1")] <- 1/(pNP1_hat*0.5*(1-pNP2_hat)*(1-pTheta_C_hat))
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A0C0")] <- 1/(pNP1_hat*0.5*pNP2_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A0D0")] <- 1/(pNP1_hat*0.5*pNP2_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B0B")] <- 1/(pNP1_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B0C1")] <- 1/(pNP1_hat*0.5*(1-pNP2_hat)*pTheta_C_hat)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B0D1")] <- 1/(pNP1_hat*0.5*(1-pNP2_hat)*(1-pTheta_C_hat))
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B0C0")] <- 1/(pNP1_hat*0.5*pNP2_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B0D0")] <- 1/(pNP1_hat*0.5*pNP2_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B1B")] <- 1/((1-pNP1_hat)*(1-pTheta_A_hat))
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B1C1")] <- 1/((1-pNP1_hat)*(1-pTheta_A_hat)*(1-pNP2_hat)*pTheta_C_hat)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B1D1")] <- 1/((1-pNP1_hat)*(1-pTheta_A_hat)*(1-pNP2_hat)*(1-pTheta_C_hat))
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B1C0")] <- 1/((1-pNP1_hat)*(1-pTheta_A_hat)*pNP2_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="B1D0")] <- 1/((1-pNP1_hat)*(1-pTheta_A_hat)*pNP2_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A1A")] <- 1/((1-pNP1_hat)*pTheta_A_hat)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A1C1")] <- 1/((1-pNP1_hat)*pTheta_A_hat*(1-pNP2_hat)*pTheta_C_hat)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A1D1")] <- 1/((1-pNP1_hat)*pTheta_A_hat*(1-pNP2_hat)*(1-pTheta_C_hat))
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A1C0")] <-1/((1-pNP1_hat)*pTheta_A_hat*pNP2_hat*0.5)
  data_output_relabel$w[which(data_output_relabel$Treatment_Path=="A1D0")] <- 1/((1-pNP1_hat)*pTheta_A_hat*pNP2_hat*0.5)
  
  # 1st dataset of responders setting T2=1, P2=0
  datareps1 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps1$T2 <- 1
  datareps1$P2 <- 0
  
  
  # 2nd dataset of responders setting T2=1, P2=1
  datareps2 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps2$T2 <- 1
  datareps2$P2 <- 1
  
  # 3rd dataset of responders setting T2=-1, P2=0 
  datareps3 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps3$T2 <- -1
  datareps3$P2 <- 0
  
  
  # 4th dataset of responders setting T2=-1, P2=1
  datareps4 <- data_output_relabel[data_output_relabel$R==1,] 
  datareps4$T2 <- -1
  datareps4$P2 <- 1
  
  
  # dataset for non-responders
  datanoresp <- data_output_relabel[data_output_relabel$R==0,]
  
  # replicated data
  replicated_dat <- rbind(datareps1,datareps2,datareps3,datareps4,datanoresp)
  
  # create wrrm analysis data
  analysis_data <- replicated_dat %>% mutate(xalpha00 = (1-P1)*(1-P2),
                                             xalpha01 = (1-P1)*(P2),
                                             xalpha10 = (P1)*(1-P2),
                                             xalpha11 = (P1)*(P2),
                                             xbeta0 = T1*(1-P1),
                                             xbeta1 = T1*P1,
                                             xtheta0 = T2*(1-P2),
                                             xtheta1 = T2*P2,
                                             xgamma = (T1*T2)) %>%
    dplyr::select(id, xalpha00, xalpha01, xalpha10, xalpha11, xbeta0, xbeta1, xtheta0, xtheta1, xgamma, Y, w, Treatment_Path) 
  
  
  return(list(data_output_relabel, analysis_data, DTR_n))
  
  
}
