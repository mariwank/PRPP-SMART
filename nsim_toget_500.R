# Functions to determine the number of simulations to run to get 500 


# Function to generate dataset 
generate_dataset <- function() {
  data <- generate_data(N=N, pNP1=pNP1, pTheta_A=pTheta_A, pNP2=pNP2, pTheta_C=pTheta_C, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
  return(data[[1]])
}

# Function to check if iteration should be skipped 
should_skip_iteration <- function(data) {
  
  # check to make sure at least three subjects per treatment path in PRPP-SMART data if not skip simulation
  trialpath_df <- data %>% dplyr::group_by(Treatment_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
  
  if (nrow(trialpath_df) < 20) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Main function to count iterations
count_iterations <- function() {
  total_datasets_generated <- 0
  iteration_number <- 0
  while (total_datasets_generated < 500) {
    iteration_number <- iteration_number + 1
    # Generate dataset
    set.seed(1000 + iteration_number)
    dataset <- generate_dataset()
    # Check if the iteration should be skipped
    if (should_skip_iteration(data=dataset)) {
      next  # Skip this iteration
    }
    
    # Increment total datasets generated
    total_datasets_generated <- total_datasets_generated + 1
  }
  
  return(iteration_number)
}