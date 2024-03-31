# PRPP-SMART
Companion code for "A Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial Design Analyzed via Weighted and Replicated Frequentist and Bayesian Methods”.

## Usage Note
The code in this repository is in active development. To view or use stable code, see the appropriate releases:
- [v1.0.0](../../releases/tag/v1.0.0): Release accompanying initial submission to _Statistical Methods in Medical Research_.

## File Descriptions
- [DataGeneration.R](DataGeneration.R): Main function used to generate data from a two-stage PRPP-SMART with binary end-of-stage outcomes.
- [Frequentist_WRRM.R](Frequentist_WRRM.R): Simulation code for frequentist weighted and replicated regression models (WRRMs) to estimate embedded dynamic treatment regimes (DTRs) under all scenarios considered in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)). 
- [Bayesian_WRRM.R](Bayesian_WRRM.R): Simulation code for Bayesian weighted and replicated regression models (WRRMs) to estimate embedded dyamic treatment regimes (DTRs) under all scenarios considered in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)). Note, Bayesian simulations are reccomened to be run on a high-performance cluster rather than a laptop. 
- [Result_Figures.R](Result_Figures.R): Code for creating figures in manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)). 

## Folder Descriptions
- The [WilliamsSavitsky2021_Functions](WilliamsSavitsky2021_Functions) folder contains the functions to implement the method outlined in Williams, M. R., and Savitsky, T. D. (2021) Uncertainty Estimation for Pseudo-Bayesian Inference Under Complex Sampling, https://doi.org/10.1111/insr.12376. Note, these functions are also available in the r package [csSampling](https://github.com/RyanHornby/csSampling). 

- The [Stan](Stan) folder contains the stan models used to perform the full and traditional Bayesian weighted and replicated regression models (WRRMs) outlined in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)). 
