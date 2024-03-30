# PRPP-SMART
Companion code for "A Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial Design Analyzed via Weighted and Replicated Frequentist and Bayesian Methods”.

## Usage Note
The code in this repository is in active development. To view or use stable code, see the appropriate releases:
- [v1.0.0](../../releases/tag/v1.0.0): Release accompanying initial submission to _Statistical Methods in Medical Research_.

## File Descriptions
- [DataGeneration.R](DataGeneration.R): Main function used to generate data from a two-stage PRPP-SMART with binary end-of-stage outcomes.
- [Frequentist_WRRM.R](Frequentist_WRRM.R): Simulation code for frequentist WRRM under all scenarios considered in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)).
- [Bayesian_WRRM.R](Bayesian_WRRM.R): Main function used to generate data from a two-stage PRPP-SMART with binary end-of-stage outcomes.
- [resultFunctions.R](resultFunctions.R): Helper functions for compiling results and tables

## Folder Descriptions
- The [WilliamsSavitsky2021_Functions](WilliamsSavitsky2021_Functions) folder contains the functions to implement the method outlined in Williams, M. R., and Savitsky, T. D. (2021) Uncertainty Estimation for Pseudo-Bayesian Inference Under Complex Sampling, https://doi.org/10.1111/insr.12376. Note, these functions are also available in the r package [csSampling](https://github.com/RyanHornby/csSampling). 

- The [Stan](Stan) folder contains the stan models used to perform the full and traditional Bayesian weighted and replicated regression models (WRRMs) outlined in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)). Simulations are designed to be run in parallel in an environment where all available cores can be dedicated to this task; as such, we recommend appropriately modifying [simulateSMART.R](simulateSMART.R) before running simulation scripts on, say, a laptop rather than a high-performance cluster.
