# PRPP-SMART
Companion code for "A Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial Design Analyzed via Weighted and Replicated Frequentist and Bayesian Methods”.

## Usage Note
The code in this repository is in active development. To view or use stable code, see the appropriate releases:
- [v1.0.0](../../releases/tag/v1.0.0): Release accompanying initial submission to _Statistical Methods in Medical Research_.

## File Descriptions
- [sampleSize.R](sampleSize.R): Contains function to compute sample size for a SMART with a longitudinal outcome in which the primary aim is to compare two embedded DTRs.
- [generateSMART.R](generateSMART.R): Main function used to generate data from a SMART.
- [resultFunctions.R](resultFunctions.R): Helper functions for compiling results and tables

The [Stan](Stan) folder contains scripts to perform simulations under all scenarios compiled in the manuscript ([ArXiv](https://arxiv.org/abs/1810.13094)). Simulations are designed to be run in parallel in an environment where all available cores can be dedicated to this task; as such, we recommend appropriately modifying [simulateSMART.R](simulateSMART.R) before running simulation scripts on, say, a laptop rather than a high-performance cluster.
