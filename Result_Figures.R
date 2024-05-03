# install/load libraries 
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
p_load(ggplot2, tidyverse, ggsci, ggpubr, gtools, utils) 

# Scenarios 
scenario_vec <- c(1,2,3)

## Load in Frequentist WRRM Results 
# path to location of frequentist WRRM results per scenario 
freq_folder_path <- "Your path here"
# Example
freq_folder_path <- "/Users/mariwank/Documents/UMICH/dissertation/WRRM_Paper/SimResults/Frequentist/Scenario"

for (i in scenario_vec){
  filenames <- list.files(paste0(freq_folder_path, i), pattern="*.csv", full.names=TRUE)
  sort(filenames)
  ldf <- lapply(filenames, read_csv)
  
  assign(paste0("ci_results_full",i, "_freq"), ldf[[1]])
  assign(paste0("ci_results_trad",i, "_freq"), ldf[[2]])
  assign(paste0("dtr_results_full",i, "_freq"), ldf[[5]])
  assign(paste0("dtr_results_trad",i, "_freq"), ldf[[6]])
  
}


## Load in Bayesian WRRM Results 
# path to location of Bayesian WRRM results per scenario 
bayes_folder_path <- "Your path here"
# Example
bayes_folder_path <- "/Users/mariwank/Documents/UMICH/dissertation/WRRM_Paper/SimResults/Bayesian/Scenario"

for (i in scenario_vec){
  filenames <- list.files(paste0(bayes_folder_path, i), pattern="*.csv", full.names=TRUE)
  sort(filenames)
  ldf <- lapply(filenames, read_csv)
  
  assign(paste0("ci_results_full",i, "_bayes"), ldf[[1]])
  assign(paste0("ci_results_trad",i, "_bayes"), ldf[[2]])
  assign(paste0("dtr_results_full",i, "_bayes"), ldf[[5]])
  assign(paste0("dtr_results_trad",i, "_bayes"), ldf[[6]])
  
}

## Create frequentist plotting dataframe 

# full frequentist WRRM, scenario 1
plot_df1_freq <- data.frame(DTR = dtr_results_full1_freq$DTR[1:4], 
                                 Bias = dtr_results_full1_freq$Bias[1:4],
                                 Sd = dtr_results_full1_freq$SE[1:4],
                                 Avg_sd = dtr_results_full1_freq$Avg_se[1:4],
                                 rMSE = dtr_results_full1_freq$rMSE[1:4],
                                 Coverage = ci_results_full1_freq$CI_Coverage[1:4],
                                 Method = "Full WRRM",
                                 Scenario = 1)

# traditional frequentist WRRM, scenario 1
plot_df1t_freq <- data.frame(DTR = dtr_results_full1_freq$DTR[1:4], 
                                  Bias = dtr_results_trad1_freq$Bias,
                                  Sd = dtr_results_trad1_freq$SE,
                                  Avg_sd = dtr_results_trad1_freq$Avg_se,
                                  rMSE = dtr_results_trad1_freq$rMSE,
                                  Coverage = ci_results_trad1_freq$CI_Coverage,
                                  Method = "Traditional WRRM",
                                  Scenario = 1)

# full frequentist WRRM, scenario 2
plot_df2_freq <- data.frame(DTR = dtr_results_full1_freq$DTR[1:4], 
                                 Bias = dtr_results_full2_freq$Bias[1:4],
                                 Sd = dtr_results_full2_freq$SE[1:4],
                                 Avg_sd = dtr_results_full2_freq$Avg_se[1:4],
                                 rMSE = dtr_results_full2_freq$rMSE[1:4],
                                 Coverage = ci_results_full2_freq$CI_Coverage[1:4],
                                 Method = "Full WRRM",
                                 Scenario = 2)

# traditional frequentist WRRM, scenario 2
plot_df2t_freq <- data.frame(DTR = dtr_results_full1_freq$DTR[1:4], 
                                  Bias = dtr_results_trad2_freq$Bias,
                                  Sd = dtr_results_trad2_freq$SE,
                                  Avg_sd = dtr_results_trad2_freq$Avg_se,
                                  rMSE = dtr_results_trad2_freq$rMSE,
                                  Coverage = ci_results_trad2_freq$CI_Coverage,
                                  Method = "Traditional WRRM",
                                  Scenario = 2)

# full frequentist WRRM, scenario 3
plot_df3_freq <- data.frame(DTR = dtr_results_full1_freq$DTR[1:4], 
                                 Bias = dtr_results_full3_freq$Bias[1:4],
                                 Sd = dtr_results_full3_freq$SE[1:4],
                                 Avg_sd = dtr_results_full3_freq$Avg_se[1:4],
                                 rMSE = dtr_results_full3_freq$rMSE[1:4],
                                 Coverage = ci_results_full3_freq$CI_Coverage[1:4],
                                 Method = "Full WRRM",
                                 Scenario = 3)

# traditional frequentist WRRM, scenario 3
plot_df3t_freq <- data.frame(DTR = dtr_results_full1_freq$DTR[1:4], 
                                  Bias = dtr_results_trad3_freq$Bias,
                                  Sd = dtr_results_trad3_freq$SE,
                                  Avg_sd = dtr_results_trad3_freq$Avg_se,
                                  rMSE = dtr_results_trad3_freq$rMSE,
                                  Coverage = ci_results_trad3_freq$CI_Coverage,
                                  Method = "Traditional WRRM",
                                  Scenario = 3)

# All scenarios full and traditional frequentist WRRM results 
plot_df_freq <- rbind(plot_df1_freq,
                           plot_df1t_freq,
                           plot_df2_freq,
                           plot_df2t_freq, 
                           plot_df3_freq,
                           plot_df3t_freq)

plot_df_freq$abs_bias <- abs(plot_df_freq$Bias)                            
plot_df_freq$Coverage <- round(plot_df_freq$Coverage,2)

## Create Bayesian plotting dataframe 

# full Bayesian WRRM, scenario 1
plot_df1_bayes <- data.frame(DTR = dtr_results_full1_bayes$DTR[1:4], 
                                  Bias = dtr_results_full1_bayes$Bias_adj[1:4],
                                  Sd = dtr_results_full1_bayes$Sd_adj[1:4],
                                  Avg_sd = dtr_results_full1_bayes$Avg_sd_adj[1:4],
                                  rMSE = dtr_results_full1_bayes$rMSE_adj[1:4],
                                  Coverage = ci_results_full1_bayes$CI_Coverage_adj[1:4],
                                  Method = "Full WRRM",
                                  Scenario = 1)

# traditional Bayesian WRRM, scenario 1
plot_df1t_bayes <- data.frame(DTR = dtr_results_full1_bayes$DTR[1:4], 
                                   Bias = dtr_results_trad1_bayes$Bias_adj,
                                   Sd = dtr_results_trad1_bayes$Sd_adj,
                                   Avg_sd = dtr_results_trad1_bayes$Avg_sd_adj,
                                   rMSE = dtr_results_trad1_bayes$rMSE_adj,
                                   Coverage = ci_results_trad1_bayes$CI_Coverage_adj,
                                   Method = "Traditional WRRM",
                                   Scenario = 1)

# full Bayesian WRRM, scenario 2
plot_df2_bayes <- data.frame(DTR = dtr_results_full1_bayes$DTR[1:4], 
                                  Bias = dtr_results_full2_bayes$Bias_adj[1:4],
                                  Sd = dtr_results_full2_bayes$Sd_adj[1:4],
                                  Avg_sd = dtr_results_full2_bayes$Avg_sd_adj[1:4],
                                  rMSE = dtr_results_full2_bayes$rMSE_adj[1:4],
                                  Coverage = ci_results_full2_bayes$CI_Coverage_adj[1:4],
                                  Method = "Full WRRM",
                                  Scenario = 2)

# traditional Bayesian WRRM, scenario 2
plot_df2t_bayes <- data.frame(DTR = dtr_results_full1_bayes$DTR[1:4], 
                                   Bias = dtr_results_trad2_bayes$Bias_adj,
                                   Sd = dtr_results_trad2_bayes$Sd_adj,
                                   Avg_sd = dtr_results_trad2_bayes$Avg_sd_adj,
                                   rMSE = dtr_results_trad2_bayes$rMSE_adj,
                                   Coverage = ci_results_trad2_bayes$CI_Coverage_adj,
                                   Method = "Traditional WRRM",
                                   Scenario = 2)

# full Bayesian WRRM, scenario 3
plot_df3_bayes <- data.frame(DTR = dtr_results_full1_bayes$DTR[1:4], 
                                  Bias = dtr_results_full3_bayes$Bias_adj[1:4],
                                  Sd = dtr_results_full3_bayes$Sd_adj[1:4],
                                  Avg_sd = dtr_results_full3_bayes$Avg_sd_adj[1:4],
                                  rMSE = dtr_results_full3_bayes$rMSE_adj[1:4],
                                  Coverage = ci_results_full3_bayes$CI_Coverage_adj[1:4],
                                  Method = "Full WRRM",
                                  Scenario = 3)

# traditional Bayesian WRRM, scenario 3
plot_df3t_bayes <- data.frame(DTR = dtr_results_full1_bayes$DTR[1:4], 
                                   Bias = dtr_results_trad3_bayes$Bias_adj,
                                   Sd = dtr_results_trad3_bayes$Sd_adj,
                                   Avg_sd = dtr_results_trad3_bayes$Avg_sd_adj,
                                   rMSE = dtr_results_trad3_bayes$rMSE_adj,
                                   Coverage = ci_results_trad3_bayes$CI_Coverage_adj,
                                   Method = "Traditional WRRM",
                                   Scenario = 3)

# All scenarios full and traditional Bayesian WRRM results 
plot_df_bayes <- rbind(plot_df1_bayes,
                            plot_df1t_bayes,
                            plot_df2_bayes,
                            plot_df2t_bayes, 
                            plot_df3_bayes,
                            plot_df3t_bayes)


plot_df_bayes$abs_bias <- abs(plot_df_bayes$Bias)                            
plot_df_bayes$Coverage <- round(plot_df_bayes$Coverage,2)


######################### PLOTS ################################

## Frequentist WRRM Indifference DTR Result Plots

# Absolute Bias
ggbarplot(plot_df_freq, x = "DTR", y = "abs_bias", fill = "Method",
          color = "white", sorting = "none",          
          palette = "grey", rotate = FALSE, position = position_dodge(0.7), 
          x.text.angle = 90, facet.by = "Scenario", ggtheme = theme_pubr(),
          main = "Absolute Bias of Indifference DTRs per Scenario", subtitle = "Frequentist Estimation",ylab = "Absolute Bias", 
          panel.labs = list(Scenario = c("Scenario 1", "Scenario 2", "Scenario 3"))) + 
  scale_x_discrete(labels = c("AAC", "AAD", "BBC", "BBD")) +
  theme(plot.title=element_text(size=14), axis.text=element_text(size=14), axis.title=element_text(size=15), legend.title=element_text(size=14), legend.text=element_text(size=14), strip.text=element_text(size=14))


# rMSE
ggbarplot(plot_df_freq, x = "DTR", y = "rMSE", fill = "Method",
          color = "white", sorting = "none",          
          palette = "grey", rotate = FALSE, position = position_dodge(0.7), 
          x.text.angle = 90, facet.by = "Scenario", ggtheme = theme_pubr(),
          main = "rMSE of Indifference DTRs per Scenario", subtitle = "Frequentist Estimation", ylab = "rMSE", 
          panel.labs = list(Scenario = c("Scenario 1", "Scenario 2", "Scenario 3"))) + 
  scale_x_discrete(labels = c("AAC", "AAD", "BBC", "BBD")) +
  theme(plot.title=element_text(size=14), axis.text=element_text(size=14), axis.title=element_text(size=15), legend.title=element_text(size=14), legend.text=element_text(size=14), strip.text=element_text(size=14))


## Bayesian WRRM Indifference DTR Result Plots

# Absolute Bias
ggbarplot(plot_df_bayes, x = "DTR", y = "abs_bias", fill = "Method",
          color = "white", sorting = "none",          
          palette = "grey", rotate = FALSE, position = position_dodge(0.7), 
          x.text.angle = 90, facet.by = "Scenario", ggtheme = theme_pubr(),
          main = "Absolute Bias of Indifference DTRs per Scenario", subtitle = "Bayesian Estimation",ylab = "Absolute Bias", 
          panel.labs = list(Scenario = c("Scenario 1", "Scenario 2", "Scenario 3"))) + 
  scale_x_discrete(labels = c("AAC", "AAD", "BBC", "BBD")) +
  theme(plot.title=element_text(size=14), axis.text=element_text(size=14), axis.title=element_text(size=15), legend.title=element_text(size=14), legend.text=element_text(size=14), strip.text=element_text(size=14))


# rMSE
ggbarplot(plot_df_bayes, x = "DTR", y = "rMSE", fill = "Method",
          color = "white", sorting = "none",          
          palette = "grey", rotate = FALSE, position = position_dodge(0.7), 
          x.text.angle = 90, facet.by = "Scenario", ggtheme = theme_pubr(),
          main = "rMSE of Indifference DTRs per Scenario", subtitle = "Bayesian Estimation", ylab = "rMSE", 
          panel.labs = list(Scenario = c("Scenario 1", "Scenario 2", "Scenario 3"))) + 
  scale_x_discrete(labels = c("AAC", "AAD", "BBC", "BBD")) +
  theme(plot.title=element_text(size=14), axis.text=element_text(size=14), axis.title=element_text(size=15), legend.title=element_text(size=14), legend.text=element_text(size=14), strip.text=element_text(size=14))

