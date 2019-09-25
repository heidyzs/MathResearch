# Math Research

Math research files used to conduct SUS Analysis.

## File Names
* Bayes_Full_Changing_Mu_Sig_Gamma.R is r script used to create algorithim using the "full Bayes" method. Analysis is included using "Full_bayes_with_sdprior_20_v2.csv"
* Bayes_Visualizations.R is r script used to create plots/summarize data for all bayes methods.
* Bayes_vs_Bootstrap_CIs.R is r script used for bayes method compared with bootstrap methods. Analysis of results and plots are also included using "bayes_vs_bootstrap_full_data_set.csv"
* BetaSimulation.R is the r script used to run boostrap CIs using a beta distribution.
* Mixture_Distribution_Changing_Means.R is r script used to create simulation using a mixture distribution and changing the means.
*MLE_Bayes.R is r script file containing bayes analysis but also using the MLE to estimate likelihood parameters. Uses truncated normal distribution to fit likelihood.
*MLE_Mixture_Distributions.R is r script containing mixture distribution analysis but using MLE to estimate parameters. Uses truncated normal distribution to fit likelihood. 
* Skew_Normal_Serial.R is r script containing bootstrap CI analysis using a skew normal distribution. Serial version that is slower.
* Skew_Normal_Vectorized_Parallel.R is r script containing bootstrap CI analysis using a skew normal distribution. Uses vectorized parallel version that is faster to run. 
