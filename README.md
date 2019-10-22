# Math Research

All important files used to conduct SUS Analysis.

## File Names

### Folders
* `Plots`: folder contains plots generated from most of the results. Includes plots that were generated from the rmarkdown files. MLE_Results.Rmd plots are not in this folder.
* `Results`: folder contains all .csv files generated from various r scripts, all files are used to generate plots in rmarkdown files. 
  * `skewnormal_results (as of 5 AUG 2019).xlsx` is the only .xlsx file that contains final results (cover percentage, width, overcoverage, and undercoveage) from analyzing bootstrap vs clt methods using a skewnormal distribution.

### Rmarkdown and HTML Files
* `10 SEP Meeting Slides.Rmd`: rmarkdown file used to build all slides and plots used for the 10 SEP meeting. (Meeting conducted with NAG)
* `10-SEP-Meeting-Slides.html`: html file of the slides generated from the rmd file.
* `23 AUG Meeting Slides.Rmd`: rmarkdown file used to build all slides and plots used for the 23 AUG meeting.
* `23-AUG-Meeting-Slides.html`: html file of the slides generated from the rmd file.
* `MLE Results.Rmd`: rmarkdown file used to build all slides and plots containing results utilizing MLE. These plots include using a skew normal distribution and using the MLE to estimate all three parameters, using a skew normal distribution and fixing gamma and only using the MLE to estimate two parameters, using a truncated normal distribution (no gamma) and using MLE to estimate the parameters, and using mixture distributions fitted to a truncated normal distribution and using MLE to estimate the parameters. 
* `MLE-Results.html`: html file of the slides generated from the rmd file.

### R Scripts
* `Bayes_Full_Changing_Mu_Sig_Gamma.R`: script used to create algorithim using the "full Bayes" method. Analysis is included using `full_bayes_with_sdprior_20_v2.csv`
* `Bayes_Visualizations.R`: script used to create plots/summarize data for all bayes methods.
* `Bayes_vs_Bootstrap_CIs.R`: script used for bayes method compared with bootstrap methods. Analysis of results and plots are also included using `bayes_vs_bootstrap_full_data_set.csv`
* `BetaSimulation.R`: r script used to run boostrap CIs using a beta distribution.
* `Mixture_Distribution_Changing_Means.R`: script used to create simulation using a mixture distribution and changing the means.
* `MLE_Bayes.R`: script containing bayes analysis but also using the MLE to estimate likelihood parameters. Uses truncated normal distribution to fit likelihood. Script can be modified to use a skew normal distribution by changing the `log.lik` function withint the script.
* `MLE_Mixture_Distributions.R`: script containing mixture distribution analysis but using MLE to estimate parameters. Uses truncated normal distribution to fit likelihood. Script can be modified to use a skew normal distribution by changing the `log.lik` function withint the script. 
* `Skew_Normal_Serial.R`: script containing bootstrap CI analysis using a skew normal distribution. Serial version that is slower.
* `Skew_Normal_Vectorized_Parallel.R`: script containing bootstrap CI analysis using a skew normal distribution. Uses vectorized parallel version that is faster to run. 


