#### Code and data for Strauss ED, Shizuka D, Holekamp KE Juvenile rank acquisition paper

Strauss ED, Shizuka D, Holekamp KE (in review) Juvenile rank acquisition is associated with fitness independent of adult rank. *Proceedings of the Royal Society, B*. 

Direct inquiries to Eli Strauss straussed@gmail.com

#### Overview
Files are numbered according to the order of their appearance in the analysis. File paths in the scripts will need to be changed in order for the code to run properly. Below is a list of files in the order they should be run to reproduce the analysis. 


* __00.define_functions.R__	Define functions used in the analysis

* __01.calculate_den_period.R__ Gather data on relevant individuals and calculate den period. 

* __02.cohortInfo.RData__	Data on all individuals used in the analysis

* __03.calculate_elo_deviance.R__	Calculate Elo deviance scores and other covariates

* __04.cub_dev_vars.RData__	Data on relevant individuals with Elo deviance score and other covariates

* __05.survival_analysis.R__	Run models of survival and LRS and generate plots and tables

* __06.rank.acquisition.RData__ Data on individuals that survive to den independence with survival covariates added

* __07.development_of_rank.R__	Assess timing of rank development for figure in supplemental materials

* __08.calculate_elo_deviance_k100.R__ Duplicate of 03, except Elo now calculated with different parameterization (k = 100)

* __09.cub_dev_vars_k100.RData__	Duplicate of 04, except Elo now calculated with different parameterization (k = 100)

* __10.survival_analysis_k100.R__	Duplicate of 05, except Elo now calculated with different parameterization (k = 100)

* __11.rank.acquisition.k100.RData__	Duplicate of 06, except Elo now calculated with different parameterization (k = 100)

* __12.development_of_rank_k100.R__	Duplicate of 07, except Elo now calculated with different parameterization (k = 100)


