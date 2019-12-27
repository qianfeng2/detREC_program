Mosaic Program to Detect Mosaic From Large Unaligned Sequences
-----------------------

### Step 1: Estimate major parameters del and eps
Scripts for this step are shown in *_nd_mosaic_est_par.sh. 
Initial del and eps are descirbed in 1nd__mosaic_est_par.sh, same setting with previous pipeline when dealing with DBLa sequences collected from 10 countries (Tonkin-Hill, unpublished).


### Step 2: Estimate recombination parameter
A grid of numbers are used to get the optimal recombination par with maximum likelihood. The script is in Sup_ghana_2.Rmd


### Step 3: Obtain the partial alignment results 
mosaic_final_submit_qian.sh offers the script for final step to get maimum likelihood path for each target sequence with previously estimated parameters. 


File Sup_ghana_3.Rmd provides the series of estimated del and eps from step 1 results, the process of getting final recombination parameter from step 2 results.


### Reference
- Tonkin-Hill, G. Q. et al. (2019). Global structure of the var genes encoding the major variant surface antigen of Plasmodium falciparum. Nat. Commun., submitted.
