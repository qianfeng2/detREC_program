
Mosaic Program to Detect Partial Alignments From Large Unaligned Sequences
-----------------------
This folder provides specific scripts of running mosaic program for a [dataset](https://github.com/qianfeng2/detREC_program/tree/master/Empirical_script/Ghana_dataset) of DBLa sequences from a cross-sectional study in Ghana. 

This entire pipeline we adopted are from an accepted paper: Tonkin-Hill et al. (2021), where they managed to handle a larger dataset. If you aim to get partial alignments results for your own large dataset, please take a look at their detailed pipeline (https://github.com/gtonkinhill/global_var_manuscript), also please cite: Tonkin-Hill G, Ruybal-Pesántez S, Tiedje KE, Rougeron V, Duffy MF, Zakeri S,et al. Evolutionary analyses of the major variant surface antigen-encoding genesreveal population structure of *Plasmodium falciparum* within and between continents. PLoS Genetics. 2021;17(2):e1009269.


If you only want to reproduce the results in our related manuscript (https://www.biorxiv.org/content/10.1101/2020.11.18.389262v1), please follow the instructions below:

### Step 1: Estimate major parameters del and eps.
We first split the data set into subset to be searched in parallel on our computing cluster. This is showed in the following picture:

<p align="center">
<img src="https://github.com/qianfeng2/detREC_program/blob/master/Empirical_script/implementation_step1.jpg" width="600" align="center">
</p>

Now we perform the Viterbi training algorithm.

Scripts for this step are shown in *_nd_mosaic_est_par.sh. Initial del and eps are descirbed in 1nd__mosaic_est_par.sh, second iteration script is 2nd_mosaic_est_par.sh, third iteration script is 3nd_mosaic_est_par.sh... The algorithm converged after 8 iterations for our dataset.

The updated del and eps are processed by estimate_transition_probs_frm_viterbi.py. 

run example:

```
python estimate_transition_probs_frm_viterbi.py --num_runs 578 --out iter1.txt --align mosaic_processed_data/results_full_1/*.txt 
```
Please note that:

(1) This python script is from Tonkin-Hill et al. (2021), and was written by Python 2. 

(2) After each iteration, you will get a log file and align.txt file per fasta file. In order for you to run this script smoothly, please put all log files and align.txt files into the same directory. Make sure the prefix of log and align files are the same. 

Sup_ghana_3.Rmd provides the series of estimated del and eps with iteration in this step.


### Step 2: Estimate recombination parameter.
A grid of numbers are used to get the optimal recombination par with maximum likelihood. The script is in Sup_ghana_3.Rmd


### Step 3: Obtain the partial alignment results.
mosaic_final_submit_qian.sh offers the script for final step to get maimum likelihood path for each target sequence with previously estimated parameters. 



### Reference
- Tonkin-Hill G, Ruybal-Pesántez S, Tiedje KE, Rougeron V, Duffy MF, Zakeri S,et al. Evolutionary analyses of the major variant surface antigen-encoding genesreveal population structure of *Plasmodium falciparum* within and between continents. PLoS Genetics. 2021;17(2):e1009269.

