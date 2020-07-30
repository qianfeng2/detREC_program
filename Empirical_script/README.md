<script type="text/javascript" async
src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?
config=TeX-MML-AM_CHTML"
</script>

Mosaic Program to Detect Partial Alignments From Large Unaligned Sequences
-----------------------

### Step 1: Estimate major parameters $$\delta$$ and $$\epsilon$$.
Scripts for this step are shown in *_nd_mosaic_est_par.sh. 
Initial del and eps are descirbed in 1nd__mosaic_est_par.sh, same setting with previous pipeline when dealing with DBLa sequences collected from 10 countries (Tonkin-Hill, unpublished). 

After each iteration, the updated del and eps are processed by customed estimate_transition_probs_frm_viterbi.py code.

run example:

```
python scripts/estimate_transition_probs_frm_viterbi_qian.py --num_runs 578 --out iter1.txt --logfiles mosaic_processed_data/results_full_1/*.log --align mosaic_processed_data/align.txt 
```

Detailed implementation is showed in the following picture:

<p align="center">
<img src="https://github.com/qianfeng2/detREC_program/blob/master/Empirical_script/actual_implementation_step1.jpg" width="600" align="center">
</p>

Sup_ghana_3.Rmd provides the series of finally estimated del and eps collected in this step.


### Step 2: Estimate recombination parameter $$\rho$$.
A grid of numbers are used to get the optimal recombination par with maximum likelihood. The script is in Sup_ghana_3.Rmd


### Step 3: Obtain the partial alignment results 
mosaic_final_submit_qian.sh offers the script for final step to get maimum likelihood path for each target sequence with previously estimated parameters. 



### Reference
- Tonkin-Hill, G. Q. et al. (2019). Global structure of the var genes encoding the major variant surface antigen of Plasmodium falciparum. Nat. Commun., submitted.
