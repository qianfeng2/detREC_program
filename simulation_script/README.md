Simulate Sequences for Evaluating Performance of Proposed Recombination Detection Program
-----------------------

### Step 1: Simuate one tree topology
One option in script msprime_tree.py is sample_size=200, providing the number of tip leaves in this tree.

#### Run Example 

```
python codes/msprime_tree.py data/simulated_tree.txt
```



### Step 2: Evolve sequences

Note that equilibrium frequencies for amino acids are from the empirical Ghana pilot DBLa datasets Protein_translateable_pilot_upper_centroids.fasta, you can customize it for your purpose.

#### Run Example 

```
python codes/simulated_seqs.py data/Protein_translateable_pilot_upper_centroids.fasta data/simulated_tree.txt data/simulated_seqs.fasta
```


### Step 3: Generate recombinant sequences
The parents for generating recombinants should be removed for following recombination detection
Four parameters are needed in this script, requirements for them are illustrated at the top of this python codes.
#### Run Example 

```
python codes/recombined_seqs.py 0 50 0 data/simulated_seqs.fasta
```

Now you have simulated a number of sequences, sequence ID starts with "r_" are all generated recombinants, starts with "seq" are non-recombinants. Afterwards you are able to proceed to the general recombination detection step with the help of [manual](https://github.com/qianfeng2/detREC_program).