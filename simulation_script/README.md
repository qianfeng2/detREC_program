Simulate Sequences for Evaluating Performance of Proposed Recombination Detection Program
-----------------------
[![Python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)](https://www.python.org/download/releases/2.7/)
[![Python 3.6](https://img.shields.io/pypi/pyversions/Django)](https://www.python.org/downloads/release/python-360/)

-- Require modules for Python 3 user:  
Bio  
pandas  
scipy  
numpy 

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

Now you have simulated a number of sequences, sequence IDs start with "r_" are all generated recombinants, start with "seq" are non-recombinants. Afterwards you are able to proceed to the general recombination detection step with the help of [manual](https://github.com/qianfeng2/detREC_program).
