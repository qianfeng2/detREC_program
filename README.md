A Program to Detect Recombinants From Unaligned Sequences
-----------------------

### About
This program is a novel approach for detecting recombinant sequences and corresponding statistical support values from unaligned biological sequences. This framework develops on the basis of the paritial alignment results from jumping hidden markov model (JHMM, or mosaic), after that, by dividing them into multiple equal-length triples, on which we use a new distance-based procedure to identify recombinant from each triple. Statistical support values calculated from Bootstrap, the bigger the better, indicating the robustness of identified recombinants.


### Required softwares
- MAFFT used to align one sequence to another two pre_aligned sequences (https://mafft.cbrc.jp/alignment/software/)
- SeqKit used to concatenate the two segments for each triple (https://bioinf.shenwei.me/seqkit/download/)
- Python  
-- Require modules for Python 2 user:  
-- Require modules for Python 3 user:

### Optional softwares (only used for simulation section)
- Msprime used to generate one arbitrary phylogenetic tree (https://msprime.readthedocs.io/en/stable/installation.html)
- Snakemake
- INDELible
- Python module  

### Required Input Files 
- Patial alignment produced by JHMM (please see [MZmosaic](https://github.com/qianfeng2/detREC_program/tree/master/MZmosaic) folder)
- Input fasta format biological sequences, maximum length for identifiers length is 15 (Zilversmit et al., 2013)

### Creating Input File

### Run Example 

```
cd /Users/fengqian/MZmosaic
./mosaic -estimate -seq input.fasta  -rec 0 -aa -tag middle_file
delta=$(grep -o 'Gap initiation: .*$' middle_file_align.txt | cut -c17-)
epsl=$(grep -o 'Gap extension:  .*$' middle_file_align.txt | cut -c17-)
./mosaic -seq input.txt -del $delta -eps $epsl -aa -tag output -grid 0.001 0.010 10 1
```

### Run Example for large number of sequences (>10000 sequences)

Instead of complete and time-consuming Baum-Welch algorithm to estimate gap open and gap extension probabilities, the slightly less accurate but much faster viterbi training algorithm has been used. Each iteration was run as on a high performing computing cluster (Helix)

Iterate until convergence:

1) Choose an initial set of parameters
2) Compute the Viterbi paths of all sequences
3) Count frequencies of events and calculate new parameters
4) Update -> 1) 
5) Stop when the major parameters del and eps change by less than 1%.

The empirical Ghana pilot DBLa dataset analyzed in manuscript contain more than 17000 sequences, the detailed code for generating partial alignment results are displayed in [Empirical_script](https://github.com/qianfeng2/detREC_program/tree/master/Empirical_script) sub folder.


### Running recombination detection program
#### Required:
- "\<filename\>": name of the file which contains your preprepared partial alignment information.
- "\<filename\>": name of the file which contains your complete biological sequences, fasta format required.
- "\<filename\>": name of the file which contains your list of identified recombinants and related statistical support value, csv format required.

Typed in above three input parameters in order.


### Output
Produces a series of files based on various stage of the implementation of recombination detection program, and places them in the directory specified by output.

- temp file folder

- complement_chunks file folder

- output csv file:
Each row records the chunk index in partial alignment result, target, db1 and db2 are three sequences ID for each triple, rec is the identified recombinant ID from this specific triple, sv is the support value.  
eg:  




### Run Example

```
python integrated_rec_det.py output_align.txt input.fasta output.csv
```

### Reference
- Martine M Zilversmit et al. "Hypervariable antigen genes in malaria have ancient roots". In: BMC evolutionary biology 13.1 (2013), p. 110.

