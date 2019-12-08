A Program to Detect Recombinants From Unaligned Sequences
-----------------------

### About
This program is a novel approach for detecting recombinant sequences and corresponding statistical support values from unaligned biological sequences. This framework develops on the basis of the paritial alignment results from jumping hidden markov model (JHMM, or mosaic), after that, by dividing them into multiple equal-length triples, on which we use a new distance-based procedure to identify recombinant from each triple. Statistical support values calculated from Bootstrap, the bigger the better, indicating the robustness of identified recombinants.


### Required softwares
- MAFFT used to align one sequence to another two pre_aligned sequences (https://mafft.cbrc.jp/alignment/software/)
- SeqKit used to concatenate the two segments for each triple (https://bioinf.shenwei.me/seqkit/download/)
- Python  
-- Require modules for Python 2 user:  
mungo (`pip install git+https://github.com/PapenfussLab/Mungo`)  
Bio  
pandas  
scipy  
numpy  
-- Require modules for Python 3 user:  
Bio  
pandas  
scipy  
numpy 

### Optional softwares (only used for simulation section)
- [Msprime](https://msprime.readthedocs.io/en/stable/installation.html) used to generate one phylogenetic tree topology
- INDELible used for simulation of biological sequences when considering indel events (Fletcher et al., 2009)
- Python module  
Pyvolve used for simulating biological sequences given a tree topology (Spielman et al. 2015) 
sklearn
statistics
csv
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) for reproducible research pipeline 

### Required Input Files 
- Patial alignment produced by JHMM (please see [MZmosaic](https://github.com/qianfeng2/detREC_program/tree/master/MZmosaic) sub folder)
- Input fasta format biological sequences, maximum length for identifiers length is 15 (Zilversmit et al., 2013)

### Creating Input File

### Run Example 

```
cd /Users/fengqian/MZmosaic
./mosaic -estimate -seq input.fasta  -rec 0 -aa -tag middle_file
delta=$(grep -o 'Gap initiation: .*$' middle_file_align.txt | cut -c17-)
epsl=$(grep -o 'Gap extension:  .*$' middle_file_align.txt | cut -c17-)
./mosaic -seq input.fasta -del $delta -eps $epsl -aa -tag output -grid 0.001 0.010 10 1
```

### Run Example for large number of sequences (>10000 sequences)

Instead of complete and time-consuming Baum-Welch algorithm to estimate gap open and gap extension probabilities, the slightly less accurate but much faster viterbi training algorithm was used. Each iteration was run as on a high performing computing cluster (Helix)

Iterate until convergence:

1) Choose an initial set of parameters
2) Compute the Viterbi paths of all sequences
3) Count frequencies of events and calculate new parameters
4) Update -> 1) 
5) Stop when the major parameters del and eps change by less than 1%.

The empirical Ghana pilot DBL\alpha dataset analyzed in manuscript contain more than 17000 sequences, the detailed code for generating partial alignment results are displayed in [Empirical_script](https://github.com/qianfeng2/detREC_program/tree/master/Empirical_script) sub folder.


### Running recombination detection program
#### Required:
- "\<filename\>": name of the file which contains your preprepared partial alignment information.
- "\<filename\>": name of the file which contains your complete biological sequences, fasta format required.
- "\<filename\>": name of the file which contains your list of identified recombinants and related statistical support value, csv format required.

Type in above three input parameters in order.


### Output
Produces a series of files based on various stage of the implementation of recombination detection program, and places them in the directory specified by output.

- temp file folder  
This folder provides all the chunks containing original triple and MAFFT processed fasta files.
- complement_chunks file folder  
This folder provides all the equal-length triples, name of each file indicates chunk index, two adjacent segment indices, and identified bkp in this triple.
- output csv file:  
Each row records the chunk index in partial alignment result, target, db1 and db2 are three sequences ID for each triple, rec is the identified recombinant ID from this specific triple, sv is the support value. For instance:  

| chunk        | target  | db1  | db2  | rec  | sv  |
| ------------|------------|------------|------------|------------|------------|
|2 | seq3|seq5|seq2|seq3|1|
|6 | seq7|seq8|seq1|seq8|0.58|
|... | ... |... |... |... |... |



### Run Example

```
python integrated_rec_det.py output_align.txt input.fasta output.csv
```
[Test_files](https://github.com/qianfeng2/detREC_program/tree/master/Test_files) sub folder, as a toy example, provides a test input.fasta and all the middle and final output files.

### Reference
- Martine M Zilversmit et al. "Hypervariable antigen genes in malaria have ancient roots". In: BMC evolutionary biology 13.1 (2013), p. 110.
- Fletcher, W., & Yang, Z. (2009). INDELible: a flexible simulator of biological sequence evolution. Molecular biology and evolution, 26(8), 1879-1888.
- Spielman, S. J., & Wilke, C. O. (2015). Pyvolve: a flexible Python module for simulating sequences along phylogenies. PloS one, 10(9), e0139047.
