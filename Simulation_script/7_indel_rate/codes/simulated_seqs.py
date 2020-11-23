#######################################################################
# Copyright (C) 2019  Qian Feng

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Usuage instruction(only suitable for python 2.7)
# This code to simulate 300 sequences based on simulated tree. output file is in the data dir, called simulated_sequences.fasta
# One input par is tree
# One output fasta file is in data dir  
# Usage example: /Users/fengqian/anaconda2/bin/python /Users/fengqian/Downloads/simulated_seqs.py /Users/fengqian/Downloads/UniMelb_shared-master/project/mosaic_data/Protein_translateable_pilot_upper_centroids.fasta /Users/fengqian/simulated_tree.txt /Users/fengqian/simulated_seqs.fasta
#######################################################################
import sys, os
import pyvolve
import glob
from mungo.fasta import FastaReader
from collections import defaultdict
input_fasta=sys.argv[1]
input_tree_txt=sys.argv[2]
output_seqs=sys.argv[3]

#f = pyvolve.ReadFrequencies("amino_acid", file = "/Users/fengqian/Downloads/UniMelb_shared-master/project/mosaic_data/Protein_translateable_pilot_upper_centroids.fasta")
#f = pyvolve.ReadFrequencies("amino_acid", file = "/data/cephfs/punim0609/qian_feng/snake_pipeline/data/Protein_translateable_pilot_upper_centroids.fasta")
f = pyvolve.ReadFrequencies("amino_acid", file = input_fasta)
frequencies = f.compute_frequencies()
my_tree_1 = pyvolve.read_tree(file=input_tree_txt)
my_model_1 = pyvolve.Model("WAG", {"state_freqs":frequencies} )
my_partition_1 = pyvolve.Partition(models = my_model_1, size = 200)
my_evolver_1 = pyvolve.Evolver(partitions = my_partition_1, tree = my_tree_1)
my_evolver_1(ratefile = None, infofile =  None, seqfile = output_seqs)

seqs = {};seq_list = [];count=0
for h,s in FastaReader(output_seqs):
    seqs["seq" + str(count)] = s
    seq_list.append("seq" + str(count))
    count+=1
##organize the seq ID name
with open(output_seqs, 'w') as outfile:
    for s in seq_list:
        outfile.write(">"+s+"\n"+seqs[s]+"\n")  
