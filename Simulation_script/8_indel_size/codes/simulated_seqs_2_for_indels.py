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
# This code to rename sequences based on simulated sequences from INDELible. output file is in the data dir, called simulated_sequences.fasta
# One input par is tree
# One output fasta file is in data dir  
# Usage example: /Users/fengqian/anaconda2/bin/python /Users/fengqian/Downloads/simulated_seqs.py /Users/fengqian/Downloads/UniMelb_shared-master/project/mosaic_data/Protein_translateable_pilot_upper_centroids.fasta /Users/fengqian/simulated_tree.txt /Users/fengqian/simulated_seqs.fasta
#######################################################################
import sys, os
import glob
from mungo.fasta import FastaReader
from collections import defaultdict
input_fasta=sys.argv[1]
output_seqs=sys.argv[2]


seqs = {};seq_list = [];count=0
for h,s in FastaReader(input_fasta):
    seqs["seq" + str(count)] = s
    seq_list.append("seq" + str(count))
    count+=1
##organize the seq ID name
with open(output_seqs, 'w') as outfile:
    for s in seq_list:
        outfile.write(">"+s+"\n"+seqs[s]+"\n")  