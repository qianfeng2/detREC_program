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

# Usuage instruction(suitable for python 3.5+)
# This code to replace the new tree into example control txt file
# two input pars: one is simulated_tree.txt, another is control.txt
# Outut is the modified control.txt
# Usage example: /Users/fengqian/anaconda2/bin/python /Users/fengqian/Downloads/simulated_seqs_1_for_indels.py /Users/fengqian/Downloads/UniMelb_shared-master/project/mosaic_data/Protein_translateable_pilot_upper_centroids.fasta /Users/fengqian/simulated_tree.txt /Users/fengqian/control.txt
#######################################################################
import sys, os
import glob
from collections import defaultdict

input_tree_txt=sys.argv[1]
input_control=sys.argv[2]

with open(input_tree_txt, 'r') as myfile:
    tree = myfile.read()
with open(input_control, 'r') as myfile:
    myfile=myfile.read()
    new=myfile.replace("(A:0.1,B:0.1);", tree)
with open(input_control, "w") as text_file:
    text_file.write(new)







