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

# Usuage instruction(only suitable for python 3+)
# This code to simulate a tree(300 leaves) based on coalescent. output txt file is in the same dir, called simulated_tree.txt
# One and only one input par is the output file name  
# Usage example: /Users/fengqian/anaconda2/envs/py3/bin/python /Users/fengqian/Downloads/msprime_tree.py /Users/fengqian/simulated_tree.txt
#######################################################################


import sys, os
import msprime
import tskit
output_tree_txt=sys.argv[1]
tree_sequence = msprime.simulate(sample_size=300, Ne=0.5)
tree = tree_sequence.first().newick(precision=10)
f = open(output_tree_txt, 'w')
f.write(tree)
f.close()
