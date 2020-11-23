
#######################################################################
# Copyright (C) 2020  Qian Feng

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

# Usuage instruction
# This code to generate jump statistics based on mosaic align result, and for the purpose of simulating changed dataset size case. 
# Input file name is ***_align.txt, output txt file is in the same dir.
# Last par: total number of sequences in output file. 
# This is used in Python 3.7.3
# Usage example: /Users/fengqian/anaconda2/envs/py3/bin/python /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/codes/recombined_seqs.py /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/20190614_newprotein_align.txt /Users/fengqian/mosaic_analysis_output.txt 
#######################################################################

import sys, os
import glob
from collections import Counter

input=sys.argv[1]
output=sys.argv[2]
input_par2=int(sys.argv[3])
with open(input, 'r') as infile:
    unique_line_index = [];Par_line_index = 0
    for i, line in enumerate(infile.readlines()):
        line = line.strip().split()
        if "Parameters" in line:#extract line index for lines that contain nothing
            Par_line_index = i
        elif "Target:" in line:
            unique_line_index.append(i) 
    startline = [m+1 for m in unique_line_index]
    #print "startline",startline
    endline = [m-3 for m in unique_line_index if m != unique_line_index[0]]+[Par_line_index-2]
    #print "endline",endline

mosaic_output=[endline[m]-startline[m] for m in range(input_par2)]
a=Counter(mosaic_output)
with open(output, 'w') as f:
    for letter in [2,3,4,5,6,7,8,9,10]:
        print ('%s : %d' % (letter, a[letter]),file=f)





