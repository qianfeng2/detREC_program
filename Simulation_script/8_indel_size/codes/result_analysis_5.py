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
# This code(used in Python 3) is to collect sv values for TP and FP results into two csv files for each simulation setting. 
# Usage example: /Users/fengqian/anaconda2/envs/py3/bin/python /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/codes/result_analysis_2.py /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/simulation_50-50 /Users/fengqian/mosaic_simulation_output2.csv
#                python /data/cephfs/punim0609/qian_feng/snake_pipeline/codes/result_analysis_2_all.py /data/cephfs/punim0609/qian_feng/snake_pipeline/simulation_50-50 /data/cephfs/punim0609/qian_feng/snake_pipeline/mosaic_simulation_output2.csv
#######################################################################

from collections import Counter
import glob
import sys
import csv
import numpy as np

inputdir=sys.argv[1]
output_TP=sys.argv[2]
output_FP=sys.argv[3]

for jump in range(100):
    dir=inputdir+"/replicate_"+str(jump)+"/"
    rec = [];sv=[];
    with open(dir+'results/result_all.csv', 'r') as f:
        reader = csv.DictReader(f, delimiter=',')  
        for n, row in enumerate(reader):
            rec.append(row['rec']) 
            sv.append(float(row['sv'])) 
    
    sv_nonrec_index=[];sv_rec_index=[]
    for i in range(len(rec)):
        if all(c in rec[i] for c in "seq"):
            sv_nonrec_index.append(i)
        else:sv_rec_index.append(i)
    sv_nonrec=[sv[c] for c in sv_nonrec_index]## "seq" in our identified recombination sv result 
    sv_rec=[sv[c] for c in sv_rec_index]## "r_" in our identified recombination sv result (true rec)
    
    with open(output_TP, "a") as f:
        writer = csv.writer(f)
        for val in sv_rec:
            writer.writerow([val])

    with open(output_FP, "a") as f:
        writer = csv.writer(f)
        for val in sv_nonrec:
            writer.writerow([val])

