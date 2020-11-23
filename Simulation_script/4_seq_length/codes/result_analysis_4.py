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
# This code is to integrate each simulated dataset's final number of triples into one csv. Generally, csv has 100 rows
# Input file name is result_all.csv from all simulation files, output is one csv file. 
# This is used in Python 3. 
# Usage example: /Users/fengqian/anaconda2/envs/py3/bin/python /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/codes/result_analysis_1.py /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/simulation_50-50 /Users/fengqian/mosaic_simulation_output1.csv
#                python /data/cephfs/punim0609/qian_feng/snake_pipeline/codes/result_analysis_1.py /data/cephfs/punim0609/qian_feng/snake_pipeline/simulation_50-50 /data/cephfs/punim0609/qian_feng/snake_pipeline/mosaic_simulation_output1.csv
#######################################################################


from collections import Counter
import glob
import csv
import sys

inputdir=sys.argv[1]
output=sys.argv[2]
for jump in range(100):
	with open(inputdir+"/replicate_"+str(jump)+"/results/result_all.csv", 'r') as infile:
		count=[];count.append(jump)
		for i, line in enumerate(infile.readlines()):
			line= line
		count.append(i)
	with open(output, "a") as f:
		writer = csv.writer(f)
		writer.writerow(count)









