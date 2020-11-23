
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

# Usuage instruction
# This code is to collect all bkp errors from 100 replicate in each specific setting.
# Input pars are two: (1) input_dir, (2) output csv file. 
# This is used in Python 3. 
# Usage example: /Users/fengqian/anaconda2/envs/py3/bin/python codes/bkp_error_collect_from_replicates.py /data/cephfs/punim0609/qian_feng/snake_pipeline/simulation_50-50/ /data/cephfs/punim0609/qian_feng/snake_pipeline/final_results/jump_once_only/simulation_50-50_output_error.csv

#######################################################################


from collections import Counter
import glob
import csv
import sys
import pandas as pd

input_dir=sys.argv[1]
output=sys.argv[2]
#input_dir= "/data/cephfs/punim0609/qian_feng/snake_pipeline/simulation_50-50/"
#output = "/data/cephfs/punim0609/qian_feng/snake_pipeline/final_results/jump_once_only/simulation_50-50_output_error.csv"
for jump in range(100):
	result_all = pd.read_csv(input_dir+"replicate_"+str(jump)+"/results/result_all_TP_bkp.csv")
	index=[jump] * len(result_all)
	result_all['replicate'] = index
	result_all.to_csv(output, header=False, mode='a',index=False)	










