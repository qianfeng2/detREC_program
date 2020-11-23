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
# This code to generate a txt file : a list of bkp errors(identified bkp - actual bkp) txt for each replicate 
# No Input pars for this file, but you should get the identified recombinant list in the results/* folder before running this script. 
# Output is in results/result_bkperror.txt
# This is used in Python 3
# Usage example: python /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/codes/bkp_error.py 
#######################################################################

from __future__ import division
import glob
import sys, os
from collections import defaultdict
import numpy as np
import random
import pandas as pd
import csv
from Bio import SeqIO




dir=os.getcwd();dir=dir+"/"
result_all = pd.read_csv(dir+"results/result_all.csv") 
input=dir+"simulated_seqs_recombined_align.txt"
input_2=dir+"data/simulated_seqs_recombined.fasta"
bkp_error = [999999999] * len(result_all)
identified_bkp = [999999999] * len(result_all)
actual_bkp = [999999999] * len(result_all)
result_all['bkp_error'] = bkp_error
result_all['identified_bkp'] = identified_bkp
result_all['actual_bkp'] = actual_bkp
with open(input, 'r') as infile:
    unique_line_index = [];Par_line_index = 0;Chunk_count = 0
    for i, line in enumerate(infile.readlines()):
        line = line.strip().split()
        if "Parameters" in line:#extract line index for lines that contain nothing
            Par_line_index = i
        elif "Target:" in line:
            unique_line_index.append(i)
            Chunk_count = Chunk_count +1
    startline = [m+1 for m in unique_line_index]
    endline = [m-3 for m in unique_line_index if m != unique_line_index[0]]+[Par_line_index-2]
mosaic_output_dbcount=[endline[m]-startline[m]-1 for m in range(Chunk_count)]





full_length_sequence = SeqIO.to_dict(SeqIO.parse(input_2, "fasta"))
for i in range(Chunk_count):
    for k in range(mosaic_output_dbcount[i]-1):
        filename1 = "temp/chunk_"+str(i)+"/s"+str((k+1))+"_mafft_down.fasta"
        filename2 = "temp/chunk_"+str(i)+"/s"+str((k+2))+"_mafft_up.fasta"        
        if os.path.isfile(filename1) and os.path.isfile(filename2): 
            records_bkp = list(SeqIO.parse(filename1, "fasta"))
            bkp=len(str(records_bkp[1].seq))
            records = list(SeqIO.parse("complement_chunks/chunk"+str(i)+"_s"+str((k+1))+str((k+2))+"_"+str(bkp)+".fasta", "fasta"))  
            #print records[0].description
            temp=result_all[(result_all["chunk"]==i) &  (result_all["target"]==records[0].description)&  (result_all["db1"]==records[1].description)&  (result_all["db2"]==records[2].description)]
            rec= [str(x) for x in temp["rec"]][-1]## get identified rec seq           
            #print records[0].description
            if "seq" in rec:continue## ignore the FP cases
            actual_bkp=  int([str(x) for x in temp["rec"]][-1].split("_")[-1])## get actual bkp for tp recombinants
            if records[2].description==rec:
                records_up = SeqIO.to_dict(SeqIO.parse(filename2,"fasta"))               
                seg = str(records_up[rec].seq).replace("-", "")
                original_sequence = str(full_length_sequence[rec].seq)
                identified_bkp=original_sequence.find(seg)
            if records[0].description==rec or records[1].description==rec:
                records_down = SeqIO.to_dict(SeqIO.parse(filename1,"fasta"))
                seg = str(records_down[rec].seq).replace("-", "")
                original_sequence = str(full_length_sequence[rec].seq)
                identified_bkp=original_sequence.find(seg)+bkp-str(records_down[rec].seq).count("-")  
            row_index= result_all[(result_all["chunk"]==i) &  (result_all["target"]==records[0].description)&  (result_all["db1"]==records[1].description)&  (result_all["db2"]==records[2].description)].index.tolist()
            result_all.loc[row_index[0],"bkp_error"] = identified_bkp-actual_bkp
            result_all.loc[row_index[0],"identified_bkp"] = identified_bkp
            result_all.loc[row_index[0],"actual_bkp"] = actual_bkp 

result_all_TP= result_all.loc[result_all['rec'].str.contains('seq')==False,:]
#bkp_error_final = result_all_TP["bkp_error"].tolist()# I include all bkp error for all TPs, even though some TP occurrence bigger than 1
#result_all_TP.loc[result_all_TP.groupby('rec')['sv'].idxmax()]
result_all_TP.to_csv(dir+"results/result_all_TP_bkp.csv",index=False)








