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
# This code to generate small segments based on mosaic align.txt and run mafft and seq concat. 
# Note that this code use all triples in jump # chunks(ont only jump once and twice chunks). 
# Input file name is ***_align.txt, output fasta file is in the complement_chunks/ dir,#=100 
# This is used in Python 2
# Usage example: /Users/fengqian/anaconda2/bin/python /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/codes/mafft_prepared_all.py /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/20190614_newprotein_align.txt
#######################################################################


from __future__ import division
import glob
from mungo.fasta import FastaReader
from collections import defaultdict
import numpy as np
import random
from collections import Counter
import os
import sys



input=sys.argv[1]
dir=os.getcwd();dir=dir+"/"
#input="/Users/fengqian/Desktop/replicate_20/simulated_seqs_recombined_align.txt"
#dir="/Users/fengqian/Desktop/replicate_20/"


### Step1: make new file folders.
path = dir+"complement_chunks";os.makedirs(path)
with open(input, 'rU') as infile:
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
for i in range(Chunk_count):
    path = dir+"temp/chunk_"+str(i);os.makedirs(path)
    for j in range(mosaic_output_dbcount[i]):
        path_detail = dir+"temp/chunk_"+str(i)+"/s"+str((j+1));os.makedirs(path_detail)        




### Step2: extract each chunk's original segment pairs. 
for i in range(Chunk_count):
    len_db=[];bkp=len_db;db_line_index=range(startline[i]+2,endline[i]+1)##note bkp include the last position of target seq
    line_index=[startline[i]]+range(startline[i]+2,endline[i]+1)
    with open(input, 'rU') as infile:
        for j, line in enumerate(infile.readlines()):           
            if j==startline[i]:
                line_target=line.strip().split();s=line_target[1];h=line_target[0]
                with open(dir+"temp/chunk_"+str(i)+"/original_chunk.txt","a+") as outfile:
                    outfile.write(line) 
            if j in db_line_index:
                line1=line.strip().split()
                h1=line1[0];s1=line1[1]
                len_db.append(len(s1))
                with open(dir+"temp/chunk_"+str(i)+"/original_chunk_db.txt","a+") as outfile:
                    outfile.write(line)
                with open(dir+"temp/chunk_"+str(i)+"/original_chunk.txt","a+") as outfile:
                    outfile.write(line)
    for m in range(1,len(len_db)):bkp[m]=bkp[(m-1)]+len_db[m]
    for k in range(mosaic_output_dbcount[i]):            
        with open(dir+"temp/chunk_"+str(i)+"/s"+str((k+1))+"/original.fasta","a+") as outfile:
            if k==0:
                with open(dir+"temp/chunk_"+str(i)+"/original_chunk_db.txt","rU") as db_infile:
                    for n, db_line in enumerate(db_infile.readlines()):
                        line1=db_line.strip().split();h1=line1[0];s1=line1[1]
                        if n==0: outfile.write(">"+h+"\n"+s[:bkp[k]]+"\n"+">"+h1+"\n"+s1+"\n")
            else: 
                with open(dir+"temp/chunk_"+str(i)+"/original_chunk_db.txt","rU") as db_infile:
                    for n, db_line in enumerate(db_infile.readlines()):
                        line1=db_line.strip().split();h1=line1[0];s1=line1[1]
                        if n==k: outfile.write(">"+h+"\n"+s[bkp[k-1]:bkp[k]]+"\n"+">"+h1+"\n"+s1+"\n")





### Step3: add the part for mafft alignment. use replace not split.
mapping_dict_seq_identifier = {}
for h,s in FastaReader(dir+"data/simulated_seqs_recombined.fasta"):
    mapping_dict_seq_identifier[h]=s
for i in range(Chunk_count):
    db_line_index=range(startline[i]+2,endline[i]+1);line_index=[startline[i]]+range(startline[i]+2,endline[i]+1)
    for k in range(mosaic_output_dbcount[i]-1):     
        with open(input, 'rU') as infile:
            for j, line in enumerate(infile.readlines()):                 
                if j == db_line_index[k]:
                    line1=line.strip().split()
                    h1=line1[0];s1=line1[1]
                    if len(s1)>=10:
                        original_seq=mapping_dict_seq_identifier[h1];s=line1[1].replace("-", "");                        
                        with open(dir+"temp/chunk_"+str(i)+"/s"+str((k+2))+"/add_up.fasta","w") as outfile:
                            if original_seq.find(s)!= -1 and original_seq.split(s)[1]!="": outfile.write(">"+h1+"\n"+original_seq.split(s)[1]+"\n")
    for m in range(1,mosaic_output_dbcount[i]):     
        with open(input, 'rU') as infile:
            for j, line in enumerate(infile.readlines()):                 
                if j == db_line_index[m]:
                    line1=line.strip().split()
                    h1=line1[0];s1=line1[1]
                    if len(s1)>=10:### only investigate the segment which longer than 10 to get reliable results.
                        original_seq=mapping_dict_seq_identifier[h1];s=line1[1].replace("-", "");
                        with open(dir+"temp/chunk_"+str(i)+"/s"+str((m))+"/add_down.fasta","w") as outfile:
                            if original_seq.find(s)!= -1 and original_seq.split(s)[-2]!="": outfile.write(">"+h1+"\n"+original_seq.split(s)[-2]+"\n")








