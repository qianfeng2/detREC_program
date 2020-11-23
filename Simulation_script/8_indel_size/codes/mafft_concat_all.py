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
# This code to run mafft and seq concat. 
# Note that this code use all triples in jump # chunks(ont only jump once and twice chunks). 
# Input file name is ***_align.txt, output fasta file is in the complement_chunks/ dir,#=100 
# This is used in Python 3
#######################################################################
from __future__ import division
import glob
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import random
from collections import Counter
import os
import sys



input=sys.argv[1]
dir=os.getcwd();dir=dir+"/"
#input="/Users/fengqian/Desktop/replicate_20/simulated_seqs_recombined_align.txt"
#dir="/Users/fengqian/Desktop/replicate_20/"

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





### Step4: let's do mafft alignment
#import subprocess
#subprocess.call(["mafft", "mafft --anysymbol --add /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/add_down.fasta --keeplength --reorder /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/original.fasta > /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/mafft_down.fasta"])
for i in range(Chunk_count):
    for k in range(mosaic_output_dbcount[i]):
        #cmd="mafft --anysymbol --add "+dir+"temp/chunk_"+str(i)+"/s"+str((k+1))+"/add_down.fasta"+" --keeplength --reorder "+dir+"temp/chunk_"+str(i)+"/s"+str((k+1))+"/original.fasta"+" > "+dir+"temp/chunk_"+str(i)+"/s"+str((k+1))+"_mafft_down.fasta"
        cmd="mafft --anysymbol --add "+"temp/chunk_"+str(i)+"/s"+str((k+1))+"/add_down.fasta"+" --keeplength --reorder "+"temp/chunk_"+str(i)+"/s"+str((k+1))+"/original.fasta"+" > "+"temp/chunk_"+str(i)+"/s"+str((k+1))+"_mafft_down.fasta"
        os.system(cmd)        
        cmd="mafft --anysymbol --add "+"temp/chunk_"+str(i)+"/s"+str((k+1))+"/add_up.fasta"+" --keeplength --reorder "+"temp/chunk_"+str(i)+"/s"+str((k+1))+"/original.fasta"+" > "+"temp/chunk_"+str(i)+"/s"+str((k+1))+"_mafft_up.fasta"
        os.system(cmd)
        cmd_del="cd "+"temp/chunk_"+str(i)+" && "+"find . -type f -empty -delete"
        os.system(cmd_del)
        
#os.system('mafft --anysymbol --add /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/add_down.fasta --keeplength --reorder /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/original.fasta > /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/mafft_down.fasta')
#os.system('mafft --anysymbol --add /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/add_up.fasta --keeplength --reorder /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/original.fasta > /Users/fengqian/Desktop/replicate_20/temp/chunk_39/s3/mafft_up.fasta')


### Step5: let's do concatenate for final algorithm
for i in range(Chunk_count):
    for k in range(mosaic_output_dbcount[i]-1):
        filename1 = "temp/chunk_"+str(i)+"/s"+str((k+1))+"_mafft_down.fasta"
        filename2 = "temp/chunk_"+str(i)+"/s"+str((k+2))+"_mafft_up.fasta"        
        if os.path.isfile(filename1) and os.path.isfile(filename2): 
            records_bkp = list(SeqIO.parse(filename1, "fasta"))
            bkp=len(str(records_bkp[1].seq))
            #for h,s in FastaReader(filename1):bkp=len(s)
            cmd="seqkit concat "+filename1+" "+filename2+" > "+"complement_chunks/chunk"+str(i)+"_s"+str((k+1))+str((k+2))+"_"+str(bkp)+".fasta"
            os.system(cmd)




