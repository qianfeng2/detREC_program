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

# Usuage instruction:
# This code to recombine sequences based on given sequences and recombinnig method. Input files are all from mafft prepared results, output csv files are in the results dir 
# no input pars, output files are in results/result_all.csv
# This is used for python 2. all.py means I use all triples in the final mosaic jump (once, twice and more) chunks.
# Usage example: /Users/fengqian/anaconda2/bin/python /Users/fengqian/Downloads/algorithm_rec_detection_all.py 

from __future__ import division
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import glob
import sys, os
from mungo.fasta import FastaReader
from collections import defaultdict
from scipy import stats
import time
import numpy as np
import random
import time
import pandas as pd
import csv

dir=os.getcwd();dir=dir+"/"

#### This function is used for returning the recombinant sequence by new method
def main_new(fastafile,bkp):
    distance_name=["ab","ac","bc"];
    seq_name=[]    
    for h,s in FastaReader(fastafile):
        seq_name.append(h)
    aln = AlignIO.read(open(fastafile), 'fasta')
    calculator = DistanceCalculator('blosum62')
    segment_1 = calculator.get_distance(aln[:, :bkp])
    segment_2 = calculator.get_distance(aln[:, bkp:])
    distance=[segment_1[seq_name[1]][0],segment_1[seq_name[2]][0],segment_1[seq_name[2]][1],segment_2[seq_name[1]][0],segment_2[seq_name[2]][0],segment_2[seq_name[2]][1]];
    #distance=[segment_1[seq_name[1]][0],segment_1[seq_name[2]][0],segment_1[seq_name[2]][1],segment_2[seq_name[1]][0],segment_2[seq_name[2]][0],segment_2[seq_name[2]][1]];
    compare_distance=[abs(distance[0]-distance[3]),abs(distance[1]-distance[4]),abs(distance[2]-distance[5])]##in order of ab,ac,bc
    temp2 = distance_name[compare_distance.index(min(compare_distance))]
    string = "abc";string = string.replace(temp2[0],"");string = string.replace(temp2[1],"")
    rec=seq_name["abc".index(string)]   
    return rec


input=dir+"complement_chunks"
output=dir+"results/result_all.csv"
headers=["chunk","target","db1","db2","rec","sv"]
with open(output, 'a+') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(headers)
for fasta_file in glob.glob(input+"/*.fasta"):
    bkp=int(fasta_file.split(".fasta")[0].split("_")[-1])
    chunk=fasta_file.split("/chunk")[1].split("_")[0]
    initial= main_new(fasta_file,bkp)
    #print initial
    mapping_dict={};seq_name=[] 
    pvalue=0;permutation=100;
    for h,s in FastaReader(fasta_file):
        mapping_dict[h]=s;seq_name.append(h);full_alignment_length=len(s)
    for m in range(permutation):
        index_1 = np.random.choice(bkp, bkp, replace=True);shuffled_sequences_1=[]
        for n in [mapping_dict[seq_name[0]][0:bkp],mapping_dict[seq_name[1]][0:bkp],mapping_dict[seq_name[2]][0:bkp]]:
            temp_list_1="";
            for j in index_1:
                temp_list_1 = temp_list_1+ n[j];
            shuffled_sequences_1.append(temp_list_1)
        s2_length= full_alignment_length-bkp
        index_2 = np.random.choice(s2_length, s2_length, replace=True);shuffled_sequences_2=[]
        for n in [mapping_dict[seq_name[0]][bkp:],mapping_dict[seq_name[1]][bkp:],mapping_dict[seq_name[2]][bkp:]]:
            temp_list_2="";
            for j in index_2:
                temp_list_2 = temp_list_2+ n[j];
            shuffled_sequences_2.append(temp_list_2)
        shuffled_sequences=[x+y for x,y in zip(shuffled_sequences_1, shuffled_sequences_2)]
        with open(dir+"temp/shuffled_sequence"+str(m)+".fasta", 'w') as outfile:
            outfile.write(">"+seq_name[0]+"\n"+shuffled_sequences[0]+"\n"+">"+seq_name[1]+"\n"+shuffled_sequences[1]+"\n"+">"+seq_name[2]+"\n"+shuffled_sequences[2]+"\n") 
        if main_new(dir+"temp/shuffled_sequence"+str(m)+".fasta",bkp)==initial:
            pvalue+=1/permutation
        os.remove(dir+"temp/shuffled_sequence"+str(m)+".fasta")
    #rec_onesetting.append(initial);pvalue_onesetting.append(pvalue)
    result=[chunk,seq_name[0],seq_name[1],seq_name[2],initial,pvalue]
    with open(output, 'a+') as outfile:
        outfile.write(",".join([str(l) for l in result]) + "\n")
