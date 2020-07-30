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
# This code to generate a list of recombinants according to proposed novel algorithm from mosaic representation results. 
# Note that this code use all triples in jump # chunks(ont only jump once and twice chunks). 
# Input are three files, 
## The first one is partial alignment from mosaic, eg.its name is ***_align.txt; 
## Second file is the original .fasta file
## Third file is the output csv file showing recombinants, its name is ***.csv
# This is used in Python >= 3.5  
# Usage example: python /Users/fengqian/Desktop/replicate_76/codes/integrated_rec_det.py /Users/fengqian/Desktop/replicate_76/simulated_seqs_recombined_align.txt /Users/fengqian/Desktop/replicate_76/integrated_rec_det_results.csv
#######################################################################



from collections import defaultdict
from collections import Counter
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import SeqIO
from scipy import stats
import numpy as np
import pandas as pd
import glob
import random
import os
import sys
import csv



input=sys.argv[1];input_fasta=sys.argv[2]
dir=os.getcwd();dir=dir+"/"
#input="/Users/fengqian/Desktop/replicate_20/simulated_seqs_recombined_align.txt"
#dir="/Users/fengqian/Desktop/replicate_20/"

# section1: mafft_prepared_all_real.py
### Step1: make new file folders.
path = dir+"complement_chunks";os.makedirs(path)
with open(input) as infile:
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
    line_index=[*[startline[i]],*range(startline[i]+2,endline[i]+1)]
    with open(input) as infile:
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
                with open(dir+"temp/chunk_"+str(i)+"/original_chunk_db.txt") as db_infile:
                    for n, db_line in enumerate(db_infile.readlines()):
                        line1=db_line.strip().split();h1=line1[0];s1=line1[1]
                        if n==0: outfile.write(">"+h+"\n"+s[:bkp[k]]+"\n"+">"+h1+"\n"+s1+"\n")
            else: 
                with open(dir+"temp/chunk_"+str(i)+"/original_chunk_db.txt") as db_infile:
                    for n, db_line in enumerate(db_infile.readlines()):
                        line1=db_line.strip().split();h1=line1[0];s1=line1[1]
                        if n==k: outfile.write(">"+h+"\n"+s[bkp[k-1]:bkp[k]]+"\n"+">"+h1+"\n"+s1+"\n")





### Step3: add the part for mafft alignment. use replace not split.

mapping_dict_seq_identifier = {}
test = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
for k,v in test.items():
    mapping_dict_seq_identifier[k]=str(v.seq)
for i in range(Chunk_count):
    db_line_index=range(startline[i]+2,endline[i]+1);line_index=[*[startline[i]],*range(startline[i]+2,endline[i]+1)] 
    for k in range(mosaic_output_dbcount[i]-1):     
        with open(input) as infile:
            for j, line in enumerate(infile.readlines()):                 
                if j == db_line_index[k]:
                    line1=line.strip().split()
                    h1=line1[0];s1=line1[1]
                    if len(s1)>=10:
                        original_seq=mapping_dict_seq_identifier[h1];s=line1[1].replace("-", "");                      
                        with open(dir+"temp/chunk_"+str(i)+"/s"+str((k+2))+"/add_up.fasta","w") as outfile:
                            if original_seq.find(s)!= -1 and original_seq.split(s)[1]!="": outfile.write(">"+h1+"\n"+original_seq.split(s)[1]+"\n")
    for m in range(1,mosaic_output_dbcount[i]):     
        with open(input) as infile:
            for j, line in enumerate(infile.readlines()):                 
                if j == db_line_index[m]:
                    line1=line.strip().split()
                    h1=line1[0];s1=line1[1]
                    if len(s1)>=10:### only investigate the segment which longer than 10 to get reliable results.
                        original_seq=mapping_dict_seq_identifier[h1];s=line1[1].replace("-", "");
                        with open(dir+"temp/chunk_"+str(i)+"/s"+str((m))+"/add_down.fasta","w") as outfile:
                            if original_seq.find(s)!= -1 and original_seq.split(s)[0]!="": outfile.write(">"+h1+"\n"+original_seq.split(s)[0]+"\n")








# section2: mafft_concat_all_real.py

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
            cmd="seqkit concat "+filename1+" "+filename2+" > "+"complement_chunks/chunk"+str(i)+"_s"+str((k+1))+str((k+2))+"_"+str(bkp)+".fasta"
            os.system(cmd)




# section3: algorithm_rec_detection_all_real.py


#### This function is used for returning the recombinant sequence by new method
def main_new(fastafile,bkp):
    distance_name=["ab","ac","bc"];
    temp = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))    
    seq_name = [*temp]
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
output=sys.argv[3]###csv file which is used for store the recombinant in each run

for fasta_file in glob.glob(input+"/*.fasta"):
    bkp=int(fasta_file.split(".fasta")[0].split("_")[-1])
    chunk=fasta_file.split("/chunk")[1].split("_")[0]
    initial= main_new(fasta_file,bkp)
    #print initial
    mapping_dict={};seq_name=[] 
    pvalue=0;permutation=100;
    test = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    seq_name = [*test]
    for k,v in test.items():
        mapping_dict[k]=str(v.seq)
    full_alignment_length=len(mapping_dict[k])
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
    result=[chunk,seq_name[0],seq_name[1],seq_name[2],initial,"{:.2f}".format(pvalue)]#headers=["chunk","target","db1","db2","rec","sv"]
    with open(output, 'a+') as outfile:
        outfile.write(",".join([str(l) for l in result]) + "\n")

