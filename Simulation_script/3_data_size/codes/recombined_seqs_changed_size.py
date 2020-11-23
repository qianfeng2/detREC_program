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
# This code to recombine sequences based on given sequences and recombinnig method. Input file name must be simulated_seqs.fasta, output fasta file is in the same dir, called simulated_seqs_recombined.fasta.
# output fasta file have different number of sequences, for the purpose of simulating changed dataset size senario. 
# Five input pars together: 
# First par is 0/1, 0:jump once to combine seqs,1:intoduce jump twice to combine seqs
# Second par: number of jump once seqs
# Third par: number of jump twice seqs
# Fourth par: total number of sequences in output file.50,150,200,250,or more
# Last input is input sequences to be recombined.
# Usage example: /Users/fengqian/anaconda2/bin/python /Users/fengqian/Downloads/recombined_seqs_changed_size.py 0 25 0 50 /Users/fengqian/Downloads/simulated_seqs.fasta
#######################################################################

from mungo.fasta import FastaReader
import sys, os
import random
input_par1=int(sys.argv[1])
input_par2=int(sys.argv[2])
input_par3=int(sys.argv[3])
input_par4=int(sys.argv[4])
input_fasta=sys.argv[5]



 
seqs = {};count=0;seq_length=0;temp_length=[]
for h,s in FastaReader(input_fasta):
    seqs[h] = s
    count+=1
    temp_length.append(len(s))

seq_length=min(temp_length)## length should always be the same without considering indel events.

if input_par1==0:## eg:set 25 jump once seqs and 25 is nonrecombinants, then input_par4==50
    need_seqs_num=2*input_par2+input_par4-input_par2
    seq_index=random.sample(range(need_seqs_num), k=2*input_par2)
    with open(input_fasta[:-20]+"simulated_seqs_recombined.fasta", 'w') as outfile:
        for count in range(input_par2):
            left_source=seq_index[2*count];right_source=seq_index[2*count+1]
            temp=min(len(seqs["seq" + str(left_source)]),len(seqs["seq" + str(right_source)]))
            bkp=random.sample(range(1,temp),1)[0]
            left_source_seq=seqs["seq" + str(left_source)][0:bkp];right_source_seq=seqs["seq" + str(right_source)][bkp:];        
            outfile.write(">"+"r_"+str(left_source)+"_"+str(right_source)+"_"+str(bkp)+"\n"+left_source_seq+right_source_seq+"\n")
        seq_index_unrecombined=list(set(range(need_seqs_num)) - set(seq_index))
        for count in seq_index_unrecombined:
            outfile.write(">"+"seq" + str(count)+"\n"+seqs["seq" + str(count)]+"\n")   
else:
####OR eg:set 0 jump once seqs and 25 jump twice seqs. Remaining 25 is nonrecombinants. 
    need_seqs_num=input_par2*2+input_par3*3+input_par4-input_par2-input_par3
    seq_index=random.sample(range(need_seqs_num), k=input_par2*2)
    seq_index_2=random.sample(list(set(range(need_seqs_num)) - set(seq_index)), k=input_par3*3)
    
    with open(input_fasta[:-20]+"simulated_seqs_recombined.fasta", 'w') as outfile:
        for count in range(input_par2):
            left_source=seq_index[2*count];right_source=seq_index[2*count+1]
            temp=min(len(seqs["seq" + str(left_source)]),len(seqs["seq" + str(right_source)]))
            bkp=random.sample(range(1,temp),1)[0]
            left_source_seq=seqs["seq" + str(left_source)][0:bkp];right_source_seq=seqs["seq" + str(right_source)][bkp:];        
            outfile.write(">"+"r_"+str(left_source)+"_"+str(right_source)+"_"+str(bkp)+"\n"+left_source_seq+right_source_seq+"\n")
        for count in range(input_par3):
            left_source=seq_index_2[3*count];middle_source=seq_index_2[3*count+1];right_source=seq_index_2[3*count+2]
            temp=min(len(seqs["seq" + str(left_source)]),len(seqs["seq" + str(middle_source)]),len(seqs["seq" + str(right_source)]))
            bkp=sorted(random.sample(range(1,temp),2))
            bkp1=bkp[0];bkp2=bkp[1]
            left_source_seq=seqs["seq" + str(left_source)][0:bkp1];middle_source_seq=seqs["seq" + str(middle_source)][bkp1:bkp2]; right_source_seq=seqs["seq" + str(right_source)][bkp2:];        
            outfile.write(">"+"r_"+str(left_source)+"_"+str(middle_source)+"_"+str(right_source)+"\n"+left_source_seq+middle_source_seq+right_source_seq+"\n")
    
        seq_index_unrecombined=list(set(range(need_seqs_num)) - set(seq_index)-set(seq_index_2))
        for count in seq_index_unrecombined:
            outfile.write(">"+"seq" + str(count)+"\n"+seqs["seq" + str(count)]+"\n")   







# In[ ]:




