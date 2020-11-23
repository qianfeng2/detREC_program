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
# This is used in Python 3.
# This code to recombine sequences based on given sequences and recombinnig strategy. Input file name must be simulated_seqs.fasta, output fasta file is in the same dir, called simulated_seqs_recombined.fasta,#=200 
# Four input pars together: 
# First par is 0/1, 0:jump once to combine seqs,1:intoduce jump twice to combine seqs
# Second par: number of jump once seqs
# Third par: number of jump twice seqs
# Last input is simulated sequences to be recombined.
# Usage example: python /Users/fengqian/Downloads/recombined_seqs.py 1 50 50 /Users/fengqian/Downloads/simulated_seqs.fasta
#######################################################################

from Bio import SeqIO
import sys, os
import random
input_par1=int(sys.argv[1])
input_par2=int(sys.argv[2])
input_par3=int(sys.argv[3])
input_fasta=sys.argv[4]

#input_par1=1
#input_par2=15
#input_par3=15
#input_fasta=dir+"data/simulated_seqs.fasta"

dir=os.getcwd();dir=dir+"/"


seqs = {};count=0;
for seq_record in SeqIO.parse(input_fasta, "fasta"):
    count+=1
    seqs[seq_record.id] = str(seq_record.seq)




if input_par1==0:## eg:set 100 jump once seqs and then remaining 100 is nonrecombinants.
    need_seqs_num=2*input_par2+200-input_par2
    seq_index=random.sample(range(need_seqs_num), k=2*input_par2)
    for count in range(input_par2):
        left_source=seq_index[2*count];right_source=seq_index[2*count+1]
        with open(dir+"data/temp.fasta", 'w') as tempfile:
            tempfile.write(">"+"seq"+ str(left_source)+"\n"+seqs["seq" + str(left_source)]+"\n"+">"+"seq"+str(right_source)+"\n"+seqs["seq" + str(right_source)]+"\n")
        cmd="mafft --auto  --inputorder data/temp.fasta > data/temp_align.fasta"
        os.system(cmd)
        records = list(SeqIO.parse("data/temp_align.fasta", "fasta"))
        aligned_len=len(str(records[1].seq))
        bkp=random.sample(range(1,aligned_len),1)[0]
        left_source_seq=str(records[0].seq)[0:bkp].replace("-","");right_source_seq=str(records[1].seq)[bkp:].replace("-","");  
        bkp=len(left_source_seq)
        with open(dir+"data/simulated_seqs_recombined.fasta", 'a+') as outfile:
            outfile.write(">"+"r_"+str(left_source)+"_"+str(right_source)+"_"+str(bkp)+"\n"+left_source_seq+right_source_seq+"\n")
    with open(dir+"data/simulated_seqs_recombined.fasta", 'a+') as outfile:                
        seq_index_unrecombined=list(set(range(need_seqs_num)) - set(seq_index))
        for count in seq_index_unrecombined:
            outfile.write(">"+"seq" + str(count)+"\n"+seqs["seq" + str(count)]+"\n")   
else:
####OR eg:set 50 jump once seqs and 50 jump twice seqs, remaining 100 is nonrecombinants. 
    need_seqs_num=input_par2*2+input_par3*3+200-input_par2-input_par3
    seq_index=random.sample(range(need_seqs_num), k=input_par2*2)
    seq_index_2=random.sample(list(set(range(need_seqs_num)) - set(seq_index)), k=input_par3*3)
    
    for count in range(input_par2):
        left_source=seq_index[2*count];right_source=seq_index[2*count+1]
        with open(input_fasta[:-20]+"temp.fasta", 'w') as tempfile:
            tempfile.write(">"+"seq"+ str(left_source)+"\n"+seqs["seq" + str(left_source)]+"\n"+">"+"seq"+str(right_source)+"\n"+seqs["seq" + str(right_source)]+"\n")
        cmd="mafft --auto  --inputorder data/temp.fasta > data/temp_align.fasta"
        os.system(cmd)
        records = list(SeqIO.parse("data/temp_align.fasta", "fasta"))
        aligned_len=len(str(records[1].seq))
        bkp=random.sample(range(1,aligned_len),1)[0]
        left_source_seq=str(records[0].seq)[0:bkp].replace("-","");right_source_seq=str(records[1].seq)[bkp:].replace("-","");  
        bkp=len(left_source_seq)
        with open(dir+"data/simulated_seqs_recombined.fasta", 'a+') as outfile:
            outfile.write(">"+"r_"+str(left_source)+"_"+str(right_source)+"_"+str(bkp)+"\n"+left_source_seq+right_source_seq+"\n")
    
    for count in range(input_par3):
        left_source=seq_index_2[3*count];middle_source=seq_index_2[3*count+1];right_source=seq_index_2[3*count+2]
        with open(dir+"data/temp.fasta", 'w') as tempfile:
            tempfile.write(">"+"seq"+ str(left_source)+"\n"+seqs["seq" + str(left_source)]+"\n"+">"+"seq"+ str(middle_source)+"\n"+seqs["seq" + str(middle_source)]+"\n"+">"+"seq"+str(right_source)+"\n"+seqs["seq" + str(right_source)]+"\n")
        cmd="mafft --auto  --inputorder data/temp.fasta > data/temp_align.fasta"
        os.system(cmd)
        records = list(SeqIO.parse("data/temp_align.fasta", "fasta"))
        aligned_len=len(str(records[1].seq))
        bkp=sorted(random.sample(range(1,aligned_len),2))
        bkp1=bkp[0];bkp2=bkp[1]
        left_source_seq=str(records[0].seq)[0:bkp1].replace("-","");middle_source_seq=str(records[1].seq)[bkp1:bkp2].replace("-",""); right_source_seq=str(records[2].seq)[bkp2:].replace("-","");
        with open(dir+"data/simulated_seqs_recombined.fasta", 'a+') as outfile:
            outfile.write(">"+"r_"+str(left_source)+"_"+str(middle_source)+"_"+str(right_source)+"\n"+left_source_seq+middle_source_seq+right_source_seq+"\n")
        
    
    with open(dir+"data/simulated_seqs_recombined.fasta", 'a+') as outfile:  
        seq_index_unrecombined=list(set(range(need_seqs_num)) - set(seq_index)-set(seq_index_2))
        for count in seq_index_unrecombined:
            outfile.write(">"+"seq" + str(count)+"\n"+seqs["seq" + str(count)]+"\n")   

os.system("rm data/temp*.fasta")




