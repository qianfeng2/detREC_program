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
# This code(used in Python 3) is to summarize final simulated rec det results into once csv file. 
# Input files are result_jump_once and result_jump_twice.csv, data/simulated_seqs_recombined.fasta
# Output csv file shows All the binary classification metric: accuracy_score,Precision,Sensitivity(TPR),Specifity,f1_score()
# Confusion matrix example: rowsum is true 50-50; colsum is prediction result.
# array([[47,  3],  (TN, FP
#      [ 8, 42]]).   FN, TP)
# Formulas for all metrics are here: https://towardsdatascience.com/hackcvilleds-4636c6c1ba53
# all.py means I use all triples in the final mosaic jump (once, twice and more) chunks.
# Usage example: /Users/fengqian/anaconda2/envs/py3/bin/python /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/codes/result_analysis_2.py /Users/fengqian/Downloads/UniMelb_shared-master/algorithm_simulation/protein/snake_pipeline/simulation_50-50 /Users/fengqian/mosaic_simulation_output2.csv
#                python /data/cephfs/punim0609/qian_feng/snake_pipeline/codes/result_analysis_2_all.py /data/cephfs/punim0609/qian_feng/snake_pipeline/simulation_50-50 /data/cephfs/punim0609/qian_feng/snake_pipeline/mosaic_simulation_output2.csv 50
# Last input is total number of sequences in each fasta file, for the purpose of simulating the changed dataset size case.
#######################################################################

from collections import Counter
import glob
import sys
import csv
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import confusion_matrix 
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from statistics import mean 
import numpy as np

inputdir=sys.argv[1]
output=sys.argv[2]
input_par=int(sys.argv[3])

for jump in range(100):
    dir=inputdir+"/replicate_"+str(jump)+"/"
    rec = [];sv=[];seqs=[]
    with open(dir+'results/result_all.csv', 'r') as f:
        reader = csv.DictReader(f, delimiter=',')  
        for n, row in enumerate(reader):
            rec.append(row['rec']) 
            sv.append(float(row['sv'])) 
    with open(dir+"data/simulated_seqs_recombined.fasta", 'r') as infile:
        for line in infile:
            if ">" in line:
                seqs.append(line.strip().split(">")[1])
    y_true=[]## 0 is nonrec; 1 is rec
    for i in range(input_par):
        if all(c in seqs[i] for c in "seq"):
            y_true.append(0)
        else:y_true.append(1)
    y_pred=[]
    for i in range(input_par):
        if seqs[i] in rec:
            y_pred.append(1)
        else:
            y_pred.append(0)     
    temp=confusion_matrix(y_true, y_pred);Precision=temp[1,1]/(temp[1,1]+temp[0,1]);Sensitivity=temp[1,1]/(temp[1,1]+temp[1,0]);Specifity=temp[0,0]/(temp[0,0]+temp[0,1])

    sv_nonrec_index=[];sv_rec_index=[]
    for i in range(len(rec)):
        if all(c in rec[i] for c in "seq"):
            sv_nonrec_index.append(i)
        else:sv_rec_index.append(i)
    sv_nonrec=mean([sv[c] for c in sv_nonrec_index]) if len(sv_nonrec_index)!=0 else 0## "seq" in our identified recombination sv result 
    sv_rec=mean([sv[c] for c in sv_rec_index]) if len(sv_rec_index)!=0 else 0## "r_" in our identified recombination sv result (true rec)
    y_scores=[0] * input_par
    for i in range(input_par):
        if seqs[i] in rec:
            temp=[j for j in range(len(rec)) if rec[j]==seqs[i]]
            y_scores[i]=mean([sv[c] for c in temp]) if len(temp)!=1 else sv[temp[0]]    
    fpr, tpr, _ = roc_curve(y_true, y_scores);roc_auc = auc(fpr, tpr)
    result=[jump,accuracy_score(y_true, y_pred),Precision,Sensitivity,Specifity,f1_score(y_true, y_pred, average='binary'),sv_nonrec,sv_rec,roc_auc]

    with open(output, "a") as f:
        writer = csv.writer(f)
        writer.writerow(result)

