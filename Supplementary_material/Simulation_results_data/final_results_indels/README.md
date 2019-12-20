Readme.txt file for Qian Feng created 12/08/2019
 
Indel process involve two key things:
(1) indel rate
(2) size distribution
First I change indel rate (0,0.04,by=0.01), size distribution is fixed NB(q=0.2,r=4),files are named with indels_0,indels_0.01,indels_0.02,indels_0.03,indels_0.04;

Second I change size distribution mean and fix variance by changing q(0.1,0.5,by=0.1) in NB distribution without changing indel rate(0.02).Files are named with indels_q_0.1,indels_q_0.2,indels_q_0.3,indels_q_0.4,indels_q_0.5. 
Note that indels_q_0.2==indel_0.02

One final thing is WAG, root tree protein length is 200, sequence involver is INDELible, recombination proportion 0.5
