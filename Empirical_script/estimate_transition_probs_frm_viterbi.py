import sys, os
import argparse
from collections import defaultdict, Counter
from itertools import izip
import math

# Amino acids - IUPAC convention
#       Abbreviation Letter Index Name
#         ---------------------------
#         ???     ? 0      Missing data
#         Ala     A 1      Alanine
#         Arg     R 2      Arginine
#         Asn     N 3      Asparagine
#         Asp     D 4      Aspartic acid (Aspartate)
#         Cys     C 5      Cysteine
#         Gln     Q 6      Glutamine
#         Glu     E 7      Glutamic acid (Glutamate)
#         Gly     G 8      Glycine
#         His     H 9      Histidine
#         Ile     I 10     Isoleucine
#         Leu     L 11     Leucine
#         Lys     K 12     Lysine
#         Met     M 13     Methionine
#         Phe     F 14     Phenylalanine
#         Pro     P 15     Proline
#         Ser     S 16     Serine
#         Thr     T 17     Threonine
#         Trp     W 18     Tryptophan
#         Tyr     Y 19     Tyrosine
#         Val     V 20     Valine
#         Asx     B 21     Aspartic acid or Asparagine
#         Glx     Z 22     Glutamine or Glutamic acid.
#         Xaa     X 23     Any amino acid.
#         TERM    * 24     termination codon

AA_LIST = ['?','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*']

def loadCounts(align_file, verbose):
  alignments = []
  isAlignment=False

  if verbose:
    print "loading file: ", align_file

  with open(align_file, 'rU') as inputfile:
    for line in inputfile:
      if "Target:" in line: #we have got to a new alignment.
        target=""
        seq=""
        isTarget=True
        isAlignment=True
        continue

      if "|" in line:
        isTarget=False

      if "_seq" in line:
        if isTarget:
          target = target + line.strip().split()[1]
        else:
          seq = seq + line.strip().split()[1]

      if not isAlignment: continue #we are not at the end of an alignment yet

      if len(line.strip().split())<1: #we have got to the end of the alignment
        assert len(target)==len(seq), "Lengths of target and seq do not match!"
        alignments.append((target,seq))
        isAlignment=False

  return alignments

def calculate_probabilities(alignments, verbose):

  if verbose:
    print "Calculating propabilities..."

  match_counts = Counter()
  aa_count = Counter()
  match_to_ins_or_del_count = 0
  insert_to_insert_count = 0
  match_to_match = 0
  insert_to_match = 0

  #First add 1 as a pseudocount (same as in zilversmit et al)
  for aa in AA_LIST:
    aa_count[aa]+=1
    match_counts[(aa,aa)] += 1
  for aa in AA_LIST:
    for bb in AA_LIST:
      if aa==bb: continue #Dont double count diagonals
      match_counts[(aa,bb)] += 1

  #Now iterate through alignments iteratively counting the various events
  prev_match=True
  for a in alignments:
    for t,s in izip(a[0],a[1]):
      if t!="-":
        aa_count[t]+=1
      if (t!="-") and (s!="-"):
        match_counts[(t,s)] += 1
        if prev_match:
          match_to_match += 1
        else:
          insert_to_match += 1
        prev_match=True
      elif prev_match:
        match_to_ins_or_del_count += 1
        prev_match=False
      else:
        insert_to_insert_count += 1
        prev_match=False

  #Now calculate the updated transition probabilities
  prob_ins = (float(match_to_ins_or_del_count)/(match_to_match+match_to_ins_or_del_count))/2.0
  prob_ext = (float(insert_to_insert_count)/(insert_to_insert_count+insert_to_match))/2.0

  #Now calculate the updated emission probabilities
  insert_emission = {}
  for aa in AA_LIST:
    insert_emission[aa] = float(aa_count[aa])/sum(aa_count.values())

  match_emission = {}
  for aa in AA_LIST:
    total_for_aa = 0
    for bb in AA_LIST:
      total_for_aa += match_counts[(aa,bb)]
    for bb in AA_LIST:
      match_emission[(aa,bb)] = float(match_counts[(aa,bb)])/total_for_aa


  return prob_ins, prob_ext, insert_emission, match_emission


def printLogProbs(prob_ins, prob_ext, insert_emission, match_emission
  , outputfile, total_log_likelihood, verbose):

  if verbose:
    print "Writing output required for mosaic to: ", outputfile

  with open(outputfile, 'w') as outfile:
    outfile.write("Total log likelihood: " + str(total_log_likelihood) + "\n\n")
    outfile.write("Insert probability (del): " + str(prob_ins) + "\n")
    outfile.write("Extension probability (eps): " + str(prob_ext) + "\n\n\n")

    outfile.write("Log emission from insert state probabilities:\n")
    for aa in AA_LIST:
      outfile.write(str(math.log(insert_emission[aa]))+",")

    outfile.write("\n\n\nLog match emission propabilities:\n")

    for aa in AA_LIST:
      outfile.write("{")
      outfile.write(",".join([str(math.log(match_emission[(aa,bb)])) for bb in AA_LIST]))
      outfile.write("},\\\n")

  return


def check_runs(number_runs, align_files, verbose):

  if verbose:
    print "Checking runs completed successfully..."
  
  log_file_ext = align_files[0].split("_run")[0]
  
  total_log_likelihood = 0.0

  for i in range(1, number_runs+1):
    logfile = log_file_ext + "_run" +str(i) + ".fasta_output.log"
    if verbose:
      print "Checking file: ", logfile
    with open(logfile, 'rU') as infile:
      finished=False
      for line in infile:
        if "Maximum log-likelihood =" in line:
          total_log_likelihood+=float(line.split("=")[1].strip())
        if "Program completed in" in line:
          finished=True
    assert finished, "Missing run " + logfile

  return total_log_likelihood




def main():

  parser = argparse.ArgumentParser(description='Calculate transition probabilities from Viterbi algorithm')

  parser.add_argument('--num_runs', dest='num_runs'
    , help="number of runs we are expecting. (error check)"
    , type=int
    , required=True)

  parser.add_argument('--out', dest='outputfile'
    , help="location of output file."
    , required=True)

  parser.add_argument('--align', nargs='+'
    , dest='align_files'
    , help='location of alignment files from Mosaic Viterbi run.')

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , default=False
    , help='print verbose output (default=False)')

  args = parser.parse_args()

  total_log_likelihood = check_runs(args.num_runs, args.align_files
    , args.verbose)

  alignments=[]
  for a in args.align_files:
    alignments = alignments + loadCounts(a, args.verbose)

  prob_ins, prob_ext, insert_emission, match_emission = calculate_probabilities(alignments
    , args.verbose)

  printLogProbs(prob_ins, prob_ext, insert_emission, match_emission
  , args.outputfile, total_log_likelihood, args.verbose)


  return


if __name__ == '__main__':
  main()