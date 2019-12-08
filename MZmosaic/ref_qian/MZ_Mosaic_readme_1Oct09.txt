*How to Run Mosaic*

_About Mosaic_:
Mosaic is a program written in C, all the commands are in C format and given at the command line from the terminal.

_Data to run on Mosaic_:
Sequence data should be amino acids and in FASTA format.  (For an exact definition of this format please see: http://en.wikipedia.org/wiki/FASTA_format .)

NOTE: make sure the sequences read-through properly when translating to amino acid sequence. A one- or two-base insertion or deletion (essentially a "missense" sequencing error) will cause the program to return a spuriously low score for the resulting alignment for that sequence.

Place all work files in the same directory as the program.

_Getting started_:
The raw components of the program are as follows:
mosaic_fb.h 
mosaic_fb.c 
seqtools.c 
seqtools.h 
tools.c 
tools.h
makefile

(These are the different algorithms and matrices that make the program work, the "makefile" includes all the information you need to compile the program so it functions.)

To compile the program, in your terminal navigate (change directory, UNIX command "cd") to the directory you have made for the component files and type the command "make" -- now you are good to go.

_Running Mosaic on Amino Acid Data_
There are two ways you will probably want to run the program.  The first time you run it with a fresh data set you will want to run it with the expectation-maximization (EM) algorithm to estimate to estimate the parameters (more info on EM here: http://en.wikipedia.org/wiki/Expectation-maximization_algorithm).  

Here is what your command looks like to run:

./mosaic -estimate -seq input_data_file.txt -aa 

Meaning: 
./mosaic - "run program called moasic"
-estimate - "use EM estimator" [leave this out if you don't need the estimation]
-seq  - "what follows is the input sequence"
input_data_file - [name of input file placed in directory]
-aa - "this is amino acid data"

you can also add:
-tag output_file_name

This will name the output alignment file anything you want, this is useful as the program's default is to output a file titled "align.txt" and it will overwrite this file every time the program is run if the name is not changed manually or the file is not moved out of the directory.

_Mosaic Output_
Output files are most easily viewed in a web browser, such as Firefox

Alignments:
The Mosaic program returns a simple, local alignment for each sequence in the file.  The sequence assessment line includes the following information:

| = amino acid match
^ = insertion in target sequence
- = deletion in target sequence
[blank] = mismatch

Scores:
The output also includes the sequence length and a calculated likelihood score of the alignment (MLlk).  Although the score is a good indicator of the probability that the alignment is correct given the data, with the higher the score the better the probability, keep in mind that the likelihoods are also dependent on the sequence length (longer sequences will always have lower scores).  Dividing the sequence length by the absolute value of the will normalize the scores:

Length/|MLlk|=normalized score

_To Run Mosaic on only a specific set of target sequences_

The Mosaic program works by making homology-based alignments using "target" sequences and matching them to "destination" sequences.  If not specifications are made otherwise, the program will use every sequence in the data set as both target and destination sequences.  Here is how to specify sequences are one group or the other:

(1) Create a data file of ALL the sequences to be used, but add to the names tags that can make them recognizable as separate groups, i.e. tags that are shared among all members of the group and are exclusive to the group.  

(2) Create a query using the -group and -target commands in the command line, this specifies the groups that can be used as either target or destination sequences, and also which group to use as the target sequences. 

Here is an example:

./mosaic -estimate -seq input_data_file.txt -aa -group 3 3D7 HB3 reich -target reich -tag output_filename

This command would use only the reich (reichenowi sequences) as the targets, thus you would only get as many alignments as their are "reich"-tagged sequences.  And the resulting alignments, consisting of destination sequences, are composed of 3D7- and HB3-tagged sequences. 
	--> If the sequences used as target sequences should be included as destination sequences as well, duplicate the sequences using the proper group tag to be specified.
	--> Be careful to make sure that all sequences in an input file contain one of the group tags, a sequence with out a specified tag (if using the -group option) will cause an error and the program will terminate.






