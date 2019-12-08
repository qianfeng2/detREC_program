/*******************************************************************
Data structures and functions for reading sequences from input files
*******************************************************************/


#include "seqtools.h"

#define DEBUG_SEQ 0

struct data * read_fasta(struct pars *my_pars, int aligned) {

	int i, clen=2, pos;
	char c;
	struct sequence *new_seq;
	struct data *ndata;
	FILE *ifp;

	printf("\n\n*** Reading sequences from input file: %s ***\n\n", my_pars->input_filename);

	ifp = fopen(my_pars->input_filename, "r");
	if (!ifp) {
		printf("\n\n***Error: Cannot open input file ***\n\n");
		exit(1);
	}

	ndata = (struct data *) malloc((size_t) sizeof(struct data));
	if (!ndata) {
		printf("\n\n*** Error: cannot allocate memory for data ***\n\n");
		exit(1);
	}
	ndata->maxlseq=0;
	ndata->nseq=0;
	ndata->type=my_pars->type;
	ndata->aligned=aligned;
	ndata->seqs = (struct sequence **) malloc((size_t) sizeof(struct sequence *));
	if (!ndata->seqs) {
		printf("\n\n*** Error: cannot allocate memory for sequences ***\n\n");
		exit(1);
	}

	if (DEBUG_SEQ) {
		printf("\nAllocated initial memory");
		fflush(stdout);
	}

	while ((c=fgetc(ifp))!=EOF) {

	/*Ignore lines starting with # as comments*/

		if (c == '#') {
			while ((c=fgetc(ifp)) != '\n');
		}
		
		if (c == '>') {
			ndata->nseq++;

			if (DEBUG_SEQ) {
				printf("\nFound new sequence (%i) - reallocating memory\n", ndata->nseq);
				fflush(stdout);
			}

			ndata->seqs = (struct sequence **) realloc(ndata->seqs, (ndata->nseq+1)*sizeof(struct sequence *));
			new_seq = (struct sequence *) malloc((size_t) sizeof(struct sequence));
			new_seq->seq = (int *) malloc((size_t) clen*sizeof(int));
			new_seq->name = (char *) malloc((size_t) MAXNAME*sizeof(char));
			new_seq->length = 0;
			new_seq->name = fgets(new_seq->name, MAXNAME, ifp);
			new_seq->group=1;
			for (i=0;i<MAXNAME;i++) {
				if (new_seq->name[i]=='\n' || new_seq->name[i]==' ' ||  new_seq->name[i]==13) new_seq->name[i]='\0';
			}
			if (DEBUG_SEQ) {
				printf("\nReading new sequence: %s", new_seq->name);
				fflush(stdout);
			}
			ndata->seqs[ndata->nseq]=new_seq;
			pos=0;
		}

		else if (c == ' ' || c =='\t' || c == '\n' || c == '*' || c == 13) ;  /*Spacing - ignore (also ignoring TERM codons)*/

		else if (ndata->type == 1) {
		
		/*Amino acids - IUPAC convention

	    Abbreviation    Amino acid name
        ------------    ---------------
        

        Ala     A 1      Alanine
        Arg     R 2      Arginine
        Asn     N 3      Asparagine
        Asp     D 4      Aspartic acid (Aspartate)
        Cys     C 5      Cysteine
        Gln     Q 6      Glutamine
        Glu     E 7      Glutamic acid (Glutamate)
        Gly     G 8      Glycine
        His     H 9      Histidine
        Ile     I 10     Isoleucine
        Leu     L 11     Leucine
        Lys     K 12     Lysine
        Met     M 13     Methionine
        Phe     F 14     Phenylalanine
        Pro     P 15     Proline
        Ser     S 16     Serine
        Thr     T 17     Threonine
        Trp     W 18     Tryptophan
        Tyr     Y 19     Tyrosine
        Val     V 20     Valine
        Asx     B 21     Aspartic acid or Asparagine
        Glx     Z 22     Glutamine or Glutamic acid.
        Xaa     X 23     Any amino acid.
        TERM      24     termination codon

		*/
			pos++;

			switch(c) {
			case '?': case '-': 
				ndata->seqs[ndata->nseq]->seq[pos]=0;
				break;
			case 'A' : case 'a':
				ndata->seqs[ndata->nseq]->seq[pos]=1;
				break;
			case 'R' : case 'r':
				ndata->seqs[ndata->nseq]->seq[pos]=2;
				break;
			case 'N' : case 'n':
				ndata->seqs[ndata->nseq]->seq[pos]=3;
				break;
			case 'D' : case 'd':
				ndata->seqs[ndata->nseq]->seq[pos]=4;
				break;
			case 'C' : case 'c':
				ndata->seqs[ndata->nseq]->seq[pos]=5;
				break;
			case 'Q' : case 'q':
				ndata->seqs[ndata->nseq]->seq[pos]=6;
				break;
			case 'E' : case 'e':
				ndata->seqs[ndata->nseq]->seq[pos]=7;
				break;
			case 'G' : case 'g':
				ndata->seqs[ndata->nseq]->seq[pos]=8;
				break;
			case 'H' : case 'h':
				ndata->seqs[ndata->nseq]->seq[pos]=9;
				break;
			case 'I' : case 'i':
				ndata->seqs[ndata->nseq]->seq[pos]=10;
				break;
			case 'L' : case 'l':
				ndata->seqs[ndata->nseq]->seq[pos]=11;
				break;
			case 'K' : case 'k':
				ndata->seqs[ndata->nseq]->seq[pos]=12;
				break;
			case 'M' : case 'm':
				ndata->seqs[ndata->nseq]->seq[pos]=13;
				break;
			case 'F' : case 'f':
				ndata->seqs[ndata->nseq]->seq[pos]=14;
				break;
			case 'P' : case 'p':
				ndata->seqs[ndata->nseq]->seq[pos]=15;
				break;
			case 'S' : case 's':
				ndata->seqs[ndata->nseq]->seq[pos]=16;
				break;
			case 'T' : case 't':
				ndata->seqs[ndata->nseq]->seq[pos]=17;
				break;
			case 'W' : case 'w':
				ndata->seqs[ndata->nseq]->seq[pos]=18;
				break;
			case 'Y' : case 'y':
				ndata->seqs[ndata->nseq]->seq[pos]=19;
				break;
			case 'V' : case 'v':
				ndata->seqs[ndata->nseq]->seq[pos]=20;
				break;
			case 'B' : case 'b':
				ndata->seqs[ndata->nseq]->seq[pos]=21;
				break;
			case 'Z' : case 'z':
				ndata->seqs[ndata->nseq]->seq[pos]=22;
				break;
			case 'X' : case 'x':
				ndata->seqs[ndata->nseq]->seq[pos]=23;
				break;
			case '.' : case '*' : /*TERM*/
				ndata->seqs[ndata->nseq]->seq[pos]=24;
				break;
			default: /*Unrecognised character*/
				printf("\n\n*** Error: unrecognised AA (%c) ***\n\n", c);
				exit(1);
			}

			ndata->seqs[ndata->nseq]->length = pos;
		}

		else { /*Nucleotide*/

			pos++;
			ndata->seqs[ndata->nseq]->length=pos;

			switch(c) {
			case '-': case '?': case 'N' : case 'n' : case 'X': case 'x':
				ndata->seqs[ndata->nseq]->seq[pos]=0;
				break;
			/*Set all IUPAC ambiguity NTs to 'missing'*/
			case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M': case 'B': case 'H': case 'D': case 'V':
				ndata->seqs[ndata->nseq]->seq[pos]=0;
				break;
			case 'r': case 'y': case 's': case 'w': case 'k': case 'm': case 'b': case 'h': case 'd': case 'v':
				ndata->seqs[ndata->nseq]->seq[pos]=0;
				break;
			case 'T': case 't': case '0': case 'U': case 'u':
				ndata->seqs[ndata->nseq]->seq[pos]=1;
				break;
			case 'C': case 'c': case '1':
				ndata->seqs[ndata->nseq]->seq[pos]=2;
				break;
			case 'A': case 'a':
				ndata->seqs[ndata->nseq]->seq[pos]=3;
				break;
			case 'G': case 'g': 
				ndata->seqs[ndata->nseq]->seq[pos]=4;
				break;
			default: /*Unrecognised character*/
				printf("\n\n*** Error: unrecognised NT (%c) ***\n\n", c);
				exit(1);
			}
		}

		/*Check to see if need to add more to sequence length*/
		if (pos>(clen-10)) {
			clen *= 2;
			ndata->seqs[ndata->nseq]->seq = (int *) realloc(ndata->seqs[ndata->nseq]->seq, clen*sizeof(int));
		}

	}

	fclose(ifp);

	return ndata;
}





double watterson(int n) {

	int i;
	double cump=1;

	for (i=2;i<n;i++) cump+=(double) 1/i;
	return cump;
}


char num2nuc(int i, int type) {

	char aa[26]={'?','A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '.'};
	char nuc[6]={'?','T', 'C', 'A', 'G'};

	if (type == 1) {
		if (i<0 || i>25) {
			printf("\n\n***Error: cannot convert number to AA ***\n\n");
			exit(0);
		}
		return aa[i];
	}
	else {
		if (i<0 || i>5) {
			printf("\n\n***Error: cannot convert number to NT ***\n\n");
			exit(0);
		}
		return nuc[i];
	}
}


void print_sequences(struct data *my_data, FILE *ofp) {

	int i, j;

	for (i=1;i<=my_data->nseq;i++) {
		fprintf(ofp,"\nSequence %i\t: %s (length = %i)\n", i, my_data->seqs[i]->name, my_data->seqs[i]->length);
		for (j=1;j<=my_data->seqs[i]->length;j++) fprintf(ofp,"%c",num2nuc(my_data->seqs[i]->seq[j], my_data->type));
	}

	fprintf(ofp,"\n\n");
}



/*Routine to print out alignemnt from pairwise routine*/

void print_pair_align(struct palign *my_align, int *sq1, int *sq2) {

	int i, lx;

	lx = my_align->length;
	
	printf("\n\nAlignment\n\n");
		for (i=1;i<=lx;i++) {
			if(my_align->global_trace[i][1]==2) printf("-");
			else printf("%c",num2nuc(sq1[my_align->global_trace[i][3]],1));
		}

		printf("\n");
		for (i=1;i<=lx;i++) {
			if (my_align->global_trace[i][1]==1) printf("|");
			else if (my_align->global_trace[i][1]==2 || my_align->global_trace[i][1]==3) printf(" ");
			else printf("x");
		}
		printf("\n");

		for (i=1;i<=lx;i++) {
			if (my_align->global_trace[i][1] == 3) printf("-");
			else printf("%c",num2nuc(sq2[my_align->global_trace[i][2]],1));
		}
		printf("\n\n");
}


void print_kalign(struct kalign *kwise, struct data *my_data, FILE *ofp) {

	int l, line=50, base=0;

	fprintf(ofp,"\n\nAlignment with recombination\n\n");

	for (l=1;l<=kwise->length;l++) {
		if (kwise->trace[1][l]==2) fprintf(ofp,"-");
		else fprintf(ofp,"%c", num2nuc(my_data->seqs[kwise->target]->seq[kwise->trace[4][l]], 1));
	}
	fprintf(ofp,"\n");
	for (l=1;l<=kwise->length;l++) {
		if (kwise->trace[1][l]==1) {
			if (my_data->seqs[kwise->target]->seq[kwise->trace[4][l]] == my_data->seqs[kwise->trace[2][l]]->seq[kwise->trace[3][l]]) fprintf(ofp,"|");
			else fprintf(ofp,".");
		}
		else if (kwise->trace[1][l]==2 || kwise->trace[1][l]==3) fprintf(ofp," ");
		else if (kwise->trace[1][l]==4) fprintf(ofp,"x");
	}
	fprintf(ofp,"\n");
	for (l=1;l<=kwise->length;l++) {
		if (kwise->trace[1][l]==3) fprintf(ofp,"-");
		else fprintf(ofp,"%c", num2nuc(my_data->seqs[kwise->trace[2][l]]->seq[kwise->trace[3][l]], 1));
	}
	fprintf(ofp,"\n");
	for (l=1;l<=kwise->length;l++) {
		fprintf(ofp,"%1i", kwise->trace[2][l]);
	}
	fprintf(ofp,"\n\n");

}



