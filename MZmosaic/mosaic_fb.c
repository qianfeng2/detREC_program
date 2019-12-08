
/*********************************************
 Program to detect mosaic gene structures
*********************************************/

#include "tools.h"
#include "seqtools.h"
#include "mosaic_fb.h"

#define DEBUG 0


main(int argc, char *argv[]) {

	int target;
	long seed = -setseed();
	struct data *my_data;
	struct pars *my_pars;
	struct matrices *my_matrices;
	FILE *ofp;
	clock_t start, end;

	printf("\n\n*** Running MOSAIC ***\n\nA program for detecting mosaic gene structures");
	printf("\n\nProgram returns (for each sequence) a maximum accuracy alignment path\n\n");

	start = clock();

	if (argc==1 || argv[1]=="-h" || argv[1]=="-?") print_help(stdout);

	my_pars = (struct pars *) get_pars(argc, argv);
	my_data = read_fasta(my_pars, 0);  

	printf("\n\nRead %i sequences\n\n", my_data->nseq);
	if (DEBUG) print_sequences(my_data, stdout);

	/*If groups identified, need to classify sequences*/
	if (my_pars->ngroups>1) {
		classify_sequences(my_pars, my_data);
	}

	/*Add null term to end of each sequence - makes for simpler calculations*/
	my_data = add_null_term(my_data);

	printf("\n\n*** Computing alignment paths for each sequence ***\n\n");

	/*Allocate memory*/
	my_matrices = (struct matrices *) allocate_matrices(my_data, my_pars);

	/*If required, estimate parameters using EM*/
	if (my_pars->estimate==1) {
		estimate_parameters_em_no_rec(my_data, my_pars, my_matrices);
	}

	/*If required, calculate likelihood for each sequence over grid of rho values and estimate MLE*/
	if (my_pars->grid) {
		calculate_llk_over_rho_grid(my_data, my_pars, my_matrices);
	}

	/*Now do alignments*/
	/*If target group defined, run each target sequence against rest*/
	my_pars->combined_llk=0.0;
	for (target=1; target<=my_data->nseq; target++) {
		if (my_data->seqs[target]->group == my_pars->target_group) {
			if (my_pars->ml) kalign_vt(my_data, my_pars, my_matrices, target);
			else kalign_fb(my_data, my_pars, my_matrices, target, 0);
			
		}
	}

	/*If estimated parameters, add them to alignment file*/
	if (my_pars->estimate==1 || my_pars->grid  || my_pars->verbose) {
		ofp = fopen(my_pars->alignment_file, "a");
		print_parameters(my_pars, ofp);
		fclose(ofp);
	}

	/*Deallocate memory for matrices*/
	deallocate_matrices(my_data, my_pars, my_matrices);

	printf("Combined log likelihood = %.3lf", my_pars->combined_llk);

	end = clock();

	printf("\n\n*** Program completed in %.3lf secs CPU time! ***\n\n", ((double) (end-start)/CLOCKS_PER_SEC));

	exit(0);

}


/*Have to add a null character to end of sequences for following algorithm*/

struct data * add_null_term(struct data *my_data) {

	int seq;

	for (seq=1;seq<=my_data->nseq;seq++) {
		my_data->seqs[seq]->seq = (int *) realloc(my_data->seqs[seq]->seq, (size_t) (my_data->seqs[seq]->length+2)*sizeof(int));
		my_data->seqs[seq]->seq[my_data->seqs[seq]->length+1]=0;
	}

	return my_data;

}


/************************************
Print parameters
************************************/

void print_parameters(struct pars *my_pars, FILE *ofp) {

	int i, j;

	fprintf(ofp,"Parameters\n\n");
	fprintf(ofp,"Gap initiation: %.5lf\n", my_pars->del);
	fprintf(ofp,"Gap extension:  %.5lf\n", my_pars->eps);
	fprintf(ofp,"Termination:    %.5lf\n", my_pars->term);
	fprintf(ofp,"Recombination:  %.5lf\n", my_pars->rho);
	fprintf(ofp,"fMatch:         %.5lf\n", my_pars->piM);

	fprintf(ofp,"\nState frequencies:\n");
	for (i=0;i<my_pars->nstate;i++) fprintf(ofp,"%c\t%.5lf\n",num2nuc(i, my_pars->type), my_pars->si[i]);

	fprintf(ofp,"\nEmission probabilities:\n\t");
	for (i=0;i<my_pars->nstate;i++) fprintf(ofp,"%c\t",num2nuc(i, my_pars->type));
	for (i=0;i<my_pars->nstate;i++) {
		fprintf(ofp,"\n%c\t",num2nuc(i, my_pars->type));
		for (j=0;j<my_pars->nstate;j++) fprintf(ofp,"%.5lf\t", my_pars->sm[i][j]);
	}

	fprintf(ofp,"\n\n");

	return;

}

/************************************
Get parameters from command line
************************************/

struct pars * get_pars(int argc, char *argv[]) {

	int i, j, k, l;
	FILE *ifp, *ofp;
	char *in_str;
	struct pars *my_pars;

	extern double emiss_gap_nt[NSTATE_NT];
	extern double emiss_match_nt[NSTATE_NT][NSTATE_NT];
	extern double emiss_gap_aa[NSTATE_AA];
	extern double emiss_match_aa[NSTATE_AA][NSTATE_AA];

	time_t calendar_time;
	struct tm *output_time;

	my_pars = (struct pars *) malloc((size_t) sizeof(struct pars));

	printf("\n\n*** Reading input parameters and initialising ***\n\n");

	/*Initialise with default values*/
	my_pars->del=(double) DEFAULT_DEL;
	my_pars->eps=(double) DEFAULT_EPS;
	my_pars->rho=(double) DEFAULT_REC;
	my_pars->term=(double) DEFAULT_TERM;
	my_pars->pmatch = (double) DEFAULT_PMATCH;
	my_pars->type = -1;

	my_pars->input_filename = (char *) malloc((size_t) MAXFILENAME*sizeof(char));
	my_pars->input_filename[0]='\0';
	my_pars->output_tag = (char *) malloc((size_t) MAXFILENAME*sizeof(char));
	my_pars->output_tag[0] = '\0';
	my_pars->posterior_file = (char *) malloc((size_t) (MAXFILENAME+8)*sizeof(char));
	my_pars->posterior_file[0]='\0';
	my_pars->alignment_file = (char *) malloc((size_t) (MAXFILENAME+8)*sizeof(char));
	my_pars->alignment_file[0]='\0';
	my_pars->parameter_file = (char *) malloc((size_t) (MAXFILENAME+8)*sizeof(char));
	my_pars->parameter_file[0]='\0';
	my_pars->llk_grid_file = (char *) malloc((size_t) (MAXFILENAME+8)*sizeof(char));
	my_pars->llk_grid_file[0]='\0';

	my_pars->ngroups=1;
	my_pars->target_group=1;

	my_pars->verbose=1;
	my_pars->emiss=1;
	my_pars->ml=1;
	my_pars->estimate=0;
	my_pars->combined_llk=0;

	my_pars->psum=0;

	my_pars->grid=0;
	my_pars->npts=3;
	my_pars->rho_low=1e-3;
	my_pars->rho_high=1e-1;
	my_pars->log_scale=1;


	for (k=1; k<argc; k++) if (strchr(argv[k], '-')) {
		
		in_str = argv[k];

		/*Indicates that data are AA*/
		if (strcmp(in_str, "-aa") == 0) {
			my_pars->type = 1;
			my_pars->nstate = NSTATE_AA;
			my_pars->sm = dmatrix(0, NSTATE_AA, 0, NSTATE_AA);
			my_pars->lsm = dmatrix(0, NSTATE_AA, 0, NSTATE_AA);
			for (i=0;i<NSTATE_AA;i++) for (j=0;j<NSTATE_AA;j++) {
				my_pars->lsm[i][j]=emiss_match_aa[i][j];
				my_pars->sm[i][j] = exp(emiss_match_aa[i][j]);
			}
			my_pars->si = dvector(0, NSTATE_AA);
			my_pars->lsi = dvector(0, NSTATE_AA);
			for (i=0;i<NSTATE_AA;i++) {
				my_pars->lsi[i] = emiss_gap_aa[i];
				my_pars->si[i] = exp(emiss_gap_aa[i]);
			}
			my_pars->nstate = NSTATE_AA; 
		}

		/*Indicates that data are NTs*/
		if (strcmp(in_str, "-nt") == 0) {
			my_pars->type = 0;
			my_pars->nstate = NSTATE_NT;
			my_pars->sm = dmatrix(0, NSTATE_NT, 0, NSTATE_NT);
			my_pars->lsm = dmatrix(0, NSTATE_NT, 0, NSTATE_NT);
			for (i=0;i<NSTATE_NT;i++) for (j=0;j<NSTATE_NT;j++) {
				my_pars->sm[i][j]=emiss_match_nt[i][j];
				my_pars->lsm[i][j] = log(emiss_match_nt[i][j]);
			}
			my_pars->si = dvector(0, NSTATE_NT);
			my_pars->lsi = dvector(0, NSTATE_NT);
			for (i=0;i<NSTATE_NT;i++) {
				my_pars->si[i]=emiss_gap_nt[i];
				my_pars->lsi[i] = log(emiss_gap_nt[i]);
			}
			my_pars->nstate = NSTATE_NT; 
		}

		/*Input filename*/
		if (strcmp(in_str, "-seq") == 0) 
			strncpy(my_pars->input_filename, argv[k+1], MAXFILENAME);

		/*Output filename tag*/
		if (strcmp(in_str, "-tag") == 0) 
			strncpy(my_pars->output_tag, argv[k+1], MAXFILENAME);

		/*Classify sequences into groups*/
		if (strcmp(in_str, "-group") == 0) {
			my_pars->ngroups = atoi(argv[k+1]);
			my_pars->group_identifiers = (char **) malloc((size_t) (my_pars->ngroups+1)*sizeof(char));
			for (l=1;l<=my_pars->ngroups; l++) {
				my_pars->group_identifiers[l] = (char *) malloc((size_t) MAXNAME*sizeof(char));
				strncpy(my_pars->group_identifiers[l], argv[k+l+1], MAXNAME);
			}
			printf("\nGroups identified: ");
			for (l=1;l<=my_pars->ngroups;l++) printf("%s\t", my_pars->group_identifiers[l]);
		}

		/*Identify a group of sequences for target*/
		if (strcmp(in_str, "-target") == 0) {
			my_pars->target_name = (char *) malloc((size_t) MAXFILENAME*sizeof(char));
			strncpy(my_pars->target_name, argv[k+1], MAXFILENAME);
			printf("\nTarget group: %s", my_pars->target_name);
			my_pars->target_group=0;
		}

		if (strcmp(in_str, "-rec") == 0) {
			my_pars->rho = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-del") == 0) {
			my_pars->del = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-eps") == 0) {
			my_pars->eps = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-term") == 0) {
			my_pars->term = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-e1") == 0) {
			my_pars->emiss=0;
		}

		if (strcmp(in_str, "-pmatch") == 0) {
			my_pars->pmatch = (double) atof(argv[k+1]);
		}

		if (strcmp(in_str, "-ma") == 0) {
			my_pars->ml = 0;
			printf("\nReturning maximum accuracy alignment\n");
		}

		if (strcmp(in_str, "-estimate") == 0) {
			my_pars->estimate=1;
		}

		if (strcmp(in_str, "-params") == 0) {
			my_pars->parameter_file = strcpy(my_pars->parameter_file, argv[k+1]);
		}

		if (strcmp(in_str, "-psum") == 0) {
			my_pars->psum=1;
		}

		if (strcmp(in_str, "-grid") == 0) {
			my_pars->grid = 1;
			my_pars->rho_low = (double) atof(argv[k+1]);
			my_pars->rho_high = (double) atof(argv[k+2]);
			my_pars->npts = (int) atoi(argv[k+3]);
			my_pars->log_scale = (int) atoi(argv[k+4]);
		}

		if (strcmp(in_str, "-brief") == 0) {
			my_pars->verbose=0;
		}

	}


	if (my_pars->type<0) {
		printf("\n\n*** Error: Need to indicate whether data is protein (-aa) or nucleotide (-nt) ***\n\n");
		exit(0);
	}


	/*Read parameters from file*/
	if (my_pars->parameter_file[0]!='\0') {
		ifp = fopen(my_pars->parameter_file, "r");
		if (!ifp) {
			printf("\n\n***Error: cannot open file with input parameters ***\n\n");
			exit(1);
		}
		read_params_from_file(my_pars, ifp);
		fclose(ifp);
	}

	
	/*Check that all parameters are >0 and sum is <1*/
	if (my_pars->rho<=0) my_pars->rho = 1e-32;
	if (my_pars->del<=0) my_pars->del = 1e-6;
	if (my_pars->eps<=0) my_pars->eps = 1e-6;
	if (my_pars->term<=0) my_pars->term = 1e-6;

	if (2*my_pars->del+my_pars->rho+my_pars->term >= 1) {
		printf("\n\n*** Error: event probabilities sum to >=1 ***\n\n");
		exit(1);
	}

	/*Quasi-stationary probabilities for match (M) and Insert (I)*/
	my_pars->piM = (double) (1-my_pars->eps)/(1-my_pars->eps+2*my_pars->del);

	/*Note that seems to be better to have this as separate parameter*/
	my_pars->piM = (double) DEFAULT_PIM;

	my_pars->piI = (double) 1-my_pars->piM;
	my_pars->mm = (double) 1-2*my_pars->del-my_pars->rho-my_pars->term;
	my_pars->gm = (double) 1-my_pars->eps-my_pars->rho-my_pars->term;
	my_pars->dm = 1-my_pars->eps;


	/*Make log versions*/
	my_pars->lpiM = log(my_pars->piM);
	my_pars->lpiI = log(my_pars->piI);
	my_pars->ldel=log(my_pars->del);
	my_pars->leps=log(my_pars->eps);
	my_pars->lrho=log(my_pars->rho);
	my_pars->lterm=log(my_pars->term);
	my_pars->lmm=log(my_pars->mm);
	my_pars->lgm=log(my_pars->gm);
	my_pars->ldm = log(my_pars->dm);



	/*Check input filename*/
	ifp = fopen(my_pars->input_filename, "r");
	if (!ifp) {
		printf("\n\n***Error: Cannot file input file (%s) ***\n\n", my_pars->input_filename);
		exit(1);
	}
	else {
		fclose(ifp);
	}

	/*Make output filenames*/
	if (my_pars->output_tag[0] != '\0') {
		my_pars->posterior_file = strncpy(my_pars->posterior_file, my_pars->output_tag, MAXFILENAME);
		my_pars->posterior_file = strcat(my_pars->posterior_file, "_post.txt");
		my_pars->alignment_file = strncpy(my_pars->alignment_file, my_pars->output_tag, MAXFILENAME);
		my_pars->alignment_file = strcat(my_pars->alignment_file, "_align.txt");
		my_pars->llk_grid_file = strncpy(my_pars->llk_grid_file, my_pars->output_tag, MAXFILENAME);
		my_pars->llk_grid_file = strcat(my_pars->llk_grid_file, "_grid.txt");

	}
	else {
		my_pars->posterior_file = strcat(my_pars->posterior_file,  "post.txt");
		my_pars->alignment_file = strcat(my_pars->alignment_file,  "align.txt");
		my_pars->llk_grid_file = strcat(my_pars->llk_grid_file, "grid.txt");
	}

	/*Check things about groups and targets*/
	/*NB with a single target, program runs as originally*/

	if (my_pars->target_group == 0) {

		for (k=1;k<=my_pars->ngroups;k++) {
			if (strcmp(my_pars->target_name, my_pars->group_identifiers[k])==0) {
				my_pars->target_group = k;

			}
		}
		if (my_pars->target_group==0) {
			printf("\n\n*** Error: Could not match target name with group labels ***\n\n");
			exit(1);
		}
		else {
			printf("\nTarget identified as group %i", my_pars->target_group);
		}
	}

	/*If match probs need to re-organise*/
	if (my_pars->pmatch>0) {
		for (i=0;i<my_pars->nstate;i++) {
			my_pars->si[i] = (double) 1/my_pars->nstate;
			my_pars->lsi[i] = (double) log(my_pars->si[i]);
		}
		for (i=0;i<my_pars->nstate;i++) for (j=0;j<=my_pars->nstate;j++) {
			if (i==j) {
				my_pars->sm[i][j] = (double) my_pars->pmatch;
				my_pars->lsm[i][j] = (double) log(my_pars->pmatch);
			}
			else {
				my_pars->sm[i][j] = (double) (1-my_pars->pmatch)/(my_pars->nstate-1);
				my_pars->lsm[i][j] = (double) log((1-my_pars->pmatch)/(my_pars->nstate-1));
			}
		}
	}

	/*If cut emission probs need to re-organise*/
	if (my_pars->emiss==0) {
		for (i=0;i<my_pars->nstate;i++) {
			my_pars->si[i] = 1.0;
			my_pars->lsi[i] = 0.0;
			for (j=0;j<=my_pars->nstate;j++) {
				my_pars->sm[i][j]=1.0;
				my_pars->lsm[i][j]=0.0;
			}
		}
	}

	/*Initialise some output files*/
	if (my_pars->psum) {
		ofp = fopen(my_pars->posterior_file, "w");
		if (my_pars->ml==0) fprintf(ofp,"#File containing summed posteriors for each target vs copy set\n");
		else fprintf(ofp,"#File containing summed ML paths for each target vs copy set\n");
		fprintf(ofp,"#Input parameters = ");
		for (i=0;i<argc;i++) fprintf(ofp,"%s ", argv[i]);
		time(&calendar_time);
		output_time = localtime(&calendar_time);
		fprintf(ofp,"\n#Created on %s\n", asctime(output_time)); 
		fprintf(ofp, "Sequence\tLength\tLogLk");
		fclose(ofp);
	}

	ofp = fopen(my_pars->alignment_file, "w");
	if (my_pars->ml) fprintf(ofp,"#File containing maximum likelihood alignment for each sequences\n");
	else fprintf(ofp,"#File containing maximum accuracy alignment for each sequences\n");
	fprintf(ofp,"#Input parameters = ");
	for (i=0;i<argc;i++) fprintf(ofp,"%s ", argv[i]);
	time(&calendar_time);
	output_time = localtime(&calendar_time);
	fprintf(ofp,"\n#Created on %s\n", asctime(output_time));
	fclose(ofp);

	/*Print out current model parameters*/
	print_parameters(my_pars, stdout);

	return my_pars;
}



/*Do kwise alignment with forward-backward*/

void kalign_fb(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target, int fonly) {

	int pos_target, pos_seq, seq, l1, l2, *s1, *s2, i, tmp_copy;
	static int first_target=1;
	double ***forward_m, ***backward_m, ***forward_i, ***backward_i, ***forward_d, ***backward_d;
	double llk_f, llk_b=SMALL, max_r, max_rn, llk_r, psum;
	FILE *ofp;

	printf("\rAligning sequence %s to rest using HMM", my_data->seqs[target]->name);
	fflush(stdout);

	tmp_copy = my_matrices->who_copy[target];
	my_matrices->who_copy[target]=0;

	l1=my_data->seqs[target]->length;
	s1=my_data->seqs[target]->seq;
	my_pars->sizeL=0.0;
				
	for (seq=1;seq<=my_data->nseq;seq++) if (seq!=target) my_pars->sizeL+=(double) my_data->seqs[seq]->length;
	my_pars->lsizeL = (double) log(my_pars->sizeL);

	if (DEBUG) printf("\nSizeL = %.0lf",my_pars->sizeL);

	/*Make foward and backward matrices*/
	forward_m = (double ***) my_matrices->m1_m;
	forward_i = (double ***) my_matrices->m1_i;
	forward_d = (double ***) my_matrices->m1_d;
	backward_m = (double ***) my_matrices->m2_m;
	backward_i = (double ***) my_matrices->m2_i;
	backward_d = (double ***) my_matrices->m2_d;


	/*Set everything to small*/
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		l2=my_data->seqs[seq]->length;
		for (pos_target=0;pos_target<=(l1+1);pos_target++) 
			for (pos_seq=0;pos_seq<=(l2+1); pos_seq++) {
				forward_m[seq][pos_seq][pos_target]=SMALL;
				forward_i[seq][pos_seq][pos_target]=SMALL;
				forward_d[seq][pos_seq][pos_target]=SMALL;
				backward_m[seq][pos_seq][pos_target]=SMALL;
				backward_i[seq][pos_seq][pos_target]=SMALL;
				backward_d[seq][pos_seq][pos_target]=SMALL;
			}
	}


	/*Initialise forward and backward matrices*/
			
	for (seq=1, max_r = SMALL;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2=my_data->seqs[seq]->length;
		s2=my_data->seqs[seq]->seq;

		for (pos_seq=1;pos_seq<=l2;pos_seq++){

			/*Forward Ms*/
			forward_m[seq][pos_seq][1] = (double) my_pars->lpiM-my_pars->lsizeL;
			forward_m[seq][pos_seq][1] += my_pars->lsm[s1[1]][s2[pos_seq]];
			forward_i[seq][pos_seq][1] = my_pars->lpiI-my_pars->lsizeL;
			forward_i[seq][pos_seq][1] += my_pars->lsi[s1[1]];
			if (pos_seq>1) {
				forward_d[seq][pos_seq][1] = (double) forward_m[seq][pos_seq-1][1] + log(my_pars->del + my_pars->eps*exp(forward_d[seq][pos_seq-1][1]-forward_m[seq][pos_seq-1][1]));
			}

			/*Find biggest value prior to recombination*/
			if (forward_m[seq][pos_seq][1] > max_r) max_r = forward_m[seq][pos_seq][1];
			if (forward_i[seq][pos_seq][1] > max_r) max_r = forward_i[seq][pos_seq][1];
			if (forward_d[seq][pos_seq][1] > max_r) max_r = forward_d[seq][pos_seq][1];

			/*Backward Ms*/
			backward_m[seq][pos_seq][l1] = my_pars->lterm;
			backward_i[seq][pos_seq][l1] = my_pars->lterm;
		}
	}



	/*Now loop forward over positions*/

	for (pos_target=2;pos_target<=l1;pos_target++) {

		/*First get contribution from recombination*/

		for (seq=1, llk_r=0.0, max_rn=SMALL; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				llk_r += (double) exp(forward_m[seq][pos_seq][pos_target-1]-max_r);
				llk_r += (double) exp(forward_i[seq][pos_seq][pos_target-1]-max_r);

			}
		}


		/*Match, Insert and Delete Matrices*/
		for (seq=1; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				/*Match*/
				forward_m[seq][pos_seq][pos_target] =  exp(forward_m[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->mm;
				forward_m[seq][pos_seq][pos_target] += exp(forward_i[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->gm;
				forward_m[seq][pos_seq][pos_target] += exp(forward_d[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->dm;
				forward_m[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho*my_pars->piM/my_pars->sizeL;
				forward_m[seq][pos_seq][pos_target] = (double) log(forward_m[seq][pos_seq][pos_target])+max_r;
				forward_m[seq][pos_seq][pos_target] += my_pars->lsm[s1[pos_target]][s2[pos_seq]];

				/*Insert*/
				forward_i[seq][pos_seq][pos_target] =  exp(forward_m[seq][pos_seq][pos_target-1]-max_r)*my_pars->del;
				forward_i[seq][pos_seq][pos_target] += exp(forward_i[seq][pos_seq][pos_target-1]-max_r)*my_pars->eps;
				forward_i[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho*my_pars->piI/my_pars->sizeL;
				forward_i[seq][pos_seq][pos_target] = (double) log(forward_i[seq][pos_seq][pos_target])+max_r;
				forward_i[seq][pos_seq][pos_target] += my_pars->lsi[s1[pos_target]];

				/*Delete: NB cannot terminate from a delete state*/
				if (pos_target<l1 && pos_seq>1) {
					forward_d[seq][pos_seq][pos_target] =  exp(forward_d[seq][pos_seq-1][pos_target]-max_r)*my_pars->eps;
					forward_d[seq][pos_seq][pos_target] += exp(forward_m[seq][pos_seq-1][pos_target]-max_r)*my_pars->del;
					forward_d[seq][pos_seq][pos_target] = (double) log(forward_d[seq][pos_seq][pos_target]) + max_r;
				}

				/*Get new max_r*/
				if (forward_m[seq][pos_seq][pos_target] > max_rn) max_rn = forward_m[seq][pos_seq][pos_target];
				if (forward_i[seq][pos_seq][pos_target] > max_rn) max_rn = forward_i[seq][pos_seq][pos_target];
			}
		}

		/*Completed position in target sequence*/
		max_r = max_rn;
	}
	/*Completed forward matrices*/

	if (DEBUG>1) print_forward_matrices(my_data, my_pars, my_matrices, target, stdout);



	/*Now calculate likelihood*/

	for (seq=1, llk_f=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2 = my_data->seqs[seq]->length;
			for (pos_seq=1;pos_seq<=l2;pos_seq++){ 
				llk_f += (double) exp(forward_m[seq][pos_seq][l1]-max_r);
				llk_f += (double) exp(forward_i[seq][pos_seq][l1]-max_r);
			}
	}
	llk_f = max_r + log(llk_f) +my_pars->lterm;
	printf("\nLog likelihood from forward algorithm  = %.5lf\n", llk_f);
	my_matrices->llk = llk_f;
	my_pars->combined_llk += llk_f;




	/*If only need forward, can return from this point*/
	if (fonly) {
		/*Reset copy state*/
		my_matrices->who_copy[target]=tmp_copy;
		return;
	}
	

	/*Now loop backward*/
	for (pos_target = l1-1, max_r = my_pars->lterm; pos_target>0; pos_target--) {

		/*First get contribution from recombination*/
		
		for (seq=1, llk_r=0.0, max_rn=SMALL; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

  				llk_r += (double) exp(backward_m[seq][pos_seq][pos_target+1]-max_r)*(my_pars->sm[s1[pos_target+1]][s2[pos_seq]])*(my_pars->piM);
				llk_r += (double) exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*(my_pars->si[s1[pos_target+1]])*(my_pars->piI);
			}
		}

			/*Match, Insert and Delete Matrices*/
		for (seq=1; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			/*Note: Have to do Delete state first and in reverse position order*/

			for (pos_seq=l2;pos_seq>0;pos_seq--) {
				
				/*Delete*/
				if (pos_seq<l2) {
					backward_d[seq][pos_seq][pos_target] =  exp(backward_d[seq][pos_seq+1][pos_target]-max_r)*my_pars->eps;
					backward_d[seq][pos_seq][pos_target] += exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]]*my_pars->dm;
					backward_d[seq][pos_seq][pos_target] = (double) log(backward_d[seq][pos_seq][pos_target]) + max_r;
				}

				/*Insert*/
				backward_i[seq][pos_seq][pos_target] =  exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*my_pars->eps*my_pars->si[s1[pos_target+1]];
				backward_i[seq][pos_seq][pos_target] += exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->gm*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]];
				backward_i[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho/my_pars->sizeL;
				backward_i[seq][pos_seq][pos_target] = (double) log(backward_i[seq][pos_seq][pos_target])+max_r;

				/*Match*/
				backward_m[seq][pos_seq][pos_target] =  exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->mm*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]];
				backward_m[seq][pos_seq][pos_target] += exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*my_pars->del*my_pars->si[s1[pos_target+1]];
				backward_m[seq][pos_seq][pos_target] += exp(backward_d[seq][pos_seq+1][pos_target]-max_r)*my_pars->del;
				backward_m[seq][pos_seq][pos_target] += (double) llk_r*my_pars->rho/my_pars->sizeL;
				backward_m[seq][pos_seq][pos_target] = (double) log(backward_m[seq][pos_seq][pos_target])+max_r;

				/*Get new max_r*/
				if (backward_m[seq][pos_seq][pos_target] > max_rn) max_rn = backward_m[seq][pos_seq][pos_target];
				if (backward_i[seq][pos_seq][pos_target] > max_rn) max_rn = backward_i[seq][pos_seq][pos_target];
			}
		}
		/*Completed position in target sequence*/
		max_r = max_rn;
	}
	/*Completed backward matrices*/

	/*Calculate backwards likelihood*/
	for (seq=1, llk_f=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;
			for (pos_seq=1;pos_seq<=l2;pos_seq++){ 
				llk_f += (double) exp(backward_m[seq][pos_seq][1]-max_r)*my_pars->sm[s1[1]][s2[pos_seq]]*my_pars->piM;
				llk_f += (double) exp(backward_i[seq][pos_seq][1]-max_r)*my_pars->si[s1[1]]*my_pars->piI;
			}
	}
	llk_f = max_r + log(llk_f) - log(my_pars->sizeL);
	printf("Log likelihood from backward algorithm = %.5lf\n\n", llk_f);



	if (DEBUG>1) print_backward_matrices(my_data, my_pars, my_matrices, target, stdout);



	/*Now calculate forward-backward probs*/

	for (seq=1;seq<=my_data->nseq;seq++) my_matrices->ppsum_match[seq]=0.0;
	for (pos_target=0;pos_target<=(l1+1);pos_target++) {
		for (i=1;i<=3;i++) my_matrices->ppsum_state[pos_target][i]=0.0;
		for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2=my_data->seqs[seq]->length;
			for (pos_seq=0;pos_seq<=(l2+1);pos_seq++) {
				forward_m[seq][pos_seq][pos_target] = (double) exp(forward_m[seq][pos_seq][pos_target]+backward_m[seq][pos_seq][pos_target]-llk_f);
				forward_i[seq][pos_seq][pos_target] = (double) exp(forward_i[seq][pos_seq][pos_target]+backward_i[seq][pos_seq][pos_target]-llk_f);
				forward_d[seq][pos_seq][pos_target] = (double) exp(forward_d[seq][pos_seq][pos_target]+backward_d[seq][pos_seq][pos_target]-llk_f);
			}
		}
	}
	

	if (DEBUG>1) print_posterior_matrices(my_data, my_pars, my_matrices, target);

	/*Print summed match posteriors to file*/
	if (my_pars->psum) {
		ofp = fopen(my_pars->posterior_file, "a");
		if (first_target) {
			for (i=1;i<=my_data->nseq;i++) if (my_data->seqs[i]->group == my_pars->target_group) {
				fprintf(ofp,"\t%s", my_data->seqs[i]->name);
			}
			fprintf(ofp,"\n");
			first_target=0;
		}
		fprintf(ofp, "%s\t%i\t%.3lf\t", my_data->seqs[target]->name, my_data->seqs[target]->length, my_matrices->llk);
		for (seq=1;seq<=my_data->nseq;seq++) {
			if (my_matrices->who_copy[seq]) {
				l2=my_data->seqs[seq]->length;
				for (pos_target=0, psum=0.0;pos_target<=l1;pos_target++) for (pos_seq=0;pos_seq<=l2;pos_seq++){
					psum += (double) forward_m[seq][pos_seq][pos_target];
				}
				fprintf(ofp,"%.3lf\t", psum);
			}
			else if (tmp_copy==1) {
				fprintf(ofp,"0.000\t");
			}
		}
		fprintf(ofp,"\n");
		fclose(ofp);
	}


	if (my_pars->verbose) print_max_acc_alignment(my_data, my_pars, my_matrices, target);


	/*Reset copy state*/
	my_matrices->who_copy[target]=tmp_copy;

	return;


}


/*Assign sequences to groups*/

void classify_sequences(struct pars *my_pars, struct data *my_data) {

	int seq, j;

	for (seq=1;seq<=my_data->nseq;seq++) {
		my_data->seqs[seq]->group=0;
		for (j=1;j<=my_pars->ngroups;j++) {
			if (strstr(my_data->seqs[seq]->name, my_pars->group_identifiers[j])) {
				my_data->seqs[seq]->group=j;
			}
		}
		if (my_data->seqs[seq]->group<1) {
			printf("\n\n***Error: Could not assign sequence %s to group ***\n\n", my_data->seqs[seq]->name);
			fflush(stdout);
			exit(1);
		}
	}

}




struct matrices * allocate_matrices(struct data *my_data, struct pars *my_pars) {

	int seq, maxl, ct=0, tmp;
	struct matrices *my_matrices;

	my_matrices = (struct matrices *) malloc((size_t) sizeof(struct matrices));

	for (seq=2, maxl=my_data->seqs[1]->length;seq<=my_data->nseq; seq++) {
		if (my_data->seqs[seq]->length>maxl) maxl = my_data->seqs[seq]->length;
	}
	my_matrices->maxl = maxl;

	my_matrices->m1_m = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m1_i = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m1_d = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m2_m = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m2_i = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));
	my_matrices->m2_d = (double ***) malloc((size_t) (my_data->nseq+1)*sizeof(double **));

	for (seq=1;seq<=my_data->nseq;seq++) if (my_pars->ngroups==1 || my_data->seqs[seq]->group != my_pars->target_group) {
		my_matrices->m1_m[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m1_i[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m1_d[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m2_m[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m2_i[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		my_matrices->m2_d[seq] = dmatrix(0, maxl+1, 0, maxl+1);
		ct++;
	}

	/*Make list of who can be copied from*/
	my_matrices->who_copy = (int *) malloc((size_t) (my_data->nseq+1)*sizeof(int));
	for (seq=1;seq<=my_data->nseq;seq++) {
		if (my_pars->ngroups==1 || my_data->seqs[seq]->group != my_pars->target_group) 
			my_matrices->who_copy[seq]=1;
		else my_matrices->who_copy[seq]=0;
	}


	/*Make paths for best alignment*/
	my_matrices->maxpath_copy = (int *) malloc((size_t) (2*maxl+1)*sizeof(int));
	my_matrices->maxpath_state = (int *) malloc((size_t) (2*maxl+1)*sizeof(int));
	my_matrices->maxpath_pos = (int *) malloc((size_t) (2*maxl+1)*sizeof(int));
	my_matrices->ppsum_match = (double *) malloc((size_t) (my_data->nseq+1)*sizeof(double));
	my_matrices->ppsum_state = (double **) dmatrix(0, maxl+1, 1, 3);

	tmp = (int) log10(maxl) + 1;
	my_matrices->tb_divisor = (double) pow(10, (double) tmp);

	printf("\n\n*** Allocated memory for matrices (%i) ***\n\n", ct);
	fflush(stdout);

	return my_matrices;
}



/*Deallocate memory for matrices*/

void deallocate_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices) {

	int i;

	for (i=1;i<=my_data->nseq; i++) {
		if (my_pars->ngroups==1 || my_data->seqs[i]->group != my_pars->target_group) {
			free_dmatrix(my_matrices->m1_m[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m1_i[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m1_d[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m2_m[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m2_i[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
			free_dmatrix(my_matrices->m2_d[i], 0, my_matrices->maxl+1, 0, my_matrices->maxl+1);
		}
	}

	free(my_matrices->m1_m);
	free(my_matrices->m1_i);
	free(my_matrices->m1_d);
	free(my_matrices->m2_m);
	free(my_matrices->m2_i);
	free(my_matrices->m2_d);

	free(my_matrices->who_copy);

	free(my_matrices->maxpath_copy);
	free(my_matrices->maxpath_state);
	free(my_matrices->maxpath_pos);
	free(my_matrices->ppsum_match);
	free_dmatrix(my_matrices->ppsum_state, 0, my_matrices->maxl+1, 1, 3);

	free(my_matrices);

	printf("\n\n*** Memory deallocated ***\n\n");

	return;
}



/* To print max accuracy alginment*/


void print_max_acc_alignment(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	int seq, pos_seq, pos_target, who_max, state_max, pos_max, l1, l2, cp, i, j;
	double max_acc, max_acc_rec, max_acc_rec_n;
	double ***post_m, ***post_i, ***post_d;
	FILE *ofp;

	/*Redirect*/
	post_m = (double ***) my_matrices->m1_m;
	post_i = (double ***) my_matrices->m1_i;
	post_d = (double ***) my_matrices->m1_d;

	/*First calculate maximum acc alignment*/

	l1 = my_data->seqs[target]->length;
	
	/*Initialise for delete state and get max_acc_rec for first position*/
	for (seq=1, max_acc_rec=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		
		l2 = my_data->seqs[seq]->length;

		for (pos_seq=1;pos_seq<=l2;pos_seq++) {

			if (post_m[seq][pos_seq][1]>max_acc_rec) max_acc_rec = post_m[seq][pos_seq][1];
			if (post_i[seq][pos_seq][1]>max_acc_rec) max_acc_rec = post_i[seq][pos_seq][1];

			/*Delete state*/
			if (post_m[seq][pos_seq-1][1] > post_d[seq][pos_seq-1][1]) {
				post_d[seq][pos_seq][1] += post_m[seq][pos_seq-1][1];
			}
			else {
				post_d[seq][pos_seq][1] += post_d[seq][pos_seq-1][1];
			}
		}
	}

	for (pos_target=2;pos_target<=l1;pos_target++) {

		for (seq=1, max_acc_rec_n=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				/*Match state*/
				if (post_d[seq][pos_seq-1][pos_target-1]>max_acc_rec) {
					post_m[seq][pos_seq][pos_target] += post_d[seq][pos_seq-1][pos_target-1];
				}
				else {
					post_m[seq][pos_seq][pos_target] += max_acc_rec;
				}

				/*Insert state*/
				post_i[seq][pos_seq][pos_target] += max_acc_rec;


				/*Delete state*/
				if (post_m[seq][pos_seq-1][pos_target] > post_d[seq][pos_seq-1][pos_target]) {
					post_d[seq][pos_seq][pos_target] += post_m[seq][pos_seq-1][pos_target];
				}
				else {
					post_d[seq][pos_seq][pos_target] += post_d[seq][pos_seq-1][pos_target];
				}
				
				/*Find new maximum for next step*/
				if (post_m[seq][pos_seq][pos_target]>max_acc_rec_n) max_acc_rec_n = post_m[seq][pos_seq][pos_target];
				if (post_i[seq][pos_seq][pos_target]>max_acc_rec_n) max_acc_rec_n = post_i[seq][pos_seq][pos_target];


			}
		}
		max_acc_rec = max_acc_rec_n;
	}

	if (DEBUG>1) {
		printf("\n\nSummed posteriors\n\n");
		print_posterior_matrices(my_data, my_pars, my_matrices, target);
	}

	/*Now work backwards finding the best alignment*/

	/*First find best value in end state*/
	cp = 2*my_matrices->maxl;
	for (seq=1, max_acc=0.0; seq<=my_data->nseq; seq++) if (my_matrices->who_copy[seq]) {
		l2 = my_data->seqs[seq]->length;
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			if (post_m[seq][pos_seq][l1]>max_acc) {
				max_acc = post_m[seq][pos_seq][l1];
				who_max = seq;
				state_max = 1;
				pos_max = pos_seq;
			}
			if (post_i[seq][pos_seq][l1]>max_acc) {
				max_acc = post_i[seq][pos_seq][l1];
				who_max = seq;
				state_max = 2;
				pos_max = pos_seq;
			}
		}
	}
	my_matrices->maxpath_copy[cp] = who_max;
	my_matrices->maxpath_state[cp] = state_max;
	my_matrices->maxpath_pos[cp] = pos_max;

	if (DEBUG) printf("\nEnding in state %i for sequence %i at position %i (%.3lf)", state_max, who_max, pos_max, max_acc);

	pos_target=l1;
	while(pos_target>=1) {

		/*Currently in delete state - can only come from delete or match state*/
		if (my_matrices->maxpath_state[cp]==3) {
			max_acc = post_d[my_matrices->maxpath_copy[cp]][my_matrices->maxpath_pos[cp]-1][pos_target];
			who_max = my_matrices->maxpath_copy[cp];
			state_max = 3;
			pos_max = my_matrices->maxpath_pos[cp]-1;
			if (post_m[who_max][my_matrices->maxpath_pos[cp]-1][pos_target]>=max_acc) {
				state_max = 1;
			}
		}

		/*If from I then just look at marginal best*/
		/*Suggests that allowing recombination from delete state would be a simplifying idea*/
		else if (my_matrices->maxpath_state[cp]==2) {

			for (seq=1, max_acc=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

				l2 = my_data->seqs[seq]->length;

				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					if (post_m[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_m[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 1;
						pos_max = pos_seq;
					}
					if (post_i[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_i[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 2;
						pos_max = pos_seq;
					}
				}
			}
		}

		/*If from M then  look at marginal best and delete state*/
		else if (my_matrices->maxpath_state[cp]==1) {

			for (seq=1, max_acc=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

				l2 = my_data->seqs[seq]->length;

				for (pos_seq=1;pos_seq<=l2;pos_seq++) {
					if (post_m[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_m[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 1;
						pos_max = pos_seq;
					}
					if (post_i[seq][pos_seq][pos_target-1]>max_acc) {
						max_acc = post_i[seq][pos_seq][pos_target-1];
						who_max = seq;
						state_max = 2;
						pos_max = pos_seq;
					}
				}
			}

			if (post_d[my_matrices->maxpath_copy[cp]][my_matrices->maxpath_pos[cp]-1][pos_target-1] > max_acc) {
				max_acc = post_d[my_matrices->maxpath_copy[cp]][my_matrices->maxpath_pos[cp]-1][pos_target-1];
				who_max = my_matrices->maxpath_copy[cp];
				state_max = 3;
				pos_max = my_matrices->maxpath_pos[cp]-1;
			}
		}

		cp--;
		if (cp<=0) {
			printf("\n\n***Error: reconstructed path longer than maximum possible ***\n\n");
			exit(1);
		}
		my_matrices->maxpath_copy[cp] = who_max;
		my_matrices->maxpath_state[cp] = state_max;
		my_matrices->maxpath_pos[cp] = pos_max;

		if (my_matrices->maxpath_state[cp+1] != 3) pos_target--;

		if (DEBUG) printf("\nNext = state %i for sequence %i at position %i (%.3lf).  Pos target = %i", state_max, who_max, pos_max, max_acc, pos_target);

	}

	if (DEBUG) {

		printf("\n\nMax Accuracy alignment\n\nWho\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_copy[i]);
		printf("\nState\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_state[i]);
		printf("\nPos\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_pos[i]);
		printf("\n\n");
	}


	ofp = fopen(my_pars->alignment_file, "a");
	fprintf(ofp,"\nTarget: %s\tLength: %i\tLlk: %.3lf\n",my_data->seqs[target]->name, my_data->seqs[target]->length, my_matrices->llk);

	/*First print target sequence*/
	cp++;
	fprintf(ofp,"%10s\t", my_data->seqs[target]->name);
	for (i=cp,pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==3) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[target]->seq[pos_target], my_data->type));
			pos_target++;
		}
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do matching*/
	fprintf(ofp,"          \t");
	for (i=cp, pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==1) {
			if (my_data->seqs[target]->seq[pos_target] == my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]]) fprintf(ofp,"|");
			else fprintf(ofp," ");
			pos_target++;
		}
		else if (my_matrices->maxpath_state[i]==2) {
			pos_target++;
			fprintf(ofp,"^");
		}
		else fprintf(ofp,"~");
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do copy tracks - switch whenever gets to new value*/
	fprintf(ofp,"%10s\t",my_data->seqs[my_matrices->maxpath_copy[cp]]->name);
	for (i=cp;i<=2*my_matrices->maxl;i++) {

		/*Check to see if need to make recombination event*/
		if (i>cp && my_matrices->maxpath_copy[i]!=my_matrices->maxpath_copy[i-1]) {
			fprintf(ofp,"\n%10s\t", my_data->seqs[my_matrices->maxpath_copy[i]]->name);
			for (j=1;j<=(i-cp);j++) fprintf(ofp," ");
		}

		if (my_matrices->maxpath_state[i]==2) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]], my_data->type));
		}
	}

	fprintf(ofp,"\n\n");
	fflush(ofp);
	fclose(ofp);


}





void print_posterior_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	int seq, pos, pos2;
	double *site_sums;

	printf("\n\nPosterior matrices for sequence %s\n", my_data->seqs[target]->name);

	site_sums = (double *) malloc((size_t) (my_data->seqs[target]->length+1)*sizeof(double));
	for (pos=1;pos<=my_data->seqs[target]->length; pos++) site_sums[pos]=0.0;

	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		printf("\n\tSequence %s: Match\n", my_data->seqs[seq]->name);
		for (pos2=1;pos2<=my_data->seqs[seq]->length;pos2++) {
			printf("\t%i\t", pos2);
			for (pos=1;pos<=my_data->seqs[target]->length;pos++) {
				printf("%.4lf\t", my_matrices->m1_m[seq][pos2][pos]);
				site_sums[pos] += my_matrices->m1_m[seq][pos2][pos];
			}
			printf("\n");
		}
		printf("\n\tSequence %s: Insert\n", my_data->seqs[seq]->name);
		for (pos2=1;pos2<=my_data->seqs[seq]->length;pos2++) {
			printf("\t%i\t", pos2);
			for (pos=1;pos<=my_data->seqs[target]->length;pos++) {
				printf("%.4lf\t", my_matrices->m1_i[seq][pos2][pos]);
				site_sums[pos] += my_matrices->m1_i[seq][pos2][pos];
			}
			printf("\n");
		}
		printf("\n\tSequence %s: Delete\n", my_data->seqs[seq]->name);
		for (pos2=1;pos2<=my_data->seqs[seq]->length;pos2++) {
			printf("\t%i\t", pos2);
			for (pos=1;pos<=my_data->seqs[target]->length;pos++) {
				printf("%.4lf\t", my_matrices->m1_d[seq][pos2][pos]);
			}
			printf("\n");
		}
	}

	printf("\n\nSite specific posterior sums\n\nPosition\tSum\n");
	for (pos=1; pos<=my_data->seqs[target]->length;pos++) {
		printf("%i\t%.4lf\n", pos, site_sums[pos]);
	}

	printf("\n\n");


	free(site_sums);

	return;


}



/*Print out forward matrices*/

void print_forward_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target, FILE *ofp) {

	int seq, pos_target, pos_seq, l1, l2, *s1, *s2;


	fprintf(ofp, "\n\n*** Forward matrices for sequence: %s ***\n\n", my_data->seqs[target]->name);

	l1 = my_data->seqs[target]->length;
	s1 = my_data->seqs[target]->seq;

	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		fprintf(ofp,"\n\nCopy_target %s\n\n", my_data->seqs[seq]->name);

		l2 = my_data->seqs[seq]->length;
		s2 = my_data->seqs[seq]->seq;

		fprintf(ofp,"Match %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++) 
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");
		
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m1_m[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nInsert %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++) 
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");
		
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m1_i[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nDelete %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++) 
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");
		
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m1_d[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");

		}
	}

}



/*Print out forward matrices*/

void print_backward_matrices(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target, FILE *ofp) {

	int seq, pos_target, pos_seq, l1, l2, *s1, *s2;


	fprintf(ofp, "\n\n*** Backward matrices for sequence: %s ***\n\n", my_data->seqs[target]->name);

	l1 = my_data->seqs[target]->length;
	s1 = my_data->seqs[target]->seq;

	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		fprintf(ofp,"\n\nCopy_target %s\n\n", my_data->seqs[seq]->name);

		l2 = my_data->seqs[seq]->length;
		s2 = my_data->seqs[seq]->seq;

		fprintf(ofp,"Match %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++) 
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");
		
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m2_m[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nInsert %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++) 
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");
		
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m2_i[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");
		}

		fprintf(ofp,"\n\nDelete %s\n\n\t", my_data->seqs[seq]->name);
		for (pos_target = 1; pos_target<=l1; pos_target++) 
			fprintf(ofp,"%c\t", num2nuc(s1[pos_target], my_data->type));
		fprintf(ofp,"\n");
		
		for (pos_seq=1;pos_seq<=l2;pos_seq++) {
			fprintf(ofp,"%c\t", num2nuc(s2[pos_seq], my_data->type));
			for (pos_target=1;pos_target<=l1;pos_target++) fprintf(ofp,"%.4lf\t", my_matrices->m2_d[seq][pos_seq][pos_target]);
			fprintf(ofp,"\n");

		}
	}

}








/*Do kwise alignment with Viterbi*/

void kalign_vt(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	int pos_target, pos_seq, seq, l1, l2, *s1, *s2, tmp_copy;
	int who_max, state_max, pos_max, who_max_n, state_max_n, pos_max_n;
	int who_next, state_next, pos_next, cp, i, j;
	double ***vt_m, ***vt_i, ***vt_d, ***tb_m, ***tb_i, ***tb_d;
	double max_r, max_rn;
	char tmp_name[MAXNAME];
	FILE *ofp;

/*	printf("\rAligning sequence %5i to rest using ML", target);
	fflush(stdout);
*/	

	tmp_copy = my_matrices->who_copy[target];
	my_matrices->who_copy[target]=0;

	l1=my_data->seqs[target]->length;
	s1=my_data->seqs[target]->seq;
	my_pars->sizeL=0.0;
	for (seq=1;seq<=my_data->nseq;seq++) if (seq!=target) my_pars->sizeL+=(double) my_data->seqs[seq]->length;
	my_pars->lsizeL = (double) log(my_pars->sizeL);

	if (DEBUG) printf("\nSizeL = %.0lf",my_pars->sizeL);

	/*Make Viterbi and traceback matrices*/
	vt_m = (double ***) my_matrices->m1_m;
	vt_i = (double ***) my_matrices->m1_i;
	vt_d = (double ***) my_matrices->m1_d;
	tb_m = (double ***) my_matrices->m2_m;
	tb_i = (double ***) my_matrices->m2_i;
	tb_d = (double ***) my_matrices->m2_d;

	/*Set everything to small*/
	
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		l2=my_data->seqs[seq]->length;
		for (pos_target=0;pos_target<=(l1+1);pos_target++) 
			for (pos_seq=0;pos_seq<=(l2+1); pos_seq++) {
				vt_m[seq][pos_seq][pos_target]=SMALL;
				vt_i[seq][pos_seq][pos_target]=SMALL;
				vt_d[seq][pos_seq][pos_target]=SMALL;
				tb_m[seq][pos_seq][pos_target]=0;
				tb_i[seq][pos_seq][pos_target]=0;
				tb_d[seq][pos_seq][pos_target]=0;
			}
	}
	

	/*Initialise Viterbi matrices*/
			
	for (seq=1, max_r = SMALL; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2=my_data->seqs[seq]->length;
		s2=my_data->seqs[seq]->seq;

		for (pos_seq=1;pos_seq<=l2;pos_seq++){

			/*Viterbi Ms*/
			vt_m[seq][pos_seq][1] = (double) my_pars->lpiM-my_pars->lsizeL;
			vt_m[seq][pos_seq][1] += my_pars->lsm[s1[1]][s2[pos_seq]];
			vt_i[seq][pos_seq][1] = my_pars->lpiI-my_pars->lsizeL;
			vt_i[seq][pos_seq][1] += my_pars->lsi[s1[1]];

			if (pos_seq>1) {
				vt_d[seq][pos_seq][1] = (double) vt_m[seq][pos_seq-1][1] + my_pars->ldel;
				tb_d[seq][pos_seq][1] = (double) seq*10 + 1 + (pos_seq-1)/my_matrices->tb_divisor;
				if ((vt_d[seq][pos_seq-1][1] + my_pars->leps) > vt_d[seq][pos_seq][1]) {
					vt_d[seq][pos_seq][1] = (double) vt_d[seq][pos_seq-1][1] + my_pars->leps;
					tb_d[seq][pos_seq][1] = (double) seq*10 + 3 + (pos_seq-1)/my_matrices->tb_divisor;
				}
			}

			/*Find biggest value prior to recombination*/
			if (vt_m[seq][pos_seq][1] > max_r) {
				max_r = vt_m[seq][pos_seq][1];
				who_max = seq;
				state_max = 1;
				pos_max = pos_seq;
			}
			if (vt_i[seq][pos_seq][1] > max_r) {
				max_r = vt_i[seq][pos_seq][1];
				who_max = seq;
				state_max = 2;
				pos_max = pos_seq;
			}
		}
	}



	/*Now loop forward over positions*/

	for (pos_target=2;pos_target<=l1;pos_target++) {

		printf("\rAligning sequence %s to rest using ML: Pos %5i", my_data->seqs[target]->name, pos_target);
		fflush(stdout);

		/*Match, Insert and Delete Matrices*/
		for (seq=1, max_rn = SMALL+max_r; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				/*Match*/
				/*Initialise with recombination*/
				vt_m[seq][pos_seq][pos_target] = max_r + my_pars->lrho + my_pars->lpiM - my_pars->lsizeL;
				tb_m[seq][pos_seq][pos_target] = (double) who_max*10 + state_max + pos_max/my_matrices->tb_divisor;

				/*Compare to MM*/
				if ((vt_m[seq][pos_seq-1][pos_target-1] + my_pars->lmm) > vt_m[seq][pos_seq][pos_target]) {
					vt_m[seq][pos_seq][pos_target] = vt_m[seq][pos_seq-1][pos_target-1] + my_pars->lmm;
					tb_m[seq][pos_seq][pos_target] = (double) seq*10 + 1 + (pos_seq-1)/my_matrices->tb_divisor;
				}
				/*Compare to IM*/
				if ((vt_i[seq][pos_seq-1][pos_target-1] + my_pars->lgm) > vt_m[seq][pos_seq][pos_target]) {
					vt_m[seq][pos_seq][pos_target] = vt_i[seq][pos_seq-1][pos_target-1] + my_pars->lgm;
					tb_m[seq][pos_seq][pos_target] = (double) seq*10 + 2 + (pos_seq-1)/my_matrices->tb_divisor;
				}

				/*Compare to DM*/
				if ((vt_d[seq][pos_seq-1][pos_target-1] + my_pars->ldm) > vt_m[seq][pos_seq][pos_target]) {
					vt_m[seq][pos_seq][pos_target] = vt_d[seq][pos_seq-1][pos_target-1] + my_pars->ldm;
					tb_m[seq][pos_seq][pos_target] = (double) seq*10 + 3 + (pos_seq-1)/my_matrices->tb_divisor;
				}

				/*Add in state match*/
				vt_m[seq][pos_seq][pos_target] += my_pars->lsm[s1[pos_target]][s2[pos_seq]];


				/*Insert*/
				/*Initialise with recombination*/
				vt_i[seq][pos_seq][pos_target] = max_r + my_pars->lrho + my_pars->lpiI - my_pars->lsizeL;
				tb_i[seq][pos_seq][pos_target] = (double) who_max*10 + state_max + pos_max/my_matrices->tb_divisor;

				/*Compare to MI*/
				if ((vt_m[seq][pos_seq][pos_target-1] + my_pars->ldel) > vt_i[seq][pos_seq][pos_target]) {
					vt_i[seq][pos_seq][pos_target] = vt_m[seq][pos_seq][pos_target-1] + my_pars->ldel;
					tb_i[seq][pos_seq][pos_target] = (double) seq*10 + 1 + (pos_seq)/my_matrices->tb_divisor;
				}
				/*Compare to II*/
				if ((vt_i[seq][pos_seq][pos_target-1] + my_pars->leps) > vt_i[seq][pos_seq][pos_target]) {
					vt_i[seq][pos_seq][pos_target] = vt_i[seq][pos_seq][pos_target-1] + my_pars->leps;
					tb_i[seq][pos_seq][pos_target] = (double) seq*10 + 2 + (pos_seq)/my_matrices->tb_divisor;
				}

				/*Add in state insert*/
				vt_i[seq][pos_seq][pos_target] += my_pars->lsi[s1[pos_target]];



				/*Delete*/
				if (pos_target< l1 && pos_seq>1) {
					/*Initialise with match*/
					vt_d[seq][pos_seq][pos_target] = vt_m[seq][pos_seq-1][pos_target] + my_pars->ldel;
					tb_d[seq][pos_seq][pos_target] = (double) seq*10 + 1 + (pos_seq-1)/my_matrices->tb_divisor;

					/*Compare to DD*/
					if ((vt_d[seq][pos_seq-1][pos_target] + my_pars->leps) > vt_d[seq][pos_seq][pos_target]) {
						vt_d[seq][pos_seq][pos_target] = vt_d[seq][pos_seq-1][pos_target] + my_pars->leps;
						tb_d[seq][pos_seq][pos_target] = (double) seq*10 + 3 + (pos_seq-1)/my_matrices->tb_divisor;
					}
				}


				/*Get new max_r, etc.*/
				if (vt_m[seq][pos_seq][pos_target] > max_rn) {
					max_rn = vt_m[seq][pos_seq][pos_target];
					who_max_n = seq;
					state_max_n = 1;
					pos_max_n = pos_seq;
				}
				if (vt_i[seq][pos_seq][pos_target] > max_rn) {
					max_rn = vt_i[seq][pos_seq][pos_target];
					who_max_n = seq;
					state_max_n = 2;
					pos_max_n = pos_seq;
				}
			}
		}

		/*Completed position in target sequence*/
		max_r = max_rn;
		who_max = who_max_n;
		state_max = state_max_n;
		pos_max = pos_max_n;
	}
	/*Completed Viterbi matrices*/

	if (DEBUG>1) {
		print_forward_matrices(my_data, my_pars, my_matrices, target, stdout);
		print_backward_matrices(my_data, my_pars, my_matrices, target, stdout);
	}


	/*Now print maximum likelihood*/
	printf("\nMaximum Log likelihood  = %.5lf\n", max_r + my_pars->lterm);
	my_matrices->llk = max_r + my_pars->lterm;
	my_pars->combined_llk += max_r + my_pars->lterm;

	/*Now reconstruct ML path*/
	cp = 2*my_matrices->maxl;
	my_matrices->maxpath_copy[cp] = who_max;
	my_matrices->maxpath_state[cp] = state_max;
	my_matrices->maxpath_pos[cp] = pos_max;

	if (DEBUG) printf("\nEnding in state %i for sequence %i at position %i", state_max, who_max, pos_max);

	pos_target=l1;
	while(pos_target>=1) {

		if (state_max==1) {
			who_next = (int) (tb_m[who_max][pos_max][pos_target])/10;
			state_next = (int) (tb_m[who_max][pos_max][pos_target]-who_next*10);
			pos_next = (int) ((tb_m[who_max][pos_max][pos_target]-who_next*10-state_next)*my_matrices->tb_divisor+1e-6);
		}
		else if (state_max==2) {
			who_next = (int) (tb_i[who_max][pos_max][pos_target])/10;
			state_next = (int) (tb_i[who_max][pos_max][pos_target]-who_next*10);
			pos_next = (int) ((tb_i[who_max][pos_max][pos_target]-who_next*10-state_next)*my_matrices->tb_divisor+1e-6);
		}
		else if (state_max==3) {
			who_next = (int) (tb_d[who_max][pos_max][pos_target])/10;
			state_next = (int) (tb_d[who_max][pos_max][pos_target]-who_next*10);
			pos_next = (int) ((tb_d[who_max][pos_max][pos_target]-who_next*10-state_next)*my_matrices->tb_divisor+1e-6);
		}

		cp--;
		if (cp<=0) {
			printf("\n\n***Error: reconstructed path longer than maximum possible ***\n\n");
			exit(1);
		}
		my_matrices->maxpath_copy[cp] = who_next;
		my_matrices->maxpath_state[cp] = state_next;
		my_matrices->maxpath_pos[cp] = pos_next;

		who_max = who_next;
		state_max = state_next;
		pos_max = pos_next;

		if (my_matrices->maxpath_state[cp+1] != 3) pos_target--;

		if (DEBUG) printf("\nNext = state %i for sequence %i at position %i.  Pos target = %i", state_max, who_max, pos_max, pos_target);

	}

	if (DEBUG) {

		printf("\n\nMax Likelihood alignment\n\nWho\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_copy[i]);
		printf("\nState\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_state[i]);
		printf("\nPos\t");
		for (i=cp;i<=2*my_matrices->maxl;i++) printf("%i ", my_matrices->maxpath_pos[i]);
		printf("\n\n");
	}


	/*Reset copy state*/
	my_matrices->who_copy[target]=tmp_copy;

	ofp = fopen(my_pars->alignment_file, "a");
	fprintf(ofp,"\nTarget: %s\tLength: %i\tMLlk: %.3lf\n",my_data->seqs[target]->name, my_data->seqs[target]->length, my_matrices->llk);

	/*First print target sequence*/
	cp++;
	strncpy(tmp_name, my_data->seqs[target]->name, 15);
	tmp_name[15]='\0';
	fprintf(ofp,"%15s\t", tmp_name);
	for (i=cp,pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==3) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[target]->seq[pos_target], my_data->type));
			pos_target++;
		}
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do matching*/
	for (i=1;i<=15;i++) fprintf(ofp," ");
	fprintf(ofp,"\t");
	for (i=cp, pos_target=1;i<=2*my_matrices->maxl;i++) {
		if (my_matrices->maxpath_state[i]==1) {
			if (my_data->seqs[target]->seq[pos_target] == my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]]) fprintf(ofp,"|");
			else fprintf(ofp," ");
			pos_target++;
		}
		else if (my_matrices->maxpath_state[i]==2) {
			pos_target++;
			fprintf(ofp,"^");
		}
		else fprintf(ofp,"~");
	}
	fprintf(ofp,"\n");
	fflush(ofp);

	/*Now do copy tracks - switch whenever gets to new value*/
	strncpy(tmp_name, my_data->seqs[my_matrices->maxpath_copy[cp]]->name, 15);
	tmp_name[15]='\0';
	fprintf(ofp,"%15s\t",tmp_name);
	for (i=cp;i<=2*my_matrices->maxl;i++) {

		/*Check to see if need to make recombination event*/
		if (i>cp && my_matrices->maxpath_copy[i]!=my_matrices->maxpath_copy[i-1]) {
			strncpy(tmp_name, my_data->seqs[my_matrices->maxpath_copy[i]]->name, 15);
			tmp_name[15]='\0';
			fprintf(ofp,"\n%15s\t", tmp_name);
			for (j=1;j<=(i-cp);j++) fprintf(ofp," ");
		}

		if (my_matrices->maxpath_state[i]==2) fprintf(ofp,"-");
		else {
			fprintf(ofp,"%c",num2nuc(my_data->seqs[my_matrices->maxpath_copy[i]]->seq[my_matrices->maxpath_pos[i]], my_data->type));
		}
	}

	fprintf(ofp,"\n\n");



	fflush(ofp);
	fclose(ofp);

	return;


}



/*To estimate parameters by EM assuming no recombination*/

void estimate_parameters_em_no_rec(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices) {

	int iteration, i, j, target, tmp_copy;
	double current_llk=SMALL, tol=10000, sum;
	extern char *states_NT[];
	extern char *states_AA[];

	printf("\n\n*** Estimating emission/transition parameters for Pr{rec}=0 ***\n\n");

	/*First calculate insert state emission probabilities from average sequence composition*/
	for (i=0;i<my_pars->nstate;i++) my_pars->si[i]=(double) 0.01;
	for (target=1;target<=my_data->nseq; target++) if (my_data->seqs[target]->group == my_pars->target_group) {
		for (i=1;i<=my_data->seqs[target]->length;i++) my_pars->si[my_data->seqs[target]->seq[i]]++;
	}
	for (i=0,sum=0;i<my_pars->nstate;i++) sum+= (double) my_pars->si[i];
	for (i=0;i<my_pars->nstate;i++) {
		my_pars->si[i] = (double) my_pars->si[i]/sum;
		my_pars->lsi[i] = (double) log(my_pars->si[i]);
	}
	printf("\n\nInsert state emission probabilities\n\n");
	for (i=0;i<my_pars->nstate;i++) {
		printf("%c\t%.3lf\t%.3lf\n", num2nuc(i,my_pars->type), my_pars->si[i], my_pars->lsi[i]);
	}

	/*Set up matrices to store expected transitions and emissions*/
	my_matrices->expected_transitions = (double **) malloc((size_t) 3*sizeof(double *));
	my_matrices->expected_emissions = (double **) malloc((size_t) (my_pars->nstate+1)*sizeof(double *));
	for (i=0;i<3;i++) my_matrices->expected_transitions[i] = (double *) malloc((size_t) 3*sizeof(double));
	for (i=0;i<=my_pars->nstate;i++) my_matrices->expected_emissions[i] = (double *) malloc((size_t) (my_pars->nstate+1)*sizeof(double));

	for (iteration=1; iteration<=(MAX_IT_EM) && tol>LLK_TOL; iteration++) {

		my_pars->combined_llk=0.0;

		/*Initialise expected matrices*/
		for (i=0;i<3;i++) for (j=0;j<3;j++) my_matrices->expected_transitions[i][j]=(double) 0.01;

		/*For emission probs, have a pseudo-count of 1 for each element to prevent non-zero values*/
		for (i=0;i<my_pars->nstate;i++) for (j=0;j<my_pars->nstate;j++) my_matrices->expected_emissions[i][j]=(double) 0.01;

		for (target=1; target<=my_data->nseq; target++) {
			if (my_data->seqs[target]->group == my_pars->target_group) {

				tmp_copy = my_matrices->who_copy[target];
				my_matrices->who_copy[target]=0;

				kalign_fb_no_rec(my_data, my_pars, my_matrices, target);
				calculate_expected_values(my_data, my_pars, my_matrices, target);

				/*Reset copy state*/
				my_matrices->who_copy[target]=tmp_copy;

			}
		}

		/*Now calculate new values of emission parameters*/
		for (i=0;i<my_pars->nstate;i++) {
			for (j=0, sum=0.0;j<my_pars->nstate;j++) sum+=my_matrices->expected_emissions[i][j];
			for (j=0;j<my_pars->nstate;j++) {
				my_pars->sm[i][j] = (double) my_matrices->expected_emissions[i][j]/sum;
				my_pars->lsm[i][j] = (double) log(my_pars->sm[i][j]);
			}
		}

		/*Update transition parameters: NB this effectively assumes the rec and term->0*/
		my_pars->del = (double) (my_matrices->expected_transitions[0][1]+my_matrices->expected_transitions[0][2])/(2*(my_matrices->expected_transitions[0][1]+my_matrices->expected_transitions[0][2]+my_matrices->expected_transitions[0][0]));
		my_pars->eps = (double) (my_matrices->expected_transitions[1][1]+my_matrices->expected_transitions[2][2])/(my_matrices->expected_transitions[1][1]+my_matrices->expected_transitions[2][2]+my_matrices->expected_transitions[1][0]+my_matrices->expected_transitions[2][0]);
		my_pars->gm = 1-my_pars->eps;
		my_pars->dm = 1-my_pars->eps;
		my_pars->mm = 1-2*my_pars->del;
		my_pars->piM = (1-my_pars->eps)/(1-my_pars->eps+2*my_pars->del);

		/*NB seems to be better to have piM as a separate parameter*/
		my_pars->piM = (double) DEFAULT_PIM;
		my_pars->piI = 1-my_pars->piM;

		my_pars->ldel = log(my_pars->del);
		my_pars->leps = log(my_pars->eps);
		my_pars->lgm = log(my_pars->gm);
		my_pars->ldm = log(my_pars->dm);
		my_pars->lmm = log(my_pars->mm);
		my_pars->lpiM = log(my_pars->piM);
		my_pars->lpiI = log(my_pars->piI);
		
		printf("\nIteration %4i: log likelihood = %.3lf", iteration, my_pars->combined_llk);
		printf("\nParameters: gap insertion = %.5lf, gap extension = %.5lf\n\n", my_pars->del, my_pars->eps);


		tol = my_pars->combined_llk-current_llk;
		current_llk = my_pars->combined_llk;

	}

	printf("\n\nEM estimation completed\n\n");
	printf("Parameters: gap insertion = %.5lf, gap extension = %.5lf\n", my_pars->del, my_pars->eps);
	printf("\n\nEmission parameters\n\n\t");
	for (i=0;i<my_pars->nstate;i++) printf("%c\t",num2nuc(i, my_pars->type));
	for (i=0;i<my_pars->nstate;i++) {
		printf("\n%c\t",num2nuc(i, my_pars->type));
		for (j=0;j<my_pars->nstate;j++) printf("%.3lf\t", my_pars->sm[i][j]);
	}

	printf("\n\n");


	/*Free up memory*/
	for (i=0;i<3;i++) free(my_matrices->expected_transitions[i]);
	for (i=0;i<=my_pars->nstate;i++) free(my_matrices->expected_emissions[i]);
	free(my_matrices->expected_transitions);
	free(my_matrices->expected_emissions);

	/*At end need to reinstate values of rec and term*/
	my_pars->mm = 1-2*my_pars->del-my_pars->term-my_pars->rho;
	my_pars->gm = 1-my_pars->eps - my_pars->rho - my_pars->term;
	if (my_pars->mm <=0 || my_pars->gm <=0) {
		printf("\n\n***Error: after EM, some parameters have been set <=0 ***\n\n");
		exit(1);
	}
	my_pars->lmm = log(my_pars->mm);
	my_pars->lgm = log(my_pars->gm);


	return;

}


/*Function to calculate FB matrices with no recombination*/

void kalign_fb_no_rec(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	int pos_target, pos_seq, seq, l1, l2, *s1, *s2;
	double ***forward_m, ***backward_m, ***forward_i, ***backward_i, ***forward_d, ***backward_d;
	double llk_f, llk_b=SMALL, max_r;

	

	l1=my_data->seqs[target]->length;
	s1=my_data->seqs[target]->seq;
	my_pars->sizeL=0.0;
	for (seq=1;seq<=my_data->nseq;seq++) if (seq!=target) my_pars->sizeL+=(double) my_data->seqs[seq]->length;
	my_pars->lsizeL = (double) log(my_pars->sizeL);

	if (DEBUG) printf("\nSizeL = %.0lf",my_pars->sizeL);

	/*Make foward and backward matrices*/
	forward_m = (double ***) my_matrices->m1_m;
	forward_i = (double ***) my_matrices->m1_i;
	forward_d = (double ***) my_matrices->m1_d;
	backward_m = (double ***) my_matrices->m2_m;
	backward_i = (double ***) my_matrices->m2_i;
	backward_d = (double ***) my_matrices->m2_d;

	/*Set everything to small*/
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		l2=my_data->seqs[seq]->length;
		for (pos_target=0;pos_target<=(l1+1);pos_target++) 
			for (pos_seq=0;pos_seq<=(l2+1); pos_seq++) {
				forward_m[seq][pos_seq][pos_target]=SMALL;
				forward_i[seq][pos_seq][pos_target]=SMALL;
				forward_d[seq][pos_seq][pos_target]=SMALL;
				backward_m[seq][pos_seq][pos_target]=SMALL;
				backward_i[seq][pos_seq][pos_target]=SMALL;
				backward_d[seq][pos_seq][pos_target]=SMALL;
			}
	}

	/*Initialise forward and backward matrices*/
			
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2=my_data->seqs[seq]->length;
		s2=my_data->seqs[seq]->seq;

		for (pos_seq=1;pos_seq<=l2;pos_seq++){

			/*Forward Ms*/
			forward_m[seq][pos_seq][1] = (double) my_pars->lpiM-my_pars->lsizeL;
			forward_m[seq][pos_seq][1] += my_pars->lsm[s1[1]][s2[pos_seq]];
			forward_i[seq][pos_seq][1] = my_pars->lpiI-my_pars->lsizeL;
			forward_i[seq][pos_seq][1] += my_pars->lsi[s1[1]];
			if (pos_seq>1) {
				forward_d[seq][pos_seq][1] = (double) forward_m[seq][pos_seq-1][1] + log(my_pars->del + my_pars->eps*exp(forward_d[seq][pos_seq-1][1]-forward_m[seq][pos_seq-1][1]));
			}

			/*Backward Ms*/
			backward_m[seq][pos_seq][l1] = my_pars->lterm;
			backward_i[seq][pos_seq][l1] = my_pars->lterm;
		}
	}



	/*Now loop forward over positions*/

	for (pos_target=2;pos_target<=l1;pos_target++) {

		/*Match, Insert and Delete Matrices*/
		for (seq=1; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			for (pos_seq=1;pos_seq<=l2;pos_seq++) {

				/*Match*/
				max_r = forward_m[seq][pos_seq-1][pos_target-1];
				if (forward_i[seq][pos_seq-1][pos_target-1]>max_r) max_r = forward_i[seq][pos_seq-1][pos_target-1];
				if (forward_d[seq][pos_seq-1][pos_target-1]>max_r) max_r = forward_d[seq][pos_seq-1][pos_target-1];
				forward_m[seq][pos_seq][pos_target] =  exp(forward_m[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->mm;
				forward_m[seq][pos_seq][pos_target] += exp(forward_i[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->gm;
				forward_m[seq][pos_seq][pos_target] += exp(forward_d[seq][pos_seq-1][pos_target-1]-max_r)*my_pars->dm;
				forward_m[seq][pos_seq][pos_target] = (double) log(forward_m[seq][pos_seq][pos_target])+max_r;
				forward_m[seq][pos_seq][pos_target] += my_pars->lsm[s1[pos_target]][s2[pos_seq]];

				/*Insert*/
				max_r = forward_m[seq][pos_seq][pos_target-1];
				if (forward_i[seq][pos_seq][pos_target-1]>max_r) max_r = forward_i[seq][pos_seq][pos_target-1];
				forward_i[seq][pos_seq][pos_target] =  exp(forward_m[seq][pos_seq][pos_target-1]-max_r)*my_pars->del;
				forward_i[seq][pos_seq][pos_target] += exp(forward_i[seq][pos_seq][pos_target-1]-max_r)*my_pars->eps;
				forward_i[seq][pos_seq][pos_target] = (double) log(forward_i[seq][pos_seq][pos_target])+max_r;
				forward_i[seq][pos_seq][pos_target] += my_pars->lsi[s1[pos_target]];

				/*Delete: NB cannot terminate from a delete state*/
				if (pos_target<l1 && pos_seq>1) {
					max_r = forward_d[seq][pos_seq-1][pos_target];
					if (forward_m[seq][pos_seq-1][pos_target]>max_r) max_r = forward_m[seq][pos_seq-1][pos_target];
					forward_d[seq][pos_seq][pos_target] =  exp(forward_d[seq][pos_seq-1][pos_target]-max_r)*my_pars->eps;
					forward_d[seq][pos_seq][pos_target] += exp(forward_m[seq][pos_seq-1][pos_target]-max_r)*my_pars->del;
					forward_d[seq][pos_seq][pos_target] = (double) log(forward_d[seq][pos_seq][pos_target]) + max_r;
				}

			}
		}
	}
	/*Completed forward matrices*/

	if (DEBUG>1) print_forward_matrices(my_data, my_pars, my_matrices, target, stdout);


	/*Now calculate likelihood*/
	for (seq=1, max_r = SMALL;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2 = my_data->seqs[seq]->length;
			for (pos_seq=1;pos_seq<=l2;pos_seq++){ 
				if (forward_m[seq][pos_seq][l1]>max_r) max_r = forward_m[seq][pos_seq][l1];
				if (forward_i[seq][pos_seq][l1]>max_r) max_r = forward_i[seq][pos_seq][l1];
			}
	}
	for (seq=1, llk_f=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2 = my_data->seqs[seq]->length;
			for (pos_seq=1;pos_seq<=l2;pos_seq++){ 
				llk_f += (double) exp(forward_m[seq][pos_seq][l1]-max_r);
				llk_f += (double) exp(forward_i[seq][pos_seq][l1]-max_r);
			}
	}
	llk_f = max_r + log(llk_f) +my_pars->lterm;
	printf("\nSequence %s\nLog likelihood from forward algorithm  = %.5lf\n", my_data->seqs[target]->name, llk_f);
	my_matrices->llk = llk_f;
	my_pars->combined_llk += llk_f;

	

	/*Now loop backward*/
	for (pos_target = l1-1; pos_target>0; pos_target--) {

		/*Match, Insert and Delete Matrices*/
		for (seq=1; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;

			/*Note: Have to do Delete state first and in reverse position order*/

			for (pos_seq=l2;pos_seq>0;pos_seq--) {
				
				/*Delete*/
				if (pos_seq<l2) {
					max_r = backward_d[seq][pos_seq+1][pos_target];
					if (backward_m[seq][pos_seq+1][pos_target+1]>max_r) max_r = backward_m[seq][pos_seq+1][pos_target+1];
					backward_d[seq][pos_seq][pos_target] =  exp(backward_d[seq][pos_seq+1][pos_target]-max_r)*my_pars->eps;
					backward_d[seq][pos_seq][pos_target] += exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]]*my_pars->dm;
					backward_d[seq][pos_seq][pos_target] = (double) log(backward_d[seq][pos_seq][pos_target]) + max_r;
				}

				/*Insert*/
				max_r = backward_i[seq][pos_seq][pos_target+1];
				if (backward_m[seq][pos_seq+1][pos_target+1]>max_r) max_r = backward_m[seq][pos_seq+1][pos_target+1];
				backward_i[seq][pos_seq][pos_target] =  exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*my_pars->eps*my_pars->si[s1[pos_target+1]];
				backward_i[seq][pos_seq][pos_target] += exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->gm*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]];
				backward_i[seq][pos_seq][pos_target] = (double) log(backward_i[seq][pos_seq][pos_target])+max_r;

				/*Match*/
				max_r = backward_m[seq][pos_seq+1][pos_target+1];
				if (backward_i[seq][pos_seq][pos_target+1]>max_r) max_r = backward_i[seq][pos_seq][pos_target+1];
				if (backward_d[seq][pos_seq+1][pos_target]>max_r) max_r = backward_d[seq][pos_seq+1][pos_target];
				backward_m[seq][pos_seq][pos_target] =  exp(backward_m[seq][pos_seq+1][pos_target+1]-max_r)*my_pars->mm*my_pars->sm[s1[pos_target+1]][s2[pos_seq+1]];
				backward_m[seq][pos_seq][pos_target] += exp(backward_i[seq][pos_seq][pos_target+1]-max_r)*my_pars->del*my_pars->si[s1[pos_target+1]];
				backward_m[seq][pos_seq][pos_target] += exp(backward_d[seq][pos_seq+1][pos_target]-max_r)*my_pars->del;
				backward_m[seq][pos_seq][pos_target] = (double) log(backward_m[seq][pos_seq][pos_target])+max_r;

			}
		}

	}
	/*Completed backward matrices*/

	/*Calculate backwards likelihood*/
	for (seq=1, max_r=SMALL; seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
		l2 = my_data->seqs[seq]->length;
		for (pos_seq=1;pos_seq<=l2;pos_seq++){ 
			if (backward_m[seq][pos_seq][1]>max_r) max_r = backward_m[seq][pos_seq][1];
			if (backward_i[seq][pos_seq][1]>max_r) max_r = backward_i[seq][pos_seq][1];
		}
	}
	for (seq=1, llk_f=0.0;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {
			l2 = my_data->seqs[seq]->length;
			s2 = my_data->seqs[seq]->seq;
			for (pos_seq=1;pos_seq<=l2;pos_seq++){ 
				llk_f += (double) exp(backward_m[seq][pos_seq][1]-max_r)*my_pars->sm[s1[1]][s2[pos_seq]]*my_pars->piM;
				llk_f += (double) exp(backward_i[seq][pos_seq][1]-max_r)*my_pars->si[s1[1]]*my_pars->piI;
			}
	}
	llk_f = max_r + log(llk_f) - log(my_pars->sizeL);
	printf("Log likelihood from backward algorithm = %.5lf\n\n", llk_f);



	if (DEBUG>1) print_backward_matrices(my_data, my_pars, my_matrices, target, stdout);


	if (DEBUG) print_posterior_matrices(my_data, my_pars, my_matrices, target);

	return;

}




/*Function to calculate expected transitions and emissions from FB matrices*/

void calculate_expected_values(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices, int target) {

	int seq, l1, l2, pos_seq, pos_target, *s1, *s2;
	double x;

	l1=my_data->seqs[target]->length;
	s1=my_data->seqs[target]->seq;

	/*First do emission probs*/
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2=my_data->seqs[seq]->length;
		s2=my_data->seqs[seq]->seq;

		for (pos_target=1;pos_target<=l1;pos_target++) for (pos_seq=1;pos_seq<=l2;pos_seq++){
			my_matrices->expected_emissions[s2[pos_seq]][s1[pos_target]] += exp(my_matrices->m1_m[seq][pos_seq][pos_target]+my_matrices->m2_m[seq][pos_seq][pos_target]-my_matrices->llk);
		}
	}

	/*Now do transitions*/
	for (seq=1;seq<=my_data->nseq;seq++) if (my_matrices->who_copy[seq]) {

		l2=my_data->seqs[seq]->length;
		s2=my_data->seqs[seq]->seq;

		for (pos_target=1;pos_target<=l1;pos_target++) for (pos_seq=1;pos_seq<=l2;pos_seq++) {

			/*Match-match*/
			x = (double) my_matrices->m1_m[seq][pos_seq-1][pos_target-1]+my_matrices->m2_m[seq][pos_seq][pos_target]+my_pars->lmm+my_pars->lsm[s2[pos_seq]][s1[pos_target]]-my_matrices->llk;
			my_matrices->expected_transitions[0][0] += exp(x);

			/*Match-insert*/
			x = (double) my_matrices->m1_m[seq][pos_seq][pos_target-1]+my_matrices->m2_i[seq][pos_seq][pos_target]+my_pars->ldel+my_pars->lsi[s1[pos_target]]-my_matrices->llk;
			my_matrices->expected_transitions[0][1] += exp(x);

			/*Match-delete*/
			x = (double) my_matrices->m1_m[seq][pos_seq-1][pos_target]+my_matrices->m2_d[seq][pos_seq][pos_target]+my_pars->ldel-my_matrices->llk;
			my_matrices->expected_transitions[0][2] += exp(x);

			/*Insert-match*/
			x = (double) my_matrices->m1_i[seq][pos_seq-1][pos_target-1]+my_matrices->m2_m[seq][pos_seq][pos_target]+my_pars->lgm+my_pars->lsm[s2[pos_seq]][s1[pos_target]]-my_matrices->llk;
			my_matrices->expected_transitions[1][0] += exp(x);

			/*Insert-insert*/
			x = (double) my_matrices->m1_i[seq][pos_seq][pos_target-1]+my_matrices->m2_i[seq][pos_seq][pos_target]+my_pars->leps+my_pars->lsi[s1[pos_target]]-my_matrices->llk;
			my_matrices->expected_transitions[1][1] += exp(x);

			/*Delete-match*/
			x = (double) my_matrices->m1_d[seq][pos_seq-1][pos_target]+my_matrices->m2_m[seq][pos_seq][pos_target]+my_pars->ldm+my_pars->lsm[s2[pos_seq]][s1[pos_target]]-my_matrices->llk;
			my_matrices->expected_transitions[2][1] += exp(x);

			/*Delete-delete*/
			x = (double) my_matrices->m1_d[seq][pos_seq-1][pos_target]+my_matrices->m2_d[seq][pos_seq][pos_target]+my_pars->leps-my_matrices->llk;
			my_matrices->expected_transitions[2][1] += exp(x);

		}

	}

	return;

}



/*To read parameters from input file*/

void read_params_from_file(struct pars *my_pars, FILE *ifp){

	int i, j;
	char *c, *line, **end_ptr;
	double x;

	printf("\n\n*** Reading model parameters from file ***\n\n");

	line = (char *) malloc((size_t) MAXLINE*sizeof(char));
	end_ptr = (char **) malloc((size_t) sizeof(char *));
	end_ptr[0] = (char *) malloc((size_t) 5*sizeof(char));
	end_ptr[0] = strcpy(end_ptr[0], "\n");

	while(!feof(ifp)) {

		line=fgets(line, MAXLINE, ifp);
		if(!line) break;

		if (c=strstr(line, "Gap initiation")) {
			c = strchr(c, ':');
			if (!c) {
				printf("\n\n***Error: incorrect input format in parameter file ***\n\n");
				exit(1);
			}
			c++;
			my_pars->del = (double) strtod(c, end_ptr);
			if (my_pars->del<=MIN_PROB) my_pars->del = MIN_PROB;

		}
		else if (c=strstr(line, "Gap extension")) {
			c = strchr(c, ':');
			if (!c) {
				printf("\n\n***Error: incorrect input format in parameter file ***\n\n");
				exit(1);
			}
			c++;
			my_pars->eps = (double) strtod(c, end_ptr);
			if (my_pars->emiss<=MIN_PROB) my_pars->eps = MIN_PROB;

		}
		else if (c=strstr(line, "Termination")) {
			c = strchr(c, ':');
			if (!c) {
				printf("\n\n***Error: incorrect input format in parameter file ***\n\n");
				exit(1);
			}
			c++;
			my_pars->term = (double) strtod(c, end_ptr);
			if (my_pars->term<=MIN_PROB) my_pars->term = MIN_PROB;

		}
		else if (c=strstr(line, "Recombination")) {
			c = strchr(c, ':');
			if (!c) {
				printf("\n\n***Error: incorrect input format in parameter file ***\n\n");
				exit(1);
			}
			c++;
			my_pars->rho = (double) strtod(c, end_ptr);
			if (my_pars->rho<=MIN_PROB) my_pars->rho = MIN_PROB;
;
		}

		else if(c=strstr(line, "State frequencies")) {
			if (!c) {
				printf("\n\n***Error: incorrect input format in parameter file ***\n\n");
				exit(1);
			}
			c++;
			for (i=0;i<my_pars->nstate;i++) {
				fgetc(ifp);
				fscanf(ifp, "%lf", &x);
				my_pars->si[i]=x;
				if (my_pars->si[i]<=MIN_PROB) my_pars->si[i]=(double) MIN_PROB;
				my_pars->lsi[i]=log(my_pars->si[i]);
				while(fgetc(ifp)!='\n') if (feof(ifp)) break;
			}
		}

		else if(c=strstr(line, "Emission probabilities")) {
			c = strchr(c, ':');
			if (!c) {
				printf("\n\n***Error: incorrect input format in parameter file ***\n\n");
				exit(1);
			}
			c++;
			/*Get past row with state names*/
			while(fgetc(ifp)!='\n');
			
			for (i=0;i<my_pars->nstate;i++) {
				fgetc(ifp);
				for (j=0;j<my_pars->nstate;j++) {
					fscanf(ifp,"%lf", &x);
					my_pars->sm[i][j]=x;
					if (my_pars->sm[i][j]<=MIN_PROB) my_pars->sm[i][j] = (double) MIN_PROB;
					my_pars->lsm[i][j]=(double) log(my_pars->sm[i][j]);
				}
				while(fgetc(ifp)!='\n') if (feof(ifp)) break;
			}
		}

	}

	free(line);

	return;

}




void calculate_llk_over_rho_grid(struct data *my_data, struct pars *my_pars, struct matrices *my_matrices) {

	int i, which_max=0, target;
	double *rho_list, **llk_grid, mllk=0;
	FILE *ofp;

	rho_list = (double *) malloc((size_t) (my_pars->npts+1)*sizeof(double));
	llk_grid = (double **) malloc((size_t) (my_data->nseq+1)*sizeof(double *));
	for (i=1;i<=my_data->nseq;i++) llk_grid[i] = (double *) malloc((size_t) (my_pars->npts+1)*sizeof(double));
	
	/*Set up rho values for grid*/
	rho_list[1] = (double) my_pars->rho_low;
	rho_list[my_pars->npts] = (double) my_pars->rho_high;
	if (rho_list[1]<1e-32) rho_list[1]=1e-32;
	if (rho_list[my_pars->npts]>0.1) rho_list[my_pars->npts]=0.1;
	if (my_pars->log_scale) {
		rho_list[1]=log(rho_list[1]);
		rho_list[my_pars->npts] = log(rho_list[my_pars->npts]);
	}
	for (i=2;i<my_pars->npts;i++) {
		rho_list[i] = (double) rho_list[1]+(rho_list[my_pars->npts]-rho_list[1])*(i-1)/(my_pars->npts-1);
	}
	if(my_pars->log_scale) {
		for (i=1;i<=my_pars->npts;i++) rho_list[i]=exp(rho_list[i]);
	}
	printf("\n\nEsimating Llk over grid for rho\nRho\tLlk\n");
	for (i=1;i<=my_pars->npts;i++) {

		/*Initialise llk*/
		my_pars->combined_llk=0.0;

		/*Set parameters*/
		my_pars->rho = rho_list[i];
		my_pars->lrho=log(my_pars->rho);
		my_pars->mm = (double) 1-2*my_pars->del-my_pars->rho-my_pars->term;
		my_pars->gm = (double) 1-my_pars->eps-my_pars->rho-my_pars->term;
		my_pars->lmm=log(my_pars->mm);
		my_pars->lgm=log(my_pars->gm);

		/*Iterate through target sequences*/
		for (target=1; target<=my_data->nseq; target++) {
			if (my_data->seqs[target]->group == my_pars->target_group) {
				kalign_fb(my_data, my_pars, my_matrices, target, 1);
				llk_grid[target][i] = (double) my_matrices->llk;
			}
		}
		printf("%.5lf\t%.3lf\n", rho_list[i], my_pars->combined_llk);
		if (which_max==0 || my_pars->combined_llk>mllk) {
			mllk = my_pars->combined_llk;
			which_max=i;
		}
	}

	printf("\n\nMax llk at rho = %.5lf\n\n", rho_list[which_max]);

	my_pars->rho = rho_list[which_max];
	my_pars->lrho=log(my_pars->rho);
	my_pars->mm = (double) 1-2*my_pars->del-my_pars->rho-my_pars->term;
	my_pars->gm = (double) 1-my_pars->eps-my_pars->rho-my_pars->term;
	my_pars->lmm=log(my_pars->mm);
	my_pars->lgm=log(my_pars->gm);

	/*Print results to file*/
	ofp = fopen(my_pars->llk_grid_file, "w");
	fprintf(ofp,"Sequence\t");
	for (i=1;i<=my_pars->npts;i++) fprintf(ofp,"%.5lf\t", rho_list[i]);
	fprintf(ofp,"\n");
	for (target=1; target<=my_data->nseq; target++) {
		if (my_data->seqs[target]->group == my_pars->target_group) {
			fprintf(ofp,"%s\t", my_data->seqs[target]->name);
			for (i=1;i<=my_pars->npts;i++) fprintf(ofp,"%.3lf\t", llk_grid[target][i]);
			fprintf(ofp,"\n");
		}
	}
	fclose(ofp);


	free(rho_list);
	for (i=1;i<=my_data->nseq;i++) free(llk_grid[i]);
	free(llk_grid);

	return;
}

void print_help(FILE *ofp) {

	fprintf(ofp,"\n\nHelp for mosaic: command line flags\n\n");

	fprintf(ofp,"\nRequired parameters\n\n");

	fprintf(ofp,"-nt\t\t\tIndicates sequences are nucleotides\n");
	fprintf(ofp,"-aa\t\t\tIndicates sequences are amino acids\n");
	fprintf(ofp,"-seq <file>\t\tName of file containing sequences in FASTA format\n");

	fprintf(ofp,"\n\nOptional parameters\n\n");
	fprintf(ofp,"-params <file>\tName of file containing input parameters\n");
	fprintf(ofp,"-rec <double>\t\tRecombination probability (default = %.4lf)\n", (double) DEFAULT_REC);
	fprintf(ofp,"-del <double>\t\tRate of moving to delete state (default = %.3lf)\n", (double) DEFAULT_DEL);
	fprintf(ofp,"-eps <double>\t\tRate of extending in delete/insert state (default = %.3lf)\n", (double) DEFAULT_EPS);
	fprintf(ofp,"-term <double>\t\tRate of moving to termination state (default = %.3lf)\n", (double) DEFAULT_TERM);
	fprintf(ofp,"-pmatch <double>\tProbability of match (default = %.2lf)\n", (double) DEFAULT_PMATCH);
	fprintf(ofp,"-tag <string>\t\tTag to attach to output files\n");
	fprintf(ofp,"-group <int> <string1> <string2> ... \tNumber of groups to split data into and identifiers\n");
	fprintf(ofp,"-target <string>\tGroup to analyse as target sequences\n");
	fprintf(ofp,"-psum           \tPrint out summed posteriors for match states in each sequence\n");
	fprintf(ofp,"-estimate\t\tEstimate parameters by EM with Pr{rec}=0\n");
	fprintf(ofp,"-grid <double> <double> <int> <flag> \tGrid for estimation of rho: values required are start, end, #pts and flag indicating whether on log scale\n");

	fprintf(ofp,"\n\nAdditional options\n\n");
/*	fprintf(ofp,"-e1\t\t\tFlag that sets emission probabilities to 1\n");*/
	fprintf(ofp,"-ma\t\t\tReturns maximum accuracy alignment rather than Viterbi\n");

	fprintf(ofp,"\n\n");


	exit(0);


}



/*End of file*/
