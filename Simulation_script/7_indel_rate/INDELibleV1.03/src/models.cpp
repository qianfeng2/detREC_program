/* 
   INDELible V1.03
    "A comprehensive and flexible simulator of molecular sequence evolution"
    Copyright (C) 2009 William Fletcher

    If using this software please cite the following reference:

    Fletcher, W. and Yang, Z. 2009. 
	"INDELible: a flexible simulator of biological sequence evolution." 
	Mol. Biol. and Evol. (in press). 
 
    If you need to contact me with any queries, questions, problems or (god-forbid) bugs
    then please go to http://abacus.gene.ucl.ac.uk/software/indelible/bug.php

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    Please visit http://www.gnu.org/licenses/ for a full copy of the license.
*/

             

#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>  //???
#include "paml.cpp"


using namespace std;

#pragma warning(disable:4786)	


vector<vector<double> > totalusermodels;

class indelmodel
{
	// all indel models conform to these parameters.
	// not all are needed for each model
	// but they must all be present for the pointer to the functions below to work.
public:
	double meansize;
	int type;
	int r;
	double q;
	double a;
	double b;
	int M;

	vector<double> usermodel;

	indelmodel()
	{
		//these settings chosen to cause a crash when indel model parameters not specified.
		type=r=M=-1;
		meansize=q=a=b=-1;		
	}
};


// functions below are all of form int double vector<double>&
// this is so generic pointer to these functions works

int choosenewzipf(int M, double a, vector<double> &z) {int u; do {u=newZipf(a);} while(u>M || u<1); return u;} //do {u=newZipf(a);} while(u>M); return u;}
int chooseoldzipf(int M, double a, vector<double> &z) {int u; do {u=oldZipf(a);} while(u>M || u<1); return u;} //do {u=oldZipf(a);} while(u>M); return u;}

	int choosenewNB(int r, double q, vector<double> &z) {return randnegbin(r,q);}
	int chooseoldNB(int r, double q, vector<double> &z) {return oldrandnegbin(r,q);}
		
	int userrandomsize(int x, double y, vector<double> &usermodel)
	{
		double myrand=mtrand1();

		for(int i=0; i<usermodel.size(); i++) if(myrand<usermodel.at(i)) return i+1;

		return usermodel.size();
	}


	int f(){return 1;}
	int f2(double x, int y) {return int(x);}



int type;	// 1 for NUCLEOTIDE, 2 for AMINOACID, 4 for CODON


extern bool controldebug;

extern void controlerrorprint2(string,string,string,string,string);


const char GeneticCodeTable[24][65]={

	// the following codes are as listed on genbank on October 2008

	// 0 - NOTHING
	"",
	// 1  - The Standard Code
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    // 2  - The Vertebrate Mitochondrial Code
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    // 3  - The Yeast Mitochondrial Code
    "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    // 4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 5  - The Invertebrate Mitochondrial Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
	// 6  - The Ciliate, Dasycladacean and Hexamita Nuclear Code
	"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 7  - deleted
	"",
	// 8  - deleted
	"",
	// 9  - The Echinoderm and Flatworm Mitochondrial Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	// 10 - The Euplotid Nuclear Code
	"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 11 - The Bacterial and Plant Plastid Code
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 12 - The Alternative Yeast Nuclear Code
	"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 13 - The Ascidian Mitochondrial Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
	// 14 - The Alternative Flatworm Mitochondrial Code
	"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	// 15 - The Blepharisma Nuclear Code
	"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 16 - The Chlorophycean Mitochondrial Code
	"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 17 - deleted
	"",
	// 18 - deleted
	"",
	// 19 - deleted
	"",
	// 20 - deleted
	"",
	// 21 - The Trematode Mitochondrial Code 
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	// 22 - The Scenedesmus obliquus mitochondrial Code 
	"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	// 23 - The Thraustochytrium Mitochondrial Code
	"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
};

vector<int> allowedcodes(int gencode)
{
	vector<int> allowedlist;  

	if(gencode==2 || gencode==6 || gencode ==14 || gencode==22 || gencode==23) allowedlist.push_back(gencode);

	else if(gencode==15 || gencode==16) {allowedlist.push_back(15); allowedlist.push_back(16);}

	else if(gencode==1 || gencode==11 || gencode==12) {allowedlist.push_back(1);allowedlist.push_back(11);allowedlist.push_back(12);}

	else {allowedlist.push_back(3);allowedlist.push_back(4);allowedlist.push_back(5);allowedlist.push_back(9);
	allowedlist.push_back(10);allowedlist.push_back(13);allowedlist.push_back(21);}

	return allowedlist;
}

vector<int> getstops(int geneticcode)
{
	// finds the stops in genetic codes listed above

	vector<int> stops; 

	for(int i=0; i<64; i++) if(GeneticCodeTable[geneticcode][i]=='*') stops.push_back(i);

	return stops;
}

void enforcestops(int geneticcode, vector<double> &basefreqs)
{
	// makes base frequencies of stop codons equal to zero

	vector<int> stops=getstops(geneticcode);

	for(int i=0; i<stops.size(); i++) basefreqs.at(stops.at(i))=0;
}

int querystops(int geneticcode, vector<double> &basefreqs)
{
	// returns non zero base frequencies at stop codons

	vector<int> stops=getstops(geneticcode);

	for(int i=0; i<stops.size(); i++) if(basefreqs.at(stops.at(i))!=0) return stops.at(i);

	return -1;
}
//////////////////////////////////////////////////////////////////////////////////////////



class model
{
	// this BIG class contains all details of substitution/indel model and all construtors etc
public:

	// parameters for indel model
	double q1,  q2,  Hx0,  s,  Himax;
	
	double insD;
	int insI;
	vector<double> insV;

	double delD;
	int delI;
	vector<double> delV;

	double delmeansize;		// used in formula for total deletion rate
	
	// pointer for indel size generation
	int (*delrandomsize)(int, double, vector<double>&) ;  
	int (*insrandomsize)(int, double, vector<double>&) ;   
	
	// test function pointers
	int (*pf)();
	int (*pf2)(double, int);

	// model name
	string name;

	int modelnumber;	//?
	int geneticcode;

	bool copiedmodel;			// whether the model is copied - defunct?
	bool codonratestrue;		// whether model has different codon position specific relative rates
	
//	vector<double> insertrates;
//	vector<double> deleterates;

	double insertrate;			// relative instantaneous rate of insertion
	double deleterate;			// relative instantaneous rate of deletion
	double indelrate;			// relative instantaneous rate of insertion and deletion
	double subrate;				// relative instantaneous rate of substitution (=1)

	int modelpos;				// position in totalmodels
	int type;					// nucleotide=1, aminoacid=2, codon=3
	int error;				
	int rootlength;				// ?
	double alpha;				// alpha for gamma models
	double pinv;			
	int ngamcat;				// number of gamma categories for discrete gamma

	double codonrates[3];		// relative substitution rates for codon positions

	bool continuousgamma;		// whether this model has continuous gamma rate heterogeneity

	int numberofsiteclasses;	// number of classes in codon sites model, or number of gamma categories too..

	int medianORmean;			// 1 = use medians, 0 = use means, to represent categories in discrete gamma rate variation

	vector<double> cumfreqs;	// cumulative frequencies for discrete gamma, or codon sites models	
	vector<double> Rrates;		// relative rates for discrete gamma categories
	vector<double> myomegas;	// different omegas for different site classes.

	vector<double> rootbasefreqs;	// base frequencies used in root sequence creation if model at root
	vector<double> basefreqs;		// stationary frequencies of model
	vector<double> insertfreqs;		// frequencies used to generateinserted sequences
	vector<double> myrates;			// vector of substitution rates "away" from a given state, used in method 2 jump chain

	vector<vector<double> > myratesvec;		//vector of different "myrates" vectors for different site classes in codon site models

	vector<int> ratevec;			// not used? ---> done in main skeleton file now.

	vector<vector<double> > Jvec;	// transition matrix of the jump chain
	vector<vector<double> > Qvec;	// transition probabilities

	vector<double> scalefactors;	// used when scaling different Qvec from different site classes

	vector<vector<vector<double> > > Jvecs;	///////   collection of Jvecs from different site classes
	vector<vector<vector<double> > > Qvecs;	///////   collection of Qvecs from different site classes

	model(int mymodelpos, int &mytype, string &myname, int &mymodelnumber, int &mygeneticcode,
	bool &mycopiedmodel, double &myinsertrate, double &mydeleterate, double &myalpha, double &mypinv, 
	int &myngamcat, double mycodonrates[], vector<double> &mybasefreqs, vector<double> &myrootbasefreqs, 
	vector<double> &myinsertfreqs, vector<double> &myparams, vector<double> &aamodel, indelmodel &insertmodel,
	indelmodel &deletemodel)
	{
		
	
		// set deletion model
		delmeansize=deletemodel.meansize;

		if(deletemodel.type == 0 || deletemodel.type == 3 || deletemodel.type == 12 ) 
		{
			delV=deletemodel.usermodel;
			delrandomsize=&userrandomsize;
		}
		else if(deletemodel.type == 2 || deletemodel.type == 13) 
		{
			//Zipfian model

			delI=deletemodel.M;   
			delD=deletemodel.a;
			
			if(deletemodel.type==2) 
			{
				double v=1, q=delD;
	
				q1=1-q; q2=1/q1;
	
				Hx0 = H(0.5,q1,q2,v) - exp( log(v) * (-q) );
	
				s = 1 - H1( H(1.5,q1,q2,v) -exp( log(v+1) * (-q) )  ,q1,q2,v );

				Himax = H( imax + 0.5 ,q1,q2,v);

				delrandomsize=&choosenewzipf;
			}
			else					delrandomsize=&chooseoldzipf;
		}
		else if(deletemodel.type == 1 || deletemodel.type == 11) 
		{
			// Negative Binomial or Geometric

			delI=deletemodel.r; 
			delD=deletemodel.q;

			if(deletemodel.type==1) delrandomsize=&choosenewNB;
			else					delrandomsize=&chooseoldNB;
		}


		// set insertion model
		if(insertmodel.type == 0 || insertmodel.type == 3 || insertmodel.type == 12 ) 
		{
			insV=insertmodel.usermodel;

			insrandomsize=&userrandomsize;
		}
		else if(insertmodel.type == 2 || insertmodel.type == 13) 
		{
			//Zipfian model

			insI=insertmodel.M; 
			insD=insertmodel.a;

			if(insertmodel.type==2) 
			{
				double v=1, q=insD;
	
				q1=1-q; q2=1/q1;
	
				Hx0 = H(0.5,q1,q2,v) - exp( log(v) * (-q) );
	
				s = 1 - H1( H(1.5,q1,q2,v) -exp( log(v+1) * (-q) )  ,q1,q2,v );

				Himax = H( imax + 0.5 ,q1,q2,v);

				insrandomsize=&choosenewzipf;
			}
			else					insrandomsize=&chooseoldzipf;
		}
		else if(insertmodel.type == 1 || insertmodel.type == 11) 
		{
			// Negative Binomial or Geometric

			insI=insertmodel.r; 
			insD=insertmodel.q;

			if(insertmodel.type==1) insrandomsize=&choosenewNB;
			else					insrandomsize=&chooseoldNB;
		}



						

		continuousgamma=false;
		numberofsiteclasses=1;

		insertrate	=myinsertrate;
		deleterate	=mydeleterate;
		subrate=1;
		//indelrate=insertrate+deleterate;
		//subrate		=1-insertrate-deleterate;

		modelpos=mymodelpos;

		error=1;
		type=mytype;
		name=myname;
		modelnumber=mymodelnumber;
		geneticcode=mygeneticcode;
		copiedmodel=mycopiedmodel;

		alpha=myalpha;
		pinv=mypinv;
		ngamcat=myngamcat;

		medianORmean=0;
	
		if(type!=3)
		{
			// set up discrete gamma and/or pinv --> actual rates are picked in SetSiteRates in main skeleton file.
			if(alpha>0) 
			{
				if(ngamcat==0) {Rrates.push_back(1); cumfreqs.push_back(1);  continuousgamma=true;  }
				else
				{
					DiscreteGamma(cumfreqs,Rrates,pinv, alpha,alpha,ngamcat,medianORmean);
				}
			}
			else
			{
				if(pinv>0) {Rrates.push_back(0); Rrates.push_back(1/(1-pinv));  cumfreqs.push_back(pinv); cumfreqs.push_back(1); }
				else {Rrates.push_back(1); cumfreqs.push_back(1); }
			} 

			numberofsiteclasses=Rrates.size();
		}

		// set up codon position specific relative substitution rates.
		// has no effect if gamma model used - that supercedes this setting.

		if( alpha==0 && (mycodonrates[0]!=1 || mycodonrates[1]!=1 || mycodonrates[2]!=1) ) codonratestrue=true; else codonratestrue=false;
		codonrates[0]=mycodonrates[0];
		codonrates[1]=mycodonrates[1];
		codonrates[2]=mycodonrates[2];

		if(type==1 && !copiedmodel)
		{
			// for DNA
			if(mybasefreqs.empty()) 
			{
				// if base frequencies are empty
				if(modelnumber%2==1 && modelnumber<16 ) //&& !copiedmodel)  
				{
					//when they shouldn't be give a warning
					stringstream med; med<<modelnumber; string modelnumberX=med.str();
					controlerrorprint2("[MODEL]",name,"basefreqs","A model with unequal base frequencies was chosen: model "+modelnumberX+"\nbut you have not specified base frequencies.  They will be equal.",""); 
				}
				// make equal frequencies either way
				makeequalfreqs(type,basefreqs);
			}
			else
			{
				// if base frequencies are given
				if(modelnumber%2==0 && modelnumber<16 ) // && !copiedmodel)  
				{
					// on a model that wants equal frequencies
					// force frequencies to be equal and give warning
					makeequalfreqs(type,basefreqs);
					stringstream med; med<<modelnumber; string modelnumberX=med.str();
					controlerrorprint2("[MODEL]",name,"basefreqs","A model with equal base frequencies was chosen: model "+modelnumberX+"\nbut you have specified base frequencies.  They will be set equal.",""); 
				}
				// otherwise set frequencies as given
				else basefreqs=mybasefreqs;  
			}
			
			
		}

		if(type==2 && !mybasefreqs.empty()) basefreqs=mybasefreqs;  // this will force +F models		

		
		if(type==3) 
		{
			if(mybasefreqs.empty()) makeequalfreqs(type,basefreqs);			// fequal frequencies

			else if(mybasefreqs.size()==4) basefreqs=fx(mybasefreqs,1);		// f1x4 frequencies
			
			else if(mybasefreqs.size()==12) basefreqs=fx(mybasefreqs,3);	// f3x4 frequencies
			
			else if(mybasefreqs.size()==64) 								// fcodon frequencies
			{
				basefreqs=mybasefreqs;			

				testmyfreqs(basefreqs,"[basefreq]");
			}
			else cout<<"INTERNAL ERROR 463"<<endl;

		}
			

		
		
		// these make the correct Q matrix for a given type and model number
				
		// make Qvec and Jvec for nucleotide models
		if(type==1) 
		{

		//	for(int y=0; y<Rrates.size(); y++)
		//	{ 
				Qvec=getDNA(name,myparams,basefreqs, modelnumber);   Qvecs.push_back(Qvec);
				getJvec(0, /*Rrates.at(y)*/name,myrates,Qvec,Jvec,basefreqs); 
				
				if(Rrates.size()!=0) Jvecs.assign(Rrates.size(),Jvec); else Jvecs.push_back(Jvec);

		//	}
		}

		// make Qvec and Jvec for amino acid models
		if(type==2) 
		{
		//	for(int y=0; y<Rrates.size(); y++)
		//	{ 
				Qvec=getAA( name,myparams,basefreqs, modelnumber, aamodel); Qvecs.push_back(Qvec);
				
				getJvec(0,name,myrates,Qvec,Jvec,basefreqs); 
				
				if(Rrates.size()!=0) Jvecs.assign(Rrates.size(),Jvec); else Jvecs.push_back(Jvec);
				
		//	}
		}

		if(type==3) 
		{
			
			/*
			(*) Codon models for variable dN/dS ratios among sites
				(com.nrate includes kappa & omega) (see also CDFdN_dS)

				NSsites          npara

				0  one-ratio     0:    one ratio for all sites
				1  neutral       1:    p0 (w0=0, w1=1)
				2  selection     3:    p0, p1, w2 (w0=0, w1=1)
				3  discrete      2K-1: p0,p1,..., and w0,w1,...
				4  freqs         K:    p's (w's are fixed)
				5  gamma         2:    alpha, beta
				6  2gamma        4:    p0, alpha1,beta1, alpha2=beta2
				7  beta          2:    p_beta, q_beta
				8  beta&w        4:    p0, p_beta, q_beta, w estimated
				9  beta&gamma    5:    p0, p_beta, q_beta, alpha, beta
			   10  beta&1+gamma  5:    p0, p_beta, q_beta, alpha, beta (1+gamma used)
			   11  beta&1>normal 5:    p0, p_beta, q_beta, mu, s    (normal truncated w>1)
			   12  0&2normal     5:    p0, p1, mu2, s1, s2
			   13  3normal       6:    p0, p1, mu2, s0, s1, s2
			   14  M8a:beta&w=1  3:    p0, p_beta, q_beta, w=1 fixed
			   15  M8a:beta&w>=1 4:    p0, p_beta, q_beta, w>=1 estimated


			   */		

			// CODON MODELS - numbered 0 to 15

			// only 3 (M3) and 14, 15 are in proper use.  M0-M13 can all be expressed as M3

			// 0 is the codon model of Goldman and Yang 1994 ; Goldman, N., and Z. Yang. 1994. A codon-based model of nucleotide substitution for protein-coding DNA sequences. Molecular Biology and Evolution 11:725-736.

			// 1 to 13 are the codon sites-models of :Yang, Z., Nielsen, R., Goldman, N. and Pedersen, A-M. K. (2000) Codon-Substitution Models for Heterogeneous Selection Pressure at Amino Acid Sites. Genetics 155: 431-439 (May 2000)

			// 14 and 15 are the empirical codon models ECM : Kosiol, C., Holmes, I. and Goldman, N. (2007) An Empirical Codon Model for Protein Sequence Evolution.  Molecular Biology and Evolution 24(7): 1464-1479.


			if(modelnumber==14) 
			{
				Qvec=getECMr();
				d(Qvec,scalefactors.at(0)); 
				Qvecs.push_back(Qvec); 
				getJvec(0,name,myrates,Qvec,Jvec,basefreqs); 
				Jvecs.push_back(Jvec); 
				cumfreqs.push_back(1);
			}
			else if(modelnumber==15) 
			{
				Qvec=getECMu();
				d(Qvec,scalefactors.at(0));
				Qvecs.push_back(Qvec); 
				getJvec(0,name,myrates,Qvec,Jvec,basefreqs); 
				Jvecs.push_back(Jvec); 
				cumfreqs.push_back(1); 
			}
			else
			{
				// all other models  follow same pattern of generation.
				// for each site class make the Qvec for that omega, and calculate the scale factor for that Qvec
				// then calculate overall scale factor from the other scale factors and the proportions
				// divide every Qvec for every site class by the overall scale factor.using d()
				// then calculate Jvecs for each site class' Qvec

				if(myparams.size()==0)
				{
					myparams.push_back(1); myparams.push_back(1);
					controlerrorprint2("[MODEL]",name,"","No kappa/omega have been defined so they have been set equal to 1.","");
				
				}
				double kappa=myparams.at(0);
				double omega;

				if(modelnumber==0) 
				{
					// (Goldman and Yang, 1994)  
					
					cumfreqs.push_back(1);
					omega=myparams.at(1);  myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega);
					d(Qvec,scalefactors.at(0));
					Qvecs.push_back(Qvec);

					getJvec(0,name,myrates,Qvec,Jvec,basefreqs);
					Jvecs.push_back(Jvec); 
				}
		
				else if(modelnumber==1)
				{
					//cout<<myparams.at(1)<<" 1 1 "<<endl;
					double p0=myparams.at(1), p1=1-p0;

					cumfreqs.push_back(p0);
					omega=myparams.at(2);   myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);					

					//cout<<1-myparams.at(1)<<" 1 2 "<<endl;
					cumfreqs.push_back(p1);
					omega=1;   myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);

					//double S=(p0*scalefactors.at(0))+(p1*scalefactors.at(1));

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) 
					{
						//ccout<<endl<<scalefactors.at(y1)<<"  "<<cumfreqs.at(y1)<<"  "<<scalefactors.at(y1)*cumfreqs.at(y1)<<endl;
					
						S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					}
					
					d(Qvecs.at(0),S); getJvec(S,name,myrates,Qvecs.at(0),Jvec,basefreqs); Jvecs.push_back(Jvec); 					
					d(Qvecs.at(1),S); getJvec(S,name,myrates,Qvecs.at(1),Jvec,basefreqs); Jvecs.push_back(Jvec); 					

				}
				else if(modelnumber==2)
				{
				//	cout<<myparams.at(1)<<" 1 1 "<<endl;
				//	cout<<myparams.at(2)<<" 1 2 "<<endl;
				//	cout<<myparams.at(3)<<" 1 3 "<<endl;
				//	cout<<1-myparams.at(1)-myparams.at(2)<<" 1 4 "<<endl;


					double p0=myparams.at(1), p1=myparams.at(2), p2=1-p0-p1;
					cumfreqs.push_back(p0);
					omega=myparams.at(3); //cout<<p0<<"  "<<omega<<endl; 
					myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);					

					cumfreqs.push_back(p1);
					omega=1;  //cout<<p1<<"  "<<omega<<endl;
 				    myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);					

					cumfreqs.push_back(p2);
					omega=myparams.at(4); // cout<<p2<<"  "<<omega<<endl;
				    myomegas.push_back(omega);
					Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
					
					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));

					//for(int t1=0; t1<Qvecs.size(); t1++) {getJvec(S,name,myrates,Qvecs.at(t1),Jvec,basefreqs); Jvecs.push_back(Jvec); 	}	
					d(Qvecs.at(0),S); getJvec(S,name,myrates,Qvecs.at(0),Jvec,basefreqs); Jvecs.push_back(Jvec); 
					d(Qvecs.at(1),S); getJvec(S,name,myrates,Qvecs.at(1),Jvec,basefreqs); Jvecs.push_back(Jvec); 
					d(Qvecs.at(2),S); getJvec(S,name,myrates,Qvecs.at(2),Jvec,basefreqs); Jvecs.push_back(Jvec); 

				}
				else if(modelnumber==3)
				{

					int mybit=myparams.size()/2; 
					double sum=0;
					for(int yf=1; yf<mybit; yf++)   
					{
						//cout<<yf<<" "<<"1"<<endl;
						omega=myparams.at(yf+mybit-1);   myomegas.push_back(omega);
						//cout<<yf<<" "<<"2 "<<myparams.at(yf+mybit-1)<<endl;
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
						
										
						//cout<<yf<<" "<<"3"<<endl;
						sum+=myparams.at(yf);
						//cout<<yf<<" "<<"4"<<endl;
						cumfreqs.push_back(myparams.at(yf));
						//cout<<yf<<" "<<"5 "<<myparams.at(yf)<<endl;
					}
	
					if(sum<=1)
					{
						//cout<<"BLAH"<<endl;
						omega=myparams.at(2*mybit-1);   myomegas.push_back(omega);
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
						
					 					
						cumfreqs.push_back(1-sum);
						//cout<<"BLAH "<<1-sum<<" "<<omega<<endl;
					}
					else cout<<"Error in sum of category frequencies in codon model 3"<<endl;

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					for(int yf0=0; yf0<Qvecs.size(); yf0++)   {d(Qvecs.at(yf0),S); getJvec(S,name,myrates,Qvecs.at(yf0),Jvec,basefreqs); Jvecs.push_back(Jvec); }
				}




				// model numbers 4 to 13 are not used any mor
				// N.B.  THERE IS NO VECTOR OF DOUBLES CALLED myomegas IN THIS SECTION _ JUST IN CASE OF CRASHES
				else if(modelnumber==4)
				{

					//as K-1 for 4, 2K-1 for 3, but is K and 2K because of kappa, difference of K
					// for modelnumber 4, K=5 means size is 5 (kappa, p0, p1, p2, p3) with p4=1-p1-p2-p3 etc
					double sum=0, mysize=myparams.size()-2;  // mysize is 3
					for(int yf=1; yf<mysize+2; yf++)   //from 1 in the list (p0 after kappa) to mysize+2=5 goes up to p3
					{
						//cout<<yf<<" "<<"1"<<endl;
						omega=(yf-1)/mysize;  //omega is 0, 1/3, 2/3. 1
						//cout<<yf<<" "<<"2 "<<endl;
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
						
							
						//cout<<yf<<" "<<"3"<<endl;
						sum+=myparams.at(yf);
						//cout<<yf<<" "<<"4"<<endl;
						cumfreqs.push_back(myparams.at(yf));
						//cout<<yf<<" "<<"5 "<<myparams.at(yf)<<"  "<<omega<<endl;

					}
					if(sum<=1)
					{
						//cout<<"BLAH"<<endl;
						omega=mysize; //omega is 3
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
						
									
						cumfreqs.push_back(1-sum);
						//cout<<"BLAH "<<1-sum<<" "<<omega<<endl;

					}
					else cout<<"Error in sum of category frequencies in codon model 4"<<endl;
				
					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					for(int yf0=0; yf0<Qvecs.size(); yf0++)   {d(Qvecs.at(yf0),S); getJvec(S,name,myrates,Qvecs.at(yf0),Jvec,basefreqs); Jvecs.push_back(Jvec); }

				}
				else 
				{
					if(modelnumber==12 || modelnumber==8)
					{
						if(modelnumber==12)
						{
							omega=0;
							Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
							
					 					
							cumfreqs.push_back(myparams.at(1));
							for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back((1-myparams.at(1))/double(ngamcat));
						}
						else
						{
							omega=myparams.at(4);
							Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
							
							 					
							cumfreqs.push_back(1-myparams.at(1));
							for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back(myparams.at(1)/double(ngamcat));
						}
					}
					else cumfreqs.assign(ngamcat,1/double(ngamcat));
 
					//for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back(1/ngamcat);
						
					vector<double> output;

					//double *mypar; mypar=new double[myparams.size()];
					
					double mypar[10]={0,0,0,0,0,0,0,0,0,0};

					//mypar[0]=0;
					for(int iu=1; iu<myparams.size(); iu++) {mypar[iu-1]=myparams.at(iu); }  

					// this function of Ziheng's calculates the discrete rates for different site classes from the model parameters

					DiscreteNSsites(mypar, ngamcat, modelnumber, output);

					for(int i=0; i<ngamcat; i++)
					{
						omega=output.at(i);
						Qvec=getCOD(name,basefreqs, modelnumber,kappa,omega); Qvecs.push_back(Qvec);
											
					}

					double S=0; for(int y1=0; y1<cumfreqs.size(); y1++) S+=(scalefactors.at(y1)*cumfreqs.at(y1));
					for(int yf0=0; yf0<Qvecs.size(); yf0++)   {d(Qvecs.at(yf0),S);  getJvec(S,name,myrates,Qvecs.at(yf0),Jvec,basefreqs); Jvecs.push_back(Jvec); }

				}//end of model>4 bracket

			} //end of solitary else bracket


			//make cumfreqs cumulative;
			double mysum=0;
			vector<double> blah=cumfreqs;

			numberofsiteclasses=cumfreqs.size();

			for(int gfv=1; gfv<numberofsiteclasses; gfv++) cumfreqs.at(gfv)+=cumfreqs.at(gfv-1);
				
			for(int f=0; f<numberofsiteclasses; f++) Rrates.push_back(1);

		}//end of type==3 bracket

		// the Qvec/Jvec must be generated before the root/insert freqs are set in case the base frequencies come from 
		// the model like in empirical sub models (codon/protein) or when using the UNREST model for DNA.

		if(myrootbasefreqs.empty()) rootbasefreqs=basefreqs; else {rootbasefreqs=myrootbasefreqs; if(type==3) testmyfreqs(rootbasefreqs,"[rootbasefreq]"); }
		
		if(myinsertfreqs.empty()) insertfreqs=basefreqs; else {insertfreqs=myinsertfreqs; if(type==3) testmyfreqs(insertfreqs,"[ibasefreq]"); }
	

	}
/////////////////////////////////////
void changeQandJ(int numcats)
{
	Jvecs.assign(numcats,Jvec);	
	Rrates.assign(numcats,1);
	Qvecs.assign(numcats,Qvec);
	myratesvec.assign(numcats,myrates);
	numberofsiteclasses=numcats;
}


//////////////////////////////////////////////////////////////////


private:


void testmyfreqs(vector<double> &basefreqs, string mycommand)
{
	// this will point an error if there is a non-zero base frequency 
	// for a stop codon depending on the genetic code used.
	// It will not give an error when a zero is put by accident
	// for a frequency not representing a stop codon.
								
	int wow=querystops(geneticcode, basefreqs);

	if(wow!=-1) 
	{
		stringstream fs,fg1; fs<<wow; fg1<<geneticcode;  string fv=fs.str(), fg=fg1.str(); 
		stringstream fs1; fs1<<basefreqs.at(wow);  string fv1=fs1.str(); 
		controlerrorprint2("[MODEL]", name, mycommand, "Base frequencies are incorrect. For genetic code chosen: "+fg+"\nThe "+fv+"th stationary frequency should be zero not "+fv1+"\notherwise you might have stop codons in the middle of a sequence.",""); 

		error=-1;
	}
}
/////////////////////////////////////
	void d(vector<vector<double> > &Q, double S) 
	{
		// multiply scale factor over Q matrix.

		int i,j,s; 
		if(type==2) s=20; else if(type==3) s=64; else s=4; 

		for(i=0; i<s; i++)  for(j=0; j<s; j++)  ( (Q.at(i)).at(j) )/=S;
	}
/////////////////////////////////////


vector<double> fx(vector<double> &basefreqs, int which)
{
	vector<double> newbasefreqs;

	int i,j,k,l=0,m=0;  if(which==3) {l=4; m=8;}


	for(i=0; i<4; i++) for(j=0; j<4; j++) for(k=0; k<4; k++) newbasefreqs.push_back(   (basefreqs.at(i)) * (basefreqs.at(j+l)) * (basefreqs.at(k+m))    );

	enforcestops(geneticcode,newbasefreqs);

	double tot;
	tot=0; for(i=0; i<64; i++) tot+=(newbasefreqs.at(i)); //	cout<<"TOT "<<tot<<endl;
	for(i=0; i<64; i++) (newbasefreqs.at(i))/=tot;
	tot=0; for(i=0; i<64; i++) tot+=(newbasefreqs.at(i)); //	cout<<"TOT "<<tot<<endl;


return newbasefreqs;


}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////




	void makeequalfreqs(int &type, vector<double> &tbasefreqs)
	{
		// makes equal stationary frequencies.
		tbasefreqs.clear();

			if(type==1)		for(int it=0; it<4; it++) tbasefreqs.push_back(0.25);
		else if(type==2)	for(int it=0; it<20; it++) tbasefreqs.push_back(0.05);
		else if(type==3) 
		{
			// for CODON the genetic code and the relevant stop codons must be considered.
			vector<int> stops=getstops(geneticcode);

			double y=stops.size();
			double x=1/(64-y);
			
			for(int it=0; it<64; it++) {tbasefreqs.push_back(x);} //tbasefreqs.assign(64,0); ?? 

			enforcestops(geneticcode,tbasefreqs);
		}
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////


	void getJvec(double S, string name, vector<double> &myrates, vector<vector<double> > Qvec, vector<vector<double> > &Jvec, vector<double> &basefreqs)
	{
		// makes Jump Chain Transition Matrix from Q matrix = only used by method 2

		int Qsize=Qvec.size();
		int Qrowsize=(Qvec.at(0)).size();
		int basesize=basefreqs.size();
		int i,j;

		vector<double> row;
		Jvec=Qvec;

		if(Qsize!=Qrowsize || Qsize!=basesize || Qrowsize!=basesize)
		{controlerrorprint2("[MODEL]",name,"getJvec","mis-match in Q matrix and base frequency size",""); error=-1;}

		// set rates from diagonal of matrix 
		myrates.clear();
		for (i=0; i<basesize; i++) {myrates.push_back((-1)*(Qvec.at(i)).at(i)); } //cout<<myrates.at(i)<<"\t"; } cout<<endl;

		//create Jump matrix
		for (i=0; i<basesize; i++)
		{		
			for (j=0; j<basesize; j++)
			{
				if(i!=j && myrates.at(i)!=0) (Jvec.at(i)).at(j) = (Qvec.at(i)).at(j)/myrates.at(i);
				else{(Jvec.at(i)).at(j)=0;}
			}
		}

		// turn into cumulative jump matrix
		for (i=0; i<basesize; i++) for (j=1; j<basesize; j++) (Jvec.at(i)).at(j)+=(Jvec.at(i)).at(j-1);

		// store rates for use in jump chain
		myratesvec.push_back(myrates);

	}

////////////////////////////////

vector<vector<double> > getDNA( string name, vector<double> nstnums, vector<double> &basefreqs, int mymodel)
{
	// function to make Q matrix for NUCLEOTIDE models

	vector<vector<double> > myQvec;

	double a=-1,b=-1,c=-1,d=-1,e=-1,f=-1,a1=-1,b1=-1,c1=-1,d1=-1,e1=-1,f1=-1;

	double nstsize=nstnums.size()+1;

	if(mymodel==0  || mymodel==1 ) 
	{
		/*
			0 is JC69: Jukes, T. H. and Cantor, C. (1969) Evolution of protein molecules. pp. 21-132 in Mammalian Protein Metabolism, ed. M. N. Munro. Academic Press, New York.

			1 is F81: Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum likelihood approach. J. Mol. Evol., 17:368-376.

		*/
		a=b=c=d=e=f=1;
	}
	if(mymodel==2  || mymodel==3 ) 
	{
		/*
			2 is K80: Kimura, M. (1980) A simple model for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. J. Mol. Evol. 16:111-120.

			3 is HKY: Hasegawa, M., H. Kishino, and T. Yano. (1985) Dating of the human-ape splitting by a molecular clock of mitchondrial DNA. J. Mol. Evol. 22:160-174.

		*/
		//cout<<nstsize<<endl; cout<<nstnums.at(0)<<"Q"<<endl; if(nstnums.size()>1) cout<<nstnums.at(1)<<"Q"<<endl;
		a=f=nstnums.at(0); if(nstsize==3) b=c=d=e=nstnums.at(1); else b=c=d=e=1;
	}
	if(mymodel==4  || mymodel==5 ) 
	{
		/*
			4 is Tn93ef, 5 is Tn93   a.k.a. TrN

			Tn93: Tamura, K. and Nei, M. (1993) Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. Mol. Bio. Evol. 10:512-526. 
		*/

		a=nstnums.at(0);   if(nstsize==4) b=c=d=e=nstnums.at(2); else b=c=d=e=1; f=nstnums.at(1);
	}
	if(mymodel==6  || mymodel==7 ) 
	{
		/*
			6 is K81, 7 is K81uf

			K81: Kimura, M. (1981) Estimation of evolutionary distances between homologous nucleotide sequences. Proc. Natl. Acad. Sci. USA. 78:454-458. 
		*/

		a=f=1; b=e=nstnums.at(0); c=d=nstnums.at(1); if(nstsize==4) a=f=nstnums.at(2); 
		//if(nstsize==4) {a=f=nstnums.at(0); b=e=nstnums.at(1); c=d=nstnums.at(2);} else {a=f=1; b=e=nstnums.at(0); c=d=nstnums.at(1);} 
	}
	if(mymodel==8  || mymodel==9 ) 
	{
		/*
			8 is TIMef, 9 is TIM
		*/

		a=nstnums.at(0);   b=e=nstnums.at(1);     c=d=nstnums.at(2); f=1; if(nstsize==5) f=nstnums.at(3);
	}
	if(mymodel==10 || mymodel==11) 
	{
		/*
			10 is TVMef, 11 is TVM
		*/
		
		//if(nstsize==6) {a=f=nstnums.at(0); b=nstnums.at(1);       c=nstnums.at(2);   d=nstnums.at(3); e=nstnums.at(4);} else {a=f=1; b=nstnums.at(0); c=nstnums.at(1); d=nstnums.at(2); e=nstnums.at(3);}
		a=f=1; b=nstnums.at(0); c=nstnums.at(1); d=nstnums.at(2); e=nstnums.at(3); if(nstsize==6) a=f=nstnums.at(4);
	}
	if(mymodel==12 || mymodel==13) 
	{
		/*
			12 is SYM : Zharkikh, A. (1994) Estimation of evolutionary distances between nucleotide sequences. J. Mol. Evol. 39:315-329.

			13 is GTR: 
				Lanave, C., G. Preparata, C. Saccone, and G. Serio. (1984) A new method for calculating evolutionary substitution rates. J. Mol. Evol. 20:86-93. 
				Tavare, S. (1986) Some probabilistic and statistical problems in the analysis of DNA sequences. Lec. Math. Life Sci. 17:57-86.
		*/
		a=nstnums.at(0);   b=nstnums.at(1);       c=nstnums.at(2);   d=nstnums.at(3); e=nstnums.at(4); f=1; if(nstsize==7) f=nstnums.at(5);
	}
				
	if(mymodel==14 || mymodel==15) 
	{
		/*
			14 is F84ef, 15 is F84
			
			Felsenstein, J. (1984) Distance Methods for inferring phylogenies: a justification.  Evolution. 38. 16-24.
		*/

		b=c=d=e=1; if(nstsize==3) b=c=d=e=nstnums.at(1);
		double Kappa = nstnums.at(0);

		double Pi_Y  = basefreqs.at(0)+basefreqs.at(1);
		double Pi_R  = basefreqs.at(2)+basefreqs.at(3);
						
		a=(1+Kappa/Pi_Y)*b;  
		f=(1+Kappa/Pi_R)*b;  
					
	}
					
/*	if(mymodel==16) 
	{
		//myinclude.at(10)=true;
		double pi;
		b=c=d=e=1; if(nstsize==4) b=c=d=e=nstnums.at(2);
		a=f=b*nstnums.at(1); 
		pi=nstnums.at(0);
                        
		basefreqs.at(0)=basefreqs.at(2)=(1-pi)/2;
		basefreqs.at(1)=basefreqs.at(3)=pi/2;
						
	}
					
		if(mymodel==17) 
*/
		if(mymodel==16)
		{
			/*
				16 is UNREST : Yang, Z. (1994) Estimating the pattern of nucleotide substitution. J Mol. Evol 39:105-111
			*/

			double denom,piT,piC,piA,piG;
			double TC,TA,TG,CT,CA,CG,AT,AC,AG,GT,GC,GA;
			
			TC=nstnums.at(0);
			TA=nstnums.at(1);
			TG=nstnums.at(2);
			CT=nstnums.at(3);
			CA=nstnums.at(4);
			CG=nstnums.at(5);
			AT=nstnums.at(6);
			AC=nstnums.at(7);
			AG=nstnums.at(8);
			GT=nstnums.at(9);  
			GC=nstnums.at(10); 
			GA=1;
			if(nstsize==13) GA=nstnums.at(11); //cout<<TA<<endl;
										
			denom=((AT*CA*GA)+(AT*CG*GA)+(AC*CT*GA)+(AT*CT*GA)+(AT*CA*GC)+(AC*CT*GC)+(AG*CT*GC)+(AT*CT*GC)+(AG*CA*GT)+(AT*CA*GT)+(AC*CG*GT)+(AG*CG*GT)+(AT*CG*GT)+(AC*CT*GT)+(AG*CT*GT)+(AT*CT*GT)+(AG*CA*TA)+(AC*CG*TA)+(AG*CG*TA)+(AG*CT*TA)+(AC*GA*TA)+(CA*GA*TA)+(CG*GA*TA)+(CT*GA*TA)+(AC*GC*TA)+(AG*GC*TA)+(CA*GC*TA)+(CT*GC*TA)+(AC*GT*TA)+(CA*GT*TA)+(CG*GT*TA)+(CT*GT*TA)+(AG*CA*TC)+(AC*CG*TC)+(AG*CG*TC)+(AT*CG*TC)+(AC*GA*TC)+(AT*GA*TC)+(CA*GA*TC)+(CG*GA*TC)+(AC*GC*TC)+(AG*GC*TC)+(AT*GC*TC)+(CA*GC*TC)+(AC*GT*TC)+(AG*GT*TC)+(AT*GT*TC)+(CA*GT*TC)+(AG*CA*TG)+(AT*CA*TG)+(AC*CG*TG)+(AG*CG*TG)+(AT*CG*TG)+(AC*CT*TG)+(AG*CT*TG)+(AT*CT*TG)+(AC*GA*TG)+(CA*GA*TG)+(CG*GA*TG)+(CT*GA*TG)+(AC*GC*TG)+(AG*GC*TG)+(AT*GC*TG)+(CA*GC*TG));
			piT=((AT*CA*GA)+(AT*CG*GA)+(AC*CT*GA)+(AT*CT*GA)+(AT*CA*GC)+(AC*CT*GC)+(AG*CT*GC)+(AT*CT*GC)+(AG*CA*GT)+(AT*CA*GT)+(AC*CG*GT)+(AG*CG*GT)+(AT*CG*GT)+(AC*CT*GT)+(AG*CT*GT)+(AT*CT*GT))/denom;
			piC=((AC*GA*TA)+(AC*GC*TA)+(AG*GC*TA)+(AC*GT*TA)+(AC*GA*TC)+(AT*GA*TC)+(AC*GC*TC)+(AG*GC*TC)+(AT*GC*TC)+(AC*GT*TC)+(AG*GT*TC)+(AT*GT*TC)+(AC*GA*TG)+(AC*GC*TG)+(AG*GC*TG)+(AT*GC*TG))/denom;
			piA=((CA*GA*TA)+(CG*GA*TA)+(CT*GA*TA)+(CA*GC*TA)+(CT*GC*TA)+(CA*GT*TA)+(CG*GT*TA)+(CT*GT*TA)+(CA*GA*TC)+(CG*GA*TC)+(CA*GC*TC)+(CA*GT*TC)+(CA*GA*TG)+(CG*GA*TG)+(CT*GA*TG)+(CA*GC*TG))/denom;
			piG=((AG*CA*TA)+(AC*CG*TA)+(AG*CG*TA)+(AG*CT*TA)+(AG*CA*TC)+(AC*CG*TC)+(AG*CG*TC)+(AT*CG*TC)+(AG*CA*TG)+(AT*CA*TG)+(AC*CG*TG)+(AG*CG*TG)+(AT*CG*TG)+(AC*CT*TG)+(AG*CT*TG)+(AT*CT*TG))/denom;
		
			/*
			               a = TC/piC;   b = TA/piA;   c = TG/piG;  
			a1 = CT/piT;                 d = CA/piA;   e = CG/piG;  
			b1 = AT/piT;  d1 = AC/piC;                 f = AG/piG;  
			c1 = GT/piT;  e1 = GC/piC;  f1 = GA/piA;  
			*/			
		               a = TC;   b = TA;  c = TG;  
			a1 = CT;             d = CA;  e = CG;  
			b1 = AT;  d1 = AC;            f = AG;  
			c1 = GT;  e1 = GC;  f1 = GA; 
			
			basefreqs.clear();
			basefreqs.push_back(piT);  basefreqs.push_back(piC);
			basefreqs.push_back(piA);  basefreqs.push_back(piG);
					
		
		}
					
		if(mymodel!=16) 
		{
			// make matrix symmetrical if it is not UNREST
			a1=a;b1=b;c1=c;d1=d;e1=e;f1=f;
					
			if(a==-1 || b==-1 || c==-1 || d==-1 || e==-1 || f==-1 || a1==-1 || b1==-1 || c1==-1 || d1==-1 || e1==-1 || f1==-1 )
			{	controlerrorprint2(	"[MODEL]", name, "INT:makeQmatrix","DNA substitution parameter routine error - one or more values = -1",""); error=-1;}			
		}
				


	// set entries of Q matrix.
	vector<double> row; myQvec.clear();
	row.clear(); row.push_back(0);	row.push_back(a);  row.push_back(b);  row.push_back(c); myQvec.push_back(row);
	row.clear(); row.push_back(a1); row.push_back(0);  row.push_back(d);  row.push_back(e); myQvec.push_back(row);
	row.clear(); row.push_back(b1); row.push_back(d1); row.push_back(0);  row.push_back(f); myQvec.push_back(row);
	row.clear(); row.push_back(c1); row.push_back(e1); row.push_back(f1); row.push_back(0);	myQvec.push_back(row);

	double sum = 0.0; int i,j;

	if(mymodel!=16)	// these steps not necessary for UNREST
	{
		// rescale stationary frequencies, to make certain they sum to 1.0 
		for (i=0; i<4; i++) sum += basefreqs.at(i); //cout<<basefreqs.at(i)<<" ";} cout<<endl;
		if(sum!=1) for (i=0; i<4; i++) basefreqs.at(i) /= sum;

		// multiply entries by stationary frequencies 
		for (i=0; i<4; i++) for (j=0; j<4; j++) (myQvec.at(i)).at(j) *= basefreqs.at(j); //base freqs are same along a column. j is column
	}
	
	// rescale, so branch lengths are in terms of expected number of substitutions per site 
	double scaler = 0.0;

	// calculate scale factor. note that Q_ii = 0 .
	for (i=0; i<4; i++) for (j=0; j<4; j++)  scaler += basefreqs.at(i) * (myQvec.at(i)).at(j);

	scaler = 1.0 / scaler;
	
	for (i=0; i<4; i++) for (j=0; j<4; j++) (myQvec.at(i)).at(j) *= scaler;
		
	//set diagonal of matrix
	for( i=0; i<4; i++) {sum=0; for(j=0; j<4; j++) sum+=(myQvec.at(i)).at(j); (myQvec.at(i)).at(i)=-sum; }  //cout<<-sum<<endl;}

	return myQvec;
}


////////////////////////////////////////////////////////////////////////////

vector<vector<double> > getAA( string name, vector<double> params, vector<double> &basefreqs, int modelnumber, vector<double> aamodel)
{
	// set Q matrix for AMINOACID substitution models

	int			i, j; //, k;

	double	diff;

	vector<double> myrates;

	double AAPi[20], AAmatrix[20][20], AAmatrixT[20][20], AAPiT[20];

	/* 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 */
	/* A R N D C Q E G H I L  K  M  F  P  S  T  W  Y  V  */ 

	bool rootTOmodelOVERRIDE=false;	
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(aamodel.size()!=0)
	{
		// USER MODEL

		for(i=0; i<20; i++) {AAmatrix[i][i]=0; AAPi[i]=0.05;}

		// AAmatrix[i][j]
		// 1,0
		// 2,0  2,1
		// 3,0  3,1  3,2

		int k=0;

		for(i=1; i<20; i++) for(j=0; j<i; j++) { AAmatrix[i][j]=AAmatrix[j][i]=aamodel.at(k); k++; }

		if(aamodel.size()!=190) for(i=0; i<20; i++) {basefreqs.push_back(aamodel.at(k)); k++;}

	}
	else if(modelnumber==0)
	{
		/*poisson*/

		for(i=0; i<20; i++){ AAPi[i]=0.05; for(j=0; j<20; j++){ if(i!=j) AAmatrix[i][j]=1; else AAmatrix[i][j]=0; }}  
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else
	if(modelnumber==1)
	{

	/* 
		jones or JTT 

		Jones, D.T., Taylor, W.R. & Thornton, J.M. (1992) The rapid generation of mutation data matrices from protein sequences. 
		Comput. Applic. Biosci. 8, 275- 282.
		
	*/
	AAmatrix[ 0][ 0] =   0; AAmatrix[ 0][ 1] =  58; AAmatrix[ 0][ 2] =  54; AAmatrix[ 0][ 3] =  81; AAmatrix[ 0][ 4] =  56; 
	AAmatrix[ 0][ 5] =  57; AAmatrix[ 0][ 6] = 105; AAmatrix[ 0][ 7] = 179; AAmatrix[ 0][ 8] =  27; AAmatrix[ 0][ 9] =  36; 
	AAmatrix[ 0][10] =  30; AAmatrix[ 0][11] =  35; AAmatrix[ 0][12] =  54; AAmatrix[ 0][13] =  15; AAmatrix[ 0][14] = 194; 
	AAmatrix[ 0][15] = 378; AAmatrix[ 0][16] = 475; AAmatrix[ 0][17] =   9; AAmatrix[ 0][18] =  11; AAmatrix[ 0][19] = 298; 
	AAmatrix[ 1][ 0] =  58; AAmatrix[ 1][ 1] =   0; AAmatrix[ 1][ 2] =  45; AAmatrix[ 1][ 3] =  16; AAmatrix[ 1][ 4] = 113; 
	AAmatrix[ 1][ 5] = 310; AAmatrix[ 1][ 6] =  29; AAmatrix[ 1][ 7] = 137; AAmatrix[ 1][ 8] = 328; AAmatrix[ 1][ 9] =  22; 
	AAmatrix[ 1][10] =  38; AAmatrix[ 1][11] = 646; AAmatrix[ 1][12] =  44; AAmatrix[ 1][13] =   5; AAmatrix[ 1][14] =  74; 
	AAmatrix[ 1][15] = 101; AAmatrix[ 1][16] =  64; AAmatrix[ 1][17] = 126; AAmatrix[ 1][18] =  20; AAmatrix[ 1][19] =  17; 
	AAmatrix[ 2][ 0] =  54; AAmatrix[ 2][ 1] =  45; AAmatrix[ 2][ 2] =   0; AAmatrix[ 2][ 3] = 528; AAmatrix[ 2][ 4] =  34; 
	AAmatrix[ 2][ 5] =  86; AAmatrix[ 2][ 6] =  58; AAmatrix[ 2][ 7] =  81; AAmatrix[ 2][ 8] = 391; AAmatrix[ 2][ 9] =  47; 
	AAmatrix[ 2][10] =  12; AAmatrix[ 2][11] = 263; AAmatrix[ 2][12] =  30; AAmatrix[ 2][13] =  10; AAmatrix[ 2][14] =  15; 
	AAmatrix[ 2][15] = 503; AAmatrix[ 2][16] = 232; AAmatrix[ 2][17] =   8; AAmatrix[ 2][18] =  70; AAmatrix[ 2][19] =  16; 
	AAmatrix[ 3][ 0] =  81; AAmatrix[ 3][ 1] =  16; AAmatrix[ 3][ 2] = 528; AAmatrix[ 3][ 3] =   0; AAmatrix[ 3][ 4] =  10; 
	AAmatrix[ 3][ 5] =  49; AAmatrix[ 3][ 6] = 767; AAmatrix[ 3][ 7] = 130; AAmatrix[ 3][ 8] = 112; AAmatrix[ 3][ 9] =  11; 
	AAmatrix[ 3][10] =   7; AAmatrix[ 3][11] =  26; AAmatrix[ 3][12] =  15; AAmatrix[ 3][13] =   4; AAmatrix[ 3][14] =  15; 
	AAmatrix[ 3][15] =  59; AAmatrix[ 3][16] =  38; AAmatrix[ 3][17] =   4; AAmatrix[ 3][18] =  46; AAmatrix[ 3][19] =  31; 
	AAmatrix[ 4][ 0] =  56; AAmatrix[ 4][ 1] = 113; AAmatrix[ 4][ 2] =  34; AAmatrix[ 4][ 3] =  10; AAmatrix[ 4][ 4] =   0; 
	AAmatrix[ 4][ 5] =   9; AAmatrix[ 4][ 6] =   5; AAmatrix[ 4][ 7] =  59; AAmatrix[ 4][ 8] =  69; AAmatrix[ 4][ 9] =  17; 
	AAmatrix[ 4][10] =  23; AAmatrix[ 4][11] =   7; AAmatrix[ 4][12] =  31; AAmatrix[ 4][13] =  78; AAmatrix[ 4][14] =  14; 
	AAmatrix[ 4][15] = 223; AAmatrix[ 4][16] =  42; AAmatrix[ 4][17] = 115; AAmatrix[ 4][18] = 209; AAmatrix[ 4][19] =  62; 
	AAmatrix[ 5][ 0] =  57; AAmatrix[ 5][ 1] = 310; AAmatrix[ 5][ 2] =  86; AAmatrix[ 5][ 3] =  49; AAmatrix[ 5][ 4] =   9; 
	AAmatrix[ 5][ 5] =   0; AAmatrix[ 5][ 6] = 323; AAmatrix[ 5][ 7] =  26; AAmatrix[ 5][ 8] = 597; AAmatrix[ 5][ 9] =   9; 
	AAmatrix[ 5][10] =  72; AAmatrix[ 5][11] = 292; AAmatrix[ 5][12] =  43; AAmatrix[ 5][13] =   4; AAmatrix[ 5][14] = 164; 
	AAmatrix[ 5][15] =  53; AAmatrix[ 5][16] =  51; AAmatrix[ 5][17] =  18; AAmatrix[ 5][18] =  24; AAmatrix[ 5][19] =  20; 
	AAmatrix[ 6][ 0] = 105; AAmatrix[ 6][ 1] =  29; AAmatrix[ 6][ 2] =  58; AAmatrix[ 6][ 3] = 767; AAmatrix[ 6][ 4] =   5; 
	AAmatrix[ 6][ 5] = 323; AAmatrix[ 6][ 6] =   0; AAmatrix[ 6][ 7] = 119; AAmatrix[ 6][ 8] =  26; AAmatrix[ 6][ 9] =  12; 
	AAmatrix[ 6][10] =   9; AAmatrix[ 6][11] = 181; AAmatrix[ 6][12] =  18; AAmatrix[ 6][13] =   5; AAmatrix[ 6][14] =  18; 
	AAmatrix[ 6][15] =  30; AAmatrix[ 6][16] =  32; AAmatrix[ 6][17] =  10; AAmatrix[ 6][18] =   7; AAmatrix[ 6][19] =  45; 
	AAmatrix[ 7][ 0] = 179; AAmatrix[ 7][ 1] = 137; AAmatrix[ 7][ 2] =  81; AAmatrix[ 7][ 3] = 130; AAmatrix[ 7][ 4] =  59; 
	AAmatrix[ 7][ 5] =  26; AAmatrix[ 7][ 6] = 119; AAmatrix[ 7][ 7] =   0; AAmatrix[ 7][ 8] =  23; AAmatrix[ 7][ 9] =   6; 
	AAmatrix[ 7][10] =   6; AAmatrix[ 7][11] =  27; AAmatrix[ 7][12] =  14; AAmatrix[ 7][13] =   5; AAmatrix[ 7][14] =  24; 
	AAmatrix[ 7][15] = 201; AAmatrix[ 7][16] =  33; AAmatrix[ 7][17] =  55; AAmatrix[ 7][18] =   8; AAmatrix[ 7][19] =  47; 
	AAmatrix[ 8][ 0] =  27; AAmatrix[ 8][ 1] = 328; AAmatrix[ 8][ 2] = 391; AAmatrix[ 8][ 3] = 112; AAmatrix[ 8][ 4] =  69; 
	AAmatrix[ 8][ 5] = 597; AAmatrix[ 8][ 6] =  26; AAmatrix[ 8][ 7] =  23; AAmatrix[ 8][ 8] =   0; AAmatrix[ 8][ 9] =  16; 
	AAmatrix[ 8][10] =  56; AAmatrix[ 8][11] =  45; AAmatrix[ 8][12] =  33; AAmatrix[ 8][13] =  40; AAmatrix[ 8][14] = 115; 
	AAmatrix[ 8][15] =  73; AAmatrix[ 8][16] =  46; AAmatrix[ 8][17] =   8; AAmatrix[ 8][18] = 573; AAmatrix[ 8][19] =  11; 
	AAmatrix[ 9][ 0] =  36; AAmatrix[ 9][ 1] =  22; AAmatrix[ 9][ 2] =  47; AAmatrix[ 9][ 3] =  11; AAmatrix[ 9][ 4] =  17; 
	AAmatrix[ 9][ 5] =   9; AAmatrix[ 9][ 6] =  12; AAmatrix[ 9][ 7] =   6; AAmatrix[ 9][ 8] =  16; AAmatrix[ 9][ 9] =   0; 
	AAmatrix[ 9][10] = 229; AAmatrix[ 9][11] =  21; AAmatrix[ 9][12] = 479; AAmatrix[ 9][13] =  89; AAmatrix[ 9][14] =  10; 
	AAmatrix[ 9][15] =  40; AAmatrix[ 9][16] = 245; AAmatrix[ 9][17] =   9; AAmatrix[ 9][18] =  32; AAmatrix[ 9][19] = 961; 
	AAmatrix[10][ 0] =  30; AAmatrix[10][ 1] =  38; AAmatrix[10][ 2] =  12; AAmatrix[10][ 3] =   7; AAmatrix[10][ 4] =  23; 
	AAmatrix[10][ 5] =  72; AAmatrix[10][ 6] =   9; AAmatrix[10][ 7] =   6; AAmatrix[10][ 8] =  56; AAmatrix[10][ 9] = 229; 
	AAmatrix[10][10] =   0; AAmatrix[10][11] =  14; AAmatrix[10][12] = 388; AAmatrix[10][13] = 248; AAmatrix[10][14] = 102; 
	AAmatrix[10][15] =  59; AAmatrix[10][16] =  25; AAmatrix[10][17] =  52; AAmatrix[10][18] =  24; AAmatrix[10][19] = 180; 
	AAmatrix[11][ 0] =  35; AAmatrix[11][ 1] = 646; AAmatrix[11][ 2] = 263; AAmatrix[11][ 3] =  26; AAmatrix[11][ 4] =   7; 
	AAmatrix[11][ 5] = 292; AAmatrix[11][ 6] = 181; AAmatrix[11][ 7] =  27; AAmatrix[11][ 8] =  45; AAmatrix[11][ 9] =  21; 
	AAmatrix[11][10] =  14; AAmatrix[11][11] =   0; AAmatrix[11][12] =  65; AAmatrix[11][13] =   4; AAmatrix[11][14] =  21; 
	AAmatrix[11][15] =  47; AAmatrix[11][16] = 103; AAmatrix[11][17] =  10; AAmatrix[11][18] =   8; AAmatrix[11][19] =  14; 
	AAmatrix[12][ 0] =  54; AAmatrix[12][ 1] =  44; AAmatrix[12][ 2] =  30; AAmatrix[12][ 3] =  15; AAmatrix[12][ 4] =  31; 
	AAmatrix[12][ 5] =  43; AAmatrix[12][ 6] =  18; AAmatrix[12][ 7] =  14; AAmatrix[12][ 8] =  33; AAmatrix[12][ 9] = 479; 
	AAmatrix[12][10] = 388; AAmatrix[12][11] =  65; AAmatrix[12][12] =   0; AAmatrix[12][13] =  43; AAmatrix[12][14] =  16; 
	AAmatrix[12][15] =  29; AAmatrix[12][16] = 226; AAmatrix[12][17] =  24; AAmatrix[12][18] =  18; AAmatrix[12][19] = 323; 
	AAmatrix[13][ 0] =  15; AAmatrix[13][ 1] =   5; AAmatrix[13][ 2] =  10; AAmatrix[13][ 3] =   4; AAmatrix[13][ 4] =  78; 
	AAmatrix[13][ 5] =   4; AAmatrix[13][ 6] =   5; AAmatrix[13][ 7] =   5; AAmatrix[13][ 8] =  40; AAmatrix[13][ 9] =  89; 
	AAmatrix[13][10] = 248; AAmatrix[13][11] =   4; AAmatrix[13][12] =  43; AAmatrix[13][13] =   0; AAmatrix[13][14] =  17; 
	AAmatrix[13][15] =  92; AAmatrix[13][16] =  12; AAmatrix[13][17] =  53; AAmatrix[13][18] = 536; AAmatrix[13][19] =  62; 
	AAmatrix[14][ 0] = 194; AAmatrix[14][ 1] =  74; AAmatrix[14][ 2] =  15; AAmatrix[14][ 3] =  15; AAmatrix[14][ 4] =  14; 
	AAmatrix[14][ 5] = 164; AAmatrix[14][ 6] =  18; AAmatrix[14][ 7] =  24; AAmatrix[14][ 8] = 115; AAmatrix[14][ 9] =  10; 
	AAmatrix[14][10] = 102; AAmatrix[14][11] =  21; AAmatrix[14][12] =  16; AAmatrix[14][13] =  17; AAmatrix[14][14] =   0; 
	AAmatrix[14][15] = 285; AAmatrix[14][16] = 118; AAmatrix[14][17] =   6; AAmatrix[14][18] =  10; AAmatrix[14][19] =  23; 
	AAmatrix[15][ 0] = 378; AAmatrix[15][ 1] = 101; AAmatrix[15][ 2] = 503; AAmatrix[15][ 3] =  59; AAmatrix[15][ 4] = 223; 
	AAmatrix[15][ 5] =  53; AAmatrix[15][ 6] =  30; AAmatrix[15][ 7] = 201; AAmatrix[15][ 8] =  73; AAmatrix[15][ 9] =  40; 
	AAmatrix[15][10] =  59; AAmatrix[15][11] =  47; AAmatrix[15][12] =  29; AAmatrix[15][13] =  92; AAmatrix[15][14] = 285; 
	AAmatrix[15][15] =   0; AAmatrix[15][16] = 477; AAmatrix[15][17] =  35; AAmatrix[15][18] =  63; AAmatrix[15][19] =  38; 
	AAmatrix[16][ 0] = 475; AAmatrix[16][ 1] =  64; AAmatrix[16][ 2] = 232; AAmatrix[16][ 3] =  38; AAmatrix[16][ 4] =  42; 
	AAmatrix[16][ 5] =  51; AAmatrix[16][ 6] =  32; AAmatrix[16][ 7] =  33; AAmatrix[16][ 8] =  46; AAmatrix[16][ 9] = 245; 
	AAmatrix[16][10] =  25; AAmatrix[16][11] = 103; AAmatrix[16][12] = 226; AAmatrix[16][13] =  12; AAmatrix[16][14] = 118; 
	AAmatrix[16][15] = 477; AAmatrix[16][16] =   0; AAmatrix[16][17] =  12; AAmatrix[16][18] =  21; AAmatrix[16][19] = 112; 
	AAmatrix[17][ 0] =   9; AAmatrix[17][ 1] = 126; AAmatrix[17][ 2] =   8; AAmatrix[17][ 3] =   4; AAmatrix[17][ 4] = 115; 
	AAmatrix[17][ 5] =  18; AAmatrix[17][ 6] =  10; AAmatrix[17][ 7] =  55; AAmatrix[17][ 8] =   8; AAmatrix[17][ 9] =   9; 
	AAmatrix[17][10] =  52; AAmatrix[17][11] =  10; AAmatrix[17][12] =  24; AAmatrix[17][13] =  53; AAmatrix[17][14] =   6; 
	AAmatrix[17][15] =  35; AAmatrix[17][16] =  12; AAmatrix[17][17] =   0; AAmatrix[17][18] =  71; AAmatrix[17][19] =  25; 
	AAmatrix[18][ 0] =  11; AAmatrix[18][ 1] =  20; AAmatrix[18][ 2] =  70; AAmatrix[18][ 3] =  46; AAmatrix[18][ 4] = 209; 
	AAmatrix[18][ 5] =  24; AAmatrix[18][ 6] =   7; AAmatrix[18][ 7] =   8; AAmatrix[18][ 8] = 573; AAmatrix[18][ 9] =  32; 
	AAmatrix[18][10] =  24; AAmatrix[18][11] =   8; AAmatrix[18][12] =  18; AAmatrix[18][13] = 536; AAmatrix[18][14] =  10; 
	AAmatrix[18][15] =  63; AAmatrix[18][16] =  21; AAmatrix[18][17] =  71; AAmatrix[18][18] =   0; AAmatrix[18][19] =  16; 
	AAmatrix[19][ 0] = 298; AAmatrix[19][ 1] =  17; AAmatrix[19][ 2] =  16; AAmatrix[19][ 3] =  31; AAmatrix[19][ 4] =  62; 
	AAmatrix[19][ 5] =  20; AAmatrix[19][ 6] =  45; AAmatrix[19][ 7] =  47; AAmatrix[19][ 8] =  11; AAmatrix[19][ 9] = 961; 
	AAmatrix[19][10] = 180; AAmatrix[19][11] =  14; AAmatrix[19][12] = 323; AAmatrix[19][13] =  62; AAmatrix[19][14] =  23; 
	AAmatrix[19][15] =  38; AAmatrix[19][16] = 112; AAmatrix[19][17] =  25; AAmatrix[19][18] =  16; AAmatrix[19][19] =   0; 

	AAPi[ 0] = 0.076748;
	AAPi[ 1] = 0.051691;
	AAPi[ 2] = 0.042645;
	AAPi[ 3] = 0.051544;
	AAPi[ 4] = 0.019803;
	AAPi[ 5] = 0.040752;
	AAPi[ 6] = 0.061830;
	AAPi[ 7] = 0.073152;
	AAPi[ 8] = 0.022944;
	AAPi[ 9] = 0.053761;
	AAPi[10] = 0.091904;
	AAPi[11] = 0.058676;
	AAPi[12] = 0.023826;
	AAPi[13] = 0.040126;
	AAPi[14] = 0.050901;
	AAPi[15] = 0.068765;
	AAPi[16] = 0.058565;
	AAPi[17] = 0.014261;
	AAPi[18] = 0.032102;
	AAPi[19] = 0.066005;

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==2)
	{

	/* 
		Jones DCMUT
		
		Kosiol, C., and Goldman, N. (2005) Different versions of the Dayhoff rate matrix. Molecular Biology and Evolution 22:193-199.
	*/
	AAmatrix[0][0]=0;			AAmatrix[0][1]=0.531678; 	AAmatrix[0][2]=0.557967; 	AAmatrix[0][3]=0.827445; 	AAmatrix[0][4]=0.574478; 	AAmatrix[0][5]=0.556725; 	AAmatrix[0][6]=1.066681; 	AAmatrix[0][7]=1.740159; 	AAmatrix[0][8]=0.21997; 	AAmatrix[0][9]=0.361684; 	AAmatrix[0][10]=0.310007; 	AAmatrix[0][11]=0.369437; 	AAmatrix[0][12]=0.469395; 	AAmatrix[0][13]=0.138293; 	AAmatrix[0][14]=1.959599; 	AAmatrix[0][15]=3.887095; 	AAmatrix[0][16]=4.582565; 	AAmatrix[0][17]=0.084329; 	AAmatrix[0][18]=0.139492; 	AAmatrix[0][19]=2.924161; 		
	AAmatrix[1][0]=0.531678; 	AAmatrix[1][1]=0;			AAmatrix[1][2]=0.451095; 	AAmatrix[1][3]=0.154899; 	AAmatrix[1][4]=1.019843; 	AAmatrix[1][5]=3.021995; 	AAmatrix[1][6]=0.318483; 	AAmatrix[1][7]=1.359652; 	AAmatrix[1][8]=3.210671; 	AAmatrix[1][9]=0.239195; 	AAmatrix[1][10]=0.372261; 	AAmatrix[1][11]=6.529255; 	AAmatrix[1][12]=0.431045; 	AAmatrix[1][13]=0.065314; 	AAmatrix[1][14]=0.710489; 	AAmatrix[1][15]=1.001551; 	AAmatrix[1][16]=0.650282; 	AAmatrix[1][17]=1.257961; 	AAmatrix[1][18]=0.235601; 	AAmatrix[1][19]=0.171995; 		
	AAmatrix[2][0]=0.557967; 	AAmatrix[2][1]=0.451095; 	AAmatrix[2][2]=0; 			AAmatrix[2][3]=5.54953; 	AAmatrix[2][4]=0.313311; 	AAmatrix[2][5]=0.768834; 	AAmatrix[2][6]=0.578115; 	AAmatrix[2][7]=0.773313; 	AAmatrix[2][8]=4.025778; 	AAmatrix[2][9]=0.491003; 	AAmatrix[2][10]=0.137289; 	AAmatrix[2][11]=2.529517; 	AAmatrix[2][12]=0.33072; 	AAmatrix[2][13]=0.073481; 	AAmatrix[2][14]=0.121804; 	AAmatrix[2][15]=5.057964; 	AAmatrix[2][16]=2.351311; 	AAmatrix[2][17]=0.0277; 	AAmatrix[2][18]=0.700693; 	AAmatrix[2][19]=0.164525; 		
	AAmatrix[3][0]=0.827445; 	AAmatrix[3][1]=0.154899; 	AAmatrix[3][2]=5.54953; 	AAmatrix[3][3]=0;			AAmatrix[3][4]=0.105625; 	AAmatrix[3][5]=0.521646; 	AAmatrix[3][6]=7.766557; 	AAmatrix[3][7]=1.272434; 	AAmatrix[3][8]=1.032342; 	AAmatrix[3][9]=0.115968; 	AAmatrix[3][10]=0.061486; 	AAmatrix[3][11]=0.282466; 	AAmatrix[3][12]=0.190001; 	AAmatrix[3][13]=0.032522; 	AAmatrix[3][14]=0.127164; 	AAmatrix[3][15]=0.589268; 	AAmatrix[3][16]=0.425159; 	AAmatrix[3][17]=0.057466; 	AAmatrix[3][18]=0.453952; 	AAmatrix[3][19]=0.315261; 		
	AAmatrix[4][0]=0.574478; 	AAmatrix[4][1]=1.019843; 	AAmatrix[4][2]=0.313311; 	AAmatrix[4][3]=0.105625; 	AAmatrix[4][4]=0; 			AAmatrix[4][5]=0.091304; 	AAmatrix[4][6]=0.053907; 	AAmatrix[4][7]=0.546389; 	AAmatrix[4][8]=0.724998; 	AAmatrix[4][9]=0.150559; 	AAmatrix[4][10]=0.164593; 	AAmatrix[4][11]=0.049009; 	AAmatrix[4][12]=0.409202; 	AAmatrix[4][13]=0.678335; 	AAmatrix[4][14]=0.123653; 	AAmatrix[4][15]=2.155331; 	AAmatrix[4][16]=0.469823; 	AAmatrix[4][17]=1.104181; 	AAmatrix[4][18]=2.114852; 	AAmatrix[4][19]=0.621323; 		
	AAmatrix[5][0]=0.556725; 	AAmatrix[5][1]=3.021995; 	AAmatrix[5][2]=0.768834; 	AAmatrix[5][3]=0.521646; 	AAmatrix[5][4]=0.091304; 	AAmatrix[5][5]=0;			AAmatrix[5][6]=3.417706; 	AAmatrix[5][7]=0.231294; 	AAmatrix[5][8]=5.68408; 	AAmatrix[5][9]=0.07827; 	AAmatrix[5][10]=0.709004; 	AAmatrix[5][11]=2.966732; 	AAmatrix[5][12]=0.456901; 	AAmatrix[5][13]=0.045683; 	AAmatrix[5][14]=1.608126; 	AAmatrix[5][15]=0.548807; 	AAmatrix[5][16]=0.523825; 	AAmatrix[5][17]=0.172206; 	AAmatrix[5][18]=0.254745; 	AAmatrix[5][19]=0.179771; 		
	AAmatrix[6][0]=1.066681; 	AAmatrix[6][1]=0.318483; 	AAmatrix[6][2]=0.578115; 	AAmatrix[6][3]=7.766557; 	AAmatrix[6][4]=0.053907; 	AAmatrix[6][5]=3.417706; 	AAmatrix[6][6]=0;			AAmatrix[6][7]=1.115632; 	AAmatrix[6][8]=0.243768; 	AAmatrix[6][9]=0.111773; 	AAmatrix[6][10]=0.097485; 	AAmatrix[6][11]=1.731684; 	AAmatrix[6][12]=0.175084; 	AAmatrix[6][13]=0.043829; 	AAmatrix[6][14]=0.191994; 	AAmatrix[6][15]=0.312449; 	AAmatrix[6][16]=0.331584; 	AAmatrix[6][17]=0.114381; 	AAmatrix[6][18]=0.063452; 	AAmatrix[6][19]=0.465271; 		
	AAmatrix[7][0]=1.740159; 	AAmatrix[7][1]=1.359652; 	AAmatrix[7][2]=0.773313; 	AAmatrix[7][3]=1.272434; 	AAmatrix[7][4]=0.546389; 	AAmatrix[7][5]=0.231294; 	AAmatrix[7][6]=1.115632; 	AAmatrix[7][7]=0; 			AAmatrix[7][8]=0.201696; 	AAmatrix[7][9]=0.053769; 	AAmatrix[7][10]=0.069492; 	AAmatrix[7][11]=0.26984; 	AAmatrix[7][12]=0.130379; 	AAmatrix[7][13]=0.050212; 	AAmatrix[7][14]=0.208081; 	AAmatrix[7][15]=1.874296; 	AAmatrix[7][16]=0.316862; 	AAmatrix[7][17]=0.54418; 	AAmatrix[7][18]=0.0525; 	AAmatrix[7][19]=0.47014; 		
	AAmatrix[8][0]=0.21997; 	AAmatrix[8][1]=3.210671; 	AAmatrix[8][2]=4.025778; 	AAmatrix[8][3]=1.032342; 	AAmatrix[8][4]=0.724998; 	AAmatrix[8][5]=5.68408; 	AAmatrix[8][6]=0.243768; 	AAmatrix[8][7]=0.201696; 	AAmatrix[8][8]=0;			AAmatrix[8][9]=0.181788; 	AAmatrix[8][10]=0.540571; 	AAmatrix[8][11]=0.525096; 	AAmatrix[8][12]=0.32966; 	AAmatrix[8][13]=0.453428; 	AAmatrix[8][14]=1.141961; 	AAmatrix[8][15]=0.743458; 	AAmatrix[8][16]=0.477355; 	AAmatrix[8][17]=0.128193; 	AAmatrix[8][18]=5.8484; 	AAmatrix[8][19]=0.121827; 		
	AAmatrix[9][0]=0.361684; 	AAmatrix[9][1]=0.239195; 	AAmatrix[9][2]=0.491003; 	AAmatrix[9][3]=0.115968; 	AAmatrix[9][4]=0.150559; 	AAmatrix[9][5]=0.07827; 	AAmatrix[9][6]=0.111773; 	AAmatrix[9][7]=0.053769; 	AAmatrix[9][8]=0.181788; 	AAmatrix[9][9]=0;			AAmatrix[9][10]=2.335139; 	AAmatrix[9][11]=0.202562; 	AAmatrix[9][12]=4.831666; 	AAmatrix[9][13]=0.77709; 	AAmatrix[9][14]=0.09858; 	AAmatrix[9][15]=0.405119; 	AAmatrix[9][16]=2.553806; 	AAmatrix[9][17]=0.13451; 	AAmatrix[9][18]=0.303445; 	AAmatrix[9][19]=9.533943; 		
	AAmatrix[10][0]=0.310007; 	AAmatrix[10][1]=0.372261; 	AAmatrix[10][2]=0.137289; 	AAmatrix[10][3]=0.061486; 	AAmatrix[10][4]=0.164593; 	AAmatrix[10][5]=0.709004; 	AAmatrix[10][6]=0.097485; 	AAmatrix[10][7]=0.069492; 	AAmatrix[10][8]=0.540571; 	AAmatrix[10][9]=2.335139; 	AAmatrix[10][10]=0; 		AAmatrix[10][11]=0.146481; 	AAmatrix[10][12]=3.856906; 	AAmatrix[10][13]=2.500294; 	AAmatrix[10][14]=1.060504; 	AAmatrix[10][15]=0.592511; 	AAmatrix[10][16]=0.272514; 	AAmatrix[10][17]=0.530324; 	AAmatrix[10][18]=0.241094; 	AAmatrix[10][19]=1.761439; 		
	AAmatrix[11][0]=0.369437; 	AAmatrix[11][1]=6.529255; 	AAmatrix[11][2]=2.529517; 	AAmatrix[11][3]=0.282466; 	AAmatrix[11][4]=0.049009; 	AAmatrix[11][5]=2.966732; 	AAmatrix[11][6]=1.731684; 	AAmatrix[11][7]=0.26984; 	AAmatrix[11][8]=0.525096; 	AAmatrix[11][9]=0.202562; 	AAmatrix[11][10]=0.146481; 	AAmatrix[11][11]=0;			AAmatrix[11][12]=0.624581; 	AAmatrix[11][13]=0.024521; 	AAmatrix[11][14]=0.216345; 	AAmatrix[11][15]=0.474478; 	AAmatrix[11][16]=0.965641; 	AAmatrix[11][17]=0.089134; 	AAmatrix[11][18]=0.087904; 	AAmatrix[11][19]=0.124066; 		
	AAmatrix[12][0]=0.469395; 	AAmatrix[12][1]=0.431045; 	AAmatrix[12][2]=0.33072; 	AAmatrix[12][3]=0.190001; 	AAmatrix[12][4]=0.409202; 	AAmatrix[12][5]=0.456901; 	AAmatrix[12][6]=0.175084; 	AAmatrix[12][7]=0.130379; 	AAmatrix[12][8]=0.32966; 	AAmatrix[12][9]=4.831666; 	AAmatrix[12][10]=3.856906; 	AAmatrix[12][11]=0.624581; 	AAmatrix[12][12]=0;			AAmatrix[12][13]=0.436181; 	AAmatrix[12][14]=0.164215; 	AAmatrix[12][15]=0.285564; 	AAmatrix[12][16]=2.114728; 	AAmatrix[12][17]=0.201334; 	AAmatrix[12][18]=0.18987; 	AAmatrix[12][19]=3.038533; 		
	AAmatrix[13][0]=0.138293; 	AAmatrix[13][1]=0.065314; 	AAmatrix[13][2]=0.073481; 	AAmatrix[13][3]=0.032522; 	AAmatrix[13][4]=0.678335; 	AAmatrix[13][5]=0.045683; 	AAmatrix[13][6]=0.043829; 	AAmatrix[13][7]=0.050212; 	AAmatrix[13][8]=0.453428; 	AAmatrix[13][9]=0.77709; 	AAmatrix[13][10]=2.500294; 	AAmatrix[13][11]=0.024521; 	AAmatrix[13][12]=0.436181; 	AAmatrix[13][13]=0;			AAmatrix[13][14]=0.148483; 	AAmatrix[13][15]=0.943971; 	AAmatrix[13][16]=0.138904; 	AAmatrix[13][17]=0.537922; 	AAmatrix[13][18]=5.484236; 	AAmatrix[13][19]=0.593478; 		
	AAmatrix[14][0]=1.959599; 	AAmatrix[14][1]=0.710489; 	AAmatrix[14][2]=0.121804; 	AAmatrix[14][3]=0.127164; 	AAmatrix[14][4]=0.123653; 	AAmatrix[14][5]=1.608126; 	AAmatrix[14][6]=0.191994; 	AAmatrix[14][7]=0.208081; 	AAmatrix[14][8]=1.141961; 	AAmatrix[14][9]=0.09858; 	AAmatrix[14][10]=1.060504; 	AAmatrix[14][11]=0.216345; 	AAmatrix[14][12]=0.164215; 	AAmatrix[14][13]=0.148483; 	AAmatrix[14][14]=0;			AAmatrix[14][15]=2.788406; 	AAmatrix[14][16]=1.176961; 	AAmatrix[14][17]=0.069965; 	AAmatrix[14][18]=0.11385; 	AAmatrix[14][19]=0.211561; 		
	AAmatrix[15][0]=3.887095; 	AAmatrix[15][1]=1.001551; 	AAmatrix[15][2]=5.057964; 	AAmatrix[15][3]=0.589268; 	AAmatrix[15][4]=2.155331; 	AAmatrix[15][5]=0.548807; 	AAmatrix[15][6]=0.312449; 	AAmatrix[15][7]=1.874296; 	AAmatrix[15][8]=0.743458; 	AAmatrix[15][9]=0.405119; 	AAmatrix[15][10]=0.592511; 	AAmatrix[15][11]=0.474478; 	AAmatrix[15][12]=0.285564; 	AAmatrix[15][13]=0.943971; 	AAmatrix[15][14]=2.788406; 	AAmatrix[15][15]=0;			AAmatrix[15][16]=4.777647; 	AAmatrix[15][17]=0.310927; 	AAmatrix[15][18]=0.628608; 	AAmatrix[15][19]=0.408532; 		
	AAmatrix[16][0]=4.582565; 	AAmatrix[16][1]=0.650282; 	AAmatrix[16][2]=2.351311; 	AAmatrix[16][3]=0.425159; 	AAmatrix[16][4]=0.469823; 	AAmatrix[16][5]=0.523825; 	AAmatrix[16][6]=0.331584; 	AAmatrix[16][7]=0.316862; 	AAmatrix[16][8]=0.477355; 	AAmatrix[16][9]=2.553806; 	AAmatrix[16][10]=0.272514; 	AAmatrix[16][11]=0.965641; 	AAmatrix[16][12]=2.114728; 	AAmatrix[16][13]=0.138904; 	AAmatrix[16][14]=1.176961; 	AAmatrix[16][15]=4.777647; 	AAmatrix[16][16]=0;			AAmatrix[16][17]=0.080556; 	AAmatrix[16][18]=0.201094; 	AAmatrix[16][19]=1.14398; 		
	AAmatrix[17][0]=0.084329; 	AAmatrix[17][1]=1.257961; 	AAmatrix[17][2]=0.0277; 	AAmatrix[17][3]=0.057466; 	AAmatrix[17][4]=1.104181; 	AAmatrix[17][5]=0.172206; 	AAmatrix[17][6]=0.114381; 	AAmatrix[17][7]=0.54418; 	AAmatrix[17][8]=0.128193; 	AAmatrix[17][9]=0.13451; 	AAmatrix[17][10]=0.530324; 	AAmatrix[17][11]=0.089134; 	AAmatrix[17][12]=0.201334; 	AAmatrix[17][13]=0.537922; 	AAmatrix[17][14]=0.069965; 	AAmatrix[17][15]=0.310927; 	AAmatrix[17][16]=0.080556; 	AAmatrix[17][17]=0;			AAmatrix[17][18]=0.747889; 	AAmatrix[17][19]=0.239697; 		
	AAmatrix[18][0]=0.139492; 	AAmatrix[18][1]=0.235601; 	AAmatrix[18][2]=0.700693; 	AAmatrix[18][3]=0.453952; 	AAmatrix[18][4]=2.114852; 	AAmatrix[18][5]=0.254745; 	AAmatrix[18][6]=0.063452; 	AAmatrix[18][7]=0.0525; 	AAmatrix[18][8]=5.8484; 	AAmatrix[18][9]=0.303445; 	AAmatrix[18][10]=0.241094; 	AAmatrix[18][11]=0.087904; 	AAmatrix[18][12]=0.18987; 	AAmatrix[18][13]=5.484236; 	AAmatrix[18][14]=0.11385; 	AAmatrix[18][15]=0.628608; 	AAmatrix[18][16]=0.201094; 	AAmatrix[18][17]=0.747889; 	AAmatrix[18][18]=0;			AAmatrix[18][19]=0.165473; 		
	AAmatrix[19][0]=2.924161; 	AAmatrix[19][1]=0.171995; 	AAmatrix[19][2]=0.164525; 	AAmatrix[19][3]=0.315261; 	AAmatrix[19][4]=0.621323; 	AAmatrix[19][5]=0.179771; 	AAmatrix[19][6]=0.465271; 	AAmatrix[19][7]=0.47014; 	AAmatrix[19][8]=0.121827; 	AAmatrix[19][9]=9.533943; 	AAmatrix[19][10]=1.761439; 	AAmatrix[19][11]=0.124066; 	AAmatrix[19][12]=3.038533; 	AAmatrix[19][13]=0.593478; 	AAmatrix[19][14]=0.211561; 	AAmatrix[19][15]=0.408532; 	AAmatrix[19][16]=1.14398; 	AAmatrix[19][17]=0.239697; 	AAmatrix[19][18]=0.165473; 	AAmatrix[19][19]=0; 		

		
	AAPi[0]=0.076862; 
	AAPi[1]=0.051057; 
	AAPi[2]=0.042546; 
	AAPi[3]=0.051269; 
	AAPi[4]=0.020279; 
	AAPi[5]=0.041061; 
	AAPi[6]=0.06182; 
	AAPi[7]=0.074714; 
	AAPi[8]=0.022983; 
	AAPi[9]=0.052569; 
	AAPi[10]=0.091111; 
	AAPi[11]=0.059498; 
	AAPi[12]=0.023414; 
	AAPi[13]=0.04053; 
	AAPi[14]=0.050532; 
	AAPi[15]=0.068225; 
	AAPi[16]=0.058518; 
	AAPi[17]=0.014336; 
	AAPi[18]=0.032303; 
	AAPi[19]=0.066374; 


	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else
	if(modelnumber==3)
	{

	/* 
		dayhoff 
		
		Dayhoff, M., Schwartz, R., & Orcutt, B. (1978) A model of evolutionary change in protein. Atlas of Protein Sequences and Structure, 5, 345--352
	*/
	AAmatrix[ 0][ 0] =   0; AAmatrix[ 0][ 1] =  27; AAmatrix[ 0][ 2] =  98; AAmatrix[ 0][ 3] = 120; AAmatrix[ 0][ 4] =  36; 
	AAmatrix[ 0][ 5] =  89; AAmatrix[ 0][ 6] = 198; AAmatrix[ 0][ 7] = 240; AAmatrix[ 0][ 8] =  23; AAmatrix[ 0][ 9] =  65; 
	AAmatrix[ 0][10] =  41; AAmatrix[ 0][11] =  26; AAmatrix[ 0][12] =  72; AAmatrix[ 0][13] =  18; AAmatrix[ 0][14] = 250; 
	AAmatrix[ 0][15] = 409; AAmatrix[ 0][16] = 371; AAmatrix[ 0][17] =   0; AAmatrix[ 0][18] =  24; AAmatrix[ 0][19] = 208; 
	AAmatrix[ 1][ 0] =  27; AAmatrix[ 1][ 1] =   0; AAmatrix[ 1][ 2] =  32; AAmatrix[ 1][ 3] =   0; AAmatrix[ 1][ 4] =  23; 
	AAmatrix[ 1][ 5] = 246; AAmatrix[ 1][ 6] =   1; AAmatrix[ 1][ 7] =   9; AAmatrix[ 1][ 8] = 240; AAmatrix[ 1][ 9] =  64; 
	AAmatrix[ 1][10] =  15; AAmatrix[ 1][11] = 464; AAmatrix[ 1][12] =  90; AAmatrix[ 1][13] =  14; AAmatrix[ 1][14] = 103; 
	AAmatrix[ 1][15] = 154; AAmatrix[ 1][16] =  26; AAmatrix[ 1][17] = 201; AAmatrix[ 1][18] =   8; AAmatrix[ 1][19] =  24; 
	AAmatrix[ 2][ 0] =  98; AAmatrix[ 2][ 1] =  32; AAmatrix[ 2][ 2] =   0; AAmatrix[ 2][ 3] = 905; AAmatrix[ 2][ 4] =   0; 
	AAmatrix[ 2][ 5] = 103; AAmatrix[ 2][ 6] = 148; AAmatrix[ 2][ 7] = 139; AAmatrix[ 2][ 8] = 535; AAmatrix[ 2][ 9] =  77; 
	AAmatrix[ 2][10] =  34; AAmatrix[ 2][11] = 318; AAmatrix[ 2][12] =   1; AAmatrix[ 2][13] =  14; AAmatrix[ 2][14] =  42; 
	AAmatrix[ 2][15] = 495; AAmatrix[ 2][16] = 229; AAmatrix[ 2][17] =  23; AAmatrix[ 2][18] =  95; AAmatrix[ 2][19] =  15; 
	AAmatrix[ 3][ 0] = 120; AAmatrix[ 3][ 1] =   0; AAmatrix[ 3][ 2] = 905; AAmatrix[ 3][ 3] =   0; AAmatrix[ 3][ 4] =   0; 
	AAmatrix[ 3][ 5] = 134; AAmatrix[ 3][ 6] = 1153; AAmatrix[ 3][ 7] = 125; AAmatrix[ 3][ 8] =  86; AAmatrix[ 3][ 9] =  24; 
	AAmatrix[ 3][10] =   0; AAmatrix[ 3][11] =  71; AAmatrix[ 3][12] =   0; AAmatrix[ 3][13] =   0; AAmatrix[ 3][14] =  13; 
	AAmatrix[ 3][15] =  95; AAmatrix[ 3][16] =  66; AAmatrix[ 3][17] =   0; AAmatrix[ 3][18] =   0; AAmatrix[ 3][19] =  18; 
	AAmatrix[ 4][ 0] =  36; AAmatrix[ 4][ 1] =  23; AAmatrix[ 4][ 2] =   0; AAmatrix[ 4][ 3] =   0; AAmatrix[ 4][ 4] =   0; 
	AAmatrix[ 4][ 5] =   0; AAmatrix[ 4][ 6] =   0; AAmatrix[ 4][ 7] =  11; AAmatrix[ 4][ 8] =  28; AAmatrix[ 4][ 9] =  44; 
	AAmatrix[ 4][10] =   0; AAmatrix[ 4][11] =   0; AAmatrix[ 4][12] =   0; AAmatrix[ 4][13] =   0; AAmatrix[ 4][14] =  19; 
	AAmatrix[ 4][15] = 161; AAmatrix[ 4][16] =  16; AAmatrix[ 4][17] =   0; AAmatrix[ 4][18] =  96; AAmatrix[ 4][19] =  49; 
	AAmatrix[ 5][ 0] =  89; AAmatrix[ 5][ 1] = 246; AAmatrix[ 5][ 2] = 103; AAmatrix[ 5][ 3] = 134; AAmatrix[ 5][ 4] =   0; 
	AAmatrix[ 5][ 5] =   0; AAmatrix[ 5][ 6] = 716; AAmatrix[ 5][ 7] =  28; AAmatrix[ 5][ 8] = 606; AAmatrix[ 5][ 9] =  18; 
	AAmatrix[ 5][10] =  73; AAmatrix[ 5][11] = 153; AAmatrix[ 5][12] = 114; AAmatrix[ 5][13] =   0; AAmatrix[ 5][14] = 153; 
	AAmatrix[ 5][15] =  56; AAmatrix[ 5][16] =  53; AAmatrix[ 5][17] =   0; AAmatrix[ 5][18] =   0; AAmatrix[ 5][19] =  35; 
	AAmatrix[ 6][ 0] = 198; AAmatrix[ 6][ 1] =   1; AAmatrix[ 6][ 2] = 148; AAmatrix[ 6][ 3] = 1153; AAmatrix[ 6][ 4] =   0; 
	AAmatrix[ 6][ 5] = 716; AAmatrix[ 6][ 6] =   0; AAmatrix[ 6][ 7] =  81; AAmatrix[ 6][ 8] =  43; AAmatrix[ 6][ 9] =  61; 
	AAmatrix[ 6][10] =  11; AAmatrix[ 6][11] =  83; AAmatrix[ 6][12] =  30; AAmatrix[ 6][13] =   0; AAmatrix[ 6][14] =  51; 
	AAmatrix[ 6][15] =  79; AAmatrix[ 6][16] =  34; AAmatrix[ 6][17] =   0; AAmatrix[ 6][18] =  22; AAmatrix[ 6][19] =  37; 
	AAmatrix[ 7][ 0] = 240; AAmatrix[ 7][ 1] =   9; AAmatrix[ 7][ 2] = 139; AAmatrix[ 7][ 3] = 125; AAmatrix[ 7][ 4] =  11; 
	AAmatrix[ 7][ 5] =  28; AAmatrix[ 7][ 6] =  81; AAmatrix[ 7][ 7] =   0; AAmatrix[ 7][ 8] =  10; AAmatrix[ 7][ 9] =   0; 
	AAmatrix[ 7][10] =   7; AAmatrix[ 7][11] =  27; AAmatrix[ 7][12] =  17; AAmatrix[ 7][13] =  15; AAmatrix[ 7][14] =  34; 
	AAmatrix[ 7][15] = 234; AAmatrix[ 7][16] =  30; AAmatrix[ 7][17] =   0; AAmatrix[ 7][18] =   0; AAmatrix[ 7][19] =  54; 
	AAmatrix[ 8][ 0] =  23; AAmatrix[ 8][ 1] = 240; AAmatrix[ 8][ 2] = 535; AAmatrix[ 8][ 3] =  86; AAmatrix[ 8][ 4] =  28; 
	AAmatrix[ 8][ 5] = 606; AAmatrix[ 8][ 6] =  43; AAmatrix[ 8][ 7] =  10; AAmatrix[ 8][ 8] =   0; AAmatrix[ 8][ 9] =   7; 
	AAmatrix[ 8][10] =  44; AAmatrix[ 8][11] =  26; AAmatrix[ 8][12] =   0; AAmatrix[ 8][13] =  48; AAmatrix[ 8][14] =  94; 
	AAmatrix[ 8][15] =  35; AAmatrix[ 8][16] =  22; AAmatrix[ 8][17] =  27; AAmatrix[ 8][18] = 127; AAmatrix[ 8][19] =  44; 
	AAmatrix[ 9][ 0] =  65; AAmatrix[ 9][ 1] =  64; AAmatrix[ 9][ 2] =  77; AAmatrix[ 9][ 3] =  24; AAmatrix[ 9][ 4] =  44; 
	AAmatrix[ 9][ 5] =  18; AAmatrix[ 9][ 6] =  61; AAmatrix[ 9][ 7] =   0; AAmatrix[ 9][ 8] =   7; AAmatrix[ 9][ 9] =   0; 
	AAmatrix[ 9][10] = 257; AAmatrix[ 9][11] =  46; AAmatrix[ 9][12] = 336; AAmatrix[ 9][13] = 196; AAmatrix[ 9][14] =  12; 
	AAmatrix[ 9][15] =  24; AAmatrix[ 9][16] = 192; AAmatrix[ 9][17] =   0; AAmatrix[ 9][18] =  37; AAmatrix[ 9][19] = 889; 
	AAmatrix[10][ 0] =  41; AAmatrix[10][ 1] =  15; AAmatrix[10][ 2] =  34; AAmatrix[10][ 3] =   0; AAmatrix[10][ 4] =   0; 
	AAmatrix[10][ 5] =  73; AAmatrix[10][ 6] =  11; AAmatrix[10][ 7] =   7; AAmatrix[10][ 8] =  44; AAmatrix[10][ 9] = 257; 
	AAmatrix[10][10] =   0; AAmatrix[10][11] =  18; AAmatrix[10][12] = 527; AAmatrix[10][13] = 157; AAmatrix[10][14] =  32; 
	AAmatrix[10][15] =  17; AAmatrix[10][16] =  33; AAmatrix[10][17] =  46; AAmatrix[10][18] =  28; AAmatrix[10][19] = 175; 
	AAmatrix[11][ 0] =  26; AAmatrix[11][ 1] = 464; AAmatrix[11][ 2] = 318; AAmatrix[11][ 3] =  71; AAmatrix[11][ 4] =   0; 
	AAmatrix[11][ 5] = 153; AAmatrix[11][ 6] =  83; AAmatrix[11][ 7] =  27; AAmatrix[11][ 8] =  26; AAmatrix[11][ 9] =  46; 
	AAmatrix[11][10] =  18; AAmatrix[11][11] =   0; AAmatrix[11][12] = 243; AAmatrix[11][13] =   0; AAmatrix[11][14] =  33; 
	AAmatrix[11][15] =  96; AAmatrix[11][16] = 136; AAmatrix[11][17] =   0; AAmatrix[11][18] =  13; AAmatrix[11][19] =  10; 
	AAmatrix[12][ 0] =  72; AAmatrix[12][ 1] =  90; AAmatrix[12][ 2] =   1; AAmatrix[12][ 3] =   0; AAmatrix[12][ 4] =   0; 
	AAmatrix[12][ 5] = 114; AAmatrix[12][ 6] =  30; AAmatrix[12][ 7] =  17; AAmatrix[12][ 8] =   0; AAmatrix[12][ 9] = 336; 
	AAmatrix[12][10] = 527; AAmatrix[12][11] = 243; AAmatrix[12][12] =   0; AAmatrix[12][13] =  92; AAmatrix[12][14] =  17; 
	AAmatrix[12][15] =  62; AAmatrix[12][16] = 104; AAmatrix[12][17] =   0; AAmatrix[12][18] =   0; AAmatrix[12][19] = 258; 
	AAmatrix[13][ 0] =  18; AAmatrix[13][ 1] =  14; AAmatrix[13][ 2] =  14; AAmatrix[13][ 3] =   0; AAmatrix[13][ 4] =   0; 
	AAmatrix[13][ 5] =   0; AAmatrix[13][ 6] =   0; AAmatrix[13][ 7] =  15; AAmatrix[13][ 8] =  48; AAmatrix[13][ 9] = 196; 
	AAmatrix[13][10] = 157; AAmatrix[13][11] =   0; AAmatrix[13][12] =  92; AAmatrix[13][13] =   0; AAmatrix[13][14] =  11; 
	AAmatrix[13][15] =  46; AAmatrix[13][16] =  13; AAmatrix[13][17] =  76; AAmatrix[13][18] = 698; AAmatrix[13][19] =  12; 
	AAmatrix[14][ 0] = 250; AAmatrix[14][ 1] = 103; AAmatrix[14][ 2] =  42; AAmatrix[14][ 3] =  13; AAmatrix[14][ 4] =  19; 
	AAmatrix[14][ 5] = 153; AAmatrix[14][ 6] =  51; AAmatrix[14][ 7] =  34; AAmatrix[14][ 8] =  94; AAmatrix[14][ 9] =  12; 
	AAmatrix[14][10] =  32; AAmatrix[14][11] =  33; AAmatrix[14][12] =  17; AAmatrix[14][13] =  11; AAmatrix[14][14] =   0; 
	AAmatrix[14][15] = 245; AAmatrix[14][16] =  78; AAmatrix[14][17] =   0; AAmatrix[14][18] =   0; AAmatrix[14][19] =  48; 
	AAmatrix[15][ 0] = 409; AAmatrix[15][ 1] = 154; AAmatrix[15][ 2] = 495; AAmatrix[15][ 3] =  95; AAmatrix[15][ 4] = 161; 
	AAmatrix[15][ 5] =  56; AAmatrix[15][ 6] =  79; AAmatrix[15][ 7] = 234; AAmatrix[15][ 8] =  35; AAmatrix[15][ 9] =  24; 
	AAmatrix[15][10] =  17; AAmatrix[15][11] =  96; AAmatrix[15][12] =  62; AAmatrix[15][13] =  46; AAmatrix[15][14] = 245; 
	AAmatrix[15][15] =   0; AAmatrix[15][16] = 550; AAmatrix[15][17] =  75; AAmatrix[15][18] =  34; AAmatrix[15][19] =  30; 
	AAmatrix[16][ 0] = 371; AAmatrix[16][ 1] =  26; AAmatrix[16][ 2] = 229; AAmatrix[16][ 3] =  66; AAmatrix[16][ 4] =  16; 
	AAmatrix[16][ 5] =  53; AAmatrix[16][ 6] =  34; AAmatrix[16][ 7] =  30; AAmatrix[16][ 8] =  22; AAmatrix[16][ 9] = 192; 
	AAmatrix[16][10] =  33; AAmatrix[16][11] = 136; AAmatrix[16][12] = 104; AAmatrix[16][13] =  13; AAmatrix[16][14] =  78; 
	AAmatrix[16][15] = 550; AAmatrix[16][16] =   0; AAmatrix[16][17] =   0; AAmatrix[16][18] =  42; AAmatrix[16][19] = 157; 
	AAmatrix[17][ 0] =   0; AAmatrix[17][ 1] = 201; AAmatrix[17][ 2] =  23; AAmatrix[17][ 3] =   0; AAmatrix[17][ 4] =   0; 
	AAmatrix[17][ 5] =   0; AAmatrix[17][ 6] =   0; AAmatrix[17][ 7] =   0; AAmatrix[17][ 8] =  27; AAmatrix[17][ 9] =   0; 
	AAmatrix[17][10] =  46; AAmatrix[17][11] =   0; AAmatrix[17][12] =   0; AAmatrix[17][13] =  76; AAmatrix[17][14] =   0; 
	AAmatrix[17][15] =  75; AAmatrix[17][16] =   0; AAmatrix[17][17] =   0; AAmatrix[17][18] =  61; AAmatrix[17][19] =   0; 
	AAmatrix[18][ 0] =  24; AAmatrix[18][ 1] =   8; AAmatrix[18][ 2] =  95; AAmatrix[18][ 3] =   0; AAmatrix[18][ 4] =  96; 
	AAmatrix[18][ 5] =   0; AAmatrix[18][ 6] =  22; AAmatrix[18][ 7] =   0; AAmatrix[18][ 8] = 127; AAmatrix[18][ 9] =  37; 
	AAmatrix[18][10] =  28; AAmatrix[18][11] =  13; AAmatrix[18][12] =   0; AAmatrix[18][13] = 698; AAmatrix[18][14] =   0; 
	AAmatrix[18][15] =  34; AAmatrix[18][16] =  42; AAmatrix[18][17] =  61; AAmatrix[18][18] =   0; AAmatrix[18][19] =  28; 
	AAmatrix[19][ 0] = 208; AAmatrix[19][ 1] =  24; AAmatrix[19][ 2] =  15; AAmatrix[19][ 3] =  18; AAmatrix[19][ 4] =  49; 
	AAmatrix[19][ 5] =  35; AAmatrix[19][ 6] =  37; AAmatrix[19][ 7] =  54; AAmatrix[19][ 8] =  44; AAmatrix[19][ 9] = 889; 
	AAmatrix[19][10] = 175; AAmatrix[19][11] =  10; AAmatrix[19][12] = 258; AAmatrix[19][13] =  12; AAmatrix[19][14] =  48; 
	AAmatrix[19][15] =  30; AAmatrix[19][16] = 157; AAmatrix[19][17] =   0; AAmatrix[19][18] =  28; AAmatrix[19][19] =   0;

	AAPi[ 0] = 0.087127;
	AAPi[ 1] = 0.040904;
	AAPi[ 2] = 0.040432;
	AAPi[ 3] = 0.046872;
	AAPi[ 4] = 0.033474;
	AAPi[ 5] = 0.038255;
	AAPi[ 6] = 0.049530;
	AAPi[ 7] = 0.088612;
	AAPi[ 8] = 0.033618;
	AAPi[ 9] = 0.036886;
	AAPi[10] = 0.085357;
	AAPi[11] = 0.080482;
	AAPi[12] = 0.014753;
	AAPi[13] = 0.039772;
	AAPi[14] = 0.050680;
	AAPi[15] = 0.069577;
	AAPi[16] = 0.058542;
	AAPi[17] = 0.010494;
	AAPi[18] = 0.029916;
	AAPi[19] = 0.064718;


	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else
	if(modelnumber==4)
	{

	/*
		DayhoffDCMUT
		
		Kosiol, C., and Goldman, N. (2005) Different versions of the Dayhoff rate matrix. Molecular Biology and Evolution 22:193-199.
	*/
	AAmatrix[0][0] = 0;			AAmatrix[0][1] = 0.267828;	AAmatrix[0][2] = 0.984474;	AAmatrix[0][3] = 1.19981;	
	AAmatrix[0][4] = 0.360016;	AAmatrix[0][5] = 0.887753;	AAmatrix[0][6] = 1.96117;	AAmatrix[0][7] = 2.38611;	
	AAmatrix[0][8] = 0.228116;	AAmatrix[0][9] = 0.653416;	AAmatrix[0][10] = 0.406431;	AAmatrix[0][11] = 0.258635;	
	AAmatrix[0][12] = 0.71784;	AAmatrix[0][13] = 0.183641;	AAmatrix[0][14] = 2.48592;	AAmatrix[0][15] = 4.05187;	
	AAmatrix[0][16] = 3.68037;	AAmatrix[0][17] = 0;		AAmatrix[0][18] = 0.244139;	AAmatrix[0][19] = 2.05956;	
	AAmatrix[1][0] = 0.267828;	AAmatrix[1][1] = 0;			AAmatrix[1][2] = 0.327059;	AAmatrix[1][3] = 0;	
	AAmatrix[1][4] = 0.232374;	AAmatrix[1][5] = 2.43994;	AAmatrix[1][6] = 0;			AAmatrix[1][7] = 0.087791;	
	AAmatrix[1][8] = 2.38315;	AAmatrix[1][9] = 0.632629;	AAmatrix[1][10] = 0.154924;	AAmatrix[1][11] = 4.61012;	
	AAmatrix[1][12] = 0.896321;	AAmatrix[1][13] = 0.136906;	AAmatrix[1][14] = 1.02831;	AAmatrix[1][15] = 1.53159;	
	AAmatrix[1][16] = 0.265745;	AAmatrix[1][17] = 2.00137;	AAmatrix[1][18] = 0.078012;	AAmatrix[1][19] = 0.240368;	
	AAmatrix[2][0] = 0.984474;	AAmatrix[2][1] = 0.327059;	AAmatrix[2][2] = 0;			AAmatrix[2][3] = 8.93151;	
	AAmatrix[2][4] = 0;			AAmatrix[2][5] = 1.02851;	AAmatrix[2][6] = 1.49341;	AAmatrix[2][7] = 1.38535;	
	AAmatrix[2][8] = 5.29002;	AAmatrix[2][9] = 0.768024;	AAmatrix[2][10] = 0.341113;	AAmatrix[2][11] = 3.14837;	
	AAmatrix[2][12] = 0;		AAmatrix[2][13] = 0.138503;	AAmatrix[2][14] = 0.419244;	AAmatrix[2][15] = 4.88589;	
	AAmatrix[2][16] = 2.2717;	AAmatrix[2][17] = 0.224968;	AAmatrix[2][18] = 0.94694;	AAmatrix[2][19] = 0.158067;	
	AAmatrix[3][0] = 1.19981;	AAmatrix[3][1] = 0;			AAmatrix[3][2] = 8.93151;	AAmatrix[3][3] = 0;	
	AAmatrix[3][4] = 0;			AAmatrix[3][5] = 1.34855;	AAmatrix[3][6] = 11.3887;	AAmatrix[3][7] = 1.24098;	
	AAmatrix[3][8] = 0.868241;	AAmatrix[3][9] = 0.239248;	AAmatrix[3][10] = 0;		AAmatrix[3][11] = 0.716913;	
	AAmatrix[3][12] = 0;		AAmatrix[3][13] = 0;		AAmatrix[3][14] = 0.13394;	AAmatrix[3][15] = 0.956097;	
	AAmatrix[3][16] = 0.66093;	AAmatrix[3][17] = 0;		AAmatrix[3][18] = 0;		AAmatrix[3][19] = 0.178316;	
	AAmatrix[4][0] = 0.360016;	AAmatrix[4][1] = 0.232374;	AAmatrix[4][2] = 0;			AAmatrix[4][3] = 0;	
	AAmatrix[4][4] = 0;			AAmatrix[4][5] = 0;			AAmatrix[4][6] = 0;			AAmatrix[4][7] = 0.107278;	
	AAmatrix[4][8] = 0.282729;	AAmatrix[4][9] = 0.438074;	AAmatrix[4][10] = 0;		AAmatrix[4][11] = 0;	
	AAmatrix[4][12] = 0;		AAmatrix[4][13] = 0;		AAmatrix[4][14] = 0.18755;	AAmatrix[4][15] = 1.59836;	
	AAmatrix[4][16] = 0.162366;	AAmatrix[4][17] = 0;		AAmatrix[4][18] = 0.953164;	AAmatrix[4][19] = 0.484678;	
	AAmatrix[5][0] = 0.887753;	AAmatrix[5][1] = 2.43994;	AAmatrix[5][2] = 1.02851;	AAmatrix[5][3] = 1.34855;	
	AAmatrix[5][4] = 0;			AAmatrix[5][5] = 0;			AAmatrix[5][6] = 7.08602;	AAmatrix[5][7] = 0.281581;	
	AAmatrix[5][8] = 6.01161;	AAmatrix[5][9] = 0.180393;	AAmatrix[5][10] = 0.730772;	AAmatrix[5][11] = 1.51908;	
	AAmatrix[5][12] = 1.1275;	AAmatrix[5][13] = 0;		AAmatrix[5][14] = 1.52619;	AAmatrix[5][15] = 0.561828;	
	AAmatrix[5][16] = 0.525651;	AAmatrix[5][17] = 0;		AAmatrix[5][18] = 0;		AAmatrix[5][19] = 0.346983;	
	AAmatrix[6][0] = 1.96117;	AAmatrix[6][1] = 0;			AAmatrix[6][2] = 1.49341;	AAmatrix[6][3] = 11.3887;	
	AAmatrix[6][4] = 0;			AAmatrix[6][5] = 7.08602;	AAmatrix[6][6] = 0;			AAmatrix[6][7] = 0.811907;	
	AAmatrix[6][8] = 0.439469;	AAmatrix[6][9] = 0.609526;	AAmatrix[6][10] = 0.11288;	AAmatrix[6][11] = 0.830078;	
	AAmatrix[6][12] = 0.304803;	AAmatrix[6][13] = 0;		AAmatrix[6][14] = 0.507003;	AAmatrix[6][15] = 0.793999;	
	AAmatrix[6][16] = 0.340156;	AAmatrix[6][17] = 0;		AAmatrix[6][18] = 0.214717;	AAmatrix[6][19] = 0.36725;	
	AAmatrix[7][0] = 2.38611;	AAmatrix[7][1] = 0.087791;	AAmatrix[7][2] = 1.38535;	AAmatrix[7][3] = 1.24098;	
	AAmatrix[7][4] = 0.107278;	AAmatrix[7][5] = 0.281581;	AAmatrix[7][6] = 0.811907;	AAmatrix[7][7] = 0;	
	AAmatrix[7][8] = 0.106802;	AAmatrix[7][9] = 0;			AAmatrix[7][10] = 0.071514;	AAmatrix[7][11] = 0.267683;	
	AAmatrix[7][12] = 0.170372;	AAmatrix[7][13] = 0.153478;	AAmatrix[7][14] = 0.347153;	AAmatrix[7][15] = 2.32224;	
	AAmatrix[7][16] = 0.306662;	AAmatrix[7][17] = 0;		AAmatrix[7][18] = 0;		AAmatrix[7][19] = 0.538165;	
	AAmatrix[8][0] = 0.228116;	AAmatrix[8][1] = 2.38315;	AAmatrix[8][2] = 5.29002;	AAmatrix[8][3] = 0.868241;	
	AAmatrix[8][4] = 0.282729;	AAmatrix[8][5] = 6.01161;	AAmatrix[8][6] = 0.439469;	AAmatrix[8][7] = 0.106802;	
	AAmatrix[8][8] = 0;			AAmatrix[8][9] = 0.076981;	AAmatrix[8][10] = 0.443504;	AAmatrix[8][11] = 0.270475;	
	AAmatrix[8][12] = 0;		AAmatrix[8][13] = 0.475927;	AAmatrix[8][14] = 0.933709;	AAmatrix[8][15] = 0.353643;	
	AAmatrix[8][16] = 0.226333;	AAmatrix[8][17] = 0.270564;	AAmatrix[8][18] = 1.2654;	AAmatrix[8][19] = 0.438715;	
	AAmatrix[9][0] = 0.653416;	AAmatrix[9][1] = 0.632629;	AAmatrix[9][2] = 0.768024;	AAmatrix[9][3] = 0.239248;	
	AAmatrix[9][4] = 0.438074;	AAmatrix[9][5] = 0.180393;	AAmatrix[9][6] = 0.609526;	AAmatrix[9][7] = 0;	
	AAmatrix[9][8] = 0.076981;	AAmatrix[9][9] = 0;			AAmatrix[9][10] = 2.55668;	AAmatrix[9][11] = 0.460857;	
	AAmatrix[9][12] = 3.33273;	AAmatrix[9][13] = 1.95195;	AAmatrix[9][14] = 0.119152;	AAmatrix[9][15] = 0.247955;	
	AAmatrix[9][16] = 1.90074;	AAmatrix[9][17] = 0;		AAmatrix[9][18] = 0.374834;	AAmatrix[9][19] = 8.81004;	
	AAmatrix[10][0] = 0.406431;	AAmatrix[10][1] = 0.154924;	AAmatrix[10][2] = 0.341113;	AAmatrix[10][3] = 0;	
	AAmatrix[10][4] = 0;		AAmatrix[10][5] = 0.730772;	AAmatrix[10][6] = 0.11288;	AAmatrix[10][7] = 0.071514;	
	AAmatrix[10][8] = 0.443504;	AAmatrix[10][9] = 2.55668;	AAmatrix[10][10] = 0;		AAmatrix[10][11] = 0.180629;	
	AAmatrix[10][12] = 5.23011;	AAmatrix[10][13] = 1.56516;	AAmatrix[10][14] = 0.316258;AAmatrix[10][15] = 0.171432;	
	AAmatrix[10][16] = 0.33109;	AAmatrix[10][17] = 0.461776;AAmatrix[10][18] = 0.286572;AAmatrix[10][19] = 1.74516;	
	AAmatrix[11][0] = 0.258635;	AAmatrix[11][1] = 4.61012;	AAmatrix[11][2] = 3.14837;	AAmatrix[11][3] = 0.716913;	
	AAmatrix[11][4] = 0;		AAmatrix[11][5] = 1.51908;	AAmatrix[11][6] = 0.830078;	AAmatrix[11][7] = 0.267683;	
	AAmatrix[11][8] = 0.270475;	AAmatrix[11][9] = 0.460857;	AAmatrix[11][10] = 0.180629;AAmatrix[11][11] = 0;	
	AAmatrix[11][12] = 2.41174;	AAmatrix[11][13] = 0;		AAmatrix[11][14] = 0.335419;AAmatrix[11][15] = 0.954557;	
	AAmatrix[11][16] = 1.3506;	AAmatrix[11][17] = 0;		AAmatrix[11][18] = 0.132142;AAmatrix[11][19] = 0.10385;	
	AAmatrix[12][0] = 0.71784;	AAmatrix[12][1] = 0.896321;	AAmatrix[12][2] = 0;		AAmatrix[12][3] = 0;	
	AAmatrix[12][4] = 0;		AAmatrix[12][5] = 1.1275;	AAmatrix[12][6] = 0.304803;	AAmatrix[12][7] = 0.170372;	
	AAmatrix[12][8] = 0;		AAmatrix[12][9] = 3.33273;	AAmatrix[12][10] = 5.23011;	AAmatrix[12][11] = 2.41174;	
	AAmatrix[12][12] = 0;		AAmatrix[12][13] = 0.92186;	AAmatrix[12][14] = 0.170205;AAmatrix[12][15] = 0.619951;	
	AAmatrix[12][16] = 1.03153;	AAmatrix[12][17] = 0;		AAmatrix[12][18] = 0;		AAmatrix[12][19] = 2.56596;	
	AAmatrix[13][0] = 0.183641;	AAmatrix[13][1] = 0.136906;	AAmatrix[13][2] = 0.138503;	AAmatrix[13][3] = 0;	
	AAmatrix[13][4] = 0;		AAmatrix[13][5] = 0;		AAmatrix[13][6] = 0;		AAmatrix[13][7] = 0.153478;	
	AAmatrix[13][8] = 0.475927;	AAmatrix[13][9] = 1.95195;	AAmatrix[13][10] = 1.56516;	AAmatrix[13][11] = 0;	
	AAmatrix[13][12] = 0.92186;	AAmatrix[13][13] = 0;		AAmatrix[13][14] = 0.110506;AAmatrix[13][15] = 0.459901;	
	AAmatrix[13][16] = 0.136655;AAmatrix[13][17] = 0.762354;AAmatrix[13][18] = 6.95263;	AAmatrix[13][19] = 0.123606;	
	AAmatrix[14][0] = 2.48592;	AAmatrix[14][1] = 1.02831;	AAmatrix[14][2] = 0.419244;	AAmatrix[14][3] = 0.13394;	
	AAmatrix[14][4] = 0.18755;	AAmatrix[14][5] = 1.52619;	AAmatrix[14][6] = 0.507003;	AAmatrix[14][7] = 0.347153;	
	AAmatrix[14][8] = 0.933709;	AAmatrix[14][9] = 0.119152;	AAmatrix[14][10] = 0.316258;AAmatrix[14][11] = 0.335419;	
	AAmatrix[14][12] = 0.170205;AAmatrix[14][13] = 0.110506;AAmatrix[14][14] = 0;		AAmatrix[14][15] = 2.4272;	
	AAmatrix[14][16] = 0.782857;AAmatrix[14][17] = 0;		AAmatrix[14][18] = 0;		AAmatrix[14][19] = 0.485026;	
	AAmatrix[15][0] = 4.05187;	AAmatrix[15][1] = 1.53159;	AAmatrix[15][2] = 4.88589;	AAmatrix[15][3] = 0.956097;	
	AAmatrix[15][4] = 1.59836;	AAmatrix[15][5] = 0.561828;	AAmatrix[15][6] = 0.793999;	AAmatrix[15][7] = 2.32224;	
	AAmatrix[15][8] = 0.353643;	AAmatrix[15][9] = 0.247955;	AAmatrix[15][10] = 0.171432;AAmatrix[15][11] = 0.954557;	
	AAmatrix[15][12] = 0.619951;AAmatrix[15][13] = 0.459901;AAmatrix[15][14] = 2.4272;	AAmatrix[15][15] = 0;	
	AAmatrix[15][16] = 5.43667;	AAmatrix[15][17] = 0.740819;AAmatrix[15][18] = 0.336289;AAmatrix[15][19] = 0.303836;	
	AAmatrix[16][0] = 3.68037;	AAmatrix[16][1] = 0.265745;	AAmatrix[16][2] = 2.2717;	AAmatrix[16][3] = 0.66093;	
	AAmatrix[16][4] = 0.162366;	AAmatrix[16][5] = 0.525651;	AAmatrix[16][6] = 0.340156;	AAmatrix[16][7] = 0.306662;	
	AAmatrix[16][8] = 0.226333;	AAmatrix[16][9] = 1.90074;	AAmatrix[16][10] = 0.33109;	AAmatrix[16][11] = 1.3506;	
	AAmatrix[16][12] = 1.03153;	AAmatrix[16][13] = 0.136655;AAmatrix[16][14] = 0.782857;AAmatrix[16][15] = 5.43667;	
	AAmatrix[16][16] = 0;		AAmatrix[16][17] = 0;		AAmatrix[16][18] = 0.417839;AAmatrix[16][19] = 1.562;	
	AAmatrix[17][0] = 0;		AAmatrix[17][1] = 2.00137;	AAmatrix[17][2] = 0.224968;	AAmatrix[17][3] = 0;	
	AAmatrix[17][4] = 0;		AAmatrix[17][5] = 0;		AAmatrix[17][6] = 0;		AAmatrix[17][7] = 0;	
	AAmatrix[17][8] = 0.270564;	AAmatrix[17][9] = 0;		AAmatrix[17][10] = 0.461776;AAmatrix[17][11] = 0;	
	AAmatrix[17][12] = 0;		AAmatrix[17][13] = 0.762354;AAmatrix[17][14] = 0;		AAmatrix[17][15] = 0.740819;	
	AAmatrix[17][16] = 0;		AAmatrix[17][17] = 0;		AAmatrix[17][18] = 0.60807;	AAmatrix[17][19] = 0;	
	AAmatrix[18][0] = 0.244139;	AAmatrix[18][1] = 0.078012;	AAmatrix[18][2] = 0.94694;	AAmatrix[18][3] = 0;	
	AAmatrix[18][4] = 0.953164;	AAmatrix[18][5] = 0;		AAmatrix[18][6] = 0.214717;	AAmatrix[18][7] = 0;	
	AAmatrix[18][8] = 1.2654;	AAmatrix[18][9] = 0.374834;	AAmatrix[18][10] = 0.286572;AAmatrix[18][11] = 0.132142;	
	AAmatrix[18][12] = 0;		AAmatrix[18][13] = 6.95263;	AAmatrix[18][14] = 0;		AAmatrix[18][15] = 0.336289;	
	AAmatrix[18][16] = 0.417839;AAmatrix[18][17] = 0.60807;	AAmatrix[18][18] = 0;		AAmatrix[18][19] = 0.279379;	
	AAmatrix[19][0] = 2.05956;	AAmatrix[19][1] = 0.240368;	AAmatrix[19][2] = 0.158067;	AAmatrix[19][3] = 0.178316;	
	AAmatrix[19][4] = 0.484678;	AAmatrix[19][5] = 0.346983;	AAmatrix[19][6] = 0.36725;	AAmatrix[19][7] = 0.538165;	
	AAmatrix[19][8] = 0.438715;	AAmatrix[19][9] = 8.81004;	AAmatrix[19][10] = 1.74516;	AAmatrix[19][11] = 0.10385;	
	AAmatrix[19][12] = 2.56596;	AAmatrix[19][13] = 0.123606;AAmatrix[19][14] = 0.485026;AAmatrix[19][15] = 0.303836;	
	AAmatrix[19][16] = 1.562;	AAmatrix[19][17] = 0;		AAmatrix[19][18] = 0.279379;AAmatrix[19][19] = 0;	

	AAPi[0] = 0.087127;
	AAPi[1] = 0.040904;
	AAPi[2] = 0.040432;
	AAPi[3] = 0.046872;
	AAPi[4] = 0.033474;
	AAPi[5] = 0.038255;
	AAPi[6] = 0.04953;
	AAPi[7] = 0.088612;
	AAPi[8] = 0.033619;
	AAPi[9] = 0.036886;
	AAPi[10] = 0.085357;
	AAPi[11] = 0.080481;
	AAPi[12] = 0.014753;
	AAPi[13] = 0.039772;
	AAPi[14] = 0.05068;
	AAPi[15] = 0.069577;
	AAPi[16] = 0.058542;
	AAPi[17] = 0.010494;
	AAPi[18] = 0.029916;
	AAPi[19] = 0.064718;


	/*Dayhoff DCMUT Test*/
	AAmatrixT[0][0] = 0;			AAmatrixT[0][1] = 26.7828;		AAmatrixT[0][2] = 98.4474;		AAmatrixT[0][3] = 119.981;	
	AAmatrixT[0][4] = 36.0016;		AAmatrixT[0][5] = 88.7753;		AAmatrixT[0][6] = 196.117;		AAmatrixT[0][7] = 238.611;	
	AAmatrixT[0][8] = 22.8116;		AAmatrixT[0][9] = 65.3416;		AAmatrixT[0][10] = 40.6431;		AAmatrixT[0][11] = 25.8635;	
	AAmatrixT[0][12] = 71.784;		AAmatrixT[0][13] = 18.3641;		AAmatrixT[0][14] = 248.592;		AAmatrixT[0][15] = 405.187;	
	AAmatrixT[0][16] = 368.036;		AAmatrixT[0][17] = 0;			AAmatrixT[0][18] = 24.4139;		AAmatrixT[0][19] = 205.956;	
	AAmatrixT[1][0] = 26.7828;		AAmatrixT[1][1] = 0;			AAmatrixT[1][2] = 32.7059;		AAmatrixT[1][3] = 0;	
	AAmatrixT[1][4] = 23.2374;		AAmatrixT[1][5] = 243.994;		AAmatrixT[1][6] = 0;			AAmatrixT[1][7] = 8.7791;	
	AAmatrixT[1][8] = 238.315;		AAmatrixT[1][9] = 63.2629;		AAmatrixT[1][10] = 15.4924;		AAmatrixT[1][11] = 461.012;	
	AAmatrixT[1][12] = 89.6321;		AAmatrixT[1][13] = 13.6906;		AAmatrixT[1][14] = 102.831;		AAmatrixT[1][15] = 153.159;	
	AAmatrixT[1][16] = 26.5745;		AAmatrixT[1][17] = 200.137;		AAmatrixT[1][18] = 7.8012;		AAmatrixT[1][19] = 24.0368;	
	AAmatrixT[2][0] = 98.4474;		AAmatrixT[2][1] = 32.7059;		AAmatrixT[2][2] = 0;			AAmatrixT[2][3] = 893.152;	
	AAmatrixT[2][4] = 0;			AAmatrixT[2][5] = 102.851;		AAmatrixT[2][6] = 149.341;		AAmatrixT[2][7] = 138.535;	
	AAmatrixT[2][8] = 529.002;		AAmatrixT[2][9] = 76.8024;		AAmatrixT[2][10] = 34.1113;		AAmatrixT[2][11] = 314.837;	
	AAmatrixT[2][12] = 0;			AAmatrixT[2][13] = 13.8503;		AAmatrixT[2][14] = 41.9244;		AAmatrixT[2][15] = 488.589;	
	AAmatrixT[2][16] = 227.17;		AAmatrixT[2][17] = 22.4968;		AAmatrixT[2][18] = 94.694;		AAmatrixT[2][19] = 15.8067;	
	AAmatrixT[3][0] = 119.981;		AAmatrixT[3][1] = 0;			AAmatrixT[3][2] = 893.152;		AAmatrixT[3][3] = 0;	
	AAmatrixT[3][4] = 0;			AAmatrixT[3][5] = 134.855;		AAmatrixT[3][6] = 1138.87;		AAmatrixT[3][7] = 124.098;	
	AAmatrixT[3][8] = 86.8241;		AAmatrixT[3][9] = 23.9248;		AAmatrixT[3][10] = 0;			AAmatrixT[3][11] = 71.6913;	
	AAmatrixT[3][12] = 0;			AAmatrixT[3][13] = 0;			AAmatrixT[3][14] = 13.394;		AAmatrixT[3][15] = 95.6097;	
	AAmatrixT[3][16] = 66.093;		AAmatrixT[3][17] = 0;			AAmatrixT[3][18] = 0;			AAmatrixT[3][19] = 17.8316;	
	AAmatrixT[4][0] = 36.0016;		AAmatrixT[4][1] = 23.2374;		AAmatrixT[4][2] = 0;			AAmatrixT[4][3] = 0;	
	AAmatrixT[4][4] = 0;			AAmatrixT[4][5] = 0;			AAmatrixT[4][6] = 0;			AAmatrixT[4][7] = 10.7278;	
	AAmatrixT[4][8] = 28.2729;		AAmatrixT[4][9] = 43.8074;		AAmatrixT[4][10] = 0;			AAmatrixT[4][11] = 0;	
	AAmatrixT[4][12] = 0;			AAmatrixT[4][13] = 0;			AAmatrixT[4][14] = 18.755;		AAmatrixT[4][15] = 159.836;	
	AAmatrixT[4][16] = 16.2366;		AAmatrixT[4][17] = 0;			AAmatrixT[4][18] = 95.3164;		AAmatrixT[4][19] = 48.4678;	
	AAmatrixT[5][0] = 88.7753;		AAmatrixT[5][1] = 243.994;		AAmatrixT[5][2] = 102.851;		AAmatrixT[5][3] = 134.855;	
	AAmatrixT[5][4] = 0;			AAmatrixT[5][5] = 0;			AAmatrixT[5][6] = 708.602;		AAmatrixT[5][7] = 28.1581;	
	AAmatrixT[5][8] = 601.161;		AAmatrixT[5][9] = 18.0393;		AAmatrixT[5][10] = 73.0772;		AAmatrixT[5][11] = 151.908;	
	AAmatrixT[5][12] = 112.75;		AAmatrixT[5][13] = 0;			AAmatrixT[5][14] = 152.619;		AAmatrixT[5][15] = 56.1828;	
	AAmatrixT[5][16] = 52.5651;		AAmatrixT[5][17] = 0;			AAmatrixT[5][18] = 0;			AAmatrixT[5][19] = 34.6983;	
	AAmatrixT[6][0] = 196.117;		AAmatrixT[6][1] = 0;			AAmatrixT[6][2] = 149.341;		AAmatrixT[6][3] = 1138.87;	
	AAmatrixT[6][4] = 0;			AAmatrixT[6][5] = 708.602;		AAmatrixT[6][6] = 0;			AAmatrixT[6][7] = 81.1907;	
	AAmatrixT[6][8] = 43.9469;		AAmatrixT[6][9] = 60.9526;		AAmatrixT[6][10] = 11.288;		AAmatrixT[6][11] = 83.0078;	
	AAmatrixT[6][12] = 30.4803;		AAmatrixT[6][13] = 0;			AAmatrixT[6][14] = 50.7003;		AAmatrixT[6][15] = 79.3999;	
	AAmatrixT[6][16] = 34.0156;		AAmatrixT[6][17] = 0;			AAmatrixT[6][18] = 21.4717;		AAmatrixT[6][19] = 36.725;	
	AAmatrixT[7][0] = 238.611;		AAmatrixT[7][1] = 8.7791;		AAmatrixT[7][2] = 138.535;		AAmatrixT[7][3] = 124.098;	
	AAmatrixT[7][4] = 10.7278;		AAmatrixT[7][5] = 28.1581;		AAmatrixT[7][6] = 81.1907;		AAmatrixT[7][7] = 0;	
	AAmatrixT[7][8] = 10.6802;		AAmatrixT[7][9] = 0;			AAmatrixT[7][10] = 7.1514;		AAmatrixT[7][11] = 26.7683;	
	AAmatrixT[7][12] = 17.0372;		AAmatrixT[7][13] = 15.3478;		AAmatrixT[7][14] = 34.7153;		AAmatrixT[7][15] = 232.224;	
	AAmatrixT[7][16] = 30.6662;		AAmatrixT[7][17] = 0;			AAmatrixT[7][18] = 0;			AAmatrixT[7][19] = 53.8165;	
	AAmatrixT[8][0] = 22.8116;		AAmatrixT[8][1] = 238.315;		AAmatrixT[8][2] = 529.002;		AAmatrixT[8][3] = 86.8241;	
	AAmatrixT[8][4] = 28.2729;		AAmatrixT[8][5] = 601.161;		AAmatrixT[8][6] = 43.9469;		AAmatrixT[8][7] = 10.6802;	
	AAmatrixT[8][8] = 0;			AAmatrixT[8][9] = 7.6981;		AAmatrixT[8][10] = 44.3504;		AAmatrixT[8][11] = 27.0475;	
	AAmatrixT[8][12] = 0;			AAmatrixT[8][13] = 47.5927;		AAmatrixT[8][14] = 93.3709;		AAmatrixT[8][15] = 35.3643;	
	AAmatrixT[8][16] = 22.6333;		AAmatrixT[8][17] = 27.0564;		AAmatrixT[8][18] = 126.54;		AAmatrixT[8][19] = 43.8715;	
	AAmatrixT[9][0] = 65.3416;		AAmatrixT[9][1] = 63.2629;		AAmatrixT[9][2] = 76.8024;		AAmatrixT[9][3] = 23.9248;	
	AAmatrixT[9][4] = 43.8074;		AAmatrixT[9][5] = 18.0393;		AAmatrixT[9][6] = 60.9526;		AAmatrixT[9][7] = 0;	
	AAmatrixT[9][8] = 7.6981;		AAmatrixT[9][9] = 0;			AAmatrixT[9][10] = 255.668;		AAmatrixT[9][11] = 46.0857;	
	AAmatrixT[9][12] = 333.273;		AAmatrixT[9][13] = 195.195;		AAmatrixT[9][14] = 11.9152;		AAmatrixT[9][15] = 24.7955;	
	AAmatrixT[9][16] = 190.074;		AAmatrixT[9][17] = 0;			AAmatrixT[9][18] = 37.4834;		AAmatrixT[9][19] = 881.004;	
	AAmatrixT[10][0] = 40.6431;		AAmatrixT[10][1] = 15.4924;		AAmatrixT[10][2] = 34.1113;		AAmatrixT[10][3] = 0;	
	AAmatrixT[10][4] = 0;			AAmatrixT[10][5] = 73.0772;		AAmatrixT[10][6] = 11.288;		AAmatrixT[10][7] = 7.1514;	
	AAmatrixT[10][8] = 44.3504;		AAmatrixT[10][9] = 255.668;		AAmatrixT[10][10] = 0;			AAmatrixT[10][11] = 18.0629;	
	AAmatrixT[10][12] = 523.011;	AAmatrixT[10][13] = 156.516;	AAmatrixT[10][14] = 31.6258;	AAmatrixT[10][15] = 17.1432;	
	AAmatrixT[10][16] = 33.109;		AAmatrixT[10][17] = 46.1776;	AAmatrixT[10][18] = 28.6572;	AAmatrixT[10][19] = 174.516;	
	AAmatrixT[11][0] = 25.8635;		AAmatrixT[11][1] = 461.012;		AAmatrixT[11][2] = 314.837;		AAmatrixT[11][3] = 71.6913;	
	AAmatrixT[11][4] = 0;			AAmatrixT[11][5] = 151.908;		AAmatrixT[11][6] = 83.0078;		AAmatrixT[11][7] = 26.7683;	
	AAmatrixT[11][8] = 27.0475;		AAmatrixT[11][9] = 46.0857;		AAmatrixT[11][10] = 18.0629;	AAmatrixT[11][11] = 0;	
	AAmatrixT[11][12] = 241.174;	AAmatrixT[11][13] = 0;			AAmatrixT[11][14] = 33.5419;	AAmatrixT[11][15] = 95.4557;	
	AAmatrixT[11][16] = 135.06;		AAmatrixT[11][17] = 0;			AAmatrixT[11][18] = 13.2142;	AAmatrixT[11][19] = 10.385;	
	AAmatrixT[12][0] = 71.784;		AAmatrixT[12][1] = 89.6321;		AAmatrixT[12][2] = 0;			AAmatrixT[12][3] = 0;	
	AAmatrixT[12][4] = 0;			AAmatrixT[12][5] = 112.75;		AAmatrixT[12][6] = 30.4803;		AAmatrixT[12][7] = 17.0372;	
	AAmatrixT[12][8] = 0;			AAmatrixT[12][9] = 333.273;		AAmatrixT[12][10] = 523.011;	AAmatrixT[12][11] = 241.174;	
	AAmatrixT[12][12] = 0;			AAmatrixT[12][13] = 92.186;		AAmatrixT[12][14] = 17.0205;	AAmatrixT[12][15] = 61.9951;	
	AAmatrixT[12][16] = 103.153;	AAmatrixT[12][17] = 0;			AAmatrixT[12][18] = 0;			AAmatrixT[12][19] = 256.596;	
	AAmatrixT[13][0] = 18.3641;		AAmatrixT[13][1] = 13.6906;		AAmatrixT[13][2] = 13.8503;		AAmatrixT[13][3] = 0;	
	AAmatrixT[13][4] = 0;			AAmatrixT[13][5] = 0;			AAmatrixT[13][6] = 0;			AAmatrixT[13][7] = 15.3478;	
	AAmatrixT[13][8] = 47.5927;		AAmatrixT[13][9] = 195.195;		AAmatrixT[13][10] = 156.516;	AAmatrixT[13][11] = 0;	
	AAmatrixT[13][12] = 92.186;		AAmatrixT[13][13] = 0;			AAmatrixT[13][14] = 11.0506;	AAmatrixT[13][15] = 45.9901;	
	AAmatrixT[13][16] = 13.6655;	AAmatrixT[13][17] = 76.2354;	AAmatrixT[13][18] = 695.263;	AAmatrixT[13][19] = 12.3606;	
	AAmatrixT[14][0] = 248.592;		AAmatrixT[14][1] = 102.831;		AAmatrixT[14][2] = 41.9244;		AAmatrixT[14][3] = 13.394;	
	AAmatrixT[14][4] = 18.755;		AAmatrixT[14][5] = 152.619;		AAmatrixT[14][6] = 50.7003;		AAmatrixT[14][7] = 34.7153;	
	AAmatrixT[14][8] = 93.3709;		AAmatrixT[14][9] = 11.9152;		AAmatrixT[14][10] = 31.6258;	AAmatrixT[14][11] = 33.5419;	
	AAmatrixT[14][12] = 17.0205;	AAmatrixT[14][13] = 11.0506;	AAmatrixT[14][14] = 0;			AAmatrixT[14][15] = 242.72;	
	AAmatrixT[14][16] = 78.2857;	AAmatrixT[14][17] = 0;			AAmatrixT[14][18] = 0;			AAmatrixT[14][19] = 48.5026;	
	AAmatrixT[15][0] = 405.187;		AAmatrixT[15][1] = 153.159;		AAmatrixT[15][2] = 488.589;		AAmatrixT[15][3] = 95.6097;	
	AAmatrixT[15][4] = 159.836;		AAmatrixT[15][5] = 56.1828;		AAmatrixT[15][6] = 79.3999;		AAmatrixT[15][7] = 232.224;	
	AAmatrixT[15][8] = 35.3643;		AAmatrixT[15][9] = 24.7955;		AAmatrixT[15][10] = 17.1432;	AAmatrixT[15][11] = 95.4557;	
	AAmatrixT[15][12] = 61.9951;	AAmatrixT[15][13] = 45.9901;	AAmatrixT[15][14] = 242.72;		AAmatrixT[15][15] = 0;	
	AAmatrixT[15][16] = 543.667;	AAmatrixT[15][17] = 74.0819;	AAmatrixT[15][18] = 33.6289;	AAmatrixT[15][19] = 30.3836;	
	AAmatrixT[16][0] = 368.036;		AAmatrixT[16][1] = 26.5745;		AAmatrixT[16][2] = 227.17;		AAmatrixT[16][3] = 66.093;	
	AAmatrixT[16][4] = 16.2366;		AAmatrixT[16][5] = 52.5651;		AAmatrixT[16][6] = 34.0156;		AAmatrixT[16][7] = 30.6662;	
	AAmatrixT[16][8] = 22.6333;		AAmatrixT[16][9] = 190.074;		AAmatrixT[16][10] = 33.109;		AAmatrixT[16][11] = 135.06;	
	AAmatrixT[16][12] = 103.153;	AAmatrixT[16][13] = 13.6655;	AAmatrixT[16][14] = 78.2857;	AAmatrixT[16][15] = 543.667;	
	AAmatrixT[16][16] = 0;			AAmatrixT[16][17] = 0;			AAmatrixT[16][18] = 41.7839;	AAmatrixT[16][19] = 156.2;	
	AAmatrixT[17][0] = 0;			AAmatrixT[17][1] = 200.137;		AAmatrixT[17][2] = 22.4968;		AAmatrixT[17][3] = 0;	
	AAmatrixT[17][4] = 0;			AAmatrixT[17][5] = 0;			AAmatrixT[17][6] = 0;			AAmatrixT[17][7] = 0;	
	AAmatrixT[17][8] = 27.0564;		AAmatrixT[17][9] = 0;			AAmatrixT[17][10] = 46.1776;	AAmatrixT[17][11] = 0;	
	AAmatrixT[17][12] = 0;			AAmatrixT[17][13] = 76.2354;	AAmatrixT[17][14] = 0;			AAmatrixT[17][15] = 74.0819;	
	AAmatrixT[17][16] = 0;			AAmatrixT[17][17] = 0;			AAmatrixT[17][18] = 60.807;		AAmatrixT[17][19] = 0;	
	AAmatrixT[18][0] = 24.4139;		AAmatrixT[18][1] = 7.8012;		AAmatrixT[18][2] = 94.694;		AAmatrixT[18][3] = 0;	
	AAmatrixT[18][4] = 95.3164;		AAmatrixT[18][5] = 0;			AAmatrixT[18][6] = 21.4717;		AAmatrixT[18][7] = 0;	
	AAmatrixT[18][8] = 126.54;		AAmatrixT[18][9] = 37.4834;		AAmatrixT[18][10] = 28.6572;	AAmatrixT[18][11] = 13.2142;	
	AAmatrixT[18][12] = 0;			AAmatrixT[18][13] = 695.263;	AAmatrixT[18][14] = 0;			AAmatrixT[18][15] = 33.6289;	
	AAmatrixT[18][16] = 41.7839;	AAmatrixT[18][17] = 60.807;		AAmatrixT[18][18] = 0;			AAmatrixT[18][19] = 27.9379;	
	AAmatrixT[19][0] = 205.956;		AAmatrixT[19][1] = 24.0368;		AAmatrixT[19][2] = 15.8067;		AAmatrixT[19][3] = 17.8316;	
	AAmatrixT[19][4] = 48.4678;		AAmatrixT[19][5] = 34.6983;		AAmatrixT[19][6] = 36.725;		AAmatrixT[19][7] = 53.8165;	
	AAmatrixT[19][8] = 43.8715;		AAmatrixT[19][9] = 881.004;		AAmatrixT[19][10] = 174.516;	AAmatrixT[19][11] = 10.385;	
	AAmatrixT[19][12] = 256.596;	AAmatrixT[19][13] = 12.3606;	AAmatrixT[19][14] = 48.5026;	AAmatrixT[19][15] = 30.3836;	
	AAmatrixT[19][16] = 156.2;		AAmatrixT[19][17] = 0;			AAmatrixT[19][18] = 27.9379;	AAmatrixT[19][19] = 0;	

	AAPiT[0] = 0.087127;
	AAPiT[1] = 0.040904;
	AAPiT[2] = 0.040432;
	AAPiT[3] = 0.046872;
	AAPiT[4] = 0.033474;
	AAPiT[5] = 0.038255;
	AAPiT[6] = 0.04953;
	AAPiT[7] = 0.088612;
	AAPiT[8] = 0.033619;
	AAPiT[9] = 0.036886;
	AAPiT[10] = 0.085357;
	AAPiT[11] = 0.080481;
	AAPiT[12] = 0.014753;
	AAPiT[13] = 0.039772;
	AAPiT[14] = 0.05068;
	AAPiT[15] = 0.069577;
	AAPiT[16] = 0.058542;
	AAPiT[17] = 0.010494;
	AAPiT[18] = 0.029916;
	AAPiT[19] = 0.064718;
	
		for(int gb1=0; gb1<20; gb1++)
		{
			if(AAPiT[gb1]!=AAPiT[gb1]) {cout<<"ERROR in base frequency translation check in Jones-DCMUT model"<<endl; gb1--;}
			
			for(int gb2=0; gb2<20; gb2++)
			{
				if(AAmatrix[gb1][gb2]-AAmatrixT[gb1][gb2]/100>0.0001) {cout<<"ERROR in substitution matrix translation check in Jones-DCMUT model"<<endl; gb2--;}
			}
		}

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else
	if(modelnumber==5)
	{

	/* 
		wag 
		
		Whelan, S. and Goldman, N. (2001)  A general empirical model of protein evolution derived from multiple protein families using a 
		maximum-likelihood approach. Molecular Biology and Evolution 18:691-699
	*/
	AAmatrix[ 0][ 0] = 0.0000000; AAmatrix[ 1][ 0] = 0.5515710; AAmatrix[ 2][ 0] = 0.5098480; AAmatrix[ 3][ 0] = 0.7389980; AAmatrix[ 4][ 0] = 1.0270400; 
	AAmatrix[ 5][ 0] = 0.9085980; AAmatrix[ 6][ 0] = 1.5828500; AAmatrix[ 7][ 0] = 1.4167200; AAmatrix[ 8][ 0] = 0.3169540; AAmatrix[ 9][ 0] = 0.1933350; 
	AAmatrix[10][ 0] = 0.3979150; AAmatrix[11][ 0] = 0.9062650; AAmatrix[12][ 0] = 0.8934960; AAmatrix[13][ 0] = 0.2104940; AAmatrix[14][ 0] = 1.4385500; 
	AAmatrix[15][ 0] = 3.3707900; AAmatrix[16][ 0] = 2.1211100; AAmatrix[17][ 0] = 0.1131330; AAmatrix[18][ 0] = 0.2407350; AAmatrix[19][ 0] = 2.0060100;
	AAmatrix[ 0][ 1] = 0.5515710; AAmatrix[ 1][ 1] = 0.0000000; AAmatrix[ 2][ 1] = 0.6353460; AAmatrix[ 3][ 1] = 0.1473040; AAmatrix[ 4][ 1] = 0.5281910;  
	AAmatrix[ 5][ 1] = 3.0355000; AAmatrix[ 6][ 1] = 0.4391570; AAmatrix[ 7][ 1] = 0.5846650; AAmatrix[ 8][ 1] = 2.1371500; AAmatrix[ 9][ 1] = 0.1869790;  
	AAmatrix[10][ 1] = 0.4976710; AAmatrix[11][ 1] = 5.3514200; AAmatrix[12][ 1] = 0.6831620; AAmatrix[13][ 1] = 0.1027110; AAmatrix[14][ 1] = 0.6794890;  
	AAmatrix[15][ 1] = 1.2241900; AAmatrix[16][ 1] = 0.5544130; AAmatrix[17][ 1] = 1.1639200; AAmatrix[18][ 1] = 0.3815330; AAmatrix[19][ 1] = 0.2518490;
	AAmatrix[ 0][ 2] = 0.5098480; AAmatrix[ 1][ 2] = 0.6353460; AAmatrix[ 2][ 2] = 0.0000000; AAmatrix[ 3][ 2] = 5.4294200; AAmatrix[ 4][ 2] = 0.2652560;  
	AAmatrix[ 5][ 2] = 1.5436400; AAmatrix[ 6][ 2] = 0.9471980; AAmatrix[ 7][ 2] = 1.1255600; AAmatrix[ 8][ 2] = 3.9562900; AAmatrix[ 9][ 2] = 0.5542360;  
	AAmatrix[10][ 2] = 0.1315280; AAmatrix[11][ 2] = 3.0120100; AAmatrix[12][ 2] = 0.1982210; AAmatrix[13][ 2] = 0.0961621; AAmatrix[14][ 2] = 0.1950810;  
	AAmatrix[15][ 2] = 3.9742300; AAmatrix[16][ 2] = 2.0300600; AAmatrix[17][ 2] = 0.0719167; AAmatrix[18][ 2] = 1.0860000; AAmatrix[19][ 2] = 0.1962460;
	AAmatrix[ 0][ 3] = 0.7389980; AAmatrix[ 1][ 3] = 0.1473040; AAmatrix[ 2][ 3] = 5.4294200; AAmatrix[ 3][ 3] = 0.0000000; AAmatrix[ 4][ 3] = 0.0302949;  
	AAmatrix[ 5][ 3] = 0.6167830; AAmatrix[ 6][ 3] = 6.1741600; AAmatrix[ 7][ 3] = 0.8655840; AAmatrix[ 8][ 3] = 0.9306760; AAmatrix[ 9][ 3] = 0.0394370;  
	AAmatrix[10][ 3] = 0.0848047; AAmatrix[11][ 3] = 0.4798550; AAmatrix[12][ 3] = 0.1037540; AAmatrix[13][ 3] = 0.0467304; AAmatrix[14][ 3] = 0.4239840;  
	AAmatrix[15][ 3] = 1.0717600; AAmatrix[16][ 3] = 0.3748660; AAmatrix[17][ 3] = 0.1297670; AAmatrix[18][ 3] = 0.3257110; AAmatrix[19][ 3] = 0.1523350;
	AAmatrix[ 0][ 4] = 1.0270400; AAmatrix[ 1][ 4] = 0.5281910; AAmatrix[ 2][ 4] = 0.2652560; AAmatrix[ 3][ 4] = 0.0302949; AAmatrix[ 4][ 4] = 0.0000000;  
	AAmatrix[ 5][ 4] = 0.0988179; AAmatrix[ 6][ 4] = 0.0213520; AAmatrix[ 7][ 4] = 0.3066740; AAmatrix[ 8][ 4] = 0.2489720; AAmatrix[ 9][ 4] = 0.1701350;  
	AAmatrix[10][ 4] = 0.3842870; AAmatrix[11][ 4] = 0.0740339; AAmatrix[12][ 4] = 0.3904820; AAmatrix[13][ 4] = 0.3980200; AAmatrix[14][ 4] = 0.1094040;  
	AAmatrix[15][ 4] = 1.4076600; AAmatrix[16][ 4] = 0.5129840; AAmatrix[17][ 4] = 0.7170700; AAmatrix[18][ 4] = 0.5438330; AAmatrix[19][ 4] = 1.0021400;
	AAmatrix[ 0][ 5] = 0.9085980; AAmatrix[ 1][ 5] = 3.0355000; AAmatrix[ 2][ 5] = 1.5436400; AAmatrix[ 3][ 5] = 0.6167830; AAmatrix[ 4][ 5] = 0.0988179;  
	AAmatrix[ 5][ 5] = 0.0000000; AAmatrix[ 6][ 5] = 5.4694700; AAmatrix[ 7][ 5] = 0.3300520; AAmatrix[ 8][ 5] = 4.2941100; AAmatrix[ 9][ 5] = 0.1139170;  
	AAmatrix[10][ 5] = 0.8694890; AAmatrix[11][ 5] = 3.8949000; AAmatrix[12][ 5] = 1.5452600; AAmatrix[13][ 5] = 0.0999208; AAmatrix[14][ 5] = 0.9333720;  
	AAmatrix[15][ 5] = 1.0288700; AAmatrix[16][ 5] = 0.8579280; AAmatrix[17][ 5] = 0.2157370; AAmatrix[18][ 5] = 0.2277100; AAmatrix[19][ 5] = 0.3012810;
	AAmatrix[ 0][ 6] = 1.5828500; AAmatrix[ 1][ 6] = 0.4391570; AAmatrix[ 2][ 6] = 0.9471980; AAmatrix[ 3][ 6] = 6.1741600; AAmatrix[ 4][ 6] = 0.0213520;  
	AAmatrix[ 5][ 6] = 5.4694700; AAmatrix[ 6][ 6] = 0.0000000; AAmatrix[ 7][ 6] = 0.5677170; AAmatrix[ 8][ 6] = 0.5700250; AAmatrix[ 9][ 6] = 0.1273950;  
	AAmatrix[10][ 6] = 0.1542630; AAmatrix[11][ 6] = 2.5844300; AAmatrix[12][ 6] = 0.3151240; AAmatrix[13][ 6] = 0.0811339; AAmatrix[14][ 6] = 0.6823550;  
	AAmatrix[15][ 6] = 0.7049390; AAmatrix[16][ 6] = 0.8227650; AAmatrix[17][ 6] = 0.1565570; AAmatrix[18][ 6] = 0.1963030; AAmatrix[19][ 6] = 0.5887310;
	AAmatrix[ 0][ 7] = 1.4167200; AAmatrix[ 1][ 7] = 0.5846650; AAmatrix[ 2][ 7] = 1.1255600; AAmatrix[ 3][ 7] = 0.8655840; AAmatrix[ 4][ 7] = 0.3066740;  
	AAmatrix[ 5][ 7] = 0.3300520; AAmatrix[ 6][ 7] = 0.5677170; AAmatrix[ 7][ 7] = 0.0000000; AAmatrix[ 8][ 7] = 0.2494100; AAmatrix[ 9][ 7] = 0.0304501;  
	AAmatrix[10][ 7] = 0.0613037; AAmatrix[11][ 7] = 0.3735580; AAmatrix[12][ 7] = 0.1741000; AAmatrix[13][ 7] = 0.0499310; AAmatrix[14][ 7] = 0.2435700;  
	AAmatrix[15][ 7] = 1.3418200; AAmatrix[16][ 7] = 0.2258330; AAmatrix[17][ 7] = 0.3369830; AAmatrix[18][ 7] = 0.1036040; AAmatrix[19][ 7] = 0.1872470;
	AAmatrix[ 0][ 8] = 0.3169540; AAmatrix[ 1][ 8] = 2.1371500; AAmatrix[ 2][ 8] = 3.9562900; AAmatrix[ 3][ 8] = 0.9306760; AAmatrix[ 4][ 8] = 0.2489720;  
	AAmatrix[ 5][ 8] = 4.2941100; AAmatrix[ 6][ 8] = 0.5700250; AAmatrix[ 7][ 8] = 0.2494100; AAmatrix[ 8][ 8] = 0.0000000; AAmatrix[ 9][ 8] = 0.1381900;  
	AAmatrix[10][ 8] = 0.4994620; AAmatrix[11][ 8] = 0.8904320; AAmatrix[12][ 8] = 0.4041410; AAmatrix[13][ 8] = 0.6793710; AAmatrix[14][ 8] = 0.6961980;  
	AAmatrix[15][ 8] = 0.7401690; AAmatrix[16][ 8] = 0.4733070; AAmatrix[17][ 8] = 0.2625690; AAmatrix[18][ 8] = 3.8734400; AAmatrix[19][ 8] = 0.1183580;
	AAmatrix[ 0][ 9] = 0.1933350; AAmatrix[ 1][ 9] = 0.1869790; AAmatrix[ 2][ 9] = 0.5542360; AAmatrix[ 3][ 9] = 0.0394370; AAmatrix[ 4][ 9] = 0.1701350;  
	AAmatrix[ 5][ 9] = 0.1139170; AAmatrix[ 6][ 9] = 0.1273950; AAmatrix[ 7][ 9] = 0.0304501; AAmatrix[ 8][ 9] = 0.1381900; AAmatrix[ 9][ 9] = 0.0000000;  
	AAmatrix[10][ 9] = 3.1709700; AAmatrix[11][ 9] = 0.3238320; AAmatrix[12][ 9] = 4.2574600; AAmatrix[13][ 9] = 1.0594700; AAmatrix[14][ 9] = 0.0999288;  
	AAmatrix[15][ 9] = 0.3194400; AAmatrix[16][ 9] = 1.4581600; AAmatrix[17][ 9] = 0.2124830; AAmatrix[18][ 9] = 0.4201700; AAmatrix[19][ 9] = 7.8213000;
	AAmatrix[ 0][10] = 0.3979150; AAmatrix[ 1][10] = 0.4976710; AAmatrix[ 2][10] = 0.1315280; AAmatrix[ 3][10] = 0.0848047; AAmatrix[ 4][10] = 0.3842870;  
	AAmatrix[ 5][10] = 0.8694890; AAmatrix[ 6][10] = 0.1542630; AAmatrix[ 7][10] = 0.0613037; AAmatrix[ 8][10] = 0.4994620; AAmatrix[ 9][10] = 3.1709700;  
	AAmatrix[10][10] = 0.0000000; AAmatrix[11][10] = 0.2575550; AAmatrix[12][10] = 4.8540200; AAmatrix[13][10] = 2.1151700; AAmatrix[14][10] = 0.4158440;  
	AAmatrix[15][10] = 0.3447390; AAmatrix[16][10] = 0.3266220; AAmatrix[17][10] = 0.6653090; AAmatrix[18][10] = 0.3986180; AAmatrix[19][10] = 1.8003400;
	AAmatrix[ 0][11] = 0.9062650; AAmatrix[ 1][11] = 5.3514200; AAmatrix[ 2][11] = 3.0120100; AAmatrix[ 3][11] = 0.4798550; AAmatrix[ 4][11] = 0.0740339;  
	AAmatrix[ 5][11] = 3.8949000; AAmatrix[ 6][11] = 2.5844300; AAmatrix[ 7][11] = 0.3735580; AAmatrix[ 8][11] = 0.8904320; AAmatrix[ 9][11] = 0.3238320;  
	AAmatrix[10][11] = 0.2575550; AAmatrix[11][11] = 0.0000000; AAmatrix[12][11] = 0.9342760; AAmatrix[13][11] = 0.0888360; AAmatrix[14][11] = 0.5568960;  
	AAmatrix[15][11] = 0.9671300; AAmatrix[16][11] = 1.3869800; AAmatrix[17][11] = 0.1375050; AAmatrix[18][11] = 0.1332640; AAmatrix[19][11] = 0.3054340;
	AAmatrix[ 0][12] = 0.8934960; AAmatrix[ 1][12] = 0.6831620; AAmatrix[ 2][12] = 0.1982210; AAmatrix[ 3][12] = 0.1037540; AAmatrix[ 4][12] = 0.3904820;  
	AAmatrix[ 5][12] = 1.5452600; AAmatrix[ 6][12] = 0.3151240; AAmatrix[ 7][12] = 0.1741000; AAmatrix[ 8][12] = 0.4041410; AAmatrix[ 9][12] = 4.2574600;  
	AAmatrix[10][12] = 4.8540200; AAmatrix[11][12] = 0.9342760; AAmatrix[12][12] = 0.0000000; AAmatrix[13][12] = 1.1906300; AAmatrix[14][12] = 0.1713290;  
	AAmatrix[15][12] = 0.4939050; AAmatrix[16][12] = 1.5161200; AAmatrix[17][12] = 0.5157060; AAmatrix[18][12] = 0.4284370; AAmatrix[19][12] = 2.0584500;
	AAmatrix[ 0][13] = 0.2104940; AAmatrix[ 1][13] = 0.1027110; AAmatrix[ 2][13] = 0.0961621; AAmatrix[ 3][13] = 0.0467304; AAmatrix[ 4][13] = 0.3980200;  
	AAmatrix[ 5][13] = 0.0999208; AAmatrix[ 6][13] = 0.0811339; AAmatrix[ 7][13] = 0.0499310; AAmatrix[ 8][13] = 0.6793710; AAmatrix[ 9][13] = 1.0594700;  
	AAmatrix[10][13] = 2.1151700; AAmatrix[11][13] = 0.0888360; AAmatrix[12][13] = 1.1906300; AAmatrix[13][13] = 0.0000000; AAmatrix[14][13] = 0.1614440;  
	AAmatrix[15][13] = 0.5459310; AAmatrix[16][13] = 0.1719030; AAmatrix[17][13] = 1.5296400; AAmatrix[18][13] = 6.4542800; AAmatrix[19][13] = 0.6498920;
	AAmatrix[ 0][14] = 1.4385500; AAmatrix[ 1][14] = 0.6794890; AAmatrix[ 2][14] = 0.1950810; AAmatrix[ 3][14] = 0.4239840; AAmatrix[ 4][14] = 0.1094040;  
	AAmatrix[ 5][14] = 0.9333720; AAmatrix[ 6][14] = 0.6823550; AAmatrix[ 7][14] = 0.2435700; AAmatrix[ 8][14] = 0.6961980; AAmatrix[ 9][14] = 0.0999288;  
	AAmatrix[10][14] = 0.4158440; AAmatrix[11][14] = 0.5568960; AAmatrix[12][14] = 0.1713290; AAmatrix[13][14] = 0.1614440; AAmatrix[14][14] = 0.0000000;  
	AAmatrix[15][14] = 1.6132800; AAmatrix[16][14] = 0.7953840; AAmatrix[17][14] = 0.1394050; AAmatrix[18][14] = 0.2160460; AAmatrix[19][14] = 0.3148870;
	AAmatrix[ 0][15] = 3.3707900; AAmatrix[ 1][15] = 1.2241900; AAmatrix[ 2][15] = 3.9742300; AAmatrix[ 3][15] = 1.0717600; AAmatrix[ 4][15] = 1.4076600;  
	AAmatrix[ 5][15] = 1.0288700; AAmatrix[ 6][15] = 0.7049390; AAmatrix[ 7][15] = 1.3418200; AAmatrix[ 8][15] = 0.7401690; AAmatrix[ 9][15] = 0.3194400;  
	AAmatrix[10][15] = 0.3447390; AAmatrix[11][15] = 0.9671300; AAmatrix[12][15] = 0.4939050; AAmatrix[13][15] = 0.5459310; AAmatrix[14][15] = 1.6132800;  
	AAmatrix[15][15] = 0.0000000; AAmatrix[16][15] = 4.3780200; AAmatrix[17][15] = 0.5237420; AAmatrix[18][15] = 0.7869930; AAmatrix[19][15] = 0.2327390;
	AAmatrix[ 0][16] = 2.1211100; AAmatrix[ 1][16] = 0.5544130; AAmatrix[ 2][16] = 2.0300600; AAmatrix[ 3][16] = 0.3748660; AAmatrix[ 4][16] = 0.5129840;  
	AAmatrix[ 5][16] = 0.8579280; AAmatrix[ 6][16] = 0.8227650; AAmatrix[ 7][16] = 0.2258330; AAmatrix[ 8][16] = 0.4733070; AAmatrix[ 9][16] = 1.4581600;  
	AAmatrix[10][16] = 0.3266220; AAmatrix[11][16] = 1.3869800; AAmatrix[12][16] = 1.5161200; AAmatrix[13][16] = 0.1719030; AAmatrix[14][16] = 0.7953840;  
	AAmatrix[15][16] = 4.3780200; AAmatrix[16][16] = 0.0000000; AAmatrix[17][16] = 0.1108640; AAmatrix[18][16] = 0.2911480; AAmatrix[19][16] = 1.3882300;
	AAmatrix[ 0][17] = 0.1131330; AAmatrix[ 1][17] = 1.1639200; AAmatrix[ 2][17] = 0.0719167; AAmatrix[ 3][17] = 0.1297670; AAmatrix[ 4][17] = 0.7170700;  
	AAmatrix[ 5][17] = 0.2157370; AAmatrix[ 6][17] = 0.1565570; AAmatrix[ 7][17] = 0.3369830; AAmatrix[ 8][17] = 0.2625690; AAmatrix[ 9][17] = 0.2124830;  
	AAmatrix[10][17] = 0.6653090; AAmatrix[11][17] = 0.1375050; AAmatrix[12][17] = 0.5157060; AAmatrix[13][17] = 1.5296400; AAmatrix[14][17] = 0.1394050;  
	AAmatrix[15][17] = 0.5237420; AAmatrix[16][17] = 0.1108640; AAmatrix[17][17] = 0.0000000; AAmatrix[18][17] = 2.4853900; AAmatrix[19][17] = 0.3653690;
	AAmatrix[ 0][18] = 0.2407350; AAmatrix[ 1][18] = 0.3815330; AAmatrix[ 2][18] = 1.0860000; AAmatrix[ 3][18] = 0.3257110; AAmatrix[ 4][18] = 0.5438330;  
	AAmatrix[ 5][18] = 0.2277100; AAmatrix[ 6][18] = 0.1963030; AAmatrix[ 7][18] = 0.1036040; AAmatrix[ 8][18] = 3.8734400; AAmatrix[ 9][18] = 0.4201700;  
	AAmatrix[10][18] = 0.3986180; AAmatrix[11][18] = 0.1332640; AAmatrix[12][18] = 0.4284370; AAmatrix[13][18] = 6.4542800; AAmatrix[14][18] = 0.2160460;  
	AAmatrix[15][18] = 0.7869930; AAmatrix[16][18] = 0.2911480; AAmatrix[17][18] = 2.4853900; AAmatrix[18][18] = 0.0000000; AAmatrix[19][18] = 0.3147300;
	AAmatrix[ 0][19] = 2.0060100; AAmatrix[ 1][19] = 0.2518490; AAmatrix[ 2][19] = 0.1962460; AAmatrix[ 3][19] = 0.1523350; AAmatrix[ 4][19] = 1.0021400;  
	AAmatrix[ 5][19] = 0.3012810; AAmatrix[ 6][19] = 0.5887310; AAmatrix[ 7][19] = 0.1872470; AAmatrix[ 8][19] = 0.1183580; AAmatrix[ 9][19] = 7.8213000;  
	AAmatrix[10][19] = 1.8003400; AAmatrix[11][19] = 0.3054340; AAmatrix[12][19] = 2.0584500; AAmatrix[13][19] = 0.6498920; AAmatrix[14][19] = 0.3148870;  
	AAmatrix[15][19] = 0.2327390; AAmatrix[16][19] = 1.3882300; AAmatrix[17][19] = 0.3653690; AAmatrix[18][19] = 0.3147300; AAmatrix[19][19] = 0.0000000;
	AAPi[ 0] = 0.08662790;
	AAPi[ 1] = 0.04397200;
	AAPi[ 2] = 0.03908940;
	AAPi[ 3] = 0.05704510;
	AAPi[ 4] = 0.01930780;
	AAPi[ 5] = 0.03672810;
	AAPi[ 6] = 0.05805890;
	AAPi[ 7] = 0.08325180;
	AAPi[ 8] = 0.02443130;
	AAPi[ 9] = 0.04846600;
	AAPi[10] = 0.08620970;
	AAPi[11] = 0.06202860;
	AAPi[12] = 0.01950273;
	AAPi[13] = 0.03843190;
	AAPi[14] = 0.04576310;
	AAPi[15] = 0.06951790;
	AAPi[16] = 0.06101270;
	AAPi[17] = 0.01438590;
	AAPi[18] = 0.03527420;
	AAPi[19] = 0.07089560;

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==6)
	{

	/* 
		mtmam 
		
		Yang, Z., Nielsen, R. and Hasegawa, M. (1998) Models of amino acid substitution and applications to Mitochondrial 
		protein evolution, Molecular Biology and Evolution 15:1600-1611.
	*/
	AAmatrix[ 0][ 0] =   0; AAmatrix[ 0][ 1] =  32; AAmatrix[ 0][ 2] =   2; AAmatrix[ 0][ 3] =  11; AAmatrix[ 0][ 4] =   0;
	AAmatrix[ 0][ 5] =   0; AAmatrix[ 0][ 6] =   0; AAmatrix[ 0][ 7] =  78; AAmatrix[ 0][ 8] =   8; AAmatrix[ 0][ 9] =  75;
	AAmatrix[ 0][10] =  21; AAmatrix[ 0][11] =   0; AAmatrix[ 0][12] =  76; AAmatrix[ 0][13] =   0; AAmatrix[ 0][14] =  53;
	AAmatrix[ 0][15] = 342; AAmatrix[ 0][16] = 681; AAmatrix[ 0][17] =   5; AAmatrix[ 0][18] =   0; AAmatrix[ 0][19] = 398;
	AAmatrix[ 1][ 0] =  32; AAmatrix[ 1][ 1] =   0; AAmatrix[ 1][ 2] =   4; AAmatrix[ 1][ 3] =   0; AAmatrix[ 1][ 4] = 186;
	AAmatrix[ 1][ 5] = 246; AAmatrix[ 1][ 6] =   0; AAmatrix[ 1][ 7] =  18; AAmatrix[ 1][ 8] = 232; AAmatrix[ 1][ 9] =   0;
	AAmatrix[ 1][10] =   6; AAmatrix[ 1][11] =  50; AAmatrix[ 1][12] =   0; AAmatrix[ 1][13] =   0; AAmatrix[ 1][14] =   9;
	AAmatrix[ 1][15] =   3; AAmatrix[ 1][16] =   0; AAmatrix[ 1][17] =  16; AAmatrix[ 1][18] =   0; AAmatrix[ 1][19] =   0;
	AAmatrix[ 2][ 0] =   2; AAmatrix[ 2][ 1] =   4; AAmatrix[ 2][ 2] =   0; AAmatrix[ 2][ 3] = 864; AAmatrix[ 2][ 4] =   0;
	AAmatrix[ 2][ 5] =   8; AAmatrix[ 2][ 6] =   0; AAmatrix[ 2][ 7] =  47; AAmatrix[ 2][ 8] = 458; AAmatrix[ 2][ 9] =  19;
	AAmatrix[ 2][10] =   0; AAmatrix[ 2][11] = 408; AAmatrix[ 2][12] =  21; AAmatrix[ 2][13] =   6; AAmatrix[ 2][14] =  33;
	AAmatrix[ 2][15] = 446; AAmatrix[ 2][16] = 110; AAmatrix[ 2][17] =   6; AAmatrix[ 2][18] = 156; AAmatrix[ 2][19] =   0;
	AAmatrix[ 3][ 0] =  11; AAmatrix[ 3][ 1] =   0; AAmatrix[ 3][ 2] = 864; AAmatrix[ 3][ 3] =   0; AAmatrix[ 3][ 4] =   0;
	AAmatrix[ 3][ 5] =  49; AAmatrix[ 3][ 6] = 569; AAmatrix[ 3][ 7] =  79; AAmatrix[ 3][ 8] =  11; AAmatrix[ 3][ 9] =   0;
	AAmatrix[ 3][10] =   0; AAmatrix[ 3][11] =   0; AAmatrix[ 3][12] =   0; AAmatrix[ 3][13] =   5; AAmatrix[ 3][14] =   2;
	AAmatrix[ 3][15] =  16; AAmatrix[ 3][16] =   0; AAmatrix[ 3][17] =   0; AAmatrix[ 3][18] =   0; AAmatrix[ 3][19] =  10;
	AAmatrix[ 4][ 0] =   0; AAmatrix[ 4][ 1] = 186; AAmatrix[ 4][ 2] =   0; AAmatrix[ 4][ 3] =   0; AAmatrix[ 4][ 4] =   0;
	AAmatrix[ 4][ 5] =   0; AAmatrix[ 4][ 6] =   0; AAmatrix[ 4][ 7] =   0; AAmatrix[ 4][ 8] = 305; AAmatrix[ 4][ 9] =  41;
	AAmatrix[ 4][10] =  27; AAmatrix[ 4][11] =   0; AAmatrix[ 4][12] =   0; AAmatrix[ 4][13] =   7; AAmatrix[ 4][14] =   0;
	AAmatrix[ 4][15] = 347; AAmatrix[ 4][16] = 114; AAmatrix[ 4][17] =  65; AAmatrix[ 4][18] = 530; AAmatrix[ 4][19] =   0;
	AAmatrix[ 5][ 0] =   0; AAmatrix[ 5][ 1] = 246; AAmatrix[ 5][ 2] =   8; AAmatrix[ 5][ 3] =  49; AAmatrix[ 5][ 4] =   0;
	AAmatrix[ 5][ 5] =   0; AAmatrix[ 5][ 6] = 274; AAmatrix[ 5][ 7] =   0; AAmatrix[ 5][ 8] = 550; AAmatrix[ 5][ 9] =   0;
	AAmatrix[ 5][10] =  20; AAmatrix[ 5][11] = 242; AAmatrix[ 5][12] =  22; AAmatrix[ 5][13] =   0; AAmatrix[ 5][14] =  51;
	AAmatrix[ 5][15] =  30; AAmatrix[ 5][16] =   0; AAmatrix[ 5][17] =   0; AAmatrix[ 5][18] =  54; AAmatrix[ 5][19] =  33;
	AAmatrix[ 6][ 0] =   0; AAmatrix[ 6][ 1] =   0; AAmatrix[ 6][ 2] =   0; AAmatrix[ 6][ 3] = 569; AAmatrix[ 6][ 4] =   0;
	AAmatrix[ 6][ 5] = 274; AAmatrix[ 6][ 6] =   0; AAmatrix[ 6][ 7] =  22; AAmatrix[ 6][ 8] =  22; AAmatrix[ 6][ 9] =   0;
	AAmatrix[ 6][10] =   0; AAmatrix[ 6][11] = 215; AAmatrix[ 6][12] =   0; AAmatrix[ 6][13] =   0; AAmatrix[ 6][14] =   0;
	AAmatrix[ 6][15] =  21; AAmatrix[ 6][16] =   4; AAmatrix[ 6][17] =   0; AAmatrix[ 6][18] =   0; AAmatrix[ 6][19] =  20;
	AAmatrix[ 7][ 0] =  78; AAmatrix[ 7][ 1] =  18; AAmatrix[ 7][ 2] =  47; AAmatrix[ 7][ 3] =  79; AAmatrix[ 7][ 4] =   0;
	AAmatrix[ 7][ 5] =   0; AAmatrix[ 7][ 6] =  22; AAmatrix[ 7][ 7] =   0; AAmatrix[ 7][ 8] =   0; AAmatrix[ 7][ 9] =   0;
	AAmatrix[ 7][10] =   0; AAmatrix[ 7][11] =   0; AAmatrix[ 7][12] =   0; AAmatrix[ 7][13] =   0; AAmatrix[ 7][14] =   0;
	AAmatrix[ 7][15] = 112; AAmatrix[ 7][16] =   0; AAmatrix[ 7][17] =   0; AAmatrix[ 7][18] =   1; AAmatrix[ 7][19] =   5;
	AAmatrix[ 8][ 0] =   8; AAmatrix[ 8][ 1] = 232; AAmatrix[ 8][ 2] = 458; AAmatrix[ 8][ 3] =  11; AAmatrix[ 8][ 4] = 305;
	AAmatrix[ 8][ 5] = 550; AAmatrix[ 8][ 6] =  22; AAmatrix[ 8][ 7] =   0; AAmatrix[ 8][ 8] =   0; AAmatrix[ 8][ 9] =   0;
	AAmatrix[ 8][10] =  26; AAmatrix[ 8][11] =   0; AAmatrix[ 8][12] =   0; AAmatrix[ 8][13] =   0; AAmatrix[ 8][14] =  53;
	AAmatrix[ 8][15] =  20; AAmatrix[ 8][16] =   1; AAmatrix[ 8][17] =   0; AAmatrix[ 8][18] =1525; AAmatrix[ 8][19] =   0;
	AAmatrix[ 9][ 0] =  75; AAmatrix[ 9][ 1] =   0; AAmatrix[ 9][ 2] =  19; AAmatrix[ 9][ 3] =   0; AAmatrix[ 9][ 4] =  41;
	AAmatrix[ 9][ 5] =   0; AAmatrix[ 9][ 6] =   0; AAmatrix[ 9][ 7] =   0; AAmatrix[ 9][ 8] =   0; AAmatrix[ 9][ 9] =   0;
	AAmatrix[ 9][10] = 232; AAmatrix[ 9][11] =   6; AAmatrix[ 9][12] = 378; AAmatrix[ 9][13] =  57; AAmatrix[ 9][14] =   5;
	AAmatrix[ 9][15] =   0; AAmatrix[ 9][16] = 360; AAmatrix[ 9][17] =   0; AAmatrix[ 9][18] =  16; AAmatrix[ 9][19] =2220;
	AAmatrix[10][ 0] =  21; AAmatrix[10][ 1] =   6; AAmatrix[10][ 2] =   0; AAmatrix[10][ 3] =   0; AAmatrix[10][ 4] =  27;
	AAmatrix[10][ 5] =  20; AAmatrix[10][ 6] =   0; AAmatrix[10][ 7] =   0; AAmatrix[10][ 8] =  26; AAmatrix[10][ 9] = 232;
	AAmatrix[10][10] =   0; AAmatrix[10][11] =   4; AAmatrix[10][12] = 609; AAmatrix[10][13] = 246; AAmatrix[10][14] =  43;
	AAmatrix[10][15] =  74; AAmatrix[10][16] =  34; AAmatrix[10][17] =  12; AAmatrix[10][18] =  25; AAmatrix[10][19] = 100;
	AAmatrix[11][ 0] =   0; AAmatrix[11][ 1] =  50; AAmatrix[11][ 2] = 408; AAmatrix[11][ 3] =   0; AAmatrix[11][ 4] =   0;
	AAmatrix[11][ 5] = 242; AAmatrix[11][ 6] = 215; AAmatrix[11][ 7] =   0; AAmatrix[11][ 8] =   0; AAmatrix[11][ 9] =   6;
	AAmatrix[11][10] =   4; AAmatrix[11][11] =   0; AAmatrix[11][12] =  59; AAmatrix[11][13] =   0; AAmatrix[11][14] =  18;
	AAmatrix[11][15] =  65; AAmatrix[11][16] =  50; AAmatrix[11][17] =   0; AAmatrix[11][18] =  67; AAmatrix[11][19] =   0;
	AAmatrix[12][ 0] =  76; AAmatrix[12][ 1] =   0; AAmatrix[12][ 2] =  21; AAmatrix[12][ 3] =   0; AAmatrix[12][ 4] =   0;
	AAmatrix[12][ 5] =  22; AAmatrix[12][ 6] =   0; AAmatrix[12][ 7] =   0; AAmatrix[12][ 8] =   0; AAmatrix[12][ 9] = 378;
	AAmatrix[12][10] = 609; AAmatrix[12][11] =  59; AAmatrix[12][12] =   0; AAmatrix[12][13] =  11; AAmatrix[12][14] =   0;
	AAmatrix[12][15] =  47; AAmatrix[12][16] = 691; AAmatrix[12][17] =  13; AAmatrix[12][18] =   0; AAmatrix[12][19] = 832;
	AAmatrix[13][ 0] =   0; AAmatrix[13][ 1] =   0; AAmatrix[13][ 2] =   6; AAmatrix[13][ 3] =   5; AAmatrix[13][ 4] =   7;
	AAmatrix[13][ 5] =   0; AAmatrix[13][ 6] =   0; AAmatrix[13][ 7] =   0; AAmatrix[13][ 8] =   0; AAmatrix[13][ 9] =  57;
	AAmatrix[13][10] = 246; AAmatrix[13][11] =   0; AAmatrix[13][12] =  11; AAmatrix[13][13] =   0; AAmatrix[13][14] =  17;
	AAmatrix[13][15] =  90; AAmatrix[13][16] =   8; AAmatrix[13][17] =   0; AAmatrix[13][18] = 682; AAmatrix[13][19] =   6;
	AAmatrix[14][ 0] =  53; AAmatrix[14][ 1] =   9; AAmatrix[14][ 2] =  33; AAmatrix[14][ 3] =   2; AAmatrix[14][ 4] =   0;
	AAmatrix[14][ 5] =  51; AAmatrix[14][ 6] =   0; AAmatrix[14][ 7] =   0; AAmatrix[14][ 8] =  53; AAmatrix[14][ 9] =   5;
	AAmatrix[14][10] =  43; AAmatrix[14][11] =  18; AAmatrix[14][12] =   0; AAmatrix[14][13] =  17; AAmatrix[14][14] =   0;
	AAmatrix[14][15] = 202; AAmatrix[14][16] =  78; AAmatrix[14][17] =   7; AAmatrix[14][18] =   8; AAmatrix[14][19] =   0;
	AAmatrix[15][ 0] = 342; AAmatrix[15][ 1] =   3; AAmatrix[15][ 2] = 446; AAmatrix[15][ 3] =  16; AAmatrix[15][ 4] = 347;
	AAmatrix[15][ 5] =  30; AAmatrix[15][ 6] =  21; AAmatrix[15][ 7] = 112; AAmatrix[15][ 8] =  20; AAmatrix[15][ 9] =   0;
	AAmatrix[15][10] =  74; AAmatrix[15][11] =  65; AAmatrix[15][12] =  47; AAmatrix[15][13] =  90; AAmatrix[15][14] = 202;
	AAmatrix[15][15] =   0; AAmatrix[15][16] = 614; AAmatrix[15][17] =  17; AAmatrix[15][18] = 107; AAmatrix[15][19] =   0;
	AAmatrix[16][ 0] = 681; AAmatrix[16][ 1] =   0; AAmatrix[16][ 2] = 110; AAmatrix[16][ 3] =   0; AAmatrix[16][ 4] = 114;
	AAmatrix[16][ 5] =   0; AAmatrix[16][ 6] =   4; AAmatrix[16][ 7] =   0; AAmatrix[16][ 8] =   1; AAmatrix[16][ 9] = 360;
	AAmatrix[16][10] =  34; AAmatrix[16][11] =  50; AAmatrix[16][12] = 691; AAmatrix[16][13] =   8; AAmatrix[16][14] =  78;
	AAmatrix[16][15] = 614; AAmatrix[16][16] =   0; AAmatrix[16][17] =   0; AAmatrix[16][18] =   0; AAmatrix[16][19] = 237;
	AAmatrix[17][ 0] =   5; AAmatrix[17][ 1] =  16; AAmatrix[17][ 2] =   6; AAmatrix[17][ 3] =   0; AAmatrix[17][ 4] =  65;
	AAmatrix[17][ 5] =   0; AAmatrix[17][ 6] =   0; AAmatrix[17][ 7] =   0; AAmatrix[17][ 8] =   0; AAmatrix[17][ 9] =   0;
	AAmatrix[17][10] =  12; AAmatrix[17][11] =   0; AAmatrix[17][12] =  13; AAmatrix[17][13] =   0; AAmatrix[17][14] =   7;
	AAmatrix[17][15] =  17; AAmatrix[17][16] =   0; AAmatrix[17][17] =   0; AAmatrix[17][18] =  14; AAmatrix[17][19] =   0;
	AAmatrix[18][ 0] =   0; AAmatrix[18][ 1] =   0; AAmatrix[18][ 2] = 156; AAmatrix[18][ 3] =   0; AAmatrix[18][ 4] = 530;
	AAmatrix[18][ 5] =  54; AAmatrix[18][ 6] =   0; AAmatrix[18][ 7] =   1; AAmatrix[18][ 8] =1525; AAmatrix[18][ 9] =  16;
	AAmatrix[18][10] =  25; AAmatrix[18][11] =  67; AAmatrix[18][12] =   0; AAmatrix[18][13] = 682; AAmatrix[18][14] =   8;
	AAmatrix[18][15] = 107; AAmatrix[18][16] =   0; AAmatrix[18][17] =  14; AAmatrix[18][18] =   0; AAmatrix[18][19] =   0;
	AAmatrix[19][ 0] = 398; AAmatrix[19][ 1] =   0; AAmatrix[19][ 2] =   0; AAmatrix[19][ 3] =  10; AAmatrix[19][ 4] =   0;
	AAmatrix[19][ 5] =  33; AAmatrix[19][ 6] =  20; AAmatrix[19][ 7] =   5; AAmatrix[19][ 8] =   0; AAmatrix[19][ 9] =2220;
	AAmatrix[19][10] = 100; AAmatrix[19][11] =   0; AAmatrix[19][12] = 832; AAmatrix[19][13] =   6; AAmatrix[19][14] =   0;
	AAmatrix[19][15] =   0; AAmatrix[19][16] = 237; AAmatrix[19][17] =   0; AAmatrix[19][18] =   0; AAmatrix[19][19] =   0;

	AAPi[ 0] = 0.0692;
	AAPi[ 1] = 0.0184;
	AAPi[ 2] = 0.0400;
	AAPi[ 3] = 0.0186;
	AAPi[ 4] = 0.0065;
	AAPi[ 5] = 0.0238;
	AAPi[ 6] = 0.0236;
	AAPi[ 7] = 0.0557;
	AAPi[ 8] = 0.0277;
	AAPi[ 9] = 0.0905;
	AAPi[10] = 0.1675;
	AAPi[11] = 0.0221;
	AAPi[12] = 0.0561;
	AAPi[13] = 0.0611;
	AAPi[14] = 0.0536;
	AAPi[15] = 0.0725;
	AAPi[16] = 0.0870;
	AAPi[17] = 0.0293;
	AAPi[18] = 0.0340;
	AAPi[19] = 0.0428;
	
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==7)
	{

	/*
		Mtart
		
		Abascal, F., D. Posada, and R. Zardoya. (2007) MtArt: A new Model of amino acid replacement for Arthropoda. Mol. Biol. Evol. 24:1-5. 
	*/
	AAmatrix[0][0] = 0;			AAmatrix[0][1] = 0.2;		AAmatrix[0][2] = 0.2;		AAmatrix[0][3] = 0.6;	
	AAmatrix[0][4] = 253.5;		AAmatrix[0][5] = 0.2;		AAmatrix[0][6] = 0.2;		AAmatrix[0][7] = 199.8;	
	AAmatrix[0][8] = 0.2;		AAmatrix[0][9] = 25.7;		AAmatrix[0][10] = 3.7;		AAmatrix[0][11] = 0.2;	
	AAmatrix[0][12] = 120.6;	AAmatrix[0][13] = 13.1;		AAmatrix[0][14] = 49.3;		AAmatrix[0][15] = 673;	
	AAmatrix[0][16] = 243.9;	AAmatrix[0][17] = 0.2;		AAmatrix[0][18] = 1.2;		AAmatrix[0][19] = 339.9;	
	AAmatrix[1][0] = 0.2;		AAmatrix[1][1] = 0;			AAmatrix[1][2] = 0.2;		AAmatrix[1][3] = 4.3;	
	AAmatrix[1][4] = 35.5;		AAmatrix[1][5] = 154;		AAmatrix[1][6] = 0.2;		AAmatrix[1][7] = 0.2;	
	AAmatrix[1][8] = 41.3;		AAmatrix[1][9] = 1.8;		AAmatrix[1][10] = 1.8;		AAmatrix[1][11] = 208.6;	
	AAmatrix[1][12] = 5.2;		AAmatrix[1][13] = 4.7;		AAmatrix[1][14] = 0.2;		AAmatrix[1][15] = 2.7;	
	AAmatrix[1][16] = 0.2;		AAmatrix[1][17] = 0.2;		AAmatrix[1][18] = 3.9;		AAmatrix[1][19] = 0.2;	
	AAmatrix[2][0] = 0.2;		AAmatrix[2][1] = 0.2;		AAmatrix[2][2] = 0;			AAmatrix[2][3] = 500.2;	
	AAmatrix[2][4] = 98.2;		AAmatrix[2][5] = 261.8;		AAmatrix[2][6] = 183;		AAmatrix[2][7] = 120.5;	
	AAmatrix[2][8] = 179.5;		AAmatrix[2][9] = 21.3;		AAmatrix[2][10] = 12.6;		AAmatrix[2][11] = 467.3;	
	AAmatrix[2][12] = 78.8;		AAmatrix[2][13] = 19.7;		AAmatrix[2][14] = 16.5;		AAmatrix[2][15] = 398.4;	
	AAmatrix[2][16] = 165.9;	AAmatrix[2][17] = 7.7;		AAmatrix[2][18] = 251.2;	AAmatrix[2][19] = 22.6;	
	AAmatrix[3][0] = 0.6;		AAmatrix[3][1] = 4.3;		AAmatrix[3][2] = 500.2;		AAmatrix[3][3] = 0;	
	AAmatrix[3][4] = 10.6;		AAmatrix[3][5] = 0.2;		AAmatrix[3][6] = 861.8;		AAmatrix[3][7] = 12.5;	
	AAmatrix[3][8] = 0.2;		AAmatrix[3][9] = 6.6;		AAmatrix[3][10] = 1.2;		AAmatrix[3][11] = 1.7;	
	AAmatrix[3][12] = 0.2;		AAmatrix[3][13] = 0.2;		AAmatrix[3][14] = 0.2;		AAmatrix[3][15] = 44.4;	
	AAmatrix[3][16] = 0.2;		AAmatrix[3][17] = 0.2;		AAmatrix[3][18] = 0.2;		AAmatrix[3][19] = 0.2;	
	AAmatrix[4][0] = 253.5;		AAmatrix[4][1] = 35.5;		AAmatrix[4][2] = 98.2;		AAmatrix[4][3] = 10.6;	
	AAmatrix[4][4] = 0;			AAmatrix[4][5] = 0.2;		AAmatrix[4][6] = 0.2;		AAmatrix[4][7] = 80.5;	
	AAmatrix[4][8] = 12.4;		AAmatrix[4][9] = 63;		AAmatrix[4][10] = 78.7;		AAmatrix[4][11] = 0.2;	
	AAmatrix[4][12] = 312.3;	AAmatrix[4][13] = 184.1;	AAmatrix[4][14] = 0.2;		AAmatrix[4][15] = 664.2;	
	AAmatrix[4][16] = 182.8;	AAmatrix[4][17] = 21.6;		AAmatrix[4][18] = 72;		AAmatrix[4][19] = 350.4;	
	AAmatrix[5][0] = 0.2;		AAmatrix[5][1] = 154;		AAmatrix[5][2] = 261.8;		AAmatrix[5][3] = 0.2;	
	AAmatrix[5][4] = 0.2;		AAmatrix[5][5] = 0;			AAmatrix[5][6] = 261.6;		AAmatrix[5][7] = 2.6;	
	AAmatrix[5][8] = 313.5;		AAmatrix[5][9] = 10.5;		AAmatrix[5][10] = 16.3;		AAmatrix[5][11] = 349.3;	
	AAmatrix[5][12] = 67.3;		AAmatrix[5][13] = 0.2;		AAmatrix[5][14] = 39.3;		AAmatrix[5][15] = 52.4;	
	AAmatrix[5][16] = 43.7;		AAmatrix[5][17] = 6.7;		AAmatrix[5][18] = 86.7;		AAmatrix[5][19] = 0.2;	
	AAmatrix[6][0] = 0.2;		AAmatrix[6][1] = 0.2;		AAmatrix[6][2] = 183;		AAmatrix[6][3] = 861.8;	
	AAmatrix[6][4] = 0.2;		AAmatrix[6][5] = 261.6;		AAmatrix[6][6] = 0;			AAmatrix[6][7] = 43.9;	
	AAmatrix[6][8] = 15.2;		AAmatrix[6][9] = 6.8;		AAmatrix[6][10] = 1.7;		AAmatrix[6][11] = 106.3;	
	AAmatrix[6][12] = 0.2;		AAmatrix[6][13] = 0.2;		AAmatrix[6][14] = 7.9;		AAmatrix[6][15] = 31.5;	
	AAmatrix[6][16] = 43.4;		AAmatrix[6][17] = 11;		AAmatrix[6][18] = 7.7;		AAmatrix[6][19] = 13.6;	
	AAmatrix[7][0] = 199.8;		AAmatrix[7][1] = 0.2;		AAmatrix[7][2] = 120.5;		AAmatrix[7][3] = 12.5;	
	AAmatrix[7][4] = 80.5;		AAmatrix[7][5] = 2.6;		AAmatrix[7][6] = 43.9;		AAmatrix[7][7] = 0;	
	AAmatrix[7][8] = 0.2;		AAmatrix[7][9] = 2.7;		AAmatrix[7][10] = 1.4;		AAmatrix[7][11] = 0.2;	
	AAmatrix[7][12] = 55.7;		AAmatrix[7][13] = 0.8;		AAmatrix[7][14] = 0.2;		AAmatrix[7][15] = 226;	
	AAmatrix[7][16] = 0.2;		AAmatrix[7][17] = 1.9;		AAmatrix[7][18] = 8.6;		AAmatrix[7][19] = 2.6;	
	AAmatrix[8][0] = 0.2;		AAmatrix[8][1] = 41.3;		AAmatrix[8][2] = 179.5;		AAmatrix[8][3] = 0.2;	
	AAmatrix[8][4] = 12.4;		AAmatrix[8][5] = 313.5;		AAmatrix[8][6] = 15.2;		AAmatrix[8][7] = 0.2;	
	AAmatrix[8][8] = 0;			AAmatrix[8][9] = 0.2;		AAmatrix[8][10] = 5.5;		AAmatrix[8][11] = 0.2;	
	AAmatrix[8][12] = 0.2;		AAmatrix[8][13] = 13.8;		AAmatrix[8][14] = 0.8;		AAmatrix[8][15] = 10.6;	
	AAmatrix[8][16] = 18.6;		AAmatrix[8][17] = 0.2;		AAmatrix[8][18] = 191.4;	AAmatrix[8][19] = 0.2;	
	AAmatrix[9][0] = 25.7;		AAmatrix[9][1] = 1.8;		AAmatrix[9][2] = 21.3;		AAmatrix[9][3] = 6.6;	
	AAmatrix[9][4] = 63;		AAmatrix[9][5] = 10.5;		AAmatrix[9][6] = 6.8;		AAmatrix[9][7] = 2.7;	
	AAmatrix[9][8] = 0.2;		AAmatrix[9][9] = 0;			AAmatrix[9][10] = 514.5;	AAmatrix[9][11] = 3.5;	
	AAmatrix[9][12] = 514.8;	AAmatrix[9][13] = 117.9;	AAmatrix[9][14] = 0.2;		AAmatrix[9][15] = 7.2;	
	AAmatrix[9][16] = 203.7;	AAmatrix[9][17] = 0.2;		AAmatrix[9][18] = 12.3;		AAmatrix[9][19] = 1854.5;	
	AAmatrix[10][0] = 3.7;		AAmatrix[10][1] = 1.8;		AAmatrix[10][2] = 12.6;		AAmatrix[10][3] = 1.2;	
	AAmatrix[10][4] = 78.7;		AAmatrix[10][5] = 16.3;		AAmatrix[10][6] = 1.7;		AAmatrix[10][7] = 1.4;	
	AAmatrix[10][8] = 5.5;		AAmatrix[10][9] = 514.5;	AAmatrix[10][10] = 0;		AAmatrix[10][11] = 3.8;	
	AAmatrix[10][12] = 885.5;	AAmatrix[10][13] = 262.6;	AAmatrix[10][14] = 12.2;	AAmatrix[10][15] = 8.2;	
	AAmatrix[10][16] = 47.8;	AAmatrix[10][17] = 21.1;	AAmatrix[10][18] = 19.8;	AAmatrix[10][19] = 84.7;	
	AAmatrix[11][0] = 0.2;		AAmatrix[11][1] = 208.6;	AAmatrix[11][2] = 467.3;	AAmatrix[11][3] = 1.7;	
	AAmatrix[11][4] = 0.2;		AAmatrix[11][5] = 349.3;	AAmatrix[11][6] = 106.3;	AAmatrix[11][7] = 0.2;	
	AAmatrix[11][8] = 0.2;		AAmatrix[11][9] = 3.5;		AAmatrix[11][10] = 3.8;		AAmatrix[11][11] = 0;	
	AAmatrix[11][12] = 105.6;	AAmatrix[11][13] = 10.7;	AAmatrix[11][14] = 16.8;	AAmatrix[11][15] = 144.2;	
	AAmatrix[11][16] = 69.5;	AAmatrix[11][17] = 16;		AAmatrix[11][18] = 117.1;	AAmatrix[11][19] = 26.1;	
	AAmatrix[12][0] = 120.6;	AAmatrix[12][1] = 5.2;		AAmatrix[12][2] = 78.8;		AAmatrix[12][3] = 0.2;	
	AAmatrix[12][4] = 312.3;	AAmatrix[12][5] = 67.3;		AAmatrix[12][6] = 0.2;		AAmatrix[12][7] = 55.7;	
	AAmatrix[12][8] = 0.2;		AAmatrix[12][9] = 514.8;	AAmatrix[12][10] = 885.5;	AAmatrix[12][11] = 105.6;	
	AAmatrix[12][12] = 0;		AAmatrix[12][13] = 321.6;	AAmatrix[12][14] = 5.3;		AAmatrix[12][15] = 111.7;	
	AAmatrix[12][16] = 288.6;	AAmatrix[12][17] = 70.7;	AAmatrix[12][18] = 70.9;	AAmatrix[12][19] = 281.3;	
	AAmatrix[13][0] = 13.1;		AAmatrix[13][1] = 4.7;		AAmatrix[13][2] = 19.7;		AAmatrix[13][3] = 0.2;	
	AAmatrix[13][4] = 184.1;	AAmatrix[13][5] = 0.2;		AAmatrix[13][6] = 0.2;		AAmatrix[13][7] = 0.8;	
	AAmatrix[13][8] = 13.8;		AAmatrix[13][9] = 117.9;	AAmatrix[13][10] = 262.6;	AAmatrix[13][11] = 10.7;	
	AAmatrix[13][12] = 321.6;	AAmatrix[13][13] = 0;		AAmatrix[13][14] = 14.6;	AAmatrix[13][15] = 36.1;	
	AAmatrix[13][16] = 13.5;	AAmatrix[13][17] = 53.7;	AAmatrix[13][18] = 791.6;	AAmatrix[13][19] = 51.9;	
	AAmatrix[14][0] = 49.3;		AAmatrix[14][1] = 0.2;		AAmatrix[14][2] = 16.5;		AAmatrix[14][3] = 0.2;	
	AAmatrix[14][4] = 0.2;		AAmatrix[14][5] = 39.3;		AAmatrix[14][6] = 7.9;		AAmatrix[14][7] = 0.2;	
	AAmatrix[14][8] = 0.8;		AAmatrix[14][9] = 0.2;		AAmatrix[14][10] = 12.2;	AAmatrix[14][11] = 16.8;	
	AAmatrix[14][12] = 5.3;		AAmatrix[14][13] = 14.6;	AAmatrix[14][14] = 0;		AAmatrix[14][15] = 86.5;	
	AAmatrix[14][16] = 46.8;	AAmatrix[14][17] = 0.2;		AAmatrix[14][18] = 18.4;	AAmatrix[14][19] = 31.7;	
	AAmatrix[15][0] = 673;		AAmatrix[15][1] = 2.7;		AAmatrix[15][2] = 398.4;	AAmatrix[15][3] = 44.4;	
	AAmatrix[15][4] = 664.2;	AAmatrix[15][5] = 52.4;		AAmatrix[15][6] = 31.5;		AAmatrix[15][7] = 226;	
	AAmatrix[15][8] = 10.6;		AAmatrix[15][9] = 7.2;		AAmatrix[15][10] = 8.2;		AAmatrix[15][11] = 144.2;	
	AAmatrix[15][12] = 111.7;	AAmatrix[15][13] = 36.1;	AAmatrix[15][14] = 86.5;	AAmatrix[15][15] = 0;	
	AAmatrix[15][16] = 660.4;	AAmatrix[15][17] = 2.4;		AAmatrix[15][18] = 30.5;	AAmatrix[15][19] = 60.6;	
	AAmatrix[16][0] = 243.9;	AAmatrix[16][1] = 0.2;		AAmatrix[16][2] = 165.9;	AAmatrix[16][3] = 0.2;	
	AAmatrix[16][4] = 182.8;	AAmatrix[16][5] = 43.7;		AAmatrix[16][6] = 43.4;		AAmatrix[16][7] = 0.2;	
	AAmatrix[16][8] = 18.6;		AAmatrix[16][9] = 203.7;	AAmatrix[16][10] = 47.8;	AAmatrix[16][11] = 69.5;	
	AAmatrix[16][12] = 288.6;	AAmatrix[16][13] = 13.5;	AAmatrix[16][14] = 46.8;	AAmatrix[16][15] = 660.4;	
	AAmatrix[16][16] = 0;		AAmatrix[16][17] = 0.2;		AAmatrix[16][18] = 46;		AAmatrix[16][19] = 544.1;	
	AAmatrix[17][0] = 0.2;		AAmatrix[17][1] = 0.2;		AAmatrix[17][2] = 7.7;		AAmatrix[17][3] = 0.2;	
	AAmatrix[17][4] = 21.6;		AAmatrix[17][5] = 6.7;		AAmatrix[17][6] = 11;		AAmatrix[17][7] = 1.9;	
	AAmatrix[17][8] = 0.2;		AAmatrix[17][9] = 0.2;		AAmatrix[17][10] = 21.1;	AAmatrix[17][11] = 16;	
	AAmatrix[17][12] = 70.7;	AAmatrix[17][13] = 53.7;	AAmatrix[17][14] = 0.2;		AAmatrix[17][15] = 2.4;	
	AAmatrix[17][16] = 0.2;		AAmatrix[17][17] = 0;		AAmatrix[17][18] = 37.7;	AAmatrix[17][19] = 0.2;	
	AAmatrix[18][0] = 1.2;		AAmatrix[18][1] = 3.9;		AAmatrix[18][2] = 251.2;	AAmatrix[18][3] = 0.2;	
	AAmatrix[18][4] = 72;		AAmatrix[18][5] = 86.7;		AAmatrix[18][6] = 7.7;		AAmatrix[18][7] = 8.6;	
	AAmatrix[18][8] = 191.4;	AAmatrix[18][9] = 12.3;		AAmatrix[18][10] = 19.8;	AAmatrix[18][11] = 117.1;	
	AAmatrix[18][12] = 70.9;	AAmatrix[18][13] = 791.6;	AAmatrix[18][14] = 18.4;	AAmatrix[18][15] = 30.5;	
	AAmatrix[18][16] = 46;		AAmatrix[18][17] = 37.7;	AAmatrix[18][18] = 0;		AAmatrix[18][19] = 1.6;	
	AAmatrix[19][0] = 339.9;	AAmatrix[19][1] = 0.2;		AAmatrix[19][2] = 22.6;		AAmatrix[19][3] = 0.2;	
	AAmatrix[19][4] = 350.4;	AAmatrix[19][5] = 0.2;		AAmatrix[19][6] = 13.6;		AAmatrix[19][7] = 2.6;	
	AAmatrix[19][8] = 0.2;		AAmatrix[19][9] = 1854.5;	AAmatrix[19][10] = 84.7;	AAmatrix[19][11] = 26.1;	
	AAmatrix[19][12] = 281.3;	AAmatrix[19][13] = 51.9;	AAmatrix[19][14] = 31.7;	AAmatrix[19][15] = 60.6;	
	AAmatrix[19][16] = 544.1;	AAmatrix[19][17] = 0.2;		AAmatrix[19][18] = 1.6;		AAmatrix[19][19] = 0;	

	AAPi[0] = 0.054116;
	AAPi[1] = 0.018227;
	AAPi[2] = 0.039903;
	AAPi[3] = 0.02016;
	AAPi[4] = 0.009709;
	AAPi[5] = 0.018781;
	AAPi[6] = 0.024289;
	AAPi[7] = 0.068183;
	AAPi[8] = 0.024518;
	AAPi[9] = 0.092639;
	AAPi[10] = 0.148658;
	AAPi[11] = 0.021718;
	AAPi[12] = 0.061453;
	AAPi[13] = 0.088668;
	AAPi[14] = 0.041826;
	AAPi[15] = 0.09103;
	AAPi[16] = 0.049194;
	AAPi[17] = 0.029786;
	AAPi[18] = 0.039443;
	AAPi[19] = 0.057701;


	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else
	if(modelnumber==8)
	{

	/* 
		mtrev24 
		
		Adachi, J. and Hasegawa, M. (1996) MOLPHY version 2.3: programs for molecular phylogenetics based on maximum likelihood.  
		Computer Science Monographs of Institute of Statistical Mathematics 28:1-150.

	*/
	AAmatrix[ 0][ 0] =   0.00; AAmatrix[ 0][ 1] =  23.18; AAmatrix[ 0][ 2] =  26.95; AAmatrix[ 0][ 3] =  17.67; AAmatrix[ 0][ 4] =  59.93;
	AAmatrix[ 0][ 5] =   1.90; AAmatrix[ 0][ 6] =   9.77; AAmatrix[ 0][ 7] = 120.71; AAmatrix[ 0][ 8] =  13.90; AAmatrix[ 0][ 9] =  96.49;
	AAmatrix[ 0][10] =  25.46; AAmatrix[ 0][11] =   8.36; AAmatrix[ 0][12] = 141.88; AAmatrix[ 0][13] =   6.37; AAmatrix[ 0][14] =  54.31;
	AAmatrix[ 0][15] = 387.86; AAmatrix[ 0][16] = 480.72; AAmatrix[ 0][17] =   1.90; AAmatrix[ 0][18] =   6.48; AAmatrix[ 0][19] = 195.06;
	AAmatrix[ 1][ 0] =  23.18; AAmatrix[ 1][ 1] =   0.00; AAmatrix[ 1][ 2] =  13.24; AAmatrix[ 1][ 3] =   1.90; AAmatrix[ 1][ 4] = 103.33;
	AAmatrix[ 1][ 5] = 220.99; AAmatrix[ 1][ 6] =   1.90; AAmatrix[ 1][ 7] =  23.03; AAmatrix[ 1][ 8] = 165.23; AAmatrix[ 1][ 9] =   1.90;
	AAmatrix[ 1][10] =  15.58; AAmatrix[ 1][11] = 141.40; AAmatrix[ 1][12] =   1.90; AAmatrix[ 1][13] =   4.69; AAmatrix[ 1][14] =  23.64;
	AAmatrix[ 1][15] =   6.04; AAmatrix[ 1][16] =   2.08; AAmatrix[ 1][17] =  21.95; AAmatrix[ 1][18] =   1.90; AAmatrix[ 1][19] =   7.64;
	AAmatrix[ 2][ 0] =  26.95; AAmatrix[ 2][ 1] =  13.24; AAmatrix[ 2][ 2] =   0.00; AAmatrix[ 2][ 3] = 794.38; AAmatrix[ 2][ 4] =  58.94;
	AAmatrix[ 2][ 5] = 173.56; AAmatrix[ 2][ 6] =  63.05; AAmatrix[ 2][ 7] =  53.30; AAmatrix[ 2][ 8] = 496.13; AAmatrix[ 2][ 9] =  27.10;
	AAmatrix[ 2][10] =  15.16; AAmatrix[ 2][11] = 608.70; AAmatrix[ 2][12] =  65.41; AAmatrix[ 2][13] =  15.20; AAmatrix[ 2][14] =  73.31;
	AAmatrix[ 2][15] = 494.39; AAmatrix[ 2][16] = 238.46; AAmatrix[ 2][17] =  10.68; AAmatrix[ 2][18] = 191.36; AAmatrix[ 2][19] =   1.90;
	AAmatrix[ 3][ 0] =  17.67; AAmatrix[ 3][ 1] =   1.90; AAmatrix[ 3][ 2] = 794.38; AAmatrix[ 3][ 3] =   0.00; AAmatrix[ 3][ 4] =   1.90;
	AAmatrix[ 3][ 5] =  55.28; AAmatrix[ 3][ 6] = 583.55; AAmatrix[ 3][ 7] =  56.77; AAmatrix[ 3][ 8] = 113.99; AAmatrix[ 3][ 9] =   4.34;
	AAmatrix[ 3][10] =   1.90; AAmatrix[ 3][11] =   2.31; AAmatrix[ 3][12] =   1.90; AAmatrix[ 3][13] =   4.98; AAmatrix[ 3][14] =  13.43;
	AAmatrix[ 3][15] =  69.02; AAmatrix[ 3][16] =  28.01; AAmatrix[ 3][17] =  19.86; AAmatrix[ 3][18] =  21.21; AAmatrix[ 3][19] =   1.90;
	AAmatrix[ 4][ 0] =  59.93; AAmatrix[ 4][ 1] = 103.33; AAmatrix[ 4][ 2] =  58.94; AAmatrix[ 4][ 3] =   1.90; AAmatrix[ 4][ 4] =   0.00;
	AAmatrix[ 4][ 5] =  75.24; AAmatrix[ 4][ 6] =   1.90; AAmatrix[ 4][ 7] =  30.71; AAmatrix[ 4][ 8] = 141.49; AAmatrix[ 4][ 9] =  62.73;
	AAmatrix[ 4][10] =  25.65; AAmatrix[ 4][11] =   1.90; AAmatrix[ 4][12] =   6.18; AAmatrix[ 4][13] =  70.80; AAmatrix[ 4][14] =  31.26;
	AAmatrix[ 4][15] = 277.05; AAmatrix[ 4][16] = 179.97; AAmatrix[ 4][17] =  33.60; AAmatrix[ 4][18] = 254.77; AAmatrix[ 4][19] =   1.90;
	AAmatrix[ 5][ 0] =   1.90; AAmatrix[ 5][ 1] = 220.99; AAmatrix[ 5][ 2] = 173.56; AAmatrix[ 5][ 3] =  55.28; AAmatrix[ 5][ 4] =  75.24;
	AAmatrix[ 5][ 5] =   0.00; AAmatrix[ 5][ 6] = 313.56; AAmatrix[ 5][ 7] =   6.75; AAmatrix[ 5][ 8] = 582.40; AAmatrix[ 5][ 9] =   8.34;
	AAmatrix[ 5][10] =  39.70; AAmatrix[ 5][11] = 465.58; AAmatrix[ 5][12] =  47.37; AAmatrix[ 5][13] =  19.11; AAmatrix[ 5][14] = 137.29;
	AAmatrix[ 5][15] =  54.11; AAmatrix[ 5][16] =  94.93; AAmatrix[ 5][17] =   1.90; AAmatrix[ 5][18] =  38.82; AAmatrix[ 5][19] =  19.00;
	AAmatrix[ 6][ 0] =   9.77; AAmatrix[ 6][ 1] =   1.90; AAmatrix[ 6][ 2] =  63.05; AAmatrix[ 6][ 3] = 583.55; AAmatrix[ 6][ 4] =   1.90;
	AAmatrix[ 6][ 5] = 313.56; AAmatrix[ 6][ 6] =   0.00; AAmatrix[ 6][ 7] =  28.28; AAmatrix[ 6][ 8] =  49.12; AAmatrix[ 6][ 9] =   3.31;
	AAmatrix[ 6][10] =   1.90; AAmatrix[ 6][11] = 313.86; AAmatrix[ 6][12] =   1.90; AAmatrix[ 6][13] =   2.67; AAmatrix[ 6][14] =  12.83;
	AAmatrix[ 6][15] =  54.71; AAmatrix[ 6][16] =  14.82; AAmatrix[ 6][17] =   1.90; AAmatrix[ 6][18] =  13.12; AAmatrix[ 6][19] =  21.14;
	AAmatrix[ 7][ 0] = 120.71; AAmatrix[ 7][ 1] =  23.03; AAmatrix[ 7][ 2] =  53.30; AAmatrix[ 7][ 3] =  56.77; AAmatrix[ 7][ 4] =  30.71;
	AAmatrix[ 7][ 5] =   6.75; AAmatrix[ 7][ 6] =  28.28; AAmatrix[ 7][ 7] =   0.00; AAmatrix[ 7][ 8] =   1.90; AAmatrix[ 7][ 9] =   5.98;
	AAmatrix[ 7][10] =   2.41; AAmatrix[ 7][11] =  22.73; AAmatrix[ 7][12] =   1.90; AAmatrix[ 7][13] =   1.90; AAmatrix[ 7][14] =   1.90;
	AAmatrix[ 7][15] = 125.93; AAmatrix[ 7][16] =  11.17; AAmatrix[ 7][17] =  10.92; AAmatrix[ 7][18] =   3.21; AAmatrix[ 7][19] =   2.53;
	AAmatrix[ 8][ 0] =  13.90; AAmatrix[ 8][ 1] = 165.23; AAmatrix[ 8][ 2] = 496.13; AAmatrix[ 8][ 3] = 113.99; AAmatrix[ 8][ 4] = 141.49;
	AAmatrix[ 8][ 5] = 582.40; AAmatrix[ 8][ 6] =  49.12; AAmatrix[ 8][ 7] =   1.90; AAmatrix[ 8][ 8] =   0.00; AAmatrix[ 8][ 9] =  12.26;
	AAmatrix[ 8][10] =  11.49; AAmatrix[ 8][11] = 127.67; AAmatrix[ 8][12] =  11.97; AAmatrix[ 8][13] =  48.16; AAmatrix[ 8][14] =  60.97;
	AAmatrix[ 8][15] =  77.46; AAmatrix[ 8][16] =  44.78; AAmatrix[ 8][17] =   7.08; AAmatrix[ 8][18] = 670.14; AAmatrix[ 8][19] =   1.90;
	AAmatrix[ 9][ 0] =  96.49; AAmatrix[ 9][ 1] =   1.90; AAmatrix[ 9][ 2] =  27.10; AAmatrix[ 9][ 3] =   4.34; AAmatrix[ 9][ 4] =  62.73;
	AAmatrix[ 9][ 5] =   8.34; AAmatrix[ 9][ 6] =   3.31; AAmatrix[ 9][ 7] =   5.98; AAmatrix[ 9][ 8] =  12.26; AAmatrix[ 9][ 9] =   0.00;
	AAmatrix[ 9][10] = 329.09; AAmatrix[ 9][11] =  19.57; AAmatrix[ 9][12] = 517.98; AAmatrix[ 9][13] =  84.67; AAmatrix[ 9][14] =  20.63;
	AAmatrix[ 9][15] =  47.70; AAmatrix[ 9][16] = 368.43; AAmatrix[ 9][17] =   1.90; AAmatrix[ 9][18] =  25.01; AAmatrix[ 9][19] =1222.94;
	AAmatrix[10][ 0] =  25.46; AAmatrix[10][ 1] =  15.58; AAmatrix[10][ 2] =  15.16; AAmatrix[10][ 3] =   1.90; AAmatrix[10][ 4] =  25.65;
	AAmatrix[10][ 5] =  39.70; AAmatrix[10][ 6] =   1.90; AAmatrix[10][ 7] =   2.41; AAmatrix[10][ 8] =  11.49; AAmatrix[10][ 9] = 329.09;
	AAmatrix[10][10] =   0.00; AAmatrix[10][11] =  14.88; AAmatrix[10][12] = 537.53; AAmatrix[10][13] = 216.06; AAmatrix[10][14] =  40.10;
	AAmatrix[10][15] =  73.61; AAmatrix[10][16] = 126.40; AAmatrix[10][17] =  32.44; AAmatrix[10][18] =  44.15; AAmatrix[10][19] =  91.67;
	AAmatrix[11][ 0] =   8.36; AAmatrix[11][ 1] = 141.40; AAmatrix[11][ 2] = 608.70; AAmatrix[11][ 3] =   2.31; AAmatrix[11][ 4] =   1.90;
	AAmatrix[11][ 5] = 465.58; AAmatrix[11][ 6] = 313.86; AAmatrix[11][ 7] =  22.73; AAmatrix[11][ 8] = 127.67; AAmatrix[11][ 9] =  19.57;
	AAmatrix[11][10] =  14.88; AAmatrix[11][11] =   0.00; AAmatrix[11][12] =  91.37; AAmatrix[11][13] =   6.44; AAmatrix[11][14] =  50.10;
	AAmatrix[11][15] = 105.79; AAmatrix[11][16] = 136.33; AAmatrix[11][17] =  24.00; AAmatrix[11][18] =  51.17; AAmatrix[11][19] =   1.90;
	AAmatrix[12][ 0] = 141.88; AAmatrix[12][ 1] =   1.90; AAmatrix[12][ 2] =  65.41; AAmatrix[12][ 3] =   1.90; AAmatrix[12][ 4] =   6.18;
	AAmatrix[12][ 5] =  47.37; AAmatrix[12][ 6] =   1.90; AAmatrix[12][ 7] =   1.90; AAmatrix[12][ 8] =  11.97; AAmatrix[12][ 9] = 517.98;
	AAmatrix[12][10] = 537.53; AAmatrix[12][11] =  91.37; AAmatrix[12][12] =   0.00; AAmatrix[12][13] =  90.82; AAmatrix[12][14] =  18.84;
	AAmatrix[12][15] = 111.16; AAmatrix[12][16] = 528.17; AAmatrix[12][17] =  21.71; AAmatrix[12][18] =  39.96; AAmatrix[12][19] = 387.54;
	AAmatrix[13][ 0] =   6.37; AAmatrix[13][ 1] =   4.69; AAmatrix[13][ 2] =  15.20; AAmatrix[13][ 3] =   4.98; AAmatrix[13][ 4] =  70.80;
	AAmatrix[13][ 5] =  19.11; AAmatrix[13][ 6] =   2.67; AAmatrix[13][ 7] =   1.90; AAmatrix[13][ 8] =  48.16; AAmatrix[13][ 9] =  84.67;
	AAmatrix[13][10] = 216.06; AAmatrix[13][11] =   6.44; AAmatrix[13][12] =  90.82; AAmatrix[13][13] =   0.00; AAmatrix[13][14] =  17.31;
	AAmatrix[13][15] =  64.29; AAmatrix[13][16] =  33.85; AAmatrix[13][17] =   7.84; AAmatrix[13][18] = 465.58; AAmatrix[13][19] =   6.35;
	AAmatrix[14][ 0] =  54.31; AAmatrix[14][ 1] =  23.64; AAmatrix[14][ 2] =  73.31; AAmatrix[14][ 3] =  13.43; AAmatrix[14][ 4] =  31.26;
	AAmatrix[14][ 5] = 137.29; AAmatrix[14][ 6] =  12.83; AAmatrix[14][ 7] =   1.90; AAmatrix[14][ 8] =  60.97; AAmatrix[14][ 9] =  20.63;
	AAmatrix[14][10] =  40.10; AAmatrix[14][11] =  50.10; AAmatrix[14][12] =  18.84; AAmatrix[14][13] =  17.31; AAmatrix[14][14] =   0.00;
	AAmatrix[14][15] = 169.90; AAmatrix[14][16] = 128.22; AAmatrix[14][17] =   4.21; AAmatrix[14][18] =  16.21; AAmatrix[14][19] =   8.23;
	AAmatrix[15][ 0] = 387.86; AAmatrix[15][ 1] =   6.04; AAmatrix[15][ 2] = 494.39; AAmatrix[15][ 3] =  69.02; AAmatrix[15][ 4] = 277.05;
	AAmatrix[15][ 5] =  54.11; AAmatrix[15][ 6] =  54.71; AAmatrix[15][ 7] = 125.93; AAmatrix[15][ 8] =  77.46; AAmatrix[15][ 9] =  47.70;
	AAmatrix[15][10] =  73.61; AAmatrix[15][11] = 105.79; AAmatrix[15][12] = 111.16; AAmatrix[15][13] =  64.29; AAmatrix[15][14] = 169.90;
	AAmatrix[15][15] =   0.00; AAmatrix[15][16] = 597.21; AAmatrix[15][17] =  38.58; AAmatrix[15][18] =  64.92; AAmatrix[15][19] =   1.90;
	AAmatrix[16][ 0] = 480.72; AAmatrix[16][ 1] =   2.08; AAmatrix[16][ 2] = 238.46; AAmatrix[16][ 3] =  28.01; AAmatrix[16][ 4] = 179.97;
	AAmatrix[16][ 5] =  94.93; AAmatrix[16][ 6] =  14.82; AAmatrix[16][ 7] =  11.17; AAmatrix[16][ 8] =  44.78; AAmatrix[16][ 9] = 368.43;
	AAmatrix[16][10] = 126.40; AAmatrix[16][11] = 136.33; AAmatrix[16][12] = 528.17; AAmatrix[16][13] =  33.85; AAmatrix[16][14] = 128.22;
	AAmatrix[16][15] = 597.21; AAmatrix[16][16] =   0.00; AAmatrix[16][17] =   9.99; AAmatrix[16][18] =  38.73; AAmatrix[16][19] = 204.54;
	AAmatrix[17][ 0] =   1.90; AAmatrix[17][ 1] =  21.95; AAmatrix[17][ 2] =  10.68; AAmatrix[17][ 3] =  19.86; AAmatrix[17][ 4] =  33.60;
	AAmatrix[17][ 5] =   1.90; AAmatrix[17][ 6] =   1.90; AAmatrix[17][ 7] =  10.92; AAmatrix[17][ 8] =   7.08; AAmatrix[17][ 9] =   1.90;
	AAmatrix[17][10] =  32.44; AAmatrix[17][11] =  24.00; AAmatrix[17][12] =  21.71; AAmatrix[17][13] =   7.84; AAmatrix[17][14] =   4.21;
	AAmatrix[17][15] =  38.58; AAmatrix[17][16] =   9.99; AAmatrix[17][17] =   0.00; AAmatrix[17][18] =  26.25; AAmatrix[17][19] =   5.37;
	AAmatrix[18][ 0] =   6.48; AAmatrix[18][ 1] =   1.90; AAmatrix[18][ 2] = 191.36; AAmatrix[18][ 3] =  21.21; AAmatrix[18][ 4] = 254.77;
	AAmatrix[18][ 5] =  38.82; AAmatrix[18][ 6] =  13.12; AAmatrix[18][ 7] =   3.21; AAmatrix[18][ 8] = 670.14; AAmatrix[18][ 9] =  25.01;
	AAmatrix[18][10] =  44.15; AAmatrix[18][11] =  51.17; AAmatrix[18][12] =  39.96; AAmatrix[18][13] = 465.58; AAmatrix[18][14] =  16.21;
	AAmatrix[18][15] =  64.92; AAmatrix[18][16] =  38.73; AAmatrix[18][17] =  26.25; AAmatrix[18][18] =   0.00; AAmatrix[18][19] =   1.90;
	AAmatrix[19][ 0] = 195.06; AAmatrix[19][ 1] =   7.64; AAmatrix[19][ 2] =   1.90; AAmatrix[19][ 3] =   1.90; AAmatrix[19][ 4] =   1.90;
	AAmatrix[19][ 5] =  19.00; AAmatrix[19][ 6] =  21.14; AAmatrix[19][ 7] =   2.53; AAmatrix[19][ 8] =   1.90; AAmatrix[19][ 9] =1222.94;
	AAmatrix[19][10] =  91.67; AAmatrix[19][11] =   1.90; AAmatrix[19][12] = 387.54; AAmatrix[19][13] =   6.35; AAmatrix[19][14] =   8.23;
	AAmatrix[19][15] =   1.90; AAmatrix[19][16] = 204.54; AAmatrix[19][17] =   5.37; AAmatrix[19][18] =   1.90; AAmatrix[19][19] =   0.00;

	AAPi[ 0] = 0.072;
	AAPi[ 1] = 0.019;
	AAPi[ 2] = 0.039;
	AAPi[ 3] = 0.019;
	AAPi[ 4] = 0.006;
	AAPi[ 5] = 0.025;
	AAPi[ 6] = 0.024;
	AAPi[ 7] = 0.056;
	AAPi[ 8] = 0.028;
	AAPi[ 9] = 0.088;
	AAPi[10] = 0.168;
	AAPi[11] = 0.023;
	AAPi[12] = 0.054;
	AAPi[13] = 0.061;
	AAPi[14] = 0.054;
	AAPi[15] = 0.072;
	AAPi[16] = 0.086;
	AAPi[17] = 0.029;
	AAPi[18] = 0.033;
	AAPi[19] = 0.043;
	
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else
	if(modelnumber==9)
	{

	/* 
		rtRev
		
		Dimmic M.W., J.S. Rest, D.P. Mindell, and D. Goldstein. (2002) RArtREV: An amino acid substitution matrix for inference 
		of retrovirus and reverse transcriptase phylogeny. Journal of Molecular Evolution 55: 65-73.
	*/
	AAmatrix[ 0][ 0] =   0; AAmatrix[ 1][ 0] =  34; AAmatrix[ 2][ 0] =  51; AAmatrix[ 3][ 0] =  10; AAmatrix[ 4][ 0] = 439;
	AAmatrix[ 5][ 0] =  32; AAmatrix[ 6][ 0] =  81; AAmatrix[ 7][ 0] = 135; AAmatrix[ 8][ 0] =  30; AAmatrix[ 9][ 0] =   1;
	AAmatrix[10][ 0] =  45; AAmatrix[11][ 0] =  38; AAmatrix[12][ 0] = 235; AAmatrix[13][ 0] =   1; AAmatrix[14][ 0] =  97;
	AAmatrix[15][ 0] = 460; AAmatrix[16][ 0] = 258; AAmatrix[17][ 0] =   5; AAmatrix[18][ 0] =  55; AAmatrix[19][ 0] = 197;
	AAmatrix[ 0][ 1] =  34; AAmatrix[ 1][ 1] =   0; AAmatrix[ 2][ 1] =  35; AAmatrix[ 3][ 1] =  30; AAmatrix[ 4][ 1] =  92;
	AAmatrix[ 5][ 1] = 221; AAmatrix[ 6][ 1] =  10; AAmatrix[ 7][ 1] =  41; AAmatrix[ 8][ 1] =  90; AAmatrix[ 9][ 1] =  24;
	AAmatrix[10][ 1] =  18; AAmatrix[11][ 1] = 593; AAmatrix[12][ 1] =  57; AAmatrix[13][ 1] =   7; AAmatrix[14][ 1] =  24;
	AAmatrix[15][ 1] = 102; AAmatrix[16][ 1] =  64; AAmatrix[17][ 1] =  13; AAmatrix[18][ 1] =  47; AAmatrix[19][ 1] =  29;
	AAmatrix[ 0][ 2] =  51; AAmatrix[ 1][ 2] =  35; AAmatrix[ 2][ 2] =   0; AAmatrix[ 3][ 2] = 384; AAmatrix[ 4][ 2] = 128;
	AAmatrix[ 5][ 2] = 236; AAmatrix[ 6][ 2] =  79; AAmatrix[ 7][ 2] =  94; AAmatrix[ 8][ 2] = 320; AAmatrix[ 9][ 2] =  35;
	AAmatrix[10][ 2] =  15; AAmatrix[11][ 2] = 123; AAmatrix[12][ 2] =   1; AAmatrix[13][ 2] =  49; AAmatrix[14][ 2] =  33;
	AAmatrix[15][ 2] = 294; AAmatrix[16][ 2] = 148; AAmatrix[17][ 2] =  16; AAmatrix[18][ 2] =  28; AAmatrix[19][ 2] =  21;
	AAmatrix[ 0][ 3] =  10; AAmatrix[ 1][ 3] =  30; AAmatrix[ 2][ 3] = 384; AAmatrix[ 3][ 3] =   0; AAmatrix[ 4][ 3] =   1;
	AAmatrix[ 5][ 3] =  78; AAmatrix[ 6][ 3] = 542; AAmatrix[ 7][ 3] =  61; AAmatrix[ 8][ 3] =  91; AAmatrix[ 9][ 3] =   1;
	AAmatrix[10][ 3] =   5; AAmatrix[11][ 3] =  20; AAmatrix[12][ 3] =   1; AAmatrix[13][ 3] =   1; AAmatrix[14][ 3] =  55;
	AAmatrix[15][ 3] = 136; AAmatrix[16][ 3] =  55; AAmatrix[17][ 3] =   1; AAmatrix[18][ 3] =   1; AAmatrix[19][ 3] =   6;
	AAmatrix[ 0][ 4] = 439; AAmatrix[ 1][ 4] =  92; AAmatrix[ 2][ 4] = 128; AAmatrix[ 3][ 4] =   1; AAmatrix[ 4][ 4] =   0;
	AAmatrix[ 5][ 4] =  70; AAmatrix[ 6][ 4] =   1; AAmatrix[ 7][ 4] =  48; AAmatrix[ 8][ 4] = 124; AAmatrix[ 9][ 4] = 104;
	AAmatrix[10][ 4] = 110; AAmatrix[11][ 4] =  16; AAmatrix[12][ 4] = 156; AAmatrix[13][ 4] =  70; AAmatrix[14][ 4] =   1;
	AAmatrix[15][ 4] =  75; AAmatrix[16][ 4] = 117; AAmatrix[17][ 4] =  55; AAmatrix[18][ 4] = 131; AAmatrix[19][ 4] = 295;
	AAmatrix[ 0][ 5] =  32; AAmatrix[ 1][ 5] = 221; AAmatrix[ 2][ 5] = 236; AAmatrix[ 3][ 5] =  78; AAmatrix[ 4][ 5] =  70;
	AAmatrix[ 5][ 5] =   0; AAmatrix[ 6][ 5] = 372; AAmatrix[ 7][ 5] =  18; AAmatrix[ 8][ 5] = 387; AAmatrix[ 9][ 5] =  33;
	AAmatrix[10][ 5] =  54; AAmatrix[11][ 5] = 309; AAmatrix[12][ 5] = 158; AAmatrix[13][ 5] =   1; AAmatrix[14][ 5] =  68;
	AAmatrix[15][ 5] = 225; AAmatrix[16][ 5] = 146; AAmatrix[17][ 5] =  10; AAmatrix[18][ 5] =  45; AAmatrix[19][ 5] =  36;
	AAmatrix[ 0][ 6] =  81; AAmatrix[ 1][ 6] =  10; AAmatrix[ 2][ 6] =  79; AAmatrix[ 3][ 6] = 542; AAmatrix[ 4][ 6] =   1;
	AAmatrix[ 5][ 6] = 372; AAmatrix[ 6][ 6] =   0; AAmatrix[ 7][ 6] =  70; AAmatrix[ 8][ 6] =  34; AAmatrix[ 9][ 6] =   1;
	AAmatrix[10][ 6] =  21; AAmatrix[11][ 6] = 141; AAmatrix[12][ 6] =   1; AAmatrix[13][ 6] =   1; AAmatrix[14][ 6] =  52;
	AAmatrix[15][ 6] =  95; AAmatrix[16][ 6] =  82; AAmatrix[17][ 6] =  17; AAmatrix[18][ 6] =   1; AAmatrix[19][ 6] =  35;
	AAmatrix[ 0][ 7] = 135; AAmatrix[ 1][ 7] =  41; AAmatrix[ 2][ 7] =  94; AAmatrix[ 3][ 7] =  61; AAmatrix[ 4][ 7] =  48;
	AAmatrix[ 5][ 7] =  18; AAmatrix[ 6][ 7] =  70; AAmatrix[ 7][ 7] =   0; AAmatrix[ 8][ 7] =  68; AAmatrix[ 9][ 7] =   1;
	AAmatrix[10][ 7] =   3; AAmatrix[11][ 7] =  30; AAmatrix[12][ 7] =  37; AAmatrix[13][ 7] =   7; AAmatrix[14][ 7] =  17;
	AAmatrix[15][ 7] = 152; AAmatrix[16][ 7] =   7; AAmatrix[17][ 7] =  23; AAmatrix[18][ 7] =  21; AAmatrix[19][ 7] =   3;
	AAmatrix[ 0][ 8] =  30; AAmatrix[ 1][ 8] =  90; AAmatrix[ 2][ 8] = 320; AAmatrix[ 3][ 8] =  91; AAmatrix[ 4][ 8] = 124;
	AAmatrix[ 5][ 8] = 387; AAmatrix[ 6][ 8] =  34; AAmatrix[ 7][ 8] =  68; AAmatrix[ 8][ 8] =   0; AAmatrix[ 9][ 8] =  34;
	AAmatrix[10][ 8] =  51; AAmatrix[11][ 8] =  76; AAmatrix[12][ 8] = 116; AAmatrix[13][ 8] = 141; AAmatrix[14][ 8] =  44;
	AAmatrix[15][ 8] = 183; AAmatrix[16][ 8] =  49; AAmatrix[17][ 8] =  48; AAmatrix[18][ 8] = 307; AAmatrix[19][ 8] =   1;
	AAmatrix[ 0][ 9] =   1; AAmatrix[ 1][ 9] =  24; AAmatrix[ 2][ 9] =  35; AAmatrix[ 3][ 9] =   1; AAmatrix[ 4][ 9] = 104;
	AAmatrix[ 5][ 9] =  33; AAmatrix[ 6][ 9] =   1; AAmatrix[ 7][ 9] =   1; AAmatrix[ 8][ 9] =  34; AAmatrix[ 9][ 9] =   0;
	AAmatrix[10][ 9] = 385; AAmatrix[11][ 9] =  34; AAmatrix[12][ 9] = 375; AAmatrix[13][ 9] =  64; AAmatrix[14][ 9] =  10;
	AAmatrix[15][ 9] =   4; AAmatrix[16][ 9] =  72; AAmatrix[17][ 9] =  39; AAmatrix[18][ 9] =  26; AAmatrix[19][ 9] =1048;
	AAmatrix[ 0][10] =  45; AAmatrix[ 1][10] =  18; AAmatrix[ 2][10] =  15; AAmatrix[ 3][10] =   5; AAmatrix[ 4][10] = 110;
	AAmatrix[ 5][10] =  54; AAmatrix[ 6][10] =  21; AAmatrix[ 7][10] =   3; AAmatrix[ 8][10] =  51; AAmatrix[ 9][10] = 385;
	AAmatrix[10][10] =   0; AAmatrix[11][10] =  23; AAmatrix[12][10] = 581; AAmatrix[13][10] = 179; AAmatrix[14][10] =  22;
	AAmatrix[15][10] =  24; AAmatrix[16][10] =  25; AAmatrix[17][10] =  47; AAmatrix[18][10] =  64; AAmatrix[19][10] = 112;
	AAmatrix[ 0][11] =  38; AAmatrix[ 1][11] = 593; AAmatrix[ 2][11] = 123; AAmatrix[ 3][11] =  20; AAmatrix[ 4][11] =  16;
	AAmatrix[ 5][11] = 309; AAmatrix[ 6][11] = 141; AAmatrix[ 7][11] =  30; AAmatrix[ 8][11] =  76; AAmatrix[ 9][11] =  34;
	AAmatrix[10][11] =  23; AAmatrix[11][11] =   0; AAmatrix[12][11] = 134; AAmatrix[13][11] =  14; AAmatrix[14][11] =  43;
	AAmatrix[15][11] =  77; AAmatrix[16][11] = 110; AAmatrix[17][11] =   6; AAmatrix[18][11] =   1; AAmatrix[19][11] =  19;
	AAmatrix[ 0][12] = 235; AAmatrix[ 1][12] =  57; AAmatrix[ 2][12] =   1; AAmatrix[ 3][12] =   1; AAmatrix[ 4][12] = 156;
	AAmatrix[ 5][12] = 158; AAmatrix[ 6][12] =   1; AAmatrix[ 7][12] =  37; AAmatrix[ 8][12] = 116; AAmatrix[ 9][12] = 375;
	AAmatrix[10][12] = 581; AAmatrix[11][12] = 134; AAmatrix[12][12] =   0; AAmatrix[13][12] = 247; AAmatrix[14][12] =   1;
	AAmatrix[15][12] =   1; AAmatrix[16][12] = 131; AAmatrix[17][12] = 111; AAmatrix[18][12] =  74; AAmatrix[19][12] = 236;
	AAmatrix[ 0][13] =   1; AAmatrix[ 1][13] =   7; AAmatrix[ 2][13] =  49; AAmatrix[ 3][13] =   1; AAmatrix[ 4][13] =  70;
	AAmatrix[ 5][13] =   1; AAmatrix[ 6][13] =   1; AAmatrix[ 7][13] =   7; AAmatrix[ 8][13] = 141; AAmatrix[ 9][13] =  64;
	AAmatrix[10][13] = 179; AAmatrix[11][13] =  14; AAmatrix[12][13] = 247; AAmatrix[13][13] =   0; AAmatrix[14][13] =  11;
	AAmatrix[15][13] =  20; AAmatrix[16][13] =  69; AAmatrix[17][13] = 182; AAmatrix[18][13] =1017; AAmatrix[19][13] =  92;
	AAmatrix[ 0][14] =  97; AAmatrix[ 1][14] =  24; AAmatrix[ 2][14] =  33; AAmatrix[ 3][14] =  55; AAmatrix[ 4][14] =   1;
	AAmatrix[ 5][14] =  68; AAmatrix[ 6][14] =  52; AAmatrix[ 7][14] =  17; AAmatrix[ 8][14] =  44; AAmatrix[ 9][14] =  10;
	AAmatrix[10][14] =  22; AAmatrix[11][14] =  43; AAmatrix[12][14] =   1; AAmatrix[13][14] =  11; AAmatrix[14][14] =   0;
	AAmatrix[15][14] = 134; AAmatrix[16][14] =  62; AAmatrix[17][14] =   9; AAmatrix[18][14] =  14; AAmatrix[19][14] =  25;
	AAmatrix[ 0][15] = 460; AAmatrix[ 1][15] = 102; AAmatrix[ 2][15] = 294; AAmatrix[ 3][15] = 136; AAmatrix[ 4][15] =  75;
	AAmatrix[ 5][15] = 225; AAmatrix[ 6][15] =  95; AAmatrix[ 7][15] = 152; AAmatrix[ 8][15] = 183; AAmatrix[ 9][15] =   4;
	AAmatrix[10][15] =  24; AAmatrix[11][15] =  77; AAmatrix[12][15] =   1; AAmatrix[13][15] =  20; AAmatrix[14][15] = 134;
	AAmatrix[15][15] =   0; AAmatrix[16][15] = 671; AAmatrix[17][15] =  14; AAmatrix[18][15] =  31; AAmatrix[19][15] =  39;
	AAmatrix[ 0][16] = 258; AAmatrix[ 1][16] =  64; AAmatrix[ 2][16] = 148; AAmatrix[ 3][16] =  55; AAmatrix[ 4][16] = 117;
	AAmatrix[ 5][16] = 146; AAmatrix[ 6][16] =  82; AAmatrix[ 7][16] =   7; AAmatrix[ 8][16] =  49; AAmatrix[ 9][16] =  72;
	AAmatrix[10][16] =  25; AAmatrix[11][16] = 110; AAmatrix[12][16] = 131; AAmatrix[13][16] =  69; AAmatrix[14][16] =  62;
	AAmatrix[15][16] = 671; AAmatrix[16][16] =   0; AAmatrix[17][16] =   1; AAmatrix[18][16] =  34; AAmatrix[19][16] = 196;
	AAmatrix[ 0][17] =   5; AAmatrix[ 1][17] =  13; AAmatrix[ 2][17] =  16; AAmatrix[ 3][17] =   1; AAmatrix[ 4][17] =  55;
	AAmatrix[ 5][17] =  10; AAmatrix[ 6][17] =  17; AAmatrix[ 7][17] =  23; AAmatrix[ 8][17] =  48; AAmatrix[ 9][17] =  39;
	AAmatrix[10][17] =  47; AAmatrix[11][17] =   6; AAmatrix[12][17] = 111; AAmatrix[13][17] = 182; AAmatrix[14][17] =   9;
	AAmatrix[15][17] =  14; AAmatrix[16][17] =   1; AAmatrix[17][17] =   0; AAmatrix[18][17] = 176; AAmatrix[19][17] =  26;
	AAmatrix[ 0][18] =  55; AAmatrix[ 1][18] =  47; AAmatrix[ 2][18] =  28; AAmatrix[ 3][18] =   1; AAmatrix[ 4][18] = 131;
	AAmatrix[ 5][18] =  45; AAmatrix[ 6][18] =   1; AAmatrix[ 7][18] =  21; AAmatrix[ 8][18] = 307; AAmatrix[ 9][18] =  26;
	AAmatrix[10][18] =  64; AAmatrix[11][18] =   1; AAmatrix[12][18] =  74; AAmatrix[13][18] =1017; AAmatrix[14][18] =  14;
	AAmatrix[15][18] =  31; AAmatrix[16][18] =  34; AAmatrix[17][18] = 176; AAmatrix[18][18] =   0; AAmatrix[19][18] =  59;
	AAmatrix[ 0][19] = 197; AAmatrix[ 1][19] =  29; AAmatrix[ 2][19] =  21; AAmatrix[ 3][19] =   6; AAmatrix[ 4][19] = 295;
	AAmatrix[ 5][19] =  36; AAmatrix[ 6][19] =  35; AAmatrix[ 7][19] =   3; AAmatrix[ 8][19] =   1; AAmatrix[ 9][19] =1048;
	AAmatrix[10][19] = 112; AAmatrix[11][19] =  19; AAmatrix[12][19] = 236; AAmatrix[13][19] =  92; AAmatrix[14][19] =  25;
	AAmatrix[15][19] =  39; AAmatrix[16][19] = 196; AAmatrix[17][19] =  26; AAmatrix[18][19] =  59; AAmatrix[19][19] =   0;
	AAPi[ 0] = 0.0646;
	AAPi[ 1] = 0.0453;
	AAPi[ 2] = 0.0376;
	AAPi[ 3] = 0.0422;
	AAPi[ 4] = 0.0114;
	AAPi[ 5] = 0.0606;
	AAPi[ 6] = 0.0607;
	AAPi[ 7] = 0.0639;
	AAPi[ 8] = 0.0273;
	AAPi[ 9] = 0.0679;
	AAPi[10] = 0.1018;
	AAPi[11] = 0.0751;
	AAPi[12] = 0.0150;
	AAPi[13] = 0.0287;
	AAPi[14] = 0.0681;
	AAPi[15] = 0.0488;
	AAPi[16] = 0.0622;
	AAPi[17] = 0.0251;
	AAPi[18] = 0.0318;
	AAPi[19] = 0.0619;
	
	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	else
	if(modelnumber==10)
	{

	/* 
		cpRev 
		
		Adachi, J., P. J. Waddell, W. Martin, and M. Hasegawa. (2000) Plastid genome phylogeny and a model of amino acid 
		substitution for proteins encoded by chloroplast DNA. Journal of Molecular Evolution 50:348-358.
	*/
	AAmatrix[ 0][ 0] =    0; AAmatrix[ 0][ 1] =  105; AAmatrix[ 0][ 2] =  227; AAmatrix[ 0][ 3] =  175; AAmatrix[ 0][ 4] =  669; 
	AAmatrix[ 0][ 5] =  157; AAmatrix[ 0][ 6] =  499; AAmatrix[ 0][ 7] =  665; AAmatrix[ 0][ 8] =   66; AAmatrix[ 0][ 9] =  145; 
	AAmatrix[ 0][10] =  197; AAmatrix[ 0][11] =  236; AAmatrix[ 0][12] =  185; AAmatrix[ 0][13] =   68; AAmatrix[ 0][14] =  490; 
	AAmatrix[ 0][15] = 2440; AAmatrix[ 0][16] = 1340; AAmatrix[ 0][17] =   14; AAmatrix[ 0][18] =   56; AAmatrix[ 0][19] =  968; 
	AAmatrix[ 1][ 0] =  105; AAmatrix[ 1][ 1] =    0; AAmatrix[ 1][ 2] =  357; AAmatrix[ 1][ 3] =   43; AAmatrix[ 1][ 4] =  823; 
	AAmatrix[ 1][ 5] = 1745; AAmatrix[ 1][ 6] =  152; AAmatrix[ 1][ 7] =  243; AAmatrix[ 1][ 8] =  715; AAmatrix[ 1][ 9] =  136; 
	AAmatrix[ 1][10] =  203; AAmatrix[ 1][11] = 4482; AAmatrix[ 1][12] =  125; AAmatrix[ 1][13] =   53; AAmatrix[ 1][14] =   87; 
	AAmatrix[ 1][15] =  385; AAmatrix[ 1][16] =  314; AAmatrix[ 1][17] =  230; AAmatrix[ 1][18] =  323; AAmatrix[ 1][19] =   92; 
	AAmatrix[ 2][ 0] =  227; AAmatrix[ 2][ 1] =  357; AAmatrix[ 2][ 2] =    0; AAmatrix[ 2][ 3] = 4435; AAmatrix[ 2][ 4] =  538; 
	AAmatrix[ 2][ 5] =  768; AAmatrix[ 2][ 6] = 1055; AAmatrix[ 2][ 7] =  653; AAmatrix[ 2][ 8] = 1405; AAmatrix[ 2][ 9] =  168; 
	AAmatrix[ 2][10] =  113; AAmatrix[ 2][11] = 2430; AAmatrix[ 2][12] =   61; AAmatrix[ 2][13] =   97; AAmatrix[ 2][14] =  173; 
	AAmatrix[ 2][15] = 2085; AAmatrix[ 2][16] = 1393; AAmatrix[ 2][17] =   40; AAmatrix[ 2][18] =  754; AAmatrix[ 2][19] =   83; 
	AAmatrix[ 3][ 0] =  175; AAmatrix[ 3][ 1] =   43; AAmatrix[ 3][ 2] = 4435; AAmatrix[ 3][ 3] =    0; AAmatrix[ 3][ 4] =   10; 
	AAmatrix[ 3][ 5] =  400; AAmatrix[ 3][ 6] = 3691; AAmatrix[ 3][ 7] =  431; AAmatrix[ 3][ 8] =  331; AAmatrix[ 3][ 9] =   10; 
	AAmatrix[ 3][10] =   10; AAmatrix[ 3][11] =  412; AAmatrix[ 3][12] =   47; AAmatrix[ 3][13] =   22; AAmatrix[ 3][14] =  170; 
	AAmatrix[ 3][15] =  590; AAmatrix[ 3][16] =  266; AAmatrix[ 3][17] =   18; AAmatrix[ 3][18] =  281; AAmatrix[ 3][19] =   75; 
	AAmatrix[ 4][ 0] =  669; AAmatrix[ 4][ 1] =  823; AAmatrix[ 4][ 2] =  538; AAmatrix[ 4][ 3] =   10; AAmatrix[ 4][ 4] =    0; 
	AAmatrix[ 4][ 5] =   10; AAmatrix[ 4][ 6] =   10; AAmatrix[ 4][ 7] =  303; AAmatrix[ 4][ 8] =  441; AAmatrix[ 4][ 9] =  280; 
	AAmatrix[ 4][10] =  396; AAmatrix[ 4][11] =   48; AAmatrix[ 4][12] =  159; AAmatrix[ 4][13] =  726; AAmatrix[ 4][14] =  285; 
	AAmatrix[ 4][15] = 2331; AAmatrix[ 4][16] =  576; AAmatrix[ 4][17] =  435; AAmatrix[ 4][18] = 1466; AAmatrix[ 4][19] =  592; 
	AAmatrix[ 5][ 0] =  157; AAmatrix[ 5][ 1] = 1745; AAmatrix[ 5][ 2] =  768; AAmatrix[ 5][ 3] =  400; AAmatrix[ 5][ 4] =   10; 
	AAmatrix[ 5][ 5] =    0; AAmatrix[ 5][ 6] = 3122; AAmatrix[ 5][ 7] =  133; AAmatrix[ 5][ 8] = 1269; AAmatrix[ 5][ 9] =   92; 
	AAmatrix[ 5][10] =  286; AAmatrix[ 5][11] = 3313; AAmatrix[ 5][12] =  202; AAmatrix[ 5][13] =   10; AAmatrix[ 5][14] =  323; 
	AAmatrix[ 5][15] =  396; AAmatrix[ 5][16] =  241; AAmatrix[ 5][17] =   53; AAmatrix[ 5][18] =  391; AAmatrix[ 5][19] =   54; 
	AAmatrix[ 6][ 0] =  499; AAmatrix[ 6][ 1] =  152; AAmatrix[ 6][ 2] = 1055; AAmatrix[ 6][ 3] = 3691; AAmatrix[ 6][ 4] =   10; 
	AAmatrix[ 6][ 5] = 3122; AAmatrix[ 6][ 6] =    0; AAmatrix[ 6][ 7] =  379; AAmatrix[ 6][ 8] =  162; AAmatrix[ 6][ 9] =  148; 
	AAmatrix[ 6][10] =   82; AAmatrix[ 6][11] = 2629; AAmatrix[ 6][12] =  113; AAmatrix[ 6][13] =  145; AAmatrix[ 6][14] =  185; 
	AAmatrix[ 6][15] =  568; AAmatrix[ 6][16] =  369; AAmatrix[ 6][17] =   63; AAmatrix[ 6][18] =  142; AAmatrix[ 6][19] =  200; 
	AAmatrix[ 7][ 0] =  665; AAmatrix[ 7][ 1] =  243; AAmatrix[ 7][ 2] =  653; AAmatrix[ 7][ 3] =  431; AAmatrix[ 7][ 4] =  303; 
	AAmatrix[ 7][ 5] =  133; AAmatrix[ 7][ 6] =  379; AAmatrix[ 7][ 7] =    0; AAmatrix[ 7][ 8] =   19; AAmatrix[ 7][ 9] =   40; 
	AAmatrix[ 7][10] =   20; AAmatrix[ 7][11] =  263; AAmatrix[ 7][12] =   21; AAmatrix[ 7][13] =   25; AAmatrix[ 7][14] =   28; 
	AAmatrix[ 7][15] =  691; AAmatrix[ 7][16] =   92; AAmatrix[ 7][17] =   82; AAmatrix[ 7][18] =   10; AAmatrix[ 7][19] =   91; 
	AAmatrix[ 8][ 0] =   66; AAmatrix[ 8][ 1] =  715; AAmatrix[ 8][ 2] = 1405; AAmatrix[ 8][ 3] =  331; AAmatrix[ 8][ 4] =  441; 
	AAmatrix[ 8][ 5] = 1269; AAmatrix[ 8][ 6] =  162; AAmatrix[ 8][ 7] =   19; AAmatrix[ 8][ 8] =    0; AAmatrix[ 8][ 9] =   29; 
	AAmatrix[ 8][10] =   66; AAmatrix[ 8][11] =  305; AAmatrix[ 8][12] =   10; AAmatrix[ 8][13] =  127; AAmatrix[ 8][14] =  152; 
	AAmatrix[ 8][15] =  303; AAmatrix[ 8][16] =   32; AAmatrix[ 8][17] =   69; AAmatrix[ 8][18] = 1971; AAmatrix[ 8][19] =   25; 
	AAmatrix[ 9][ 0] =  145; AAmatrix[ 9][ 1] =  136; AAmatrix[ 9][ 2] =  168; AAmatrix[ 9][ 3] =   10; AAmatrix[ 9][ 4] =  280; 
	AAmatrix[ 9][ 5] =   92; AAmatrix[ 9][ 6] =  148; AAmatrix[ 9][ 7] =   40; AAmatrix[ 9][ 8] =   29; AAmatrix[ 9][ 9] =    0; 
	AAmatrix[ 9][10] = 1745; AAmatrix[ 9][11] =  345; AAmatrix[ 9][12] = 1772; AAmatrix[ 9][13] =  454; AAmatrix[ 9][14] =  117; 
	AAmatrix[ 9][15] =  216; AAmatrix[ 9][16] = 1040; AAmatrix[ 9][17] =   42; AAmatrix[ 9][18] =   89; AAmatrix[ 9][19] = 4797; 
	AAmatrix[10][ 0] =  197; AAmatrix[10][ 1] =  203; AAmatrix[10][ 2] =  113; AAmatrix[10][ 3] =   10; AAmatrix[10][ 4] =  396; 
	AAmatrix[10][ 5] =  286; AAmatrix[10][ 6] =   82; AAmatrix[10][ 7] =   20; AAmatrix[10][ 8] =   66; AAmatrix[10][ 9] = 1745; 
	AAmatrix[10][10] =    0; AAmatrix[10][11] =  218; AAmatrix[10][12] = 1351; AAmatrix[10][13] = 1268; AAmatrix[10][14] =  219; 
	AAmatrix[10][15] =  516; AAmatrix[10][16] =  156; AAmatrix[10][17] =  159; AAmatrix[10][18] =  189; AAmatrix[10][19] =  865; 
	AAmatrix[11][ 0] =  236; AAmatrix[11][ 1] = 4482; AAmatrix[11][ 2] = 2430; AAmatrix[11][ 3] =  412; AAmatrix[11][ 4] =   48; 
	AAmatrix[11][ 5] = 3313; AAmatrix[11][ 6] = 2629; AAmatrix[11][ 7] =  263; AAmatrix[11][ 8] =  305; AAmatrix[11][ 9] =  345; 
	AAmatrix[11][10] =  218; AAmatrix[11][11] =    0; AAmatrix[11][12] =  193; AAmatrix[11][13] =   72; AAmatrix[11][14] =  302; 
	AAmatrix[11][15] =  868; AAmatrix[11][16] =  918; AAmatrix[11][17] =   10; AAmatrix[11][18] =  247; AAmatrix[11][19] =  249; 
	AAmatrix[12][ 0] =  185; AAmatrix[12][ 1] =  125; AAmatrix[12][ 2] =   61; AAmatrix[12][ 3] =   47; AAmatrix[12][ 4] =  159; 
	AAmatrix[12][ 5] =  202; AAmatrix[12][ 6] =  113; AAmatrix[12][ 7] =   21; AAmatrix[12][ 8] =   10; AAmatrix[12][ 9] = 1772; 
	AAmatrix[12][10] = 1351; AAmatrix[12][11] =  193; AAmatrix[12][12] =    0; AAmatrix[12][13] =  327; AAmatrix[12][14] =  100; 
	AAmatrix[12][15] =   93; AAmatrix[12][16] =  645; AAmatrix[12][17] =   86; AAmatrix[12][18] =  215; AAmatrix[12][19] =  475; 
	AAmatrix[13][ 0] =   68; AAmatrix[13][ 1] =   53; AAmatrix[13][ 2] =   97; AAmatrix[13][ 3] =   22; AAmatrix[13][ 4] =  726; 
	AAmatrix[13][ 5] =   10; AAmatrix[13][ 6] =  145; AAmatrix[13][ 7] =   25; AAmatrix[13][ 8] =  127; AAmatrix[13][ 9] =  454; 
	AAmatrix[13][10] = 1268; AAmatrix[13][11] =   72; AAmatrix[13][12] =  327; AAmatrix[13][13] =    0; AAmatrix[13][14] =   43; 
	AAmatrix[13][15] =  487; AAmatrix[13][16] =  148; AAmatrix[13][17] =  468; AAmatrix[13][18] = 2370; AAmatrix[13][19] =  317; 
	AAmatrix[14][ 0] =  490; AAmatrix[14][ 1] =   87; AAmatrix[14][ 2] =  173; AAmatrix[14][ 3] =  170; AAmatrix[14][ 4] =  285; 
	AAmatrix[14][ 5] =  323; AAmatrix[14][ 6] =  185; AAmatrix[14][ 7] =   28; AAmatrix[14][ 8] =  152; AAmatrix[14][ 9] =  117; 
	AAmatrix[14][10] =  219; AAmatrix[14][11] =  302; AAmatrix[14][12] =  100; AAmatrix[14][13] =   43; AAmatrix[14][14] =    0; 
	AAmatrix[14][15] = 1202; AAmatrix[14][16] =  260; AAmatrix[14][17] =   49; AAmatrix[14][18] =   97; AAmatrix[14][19] =  122; 
	AAmatrix[15][ 0] = 2440; AAmatrix[15][ 1] =  385; AAmatrix[15][ 2] = 2085; AAmatrix[15][ 3] =  590; AAmatrix[15][ 4] = 2331; 
	AAmatrix[15][ 5] =  396; AAmatrix[15][ 6] =  568; AAmatrix[15][ 7] =  691; AAmatrix[15][ 8] =  303; AAmatrix[15][ 9] =  216; 
	AAmatrix[15][10] =  516; AAmatrix[15][11] =  868; AAmatrix[15][12] =   93; AAmatrix[15][13] =  487; AAmatrix[15][14] = 1202; 
	AAmatrix[15][15] =    0; AAmatrix[15][16] = 2151; AAmatrix[15][17] =   73; AAmatrix[15][18] =  522; AAmatrix[15][19] =  167; 
	AAmatrix[16][ 0] = 1340; AAmatrix[16][ 1] =  314; AAmatrix[16][ 2] = 1393; AAmatrix[16][ 3] =  266; AAmatrix[16][ 4] =  576; 
	AAmatrix[16][ 5] =  241; AAmatrix[16][ 6] =  369; AAmatrix[16][ 7] =   92; AAmatrix[16][ 8] =   32; AAmatrix[16][ 9] = 1040; 
	AAmatrix[16][10] =  156; AAmatrix[16][11] =  918; AAmatrix[16][12] =  645; AAmatrix[16][13] =  148; AAmatrix[16][14] =  260; 
	AAmatrix[16][15] = 2151; AAmatrix[16][16] =    0; AAmatrix[16][17] =   29; AAmatrix[16][18] =   71; AAmatrix[16][19] =  760; 
	AAmatrix[17][ 0] =   14; AAmatrix[17][ 1] =  230; AAmatrix[17][ 2] =   40; AAmatrix[17][ 3] =   18; AAmatrix[17][ 4] =  435; 
	AAmatrix[17][ 5] =   53; AAmatrix[17][ 6] =   63; AAmatrix[17][ 7] =   82; AAmatrix[17][ 8] =   69; AAmatrix[17][ 9] =   42; 
	AAmatrix[17][10] =  159; AAmatrix[17][11] =   10; AAmatrix[17][12] =   86; AAmatrix[17][13] =  468; AAmatrix[17][14] =   49; 
	AAmatrix[17][15] =   73; AAmatrix[17][16] =   29; AAmatrix[17][17] =    0; AAmatrix[17][18] =  346; AAmatrix[17][19] =   10; 
	AAmatrix[18][ 0] =   56; AAmatrix[18][ 1] =  323; AAmatrix[18][ 2] =  754; AAmatrix[18][ 3] =  281; AAmatrix[18][ 4] = 1466; 
	AAmatrix[18][ 5] =  391; AAmatrix[18][ 6] =  142; AAmatrix[18][ 7] =   10; AAmatrix[18][ 8] = 1971; AAmatrix[18][ 9] =   89; 
	AAmatrix[18][10] =  189; AAmatrix[18][11] =  247; AAmatrix[18][12] =  215; AAmatrix[18][13] = 2370; AAmatrix[18][14] =   97; 
	AAmatrix[18][15] =  522; AAmatrix[18][16] =   71; AAmatrix[18][17] =  346; AAmatrix[18][18] =    0; AAmatrix[18][19] =  119; 
	AAmatrix[19][ 0] =  968; AAmatrix[19][ 1] =   92; AAmatrix[19][ 2] =   83; AAmatrix[19][ 3] =   75; AAmatrix[19][ 4] =  592; 
	AAmatrix[19][ 5] =   54; AAmatrix[19][ 6] =  200; AAmatrix[19][ 7] =   91; AAmatrix[19][ 8] =   25; AAmatrix[19][ 9] = 4797; 
	AAmatrix[19][10] =  865; AAmatrix[19][11] =  249; AAmatrix[19][12] =  475; AAmatrix[19][13] =  317; AAmatrix[19][14] =  122; 
	AAmatrix[19][15] =  167; AAmatrix[19][16] =  760; AAmatrix[19][17] =   10; AAmatrix[19][18] =  119; AAmatrix[19][19] =    0; 

	AAPi[0] = 0.076;
	AAPi[1] = 0.062;
	AAPi[2] = 0.041;
	AAPi[3] = 0.037;
	AAPi[4] = 0.009;
	AAPi[5] = 0.038;
	AAPi[6] = 0.049;
	AAPi[7] = 0.084;
	AAPi[8] = 0.025;
	AAPi[9] = 0.081;
	AAPi[10] = 0.101;
	AAPi[11] = 0.050;
	AAPi[12] = 0.022;
	AAPi[13] = 0.051;
	AAPi[14] = 0.043;
	AAPi[15] = 0.062;
	AAPi[16] = 0.054;
	AAPi[17] = 0.018;
	AAPi[18] = 0.031;
	AAPi[19] = 0.066;

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==11)
	{

	/* 
		VT model
		
		Müller T, Vingron M (2000). Modeling amino acid replacement. J Computat Biol 7(6), 761-776
	*/
	AAmatrix[ 0][ 0] = 0.000000; AAmatrix[ 0][ 1] = 0.233108; AAmatrix[ 0][ 2] = 0.199097; AAmatrix[ 0][ 3] = 0.265145; AAmatrix[ 0][ 4] = 0.227333; 
	AAmatrix[ 0][ 5] = 0.310084; AAmatrix[ 0][ 6] = 0.567957; AAmatrix[ 0][ 7] = 0.876213; AAmatrix[ 0][ 8] = 0.078692; AAmatrix[ 0][ 9] = 0.222972; 
	AAmatrix[ 0][10] = 0.424630; AAmatrix[ 0][11] = 0.393245; AAmatrix[ 0][12] = 0.211550; AAmatrix[ 0][13] = 0.116646; AAmatrix[ 0][14] = 0.399143; 
	AAmatrix[ 0][15] = 1.817198; AAmatrix[ 0][16] = 0.877877; AAmatrix[ 0][17] = 0.030309; AAmatrix[ 0][18] = 0.087061; AAmatrix[ 0][19] = 1.230985; 
	AAmatrix[ 1][ 0] = 0.233108; AAmatrix[ 1][ 1] = 0.000000; AAmatrix[ 1][ 2] = 0.210797; AAmatrix[ 1][ 3] = 0.105191; AAmatrix[ 1][ 4] = 0.031726; 
	AAmatrix[ 1][ 5] = 0.493763; AAmatrix[ 1][ 6] = 0.255240; AAmatrix[ 1][ 7] = 0.156945; AAmatrix[ 1][ 8] = 0.213164; AAmatrix[ 1][ 9] = 0.081510; 
	AAmatrix[ 1][10] = 0.192364; AAmatrix[ 1][11] = 1.755838; AAmatrix[ 1][12] = 0.087930; AAmatrix[ 1][13] = 0.042569; AAmatrix[ 1][14] = 0.128480; 
	AAmatrix[ 1][15] = 0.292327; AAmatrix[ 1][16] = 0.204109; AAmatrix[ 1][17] = 0.046417; AAmatrix[ 1][18] = 0.097010; AAmatrix[ 1][19] = 0.113146; 
	AAmatrix[ 2][ 0] = 0.199097; AAmatrix[ 2][ 1] = 0.210797; AAmatrix[ 2][ 2] = 0.000000; AAmatrix[ 2][ 3] = 0.883422; AAmatrix[ 2][ 4] = 0.027495; 
	AAmatrix[ 2][ 5] = 0.275700; AAmatrix[ 2][ 6] = 0.270417; AAmatrix[ 2][ 7] = 0.362028; AAmatrix[ 2][ 8] = 0.290006; AAmatrix[ 2][ 9] = 0.087225; 
	AAmatrix[ 2][10] = 0.069245; AAmatrix[ 2][11] = 0.503060; AAmatrix[ 2][12] = 0.057420; AAmatrix[ 2][13] = 0.039769; AAmatrix[ 2][14] = 0.083956; 
	AAmatrix[ 2][15] = 0.847049; AAmatrix[ 2][16] = 0.471268; AAmatrix[ 2][17] = 0.010459; AAmatrix[ 2][18] = 0.093268; AAmatrix[ 2][19] = 0.049824; 
	AAmatrix[ 3][ 0] = 0.265145; AAmatrix[ 3][ 1] = 0.105191; AAmatrix[ 3][ 2] = 0.883422; AAmatrix[ 3][ 3] = 0.000000; AAmatrix[ 3][ 4] = 0.010313; 
	AAmatrix[ 3][ 5] = 0.205842; AAmatrix[ 3][ 6] = 1.599461; AAmatrix[ 3][ 7] = 0.311718; AAmatrix[ 3][ 8] = 0.134252; AAmatrix[ 3][ 9] = 0.011720; 
	AAmatrix[ 3][10] = 0.060863; AAmatrix[ 3][11] = 0.261101; AAmatrix[ 3][12] = 0.012182; AAmatrix[ 3][13] = 0.016577; AAmatrix[ 3][14] = 0.160063; 
	AAmatrix[ 3][15] = 0.461519; AAmatrix[ 3][16] = 0.178197; AAmatrix[ 3][17] = 0.011393; AAmatrix[ 3][18] = 0.051664; AAmatrix[ 3][19] = 0.048769; 
	AAmatrix[ 4][ 0] = 0.227333; AAmatrix[ 4][ 1] = 0.031726; AAmatrix[ 4][ 2] = 0.027495; AAmatrix[ 4][ 3] = 0.010313; AAmatrix[ 4][ 4] = 0.000000; 
	AAmatrix[ 4][ 5] = 0.004315; AAmatrix[ 4][ 6] = 0.005321; AAmatrix[ 4][ 7] = 0.050876; AAmatrix[ 4][ 8] = 0.016695; AAmatrix[ 4][ 9] = 0.046398; 
	AAmatrix[ 4][10] = 0.091709; AAmatrix[ 4][11] = 0.004067; AAmatrix[ 4][12] = 0.023690; AAmatrix[ 4][13] = 0.051127; AAmatrix[ 4][14] = 0.011137; 
	AAmatrix[ 4][15] = 0.175270; AAmatrix[ 4][16] = 0.079511; AAmatrix[ 4][17] = 0.007732; AAmatrix[ 4][18] = 0.042823; AAmatrix[ 4][19] = 0.163831; 
	AAmatrix[ 5][ 0] = 0.310084; AAmatrix[ 5][ 1] = 0.493763; AAmatrix[ 5][ 2] = 0.275700; AAmatrix[ 5][ 3] = 0.205842; AAmatrix[ 5][ 4] = 0.004315; 
	AAmatrix[ 5][ 5] = 0.000000; AAmatrix[ 5][ 6] = 0.960976; AAmatrix[ 5][ 7] = 0.128660; AAmatrix[ 5][ 8] = 0.315521; AAmatrix[ 5][ 9] = 0.054602; 
	AAmatrix[ 5][10] = 0.243530; AAmatrix[ 5][11] = 0.738208; AAmatrix[ 5][12] = 0.120801; AAmatrix[ 5][13] = 0.026235; AAmatrix[ 5][14] = 0.156570; 
	AAmatrix[ 5][15] = 0.358017; AAmatrix[ 5][16] = 0.248992; AAmatrix[ 5][17] = 0.021248; AAmatrix[ 5][18] = 0.062544; AAmatrix[ 5][19] = 0.112027; 
	AAmatrix[ 6][ 0] = 0.567957; AAmatrix[ 6][ 1] = 0.255240; AAmatrix[ 6][ 2] = 0.270417; AAmatrix[ 6][ 3] = 1.599461; AAmatrix[ 6][ 4] = 0.005321; 
	AAmatrix[ 6][ 5] = 0.960976; AAmatrix[ 6][ 6] = 0.000000; AAmatrix[ 6][ 7] = 0.250447; AAmatrix[ 6][ 8] = 0.104458; AAmatrix[ 6][ 9] = 0.046589; 
	AAmatrix[ 6][10] = 0.151924; AAmatrix[ 6][11] = 0.888630; AAmatrix[ 6][12] = 0.058643; AAmatrix[ 6][13] = 0.028168; AAmatrix[ 6][14] = 0.205134; 
	AAmatrix[ 6][15] = 0.406035; AAmatrix[ 6][16] = 0.321028; AAmatrix[ 6][17] = 0.018844; AAmatrix[ 6][18] = 0.055200; AAmatrix[ 6][19] = 0.205868; 
	AAmatrix[ 7][ 0] = 0.876213; AAmatrix[ 7][ 1] = 0.156945; AAmatrix[ 7][ 2] = 0.362028; AAmatrix[ 7][ 3] = 0.311718; AAmatrix[ 7][ 4] = 0.050876; 
	AAmatrix[ 7][ 5] = 0.128660; AAmatrix[ 7][ 6] = 0.250447; AAmatrix[ 7][ 7] = 0.000000; AAmatrix[ 7][ 8] = 0.058131; AAmatrix[ 7][ 9] = 0.051089; 
	AAmatrix[ 7][10] = 0.087056; AAmatrix[ 7][11] = 0.193243; AAmatrix[ 7][12] = 0.046560; AAmatrix[ 7][13] = 0.050143; AAmatrix[ 7][14] = 0.124492; 
	AAmatrix[ 7][15] = 0.612843; AAmatrix[ 7][16] = 0.136266; AAmatrix[ 7][17] = 0.023990; AAmatrix[ 7][18] = 0.037568; AAmatrix[ 7][19] = 0.082579; 
	AAmatrix[ 8][ 0] = 0.078692; AAmatrix[ 8][ 1] = 0.213164; AAmatrix[ 8][ 2] = 0.290006; AAmatrix[ 8][ 3] = 0.134252; AAmatrix[ 8][ 4] = 0.016695; 
	AAmatrix[ 8][ 5] = 0.315521; AAmatrix[ 8][ 6] = 0.104458; AAmatrix[ 8][ 7] = 0.058131; AAmatrix[ 8][ 8] = 0.000000; AAmatrix[ 8][ 9] = 0.020039; 
	AAmatrix[ 8][10] = 0.103552; AAmatrix[ 8][11] = 0.153323; AAmatrix[ 8][12] = 0.021157; AAmatrix[ 8][13] = 0.079807; AAmatrix[ 8][14] = 0.078892; 
	AAmatrix[ 8][15] = 0.167406; AAmatrix[ 8][16] = 0.101117; AAmatrix[ 8][17] = 0.020009; AAmatrix[ 8][18] = 0.286027; AAmatrix[ 8][19] = 0.068575; 
	AAmatrix[ 9][ 0] = 0.222972; AAmatrix[ 9][ 1] = 0.081510; AAmatrix[ 9][ 2] = 0.087225; AAmatrix[ 9][ 3] = 0.011720; AAmatrix[ 9][ 4] = 0.046398; 
	AAmatrix[ 9][ 5] = 0.054602; AAmatrix[ 9][ 6] = 0.046589; AAmatrix[ 9][ 7] = 0.051089; AAmatrix[ 9][ 8] = 0.020039; AAmatrix[ 9][ 9] = 0.000000; 
	AAmatrix[ 9][10] = 2.089890; AAmatrix[ 9][11] = 0.093181; AAmatrix[ 9][12] = 0.493845; AAmatrix[ 9][13] = 0.321020; AAmatrix[ 9][14] = 0.054797; 
	AAmatrix[ 9][15] = 0.081567; AAmatrix[ 9][16] = 0.376588; AAmatrix[ 9][17] = 0.034954; AAmatrix[ 9][18] = 0.086237; AAmatrix[ 9][19] = 3.654430; 
	AAmatrix[10][ 0] = 0.424630; AAmatrix[10][ 1] = 0.192364; AAmatrix[10][ 2] = 0.069245; AAmatrix[10][ 3] = 0.060863; AAmatrix[10][ 4] = 0.091709; 
	AAmatrix[10][ 5] = 0.243530; AAmatrix[10][ 6] = 0.151924; AAmatrix[10][ 7] = 0.087056; AAmatrix[10][ 8] = 0.103552; AAmatrix[10][ 9] = 2.089890; 
	AAmatrix[10][10] = 0.000000; AAmatrix[10][11] = 0.201204; AAmatrix[10][12] = 1.105667; AAmatrix[10][13] = 0.946499; AAmatrix[10][14] = 0.169784; 
	AAmatrix[10][15] = 0.214977; AAmatrix[10][16] = 0.243227; AAmatrix[10][17] = 0.083439; AAmatrix[10][18] = 0.189842; AAmatrix[10][19] = 1.337571; 
	AAmatrix[11][ 0] = 0.393245; AAmatrix[11][ 1] = 1.755838; AAmatrix[11][ 2] = 0.503060; AAmatrix[11][ 3] = 0.261101; AAmatrix[11][ 4] = 0.004067; 
	AAmatrix[11][ 5] = 0.738208; AAmatrix[11][ 6] = 0.888630; AAmatrix[11][ 7] = 0.193243; AAmatrix[11][ 8] = 0.153323; AAmatrix[11][ 9] = 0.093181; 
	AAmatrix[11][10] = 0.201204; AAmatrix[11][11] = 0.000000; AAmatrix[11][12] = 0.096474; AAmatrix[11][13] = 0.038261; AAmatrix[11][14] = 0.212302; 
	AAmatrix[11][15] = 0.400072; AAmatrix[11][16] = 0.446646; AAmatrix[11][17] = 0.023321; AAmatrix[11][18] = 0.068689; AAmatrix[11][19] = 0.144587; 
	AAmatrix[12][ 0] = 0.211550; AAmatrix[12][ 1] = 0.087930; AAmatrix[12][ 2] = 0.057420; AAmatrix[12][ 3] = 0.012182; AAmatrix[12][ 4] = 0.023690; 
	AAmatrix[12][ 5] = 0.120801; AAmatrix[12][ 6] = 0.058643; AAmatrix[12][ 7] = 0.046560; AAmatrix[12][ 8] = 0.021157; AAmatrix[12][ 9] = 0.493845; 
	AAmatrix[12][10] = 1.105667; AAmatrix[12][11] = 0.096474; AAmatrix[12][12] = 0.000000; AAmatrix[12][13] = 0.173052; AAmatrix[12][14] = 0.010363; 
	AAmatrix[12][15] = 0.090515; AAmatrix[12][16] = 0.184609; AAmatrix[12][17] = 0.022019; AAmatrix[12][18] = 0.073223; AAmatrix[12][19] = 0.307309; 
	AAmatrix[13][ 0] = 0.116646; AAmatrix[13][ 1] = 0.042569; AAmatrix[13][ 2] = 0.039769; AAmatrix[13][ 3] = 0.016577; AAmatrix[13][ 4] = 0.051127; 
	AAmatrix[13][ 5] = 0.026235; AAmatrix[13][ 6] = 0.028168; AAmatrix[13][ 7] = 0.050143; AAmatrix[13][ 8] = 0.079807; AAmatrix[13][ 9] = 0.321020; 
	AAmatrix[13][10] = 0.946499; AAmatrix[13][11] = 0.038261; AAmatrix[13][12] = 0.173052; AAmatrix[13][13] = 0.000000; AAmatrix[13][14] = 0.042564; 
	AAmatrix[13][15] = 0.138119; AAmatrix[13][16] = 0.085870; AAmatrix[13][17] = 0.128050; AAmatrix[13][18] = 0.898663; AAmatrix[13][19] = 0.247329; 
	AAmatrix[14][ 0] = 0.399143; AAmatrix[14][ 1] = 0.128480; AAmatrix[14][ 2] = 0.083956; AAmatrix[14][ 3] = 0.160063; AAmatrix[14][ 4] = 0.011137; 
	AAmatrix[14][ 5] = 0.156570; AAmatrix[14][ 6] = 0.205134; AAmatrix[14][ 7] = 0.124492; AAmatrix[14][ 8] = 0.078892; AAmatrix[14][ 9] = 0.054797; 
	AAmatrix[14][10] = 0.169784; AAmatrix[14][11] = 0.212302; AAmatrix[14][12] = 0.010363; AAmatrix[14][13] = 0.042564; AAmatrix[14][14] = 0.000000; 
	AAmatrix[14][15] = 0.430431; AAmatrix[14][16] = 0.207143; AAmatrix[14][17] = 0.014584; AAmatrix[14][18] = 0.032043; AAmatrix[14][19] = 0.129315; 
	AAmatrix[15][ 0] = 1.817198; AAmatrix[15][ 1] = 0.292327; AAmatrix[15][ 2] = 0.847049; AAmatrix[15][ 3] = 0.461519; AAmatrix[15][ 4] = 0.175270; 
	AAmatrix[15][ 5] = 0.358017; AAmatrix[15][ 6] = 0.406035; AAmatrix[15][ 7] = 0.612843; AAmatrix[15][ 8] = 0.167406; AAmatrix[15][ 9] = 0.081567; 
	AAmatrix[15][10] = 0.214977; AAmatrix[15][11] = 0.400072; AAmatrix[15][12] = 0.090515; AAmatrix[15][13] = 0.138119; AAmatrix[15][14] = 0.430431; 
	AAmatrix[15][15] = 0.000000; AAmatrix[15][16] = 1.767766; AAmatrix[15][17] = 0.035933; AAmatrix[15][18] = 0.121979; AAmatrix[15][19] = 0.127700; 
	AAmatrix[16][ 0] = 0.877877; AAmatrix[16][ 1] = 0.204109; AAmatrix[16][ 2] = 0.471268; AAmatrix[16][ 3] = 0.178197; AAmatrix[16][ 4] = 0.079511; 
	AAmatrix[16][ 5] = 0.248992; AAmatrix[16][ 6] = 0.321028; AAmatrix[16][ 7] = 0.136266; AAmatrix[16][ 8] = 0.101117; AAmatrix[16][ 9] = 0.376588; 
	AAmatrix[16][10] = 0.243227; AAmatrix[16][11] = 0.446646; AAmatrix[16][12] = 0.184609; AAmatrix[16][13] = 0.085870; AAmatrix[16][14] = 0.207143; 
	AAmatrix[16][15] = 1.767766; AAmatrix[16][16] = 0.000000; AAmatrix[16][17] = 0.020437; AAmatrix[16][18] = 0.094617; AAmatrix[16][19] = 0.740372; 
	AAmatrix[17][ 0] = 0.030309; AAmatrix[17][ 1] = 0.046417; AAmatrix[17][ 2] = 0.010459; AAmatrix[17][ 3] = 0.011393; AAmatrix[17][ 4] = 0.007732; 
	AAmatrix[17][ 5] = 0.021248; AAmatrix[17][ 6] = 0.018844; AAmatrix[17][ 7] = 0.023990; AAmatrix[17][ 8] = 0.020009; AAmatrix[17][ 9] = 0.034954; 
	AAmatrix[17][10] = 0.083439; AAmatrix[17][11] = 0.023321; AAmatrix[17][12] = 0.022019; AAmatrix[17][13] = 0.128050; AAmatrix[17][14] = 0.014584; 
	AAmatrix[17][15] = 0.035933; AAmatrix[17][16] = 0.020437; AAmatrix[17][17] = 0.000000; AAmatrix[17][18] = 0.124746; AAmatrix[17][19] = 0.022134; 
	AAmatrix[18][ 0] = 0.087061; AAmatrix[18][ 1] = 0.097010; AAmatrix[18][ 2] = 0.093268; AAmatrix[18][ 3] = 0.051664; AAmatrix[18][ 4] = 0.042823; 
	AAmatrix[18][ 5] = 0.062544; AAmatrix[18][ 6] = 0.055200; AAmatrix[18][ 7] = 0.037568; AAmatrix[18][ 8] = 0.286027; AAmatrix[18][ 9] = 0.086237; 
	AAmatrix[18][10] = 0.189842; AAmatrix[18][11] = 0.068689; AAmatrix[18][12] = 0.073223; AAmatrix[18][13] = 0.898663; AAmatrix[18][14] = 0.032043; 
	AAmatrix[18][15] = 0.121979; AAmatrix[18][16] = 0.094617; AAmatrix[18][17] = 0.124746; AAmatrix[18][18] = 0.000000; AAmatrix[18][19] = 0.125733; 
	AAmatrix[19][ 0] = 1.230985; AAmatrix[19][ 1] = 0.113146; AAmatrix[19][ 2] = 0.049824; AAmatrix[19][ 3] = 0.048769; AAmatrix[19][ 4] = 0.163831; 
	AAmatrix[19][ 5] = 0.112027; AAmatrix[19][ 6] = 0.205868; AAmatrix[19][ 7] = 0.082579; AAmatrix[19][ 8] = 0.068575; AAmatrix[19][ 9] = 3.654430; 
	AAmatrix[19][10] = 1.337571; AAmatrix[19][11] = 0.144587; AAmatrix[19][12] = 0.307309; AAmatrix[19][13] = 0.247329; AAmatrix[19][14] = 0.129315; 
	AAmatrix[19][15] = 0.127700; AAmatrix[19][16] = 0.740372; AAmatrix[19][17] = 0.022134; AAmatrix[19][18] = 0.125733; AAmatrix[19][19] = 0.000000; 

	AAPi[ 0] = 0.078837;
	AAPi[ 1] = 0.051238;
	AAPi[ 2] = 0.042313;
	AAPi[ 3] = 0.053066;
	AAPi[ 4] = 0.015175;
	AAPi[ 5] = 0.036713;
	AAPi[ 6] = 0.061924;
	AAPi[ 7] = 0.070852;
	AAPi[ 8] = 0.023082;
	AAPi[ 9] = 0.062056;
	AAPi[10] = 0.096371;
	AAPi[11] = 0.057324;
	AAPi[12] = 0.023771;
	AAPi[13] = 0.043296;
	AAPi[14] = 0.043911;
	AAPi[15] = 0.063403;
	AAPi[16] = 0.055897;
	AAPi[17] = 0.013272;
	AAPi[18] = 0.034399;
	AAPi[19] = 0.073101;

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==12)
	{

	/* 
		Blosum62 
		
		Henikoff, S., and J. G. Henikoff. (1992) Amino acid substitution matrices from protein blocks. Proc. Natl. Acad. Sci., U.S.A. 89:10915-10919.    

	*/
	AAmatrix[ 0][ 0] = 0.000000000000; AAmatrix[ 0][ 1] = 0.735790389698; AAmatrix[ 0][ 2] = 0.485391055466; AAmatrix[ 0][ 3] = 0.543161820899; AAmatrix[ 0][ 4] = 1.459995310470; 
	AAmatrix[ 0][ 5] = 1.199705704602; AAmatrix[ 0][ 6] = 1.170949042800; AAmatrix[ 0][ 7] = 1.955883574960; AAmatrix[ 0][ 8] = 0.716241444998; AAmatrix[ 0][ 9] = 0.605899003687; 
	AAmatrix[ 0][10] = 0.800016530518; AAmatrix[ 0][11] = 1.295201266783; AAmatrix[ 0][12] = 1.253758266664; AAmatrix[ 0][13] = 0.492964679748; AAmatrix[ 0][14] = 1.173275900924; 
	AAmatrix[ 0][15] = 4.325092687057; AAmatrix[ 0][16] = 1.729178019485; AAmatrix[ 0][17] = 0.465839367725; AAmatrix[ 0][18] = 0.718206697586; AAmatrix[ 0][19] = 2.187774522005; 
	AAmatrix[ 1][ 0] = 0.735790389698; AAmatrix[ 1][ 1] = 0.000000000000; AAmatrix[ 1][ 2] = 1.297446705134; AAmatrix[ 1][ 3] = 0.500964408555; AAmatrix[ 1][ 4] = 0.227826574209; 
	AAmatrix[ 1][ 5] = 3.020833610064; AAmatrix[ 1][ 6] = 1.360574190420; AAmatrix[ 1][ 7] = 0.418763308518; AAmatrix[ 1][ 8] = 1.456141166336; AAmatrix[ 1][ 9] = 0.232036445142; 
	AAmatrix[ 1][10] = 0.622711669692; AAmatrix[ 1][11] = 5.411115141489; AAmatrix[ 1][12] = 0.983692987457; AAmatrix[ 1][13] = 0.371644693209; AAmatrix[ 1][14] = 0.448133661718; 
	AAmatrix[ 1][15] = 1.122783104210; AAmatrix[ 1][16] = 0.914665954563; AAmatrix[ 1][17] = 0.426382310122; AAmatrix[ 1][18] = 0.720517441216; AAmatrix[ 1][19] = 0.438388343772; 
	AAmatrix[ 2][ 0] = 0.485391055466; AAmatrix[ 2][ 1] = 1.297446705134; AAmatrix[ 2][ 2] = 0.000000000000; AAmatrix[ 2][ 3] = 3.180100048216; AAmatrix[ 2][ 4] = 0.397358949897; 
	AAmatrix[ 2][ 5] = 1.839216146992; AAmatrix[ 2][ 6] = 1.240488508640; AAmatrix[ 2][ 7] = 1.355872344485; AAmatrix[ 2][ 8] = 2.414501434208; AAmatrix[ 2][ 9] = 0.283017326278; 
	AAmatrix[ 2][10] = 0.211888159615; AAmatrix[ 2][11] = 1.593137043457; AAmatrix[ 2][12] = 0.648441278787; AAmatrix[ 2][13] = 0.354861249223; AAmatrix[ 2][14] = 0.494887043702; 
	AAmatrix[ 2][15] = 2.904101656456; AAmatrix[ 2][16] = 1.898173634533; AAmatrix[ 2][17] = 0.191482046247; AAmatrix[ 2][18] = 0.538222519037; AAmatrix[ 2][19] = 0.312858797993; 
	AAmatrix[ 3][ 0] = 0.543161820899; AAmatrix[ 3][ 1] = 0.500964408555; AAmatrix[ 3][ 2] = 3.180100048216; AAmatrix[ 3][ 3] = 0.000000000000; AAmatrix[ 3][ 4] = 0.240836614802; 
	AAmatrix[ 3][ 5] = 1.190945703396; AAmatrix[ 3][ 6] = 3.761625208368; AAmatrix[ 3][ 7] = 0.798473248968; AAmatrix[ 3][ 8] = 0.778142664022; AAmatrix[ 3][ 9] = 0.418555732462; 
	AAmatrix[ 3][10] = 0.218131577594; AAmatrix[ 3][11] = 1.032447924952; AAmatrix[ 3][12] = 0.222621897958; AAmatrix[ 3][13] = 0.281730694207; AAmatrix[ 3][14] = 0.730628272998; 
	AAmatrix[ 3][15] = 1.582754142065; AAmatrix[ 3][16] = 0.934187509431; AAmatrix[ 3][17] = 0.145345046279; AAmatrix[ 3][18] = 0.261422208965; AAmatrix[ 3][19] = 0.258129289418; 
	AAmatrix[ 4][ 0] = 1.459995310470; AAmatrix[ 4][ 1] = 0.227826574209; AAmatrix[ 4][ 2] = 0.397358949897; AAmatrix[ 4][ 3] = 0.240836614802; AAmatrix[ 4][ 4] = 0.000000000000; 
	AAmatrix[ 4][ 5] = 0.329801504630; AAmatrix[ 4][ 6] = 0.140748891814; AAmatrix[ 4][ 7] = 0.418203192284; AAmatrix[ 4][ 8] = 0.354058109831; AAmatrix[ 4][ 9] = 0.774894022794; 
	AAmatrix[ 4][10] = 0.831842640142; AAmatrix[ 4][11] = 0.285078800906; AAmatrix[ 4][12] = 0.767688823480; AAmatrix[ 4][13] = 0.441337471187; AAmatrix[ 4][14] = 0.356008498769; 
	AAmatrix[ 4][15] = 1.197188415094; AAmatrix[ 4][16] = 1.119831358516; AAmatrix[ 4][17] = 0.527664418872; AAmatrix[ 4][18] = 0.470237733696; AAmatrix[ 4][19] = 1.116352478606; 
	AAmatrix[ 5][ 0] = 1.199705704602; AAmatrix[ 5][ 1] = 3.020833610064; AAmatrix[ 5][ 2] = 1.839216146992; AAmatrix[ 5][ 3] = 1.190945703396; AAmatrix[ 5][ 4] = 0.329801504630; 
	AAmatrix[ 5][ 5] = 0.000000000000; AAmatrix[ 5][ 6] = 5.528919177928; AAmatrix[ 5][ 7] = 0.609846305383; AAmatrix[ 5][ 8] = 2.435341131140; AAmatrix[ 5][ 9] = 0.236202451204; 
	AAmatrix[ 5][10] = 0.580737093181; AAmatrix[ 5][11] = 3.945277674515; AAmatrix[ 5][12] = 2.494896077113; AAmatrix[ 5][13] = 0.144356959750; AAmatrix[ 5][14] = 0.858570575674; 
	AAmatrix[ 5][15] = 1.934870924596; AAmatrix[ 5][16] = 1.277480294596; AAmatrix[ 5][17] = 0.758653808642; AAmatrix[ 5][18] = 0.958989742850; AAmatrix[ 5][19] = 0.530785790125; 
	AAmatrix[ 6][ 0] = 1.170949042800; AAmatrix[ 6][ 1] = 1.360574190420; AAmatrix[ 6][ 2] = 1.240488508640; AAmatrix[ 6][ 3] = 3.761625208368; AAmatrix[ 6][ 4] = 0.140748891814; 
	AAmatrix[ 6][ 5] = 5.528919177928; AAmatrix[ 6][ 6] = 0.000000000000; AAmatrix[ 6][ 7] = 0.423579992176; AAmatrix[ 6][ 8] = 1.626891056982; AAmatrix[ 6][ 9] = 0.186848046932; 
	AAmatrix[ 6][10] = 0.372625175087; AAmatrix[ 6][11] = 2.802427151679; AAmatrix[ 6][12] = 0.555415397470; AAmatrix[ 6][13] = 0.291409084165; AAmatrix[ 6][14] = 0.926563934846; 
	AAmatrix[ 6][15] = 1.769893238937; AAmatrix[ 6][16] = 1.071097236007; AAmatrix[ 6][17] = 0.407635648938; AAmatrix[ 6][18] = 0.596719300346; AAmatrix[ 6][19] = 0.524253846338; 
	AAmatrix[ 7][ 0] = 1.955883574960; AAmatrix[ 7][ 1] = 0.418763308518; AAmatrix[ 7][ 2] = 1.355872344485; AAmatrix[ 7][ 3] = 0.798473248968; AAmatrix[ 7][ 4] = 0.418203192284; 
	AAmatrix[ 7][ 5] = 0.609846305383; AAmatrix[ 7][ 6] = 0.423579992176; AAmatrix[ 7][ 7] = 0.000000000000; AAmatrix[ 7][ 8] = 0.539859124954; AAmatrix[ 7][ 9] = 0.189296292376; 
	AAmatrix[ 7][10] = 0.217721159236; AAmatrix[ 7][11] = 0.752042440303; AAmatrix[ 7][12] = 0.459436173579; AAmatrix[ 7][13] = 0.368166464453; AAmatrix[ 7][14] = 0.504086599527; 
	AAmatrix[ 7][15] = 1.509326253224; AAmatrix[ 7][16] = 0.641436011405; AAmatrix[ 7][17] = 0.508358924638; AAmatrix[ 7][18] = 0.308055737035; AAmatrix[ 7][19] = 0.253340790190; 
	AAmatrix[ 8][ 0] = 0.716241444998; AAmatrix[ 8][ 1] = 1.456141166336; AAmatrix[ 8][ 2] = 2.414501434208; AAmatrix[ 8][ 3] = 0.778142664022; AAmatrix[ 8][ 4] = 0.354058109831; 
	AAmatrix[ 8][ 5] = 2.435341131140; AAmatrix[ 8][ 6] = 1.626891056982; AAmatrix[ 8][ 7] = 0.539859124954; AAmatrix[ 8][ 8] = 0.000000000000; AAmatrix[ 8][ 9] = 0.252718447885; 
	AAmatrix[ 8][10] = 0.348072209797; AAmatrix[ 8][11] = 1.022507035889; AAmatrix[ 8][12] = 0.984311525359; AAmatrix[ 8][13] = 0.714533703928; AAmatrix[ 8][14] = 0.527007339151; 
	AAmatrix[ 8][15] = 1.117029762910; AAmatrix[ 8][16] = 0.585407090225; AAmatrix[ 8][17] = 0.301248600780; AAmatrix[ 8][18] = 4.218953969389; AAmatrix[ 8][19] = 0.201555971750; 
	AAmatrix[ 9][ 0] = 0.605899003687; AAmatrix[ 9][ 1] = 0.232036445142; AAmatrix[ 9][ 2] = 0.283017326278; AAmatrix[ 9][ 3] = 0.418555732462; AAmatrix[ 9][ 4] = 0.774894022794; 
	AAmatrix[ 9][ 5] = 0.236202451204; AAmatrix[ 9][ 6] = 0.186848046932; AAmatrix[ 9][ 7] = 0.189296292376; AAmatrix[ 9][ 8] = 0.252718447885; AAmatrix[ 9][ 9] = 0.000000000000; 
	AAmatrix[ 9][10] = 3.890963773304; AAmatrix[ 9][11] = 0.406193586642; AAmatrix[ 9][12] = 3.364797763104; AAmatrix[ 9][13] = 1.517359325954; AAmatrix[ 9][14] = 0.388355409206; 
	AAmatrix[ 9][15] = 0.357544412460; AAmatrix[ 9][16] = 1.179091197260; AAmatrix[ 9][17] = 0.341985787540; AAmatrix[ 9][18] = 0.674617093228; AAmatrix[ 9][19] = 8.311839405458; 
	AAmatrix[10][ 0] = 0.800016530518; AAmatrix[10][ 1] = 0.622711669692; AAmatrix[10][ 2] = 0.211888159615; AAmatrix[10][ 3] = 0.218131577594; AAmatrix[10][ 4] = 0.831842640142; 
	AAmatrix[10][ 5] = 0.580737093181; AAmatrix[10][ 6] = 0.372625175087; AAmatrix[10][ 7] = 0.217721159236; AAmatrix[10][ 8] = 0.348072209797; AAmatrix[10][ 9] = 3.890963773304; 
	AAmatrix[10][10] = 0.000000000000; AAmatrix[10][11] = 0.445570274261; AAmatrix[10][12] = 6.030559379572; AAmatrix[10][13] = 2.064839703237; AAmatrix[10][14] = 0.374555687471; 
	AAmatrix[10][15] = 0.352969184527; AAmatrix[10][16] = 0.915259857694; AAmatrix[10][17] = 0.691474634600; AAmatrix[10][18] = 0.811245856323; AAmatrix[10][19] = 2.231405688913; 
	AAmatrix[11][ 0] = 1.295201266783; AAmatrix[11][ 1] = 5.411115141489; AAmatrix[11][ 2] = 1.593137043457; AAmatrix[11][ 3] = 1.032447924952; AAmatrix[11][ 4] = 0.285078800906; 
	AAmatrix[11][ 5] = 3.945277674515; AAmatrix[11][ 6] = 2.802427151679; AAmatrix[11][ 7] = 0.752042440303; AAmatrix[11][ 8] = 1.022507035889; AAmatrix[11][ 9] = 0.406193586642; 
	AAmatrix[11][10] = 0.445570274261; AAmatrix[11][11] = 0.000000000000; AAmatrix[11][12] = 1.073061184332; AAmatrix[11][13] = 0.266924750511; AAmatrix[11][14] = 1.047383450722; 
	AAmatrix[11][15] = 1.752165917819; AAmatrix[11][16] = 1.303875200799; AAmatrix[11][17] = 0.332243040634; AAmatrix[11][18] = 0.717993486900; AAmatrix[11][19] = 0.498138475304; 
	AAmatrix[12][ 0] = 1.253758266664; AAmatrix[12][ 1] = 0.983692987457; AAmatrix[12][ 2] = 0.648441278787; AAmatrix[12][ 3] = 0.222621897958; AAmatrix[12][ 4] = 0.767688823480; 
	AAmatrix[12][ 5] = 2.494896077113; AAmatrix[12][ 6] = 0.555415397470; AAmatrix[12][ 7] = 0.459436173579; AAmatrix[12][ 8] = 0.984311525359; AAmatrix[12][ 9] = 3.364797763104; 
	AAmatrix[12][10] = 6.030559379572; AAmatrix[12][11] = 1.073061184332; AAmatrix[12][12] = 0.000000000000; AAmatrix[12][13] = 1.773855168830; AAmatrix[12][14] = 0.454123625103; 
	AAmatrix[12][15] = 0.918723415746; AAmatrix[12][16] = 1.488548053722; AAmatrix[12][17] = 0.888101098152; AAmatrix[12][18] = 0.951682162246; AAmatrix[12][19] = 2.575850755315; 
	AAmatrix[13][ 0] = 0.492964679748; AAmatrix[13][ 1] = 0.371644693209; AAmatrix[13][ 2] = 0.354861249223; AAmatrix[13][ 3] = 0.281730694207; AAmatrix[13][ 4] = 0.441337471187; 
	AAmatrix[13][ 5] = 0.144356959750; AAmatrix[13][ 6] = 0.291409084165; AAmatrix[13][ 7] = 0.368166464453; AAmatrix[13][ 8] = 0.714533703928; AAmatrix[13][ 9] = 1.517359325954; 
	AAmatrix[13][10] = 2.064839703237; AAmatrix[13][11] = 0.266924750511; AAmatrix[13][12] = 1.773855168830; AAmatrix[13][13] = 0.000000000000; AAmatrix[13][14] = 0.233597909629; 
	AAmatrix[13][15] = 0.540027644824; AAmatrix[13][16] = 0.488206118793; AAmatrix[13][17] = 2.074324893497; AAmatrix[13][18] = 6.747260430801; AAmatrix[13][19] = 0.838119610178; 
	AAmatrix[14][ 0] = 1.173275900924; AAmatrix[14][ 1] = 0.448133661718; AAmatrix[14][ 2] = 0.494887043702; AAmatrix[14][ 3] = 0.730628272998; AAmatrix[14][ 4] = 0.356008498769; 
	AAmatrix[14][ 5] = 0.858570575674; AAmatrix[14][ 6] = 0.926563934846; AAmatrix[14][ 7] = 0.504086599527; AAmatrix[14][ 8] = 0.527007339151; AAmatrix[14][ 9] = 0.388355409206; 
	AAmatrix[14][10] = 0.374555687471; AAmatrix[14][11] = 1.047383450722; AAmatrix[14][12] = 0.454123625103; AAmatrix[14][13] = 0.233597909629; AAmatrix[14][14] = 0.000000000000; 
	AAmatrix[14][15] = 1.169129577716; AAmatrix[14][16] = 1.005451683149; AAmatrix[14][17] = 0.252214830027; AAmatrix[14][18] = 0.369405319355; AAmatrix[14][19] = 0.496908410676; 
	AAmatrix[15][ 0] = 4.325092687057; AAmatrix[15][ 1] = 1.122783104210; AAmatrix[15][ 2] = 2.904101656456; AAmatrix[15][ 3] = 1.582754142065; AAmatrix[15][ 4] = 1.197188415094; 
	AAmatrix[15][ 5] = 1.934870924596; AAmatrix[15][ 6] = 1.769893238937; AAmatrix[15][ 7] = 1.509326253224; AAmatrix[15][ 8] = 1.117029762910; AAmatrix[15][ 9] = 0.357544412460; 
	AAmatrix[15][10] = 0.352969184527; AAmatrix[15][11] = 1.752165917819; AAmatrix[15][12] = 0.918723415746; AAmatrix[15][13] = 0.540027644824; AAmatrix[15][14] = 1.169129577716; 
	AAmatrix[15][15] = 0.000000000000; AAmatrix[15][16] = 5.151556292270; AAmatrix[15][17] = 0.387925622098; AAmatrix[15][18] = 0.796751520761; AAmatrix[15][19] = 0.561925457442; 
	AAmatrix[16][ 0] = 1.729178019485; AAmatrix[16][ 1] = 0.914665954563; AAmatrix[16][ 2] = 1.898173634533; AAmatrix[16][ 3] = 0.934187509431; AAmatrix[16][ 4] = 1.119831358516; 
	AAmatrix[16][ 5] = 1.277480294596; AAmatrix[16][ 6] = 1.071097236007; AAmatrix[16][ 7] = 0.641436011405; AAmatrix[16][ 8] = 0.585407090225; AAmatrix[16][ 9] = 1.179091197260; 
	AAmatrix[16][10] = 0.915259857694; AAmatrix[16][11] = 1.303875200799; AAmatrix[16][12] = 1.488548053722; AAmatrix[16][13] = 0.488206118793; AAmatrix[16][14] = 1.005451683149; 
	AAmatrix[16][15] = 5.151556292270; AAmatrix[16][16] = 0.000000000000; AAmatrix[16][17] = 0.513128126891; AAmatrix[16][18] = 0.801010243199; AAmatrix[16][19] = 2.253074051176; 
	AAmatrix[17][ 0] = 0.465839367725; AAmatrix[17][ 1] = 0.426382310122; AAmatrix[17][ 2] = 0.191482046247; AAmatrix[17][ 3] = 0.145345046279; AAmatrix[17][ 4] = 0.527664418872; 
	AAmatrix[17][ 5] = 0.758653808642; AAmatrix[17][ 6] = 0.407635648938; AAmatrix[17][ 7] = 0.508358924638; AAmatrix[17][ 8] = 0.301248600780; AAmatrix[17][ 9] = 0.341985787540; 
	AAmatrix[17][10] = 0.691474634600; AAmatrix[17][11] = 0.332243040634; AAmatrix[17][12] = 0.888101098152; AAmatrix[17][13] = 2.074324893497; AAmatrix[17][14] = 0.252214830027; 
	AAmatrix[17][15] = 0.387925622098; AAmatrix[17][16] = 0.513128126891; AAmatrix[17][17] = 0.000000000000; AAmatrix[17][18] = 4.054419006558; AAmatrix[17][19] = 0.266508731426; 
	AAmatrix[18][ 0] = 0.718206697586; AAmatrix[18][ 1] = 0.720517441216; AAmatrix[18][ 2] = 0.538222519037; AAmatrix[18][ 3] = 0.261422208965; AAmatrix[18][ 4] = 0.470237733696; 
	AAmatrix[18][ 5] = 0.958989742850; AAmatrix[18][ 6] = 0.596719300346; AAmatrix[18][ 7] = 0.308055737035; AAmatrix[18][ 8] = 4.218953969389; AAmatrix[18][ 9] = 0.674617093228; 
	AAmatrix[18][10] = 0.811245856323; AAmatrix[18][11] = 0.717993486900; AAmatrix[18][12] = 0.951682162246; AAmatrix[18][13] = 6.747260430801; AAmatrix[18][14] = 0.369405319355; 
	AAmatrix[18][15] = 0.796751520761; AAmatrix[18][16] = 0.801010243199; AAmatrix[18][17] = 4.054419006558; AAmatrix[18][18] = 0.000000000000; AAmatrix[18][19] = 1.000000000000; 
	AAmatrix[19][ 0] = 2.187774522005; AAmatrix[19][ 1] = 0.438388343772; AAmatrix[19][ 2] = 0.312858797993; AAmatrix[19][ 3] = 0.258129289418; AAmatrix[19][ 4] = 1.116352478606; 
	AAmatrix[19][ 5] = 0.530785790125; AAmatrix[19][ 6] = 0.524253846338; AAmatrix[19][ 7] = 0.253340790190; AAmatrix[19][ 8] = 0.201555971750; AAmatrix[19][ 9] = 8.311839405458; 
	AAmatrix[19][10] = 2.231405688913; AAmatrix[19][11] = 0.498138475304; AAmatrix[19][12] = 2.575850755315; AAmatrix[19][13] = 0.838119610178; AAmatrix[19][14] = 0.496908410676; 
	AAmatrix[19][15] = 0.561925457442; AAmatrix[19][16] = 2.253074051176; AAmatrix[19][17] = 0.266508731426; AAmatrix[19][18] = 1.000000000000; AAmatrix[19][19] = 0.000000000000; 	

	AAPi[ 0] = 0.074; 
	AAPi[ 1] = 0.052; 
	AAPi[ 2] = 0.045; 
	AAPi[ 3] = 0.054;
	AAPi[ 4] = 0.025; 
	AAPi[ 5] = 0.034; 
	AAPi[ 6] = 0.054; 
	AAPi[ 7] = 0.074;
	AAPi[ 8] = 0.026; 
	AAPi[ 9] = 0.068; 
	AAPi[10] = 0.099; 
	AAPi[11] = 0.058;
	AAPi[12] = 0.025; 
	AAPi[13] = 0.047; 
	AAPi[14] = 0.039; 
	AAPi[15] = 0.057;
	AAPi[16] = 0.051; 
	AAPi[17] = 0.013; 
	AAPi[18] = 0.032; 
	AAPi[19] = 0.073;

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==13)
	{

	/* 
		LG MODEL
		
		Le, S.Q. and Gascuel, O. (2008) LG: An Improved, General Amino-Acid Replacement Matrix. Molecular Biology and Evolution 25(7), 1307-1320 
	*/


	AAmatrix[0][1]= AAmatrix[1][0]=0.425093;
	AAmatrix[0][2]= AAmatrix[2][0]=0.276818;  AAmatrix[1][2]= AAmatrix[2][1]=0.751878;
	AAmatrix[0][3]= AAmatrix[3][0]=0.395144;  AAmatrix[1][3]= AAmatrix[3][1]=0.123954;  AAmatrix[2][3]= AAmatrix[3][2]=5.076149;
	AAmatrix[0][4]= AAmatrix[4][0]=2.489084;  AAmatrix[1][4]= AAmatrix[4][1]=0.534551;  AAmatrix[2][4]= AAmatrix[4][2]=0.528768;  AAmatrix[3][4]= AAmatrix[4][3]=0.062556;
	AAmatrix[0][5]= AAmatrix[5][0]=0.969894;  AAmatrix[1][5]= AAmatrix[5][1]=2.807908;  AAmatrix[2][5]= AAmatrix[5][2]=1.695752;  AAmatrix[3][5]= AAmatrix[5][3]=0.523386;  AAmatrix[4][5]= AAmatrix[5][4]=0.084808;
	AAmatrix[0][6]= AAmatrix[6][0]=1.038545;  AAmatrix[1][6]= AAmatrix[6][1]=0.363970;  AAmatrix[2][6]= AAmatrix[6][2]=0.541712;  AAmatrix[3][6]= AAmatrix[6][3]=5.243870;  AAmatrix[4][6]= AAmatrix[6][4]=0.003499;  AAmatrix[5][6]= AAmatrix[6][5]=4.128591;
	AAmatrix[0][7]= AAmatrix[7][0]=2.066040;  AAmatrix[1][7]= AAmatrix[7][1]=0.390192;  AAmatrix[2][7]= AAmatrix[7][2]=1.437645;  AAmatrix[3][7]= AAmatrix[7][3]=0.844926;  AAmatrix[4][7]= AAmatrix[7][4]=0.569265;  AAmatrix[5][7]= AAmatrix[7][5]=0.267959;  AAmatrix[6][7]= AAmatrix[7][6]=0.348847;
	AAmatrix[0][8]= AAmatrix[8][0]=0.358858;  AAmatrix[1][8]= AAmatrix[8][1]=2.426601;  AAmatrix[2][8]= AAmatrix[8][2]=4.509238;  AAmatrix[3][8]= AAmatrix[8][3]=0.927114;  AAmatrix[4][8]= AAmatrix[8][4]=0.640543;  AAmatrix[5][8]= AAmatrix[8][5]=4.813505;  AAmatrix[6][8]= AAmatrix[8][6]=0.423881;  AAmatrix[7][8]= AAmatrix[8][7]=0.311484;
	AAmatrix[0][9]= AAmatrix[9][0]=0.149830;  AAmatrix[1][9]= AAmatrix[9][1]=0.126991;  AAmatrix[2][9]= AAmatrix[9][2]=0.191503;  AAmatrix[3][9]= AAmatrix[9][3]=0.010690;  AAmatrix[4][9]= AAmatrix[9][4]=0.320627;  AAmatrix[5][9]= AAmatrix[9][5]=0.072854;  AAmatrix[6][9]= AAmatrix[9][6]=0.044265;  AAmatrix[7][9]= AAmatrix[9][7]=0.008705;  AAmatrix[8][9]= AAmatrix[9][8]=0.108882;
	AAmatrix[0][10]=AAmatrix[10][0]=0.395337; AAmatrix[1][10]=AAmatrix[10][1]=0.301848; AAmatrix[2][10]=AAmatrix[10][2]=0.068427; AAmatrix[3][10]=AAmatrix[10][3]=0.015076; AAmatrix[4][10]=AAmatrix[10][4]=0.594007; AAmatrix[5][10]=AAmatrix[10][5]=0.582457; AAmatrix[6][10]=AAmatrix[10][6]=0.069673; AAmatrix[7][10]=AAmatrix[10][7]=0.044261; AAmatrix[8][10]=AAmatrix[10][8]=0.366317; AAmatrix[9][10]=AAmatrix[10][9]=4.145067;
	AAmatrix[0][11]=AAmatrix[11][0]=0.536518; AAmatrix[1][11]=AAmatrix[11][1]=6.326067; AAmatrix[2][11]=AAmatrix[11][2]=2.145078; AAmatrix[3][11]=AAmatrix[11][3]=0.282959; AAmatrix[4][11]=AAmatrix[11][4]=0.013266; AAmatrix[5][11]=AAmatrix[11][5]=3.234294; AAmatrix[6][11]=AAmatrix[11][6]=1.807177; AAmatrix[7][11]=AAmatrix[11][7]=0.296636; AAmatrix[8][11]=AAmatrix[11][8]=0.697264; AAmatrix[9][11]=AAmatrix[11][9]=0.159069;  AAmatrix[10][11]=AAmatrix[11][10]=0.137500;
	AAmatrix[0][12]=AAmatrix[12][0]=1.124035; AAmatrix[1][12]=AAmatrix[12][1]=0.484133; AAmatrix[2][12]=AAmatrix[12][2]=0.371004; AAmatrix[3][12]=AAmatrix[12][3]=0.025548; AAmatrix[4][12]=AAmatrix[12][4]=0.893680; AAmatrix[5][12]=AAmatrix[12][5]=1.672569; AAmatrix[6][12]=AAmatrix[12][6]=0.173735; AAmatrix[7][12]=AAmatrix[12][7]=0.139538; AAmatrix[8][12]=AAmatrix[12][8]=0.442472; AAmatrix[9][12]=AAmatrix[12][9]=4.273607;  AAmatrix[10][12]=AAmatrix[12][10]=6.312358; AAmatrix[11][12]=AAmatrix[12][11]=0.656604;
	AAmatrix[0][13]=AAmatrix[13][0]=0.253701; AAmatrix[1][13]=AAmatrix[13][1]=0.052722; AAmatrix[2][13]=AAmatrix[13][2]=0.089525; AAmatrix[3][13]=AAmatrix[13][3]=0.017416; AAmatrix[4][13]=AAmatrix[13][4]=1.105251; AAmatrix[5][13]=AAmatrix[13][5]=0.035855; AAmatrix[6][13]=AAmatrix[13][6]=0.018811; AAmatrix[7][13]=AAmatrix[13][7]=0.089586; AAmatrix[8][13]=AAmatrix[13][8]=0.682139; AAmatrix[9][13]=AAmatrix[13][9]=1.112727;  AAmatrix[10][13]=AAmatrix[13][10]=2.592692; AAmatrix[11][13]=AAmatrix[13][11]=0.023918; AAmatrix[12][13]=AAmatrix[13][12]=1.798853;
	AAmatrix[0][14]=AAmatrix[14][0]=1.177651; AAmatrix[1][14]=AAmatrix[14][1]=0.332533; AAmatrix[2][14]=AAmatrix[14][2]=0.161787; AAmatrix[3][14]=AAmatrix[14][3]=0.394456; AAmatrix[4][14]=AAmatrix[14][4]=0.075382; AAmatrix[5][14]=AAmatrix[14][5]=0.624294; AAmatrix[6][14]=AAmatrix[14][6]=0.419409; AAmatrix[7][14]=AAmatrix[14][7]=0.196961; AAmatrix[8][14]=AAmatrix[14][8]=0.508851; AAmatrix[9][14]=AAmatrix[14][9]=0.078281;  AAmatrix[10][14]=AAmatrix[14][10]=0.249060; AAmatrix[11][14]=AAmatrix[14][11]=0.390322; AAmatrix[12][14]=AAmatrix[14][12]=0.099849; AAmatrix[13][14]=AAmatrix[14][13]=0.094464;
	AAmatrix[0][15]=AAmatrix[15][0]=4.727182; AAmatrix[1][15]=AAmatrix[15][1]=0.858151; AAmatrix[2][15]=AAmatrix[15][2]=4.008358; AAmatrix[3][15]=AAmatrix[15][3]=1.240275; AAmatrix[4][15]=AAmatrix[15][4]=2.784478; AAmatrix[5][15]=AAmatrix[15][5]=1.223828; AAmatrix[6][15]=AAmatrix[15][6]=0.611973; AAmatrix[7][15]=AAmatrix[15][7]=1.739990; AAmatrix[8][15]=AAmatrix[15][8]=0.990012; AAmatrix[9][15]=AAmatrix[15][9]=0.064105;  AAmatrix[10][15]=AAmatrix[15][10]=0.182287; AAmatrix[11][15]=AAmatrix[15][11]=0.748683; AAmatrix[12][15]=AAmatrix[15][12]=0.346960; AAmatrix[13][15]=AAmatrix[15][13]=0.361819; AAmatrix[14][15]=AAmatrix[15][14]=1.338132;
	AAmatrix[0][16]=AAmatrix[16][0]=2.139501; AAmatrix[1][16]=AAmatrix[16][1]=0.578987; AAmatrix[2][16]=AAmatrix[16][2]=2.000679; AAmatrix[3][16]=AAmatrix[16][3]=0.425860; AAmatrix[4][16]=AAmatrix[16][4]=1.143480; AAmatrix[5][16]=AAmatrix[16][5]=1.080136; AAmatrix[6][16]=AAmatrix[16][6]=0.604545; AAmatrix[7][16]=AAmatrix[16][7]=0.129836; AAmatrix[8][16]=AAmatrix[16][8]=0.584262; AAmatrix[9][16]=AAmatrix[16][9]=1.033739;  AAmatrix[10][16]=AAmatrix[16][10]=0.302936; AAmatrix[11][16]=AAmatrix[16][11]=1.136863; AAmatrix[12][16]=AAmatrix[16][12]=2.020366; AAmatrix[13][16]=AAmatrix[16][13]=0.165001; AAmatrix[14][16]=AAmatrix[16][14]=0.571468; AAmatrix[15][16]=AAmatrix[16][15]=6.472279;
	AAmatrix[0][17]=AAmatrix[17][0]=0.180717; AAmatrix[1][17]=AAmatrix[17][1]=0.593607; AAmatrix[2][17]=AAmatrix[17][2]=0.045376; AAmatrix[3][17]=AAmatrix[17][3]=0.029890; AAmatrix[4][17]=AAmatrix[17][4]=0.670128; AAmatrix[5][17]=AAmatrix[17][5]=0.236199; AAmatrix[6][17]=AAmatrix[17][6]=0.077852; AAmatrix[7][17]=AAmatrix[17][7]=0.268491; AAmatrix[8][17]=AAmatrix[17][8]=0.597054; AAmatrix[9][17]=AAmatrix[17][9]=0.111660;  AAmatrix[10][17]=AAmatrix[17][10]=0.619632; AAmatrix[11][17]=AAmatrix[17][11]=0.049906; AAmatrix[12][17]=AAmatrix[17][12]=0.696175; AAmatrix[13][17]=AAmatrix[17][13]=2.457121; AAmatrix[14][17]=AAmatrix[17][14]=0.095131; AAmatrix[15][17]=AAmatrix[17][15]=0.248862; AAmatrix[16][17]=AAmatrix[17][16]=0.140825;
	AAmatrix[0][18]=AAmatrix[18][0]=0.218959; AAmatrix[1][18]=AAmatrix[18][1]=0.314440; AAmatrix[2][18]=AAmatrix[18][2]=0.612025; AAmatrix[3][18]=AAmatrix[18][3]=0.135107; AAmatrix[4][18]=AAmatrix[18][4]=1.165532; AAmatrix[5][18]=AAmatrix[18][5]=0.257336; AAmatrix[6][18]=AAmatrix[18][6]=0.120037; AAmatrix[7][18]=AAmatrix[18][7]=0.054679; AAmatrix[8][18]=AAmatrix[18][8]=5.306834; AAmatrix[9][18]=AAmatrix[18][9]=0.232523;  AAmatrix[10][18]=AAmatrix[18][10]=0.299648; AAmatrix[11][18]=AAmatrix[18][11]=0.131932; AAmatrix[12][18]=AAmatrix[18][12]=0.481306; AAmatrix[13][18]=AAmatrix[18][13]=7.803902; AAmatrix[14][18]=AAmatrix[18][14]=0.089613; AAmatrix[15][18]=AAmatrix[18][15]=0.400547; AAmatrix[16][18]=AAmatrix[18][16]=0.245841; AAmatrix[17][18]=AAmatrix[18][17]=3.151815;
	AAmatrix[0][19]=AAmatrix[19][0]=2.547870; AAmatrix[1][19]=AAmatrix[19][1]=0.170887; AAmatrix[2][19]=AAmatrix[19][2]=0.083688; AAmatrix[3][19]=AAmatrix[19][3]=0.037967; AAmatrix[4][19]=AAmatrix[19][4]=1.959291; AAmatrix[5][19]=AAmatrix[19][5]=0.210332; AAmatrix[6][19]=AAmatrix[19][6]=0.245034; AAmatrix[7][19]=AAmatrix[19][7]=0.076701; AAmatrix[8][19]=AAmatrix[19][8]=0.119013; AAmatrix[9][19]=AAmatrix[19][9]=10.649107; AAmatrix[10][19]=AAmatrix[19][10]=1.702745; AAmatrix[11][19]=AAmatrix[19][11]=0.185202; AAmatrix[12][19]=AAmatrix[19][12]=1.898718; AAmatrix[13][19]=AAmatrix[19][13]=0.654683; AAmatrix[14][19]=AAmatrix[19][14]=0.296501; AAmatrix[15][19]=AAmatrix[19][15]=0.098369; AAmatrix[16][19]=AAmatrix[19][16]=2.188158; AAmatrix[17][19]=AAmatrix[19][17]=0.189510; AAmatrix[18][19]=AAmatrix[19][18]=0.249313;

	AAPi[0]=0.079066; AAPi[1]=0.055941; AAPi[2]=0.041977; AAPi[3]=0.053052; AAPi[4]=0.012937; AAPi[5]=0.040767; AAPi[6]=0.071586; AAPi[7]=0.057337; AAPi[8]=0.022355; AAPi[9]=0.062157; AAPi[10]=0.099081; AAPi[11]=0.064600; AAPi[12]=0.022951; AAPi[13]=0.042302; AAPi[14]=0.044040; AAPi[15]=0.061197; AAPi[16]=0.053287; AAPi[17]=0.012066; AAPi[18]=0.034155; AAPi[19]=0.069147;

	/*LG TEST*/
	/* // old version of matrix
	AAmatrixT[0][0] = 0;			AAmatrixT[0][1] = 0.449682;		AAmatrixT[0][2] = 0.267582;		AAmatrixT[0][3] = 0.401081;	
	AAmatrixT[0][4] = 2.31284;		AAmatrixT[0][5] = 0.944706;		AAmatrixT[0][6] = 1.16436;		AAmatrixT[0][7] = 2.10185;	
	AAmatrixT[0][8] = 0.341479;		AAmatrixT[0][9] = 0.122945;		AAmatrixT[0][10] = 0.391826;	AAmatrixT[0][11] = 0.556137;	
	AAmatrixT[0][12] = 1.0503;		AAmatrixT[0][13] = 0.237746;	AAmatrixT[0][14] = 1.23291;		AAmatrixT[0][15] = 4.65523;	
	AAmatrixT[0][16] = 1.98643;		AAmatrixT[0][17] = 0.179433;	AAmatrixT[0][18] = 0.223517;	AAmatrixT[0][19] = 2.36882;	
	AAmatrixT[1][0] = 0.449682;		AAmatrixT[1][1] = 0;			AAmatrixT[1][2] = 0.827348;		AAmatrixT[1][3] = 0.132811;	
	AAmatrixT[1][4] = 0.552587;		AAmatrixT[1][5] = 3.10941;		AAmatrixT[1][6] = 0.442407;		AAmatrixT[1][7] = 0.44398;	
	AAmatrixT[1][8] = 2.65765;		AAmatrixT[1][9] = 0.134451;		AAmatrixT[1][10] = 0.33036;		AAmatrixT[1][11] = 7.11437;	
	AAmatrixT[1][12] = 0.477124;	AAmatrixT[1][13] = 0.055544;	AAmatrixT[1][14] = 0.404818;	AAmatrixT[1][15] = 0.897892;	
	AAmatrixT[1][16] = 0.579784;	AAmatrixT[1][17] = 0.701255;	AAmatrixT[1][18] = 0.342216;	AAmatrixT[1][19] = 0.173721;	
	AAmatrixT[2][0] = 0.267582;		AAmatrixT[2][1] = 0.827348;		AAmatrixT[2][2] = 0;			AAmatrixT[2][3] = 5.921;	
	AAmatrixT[2][4] = 0.522133;		AAmatrixT[2][5] = 1.87744;		AAmatrixT[2][6] = 0.599223;		AAmatrixT[2][7] = 1.56619;	
	AAmatrixT[2][8] = 4.88956;		AAmatrixT[2][9] = 0.216069;		AAmatrixT[2][10] = 0.075149;	AAmatrixT[2][11] = 2.46334;	
	AAmatrixT[2][12] = 0.370061;	AAmatrixT[2][13] = 0.090929;	AAmatrixT[2][14] = 0.19063;		AAmatrixT[2][15] = 4.29942;	
	AAmatrixT[2][16] = 2.06149;		AAmatrixT[2][17] = 0.054722;	AAmatrixT[2][18] = 0.658002;	AAmatrixT[2][19] = 0.088856;	
	AAmatrixT[3][0] = 0.401081;		AAmatrixT[3][1] = 0.132811;		AAmatrixT[3][2] = 5.921;		AAmatrixT[3][3] = 0;	
	AAmatrixT[3][4] = 0.056428;		AAmatrixT[3][5] = 0.498202;		AAmatrixT[3][6] = 6.37423;		AAmatrixT[3][7] = 0.922928;	
	AAmatrixT[3][8] = 0.982202;		AAmatrixT[3][9] = 0.010922;		AAmatrixT[3][10] = 0.017176;	AAmatrixT[3][11] = 0.278545;	
	AAmatrixT[3][12] = 0.022762;	AAmatrixT[3][13] = 0.017714;	AAmatrixT[3][14] = 0.449817;	AAmatrixT[3][15] = 1.26822;	
	AAmatrixT[3][16] = 0.405969;	AAmatrixT[3][17] = 0.046559;	AAmatrixT[3][18] = 0.147235;	AAmatrixT[3][19] = 0.03872;	
	AAmatrixT[4][0] = 2.31284;		AAmatrixT[4][1] = 0.552587;		AAmatrixT[4][2] = 0.522133;		AAmatrixT[4][3] = 0.056428;	
	AAmatrixT[4][4] = 0;			AAmatrixT[4][5] = 0.080602;		AAmatrixT[4][6] = 0.00133;		AAmatrixT[4][7] = 0.529114;	
	AAmatrixT[4][8] = 0.593147;		AAmatrixT[4][9] = 0.262931;		AAmatrixT[4][10] = 0.541544;	AAmatrixT[4][11] = 0.003892;	
	AAmatrixT[4][12] = 0.773189;	AAmatrixT[4][13] = 0.950511;	AAmatrixT[4][14] = 0.076565;	AAmatrixT[4][15] = 2.60597;	
	AAmatrixT[4][16] = 0.993542;	AAmatrixT[4][17] = 0.659458;	AAmatrixT[4][18] = 1.09531;		AAmatrixT[4][19] = 1.74588;	
	AAmatrixT[5][0] = 0.944706;		AAmatrixT[5][1] = 3.10941;		AAmatrixT[5][2] = 1.87744;		AAmatrixT[5][3] = 0.498202;	
	AAmatrixT[5][4] = 0.080602;		AAmatrixT[5][5] = 0;			AAmatrixT[5][6] = 4.7998;		AAmatrixT[5][7] = 0.279365;	
	AAmatrixT[5][8] = 5.178;		AAmatrixT[5][9] = 0.073719;		AAmatrixT[5][10] = 0.61329;		AAmatrixT[5][11] = 3.46677;	
	AAmatrixT[5][12] = 1.65667;		AAmatrixT[5][13] = 0.033627;	AAmatrixT[5][14] = 0.69839;		AAmatrixT[5][15] = 1.2058;	
	AAmatrixT[5][16] = 1.02734;		AAmatrixT[5][17] = 0.249044;	AAmatrixT[5][18] = 0.244886;	AAmatrixT[5][19] = 0.204644;	
	AAmatrixT[6][0] = 1.16436;		AAmatrixT[6][1] = 0.442407;		AAmatrixT[6][2] = 0.599223;		AAmatrixT[6][3] = 6.37423;	
	AAmatrixT[6][4] = 0.00133;		AAmatrixT[6][5] = 4.7998;		AAmatrixT[6][6] = 0;			AAmatrixT[6][7] = 0.407773;	
	AAmatrixT[6][8] = 0.458209;		AAmatrixT[6][9] = 0.056153;		AAmatrixT[6][10] = 0.086633;	AAmatrixT[6][11] = 2.16893;	
	AAmatrixT[6][12] = 0.183748;	AAmatrixT[6][13] = 0.024362;	AAmatrixT[6][14] = 0.523437;	AAmatrixT[6][15] = 0.667092;	
	AAmatrixT[6][16] = 0.659097;	AAmatrixT[6][17] = 0.099542;	AAmatrixT[6][18] = 0.140547;	AAmatrixT[6][19] = 0.278624;	
	AAmatrixT[7][0] = 2.10185;		AAmatrixT[7][1] = 0.44398;		AAmatrixT[7][2] = 1.56619;		AAmatrixT[7][3] = 0.922928;	
	AAmatrixT[7][4] = 0.529114;		AAmatrixT[7][5] = 0.279365;		AAmatrixT[7][6] = 0.407773;		AAmatrixT[7][7] = 0;	
	AAmatrixT[7][8] = 0.30432;		AAmatrixT[7][9] = 0.008454;		AAmatrixT[7][10] = 0.047556;	AAmatrixT[7][11] = 0.313114;	
	AAmatrixT[7][12] = 0.137976;	AAmatrixT[7][13] = 0.080743;	AAmatrixT[7][14] = 0.226307;	AAmatrixT[7][15] = 1.78478;	
	AAmatrixT[7][16] = 0.114336;	AAmatrixT[7][17] = 0.292882;	AAmatrixT[7][18] = 0.056885;	AAmatrixT[7][19] = 0.075577;	
	AAmatrixT[8][0] = 0.341479;		AAmatrixT[8][1] = 2.65765;		AAmatrixT[8][2] = 4.88956;		AAmatrixT[8][3] = 0.982202;	
	AAmatrixT[8][4] = 0.593147;		AAmatrixT[8][5] = 5.178;		AAmatrixT[8][6] = 0.458209;		AAmatrixT[8][7] = 0.30432;	
	AAmatrixT[8][8] = 0;			AAmatrixT[8][9] = 0.106232;		AAmatrixT[8][10] = 0.363554;	AAmatrixT[8][11] = 0.682564;	
	AAmatrixT[8][12] = 0.395265;	AAmatrixT[8][13] = 0.616582;	AAmatrixT[8][14] = 0.545492;	AAmatrixT[8][15] = 0.947402;	
	AAmatrixT[8][16] = 0.526423;	AAmatrixT[8][17] = 0.559689;	AAmatrixT[8][18] = 5.44623;		AAmatrixT[8][19] = 0.108961;	
	AAmatrixT[9][0] = 0.122945;		AAmatrixT[9][1] = 0.134451;		AAmatrixT[9][2] = 0.216069;		AAmatrixT[9][3] = 0.010922;	
	AAmatrixT[9][4] = 0.262931;		AAmatrixT[9][5] = 0.073719;		AAmatrixT[9][6] = 0.056153;		AAmatrixT[9][7] = 0.008454;	
	AAmatrixT[9][8] = 0.106232;		AAmatrixT[9][9] = 0;			AAmatrixT[9][10] = 3.80151;		AAmatrixT[9][11] = 0.173179;	
	AAmatrixT[9][12] = 3.84902;		AAmatrixT[9][13] = 1.02066;		AAmatrixT[9][14] = 0.086269;	AAmatrixT[9][15] = 0.063251;	
	AAmatrixT[9][16] = 0.992803;	AAmatrixT[9][17] = 0.121839;	AAmatrixT[9][18] = 0.238891;	AAmatrixT[9][19] = 9.41677;	
	AAmatrixT[10][0] = 0.391826;	AAmatrixT[10][1] = 0.33036;		AAmatrixT[10][2] = 0.075149;	AAmatrixT[10][3] = 0.017176;	
	AAmatrixT[10][4] = 0.541544;	AAmatrixT[10][5] = 0.61329;		AAmatrixT[10][6] = 0.086633;	AAmatrixT[10][7] = 0.047556;	
	AAmatrixT[10][8] = 0.363554;	AAmatrixT[10][9] = 3.80151;		AAmatrixT[10][10] = 0;			AAmatrixT[10][11] = 0.145273;	
	AAmatrixT[10][12] = 5.83627;	AAmatrixT[10][13] = 2.42627;	AAmatrixT[10][14] = 0.265077;	AAmatrixT[10][15] = 0.184361;	
	AAmatrixT[10][16] = 0.286481;	AAmatrixT[10][17] = 0.649934;	AAmatrixT[10][18] = 0.292232;	AAmatrixT[10][19] = 1.51964;	
	AAmatrixT[11][0] = 0.556137;	AAmatrixT[11][1] = 7.11437;		AAmatrixT[11][2] = 2.46334;		AAmatrixT[11][3] = 0.278545;	
	AAmatrixT[11][4] = 0.003892;	AAmatrixT[11][5] = 3.46677;		AAmatrixT[11][6] = 2.16893;		AAmatrixT[11][7] = 0.313114;	
	AAmatrixT[11][8] = 0.682564;	AAmatrixT[11][9] = 0.173179;	AAmatrixT[11][10] = 0.145273;	AAmatrixT[11][11] = 0;	
	AAmatrixT[11][12] = 0.672252;	AAmatrixT[11][13] = 0.026721;	AAmatrixT[11][14] = 0.445474;	AAmatrixT[11][15] = 0.755746;	
	AAmatrixT[11][16] = 1.15218;	AAmatrixT[11][17] = 0.047995;	AAmatrixT[11][18] = 0.138336;	AAmatrixT[11][19] = 0.184432;	
	AAmatrixT[12][0] = 1.0503;		AAmatrixT[12][1] = 0.477124;	AAmatrixT[12][2] = 0.370061;	AAmatrixT[12][3] = 0.022762;	
	AAmatrixT[12][4] = 0.773189;	AAmatrixT[12][5] = 1.65667;		AAmatrixT[12][6] = 0.183748;	AAmatrixT[12][7] = 0.137976;	
	AAmatrixT[12][8] = 0.395265;	AAmatrixT[12][9] = 3.84902;		AAmatrixT[12][10] = 5.83627;	AAmatrixT[12][11] = 0.672252;	
	AAmatrixT[12][12] = 0;			AAmatrixT[12][13] = 1.62617;	AAmatrixT[12][14] = 0.096861;	AAmatrixT[12][15] = 0.319101;	
	AAmatrixT[12][16] = 1.86695;	AAmatrixT[12][17] = 0.660667;	AAmatrixT[12][18] = 0.436403;	AAmatrixT[12][19] = 1.59505;	
	AAmatrixT[13][0] = 0.237746;	AAmatrixT[13][1] = 0.055544;	AAmatrixT[13][2] = 0.090929;	AAmatrixT[13][3] = 0.017714;	
	AAmatrixT[13][4] = 0.950511;	AAmatrixT[13][5] = 0.033627;	AAmatrixT[13][6] = 0.024362;	AAmatrixT[13][7] = 0.080743;	
	AAmatrixT[13][8] = 0.616582;	AAmatrixT[13][9] = 1.02066;		AAmatrixT[13][10] = 2.42627;	AAmatrixT[13][11] = 0.026721;	
	AAmatrixT[13][12] = 1.62617;	AAmatrixT[13][13] = 0;			AAmatrixT[13][14] = 0.104849;	AAmatrixT[13][15] = 0.355654;	
	AAmatrixT[13][16] = 0.145526;	AAmatrixT[13][17] = 2.42582;	AAmatrixT[13][18] = 7.59878;	AAmatrixT[13][19] = 0.578417;	
	AAmatrixT[14][0] = 1.23291;		AAmatrixT[14][1] = 0.404818;	AAmatrixT[14][2] = 0.19063;		AAmatrixT[14][3] = 0.449817;	
	AAmatrixT[14][4] = 0.076565;	AAmatrixT[14][5] = 0.69839;		AAmatrixT[14][6] = 0.523437;	AAmatrixT[14][7] = 0.226307;	
	AAmatrixT[14][8] = 0.545492;	AAmatrixT[14][9] = 0.086269;	AAmatrixT[14][10] = 0.265077;	AAmatrixT[14][11] = 0.445474;	
	AAmatrixT[14][12] = 0.096861;	AAmatrixT[14][13] = 0.104849;	AAmatrixT[14][14] = 0;			AAmatrixT[14][15] = 1.42481;	
	AAmatrixT[14][16] = 0.592443;	AAmatrixT[14][17] = 0.118287;	AAmatrixT[14][18] = 0.109774;	AAmatrixT[14][19] = 0.302548;	
	AAmatrixT[15][0] = 4.65523;		AAmatrixT[15][1] = 0.897892;	AAmatrixT[15][2] = 4.29942;		AAmatrixT[15][3] = 1.26822;	
	AAmatrixT[15][4] = 2.60597;		AAmatrixT[15][5] = 1.2058;		AAmatrixT[15][6] = 0.667092;	AAmatrixT[15][7] = 1.78478;	
	AAmatrixT[15][8] = 0.947402;	AAmatrixT[15][9] = 0.063251;	AAmatrixT[15][10] = 0.184361;	AAmatrixT[15][11] = 0.755746;	
	AAmatrixT[15][12] = 0.319101;	AAmatrixT[15][13] = 0.355654;	AAmatrixT[15][14] = 1.42481;	AAmatrixT[15][15] = 0;	
	AAmatrixT[15][16] = 6.26607;	AAmatrixT[15][17] = 0.267487;	AAmatrixT[15][18] = 0.407468;	AAmatrixT[15][19] = 0.062285;	
	AAmatrixT[16][0] = 1.98643;		AAmatrixT[16][1] = 0.579784;	AAmatrixT[16][2] = 2.06149;		AAmatrixT[16][3] = 0.405969;	
	AAmatrixT[16][4] = 0.993542;	AAmatrixT[16][5] = 1.02734;		AAmatrixT[16][6] = 0.659097;	AAmatrixT[16][7] = 0.114336;	
	AAmatrixT[16][8] = 0.526423;	AAmatrixT[16][9] = 0.992803;	AAmatrixT[16][10] = 0.286481;	AAmatrixT[16][11] = 1.15218;	
	AAmatrixT[16][12] = 1.86695;	AAmatrixT[16][13] = 0.145526;	AAmatrixT[16][14] = 0.592443;	AAmatrixT[16][15] = 6.26607;	
	AAmatrixT[16][16] = 0;			AAmatrixT[16][17] = 0.144967;	AAmatrixT[16][18] = 0.236493;	AAmatrixT[16][19] = 1.94732;	
	AAmatrixT[17][0] = 0.179433;	AAmatrixT[17][1] = 0.701255;	AAmatrixT[17][2] = 0.054722;	AAmatrixT[17][3] = 0.046559;	
	AAmatrixT[17][4] = 0.659458;	AAmatrixT[17][5] = 0.249044;	AAmatrixT[17][6] = 0.099542;	AAmatrixT[17][7] = 0.292882;	
	AAmatrixT[17][8] = 0.559689;	AAmatrixT[17][9] = 0.121839;	AAmatrixT[17][10] = 0.649934;	AAmatrixT[17][11] = 0.047995;	
	AAmatrixT[17][12] = 0.660667;	AAmatrixT[17][13] = 2.42582;	AAmatrixT[17][14] = 0.118287;	AAmatrixT[17][15] = 0.267487;	
	AAmatrixT[17][16] = 0.144967;	AAmatrixT[17][17] = 0;			AAmatrixT[17][18] = 3.34452;	AAmatrixT[17][19] = 0.201078;	
	AAmatrixT[18][0] = 0.223517;	AAmatrixT[18][1] = 0.342216;	AAmatrixT[18][2] = 0.658002;	AAmatrixT[18][3] = 0.147235;	
	AAmatrixT[18][4] = 1.09531;		AAmatrixT[18][5] = 0.244886;	AAmatrixT[18][6] = 0.140547;	AAmatrixT[18][7] = 0.056885;	
	AAmatrixT[18][8] = 5.44623;		AAmatrixT[18][9] = 0.238891;	AAmatrixT[18][10] = 0.292232;	AAmatrixT[18][11] = 0.138336;	
	AAmatrixT[18][12] = 0.436403;	AAmatrixT[18][13] = 7.59878;	AAmatrixT[18][14] = 0.109774;	AAmatrixT[18][15] = 0.407468;	
	AAmatrixT[18][16] = 0.236493;	AAmatrixT[18][17] = 3.34452;	AAmatrixT[18][18] = 0;			AAmatrixT[18][19] = 0.235819;	
	AAmatrixT[19][0] = 2.36882;		AAmatrixT[19][1] = 0.173721;	AAmatrixT[19][2] = 0.088856;	AAmatrixT[19][3] = 0.03872;	
	AAmatrixT[19][4] = 1.74588;		AAmatrixT[19][5] = 0.204644;	AAmatrixT[19][6] = 0.278624;	AAmatrixT[19][7] = 0.075577;	
	AAmatrixT[19][8] = 0.108961;	AAmatrixT[19][9] = 9.41677;		AAmatrixT[19][10] = 1.51964;	AAmatrixT[19][11] = 0.184432;	
	AAmatrixT[19][12] = 1.59505;	AAmatrixT[19][13] = 0.578417;	AAmatrixT[19][14] = 0.302548;	AAmatrixT[19][15] = 0.062285;	
	AAmatrixT[19][16] = 1.94732;	AAmatrixT[19][17] = 0.201078;	AAmatrixT[19][18] = 0.235819;	AAmatrixT[19][19] = 0;	

	AAPiT[0] = 0.079611;
	AAPiT[1] = 0.053191;
	AAPiT[2] = 0.039948;
	AAPiT[3] = 0.050634;
	AAPiT[4] = 0.01359;
	AAPiT[5] = 0.038611;
	AAPiT[6] = 0.066539;
	AAPiT[7] = 0.059913;
	AAPiT[8] = 0.021738;
	AAPiT[9] = 0.063589;
	AAPiT[10] = 0.105134;
	AAPiT[11] = 0.061845;
	AAPiT[12] = 0.02299;
	AAPiT[13] = 0.044365;
	AAPiT[14] = 0.044909;
	AAPiT[15] = 0.059477;
	AAPiT[16] = 0.054114;
	AAPiT[17] = 0.012588;
	AAPiT[18] = 0.035709;
	AAPiT[19] = 0.071505;
	
		for(int gb1=0; gb1<20; gb1++)
		{
			if(AAPiT[gb1]!=AAPiT[gb1]) {cout<<"ERROR in base frequency translation check in LG model"<<endl; gb1--;}
			
			for(int gb2=0; gb2<20; gb2++)
			{
				if(AAmatrix[gb1][gb2]!=AAmatrixT[gb1][gb2]) {cout<<"ERROR in substitution matrix translation check in LG model"<<endl; gb2--;}
			}
		}
	*/

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==14)
	{

	
	/*
		HIVb
		
		Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL. (2007) HIV-Specific Probabilistic Models 
		of Protein Evolution. PLoS ONE. Jun 6;2:e503.
	*/
	AAmatrix[0][0] = 0;				AAmatrix[0][1] = 0.307507;		AAmatrix[0][2] = 0.005;			AAmatrix[0][3] = 1.45504;	
	AAmatrix[0][4] = 0.123758;		AAmatrix[0][5] = 0.0551128;		AAmatrix[0][6] = 1.48135;		AAmatrix[0][7] = 2.13536;	
	AAmatrix[0][8] = 0.0847613;		AAmatrix[0][9] = 0.005;			AAmatrix[0][10] = 0.215256;		AAmatrix[0][11] = 0.005;	
	AAmatrix[0][12] = 0.0186643;	AAmatrix[0][13] = 0.0141269;	AAmatrix[0][14] = 2.12217;		AAmatrix[0][15] = 2.46633;	
	AAmatrix[0][16] = 15.9183;		AAmatrix[0][17] = 0.005;		AAmatrix[0][18] = 0.005;		AAmatrix[0][19] = 7.61428;	
	AAmatrix[1][0] = 0.307507;		AAmatrix[1][1] = 0;				AAmatrix[1][2] = 0.295543;		AAmatrix[1][3] = 0.005;	
	AAmatrix[1][4] = 0.351721;		AAmatrix[1][5] = 3.4215;		AAmatrix[1][6] = 0.0749218;		AAmatrix[1][7] = 3.65345;	
	AAmatrix[1][8] = 9.04044;		AAmatrix[1][9] = 0.677289;		AAmatrix[1][10] = 0.701427;		AAmatrix[1][11] = 20.45;	
	AAmatrix[1][12] = 2.51394;		AAmatrix[1][13] = 0.005;		AAmatrix[1][14] = 1.28355;		AAmatrix[1][15] = 3.4791;	
	AAmatrix[1][16] = 2.86868;		AAmatrix[1][17] = 0.991338;		AAmatrix[1][18] = 0.00991826;	AAmatrix[1][19] = 0.0812454;	
	AAmatrix[2][0] = 0.005;			AAmatrix[2][1] = 0.295543;		AAmatrix[2][2] = 0;				AAmatrix[2][3] = 17.6612;	
	AAmatrix[2][4] = 0.0860642;		AAmatrix[2][5] = 0.672052;		AAmatrix[2][6] = 0.0792633;		AAmatrix[2][7] = 0.323401;	
	AAmatrix[2][8] = 7.64585;		AAmatrix[2][9] = 0.680565;		AAmatrix[2][10] = 0.005;		AAmatrix[2][11] = 7.90443;	
	AAmatrix[2][12] = 0.005;		AAmatrix[2][13] = 0.005;		AAmatrix[2][14] = 0.00739578;	AAmatrix[2][15] = 13.1447;	
	AAmatrix[2][16] = 6.88667;		AAmatrix[2][17] = 0.005;		AAmatrix[2][18] = 1.76417;		AAmatrix[2][19] = 0.026656;	
	AAmatrix[3][0] = 1.45504;		AAmatrix[3][1] = 0.005;			AAmatrix[3][2] = 17.6612;		AAmatrix[3][3] = 0;	
	AAmatrix[3][4] = 0.005;			AAmatrix[3][5] = 0.005;			AAmatrix[3][6] = 10.5872;		AAmatrix[3][7] = 2.83806;	
	AAmatrix[3][8] = 1.9169;		AAmatrix[3][9] = 0.0176792;		AAmatrix[3][10] = 0.00876048;	AAmatrix[3][11] = 0.005;	
	AAmatrix[3][12] = 0.005;		AAmatrix[3][13] = 0.005;		AAmatrix[3][14] = 0.0342658;	AAmatrix[3][15] = 0.52823;	
	AAmatrix[3][16] = 0.274724;		AAmatrix[3][17] = 0.005;		AAmatrix[3][18] = 0.674653;		AAmatrix[3][19] = 1.04793;	
	AAmatrix[4][0] = 0.123758;		AAmatrix[4][1] = 0.351721;		AAmatrix[4][2] = 0.0860642;		AAmatrix[4][3] = 0.005;	
	AAmatrix[4][4] = 0;				AAmatrix[4][5] = 0.005;			AAmatrix[4][6] = 0.005;			AAmatrix[4][7] = 0.897871;	
	AAmatrix[4][8] = 0.240073;		AAmatrix[4][9] = 0.005;			AAmatrix[4][10] = 0.129777;		AAmatrix[4][11] = 0.005;	
	AAmatrix[4][12] = 0.005;		AAmatrix[4][13] = 9.29815;		AAmatrix[4][14] = 0.005;		AAmatrix[4][15] = 4.69314;	
	AAmatrix[4][16] = 0.739969;		AAmatrix[4][17] = 2.63277;		AAmatrix[4][18] = 7.57932;		AAmatrix[4][19] = 0.420027;	
	AAmatrix[5][0] = 0.0551128;		AAmatrix[5][1] = 3.4215;		AAmatrix[5][2] = 0.672052;		AAmatrix[5][3] = 0.005;	
	AAmatrix[5][4] = 0.005;			AAmatrix[5][5] = 0;				AAmatrix[5][6] = 2.5602;		AAmatrix[5][7] = 0.0619137;	
	AAmatrix[5][8] = 7.05545;		AAmatrix[5][9] = 0.005;			AAmatrix[5][10] = 1.49456;		AAmatrix[5][11] = 6.54737;	
	AAmatrix[5][12] = 0.303676;		AAmatrix[5][13] = 0.005;		AAmatrix[5][14] = 4.47211;		AAmatrix[5][15] = 0.116311;	
	AAmatrix[5][16] = 0.243589;		AAmatrix[5][17] = 0.026656;		AAmatrix[5][18] = 0.113033;		AAmatrix[5][19] = 0.0209153;	
	AAmatrix[6][0] = 1.48135;		AAmatrix[6][1] = 0.0749218;		AAmatrix[6][2] = 0.0792633;		AAmatrix[6][3] = 10.5872;	
	AAmatrix[6][4] = 0.005;			AAmatrix[6][5] = 2.5602;		AAmatrix[6][6] = 0;				AAmatrix[6][7] = 3.92775;	
	AAmatrix[6][8] = 0.11974;		AAmatrix[6][9] = 0.00609079;	AAmatrix[6][10] = 0.005;		AAmatrix[6][11] = 4.61482;	
	AAmatrix[6][12] = 0.175789;		AAmatrix[6][13] = 0.005;		AAmatrix[6][14] = 0.0120226;	AAmatrix[6][15] = 0.005;	
	AAmatrix[6][16] = 0.289774;		AAmatrix[6][17] = 0.005;		AAmatrix[6][18] = 0.0792633;	AAmatrix[6][19] = 1.02847;	
	AAmatrix[7][0] = 2.13536;		AAmatrix[7][1] = 3.65345;		AAmatrix[7][2] = 0.323401;		AAmatrix[7][3] = 2.83806;	
	AAmatrix[7][4] = 0.897871;		AAmatrix[7][5] = 0.0619137;		AAmatrix[7][6] = 3.92775;		AAmatrix[7][7] = 0;	
	AAmatrix[7][8] = 0.005;			AAmatrix[7][9] = 0.005;			AAmatrix[7][10] = 0.005;		AAmatrix[7][11] = 0.521705;	
	AAmatrix[7][12] = 0.005;		AAmatrix[7][13] = 0.291561;		AAmatrix[7][14] = 0.005;		AAmatrix[7][15] = 4.38041;	
	AAmatrix[7][16] = 0.369615;		AAmatrix[7][17] = 1.21674;		AAmatrix[7][18] = 0.005;		AAmatrix[7][19] = 0.953155;	
	AAmatrix[8][0] = 0.0847613;		AAmatrix[8][1] = 9.04044;		AAmatrix[8][2] = 7.64585;		AAmatrix[8][3] = 1.9169;	
	AAmatrix[8][4] = 0.240073;		AAmatrix[8][5] = 7.05545;		AAmatrix[8][6] = 0.11974;		AAmatrix[8][7] = 0.005;	
	AAmatrix[8][8] = 0;				AAmatrix[8][9] = 0.103111;		AAmatrix[8][10] = 1.74171;		AAmatrix[8][11] = 0.005;	
	AAmatrix[8][12] = 0.005;		AAmatrix[8][13] = 0.145558;		AAmatrix[8][14] = 2.45318;		AAmatrix[8][15] = 0.382747;	
	AAmatrix[8][16] = 0.711594;		AAmatrix[8][17] = 0.0695179;	AAmatrix[8][18] = 18.6943;		AAmatrix[8][19] = 0.005;	
	AAmatrix[9][0] = 0.005;			AAmatrix[9][1] = 0.677289;		AAmatrix[9][2] = 0.680565;		AAmatrix[9][3] = 0.0176792;	
	AAmatrix[9][4] = 0.005;			AAmatrix[9][5] = 0.005;			AAmatrix[9][6] = 0.00609079;	AAmatrix[9][7] = 0.005;	
	AAmatrix[9][8] = 0.103111;		AAmatrix[9][9] = 0;				AAmatrix[9][10] = 5.95879;		AAmatrix[9][11] = 0.322319;	
	AAmatrix[9][12] = 11.2065;		AAmatrix[9][13] = 3.39836;		AAmatrix[9][14] = 0.0410593;	AAmatrix[9][15] = 1.21803;	
	AAmatrix[9][16] = 8.61217;		AAmatrix[9][17] = 0.005;		AAmatrix[9][18] = 0.148168;		AAmatrix[9][19] = 17.7389;	
	AAmatrix[10][0] = 0.215256;		AAmatrix[10][1] = 0.701427;		AAmatrix[10][2] = 0.005;		AAmatrix[10][3] = 0.00876048;	
	AAmatrix[10][4] = 0.129777;		AAmatrix[10][5] = 1.49456;		AAmatrix[10][6] = 0.005;		AAmatrix[10][7] = 0.005;	
	AAmatrix[10][8] = 1.74171;		AAmatrix[10][9] = 5.95879;		AAmatrix[10][10] = 0;			AAmatrix[10][11] = 0.0814995;	
	AAmatrix[10][12] = 5.31961;		AAmatrix[10][13] = 8.52484;		AAmatrix[10][14] = 2.07757;		AAmatrix[10][15] = 0.927656;	
	AAmatrix[10][16] = 0.0437673;	AAmatrix[10][17] = 0.748843;	AAmatrix[10][18] = 0.111986;	AAmatrix[10][19] = 1.41036;	
	AAmatrix[11][0] = 0.005;		AAmatrix[11][1] = 20.45;		AAmatrix[11][2] = 7.90443;		AAmatrix[11][3] = 0.005;	
	AAmatrix[11][4] = 0.005;		AAmatrix[11][5] = 6.54737;		AAmatrix[11][6] = 4.61482;		AAmatrix[11][7] = 0.521705;	
	AAmatrix[11][8] = 0.005;		AAmatrix[11][9] = 0.322319;		AAmatrix[11][10] = 0.0814995;	AAmatrix[11][11] = 0;	
	AAmatrix[11][12] = 1.28246;		AAmatrix[11][13] = 0.0342658;	AAmatrix[11][14] = 0.0313862;	AAmatrix[11][15] = 0.504111;	
	AAmatrix[11][16] = 4.67142;		AAmatrix[11][17] = 0.005;		AAmatrix[11][18] = 0.005;		AAmatrix[11][19] = 0.265829;	
	AAmatrix[12][0] = 0.0186643;	AAmatrix[12][1] = 2.51394;		AAmatrix[12][2] = 0.005;		AAmatrix[12][3] = 0.005;	
	AAmatrix[12][4] = 0.005;		AAmatrix[12][5] = 0.303676;		AAmatrix[12][6] = 0.175789;		AAmatrix[12][7] = 0.005;	
	AAmatrix[12][8] = 0.005;		AAmatrix[12][9] = 11.2065;		AAmatrix[12][10] = 5.31961;		AAmatrix[12][11] = 1.28246;	
	AAmatrix[12][12] = 0;			AAmatrix[12][13] = 0.188025;	AAmatrix[12][14] = 0.005;		AAmatrix[12][15] = 0.005;	
	AAmatrix[12][16] = 4.94026;		AAmatrix[12][17] = 0.089078;	AAmatrix[12][18] = 0.005;		AAmatrix[12][19] = 6.8532;	
	AAmatrix[13][0] = 0.0141269;	AAmatrix[13][1] = 0.005;		AAmatrix[13][2] = 0.005;		AAmatrix[13][3] = 0.005;	
	AAmatrix[13][4] = 9.29815;		AAmatrix[13][5] = 0.005;		AAmatrix[13][6] = 0.005;		AAmatrix[13][7] = 0.291561;	
	AAmatrix[13][8] = 0.145558;		AAmatrix[13][9] = 3.39836;		AAmatrix[13][10] = 8.52484;		AAmatrix[13][11] = 0.0342658;	
	AAmatrix[13][12] = 0.188025;	AAmatrix[13][13] = 0;			AAmatrix[13][14] = 0.005;		AAmatrix[13][15] = 0.956472;	
	AAmatrix[13][16] = 0.0141269;	AAmatrix[13][17] = 0.829343;	AAmatrix[13][18] = 15.34;		AAmatrix[13][19] = 0.723274;	
	AAmatrix[14][0] = 2.12217;		AAmatrix[14][1] = 1.28355;		AAmatrix[14][2] = 0.00739578;	AAmatrix[14][3] = 0.0342658;	
	AAmatrix[14][4] = 0.005;		AAmatrix[14][5] = 4.47211;		AAmatrix[14][6] = 0.0120226;	AAmatrix[14][7] = 0.005;	
	AAmatrix[14][8] = 2.45318;		AAmatrix[14][9] = 0.0410593;	AAmatrix[14][10] = 2.07757;		AAmatrix[14][11] = 0.0313862;	
	AAmatrix[14][12] = 0.005;		AAmatrix[14][13] = 0.005;		AAmatrix[14][14] = 0;			AAmatrix[14][15] = 5.37762;	
	AAmatrix[14][16] = 2.01417;		AAmatrix[14][17] = 0.0444506;	AAmatrix[14][18] = 0.0304381;	AAmatrix[14][19] = 0.005;	
	AAmatrix[15][0] = 2.46633;		AAmatrix[15][1] = 3.4791;		AAmatrix[15][2] = 13.1447;		AAmatrix[15][3] = 0.52823;	
	AAmatrix[15][4] = 4.69314;		AAmatrix[15][5] = 0.116311;		AAmatrix[15][6] = 0.005;		AAmatrix[15][7] = 4.38041;	
	AAmatrix[15][8] = 0.382747;		AAmatrix[15][9] = 1.21803;		AAmatrix[15][10] = 0.927656;	AAmatrix[15][11] = 0.504111;	
	AAmatrix[15][12] = 0.005;		AAmatrix[15][13] = 0.956472;	AAmatrix[15][14] = 5.37762;		AAmatrix[15][15] = 0;	
	AAmatrix[15][16] = 8.93107;		AAmatrix[15][17] = 0.0248728;	AAmatrix[15][18] = 0.648024;	AAmatrix[15][19] = 0.0749218;	
	AAmatrix[16][0] = 15.9183;		AAmatrix[16][1] = 2.86868;		AAmatrix[16][2] = 6.88667;		AAmatrix[16][3] = 0.274724;	
	AAmatrix[16][4] = 0.739969;		AAmatrix[16][5] = 0.243589;		AAmatrix[16][6] = 0.289774;		AAmatrix[16][7] = 0.369615;	
	AAmatrix[16][8] = 0.711594;		AAmatrix[16][9] = 8.61217;		AAmatrix[16][10] = 0.0437673;	AAmatrix[16][11] = 4.67142;	
	AAmatrix[16][12] = 4.94026;		AAmatrix[16][13] = 0.0141269;	AAmatrix[16][14] = 2.01417;		AAmatrix[16][15] = 8.93107;	
	AAmatrix[16][16] = 0;			AAmatrix[16][17] = 0.005;		AAmatrix[16][18] = 0.105652;	AAmatrix[16][19] = 0.709226;	
	AAmatrix[17][0] = 0.005;		AAmatrix[17][1] = 0.991338;		AAmatrix[17][2] = 0.005;		AAmatrix[17][3] = 0.005;	
	AAmatrix[17][4] = 2.63277;		AAmatrix[17][5] = 0.026656;		AAmatrix[17][6] = 0.005;		AAmatrix[17][7] = 1.21674;	
	AAmatrix[17][8] = 0.0695179;	AAmatrix[17][9] = 0.005;		AAmatrix[17][10] = 0.748843;	AAmatrix[17][11] = 0.005;	
	AAmatrix[17][12] = 0.089078;	AAmatrix[17][13] = 0.829343;	AAmatrix[17][14] = 0.0444506;	AAmatrix[17][15] = 0.0248728;	
	AAmatrix[17][16] = 0.005;		AAmatrix[17][17] = 0;			AAmatrix[17][18] = 1.28022;		AAmatrix[17][19] = 0.005;	
	AAmatrix[18][0] = 0.005;		AAmatrix[18][1] = 0.00991826;	AAmatrix[18][2] = 1.76417;		AAmatrix[18][3] = 0.674653;	
	AAmatrix[18][4] = 7.57932;		AAmatrix[18][5] = 0.113033;		AAmatrix[18][6] = 0.0792633;	AAmatrix[18][7] = 0.005;	
	AAmatrix[18][8] = 18.6943;		AAmatrix[18][9] = 0.148168;		AAmatrix[18][10] = 0.111986;	AAmatrix[18][11] = 0.005;	
	AAmatrix[18][12] = 0.005;		AAmatrix[18][13] = 15.34;		AAmatrix[18][14] = 0.0304381;	AAmatrix[18][15] = 0.648024;	
	AAmatrix[18][16] = 0.105652;	AAmatrix[18][17] = 1.28022;		AAmatrix[18][18] = 0;			AAmatrix[18][19] = 0.0410593;	
	AAmatrix[19][0] = 7.61428;		AAmatrix[19][1] = 0.0812454;	AAmatrix[19][2] = 0.026656;		AAmatrix[19][3] = 1.04793;	
	AAmatrix[19][4] = 0.420027;		AAmatrix[19][5] = 0.0209153;	AAmatrix[19][6] = 1.02847;		AAmatrix[19][7] = 0.953155;	
	AAmatrix[19][8] = 0.005;		AAmatrix[19][9] = 17.7389;		AAmatrix[19][10] = 1.41036;		AAmatrix[19][11] = 0.265829;	
	AAmatrix[19][12] = 6.8532;		AAmatrix[19][13] = 0.723274;	AAmatrix[19][14] = 0.005;		AAmatrix[19][15] = 0.0749218;	
	AAmatrix[19][16] = 0.709226;	AAmatrix[19][17] = 0.005;		AAmatrix[19][18] = 0.0410593;	AAmatrix[19][19] = 0;	

	AAPi[0] = 0.0604902;
	AAPi[1] = 0.0660397;
	AAPi[2] = 0.0441278;
	AAPi[3] = 0.042109;
	AAPi[4] = 0.0200759;
	AAPi[5] = 0.0536065;
	AAPi[6] = 0.0715674;
	AAPi[7] = 0.0723082;
	AAPi[8] = 0.0222939;
	AAPi[9] = 0.0697306;
	AAPi[10] = 0.0988511;
	AAPi[11] = 0.0569682;
	AAPi[12] = 0.0197683;
	AAPi[13] = 0.0288094;
	AAPi[14] = 0.0460253;
	AAPi[15] = 0.0506043;
	AAPi[16] = 0.0536368;
	AAPi[17] = 0.0330116;
	AAPi[18] = 0.0283502;
	AAPi[19] = 0.0616252;

	}
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else
	if(modelnumber==15)
	{

	/*
		HIVw
		
		Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL. (2007) HIV-Specific Probabilistic Models 
		of Protein Evolution. PLoS ONE. Jun 6;2:e503.
	*/
	AAmatrix[0][0] = 0;				AAmatrix[0][1] = 0.0744808;		AAmatrix[0][2] = 0.617509;		AAmatrix[0][3] = 4.43521;	
	AAmatrix[0][4] = 0.167653;		AAmatrix[0][5] = 0.005;			AAmatrix[0][6] = 5.56325;		AAmatrix[0][7] = 1.8685;	
	AAmatrix[0][8] = 0.005;			AAmatrix[0][9] = 0.005;			AAmatrix[0][10] = 0.16024;		AAmatrix[0][11] = 0.592784;	
	AAmatrix[0][12] = 0.005;		AAmatrix[0][13] = 0.597923;		AAmatrix[0][14] = 1.00981;		AAmatrix[0][15] = 8.5942;	
	AAmatrix[0][16] = 24.1422;		AAmatrix[0][17] = 0.005;		AAmatrix[0][18] = 0.005;		AAmatrix[0][19] = 24.8094;	
	AAmatrix[1][0] = 0.0744808;		AAmatrix[1][1] = 0;				AAmatrix[1][2] = 0.16024;		AAmatrix[1][3] = 0.0674539;	
	AAmatrix[1][4] = 2.86364;		AAmatrix[1][5] = 10.6746;		AAmatrix[1][6] = 0.0251632;		AAmatrix[1][7] = 13.4379;	
	AAmatrix[1][8] = 6.84405;		AAmatrix[1][9] = 1.34069;		AAmatrix[1][10] = 0.586757;		AAmatrix[1][11] = 39.8897;	
	AAmatrix[1][12] = 3.28652;		AAmatrix[1][13] = 0.005;		AAmatrix[1][14] = 0.404723;		AAmatrix[1][15] = 8.35024;	
	AAmatrix[1][16] = 0.928203;		AAmatrix[1][17] = 5.96564;		AAmatrix[1][18] = 0.005;		AAmatrix[1][19] = 0.279425;	
	AAmatrix[2][0] = 0.617509;		AAmatrix[2][1] = 0.16024;		AAmatrix[2][2] = 0;				AAmatrix[2][3] = 29.4087;	
	AAmatrix[2][4] = 0.0604932;		AAmatrix[2][5] = 0.342068;		AAmatrix[2][6] = 0.201526;		AAmatrix[2][7] = 0.0604932;	
	AAmatrix[2][8] = 8.59876;		AAmatrix[2][9] = 0.987028;		AAmatrix[2][10] = 0.005;		AAmatrix[2][11] = 10.6655;	
	AAmatrix[2][12] = 0.201526;		AAmatrix[2][13] = 0.005;		AAmatrix[2][14] = 0.344848;		AAmatrix[2][15] = 14.5699;	
	AAmatrix[2][16] = 4.54206;		AAmatrix[2][17] = 0.005;		AAmatrix[2][18] = 5.06475;		AAmatrix[2][19] = 0.0744808;	
	AAmatrix[3][0] = 4.43521;		AAmatrix[3][1] = 0.0674539;		AAmatrix[3][2] = 29.4087;		AAmatrix[3][3] = 0;	
	AAmatrix[3][4] = 0.005;			AAmatrix[3][5] = 0.005;			AAmatrix[3][6] = 12.1233;		AAmatrix[3][7] = 10.3969;	
	AAmatrix[3][8] = 2.31779;		AAmatrix[3][9] = 0.145124;		AAmatrix[3][10] = 0.005;		AAmatrix[3][11] = 0.894313;	
	AAmatrix[3][12] = 0.005;		AAmatrix[3][13] = 0.005;		AAmatrix[3][14] = 0.005;		AAmatrix[3][15] = 0.427881;	
	AAmatrix[3][16] = 0.630395;		AAmatrix[3][17] = 0.005;		AAmatrix[3][18] = 2.28154;		AAmatrix[3][19] = 2.91786;	
	AAmatrix[4][0] = 0.167653;		AAmatrix[4][1] = 2.86364;		AAmatrix[4][2] = 0.0604932;		AAmatrix[4][3] = 0.005;	
	AAmatrix[4][4] = 0;				AAmatrix[4][5] = 0.005;			AAmatrix[4][6] = 0.005;			AAmatrix[4][7] = 0.0489798;	
	AAmatrix[4][8] = 0.005;			AAmatrix[4][9] = 0.005;			AAmatrix[4][10] = 0.005;		AAmatrix[4][11] = 0.005;	
	AAmatrix[4][12] = 0.005;		AAmatrix[4][13] = 0.362959;		AAmatrix[4][14] = 0.005;		AAmatrix[4][15] = 1.12195;	
	AAmatrix[4][16] = 0.005;		AAmatrix[4][17] = 5.49894;		AAmatrix[4][18] = 8.34835;		AAmatrix[4][19] = 0.005;	
	AAmatrix[5][0] = 0.005;			AAmatrix[5][1] = 10.6746;		AAmatrix[5][2] = 0.342068;		AAmatrix[5][3] = 0.005;	
	AAmatrix[5][4] = 0.005;			AAmatrix[5][5] = 0;				AAmatrix[5][6] = 3.20656;		AAmatrix[5][7] = 0.0604932;	
	AAmatrix[5][8] = 18.5465;		AAmatrix[5][9] = 0.0342252;		AAmatrix[5][10] = 2.89048;		AAmatrix[5][11] = 13.0705;	
	AAmatrix[5][12] = 0.005;		AAmatrix[5][13] = 0.005;		AAmatrix[5][14] = 3.04502;		AAmatrix[5][15] = 0.16024;	
	AAmatrix[5][16] = 0.203091;		AAmatrix[5][17] = 0.0443298;	AAmatrix[5][18] = 0.005;		AAmatrix[5][19] = 0.005;	
	AAmatrix[6][0] = 5.56325;		AAmatrix[6][1] = 0.0251632;		AAmatrix[6][2] = 0.201526;		AAmatrix[6][3] = 12.1233;	
	AAmatrix[6][4] = 0.005;			AAmatrix[6][5] = 3.20656;		AAmatrix[6][6] = 0;				AAmatrix[6][7] = 14.7801;	
	AAmatrix[6][8] = 0.005;			AAmatrix[6][9] = 0.0390512;		AAmatrix[6][10] = 0.129839;		AAmatrix[6][11] = 23.9626;	
	AAmatrix[6][12] = 0.005;		AAmatrix[6][13] = 0.005;		AAmatrix[6][14] = 0.005;		AAmatrix[6][15] = 0.005;	
	AAmatrix[6][16] = 0.458743;		AAmatrix[6][17] = 0.005;		AAmatrix[6][18] = 0.005;		AAmatrix[6][19] = 2.19952;	
	AAmatrix[7][0] = 1.8685;		AAmatrix[7][1] = 13.4379;		AAmatrix[7][2] = 0.0604932;		AAmatrix[7][3] = 10.3969;	
	AAmatrix[7][4] = 0.0489798;		AAmatrix[7][5] = 0.0604932;		AAmatrix[7][6] = 14.7801;		AAmatrix[7][7] = 0;	
	AAmatrix[7][8] = 0.005;			AAmatrix[7][9] = 0.005;			AAmatrix[7][10] = 0.0489798;	AAmatrix[7][11] = 0.279425;	
	AAmatrix[7][12] = 0.0489798;	AAmatrix[7][13] = 0.005;		AAmatrix[7][14] = 0.005;		AAmatrix[7][15] = 6.27966;	
	AAmatrix[7][16] = 0.0489798;	AAmatrix[7][17] = 2.8258;		AAmatrix[7][18] = 0.005;		AAmatrix[7][19] = 2.79622;	
	AAmatrix[8][0] = 0.005;			AAmatrix[8][1] = 6.84405;		AAmatrix[8][2] = 8.59876;		AAmatrix[8][3] = 2.31779;	
	AAmatrix[8][4] = 0.005;			AAmatrix[8][5] = 18.5465;		AAmatrix[8][6] = 0.005;			AAmatrix[8][7] = 0.005;	
	AAmatrix[8][8] = 0;				AAmatrix[8][9] = 0.005;			AAmatrix[8][10] = 1.76382;		AAmatrix[8][11] = 0.22406;	
	AAmatrix[8][12] = 0.005;		AAmatrix[8][13] = 0.005;		AAmatrix[8][14] = 13.9444;		AAmatrix[8][15] = 0.725157;	
	AAmatrix[8][16] = 0.95956;		AAmatrix[8][17] = 0.005;		AAmatrix[8][18] = 47.4889;		AAmatrix[8][19] = 0.827479;	
	AAmatrix[9][0] = 0.005;			AAmatrix[9][1] = 1.34069;		AAmatrix[9][2] = 0.987028;		AAmatrix[9][3] = 0.145124;	
	AAmatrix[9][4] = 0.005;			AAmatrix[9][5] = 0.0342252;		AAmatrix[9][6] = 0.0390512;		AAmatrix[9][7] = 0.005;	
	AAmatrix[9][8] = 0.005;			AAmatrix[9][9] = 0;				AAmatrix[9][10] = 9.10246;		AAmatrix[9][11] = 0.817481;	
	AAmatrix[9][12] = 17.3064;		AAmatrix[9][13] = 1.48288;		AAmatrix[9][14] = 0.005;		AAmatrix[9][15] = 0.740091;	
	AAmatrix[9][16] = 9.36345;		AAmatrix[9][17] = 0.005;		AAmatrix[9][18] = 0.114512;		AAmatrix[9][19] = 24.8231;	
	AAmatrix[10][0] = 0.16024;		AAmatrix[10][1] = 0.586757;		AAmatrix[10][2] = 0.005;		AAmatrix[10][3] = 0.005;	
	AAmatrix[10][4] = 0.005;		AAmatrix[10][5] = 2.89048;		AAmatrix[10][6] = 0.129839;		AAmatrix[10][7] = 0.0489798;	
	AAmatrix[10][8] = 1.76382;		AAmatrix[10][9] = 9.10246;		AAmatrix[10][10] = 0;			AAmatrix[10][11] = 0.005;	
	AAmatrix[10][12] = 11.3839;		AAmatrix[10][13] = 7.48781;		AAmatrix[10][14] = 9.83095;		AAmatrix[10][15] = 6.14396;	
	AAmatrix[10][16] = 0.005;		AAmatrix[10][17] = 1.37031;		AAmatrix[10][18] = 0.005;		AAmatrix[10][19] = 2.95344;	
	AAmatrix[11][0] = 0.592784;		AAmatrix[11][1] = 39.8897;		AAmatrix[11][2] = 10.6655;		AAmatrix[11][3] = 0.894313;	
	AAmatrix[11][4] = 0.005;		AAmatrix[11][5] = 13.0705;		AAmatrix[11][6] = 23.9626;		AAmatrix[11][7] = 0.279425;	
	AAmatrix[11][8] = 0.22406;		AAmatrix[11][9] = 0.817481;		AAmatrix[11][10] = 0.005;		AAmatrix[11][11] = 0;	
	AAmatrix[11][12] = 4.09564;		AAmatrix[11][13] = 0.005;		AAmatrix[11][14] = 0.111928;	AAmatrix[11][15] = 0.005;	
	AAmatrix[11][16] = 4.04802;		AAmatrix[11][17] = 0.005;		AAmatrix[11][18] = 0.005;		AAmatrix[11][19] = 0.128065;	
	AAmatrix[12][0] = 0.005;		AAmatrix[12][1] = 3.28652;		AAmatrix[12][2] = 0.201526;		AAmatrix[12][3] = 0.005;	
	AAmatrix[12][4] = 0.005;		AAmatrix[12][5] = 0.005;		AAmatrix[12][6] = 0.005;		AAmatrix[12][7] = 0.0489798;	
	AAmatrix[12][8] = 0.005;		AAmatrix[12][9] = 17.3064;		AAmatrix[12][10] = 11.3839;		AAmatrix[12][11] = 4.09564;	
	AAmatrix[12][12] = 0;			AAmatrix[12][13] = 0.005;		AAmatrix[12][14] = 0.005;		AAmatrix[12][15] = 0.392575;	
	AAmatrix[12][16] = 7.41313;		AAmatrix[12][17] = 0.005;		AAmatrix[12][18] = 0.579198;	AAmatrix[12][19] = 14.7683;	
	AAmatrix[13][0] = 0.597923;		AAmatrix[13][1] = 0.005;		AAmatrix[13][2] = 0.005;		AAmatrix[13][3] = 0.005;	
	AAmatrix[13][4] = 0.362959;		AAmatrix[13][5] = 0.005;		AAmatrix[13][6] = 0.005;		AAmatrix[13][7] = 0.005;	
	AAmatrix[13][8] = 0.005;		AAmatrix[13][9] = 1.48288;		AAmatrix[13][10] = 7.48781;		AAmatrix[13][11] = 0.005;	
	AAmatrix[13][12] = 0.005;		AAmatrix[13][13] = 0;			AAmatrix[13][14] = 0.0342252;	AAmatrix[13][15] = 4.27939;	
	AAmatrix[13][16] = 0.114512;	AAmatrix[13][17] = 0.005;		AAmatrix[13][18] = 4.12728;		AAmatrix[13][19] = 2.28;	
	AAmatrix[14][0] = 1.00981;		AAmatrix[14][1] = 0.404723;		AAmatrix[14][2] = 0.344848;		AAmatrix[14][3] = 0.005;	
	AAmatrix[14][4] = 0.005;		AAmatrix[14][5] = 3.04502;		AAmatrix[14][6] = 0.005;		AAmatrix[14][7] = 0.005;	
	AAmatrix[14][8] = 13.9444;		AAmatrix[14][9] = 0.005;		AAmatrix[14][10] = 9.83095;		AAmatrix[14][11] = 0.111928;	
	AAmatrix[14][12] = 0.005;		AAmatrix[14][13] = 0.0342252;	AAmatrix[14][14] = 0;			AAmatrix[14][15] = 14.249;	
	AAmatrix[14][16] = 4.33701;		AAmatrix[14][17] = 0.005;		AAmatrix[14][18] = 0.005;		AAmatrix[14][19] = 0.005;	
	AAmatrix[15][0] = 8.5942;		AAmatrix[15][1] = 8.35024;		AAmatrix[15][2] = 14.5699;		AAmatrix[15][3] = 0.427881;	
	AAmatrix[15][4] = 1.12195;		AAmatrix[15][5] = 0.16024;		AAmatrix[15][6] = 0.005;		AAmatrix[15][7] = 6.27966;	
	AAmatrix[15][8] = 0.725157;		AAmatrix[15][9] = 0.740091;		AAmatrix[15][10] = 6.14396;		AAmatrix[15][11] = 0.005;	
	AAmatrix[15][12] = 0.392575;	AAmatrix[15][13] = 4.27939;		AAmatrix[15][14] = 14.249;		AAmatrix[15][15] = 0;	
	AAmatrix[15][16] = 6.34079;		AAmatrix[15][17] = 1.10156;		AAmatrix[15][18] = 0.933142;	AAmatrix[15][19] = 0.862637;	
	AAmatrix[16][0] = 24.1422;		AAmatrix[16][1] = 0.928203;		AAmatrix[16][2] = 4.54206;		AAmatrix[16][3] = 0.630395;	
	AAmatrix[16][4] = 0.005;		AAmatrix[16][5] = 0.203091;		AAmatrix[16][6] = 0.458743;		AAmatrix[16][7] = 0.0489798;	
	AAmatrix[16][8] = 0.95956;		AAmatrix[16][9] = 9.36345;		AAmatrix[16][10] = 0.005;		AAmatrix[16][11] = 4.04802;	
	AAmatrix[16][12] = 7.41313;		AAmatrix[16][13] = 0.114512;	AAmatrix[16][14] = 4.33701;		AAmatrix[16][15] = 6.34079;	
	AAmatrix[16][16] = 0;			AAmatrix[16][17] = 0.005;		AAmatrix[16][18] = 0.490608;	AAmatrix[16][19] = 0.005;	
	AAmatrix[17][0] = 0.005;		AAmatrix[17][1] = 5.96564;		AAmatrix[17][2] = 0.005;		AAmatrix[17][3] = 0.005;	
	AAmatrix[17][4] = 5.49894;		AAmatrix[17][5] = 0.0443298;	AAmatrix[17][6] = 0.005;		AAmatrix[17][7] = 2.8258;	
	AAmatrix[17][8] = 0.005;		AAmatrix[17][9] = 0.005;		AAmatrix[17][10] = 1.37031;		AAmatrix[17][11] = 0.005;	
	AAmatrix[17][12] = 0.005;		AAmatrix[17][13] = 0.005;		AAmatrix[17][14] = 0.005;		AAmatrix[17][15] = 1.10156;	
	AAmatrix[17][16] = 0.005;		AAmatrix[17][17] = 0;			AAmatrix[17][18] = 0.005;		AAmatrix[17][19] = 0.005;	
	AAmatrix[18][0] = 0.005;		AAmatrix[18][1] = 0.005;		AAmatrix[18][2] = 5.06475;		AAmatrix[18][3] = 2.28154;	
	AAmatrix[18][4] = 8.34835;		AAmatrix[18][5] = 0.005;		AAmatrix[18][6] = 0.005;		AAmatrix[18][7] = 0.005;	
	AAmatrix[18][8] = 47.4889;		AAmatrix[18][9] = 0.114512;		AAmatrix[18][10] = 0.005;		AAmatrix[18][11] = 0.005;	
	AAmatrix[18][12] = 0.579198;	AAmatrix[18][13] = 4.12728;		AAmatrix[18][14] = 0.005;		AAmatrix[18][15] = 0.933142;	
	AAmatrix[18][16] = 0.490608;	AAmatrix[18][17] = 0.005;		AAmatrix[18][18] = 0;			AAmatrix[18][19] = 1.35482;	
	AAmatrix[19][0] = 24.8094;		AAmatrix[19][1] = 0.279425;		AAmatrix[19][2] = 0.0744808;	AAmatrix[19][3] = 2.91786;	
	AAmatrix[19][4] = 0.005;		AAmatrix[19][5] = 0.005;		AAmatrix[19][6] = 2.19952;		AAmatrix[19][7] = 2.79622;	
	AAmatrix[19][8] = 0.827479;		AAmatrix[19][9] = 24.8231;		AAmatrix[19][10] = 2.95344;		AAmatrix[19][11] = 0.128065;	
	AAmatrix[19][12] = 14.7683;		AAmatrix[19][13] = 2.28;		AAmatrix[19][14] = 0.005;		AAmatrix[19][15] = 0.862637;	
	AAmatrix[19][16] = 0.005;		AAmatrix[19][17] = 0.005;		AAmatrix[19][18] = 1.35482;		AAmatrix[19][19] = 0;	

	AAPi[0] = 0.0377494;
	AAPi[1] = 0.057321;
	AAPi[2] = 0.0891129;
	AAPi[3] = 0.0342034;
	AAPi[4] = 0.0240105;
	AAPi[5] = 0.0437824;
	AAPi[6] = 0.0618606;
	AAPi[7] = 0.0838496;
	AAPi[8] = 0.0156076;
	AAPi[9] = 0.0983641;
	AAPi[10] = 0.0577867;
	AAPi[11] = 0.0641682;
	AAPi[12] = 0.0158419;
	AAPi[13] = 0.0422741;
	AAPi[14] = 0.0458601;
	AAPi[15] = 0.0550846;
	AAPi[16] = 0.0813774;
	AAPi[17] = 0.019597;
	AAPi[18] = 0.0205847;
	AAPi[19] = 0.0515639;


	}

	else
	if(modelnumber==16 )
	{
		// not used any more - see beginning of function

		/*USER DEFINED SUBSTITUTION MODEL*/
	}

	else {cout<<"ERROR IN AMINO ACID MODEL NUMBER IN SETAARATES"<<endl; }


	// now, check that the matrices are symmetrical 
	for (i=0; i<20; i++)
		{
		for (j=i+1; j<20; j++)
			{
			diff = AAmatrix[i][j] - AAmatrix[j][i];
			if (diff < 0.0)
				diff = -diff;
			if (diff > 0.001)
				{
					string errorout="ERROR: ";

				
					if(modelnumber==0)  errorout+="Poisson model ";
					else if(modelnumber==1)  errorout+="Jones model ";
					else if(modelnumber==2)  errorout+="Jones - DCMUT model ";
					else if(modelnumber==3)  errorout+="Dayhoff model ";
					else if(modelnumber==4)  errorout+="Dayhoff - DCMUT model ";
					else if(modelnumber==5)  errorout+="WAG model ";
					else if(modelnumber==6)  errorout+="Mtmam model ";
					else if(modelnumber==7)  errorout+="Mtart model ";
					else if(modelnumber==8)  errorout+="MtREV model ";
					else if(modelnumber==9)  errorout+="rtREV model ";
					else if(modelnumber==10) errorout+="cpREV model ";
					else if(modelnumber==11) errorout+="Vt model ";
					else if(modelnumber==12) errorout+="Blosum model ";
					else if(modelnumber==13) errorout+="LG model ";
					else if(modelnumber==14) errorout+="HIVb model ";
					else if(modelnumber==15) errorout+="HIVw model ";
					//else if(modelnumber==16) errorout+="User-defined PROTEIN model ";
				

					errorout+="is not symmetrical.";

					controlerrorprint2("[MODEL]",name,"getAA(name)",errorout,"");
					error=-1;
				}
			}
		}


	vector<double> row;
	vector<vector<double> > blah;


	bool dobase=false;
	// if no base frequencies have been specified prior to this command
	// then the model base frequencies will come from the empirical
	// substitution model
	if(basefreqs.empty()) dobase=true;

	for(int gj1=0; gj1<20; gj1++)
	{
		if(dobase) basefreqs.push_back(AAPi[gj1]);
		row.clear();
		for(int gj2=0; gj2<20; gj2++) row.push_back(AAmatrix[gj1][gj2]);
		blah.push_back(row);
	}


	// rescale stationary frequencies, to make certain they sum to 1.0 
	double sum = 0.0;
	for (i=0; i<20; i++) sum += basefreqs.at(i); //cout<<basefreqs.at(i)<<" ";} cout<<endl;
	if(sum!=1) for (i=0; i<20; i++) basefreqs.at(i) /= sum;

	// multiply entries by stationary frequencies 
	for (i=0; i<20; i++) for (j=0; j<20; j++) (blah.at(i)).at(j) *= basefreqs.at(j);
			
	// rescale, so branch lengths are in terms of expected number of substitutions per site 
	double scaler = 0.0;
	for (i=0; i<20; i++)
	{
		for (j=i+1; j<20; j++)
		{
			scaler += basefreqs.at(i) * (blah.at(i)).at(j);
			scaler += basefreqs.at(j) * (blah.at(j)).at(i);
		}
	}
	scaler = 1.0 / scaler;
	

	for (i=0; i<20; i++) for (j=0; j<20; j++) (blah.at(i)).at(j) *= scaler;

	//set diagonal of matrix
	for( i=0; i<20; i++) {sum=0; for(j=0; j<20; j++) if(i!=j) sum+=(blah.at(i)).at(j); (blah.at(i)).at(i)=-sum; }

//	for(int p1=0; p1<20; p1++)	for(int p2=0; p2<20; p2++) ccout<<(blah.at(p1)).at(p2)<<"\t"<<AAmatrix[p1][p2]<<endl;

//	for(int p1=0; p1<20; p1++)	for(int p2=0; p2<20; p2++) if(AAmatrix[p1][p2]-(blah.at(p1)).at(p2)!=0) cout<<AAmatrix[p1][p2]-(blah.at(p1)).at(p2)<<"  "<<AAmatrix[p1][p2]<<" "<<(blah.at(p1)).at(p2)<<" "<<endl;

	return blah;
}
////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////
vector<vector<double> > getECMr()
{
	/*
		Kosiol, C., Holmes, I. and Goldman, N. (2007) An Empirical Codon Model for Protein Sequence Evolution.  Molecular Biology and Evolution 24(7): 1464-1479.

		Empirical codon substitution model where evolution is restricted to single nucleotide substitutions
	*/
	vector<vector<double> > myQvec;
	vector<double> row;

	double ECMrestQ[64][64]=
		{{0,11.192024,1.31561,5.427076,1.658051,0,0,0,6.610894,0,0,0,3.347364,0,0,0,2.090751,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.682029,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.261813,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{11.192024,0,0.010896,4.756288,0,1.913571,0,0,0,5.17793,0,0,0,1.558523,0,0,0,2.266373,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.786043,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.923392,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{1.31561,0.010896,0,24.748755,0,0,2.952332,0,0,0,0,0,0,0,0,0,0,0,75.752638,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10.116588,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.36272,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{5.427076,4.756288,24.748755,0,0,0,0,8.126914,0,0,0,0,0,0,0,5.369644,0,0,0,20.877218,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.911096,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.022101,0,0,0,0,0,0,0,0,0,0,0,0}
		,{1.658051,0,0,0,0,13.889102,44.407955,17.057443,2.206054,0,0,0,6.191481,0,0,0,0,0,0,0,1.769355,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,38.2291,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.540052,0,0,0,0,0,0,0,0,0,0,0}
		,{0,1.913571,0,0,13.889102,0,13.681751,65.097021,0,5.615472,0,0,0,9.339206,0,0,0,0,0,0,0,2.704601,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15.793595,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.675838,0,0,0,0,0,0,0,0,0,0}
		,{0,0,2.952332,0,44.407955,13.681751,0,12.991861,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.312811,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.033932,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9.880382,0,0,0,0,0,0,0,0,0}
		,{0,0,0,8.126914,17.057443,65.097021,12.991861,0,0,0,0,0,0,0,0,4.662001,0,0,0,0,0,0,0,1.30348,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,17.103904,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,21.863158,0,0,0,0,0,0,0,0}
		,{6.610894,0,0,0,2.206054,0,0,0,0,19.942818,0,0,0.582084,0,0,0,0,0,0,0,0,0,0,0,3.444964,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.245405,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.367553,0,0,0,0,0,0,0}
		,{0,5.17793,0,0,0,5.615472,0,0,19.942818,0,0,0,0,0.144278,0,0,0,0,0,0,0,0,0,0,0,7.087801,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.228361,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.294702,0,0,0,0,0,0}
		,{3.347364,0,0,0,6.191481,0,0,0,0.582084,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,1.558523,0,0,0,9.339206,0,0,0,0.144278,44.777964,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,5.369644,0,0,0,4.662001,0,0,0.677177,0.073268,0,44.777964,0,0.677177,0,0,0,0,0,0,0,0,0,0,0,0,1.617091,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.719489,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.104727,0,0,0}
		,{2.090751,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.073268,0,0,0,0,0,0,0,0,0,0,0,0,0,1.632945,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.047654,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.552851,0,0}
		,{0,2.266373,0,0,0,0,0,0,0,0,0,0,0,8.905484,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,75.752638,0,0,0,0,0,0,0,0,0,0,56.803876,7.811205,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.026086,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.825082,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.810856}
		,{0,0,0,20.877218,0,0,0,0,0,0,0,0,0,8.432339,22.078564,5.650116,0,8.905484,56.803876,8.432339,1.263838,0,0,0,1.583116,0,0,0,0.779565,0,0,0,4.301225,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.473623,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,1.769355,0,0,0,0,0,0,0,0,1.263838,0,0,0,0,7.811205,22.078564,0,1.389735,0,0,0,3.230751,0,0,0,2.25077,0,0,0,6.381841,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.914972,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,2.704601,0,0,0,0,0,0,0,0,1.389735,0,0,17.461627,0,5.650116,0,0,1.39368,0,0,0,7.419058,0,0,0,3.023939,0,0,0,5.134459,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.737489,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,3.312811,0,0,0,0,0,0,0,0,1.39368,0,35.480963,12.053827,0,0,0,0,0.477616,0,0,0,1.81254,0,0,0,1.462945,0,0,0,2.570123,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.164805,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,1.30348,0,0,0,0,0,0,0,0,0.477616,8.407091,28.557939,11.295213,0,17.461627,35.480963,8.407091,1.021682,0,0,0,0.334165,0,0,0,0,0,0,0,6.578976,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.159185,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,3.444964,0,0,0,0,1.583116,0,0,0,1.021682,0,0,0,0,12.053827,28.557939,0,3.774544,0,0,0,1.699302,0,0,0,0,0,0,0,1.43455,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.120189,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,7.087801,0,0,0,0,3.230751,0,0,0,3.774544,0,0,28.08616,0,11.295213,0,0,5.381868,0,0,0,1.693662,0,0,0,0,0,0,0,0.925575,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.923972,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.419058,0,0,0,5.381868,0,3.44038,1.918904,0,0,0,0,1.794388,0,0,0,3.144296,0,0,0,0,0,0,0,1.23845,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.034856,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.81254,0,0,0,1.794388,1.086327,5.369463,14.959151,0,28.08616,3.44038,1.086327,3.019726,0,0,0,0,0,0,0,0,0,0,0,5.004762,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.383706,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,1.617091,0,0,0.779565,0,0,0,0.334165,0,0,0,3.019726,0,0,0,0,1.918904,5.369463,0,7.016899,0,0,0,0,0,0,0,0,0,0,0,4.105602,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.006827,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,1.632945,0,0,2.25077,0,0,0,1.699302,0,0,0,7.016899,0,0,14.603857,0,14.959151,0,0,6.415757,0,0,0,0,0,0,0,0,0,0,0,6.404436,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.485369,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.023939,0,0,0,1.693662,0,0,0,6.415757,0,99.459951,14.930266,0,0,0,0,19.920977,0,0,0,0,0,0,0,0,0,0,0,2.715692,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.686782,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,3.026086,0,0,0,1.462945,0,0,0,3.144296,0,0,0,19.920977,30.80475,79.48373,13.919752,0,14.603857,99.459951,30.80475,0,0,0,0,0,0,0,0,0,0,0,0,3.834666,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04115,0,0,0}
		,{1.682029,0,0,0,0,0,0,0,0,0,0,0,0,4.301225,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14.930266,79.48373,0,0,0,0,0,0,0,0,0,0,0,0,0,5.03363,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.25247,0,0}
		,{0,0.786043,0,0,0,0,0,0,0,0,0,0,0,0,6.381841,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10.140728,0,13.919752,0,0,0,0,0,0,0,0,0,0,0,0,0,0,92.372238,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.091041,0}
		,{0,0,10.116588,0,0,0,0,0,0,0,0,0,0,0,0,5.134459,0,0,0,0,0,0,0,0,0,0,0,0,0,18.2989,4.623936,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,133.296291,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.803738}
		,{0,0,0,7.911096,0,0,0,0,0,0,0,0,0,0,0,0,2.570123,0,0,0,0,0,0,0,0,0,0,0,0,1.281784,1.303951,2.082128,0,10.140728,18.2989,1.281784,2.801564,0,0,0,0.501054,0,0,0,0.578118,0,0,0,7.096281,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,38.2291,0,0,0,0,0,0,0,0,0,0,0,0,6.578976,0,0,0,0,0,0,0,0,0,0,0,2.801564,0,0,0,0,4.623936,1.303951,0,2.231468,0,0,0,0.292691,0,0,0,0.437779,0,0,0,10.137337,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,15.793595,0,0,0,0,0,0,0,0,0,0,0,0,1.43455,0,0,0,0,0,0,0,0,0,0,0,2.231468,0,0,6.03574,0,2.082128,0,0,4.96235,0,0,0,2.64762,0,0,0,1.903175,0,0,0,25.294298,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,6.033932,0,0,0,0,0,0,0,0,0,0,0,0,0.925575,0,0,0,0,0,0,0,0,0,0,0,4.96235,0,28.307876,6.967655,0,0,0,0,8.155285,0,0,0,0.975074,0,0,0,2.231662,0,0,0,2.078444,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,17.103904,0,0,0,0,0,0,0,0,0,0,0,0,1.23845,0,0,0,0,0,0,0,0,0,0,0,8.155285,19.578982,38.414969,12.678802,0,6.03574,28.307876,19.578982,11.715476,0,0,0,39.399322,0,0,0,0,0,0,0,5.107629,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,2.245405,0,0,0,0,0,0,0,0,0,0,0,0,5.004762,0,0,0,0,0,0,0,0.501054,0,0,0,11.715476,0,0,0,0,6.967655,38.414969,0,2.13474,0,0,0,21.337943,0,0,0,0,0,0,0,2.312255,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0.228361,0,0,0,0,0,0,0,0,0,0,0,0,4.105602,0,0,0,0,0,0,0,0.292691,0,0,0,2.13474,0,0,13.863648,0,12.678802,0,0,3.91936,0,0,0,0.754055,0,0,0,0,0,0,0,3.064069,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.404436,0,0,0,0,0,0,0,2.64762,0,0,0,3.91936,0,4.929483,0.366267,0,0,0,0,5.869857,0,0,0,22.577271,0,0,0,0,0,0,0,25.461549,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.715692,0,0,0,0,0,0,0,0.975074,0,0,0,5.869857,1.010212,0.982893,10.762877,0,13.863648,4.929483,1.010212,16.623529,0,0,0,0,0,0,0,0,0,0,0,6.091654,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,4.719489,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.834666,0,0,0,0.578118,0,0,0,39.399322,0,0,0,16.623529,0,0,0,0,0.366267,0.982893,0,7.784768,0,0,0,0,0,0,0,0,0,0,0,3.686074,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,2.047654,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.03363,0,0,0,0.437779,0,0,0,21.337943,0,0,0,7.784768,0,0,26.637668,0,10.762877,0,0,8.423762,0,0,0,0,0,0,0,0,0,0,0,2.811398,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,92.372238,0,0,0,1.903175,0,0,0,0.754055,0,0,0,8.423762,0,1.792245,0.1209,0,0,0,0,21.000358,0,0,0,0,0,0,0,0,0,0,0,2.090641,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0.825082,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,133.296291,0,0,0,2.231662,0,0,0,22.577271,0,0,0,21.000358,3.324581,6.01197,36.292705,0,26.637668,1.792245,3.324581,0,0,0,0,0,0,0,0,0,0,0,0,10.59078,0,0,0}
		,{2.261813,0,0,0,0,0,0,0,0,0,0,0,0,2.473623,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.096281,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1209,6.01197,0,0,0,0,0,0,0,0,0,0,0,0,0,4.285543,0,0}
		,{0,1.923392,0,0,0,0,0,0,0,0,0,0,0,0,5.914972,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10.137337,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.669955,0,36.292705,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.43241,0}
		,{0,0,2.36272,0,0,0,0,0,0,0,0,0,0,0,0,3.737489,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,25.294298,0,0,0,0,0,0,0,0,0,0,0,0,0,26.045078,3.531461,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.388514}
		,{0,0,0,2.022101,0,0,0,0,0,0,0,0,0,0,0,0,2.164805,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.078444,0,0,0,0,0,0,0,0,0,0,0,0,8.901167,21.657664,11.898141,0,6.669955,26.045078,8.901167,3.682092,0,0,0,0.352915,0,0,0,0.503385,0,0,0}
		,{0,0,0,0,5.540052,0,0,0,0,0,0,0,0,0,0,0,0,1.159185,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.107629,0,0,0,0,0,0,0,0,0,0,0,3.682092,0,0,0,0,3.531461,21.657664,0,4.308415,0,0,0,0.208522,0,0,0,0.542717,0,0}
		,{0,0,0,0,0,7.675838,0,0,0,0,0,0,0,0,0,0,0,0,3.120189,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.312255,0,0,0,0,0,0,0,0,0,0,0,4.308415,0,0,6.516319,0,11.898141,0,0,6.291148,0,0,0,1.277861,0,0,0,0.702411,0}
		,{0,0,0,0,0,0,9.880382,0,0,0,0,0,0,0,0,0,0,0,0,2.923972,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.064069,0,0,0,0,0,0,0,0,0,0,0,6.291148,0,21.910225,5.090423,0,0,0,0,6.166554,0,0,0,0.476105,0,0,0,0.302501}
		,{0,0,0,0,0,0,0,21.863158,0,0,0,0,0,0,0,0,0,0,0,0,6.034856,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,25.461549,0,0,0,0,0,0,0,0,0,0,0,6.166554,5.512586,20.715347,9.529141,0,6.516319,21.910225,5.512586,0.693026,0,0,0,1.541379,0,0,0}
		,{0,0,0,0,0,0,0,0,0.367553,0,0,0,0,0,0,0,0,0,0,0,0,0.383706,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.091654,0,0,0,0,0,0,0,0.352915,0,0,0,0.693026,0,0,0,0,5.090423,20.715347,0,1.866565,0,0,0,2.303487,0,0}
		,{0,0,0,0,0,0,0,0,0,0.294702,0,0,0,0,0,0,0,0,0,0,0,0,3.006827,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.686074,0,0,0,0,0,0,0,0.208522,0,0,0,1.866565,0,0,10.605899,0,9.529141,0,0,2.774445,0,0,0,2.985093,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.485369,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.811398,0,0,0,0,0,0,0,1.277861,0,0,0,2.774445,0,2.71061,0.650088,0,0,0,0,9.441919,0,0,0,6.644971}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.686782,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2.090641,0,0,0,0,0,0,0,0.476105,0,0,0,9.441919,1.296294,3.779053,10.15357,0,10.605899,2.71061,1.296294,1.042624,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,1.104727,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.04115,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10.59078,0,0,0,0.503385,0,0,0,1.541379,0,0,0,1.042624,0,0,0,0,0.650088,3.779053,0,1.561629,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0.552851,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.25247,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.285543,0,0,0,0.542717,0,0,0,2.303487,0,0,0,1.561629,0,0,9.48852,0,10.15357,0,0,0.874,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.091041,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.43241,0,0,0,0.702411,0,0,0,2.985093,0,0,0,0.874,0,20.5181,4.120953,0,0,0,0,1.39381}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0.810856,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.803738,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5.388514,0,0,0,0.302501,0,0,0,6.644971,0,0,0,1.39381,13.246936,18.064826,19.084271,0,9.48852,20.5181,13.246936}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.120953,18.064826}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,19.084271}
		,{0.022103,0,0.021383,0,0.016387,0,0.015425,0,0.01188,0,0.011131,0,0.00975,0,0.008956,0,0.015965,0,0.015782,0,0.006025,0,0.007029,0,0.01188,0,0.014467,0,0.017386,0,0.0076,0,0.028839,0,0.010007,0,0.0101,0,0.010642,0,0.011843,0,0.011097,0,0.011703,0,0.016076,0,0.020211,0,0.008311,0,0.014148,0,0.0048,0,0.007837,0,0.025576,0,0.023441,0,0.013551,0}};


		double ECMrestPi[64]=
		{0.022103,0.021383,0.016387,0.015425,0.01188,0.011131,0.00975,0.008956,0.015965,0.015782,0,0,0.006025,0.007029,0,0.01188,0.014467,0.017386,0.0076,0.028839,0.010007,0.0101,0.010642,0.011843,0.011097,0.011703,0.016076,0.020211,0.008311,0.014148,0.0048,0.007837,0.025576,0.023441,0.013551,0.020102,0.013424,0.020201,0.015528,0.012142,0.023006,0.020171,0.030001,0.026344,0.010142,0.011679,0.010372,0.008195,0.019047,0.018938,0.010901,0.022747,0.019005,0.028307,0.015908,0.018853,0.028198,0.024532,0.033223,0.031878,0.016852,0.022982,0.015796,0.010191};


		if(basefreqs.empty()) for(int gv3=0; gv3<64; gv3++) basefreqs.push_back(ECMrestPi[gv3]); 
			
		for(int gv1=0; gv1<64; gv1++) {row.clear(); {for(int gv2=0; gv2<64; gv2++) row.push_back(ECMrestQ[gv1][gv2]);} myQvec.push_back(row);}

		return myQvec;
}

vector<vector<double> > getECMu()
{
	
	/*
		Kosiol, C., Holmes, I. and Goldman, N. (2007) An Empirical Codon Model for Protein Sequence Evolution.  Molecular Biology and Evolution 24(7): 1464-1479.

		Empirical codon substitution model where doublet and triplet substitutions are allowed as well 
	*/
	vector<vector<double> > myQvec;
	vector<double> row;
		double ECMunrestQ[64][64]=
		{{0,16.011531,2.395822,1.204356,0.773935,0.030074,0.27809,0.034137,4.317981,0.481042,0,0,0.733587,0.045951,0,0.786871,2.016257,0.083684,1.036474,0.07355,0.324368,0.00014,0.025217,0.004063,0.638696,0.000467,0.12632,0.00976,0.143832,0,0.030121,0.000934,1.119411,0.056038,1.075187,0.67937,0.497293,0.000093,0.079948,0.008032,0.263567,0.000093,0.148268,0.013122,0.225554,0.000093,0.140516,0.08331,0.86397,0.026338,0.565566,0.068927,0.483563,0.000093,0.109975,0.005884,0.064397,0,0.084945,0.010974,0.164659,0,0.113991,0.018773}
		,{16.011531,0,0.151858,0.675537,0.052602,0.656004,0.056677,0.198277,0.503397,4.483501,0,0,0.076912,0.56162,0,1.183337,0.207692,2.30611,0.198558,1.341144,0.000141,0.285635,0.006558,0.079161,0.001312,0.761771,0.016628,0.107218,0.000094,0.200806,0.01602,0.027355,0.059956,1.006045,0.064968,0.800602,0.000141,0.371541,0.010539,0.036489,0.000094,0.143614,0.005996,0.043609,0.000047,0.166706,0.0126,0.056771,0.065905,0.946136,0.047403,0.371072,0.000234,0.512203,0.016535,0.074851,0,0.055366,0.006464,0.034428,0.000141,0.172889,0.018315,0.032039}
		,{2.395822,0.151858,0,18.541946,0.249707,0.011609,1.184813,0.010188,0.798582,0.033529,0,0,0.645571,0.040012,0,0.271072,12.035723,1.373823,27.219895,1.045943,0.001358,0,0.261069,0,0.026551,0.000123,0.576476,0,0.000741,0,0.136647,0,2.130663,0.042112,5.159075,1.418466,0.012473,0.000062,0.871568,0.000062,0.271628,0.00284,0.61266,0.000062,0.010929,0.000432,0.423774,0.000247,0.748196,0.00599,2.299543,0.024699,0.008892,0,1.041312,0,0.042112,0.000062,0.287373,0,0.000741,0,0.201112,0}
		,{1.204356,0.675537,18.541946,0,0.27499,0.158873,0.611887,0.694091,0.337279,0.177833,0,0,0.395942,0.240632,0,0.632947,11.161511,5.651603,16.560966,12.455337,0.003499,0.000382,0.005535,0.112999,0.040275,0.002163,0.007508,0.250748,0.003054,0.000064,0.001527,0.127696,1.292935,0.478019,1.065537,3.062807,0.015652,0.00299,0.011134,0.586818,0.077878,0.060699,0.004963,0.33378,0.009289,0.00579,0.001718,0.506841,0.529619,0.140867,0.762425,1.108802,0.03563,0.001654,0.019406,0.220908,0.038557,0.006808,0.005472,0.159955,0.003881,0.000191,0.001082,0.175861}
		,{0.773935,0.052602,0.249707,0.27499,0,23.65509,35.921779,11.510965,0.688169,0.069588,0,0,1.811753,0.138244,0,0.069758,0.277929,0.000085,0.000678,0,2.846677,0.101204,0.487542,0.021444,1.253945,0.000593,0.508308,0.049246,0.660622,0,0.006103,0.000085,0.172403,0,0.000424,0.093491,4.693944,0.196983,2.018483,0.361164,1.773102,0.00339,0.340991,0.156468,12.169045,0.094762,0.04789,0.063994,0.563995,0.000085,0.001356,0.000254,3.994417,0.111968,1.93135,0.103323,1.120532,0.000254,0.330481,0.04238,0.976185,0.000085,0.012121,0.002797}
		,{0.030074,0.656004,0.011609,0.158873,23.65509,0,15.982573,35.359077,0.047115,0.524116,0,0,0.343463,1.323765,0,0.081312,0.000186,0.342813,0.000186,0.001022,0.196358,2.487136,0.138742,0.371063,0.002137,1.144692,0.080383,0.423382,0.001208,0.685812,0.004089,0.004832,0,0.115975,0,0.042282,0.487317,3.82958,0.409721,1.562591,0.000929,1.118208,0.020909,0.25165,0.083636,10.892976,0.002044,0.028715,0.000186,0.396897,0.000279,0.001394,0.771771,3.45257,0.5585,1.262618,0.003717,1.023142,0.08568,0.283432,0.001951,0.880032,0.001951,0.002974}
		,{0.27809,0.056677,1.184813,0.611887,35.921779,15.982573,0,17.424222,0.341791,0.070809,0,0,0.75198,0.121937,0,0.195833,0.000289,0.000096,0.496046,0,0.544474,0.072352,3.121656,0.064924,0.128111,0.01447,2.066955,0.002122,0.000579,0,0.557015,0,0.000386,0.000096,0.403435,0.246094,2.297807,0.150107,5.709933,0.516594,0.872084,0.187826,1.496628,0.004438,5.964323,1.049877,1.094736,0.009261,0.002219,0.000096,0.81372,0.000772,1.825878,0.084218,4.380679,0.150589,0.348448,0.007428,1.265487,0.001061,0.011673,0.000289,1.720919,0.003376}
		,{0.034137,0.198277,0.010188,0.694091,11.510965,35.359077,17.424222,0,0.058136,0.213967,0,0,0.143447,0.493179,0,0.410046,0,0.000344,0.000115,0.266943,0.078776,0.43252,0.151589,2.075226,0.07373,0.114551,0.002179,1.519211,0.00172,0.032106,0.003211,1.903571,0,0.000344,0,0.527005,0.199748,2.395833,0.349846,4.919174,0.040706,0.725836,0.000459,0.768607,0.681575,9.818281,0.000115,0.514392,0.000115,0.001261,0.000115,0.451899,0.250545,1.978448,0.505677,4.658653,0.117533,0.670108,0.002179,1.029128,0.007109,0.356038,0.00172,2.163175}
		,{4.317981,0.503397,0.798582,0.337279,0.688169,0.047115,0.341791,0.058136,0,24.177765,0,0,0.822999,0.068342,0,1.140051,0.485469,0.000116,0.01665,0.004815,0.337879,0.000116,0.03214,0.004177,3.088481,0.265766,0.486281,0.09207,0.534375,0.000116,0.043917,0.003713,0.352731,0.000116,0.573013,0.294368,0.599932,0.000116,0.208041,0.042119,0.74787,0.046586,0.467136,0.072867,0.50647,0.000116,0.424321,0.200499,0.571505,0.000116,0.236179,0.00963,0.370309,0.000058,0.176829,0.027035,0.223763,0.010037,0.257122,0.042815,0.13094,0.000058,0.082323,0.007948}
		,{0.481042,4.483501,0.033529,0.177833,0.069588,0.524116,0.070809,0.213967,24.177765,0,0,0,0.05486,0.628438,0,1.421996,0.000299,0.622089,0.011978,0.308859,0.000479,0.310416,0.025873,0.037372,0.340541,3.193996,0.079236,0.332396,0.001377,0.604541,0.051686,0.081931,0.00018,0.255975,0.025454,0.300354,0.000599,0.393545,0.053363,0.177757,0.042762,0.487814,0.042942,0.140864,0.00018,0.34689,0.060909,0.26951,0.000359,0.582801,0.06043,0.116309,0.000419,0.380788,0.034737,0.106187,0.015452,0.184704,0.043721,0.136432,0.00012,0.127388,0.029826,0.014314}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0.733587,0.076912,0.645571,0.395942,1.811753,0.343463,0.75198,0.143447,0.822999,0.05486,0,0,0,56.838378,0,0.264556,0.54324,0.000466,0.020649,0.002639,0.239715,0,0.002795,0.000155,0.634065,0.000155,0.163174,0.05791,0.726908,0.012886,0.232728,0.023909,0.431456,0.000311,0.069555,0.298091,1.089585,0.000776,0.415775,0.053874,0.360038,0.000932,0.099519,0.011489,2.007768,0.248099,0.144388,0.122186,1.598824,0.001087,0.754233,0.126844,1.899865,0.000466,0.806554,0.028567,0.099985,0,0.028878,0.014439,0.420899,0.007608,0.197641,0.105884}
		,{0.045951,0.56162,0.040012,0.240632,0.138244,1.323765,0.121937,0.493179,0.068342,0.628438,0,0,56.838378,0,0,0.210115,0.000674,0.674176,0.021578,0.265948,0.00027,0.215779,0.01025,0.004585,0.001483,0.483076,0.032232,0.105597,0.077815,0.516927,0.16615,0.183143,0.000405,0.309643,0.012138,0.324613,0.001483,0.806071,0.061227,0.173298,0.000135,0.325422,0.004316,0.034794,0.181794,1.372357,0.030883,0.070533,0.004316,1.273908,0.086986,0.530278,0.011733,1.514097,0.297371,0.586111,0.000135,0.071612,0.003641,0.013216,0.045044,0.309374,0.061497,0.183952}
		,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
		,{0.786871,1.183337,0.271072,0.632947,0.069758,0.081312,0.195833,0.410046,1.140051,1.421996,0,0,0.264556,0.210115,0,0,0.010122,0.113701,0.017106,0.504866,0.061833,0.032564,0.070308,0.215788,0.195073,0.369273,0.055163,0.24749,0.019696,0.176476,0.146501,1.13591,0.078312,0.136849,0.170041,0.321642,0.035939,0.037822,0.060421,0.362681,0.074859,0.053908,0.051005,0.105226,0.139046,0.167138,0.145245,0.496706,0.038763,0.092044,0.058616,0.277072,0.03586,0.044178,0.031466,0.446015,0.028249,0.066384,0.009966,0.137634,0.039313,0.105305,0.073682,0.381671}
		,{2.016257,0.207692,12.035723,11.161511,0.277929,0.000186,0.000289,0,0.485469,0.000299,0,0,0.54324,0.000674,0,0.010122,0,15.874441,21.437257,4.802017,0.822643,0.026571,0.065669,0.007978,0.664866,0.058614,0.529045,0.079119,0.663877,0.016022,0.424607,0.039428,3.330793,0.39019,1.260239,1.22002,0.831215,0.000396,0.000857,0.000264,0.25938,0.00033,0.060922,0.039362,0.322807,0.00033,0.004747,0.00956,1.89715,0.105361,0.509595,0.091844,0.879741,0.000132,0.002374,0.000066,0.129492,0.000066,0.03956,0.00422,0.169777,0,0.00033,0.000066}
		,{0.083684,2.30611,1.373823,5.651603,0.000085,0.342813,0.000096,0.000344,0.000116,0.622089,0,0,0.000466,0.674176,0,0.113701,15.874441,0,8.808275,15.484088,0.036254,0.648769,0.016609,0.118229,0.057328,0.617694,0.071794,0.422671,0.068758,0.544828,0.112395,0.640495,0.18401,3.765697,0.136148,1.434579,0.00006,0.89749,0.000119,0.000595,0,0.18788,0.002143,0.022384,0.00006,0.300513,0.00006,0.000952,0.072866,2.720153,0.03173,0.477439,0.000417,1.191514,0.000357,0.000893,0,0.135255,0.002679,0.010061,0.00006,0.240505,0.00006,0.000119}
		,{1.036474,0.198558,27.219895,16.560966,0.000678,0.000186,0.496046,0.000115,0.01665,0.011978,0,0,0.020649,0.021578,0,0.017106,21.437257,8.808275,0,8.319767,0.181411,0.040087,1.07379,0.016442,0.438648,0.059927,1.205738,0.105449,0.134394,0.005436,0.918198,0.040902,1.328581,0.203017,4.40061,1.635281,0.004348,0.000679,0.94456,0.000408,0.019568,0.01454,0.43362,0.000679,0.009512,0.002718,0.409432,0.000951,0.555104,0.043756,2.047159,0.173122,0.004077,0.000136,0.741135,0,0.012366,0.001359,0.313495,0.000136,0.000272,0,0.165784,0}
		,{0.07355,1.341144,1.045943,12.455337,0,0.001022,0,0.266943,0.004815,0.308859,0,0,0.002639,0.265948,0,0.504866,4.802017,15.484088,8.319767,0,0.014388,0.149771,0.040917,0.495176,0.044742,0.330036,0.033372,0.703795,0.015019,0.563956,0.041969,0.794366,0.089308,2.469249,0.048882,2.236557,0.00007,0.073657,0.00007,0.733202,0.00014,0.034916,0.000035,0.133032,0.000105,0.0173,0,0.04032,0.01158,1.08517,0.020529,2.26249,0.000772,0.103766,0.000351,1.524024,0.000491,0.015686,0.00014,0.1763,0.000175,0.047268,0.00007,0.185038}
		,{0.324368,0.000141,0.001358,0.003499,2.846677,0.196358,0.544474,0.078776,0.337879,0.000479,0,0,0.239715,0.00027,0,0.061833,0.822643,0.036254,0.181411,0.014388,0,23.496083,40.922701,10.291826,0.775254,0.061583,0.435109,0.107434,0.500433,0.002398,0.352807,0.00988,0.292855,0.000096,0.002014,0.081631,1.050363,0.044125,0.368538,0.079712,0.78734,0.000959,0.247195,0.220816,1.035015,0.001247,0.011511,0.060432,0.708395,0.000192,0.001151,0.000384,1.272426,0.022062,0.335924,0.014101,0.661776,0.000096,0.184749,0.034437,0.418802,0.000096,0.003549,0.001918}
		,{0.00014,0.285635,0,0.000382,0.101204,2.487136,0.072352,0.43252,0.000116,0.310416,0,0,0,0.215779,0,0.032564,0.026571,0.648769,0.040087,0.149771,23.496083,0,15.426733,33.45378,0.091276,0.730306,0.074846,0.529594,0.124232,0.563799,0.154017,0.897101,0.000096,0.270274,0,0.008455,0.053805,0.870967,0.02623,0.435243,0.000192,0.532092,0.012394,0.303613,0.000288,1.337629,0.000384,0.033244,0.000192,0.76009,0.000192,0.001537,0.073021,1.394892,0.069754,0.417565,0.000769,0.976175,0.105112,0.294294,0.000288,0.636916,0.000384,0.001441}
		,{0.025217,0.006558,0.261069,0.005535,0.487542,0.138742,3.121656,0.151589,0.03214,0.025873,0,0,0.002795,0.01025,0,0.070308,0.065669,0.016609,1.07379,0.040917,40.922701,15.426733,0,15.127582,0.286252,0.089835,1.05204,0.184327,0.063413,0.001702,0.626603,0.0103,0.002597,0.000448,0.480521,0.042006,0.345545,0.027138,1.018005,0.083475,0.096104,0.074161,1.042457,0.012002,0.021048,0.005195,0.493508,0.009225,0.009225,0.000537,0.732111,0.000717,0.346978,0.013077,1.236457,0.017824,0.147873,0.003672,0.890822,0.001791,0.002508,0.00009,0.556204,0.001254}
		,{0.004063,0.079161,0,0.112999,0.021444,0.371063,0.064924,2.075226,0.004177,0.037372,0,0,0.000155,0.004585,0,0.215788,0.007978,0.118229,0.016442,0.495176,10.291826,33.45378,15.127582,0,0.054021,0.364129,0.063366,0.715716,0.044676,0.798346,0.091073,1.164525,0.000246,0.021723,0,0.193459,0.011476,0.527094,0.015001,0.962295,0.002705,0.095418,0.00041,0.561523,0.001066,0.104681,0,0.273876,0.000246,0.031642,0.000082,0.602428,0.013526,0.971148,0.087384,1.950083,0.03156,0.644235,0.00541,0.99033,0.001312,0.395771,0.000164,0.703092}
		,{0.638696,0.001312,0.026551,0.040275,1.253945,0.002137,0.128111,0.07373,3.088481,0.340541,0,0,0.634065,0.001483,0,0.195073,0.664866,0.057328,0.438648,0.044742,0.775254,0.091276,0.286252,0.054021,0,38.685701,2.473439,1.106179,2.460976,0.170088,1.35386,0.316372,0.193328,0.000469,0.040109,0.323588,0.898794,0.001125,0.082373,0.076938,2.691226,0.46097,0.899262,0.324525,1.293228,0.00506,0.411115,0.140943,0.338207,0.000375,0.019586,0.032237,0.42058,0.000562,0.039265,0.080124,0.746792,0.100928,0.452442,0.159217,0.388156,0.000843,0.097554,0.08406}
		,{0.000467,0.761771,0.000123,0.002163,0.000593,1.144692,0.01447,0.114551,0.265766,3.193996,0,0,0.000155,0.483076,0,0.369273,0.058614,0.617694,0.059927,0.330036,0.061583,0.730306,0.089835,0.364129,38.685701,0,0.751904,2.503268,0.277265,2.478358,0.526904,2.20843,0,0.127899,0.000272,0.163406,0.000272,0.969387,0.02192,0.10308,0.188587,2.203539,0.143479,0.571469,0.000453,1.282522,0.025544,0.169022,0.000091,0.491034,0.002808,0.044112,0.000091,0.696016,0.004982,0.190037,0.046739,0.975727,0.106069,0.566034,0.000091,0.566759,0.004982,0.053714}
		,{0.12632,0.016628,0.576476,0.007508,0.508308,0.080383,2.066955,0.002179,0.486281,0.079236,0,0,0.163174,0.032232,0,0.055163,0.529045,0.071794,1.205738,0.033372,0.435109,0.074846,1.05204,0.063366,2.473439,0.751904,0,17.923045,1.164262,0.14894,4.72584,0.299978,0.078926,0.010543,0.390087,0.443617,0.374419,0.066519,1.660885,0.002854,1.759732,0.377447,3.215817,0.461383,1.177718,0.232701,1.24214,0.003786,0.150804,0.006757,0.606945,0.000466,0.27295,0.018814,1.274118,0.001165,0.706782,0.121389,3.081614,0.343314,0.042812,0.016193,0.551493,0.003029}
		,{0.00976,0.107218,0,0.250748,0.049246,0.423382,0.002122,1.519211,0.09207,0.332396,0,0,0.05791,0.105597,0,0.24749,0.079119,0.422671,0.105449,0.703795,0.107434,0.529594,0.184327,0.715716,1.106179,2.503268,17.923045,0,0.340811,2.029914,0.61732,4.718199,0.003859,0.105885,0.000048,0.834976,0.029088,0.617464,0.001302,1.766961,0.206851,0.985145,0.285384,2.285052,0.083261,1.418383,0.000289,1.148387,0.0096,0.06932,0.000145,0.475496,0.040907,0.634927,0.002219,1.544626,0.130873,0.928319,0.536567,3.036767,0.003377,0.336277,0.000289,0.634203}
		,{0.143832,0.000094,0.000741,0.003054,0.660622,0.001208,0.000579,0.00172,0.534375,0.001377,0,0,0.726908,0.077815,0,0.019696,0.663877,0.068758,0.134394,0.015019,0.500433,0.124232,0.063413,0.044676,2.460976,0.277265,1.164262,0.340811,0,27.244097,39.595443,12.868484,0.076434,0.000118,0,0.028736,0.344601,0.00106,0.000589,0.004593,0.682254,0.00318,1.769879,0.560831,0.772349,0.387941,17.450524,6.798629,0.336828,0.000589,0.000353,0.001413,0.324227,0.000236,0.000589,0.001531,0.162525,0.000236,0.034978,0.007891,0.241197,0.021435,0.01531,0.043222}
		,{0,0.200806,0,0.000064,0,0.685812,0,0.032106,0.000116,0.604541,0,0,0.012886,0.516927,0,0.176476,0.016022,0.544828,0.005436,0.563956,0.002398,0.563799,0.001702,0.798346,0.170088,2.478358,0.14894,2.029914,27.244097,0,12.677657,35.563093,0,0.238839,0,0.029786,0,1.481721,0,0.014243,0.000068,0.863191,0.296358,2.721043,0.018146,1.320875,1.113671,4.087042,0,0.648523,0,0.003287,0,1.610248,0.000068,0.083744,0,0.915505,0.025678,0.528692,0.004656,0.676049,0.000753,0.097165}
		,{0.030121,0.01602,0.136647,0.001527,0.006103,0.004089,0.557015,0.003211,0.043917,0.051686,0,0,0.232728,0.16615,0,0.146501,0.424607,0.112395,0.918198,0.041969,0.352807,0.154017,0.626603,0.091073,1.35386,0.526904,4.72584,0.61732,39.595443,12.677657,0,30.574631,0.000416,0.001248,0.121855,0.015596,0.00104,0.002079,0.360575,0.002911,0.009981,0.016636,3.06551,0.034519,0.530881,0.354545,31.949764,15.287419,0.003743,0.001871,0.255147,0.000832,0.001871,0.000416,0.286547,0.000624,0.00707,0.009981,0.440217,0.00104,0.042005,0.008942,0.247245,0.143481}
		,{0.000934,0.027355,0,0.127696,0.000085,0.004832,0,1.903571,0.003713,0.081931,0,0,0.023909,0.183143,0,1.13591,0.039428,0.640495,0.040902,0.794366,0.00988,0.897101,0.0103,1.164525,0.316372,2.20843,0.299978,4.718199,12.868484,35.563093,30.574631,0,0,0.003064,0,0.40868,0,0.017774,0,2.985421,0.000981,0.065212,0.011032,3.8322,0.025006,1.360141,2.418859,18.531553,0,0.005026,0,0.701399,0,0.020471,0.000123,3.409178,0.000368,0.150527,0.000858,2.171984,0.011768,0.703728,0.010419,0.590833}
		,{1.119411,0.059956,2.130663,1.292935,0.172403,0,0.000386,0,0.352731,0.00018,0,0,0.431456,0.000405,0,0.078312,3.330793,0.18401,1.328581,0.089308,0.292855,0.000096,0.002597,0.000246,0.193328,0,0.078926,0.003859,0.076434,0,0.000416,0,0,13.60931,16.415611,1.155098,1.266654,0.034573,0.296014,0.090674,0.239388,0.025405,0.221536,0.04144,0.359183,0.03114,0.116039,0.093542,6.546163,0.458622,2.581654,0.953353,0.536332,0.00004,0.002262,0.000081,0.066966,0,0.038612,0.003312,0.069995,0.000283,0.000283,0.000081}
		,{0.056038,1.006045,0.042112,0.478019,0,0.115975,0.000096,0.000344,0.000116,0.255975,0,0,0.000311,0.309643,0,0.136849,0.39019,3.765697,0.203017,2.469249,0.000096,0.270274,0.000448,0.021723,0.000469,0.127899,0.010543,0.105885,0.000118,0.238839,0.001248,0.003064,13.60931,0,5.784672,1.428293,0.075878,1.066285,0.048902,0.311759,0.014351,0.175491,0.02302,0.116405,0.018727,0.238533,0.012331,0.062369,0.575921,7.487816,0.313779,3.487342,0.000253,0.536362,0.000589,0.004629,0.000042,0.032447,0.009174,0.031984,0,0.055425,0.000084,0.000295}
		,{1.075187,0.064968,5.159075,1.065537,0.000424,0,0.403435,0,0.573013,0.025454,0,0,0.069555,0.012138,0,0.170041,1.260239,0.136148,4.40061,0.048882,0.002014,0,0.480521,0,0.040109,0.000272,0.390087,0.000048,0,0,0.121855,0,16.415611,5.784672,0,2.230691,0.351882,0.016701,2.260424,0.154441,0.256283,0.020837,0.667397,0.056267,0.269862,0.021539,0.780008,0.317466,2.577578,0.154207,11.271062,0.911193,0.001795,0,0.983926,0.000078,0.001171,0,0.361403,0.000078,0.000156,0,0.194319,0.000078}
		,{0.67937,0.800602,1.418466,3.062807,0.093491,0.042282,0.246094,0.527005,0.294368,0.300354,0,0,0.298091,0.324613,0,0.321642,1.22002,1.434579,1.635281,2.236557,0.081631,0.008455,0.042006,0.193459,0.323588,0.163406,0.443617,0.834976,0.028736,0.029786,0.015596,0.40868,1.155098,1.428293,2.230691,0,0.419831,0.433759,0.853779,1.376727,0.208924,0.21917,0.275355,0.497593,0.306375,0.304581,0.305714,0.90581,0.430124,0.350473,0.569076,1.399673,0.451134,0.169264,0.517896,0.837302,0.059065,0.011379,0.033994,0.262465,0.027479,0.018603,0.037724,0.410199}
		,{0.497293,0.000141,0.012473,0.015652,4.693944,0.487317,2.297807,0.199748,0.599932,0.000599,0,0,1.089585,0.001483,0,0.035939,0.831215,0.00006,0.004348,0.00007,1.050363,0.053805,0.345545,0.011476,0.898794,0.000272,0.374419,0.029088,0.344601,0,0.00104,0,1.266654,0.075878,0.351882,0.419831,0,13.991583,31.915858,12.116657,2.057449,0.296066,0.878163,0.291009,4.943439,0.622868,0.432918,0.309656,2.898791,0.001027,0.009008,0.003635,2.359362,0.038559,0.381796,0.023862,0.928969,0.000158,0.251423,0.033581,0.380374,0,0.002449,0.000553}
		,{0.000093,0.371541,0.000062,0.00299,0.196983,3.82958,0.150107,2.395833,0.000116,0.393545,0,0,0.000776,0.806071,0,0.037822,0.000396,0.89749,0.000679,0.073657,0.044125,0.870967,0.027138,0.527094,0.001125,0.969387,0.066519,0.617464,0.00106,1.481721,0.002079,0.017774,0.034573,1.066285,0.016701,0.433759,13.991583,0,8.373639,28.470047,0.067554,1.346385,0.089476,0.623366,0.275865,4.699375,0.021698,0.140701,0.00302,2.862216,0.001957,0.089643,0.120121,2.159328,0.077341,0.728891,0.000559,1.013424,0.109664,0.360196,0.000112,0.518903,0.000112,0.000447}
		,{0.079948,0.010539,0.871568,0.011134,2.018483,0.409721,5.709933,0.349846,0.208041,0.053363,0,0,0.415775,0.061227,0,0.060421,0.000857,0.000119,0.94456,0.00007,0.368538,0.02623,1.018005,0.015001,0.082373,0.02192,1.660885,0.001302,0.000589,0,0.360575,0,0.296014,0.048902,2.260424,0.853779,31.915858,8.373639,0,19.459275,1.243753,0.259909,1.523251,0.256174,2.397415,0.441084,1.316696,0.511968,0.05625,0.002515,3.879355,0.011433,0.324162,0.019665,2.735831,0.049848,0.092988,0.003354,1.164866,0.000838,0.000534,0,0.46677,0.00061}
		,{0.008032,0.036489,0.000062,0.586818,0.361164,1.562591,0.516594,4.919174,0.042119,0.177757,0,0,0.053874,0.173298,0,0.362681,0.000264,0.000595,0.000408,0.733202,0.079712,0.435243,0.083475,0.962295,0.076938,0.10308,0.002854,1.766961,0.004593,0.014243,0.002911,2.985421,0.090674,0.311759,0.154441,1.376727,12.116657,28.470047,19.459275,0,0.224397,0.822133,0.199589,1.144639,0.563566,2.871848,0.087905,0.765495,0.004868,0.010298,0.002528,3.0925,0.108501,0.498504,0.318574,2.866325,0.042595,0.095207,0.003464,1.447392,0.000374,0.006459,0.000187,0.716441}
		,{0.263567,0.000094,0.271628,0.077878,1.773102,0.000929,0.872084,0.040706,0.74787,0.042762,0,0,0.360038,0.000135,0,0.074859,0.25938,0,0.019568,0.00014,0.78734,0.000192,0.096104,0.002705,2.691226,0.188587,1.759732,0.206851,0.682254,0.000068,0.009981,0.000981,0.239388,0.014351,0.256283,0.208924,2.057449,0.067554,1.243753,0.224397,0,17.634677,2.075154,0.524647,4.971507,0.643789,0.93684,0.454347,0.318595,0,0.008844,0.000628,0.64316,0.000045,0.08462,0.003771,3.529593,0.167041,0.975582,0.149578,1.322234,0.001122,0.909861,0.194964}
		,{0.000093,0.143614,0.00284,0.060699,0.00339,1.118208,0.187826,0.725836,0.046586,0.487814,0,0,0.000932,0.325422,0,0.053908,0.00033,0.18788,0.01454,0.034916,0.000959,0.532092,0.074161,0.095418,0.46097,2.203539,0.377447,0.985145,0.00318,0.863191,0.016636,0.065212,0.025405,0.175491,0.020837,0.21917,0.296066,1.346385,0.259909,0.822133,17.634677,0,0.413957,1.038682,0.586685,4.127466,0.273855,0.600415,0.000295,0.215492,0.002362,0.013582,0.001329,0.565968,0.041435,0.068501,0.371685,2.729647,0.193544,0.372719,0.005905,1.110726,0.2804,0.293884}
		,{0.148268,0.005996,0.61266,0.004963,0.340991,0.020909,1.496628,0.000459,0.467136,0.042942,0,0,0.099519,0.004316,0,0.051005,0.060922,0.002143,0.43362,0.000035,0.247195,0.012394,1.042457,0.00041,0.899262,0.143479,3.215817,0.285384,1.769879,0.296358,3.06551,0.011032,0.221536,0.02302,0.667397,0.275355,0.878163,0.089476,1.523251,0.199589,2.075154,0.413957,0,12.931524,1.29386,0.334224,5.815294,1.868194,0.20148,0.00193,0.708204,0.000193,0.216244,0.0055,0.923644,0.000482,0.604859,0.053168,2.258321,0.159248,0.04873,0.002863,0.713961,0.001158}
		,{0.013122,0.043609,0.000062,0.33378,0.156468,0.25165,0.004438,0.768607,0.072867,0.140864,0,0,0.011489,0.034794,0,0.105226,0.039362,0.022384,0.000679,0.133032,0.220816,0.303613,0.012002,0.561523,0.324525,0.571469,0.461383,2.285052,0.560831,2.721043,0.034519,3.8322,0.04144,0.116405,0.056267,0.497593,0.291009,0.623366,0.256174,1.144639,0.524647,1.038682,12.931524,0,0.389004,0.928876,1.197614,7.316623,0.072138,0.056145,0.000524,0.311885,0.268662,0.373948,0.004382,0.759132,0.188097,0.426684,0.308851,1.563846,0.021649,0.176224,0.00176,0.744}
		,{0.225554,0.000047,0.010929,0.009289,12.169045,0.083636,5.964323,0.681575,0.50647,0.00018,0,0,2.007768,0.181794,0,0.139046,0.322807,0.00006,0.009512,0.000105,1.035015,0.000288,0.021048,0.001066,1.293228,0.000453,1.177718,0.083261,0.772349,0.018146,0.530881,0.025006,0.359183,0.018727,0.269862,0.306375,4.943439,0.275865,2.397415,0.563566,4.971507,0.586685,1.29386,0.389004,0,28.579806,1.644621,1.477696,0.501456,0.000097,0.006305,0.001358,2.323091,0.000485,0.126382,0.006402,1.702817,0.000388,0.832592,0.129098,2.382451,0.054025,1.179053,0.684968}
		,{0.000093,0.166706,0.000432,0.00579,0.094762,10.892976,1.049877,9.818281,0.000116,0.34689,0,0,0.248099,1.372357,0,0.167138,0.00033,0.300513,0.002718,0.0173,0.001247,1.337629,0.005195,0.104681,0.00506,1.282522,0.232701,1.418383,0.387941,1.320875,0.354545,1.360141,0.03114,0.238533,0.021539,0.304581,0.622868,4.699375,0.441084,2.871848,0.643789,4.127466,0.334224,0.928876,28.579806,0,0.403913,1.28699,0.000219,0.381287,0.000657,0.00416,0.005328,2.214735,0.063353,0.200205,0.012481,2.005334,0.308372,0.822643,0.326035,2.392606,0.298738,1.149846}
		,{0.140516,0.0126,0.423774,0.001718,0.04789,0.002044,1.094736,0.000115,0.424321,0.060909,0,0,0.144388,0.030883,0,0.145245,0.004747,0.00006,0.409432,0,0.011511,0.000384,0.493508,0,0.411115,0.025544,1.24214,0.000289,17.450524,1.113671,31.949764,2.418859,0.116039,0.012331,0.780008,0.305714,0.432918,0.021698,1.316696,0.087905,0.93684,0.273855,5.815294,1.197614,1.644621,0.403913,0,43.916187,0.032116,0.000205,0.560642,0.000103,0.013852,0,0.461729,0,0.030474,0.000718,0.668173,0.00041,0.037657,0.000821,0.938439,0.069567}
		,{0.08331,0.056771,0.000247,0.506841,0.063994,0.028715,0.009261,0.514392,0.200499,0.26951,0,0,0.122186,0.070533,0,0.496706,0.00956,0.000952,0.000951,0.04032,0.060432,0.033244,0.009225,0.273876,0.140943,0.169022,0.003786,1.148387,6.798629,4.087042,15.287419,18.531553,0.093542,0.062369,0.317466,0.90581,0.309656,0.140701,0.511968,0.765495,0.454347,0.600415,1.868194,7.316623,1.477696,1.28699,43.916187,0,0.051709,0.002062,0.002062,0.314379,0.045374,0.001768,0.00442,0.187832,0.015763,0.008986,0.00442,1.19579,0.047437,0.012227,0.165587,1.558784}
		,{0.86397,0.065905,0.748196,0.529619,0.563995,0.000186,0.002219,0.000115,0.571505,0.000359,0,0,1.598824,0.004316,0,0.038763,1.89715,0.072866,0.555104,0.01158,0.708395,0.000192,0.009225,0.000246,0.338207,0.000091,0.150804,0.0096,0.336828,0,0.003743,0,6.546163,0.575921,2.577578,0.430124,2.898791,0.00302,0.05625,0.004868,0.318595,0.000295,0.20148,0.072138,0.501456,0.000219,0.032116,0.051709,0,10.956917,25.313949,8.832621,2.600428,0.046583,0.718719,0.054049,0.153418,0.004101,0.276499,0.049842,0.164143,0.00205,0.080337,0.032177}
		,{0.026338,0.946136,0.00599,0.140867,0.000085,0.396897,0.000096,0.001261,0.000116,0.582801,0,0,0.001087,1.273908,0,0.092044,0.105361,2.720153,0.043756,1.08517,0.000192,0.76009,0.000537,0.031642,0.000375,0.491034,0.006757,0.06932,0.000589,0.648523,0.001871,0.005026,0.458622,7.487816,0.154207,0.350473,0.001027,2.862216,0.002515,0.010298,0,0.215492,0.00193,0.056145,0.000097,0.381287,0.000205,0.002062,10.956917,0,5.637509,18.744445,0.131929,2.620833,0.092405,0.968351,0.007112,0.119062,0.042565,0.245019,0.016776,0.201477,0.009773,0.064227}
		,{0.565566,0.047403,2.299543,0.762425,0.001356,0.000279,0.81372,0.000115,0.236179,0.06043,0,0,0.754233,0.086986,0,0.058616,0.509595,0.03173,2.047159,0.020529,0.001151,0.000192,0.732111,0.000082,0.019586,0.002808,0.606945,0.000145,0.000353,0,0.255147,0,2.581654,0.313779,11.271062,0.569076,0.009008,0.001957,3.879355,0.002528,0.008844,0.002362,0.708204,0.000524,0.006305,0.000657,0.560642,0.002062,25.313949,5.637509,0,13.945647,0.662578,0.028569,3.415722,0.081861,0.078381,0.006776,0.469281,0.053017,0.072521,0.001557,0.324696,0.074536}
		,{0.068927,0.371072,0.024699,1.108802,0.000254,0.001394,0.000772,0.451899,0.00963,0.116309,0,0,0.126844,0.530278,0,0.277072,0.091844,0.477439,0.173122,2.26249,0.000384,0.001537,0.000717,0.602428,0.032237,0.044112,0.000466,0.475496,0.001413,0.003287,0.000832,0.701399,0.953353,3.487342,0.911193,1.399673,0.003635,0.089643,0.011433,3.0925,0.000628,0.013582,0.000193,0.311885,0.001358,0.00416,0.000103,0.314379,8.832621,18.744445,13.945647,0,0.152215,0.682579,0.415718,2.211488,0.011491,0.04128,0.055025,0.362328,0.024883,0.051048,0.016839,0.276276}
		,{0.483563,0.000234,0.008892,0.03563,3.994417,0.771771,1.825878,0.250545,0.370309,0.000419,0,0,1.899865,0.011733,0,0.03586,0.879741,0.000417,0.004077,0.000772,1.272426,0.073021,0.346978,0.013526,0.42058,0.000091,0.27295,0.040907,0.324227,0,0.001871,0,0.536332,0.000253,0.001795,0.451134,2.359362,0.120121,0.324162,0.108501,0.64316,0.001329,0.216244,0.268662,2.323091,0.005328,0.013852,0.045374,2.600428,0.131929,0.662578,0.152215,0,9.612709,24.400553,5.140068,0.396521,0.018617,0.502355,0.106257,1.572808,0.022214,0.658541,0.238907}
		,{0.000093,0.512203,0,0.001654,0.111968,3.45257,0.084218,1.978448,0.000058,0.380788,0,0,0.000466,1.514097,0,0.044178,0.000132,1.191514,0.000136,0.103766,0.022062,1.394892,0.013077,0.971148,0.000562,0.696016,0.018814,0.634927,0.000236,1.610248,0.000416,0.020471,0.00004,0.536362,0,0.169264,0.038559,2.159328,0.019665,0.498504,0.000045,0.565968,0.0055,0.373948,0.000485,2.214735,0,0.001768,0.046583,2.620833,0.028569,0.682579,9.612709,0,6.74656,19.373137,0.01514,0.802516,0.140546,0.938586,0.086923,1.797671,0.036022,0.496552}
		,{0.109975,0.016535,1.041312,0.019406,1.93135,0.5585,4.380679,0.505677,0.176829,0.034737,0,0,0.806554,0.297371,0,0.031466,0.002374,0.000357,0.741135,0.000351,0.335924,0.069754,1.236457,0.087384,0.039265,0.004982,1.274118,0.002219,0.000589,0.000068,0.286547,0.000123,0.002262,0.000589,0.983926,0.517896,0.381796,0.077341,2.735831,0.318574,0.08462,0.041435,0.923644,0.004382,0.126382,0.063353,0.461729,0.00442,0.718719,0.092405,3.415722,0.415718,24.400553,6.74656,0,11.561124,0.18909,0.027912,0.905488,0.157605,0.585071,0.027973,1.693998,0.672077}
		,{0.005884,0.074851,0,0.220908,0.103323,1.262618,0.150589,4.658653,0.027035,0.106187,0,0,0.028567,0.586111,0,0.446015,0.000066,0.000893,0,1.524024,0.014101,0.417565,0.017824,1.950083,0.080124,0.190037,0.001165,1.544626,0.001531,0.083744,0.000624,3.409178,0.000081,0.004629,0.000078,0.837302,0.023862,0.728891,0.049848,2.866325,0.003771,0.068501,0.000482,0.759132,0.006402,0.200205,0,0.187832,0.054049,0.968351,0.081861,2.211488,5.140068,19.373137,11.561124,0,0.043198,0.702594,0.227527,1.251589,0.083552,1.398079,0.046588,1.526141}
		,{0.064397,0,0.042112,0.038557,1.120532,0.003717,0.348448,0.117533,0.223763,0.015452,0,0,0.099985,0.000135,0,0.028249,0.129492,0,0.012366,0.000491,0.661776,0.000769,0.147873,0.03156,0.746792,0.046739,0.706782,0.130873,0.162525,0,0.00707,0.000368,0.066966,0.000042,0.001171,0.059065,0.928969,0.000559,0.092988,0.042595,3.529593,0.371685,0.604859,0.188097,1.702817,0.012481,0.030474,0.015763,0.153418,0.007112,0.078381,0.011491,0.396521,0.01514,0.18909,0.043198,0,14.214694,2.738552,1.091224,0.629243,0.037461,0.375097,0.235747}
		,{0,0.055366,0.000062,0.006808,0.000254,1.023142,0.007428,0.670108,0.010037,0.184704,0,0,0,0.071612,0,0.066384,0.000066,0.135255,0.001359,0.015686,0.000096,0.976175,0.003672,0.644235,0.100928,0.975727,0.121389,0.928319,0.000236,0.915505,0.009981,0.150527,0,0.032447,0,0.011379,0.000158,1.013424,0.003354,0.095207,0.167041,2.729647,0.053168,0.426684,0.000388,2.005334,0.000718,0.008986,0.004101,0.119062,0.006776,0.04128,0.018617,0.802516,0.027912,0.702594,14.214694,0,0.892903,3.195698,0.03512,1.228004,0.067431,0.403521}
		,{0.084945,0.006464,0.287373,0.005472,0.330481,0.08568,1.265487,0.002179,0.257122,0.043721,0,0,0.028878,0.003641,0,0.009966,0.03956,0.002679,0.313495,0.00014,0.184749,0.105112,0.890822,0.00541,0.452442,0.106069,3.081614,0.536567,0.034978,0.025678,0.440217,0.000858,0.038612,0.009174,0.361403,0.033994,0.251423,0.109664,1.164866,0.003464,0.975582,0.193544,2.258321,0.308851,0.832592,0.308372,0.668173,0.00442,0.276499,0.042565,0.469281,0.055025,0.502355,0.140546,0.905488,0.227527,2.738552,0.892903,0,12.984714,0.089148,0.030585,0.639125,0.136937}
		,{0.010974,0.034428,0,0.159955,0.04238,0.283432,0.001061,1.029128,0.042815,0.136432,0,0,0.014439,0.013216,0,0.137634,0.00422,0.010061,0.000136,0.1763,0.034437,0.294294,0.001791,0.99033,0.159217,0.566034,0.343314,3.036767,0.007891,0.528692,0.00104,2.171984,0.003312,0.031984,0.000078,0.262465,0.033581,0.360196,0.000838,1.447392,0.149578,0.372719,0.159248,1.563846,0.129098,0.822643,0.00041,1.19579,0.049842,0.245019,0.053017,0.362328,0.106257,0.938586,0.157605,1.251589,1.091224,3.195698,12.984714,0,0.030223,0.239725,0.053748,0.968146}
		,{0.164659,0.000141,0.000741,0.003881,0.976185,0.001951,0.011673,0.007109,0.13094,0.00012,0,0,0.420899,0.045044,0,0.039313,0.169777,0.00006,0.000272,0.000175,0.418802,0.000288,0.002508,0.001312,0.388156,0.000091,0.042812,0.003377,0.241197,0.004656,0.042005,0.011768,0.069995,0,0.000156,0.027479,0.380374,0.000112,0.000534,0.000374,1.322234,0.005905,0.04873,0.021649,2.382451,0.326035,0.037657,0.047437,0.164143,0.016776,0.072521,0.024883,1.572808,0.086923,0.585071,0.083552,0.629243,0.03512,0.089148,0.030223,0,13.93595,21.171295,13.981617}
		,{0,0.172889,0,0.000191,0.000085,0.880032,0.000289,0.356038,0.000058,0.127388,0,0,0.007608,0.309374,0,0.105305,0,0.240505,0,0.047268,0.000096,0.636916,0.00009,0.395771,0.000843,0.566759,0.016193,0.336277,0.021435,0.676049,0.008942,0.703728,0.000283,0.055425,0,0.018603,0,0.518903,0,0.006459,0.001122,1.110726,0.002863,0.176224,0.054025,2.392606,0.000821,0.012227,0.00205,0.201477,0.001557,0.051048,0.022214,1.797671,0.027973,1.398079,0.037461,1.228004,0.030585,0.239725,13.93595,0,5.214689,18.675227}
		,{0.113991,0.018315,0.201112,0.001082,0.012121,0.001951,1.720919,0.00172,0.082323,0.029826,0,0,0.197641,0.061497,0,0.073682,0.00033,0.00006,0.165784,0.00007,0.003549,0.000384,0.556204,0.000164,0.097554,0.004982,0.551493,0.000289,0.01531,0.000753,0.247245,0.010419,0.000283,0.000084,0.194319,0.037724,0.002449,0.000112,0.46677,0.000187,0.909861,0.2804,0.713961,0.00176,1.179053,0.298738,0.938439,0.165587,0.080337,0.009773,0.324696,0.016839,0.658541,0.036022,1.693998,0.046588,0.375097,0.067431,0.639125,0.053748,21.171295,5.214689,0,25.64086}
		,{0.018773,0.032039,0,0.175861,0.002797,0.002974,0.003376,2.163175,0.007948,0.014314,0,0,0.105884,0.183952,0,0.381671,0.000066,0.000119,0,0.185038,0.001918,0.001441,0.001254,0.703092,0.08406,0.053714,0.003029,0.634203,0.043222,0.097165,0.143481,0.590833,0.000081,0.000295,0.000078,0.410199,0.000553,0.000447,0.00061,0.716441,0.194964,0.293884,0.001158,0.744,0.684968,1.149846,0.069567,1.558784,0.032177,0.064227,0.074536,0.276276,0.238907,0.496552,0.672077,1.526141,0.235747,0.403521,0.136937,0.968146,13.981617,18.675227,25.64086,0}};


		double ECMunrestPi[64]=
		{0.021414,0.021349,0.016195,0.015717,0.011798,0.010761,0.010366,0.008721,0.017237,0.016697,0,0,0.006441,0.007415,0,0.012744,0.015167,0.016798,0.007359,0.028497,0.010425,0.010408,0.011165,0.012199,0.010671,0.01104,0.017168,0.02073,0.008491,0.014604,0.004809,0.008158,0.024759,0.023762,0.012814,0.02118,0.012656,0.017882,0.01312,0.010682,0.022276,0.020321,0.03109,0.026699,0.01031,0.013701,0.009746,0.006788,0.01902,0.018419,0.010921,0.022626,0.018907,0.026817,0.016516,0.018288,0.02859,0.025285,0.034527,0.030606,0.016883,0.023659,0.016386,0.010223};

		if(basefreqs.empty()) for(int gv3=0; gv3<64; gv3++) basefreqs.push_back(ECMunrestPi[gv3]); 
			
		for(int gv1=0; gv1<64; gv1++) {row.clear(); {for(int gv2=0; gv2<64; gv2++) row.push_back(ECMunrestQ[gv1][gv2]);} myQvec.push_back(row);}

		return myQvec;
}

vector<vector<double> > getCOD(string &name, vector<double> &basefreqs, int mymodel, double kappa, double omega)
{
	//	make Q matrix for codon models

	vector<vector<double> > myQvec;
	vector<double> row;


		for(int xh1=0; xh1<4; xh1++)
		{for(int xh2=0; xh2<4; xh2++)
		{for(int xh3=0; xh3<4; xh3++)
		{

		for(int yh1=0; yh1<4; yh1++)
		{for(int yh2=0; yh2<4; yh2++)
		{for(int yh3=0; yh3<4; yh3++)
		{
			int xcodonNumber  = (xh1<<4)+(xh2<<2)+xh3; //getCodonNumber2(xh1,xh2,xh3);
			int ycodonNumber  = (yh1<<4)+(yh2<<2)+yh3; //getCodonNumber2(yh1,yh2,yh3);
 
			double matrixrate=0;
 
			char xcodon=GeneticCodeTable[geneticcode][xcodonNumber];
			char ycodon=GeneticCodeTable[geneticcode][ycodonNumber];

		    if(xcodon!='*' && ycodon!='*')
			{				
				if(yh1==xh1 && yh2==xh2 && yh3==xh3) matrixrate=0;

				else 
				if(yh1==xh1 && yh2==xh2)
				{
					if(xcodon!=ycodon) matrixrate=omega; else matrixrate=1; 
//	 				if( (yh3==0 && xh3==1) || (yh3==1 && xh3==0) || (yh3==3 && xh3==2) || (yh3==2 && xh3==3) )  matrixrate*=kappa;
					if((yh3+xh3-1)*(yh3+xh3-5)==0)  matrixrate*=kappa;  
				}
				else if(yh1==xh1 && yh3==xh3)
				{
					if(xcodon!=ycodon) matrixrate=omega; else matrixrate=1; 
//					if( (yh2==0 && xh2==1) || (yh2==1 && xh2==0) || (yh2==3 && xh2==2) || (yh2==2 && xh2==3) )  matrixrate*=kappa;
					if((yh2+xh2-1)*(yh2+xh2-5)==0)  matrixrate*=kappa;  
				}
				else if(yh2==xh2 && yh3==xh3) 
				{
					if(xcodon!=ycodon) matrixrate=omega; else matrixrate=1; 
//	 				if( (yh1==0 && xh1==1) || (yh1==1 && xh1==0) || (yh1==3 && xh1==2) || (yh1==2 && xh1==3) )  matrixrate*=kappa;
					if((yh1+xh1-1)*(yh1+xh1-5)==0)  matrixrate*=kappa;  
				}
			}
			//else cout<<xcodon<<" "<<ycodon<<"     ";   
			
			row.push_back(matrixrate);
			}}} // end of second for triple

		myQvec.push_back(row);
		row.clear();
		}}} // end of first for triple
		

		



		/*
		// checks that the matrices are symmetrical 
		for ( i=0; i<64; i++)
		{						CPi[i]=mybasefreqs.at(i);
		for (j=i+1; j<64; j++)
		{ 
			diff = Cmatrix[i][j] - Cmatrix[j][i];
			if (diff < 0.0)
				diff = -diff;
			if (diff > 0.001)  
				if(ECMrest)cout<<" ERROR: ECMrestCodon model is not symmetrical before frequencies.\n";
				else if(ECMunrest)cout<<" ERROR: ECMunrest Codon model is not symmetrical before frequencies.\n";
				else cout<<" ERROR: User Codon model is not symmetrical before frequencies.\n";
				
				
			}
		}
		*/

		int i,j;
		// rescale stationary frequencies, to make certain they sum to 1.0 
		double sum = 0.0;
		for (i=0; i<64; i++) sum += basefreqs.at(i); //cout<<basefreqs.at(i)<<" ";} cout<<endl;
		if(sum!=1) for (i=0; i<64; i++) basefreqs.at(i) /= sum;
			
		
		// multiply entries by stationary frequencies 
		for (i=0; i<64; i++) for (j=0; j<64; j++) (myQvec.at(i)).at(j) *= basefreqs.at(j);
		

		// rescale, so branch lengths are in terms of expected number of substitutions per site 
		double scaler = 0.0;
		for ( i=0; i<64; i++)
		{
			for ( j=i+1; j<64; j++)
			{
				scaler += basefreqs.at(i) * (myQvec.at(i)).at(j);
				scaler += basefreqs.at(j) * (myQvec.at(j)).at(i);
			}
		}

		//set diagonal of matrix
		for( i=0; i<64; i++) {sum=0; for(j=0; j<64; j++) if(i!=j) sum+=(myQvec.at(i)).at(j); (myQvec.at(i)).at(i)=-sum; }
		
		scalefactors.push_back(scaler);

	return myQvec;

}
//////////////////////////////////////////////////////////////////////////////////////////




}; // end of modelclass

///////////////////
/*
double main(int argc, char* argv[])
{	

	return 0;
}


*/

