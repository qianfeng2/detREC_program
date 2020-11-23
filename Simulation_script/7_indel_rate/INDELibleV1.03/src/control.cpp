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


#include <iostream>
#include <iomanip> 
#include <fstream> 
#include <string>
#include <vector> 
#include <math.h>
#include <sstream>  
#include <algorithm>
#include <map> 
#include "models.cpp"

using namespace std;

bool controldebug=false;   //"false" will cease everything when a -1 error occurs, "true" will continue on - for debugging (but may lead to crashes)

extern string masterfilename;//="control.txt";

bool oldmethod;				// false for method 1, true for waiting time method

bool breakonerror=true;   // "true" will allow "return -1;" statements after an error.
                          // "false" will continue on and probably make the program crash but is useful
                          // for checking errors in the control file programs.

#pragma warning(disable:4786)

bool globalseed=false;				// whether one seed sets both random trees, and sequence generation.

bool insertaslowercase   = false;    // prints inserted sites as lower case and core sites as uppercase

bool fixtrueproportions = false;   //  if true then the proportions of codon sites models will be true rather than random.

bool markdeletedinsertions = false; //
bool fileperrep          = false;	// prints all replicates in same file if false, otherwise in their own files.

bool printrates			 = false;	// whether to print out relative rate information or not.
bool printallrates			 = false;	// whether to print out relative rate information or not for every branch

bool printcodonsasDNA=true;		// whether to print codons as DNA bases or not (true for DNA/false for amino acids)


bool printinsert		 = false;	// whether to print indel statistics
bool ancestralprint		 = false;	// whether to print ancestral sequences to another file
bool ancestralfile		 = true;	// whether ancestral sequences go in same or different file.  true for different file
//bool replicatesfile		 = true;	// whether replicate datasets should be printed in the same file or into separate files.
bool phylipnametruncate	 = false;	// truncates taxa names to 10 characters to comply with PHYLIP format
bool guidetreebinary	 = true;	// whether guide tree contains polytomies
bool guidetreerooted	 = true;	// whether guide tree is unrooted or rooted

string phylipextension	="phy";		// used for output file filename extension
string fastaextension	="fas";		// used for output file filename extension
string nexusextension	="nex";		// used for output file filename extension

int outputtype			= 2;		// 1 for FASTA, 2 for PHYLIP, 3 for NEXUS
int guidetreetype		= 0;		// used for logfile


// makes a time string to append to output file names.
string makesimtime()
{
	time_t stringtime;

	struct tm *mytime; time(&stringtime);
	
	mytime=localtime(&stringtime); 

	string s1="_",s2=asctime(mytime); //hexToAsci(%3A)
	for(int i=0; i<s2.size()-1; i++) {if(s2[i]==':') s1+="."; else s1+=s2[i]; }
	for(int j=0; j<s1.size(); j++) if(s1[j]==' ') s1[j]='_';
	return s1;
}

string simtime; //=makesimtime();   //disabled time stamp so that output is over-written.

ofstream* LOG=new ofstream;		// used to print out information from simulations


// need "total" vectors for models, sites, branches, and modelnames, sitenames, branchesnames etc

vector<char> originalcontrol;  // used to store control file used.

vector<string>		totalmodelnames;
vector<string>		totalsitenames;
vector<string>		totalbranchesnames;
vector<string>		totaltreenames;
vector<string>		totalpartitionnames;
vector<model>		totalmodels;


// blocktype e.g. [MODEL], blockname #de1#  (* ??? *), commandname= [submodel], instring error print, myline last input
void controlerrorprint2(string blocktype, string blockname, string commandname, string instring, string myline)
{
	// this functions provides a framework for a standard output of an error to the screen and log file
	// mostly it is just formatting the border of the error box and white space etc

//	cout<<"\r                                                                      \n ERROR occurred in "<<blocktype<<" block "<<blockname<<".\n See below for details or consult LOG.txt                       ";
	cout<<"\r                                                                      \n";

	vector<string> toprint;
	
	string startline="ERROR in "+blocktype+" block ";
	if(blockname!="")   startline+=blockname;
	if(commandname!="") startline+=" in command "+commandname;

//	stringstream fd1; fd1<<linecount; string fd1s=fd1.str();
//	startline+=fd1s; startline+=" of control file:";
	
	toprint.push_back(startline);

	string tempstring;
	char c;
	int themaxsize=startline.size();
	
	for(int j0=0; j0<instring.size(); j0++) 
	{
		c=instring[j0]; 
		if(c=='\n') 
		{
			toprint.push_back(tempstring);
			if(themaxsize<tempstring.size()) themaxsize=tempstring.size();
			tempstring="";
		} 
		else tempstring+=c;
	} 
	toprint.push_back(tempstring);if(themaxsize<tempstring.size()) themaxsize=tempstring.size();
	
	if(myline!="") {string endline="Last Input read was: "; endline+=myline; if(themaxsize<endline.size()) themaxsize=endline.size(); toprint.push_back(endline);}

	cout<<endl<<endl; (*LOG)<<endl;

	for(int i0=0; i0<toprint.size(); i0++)
	{
		string tempstring2=toprint.at(i0);
		for(int h1=tempstring2.size(); h1<themaxsize; h1++) tempstring2+=" ";
		toprint.at(i0)=tempstring2;
	}

	cout<<endl<<" +";  (*LOG)<<endl<<"+";  
	for(int i1=0; i1<themaxsize+2; i1++) {cout<<"-"; (*LOG)<<"-";}
	cout<<"+"<<endl;	(*LOG)<<"+"<<endl;
	
	for(int i2=0; i2<toprint.size(); i2++) {cout<<" | "<<toprint.at(i2)<<" |"<<endl; (*LOG)<<"| "<<toprint.at(i2)<<" |"<<endl;}
	
	cout<<" +"; (*LOG)<<"+"; 
	for(int i3=0; i3<themaxsize+2; i3++) {cout<<"-"; (*LOG)<<"-";}
	cout<<"+"<<endl<<endl;	(*LOG)<<"+"<<endl<<endl;

	
}


int checkthere(string mytest,vector<string> mylist)
{
	// returns -1 if mytest is not in mylist otherwise returns position in list
	int blah=-1;
	for(int h1=0; h1<mylist.size(); h1++) if(mylist.at(h1)==mytest) {blah=h1; break;}
	return blah;
}

int checkthereN(int mytest,vector<int> mylist)
{
	// returns -1 if mytest is not in mylist otherwise returns position in list
	int blah=-1;
	for(int h1=0; h1<mylist.size(); h1++) if(mylist.at(h1)==mytest) {blah=h1; break;}
	return blah;
}


bool AinB(char testchar, string referencestring)
{
	// returns true if testchar is in referencestring
	// returns false if testchar is not in referencestring

	bool answer1=false;

	for(int i=0; i<referencestring.size(); i++) if(testchar==referencestring[i]) {answer1=true; break;}

	return answer1;
}

bool AonlyfromB(string teststring, string referencestring)
{
	// checks test string to make sure it is only made up of 
	// characters from referencestring. Anything else gives false.
	// returns false if teststring contains illegal characters

	bool answer;

	for(int i=0; i<teststring.size(); i++)
	{
		answer=AinB(teststring[i],referencestring);
		if(!answer) break;
	}

return answer;

}


bool allAinB(string teststring, string referencestring)
{
	// checks test string to make sure it is only made up of 
	// characters from referencestring. Anything else gives false.
	// returns false if teststring contains illegal characters

	bool answer;

	for(int i=0; i<teststring.size(); i++)
	{
		answer=AinB(teststring[i],referencestring);
		if(!answer) break;
	}

return answer;

}


bool noneAinB(string teststring, string referencestring)
{
	// checks test string to make sure it is only made up of 
	// characters from referencestring. Anything else gives false.
	// returns false if teststring contains illegal characters

	bool answer;

	for(int i=0; i<teststring.size(); i++)
	{
		answer=AinB(teststring[i],referencestring);
		if(answer) break;
	}

return !answer;

}


////////////////////////////////////////////////

class branchclass
{
	// this class provides the framework for simulations that use different models on different branches in the guide tree.

public:
	string treewithrootmodel;         // original tree given in control file, has a model at the root.
	string tree;                      // tree minus model at root.
	string baretree;                  // tree expressed as only a pattern of parentheses and commas
	string name;					  // name of branch class
	vector<string> modelnames;		  // names of models used, given in same order as modelpositions
	//vector<model> branchmodels;
	vector<int> modelpositions;       // position of models used in totalmodels.
	vector<int> allmodelsused;		  // list of all models used in branch class (listed by model number)
	vector<double> rootbasefreqs;     // base frequencies of model defined at root of branch class tree.
	double insertrate;                // insertion rate at root?
	double deleterate;                // deletion rate at root?
	int rootmodelpos;                 // position of root model in totalmodels.

	bool geneticcodefixed;			  // whether to allow geneticcode to change on trees.
	int geneticcode;                  // genetic code used (for codon simulations) must be the same across entire tree
	int error;                        // whether branch class successfully set itself up without errors.

	int numcats;					  // number of categories in discrete gamma model or in site/branch-site models
	                                  //    must be the same across the whole tree.

	bool changecats;				  // says whether to rescale root cats and therefore check all for change.
	
	vector<double> catprops;          // relative proportion of categories (must be same in different models).

	branchclass(string &ttree, string &myname, vector<string> &totalmodelnames, bool iscodefixed)
	{
		// constructor
		geneticcodefixed=iscodefixed;

		error=1;
		treewithrootmodel=ttree;
		for(int fd=0; fd<treewithrootmodel.size(); fd++) {char c=treewithrootmodel[fd]; if(c=='(' || c==')' || c==',' || c==';') baretree+=c;}
		name=myname;
		
		if(error!=-1) error=getrootfreqs(treewithrootmodel);
		if(error!=-1) error=buildbranches(tree,myname,modelnames);
			
		if(error!=-1) error=getbranchmodels(modelnames,totalmodelnames,modelpositions); //branchmodels);

		if(error!=-1 /*&& geneticcodefixed*/) error=testgeneticcode2(modelpositions);

	}

private:


	int getrootfreqs(string &testtree2)
	{
		// this function takes the original tree with root model and strips off the defined model at root.
		// it then checks that that model exists, and collects information about the model.
			
		int testtreesize=testtree2.size()-1;
		string testtree, rootmodel;
		int pospos;

		for(int pp=testtreesize; pp>0; pp--) if(testtree2[pp]==')') {pospos=pp; break;}

		for(int pp2=0; pp2<testtreesize; pp2++) 
		{
			char c1=testtree2[pp2];
			if(pp2>pospos) 
			{
				if(pp2!=pospos+1) rootmodel+=c1; 
			}
			else testtree+=c1;
		}
 		testtree+=';';

		if(rootmodel=="") {controlerrorprint2("[BRANCHES]", name ,"","You haven't named a model for the root of this branch class.\nThe branch class tree should be of the form: (...)#m;\nwhere #m represents the (previously defined) root model m.",""	); return -1;} 
		else
		{
			int ans=checkthere(rootmodel,totalmodelnames);

			if(ans==-1) {controlerrorprint2("[BRANCHES]", name ,"","The model "+rootmodel+" in this branches block has not been defined",""	); return -1; }
			else 
			{
				modelpositions.push_back(ans);
				model* n=&(totalmodels.at(ans));

			//	if(type==3) 
			//	{
					catprops=(*n).cumfreqs;
					numcats=catprops.size();
					if(numcats==1) changecats=true; else changecats=false;
			//	}
				
				insertrate=(*n).insertrate;
				deleterate=(*n).deleterate;
				rootbasefreqs=(*n).rootbasefreqs;
				rootmodelpos=(*n).modelpos;

			}
		}

		tree=testtree;

		return 0;
	}	

///////////////////////////

	int testgeneticcode(vector<int> modelpositions)
	{
		// this function simply tests that the genetic code defined in each model on each
		// branch in the guide tree is the same and throws an error if they are not.

		model* amodel=&(totalmodels.at(modelpositions.at(0)));
		int mygeneticcode=(*amodel).geneticcode;
		int lastgeneticcode;
		string n2;

		string n1=totalmodelnames.at(modelpositions.at(0));

		geneticcode=mygeneticcode;
		for(int ds=1; ds<modelpositions.size(); ds++)
		{
			lastgeneticcode=mygeneticcode;
			n2=n1;

			amodel=&(totalmodels.at(modelpositions.at(ds)));
			mygeneticcode=(*amodel).geneticcode;

			string n2=totalmodelnames.at(modelpositions.at(ds));

			stringstream dd1; dd1<<lastgeneticcode;   string d1=dd1.str();
			stringstream dd2; dd2<<mygeneticcode;     string d2=dd2.str();
			
			if(mygeneticcode!=lastgeneticcode) {controlerrorprint2("[BRANCHES]", name ,"","All models in a branches block should use the same genetic code.\nmodel "+n1+" has genetic code "+d1+"\nmodel "+n2+" has genetic code "+d2+"\nIf you want to allow the genetic code to change then use the command\n[BRANCHES*] instead of [BRANCHES] to turn off this safeguard.",""	); {if(breakonerror) return -1;} } 
		}

		return 0;
	}
///////////////////////////
	
	int testgeneticcode2(vector<int> modelpositions)
	{
		// this function simply tests that the genetic code defined in each model on each
		// branch in the guide tree is the same and throws an error if they are not.

		model* amodel=&(totalmodels.at(modelpositions.at(0)));
		int rootgeneticcode=(*amodel).geneticcode;
		vector<int> myallowedlist=allowedcodes(rootgeneticcode);
		int lastgeneticcode;
		
		stringstream dd2; dd2<<rootgeneticcode;     string d2=dd2.str();
			
		string n1=totalmodelnames.at(modelpositions.at(0)), n2;

		geneticcode=rootgeneticcode;
		for(int ds=1; ds<modelpositions.size(); ds++)
		{
			amodel=&(totalmodels.at(modelpositions.at(ds)));
			lastgeneticcode=(*amodel).geneticcode;

			string n2=totalmodelnames.at(modelpositions.at(ds));

			stringstream dd1; dd1<<lastgeneticcode;   string d1=dd1.str();

			bool codeerror=true;
			for(int qwe=0; qwe<myallowedlist.size(); qwe++)
			{
				if(lastgeneticcode==myallowedlist.at(qwe)) {codeerror=false; break;}
			}
			if(codeerror) {controlerrorprint2("[BRANCHES]", name ,"","All models in a branches block should use genetic codes that share \nthe same stop codons.\nmodel "+n1+" has genetic code "+d1+"\nmodel "+n2+" has genetic code "+d2+"\nThese codes do not share all stop codons and so cannot be combined. ",""	); {if(breakonerror) return -1;} } 	
		}

		return 0;
	}
///////////////////////////
	int getbranchmodels(vector<string> &modelnames,vector<string> &totalmodelnames, vector<int> &modelpositions) //vector<model> &branchmodels)
	{	
		// this function tests all the models found by buildbranches function. if the models do not exist, or if they
		// have the wrong number of categories, or the categories have different proportions then an error is thrown.

		for(int hgg=0; hgg<modelnames.size(); hgg++)
		{
			string thisname=modelnames.at(hgg);

			int found=checkthere(thisname,totalmodelnames);
			if(found==-1) {controlerrorprint2("[BRANCHES]", name ,"","The model "+thisname+" in this branches block has not been defined",""	); return -1;}
			else 
			{
	//			if(type==3) 
	//			{
					model *k=&(totalmodels.at(found));
			
					int x=((*k).cumfreqs).size();
				
					if( numcats!=1 && x==1 && type!=3)
					{
						if(numcats!=0)
						{
							totalmodels.push_back(totalmodels.at(found));
						
							model *k2=&(totalmodels.back());

							totalmodelnames.push_back(  ((*k2).name)+"_with_all_classes_equal"  );
					
							//for(int pq=0; pq<((*k2).cumfreqs).size(); pq++) cout<<"1 "<<pq<<" "<<((*k2).cumfreqs).at(pq)<<endl;
							(*k2).cumfreqs=(totalmodels.at(rootmodelpos)).cumfreqs;
							//for(int pq=0; pq<((*k2).cumfreqs).size(); pq++) cout<<"2 "<<pq<<" "<<((*k2).cumfreqs).at(pq)<<endl;

							(*k2).changeQandJ(numcats);

							modelnames.at(hgg)=(*k2).name;

							found=totalmodels.size()-1;
						}
					}
					else
					{
				

						if( numcats!=x )  
						{

							stringstream s1,s2; s1<<x; s2<<numcats; string ss1=s1.str(), ss2=s2.str(); 
							string nn1=totalmodelnames.at(found); string nn2=totalmodelnames.at(rootmodelpos);
							controlerrorprint2("[BRANCHES]", name ,"","Number of site/rate classes in model "+nn1+" is "+ss1+"\nNumber of site/rate classes in root model "+nn2+" is "+ss2+"\nThese should be equal.",""); 
								
							return -1;
						}


						
						for(int op=0; op<numcats; op++) 
						{
							// cout<<"\t"<<catprops.at(op)<<"\t"<<((*k).cumfreqs).at(op)<<endl;

							if( catprops.at(op)!=((*k).cumfreqs).at(op) )  
							{
								stringstream s1,s2,s3; s1<<((*k).cumfreqs).at(op); s2<<catprops.at(op); s3<<op+1; 
								string ss1=s1.str(), ss2=s2.str(), op2=s3.str(); 
								string nn1=totalmodelnames.at(found); string nn2=totalmodelnames.at(rootmodelpos);
								controlerrorprint2("[BRANCHES]", name ,"","Site/rate class "+op2+" in model "+nn1+" has proportion "+ss1+"\nSite/rate class "+op2+" in root model "+nn2+" has proportion "+ss2+"\nThese should be equal.",""); 
									
								return -1;
							}

						}

					}
			//	}
				
				modelpositions.push_back(found);
				int minicheck=checkthereN(found,allmodelsused);
				if(minicheck==-1) allmodelsused.push_back(found);
			
			}
			
			//cout<<"model found "<<thisname<<" "<<found<<endl;
		}

		return 0;

	}
	///////////////

		int buildbranches(string &testtree, string &branchesname, vector<string> &modelnames)
		{
			// this function navigates the branch class guide tree defined in the control file
			// and extracts model names subject to further testing, whether they exist and so on...


			char c1='Q';			

			int bracketcount=0;
			bool istherehash=false;
			for(int i=0; i<testtree.size(); i++)
			{
				c1=testtree[i];

				if(c1=='#') istherehash=true;
				if(c1=='(') bracketcount++;
				if(c1==')') bracketcount--;
			}
			if(bracketcount!=0) {controlerrorprint2("[BRANCHES]", name ,"","Number of parentheses in tree for this branch class do not match",""); return -1;}

			//cout<<"Number of parentheses in tree for branch class "<<branchesname<<" do not match"<<endl;
			//if(c1!=';') cout<<"Tree for branch class "<<branchesname<<" must be in newick format and end in a ;"<<endl;
			
			if(testtree.size()>0 && testtree[testtree.size()-1]==';')
			{
				bracketcount=1;
				string rfstring="";
				

				for(int rf=1; rf<testtree.size(); rf++)
				{
					c1=testtree[rf];
					if(c1=='(') bracketcount++;
					if(c1==')') bracketcount--;
			//		if(bracketcount==0 ||(bracketcount==1 && (c1==','||c1==')')) )
					if((bracketcount==1 && c1==',') || bracketcount==0) 
					{
						buildbranches(rfstring,branchesname,modelnames);
						rfstring="";
					}
					else if(c1!=' ') rfstring+=c1;

				}
			}
			else if(istherehash)
			{
				bracketcount=1;
				string teststring, modelnamestring;
				vector<string>remaining;

				bool modelnamewrite=false;

				if(testtree[0]!='(') 
				{
					for(int i1=0; i1<testtree.size(); i1++)
					{
						///////////////////////////////////////////////////////////////////
						c1=testtree[i1]; 
						
						//if(modelnamewrite) cout<<"Q"<<c1<<"Q ";
						//if((c1==','||c1==')') && modelnamewrite) {modelnames.push_back(modelnamestring); break;}
						//if(c1=='#' && modelnamewrite) {modelnames.push_back(modelnamestring); break;}
			
						if(modelnamewrite) {modelnamestring+=c1;}		// cout<<"X "<<modelnamestring<<endl;}
						if(c1=='#' && !modelnamewrite) modelnamewrite=true;
					}
					modelnames.push_back(modelnamestring); 
	
				}
				else
				{
					for(int i=1; i<testtree.size(); i++)
					{
						c1=testtree[i];

						if(c1=='(') bracketcount++;
						if(c1==')') bracketcount--;

				//		if(bracketcount==0 ||(bracketcount==1 && (c1==','||c1==')')) )			
						if((bracketcount==1 && c1==',') || bracketcount==0) 
						{
							//	cout<<"AA "<<teststring<<endl;
							remaining.push_back(teststring);
							teststring="";
							if(bracketcount==0) bracketcount=-99;
						}
						else if(bracketcount==-99)
						{
								//if((c1==','||c1==')') && modelnamewrite) { modelnames.push_back(modelnamestring); break;}
								//if(c1=='#' && modelnamewrite) {modelnames.push_back(modelnamestring); break;}
								
								//if(modelnamewrite) cout<<"Q2"<<c1<<"Q ";

							if(modelnamewrite) {modelnamestring+=c1; }		//cout<<"X2 "<<modelnamestring<<endl;}
							if(c1=='#' && !modelnamewrite) modelnamewrite=true;
						}
						else teststring+=c1;
					}
					modelnames.push_back(modelnamestring);
				}
				for(int kd=0; kd<remaining.size(); kd++) buildbranches(remaining.at(kd),branchesname,modelnames);
					
				

			}

			return 0;
		}	
};
///////////////////////////////////////////////////////////

vector<branchclass> totalbranches;			// storage for branch classes




class siteclass
{
	// site class is now disabled.  It allows the use of different models in different proportions along a sequence.

public:

	vector<int> allmodelsused;
	vector<bool>   trueformodel;  //true for model, false for branch class, in mbnames, mbpositions
	vector<int>    mbpositions;
	vector<string> mbnames;
	vector<double> props;  // cumulative props remember!
	vector<vector<double> > trootbasefreqs;

	int geneticcode;
	int error;
	string mbtree;
	string name;

	bool therearemodels;
	bool therearebranches;
	
//	siteclasscopy(site copysite)
//	{
//		//??????
//	}

//	createsiteprops()
//	{
//		//????
//	}
	
	siteclass(string sitename, vector<string> &mynames, vector<double> &myprops, vector<string> &totalmodelnames, vector<string> &totalbranchnames)
	{
		error=1;
		therearemodels=false;
		therearebranches=false;
		vector<string> mbtrees;

		name=sitename;
		props=myprops;
		mbnames=mynames;
	
		if(error!=-1) error=sortoutnames(mbtrees,totalmodelnames,totalbranchnames);
		if(error!=-1) error=sortouttrees(mbtrees);
		if(error!=-1) error=testgeneticcode();

	}


///////////////////////////////////
private:

	
	int testgeneticcode()
	{
		int lastgeneticcode;
		string modelorbranch1="model";
		string modelorbranch2="model";
		string n1,n2;

		model* amodel;

		branchclass* abranchclass;

		if(trueformodel.at(0))
		{	
			amodel=&(totalmodels.at(mbpositions.at(0)));
			geneticcode=(*amodel).geneticcode;
		}
		else
		{	
			abranchclass=&(totalbranches.at(mbpositions.at(0)));
			geneticcode=(*abranchclass).geneticcode;
			
			modelorbranch1="branch class";
			modelorbranch2="branch class";
		}


		for(int ds=0; ds<mbpositions.size(); ds++)
		{
			lastgeneticcode=geneticcode;
			n2=n1;
			int pos=mbpositions.at(ds);
			modelorbranch2=modelorbranch1;
			if(trueformodel.at(ds)) 
			{
				modelorbranch1="model"; 
				n1=totalmodelnames.at(pos);
				amodel=&(totalmodels.at(pos));
				geneticcode=(*amodel).geneticcode;
			}
			else 
			{
				modelorbranch1="branch class"; 
				n1=totalbranchesnames.at(pos);
				abranchclass=&(totalbranches.at(pos));
				geneticcode=(*abranchclass).geneticcode;
			}

			
			
			stringstream dd1; dd1<<lastgeneticcode; string d1=dd1.str();
			stringstream dd2; dd2<<geneticcode;     string d2=dd2.str();
		
			if(geneticcode!=lastgeneticcode) {controlerrorprint2("[SITES]", name ,"","All models/branches in a sites block should use the same genetic code.\n"+modelorbranch1+" "+n1+" has genetic code "+d1+"\n"+modelorbranch2+" "+n2+" has genetic code "+d2,""	); {if(breakonerror) return -1;} } 
			
		}
		return 0;
	}

	int sortoutnames(vector<string> &mbtrees, vector<string> &totalmodelnames, vector<string> &totalbranchnames)
	{
		for(int gh=0; gh<mbnames.size(); gh++)
		{
			bool modelorbranch=true;
			string mbname=mbnames.at(gh);
			string tree,stree;

			int ans=checkthere(mbname,totalmodelnames);
			
			//cout<<"ans "<<ans<<endl;
			if(ans==-1) 
			{
				ans=checkthere(mbname,totalbranchnames); 

				
			//	cout<<"\tans "<<ans<<endl;
				modelorbranch=false;
				
				if(ans==-1) {controlerrorprint2("[SITES]", name, "", "You specified \""+mbname+"\" in this sites class.\nNo model/branch class with this name has been defined yet.",""); {if(breakonerror) return -1;} }
				else 
				{
					branchclass* b=&(totalbranches.at(ans));
					
					trootbasefreqs.push_back( (*b).rootbasefreqs );
		
					tree=(*b).tree; 
					therearebranches=true;
					
					vector<int> modpos=((*b).allmodelsused);
					for(int o1=0; o1<modpos.size(); o1++)
					{
						int y=modpos.at(o1);
						int minicheck=checkthereN(y,allmodelsused);
						if(minicheck==-1) allmodelsused.push_back(y);	
					}

					

				}

				for(int yf=0; yf<tree.size(); yf++) {char c=tree[yf]; if(c=='(' || c==')' || c==',' || c==';') stree+=c;}

				if(mbtree=="") mbtree=stree;

			}
			else 
			{
				therearemodels=true;
				int minicheck=checkthereN(ans,allmodelsused);
				if(minicheck==-1) allmodelsused.push_back(ans);
				model* m=&(totalmodels.at(ans));
				trootbasefreqs.push_back( (*m).rootbasefreqs );
			}

			trueformodel.push_back(modelorbranch);
			mbpositions.push_back(ans);
			
			mbtrees.push_back(stree);
		}

		return 0;
	}
	/////////////////
	int sortouttrees(vector<string> &mbtrees)
	{

		int lastpos=-1;
		for(int fc=0; fc<mbtrees.size(); fc++)
		{
			string testtree=mbtrees.at(fc);

			if(testtree!="")
			{
				if(lastpos==-1) lastpos=fc;
				else if(mbtree!=testtree)
				{
					string s1=mbnames.at(lastpos), s2=mbnames.at(fc), t1=mbtree, t2=testtree, u1;
					int s1s=s1.size(), s2s=s2.size();
						
					if(s1s>s2s) {for(int fv=0; fv<s1s-s2s; fv++) s2+=" ";} 
					else        {for(int fv=0; fv<s2s-s1s; fv++) s1+=" ";}

					controlerrorprint2("[SITES]", name, "", "Two branch classes in this site class have different guide trees.\n1) "+s1+"  "+t1+"\n2) "+s2+"  "+t2+"\nTrees must have their taxa/branches written in identical order.",""); 
					{if(breakonerror) return -1;} 
				}
			}
		}
		return 0;
	}
	////////
};
///////////////////////////////////////////////////////////
vector<siteclass>	totalsites;	// storage for site classes



///////////////
int teststring(int lastresult,  string com, vector<string> &allowed, string blocktype, string allowedstring, string blockname)
{
	// test function used by many other functions for testing "chunks" in the control file until they are processed.

	// returns -1 if com not a statement beginning with [
	// reurns x if com is the xth command in allowed
	// returns x if currently on command x and input is not a new command or non-decimal number.
	// returns -2 if com begins with [ but is not in allowed
	// returns -2 if com contains characters not in allowedstring.
	int myresult=-1;
	
	if(com[0]=='[') 
	{	
		for(int k=0; k<allowed.size(); k++) if(com==allowed.at(k)) {myresult=k; break;}

		if(myresult==-1) 
		{
			controlerrorprint2(blocktype, blockname,"" ,"Unknown command.",com); 
		
			myresult=-2;
		}
	}
	else
	{
		if(!AonlyfromB(com,allowedstring)) 
		{
			
			//	cout<<com<<"  "<<allowedstring<<endl;
	//		if(allowedstring!="0123456789.") myresult=lastresult;
			
			myresult=-2; 

			string ds; if(lastresult!=-999) ds=allowed.at(lastresult);
			
//			cout<<blocktype<<" "<<com<<" teststring func 1"<<endl;
			if(blocktype=="[TREE]" && ds=="[user]") myresult=lastresult;
			else 
			{
				ifstream gf1; gf1.open(com.c_str()); 
				if(gf1.good()) myresult=lastresult;	   // i.e. if it is a filename and the file exists and can be opened.
					else controlerrorprint2(blocktype, blockname,ds ,"Input value must be a positive decimal number or a command [..]\nIf this was meant to be a filename no such file exists.",com);
			}

		} else myresult=lastresult;
	}
//	cout<<"myresult "<<myresult<<endl;
	return myresult;
}
////////////////////////////////////////////////////////////////////////////////////

int dealwithsettings(vector<string> &block)	
{
	// this function deals with [SETTINGS] blocks in the control file.

	//for(int i=0; i<block.size(); i++) cout<<"settings "<<block.at(i)<<endl;
	
	int size=block.size();

	if(size==0) {controlerrorprint2("[SETTINGS]","","","The settings block seems to be empty.",""); {if(breakonerror) return -1;} }

	if(size%2!=0) {controlerrorprint2("[SETTINGS]","","","The settings block must contain an even number of statements:\n[COMMAND1] value1\n[COMMAND2] value2\nand so on....",""); {if(breakonerror) return -1;} }

	string s1,s2;
	string scommands[15]={"[printindelstatistics]","[ancestralprint]","[phylipextension]","[nexusextension]","[fastaextension]","[output]","[randomseed]","[printrates]","[insertaslowercase]","[printcodonsasaminoacids]","[fileperrep]","[markdeletedinsertions]","[fixproportions]","[printallrates]","[globalseed]"};
	vector<string> commands(scommands,scommands+15); 

	int idum2;
	int commanderror;
	
	for(int j=1; j<size; j=j+2)
	{
		s1=block.at(j-1);		
		s2=block.at(j);
		if(!allAinB(s2,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[SETTINGS]", "", "", "Values after commands should only contain ,.;()-_ \nand alpha-numeric characters.",""); {if(breakonerror) return -1;} }

		commanderror=-1;
		
		if(s1[0]!='[' || s1[s1.size()-1]!=']') {controlerrorprint2("[SETTINGS]","","","Was expecting a command in square brackets e.g. [COMMAND]\nbut found: "+s1,""); {if(breakonerror) return -1;} }
		else commanderror=checkthere(s1,commands);

		if(commanderror==-1) {controlerrorprint2("[SETTINGS]","","","The command "+s1+" is not valid in a setting block.\nPlease consult manual or example control file.",""); {if(breakonerror) return -1;} }
		
		else if(commanderror==0)
		{
			//printindelstatistics
			//myinclude.at(25)=true;

			if(s2=="TRUE") printinsert=true; else if(s2=="FALSE") printinsert=false; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==1)
		{
			//ancestralprint
			//myinclude.at(26)=true;

			     if(s2=="SAME") {ancestralprint=true; ancestralfile=false;} 
			else if(s2=="NEW")  {ancestralprint=true; ancestralfile=true;} 
			else if(s2=="FALSE") ancestralprint=false; 
			else {controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be SAME, NEW or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==2) phylipextension=s2;			
		
		else if(commanderror==3) nexusextension=s2;			
		
		else if(commanderror==4) fastaextension=s2;			

		else if(commanderror==5)
		{
			//output
			//myinclude.at(22)=true;
			
			     if(s2=="FASTA"   ) {outputtype=1; phylipnametruncate=false;}
  		    else if(s2=="PHYLIP"  ) {outputtype=2; phylipnametruncate=false;}
			else if(s2=="PHYLIPT" ) {outputtype=2; phylipnametruncate=true; }
			else if(s2=="NEXUS"   ) {outputtype=3; phylipnametruncate=false;}
			else {controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be FASTA, PHYLIP, PHYLIPT or NEXUS for this command.",""); {if(breakonerror) return -1;} }					
		}
		else if(commanderror==6)
		{
			//randomseed

			if(!allAinB(s2,"1234567890-")) { controlerrorprint2("[SETTINGS]", "",scommands[commanderror], "This value must be an integer.",""); {if(breakonerror) return -1;} }

			idum=atoi(s2.c_str()); if(idum<0) idum=-idum; 	mtrand1.seed(idum);	
		}
		else if(commanderror==7)
		{
			//printindelstatistics
			//myinclude.at(25)=true;

			if(s2=="TRUE") printrates=true; else if(s2=="FALSE") printrates=false; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==8)
		{

			if(s2=="TRUE") insertaslowercase=true; else if(s2=="FALSE") insertaslowercase=false; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==9)
		{
			// command is "[printcodonsasaminoacids]"
			if(s2=="TRUE") printcodonsasDNA=false; else if(s2=="FALSE") printcodonsasDNA=true; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==10)
		{
			// command is "[fileperrep]"
			if(s2=="TRUE") fileperrep=true; else if(s2=="FALSE") fileperrep=false; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==11)
		{
			// command is "[markdeletedinsertions]"
			if(s2=="TRUE") markdeletedinsertions=true; else if(s2=="FALSE") markdeletedinsertions=false; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==12)
		{
			// command is "[fixproportions]"
			if(s2=="TRUE") fixtrueproportions=true; else if(s2=="FALSE") fixtrueproportions=false; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==13)
		{
			//printallrates


			if(s2=="TRUE") printallrates=true; else if(s2=="FALSE") printallrates=false; else
			{controlerrorprint2("[SETTINGS]","",scommands[commanderror],"Value must be TRUE or FALSE for this command.",""); {if(breakonerror) return -1;} }
		}
		else if(commanderror==14)
		{
			//randomseed

			globalseed=true;
			if(!allAinB(s2,"1234567890-")) { controlerrorprint2("[SETTINGS]", "",scommands[commanderror], "This value must be an integer.",""); {if(breakonerror) return -1;} }

			idum2=atoi(s2.c_str()); if(idum2<0) idum2=-idum2; 	
		}

	}

	if(globalseed) {idum=idum2; mtrand1.seed(idum);	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////


int dealwithmodel(vector<string> &block)		
{
	// this is a big function that deals with setting up of substitution and indel models.  What is not tested or pre-processed here is passed to the model class constructor.
	
	vector<double> insertrates, deleterates;
	if(!oldmethod) {insertrates.push_back(0); deleterates.push_back(0);}
	else 
	{
		int tempsize; 
		if(type==1) tempsize=4; else if(type==2) tempsize=20; else if(type==3) tempsize=64; else cout<<"tempsize and type error in dealwithmodel"<<endl;
		for(int jp=0; jp<tempsize; jp++) {insertrates.push_back(0); deleterates.push_back(0);}
	}
 
	vector<double> aamodel;
	string name="DEFAULT";
	int modelnumber=0;
	int geneticcode=1;
	bool copiedmodel=false;
	double insertrate=0;
	double deleterate=0;
	indelmodel insertmodel=indelmodel();
	indelmodel deletemodel=indelmodel();
	double alpha=0;
	double pinv=0;
	int ngamcat=5;
	double codonrates[3]; codonrates[0]=codonrates[1]=codonrates[2]=1;
	vector<double> basefreqs, rootbasefreqs, insertfreqs, params;

	string commandsarray[143]={"[insertmodel]","[deletemodel]","[indelmodel]","[insertrate]","[deleterate]","[indelrate]",
		"[geneticcode]", "[submodel]","[rates]","[Crates]","[statefreq]","[istatefreq]","[rootstatefreq]","[XXXrootseq]",
		"[testindeldistributions]"};  // acceptable list of commands
	
	// what about tree creation etc?

	vector<string> commands(commandsarray,commandsarray+13);
	string mymodelname;
			
//	for(int i=0; i<block.size(); i++) cout<<"model "<<block.at(i)<<endl;


	string hg=block.at(0), hg1;
//	if(hg[0]!='#' || hg[hg.size()-1]!='#') { controlerrorprint2("[MODEL]", "?", "", "First statement in a [MODEL] block must be a model name statement in the form #modelname#",""); {if(breakonerror) return -1;} }
	if(!allAinB(hg,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[MODEL]", hg, "", "First statement in a [MODEL] block must be a model name statement.\nThe name should only contain  ,.;()-_ and alpha-numeric characters.",hg); {if(breakonerror) return -1;} }
	else 
	{
		for(int i0=0; i0<hg.size(); i0++) if(hg[i0]!=' ') mymodelname+=hg[i0];
	//	cout<<"MODEL NAME "<<mymodelname<<endl;

		if(checkthere(mymodelname,totalmodelnames)!=-1)
		{ controlerrorprint2("[MODEL]", mymodelname, "", "A model with this name has been specified twice.",""); {if(breakonerror) return -1;} }
		else {
				totalmodelnames.push_back(mymodelname); 
				name=mymodelname; 
			}
	
		vector<string> tempvec;
		int blocksize=block.size()-1;
		
		//bool waitforcommand=false;
		
		hg="";
		int lasttest=-999;
		int mytest=-999;
		int astart=1;
		if(blocksize+1>astart) hg=block.at(astart);
		if(hg[0]=='#')
		{ 
			// this if bracket is concerned with automatically copying from model to model.
			// feature will be disabled in final program - so has not been checked after Dec 2008

			mymodelname="";
			for(int i0=1; i0<hg.size(); i0++) if(hg[i0]!=' ') mymodelname+=hg[i0];
		//	cout<<"MODEL NAME 2"<<mymodelname<<endl;
	
			if(mymodelname==name) {controlerrorprint2("[MODEL]", mymodelname, "", "You cannot copy a model within it's own definition statement!\nThe second model name should be different from "+mymodelname,""); {if(breakonerror) return -1;} }
	
			int ans=checkthere(mymodelname,totalmodelnames);
			if(ans==-1) {controlerrorprint2("[MODEL]", name, "", "Cannot copy model "+mymodelname+".\nNo model with this name has been defined yet.",""); {if(breakonerror) return -1;} }
			else
			{
				model* m=&(totalmodels.at(ans));

				modelnumber=(*m).modelnumber;
				geneticcode=(*m).geneticcode;
				copiedmodel=true;
				insertrate=(*m).insertrate;
				deleterate=(*m).deleterate;
				
				//insertrates=(*m).insertrates;
				//deleterates=(*m).deleterates;

				//inlength=(*m).inlength;
				//dellength=(*m).dellength;
				//rootlength=(*m).rootlength;
				alpha=(*m).alpha;
				pinv=(*m).pinv;
				ngamcat=(*m).ngamcat;
				codonrates[0]=(*m).codonrates[0];
				codonrates[1]=(*m).codonrates[1];
				codonrates[2]=(*m).codonrates[2];

				//	vector<double> params;  // ??

				rootbasefreqs=(*m).rootbasefreqs;
				basefreqs=(*m).basefreqs;
				insertfreqs=(*m).insertfreqs;
				astart++;
				//	vector<double> ratevec; // ??
			}

		}

		string oldhg="Z";
		for(int i1=astart; i1<blocksize+1; i1++) 
		{
			hg=block.at(i1); 
			lasttest=mytest;
			string myallowed="0123456789.";
			
			if(oldhg=="[insertmodel]" || oldhg=="[deletemodel]" || oldhg=="[indelmodel]"  ) myallowed="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890."; //"NBUSEROLDPOWCUMDAWGLAV";

			mytest=teststring(lasttest, hg, commands, "[MODEL]",myallowed,mymodelname);
			
			oldhg=hg;

			if(mytest==-2) {if(breakonerror) return -1;} 
			else 
			{//a1
			 //			cout<<"H2 "<<hg<<" "<<mytest<<" "<<lasttest<<"  "<<blocksize<<"  "<<astart<<endl;
				
			if(i1==blocksize )  {tempvec.push_back(hg); }//cout<<"H1 "<<hg<<endl;}
				
			if(mytest==lasttest && i1!=blocksize )  {tempvec.push_back(hg); }// cout<<"H1 "<<hg<<endl;}
				else if(i1!=astart)
				{//a2

					if(tempvec.size()==0) 
					{
						controlerrorprint2("[MODEL]", name, commands.at(lasttest), "No values were found after this command in the control file.","");
						{if(breakonerror) return -1;} 
					}
					else if(lasttest>-1 && lasttest<3)
					{
						// indel models
		
						// lasttest = 0  for insertion length model
						//            1  for deletion  length model
						//            2 for same length model for both insertions and deletions

						bool noindelerror=true;
						string myin=tempvec.at(0);
						int tempsize=tempvec.size();

						int myindelmodelint;
						double mr, mM;
						double mq, ma, mb;

						double entryd=0, totald=0, meand=0;
						int thecount=0;
						
						vector<double> usermodel;
	
						// 0 USER
						// 1 NB
						// 2 POW
						// 3 LAV

						// 11 NBOLD - NB calculated using p GEO
						// 12 POWCUM - power law calculated using cumulative method used for USER and LAV
						// 13 POWDAWG - power law calculated using zipf function from DAWG

						// User Model
						if(myin=="USER") 
						{
							myindelmodelint=0;

							if(tempsize!=2) noindelerror=false;
														
							if(noindelerror)
							{
								myin=tempvec.at(1); 

								ifstream if1; if1.open(myin.c_str());
								
								double countcount=0;
								if(!(if1.good()))
								{
									controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Could not find a file named:\n"+myin,"");
									mytest=-2; {if(breakonerror) return -1;}
								}
								char c=if1.get();
								string s;

								bool writeon=false;
								while(if1.good())
								{ 
									if(!AinB(c,"0123456789. \n\t\r")) 
									{
										string y; y+=c;
										controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Illegal character "+y+" found in file named:\n"+myin+"\nAllowed characters are 0.123456789\nSpaces, horizontal tabs, and new lines are also allowed","");
										mytest=-2; {if(breakonerror) return -1;}
									}
									if(AinB(c,"0123456789.")) writeon=true; 

									if(writeon && !AinB(c,"0123456789.")) 
									{
										countcount++; 
										writeon=false; 
										entryd=atof(s.c_str());
										meand+=(countcount*entryd); 
										totald+=entryd; 
										usermodel.push_back(totald); 
										s=""; 
									}

									if(writeon) s+=c;

									c=if1.get();
								}



								if(s!="") 
								{
									countcount++; 
									entryd=atof(s.c_str()); 
									meand+=(countcount*entryd); 
									totald+=entryd; 
									usermodel.push_back(totald); 
								}
								
								if(usermodel.empty()) 							
								{
									controlerrorprint2("[MODEL]", name, commands.at(lasttest), "The user-defined indel length model must contain at least one value.!","");
									mytest=-2; {if(breakonerror) return -1;}
								}

								double diff=totald-1;
								if(diff<0) diff=-diff;
								if(diff>0.00001) 
								{
									stringstream df; df<<totald; string ff=df.str();
									controlerrorprint2("[MODEL]", name, commands.at(lasttest), "User-defined indel model does not sum to 1 and so has been rescaled.\n(It summed to "+ff+") ","");

									for(int i=0; i<usermodel.size(); i++) {(usermodel.at(i))/=totald; }
								}
	
								meand/=totald;
							}

						} // end of User Model

						// Negative Binomial Distribution
						else if(myin=="NB" || myin=="NBOLD") 
						{
							if(myin=="NB")	myindelmodelint=1;	// used in program
							else			myindelmodelint=11;	// alternative parameterisation for testing
										
							if(tempsize!=3) noindelerror=false;
														
							if(noindelerror)
							{
								myin=tempvec.at(1); 
								if(!allAinB(myin,"0123456789.")) noindelerror=false;
								else mq=atof(myin.c_str());
								if(mq<0 || mq>1) noindelerror=false;
							}
								
							if(noindelerror)
							{
								myin=tempvec.at(2); 
								if(!allAinB(myin,"0123456789")) noindelerror=false;
								else mr=atoi(myin.c_str());
								if(mr==0) noindelerror=false;
							}

							if(!noindelerror) 							
							{	
								controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Expecting 2 values after \"NB\" when using Negative Binomial Model.\nThe first should be a proportion 0 <= q <= 1.\nThe second should be an integer r > 0.","");
								mytest=-2; {if(breakonerror) return -1;} 
							}

							meand =  1.0+mr*mq/(1.0-mq);

						} // end of Negative Binomial Distribution


						// Zipfian Distribution (Power Law)
						else if(myin=="POW" || myin=="POWDAWG" || myin=="POWCUM") 
						{
							if(myin=="POW")     myindelmodelint=2;  // Fast Zipf

							if(myin=="POWDAWG") myindelmodelint=13; // Slower Zipf from DAWG - here for checking fast zipf!

							if(myin=="POWCUM")  myindelmodelint=12; // cumulative method as for USER and LAV used to check both above!

							if(tempsize!=2 && tempsize!=3) noindelerror=false;
														
							if(noindelerror)
							{
								myin=tempvec.at(1); 
								if(!allAinB(myin,"0123456789.")) noindelerror=false;
								else ma=atof(myin.c_str());
								if(ma<=1) noindelerror=false;
							}

							if(tempsize==3 && noindelerror)
							{
								myin=tempvec.at(2); 
								if(!allAinB(myin,"0123456789")) noindelerror=false;
								else mM=atoi(myin.c_str());
								if(mM<=1) noindelerror=false;
							}
							else mM=2147483647; // largest integer value possible. //mM=pow(10,12); //mM=~0-1; // largest known genome is 132 pg  (1pg is ~ 978 MB, so 132pg is ~ 129 GB, i.e. 10^12 is bigger than any genome!)
							
							if(!noindelerror) 							
							{	
								controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Expecting 1 or 2 values after \"POW\" when using Zipfian Power Model.\nThis obligatory first value should be a decimal number a > 1.\nThe second optional value should be an integer M > 1.","");
								mytest=-2; {if(breakonerror) return -1;} 
							}

							double blah;
							if(myindelmodelint==12)
							{
								for(double u=1; u<mM+1; u++) {entryd=pow(u,-1*ma); totald+=entryd; meand+=( entryd * u );usermodel.push_back(totald); }
		
								for(int i=0; i<usermodel.size(); i++) {(usermodel.at(i))/=totald;  }
							}
							else
							{
								//cout<<"mM is "<<mM<<endl;
								double last=2, diff, bit=pow(10,-15);
								for(double u=1; u<mM+1; u++) {entryd=pow(u,-1*ma); meand+=(entryd*u); totald+=entryd; diff=last-entryd; last=entryd;  if(diff<bit) {blah=u; break;}}	
							}
							//cout<<"BLAH  "<<blah<<endl;
							meand/=totald;
							//cout<<meand<<"  "<<totald<<endl;
	
						} // end of Zipfian Distribution (Power Law)

						// Lavalette Distribution
						else  if(myin=="LAV") 
						{					
							/*
							Lavalette Distribution References

							Lavalette distribution (Lavalette, 1996; Popescu et al., 1997; Popescu, 2003).  
							Lavalette D. (1996) Facteur d’impact: impartialité ou impuissance?, Internal Report, INSERM U350, Institut Curie - Recherche, Bât. 112, Centre Universitaire, 91405 Orsay, France (November 1996)
							Popescu, I.I. (1997) On the Lavalette Ranking Law (with M. Ganciu, M. C. Penache, and D. Penache), Romanian Reports in Physics, 49, 3-27
							Popescu, I.I. (2003) On a Zipf’s law extension to impact factors. Glottometrics, 6. 83-93.
							*/

							myindelmodelint=3;

							if(tempsize!=3) noindelerror=false;
														
							if(noindelerror)
							{
								myin=tempvec.at(1); 
								if(!allAinB(myin,"0123456789.")) noindelerror=false;
								else mb=atof(myin.c_str());
								if(mb<=1) noindelerror=false;
							}

							if(noindelerror)
							{
								myin=tempvec.at(2); 
								if(!allAinB(myin,"0123456789")) noindelerror=false;
								else mM=atoi(myin.c_str());
								if(mM<=1) noindelerror=false;
							}

							if(!noindelerror) 							
							{	
								controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Expecting 2 values after \"LAV\" when using Lavalette Model.\nThe first should be a decimal b > 1.\nThe second (maximum length) should be an integer M>1.","");
								mytest=-2; {if(breakonerror) return -1;} 
							}

							
							for(double u=1; u<mM+1; u++) {entryd=pow(mM*u/(mM-u+1),-1*mb); meand+=(u*entryd); totald+=entryd; usermodel.push_back(totald); }
		
							for(int i=0; i<usermodel.size(); i++) {(usermodel.at(i))/=totald; }
							meand/=totald;


						} // end of Lavalette Distribution

						else
						{	
							controlerrorprint2("[MODEL]", name, commands.at(lasttest), "The first value in this command specifies the insertion or deletion model.\nPossible values are:\n\n     NB = Negative Binomial Distribution\n    POW = Zipfian Distribution (Power Law)\n    LAV = Lavalette Distribution\n   USER = explicit, user-defined indel length distribution\n",myin);
							mytest=-2; {if(breakonerror) return -1;} 
						}
		
						if(lasttest==0 || lasttest==2) 
						{
							insertmodel.type=myindelmodelint;
							insertmodel.r=int(mr);
							insertmodel.q=mq;
							insertmodel.a=ma;
							insertmodel.b=mb;
							insertmodel.M=int(mM);			//cout<<"mM before "<<mM<<"  and after "<<int(mM)<<endl;
							insertmodel.meansize=meand;
							insertmodel.usermodel=usermodel;
						}
						if(lasttest==1 || lasttest==2) 
						{
							deletemodel.type=myindelmodelint;
							deletemodel.r=int(mr);
							deletemodel.q=mq;
							deletemodel.a=ma;
							deletemodel.b=mb;
							deletemodel.M=int(mM);			//cout<<"mM before "<<mM<<"  and after "<<int(mM)<<endl;
							deletemodel.meansize=meand;
							deletemodel.usermodel=usermodel;
						}
						
					} // end of indel model stuff

					else if(lasttest>2 && lasttest<7) 
					{
						// last test = 3  for insertion rate
						// last test = 4  for deletion rate
						// last test = 5  for insertion AND deletion rate
						// last test = 6  for genetic code

						if(tempvec.size()>1) 
						{	
							controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Expecting 1 value after this command but found "+tempvec.at(0)+" then "+tempvec.at(1)+".","");
							mytest=-2; {if(breakonerror) return -1;} 
						}

						string myin=tempvec.at(0);
			
						if(lasttest==6)
						{
							// genetic code
							
							if(type!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "This command only has effect in CODON simulations.","");   {if(breakonerror) return -1;} }  

							if(AinB('.',myin)) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Value for this command should not be decimal, was expecting integer.",myin);   {if(breakonerror) return -1;} }  

							else
							{
								int myout=atoi(myin.c_str());
									
								if(  (myout>0 && myout<7) || (myout>8 && myout<17) || (myout>20 && myout<24)  ) geneticcode=myout; 
								
								else {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Value for this command should be an integer 1-6, 9-16, or 21-23\nThese correspond to the entries on GenBank (October 2008).",myin);   {if(breakonerror) return -1;} }
							}
						}
						else 
						{ 
							double myout=atof(myin.c_str());

							if(myout<0) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Value for this command should be greater than or equal to zero.",myin);   {if(breakonerror) return -1;} }  

							if(lasttest==3) insertrate=myout; 
							
							else if(lasttest==4) deleterate=myout;

							else insertrate=deleterate=myout;
						}
							
					} 
					else if(lasttest==7)
					{
						// lasttest == 7 is the substitution model  -----> big if bracket

						string myin=tempvec.at(0);

						if(!allAinB(myin,"0123456789.") )
						{
							ifstream ig1; ig1.open(myin.c_str()); 

							if(type==2)
							{
								if(!(ig1.good())) {controlerrorprint2("[MODEL]", name, "[submodel]", myin+"\nis not a model name or number and no file of this name exists.\nFor protein models, the entry after a [submodel] command must be\nan integer or a filename.",""); {if(breakonerror) return -1;} }
								
								modelnumber=-1;
								char c2; string sss;

								c2=ig1.get(); aamodel.clear(); 
								bool writeon=false;
								while(ig1.good())
								{	
									if(AinB(c2,"1234567890. \t\n\r"))
									{
										if(AinB(c2," \t\n\r")) writeon=false; else writeon=true;
										
										if(writeon) sss+=c2;
										else if(sss!="") {aamodel.push_back(  atof(sss.c_str())  ); sss="";}
									}
									else {controlerrorprint2("[MODEL]", name, "[submodel]", "Reading user defined protein substitution model from file:\n"+myin+"\n...but have found illegal character \""+c2+"\"\nEntries in this file can only be decimal numbers separated by\nspaces, horizontal tabs and new lines.",""); {if(breakonerror) return -1;} }
								
									c2=ig1.get(); 
								}
								int wrongsize2=aamodel.size(); stringstream dd; dd<<wrongsize2; string wrongsize=dd.str();

								if(wrongsize2!=210) {controlerrorprint2("[MODEL]", name, "[submodel]", "Reading user defined protein substitution model from file:\n"+myin+"\n...but have found "+wrongsize+" entries. Was expecting 210 entries.\nThe first 190 entries are from the substitution matrix.\nThe final 20 values are for the stationary amino-acid frequencies.\nThe substititution model must be symmetric and the file must only contain\nthe lower triangular half of the substitution matrix.",""); {if(breakonerror) return -1;} }
							}
							else if(type==1) {controlerrorprint2("[MODEL]", name, "[submodel]", myin+" is not an integer model number.\nThe entry after a [submodel] command must be an integer.",""); {if(breakonerror) return -1;} }
		
						}
						else
						{

							int mymodel=atoi(myin.c_str());

							int b=17; string typestring="NUCLEOTIDE"; string bnum="16";  //for type 1
							if(type==2) {typestring="AMINOACID";}
							else if(type==3) {b=16; bnum="15"; typestring="CODON";}

							if(type!=3)
							{
								if(mymodel>-1 && mymodel<b && !AinB('.',myin)) {modelnumber=mymodel;} 
								else {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "First value for this command represents the model number.\nThis should be a integer between 0 and "+bnum+" inclusive\nwhen type is set to "+typestring+".",myin);   {if(breakonerror) return -1;} }
							}
							else
							{
								if(mymodel!=14 && mymodel!=15) modelnumber=mymodel=3;
							}
							int nstsize=tempvec.size();
							
							// old way
							//for(int hy=1; hy<nstsize; hy++) {params.push_back(atof((tempvec.at(hy)).c_str())); }    //cout<<"WER "<<atof((tempvec.at(hy)).c_str())<<endl;}

							if(type==3) for(int hy=0; hy<nstsize; hy++) {params.push_back(atof((tempvec.at(hy)).c_str())); }    //cout<<"WER "<<atof((tempvec.at(hy)).c_str())<<endl;}
							else        for(int hy=1; hy<nstsize; hy++) {params.push_back(atof((tempvec.at(hy)).c_str())); }    //cout<<"WER "<<atof((tempvec.at(hy)).c_str())<<endl;}

							if(type==2) 
							{
								if(nstsize>1) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Value for this command should be a single integer between 0 and "+bnum+"\nor a model/file name when type is set to AMINOACID.\nBut found "+tempvec.at(0)+" then "+tempvec.at(1)+".","");   {if(breakonerror) return -1;} }


							}
							
							else if(type==1)
							{
								stringstream fr; fr<<nstsize-1; string fh=fr.str();

	 								 if(mymodel==0  && nstsize!=1)                 {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as JC69. \nNot expecting any substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==1  && nstsize!=1)                 {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as F81. \nNot expecting any substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==2  && nstsize!=2 && nstsize!=3)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as K80. \nExpecting 1 or 2 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==3  && nstsize!=2 && nstsize!=3)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as HKY85. \nExpecting 1 or 2 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==4  && nstsize!=3 && nstsize!=4)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as TN93ef. \nExpecting 2 or 3 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==5  && nstsize!=3 && nstsize!=4)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as TN93. \nExpecting 2 or 3 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==6  && nstsize!=3 && nstsize!=4)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as K81. \nExpecting 2 or 3 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==7  && nstsize!=3 && nstsize!=4)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as K81uf. \nExpecting 2 or 3 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==8  && nstsize!=4 && nstsize!=5)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as TIMef. \nExpecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==9  && nstsize!=4 && nstsize!=5)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as TIM. \nExpecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==10 && nstsize!=5 && nstsize!=6)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as TVMef. \nExpecting 4 or 5 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==11 && nstsize!=5 && nstsize!=6)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as TVM. \nExpecting 4 or 5 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==12 && nstsize!=6 && nstsize!=7)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as SYM. \nExpecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==13 && nstsize!=6 && nstsize!=7)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as GTR. \nExpecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								
								else if(mymodel==14 && nstsize!=2 && nstsize!=3)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as F84ef. \nExpecting 1 or 2 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==15 && nstsize!=2 && nstsize!=3)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as F84. \nExpecting 1 or 2 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
							//	else if(mymodel==16 && nstsize!=3 && nstsize!=4)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as T92. \nExpecting 2 or 3 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
							//	else if(mymodel==17 && nstsize!=12 && nstsize!=13) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as UNREST. \nExpecting 11 or 12 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==16 && nstsize!=12 && nstsize!=13) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as UNREST. \nExpecting 11 or 12 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
							}

							else if(type==3)
							{
								ngamcat=1;
								//nstsize--;
								stringstream fr; fr<<nstsize; string fh=fr.str();

	 								 if(mymodel==0  && nstsize!=2) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M0. \nWas expecting 2 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==1  && nstsize!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M1. \nWas expecting 3 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==2  && nstsize!=5) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M2. \nWas expecting 5 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								

								/*
								else if(mymodel==5 ) 
								{
									if(nstsize==4) 
									{
										if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 4 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} }
										else ngamcat=atoi((tempvec.back()).c_str());
									}
									else 
									{
										if(nstsize!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } 
										else ngamcat=0;
									}
								}

								
								else if(mymodel==5 ) {if(nstsize==4){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 4 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==6 ) {if(nstsize==6){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 6 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=5) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==7 ) {if(nstsize==4){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 4 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==8 ) {if(nstsize==6){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 6 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=5) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==9 ) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==10) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==11) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==12) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								else if(mymodel==13) {if(nstsize==8){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWhen specifying 8 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=7) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 7 or 8 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=0;}}
								/*
								else if(mymodel==5  && nstsize!=3 && nstsize!=4) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==6  && nstsize!=5 && nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M6. \nWas expecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==7  && nstsize!=3 && nstsize!=4) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M7. \nWas expecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==8  && nstsize!=5 && nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M8. \nWas expecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==9  && nstsize!=6 && nstsize!=7) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M9. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==10 && nstsize!=6 && nstsize!=7) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M10.\nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==11 && nstsize!=6 && nstsize!=7) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M11.\nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==12 && nstsize!=6 && nstsize!=7) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M12.\nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==13 && nstsize!=7 && nstsize!=8) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M13.\nWas expecting 7 or 8 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								*/

								else if(mymodel==5 ) {if(nstsize==4){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5.\nWhen specifying 4 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==6 ) {if(nstsize==6){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M6.\nWhen specifying 6 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=5) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==7 ) {if(nstsize==4){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M7.\nWhen specifying 4 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 3 or 4 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==8 ) {if(nstsize==6){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M8.\nWhen specifying 6 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=5) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 5 or 6 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==9 ) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M9.\nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==10) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M10. \nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==11) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M11. \nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==12) {if(nstsize==7){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M12. \nWhen specifying 7 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=6) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 6 or 7 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}
								else if(mymodel==13) {if(nstsize==8){ if(AinB('.',tempvec.back())) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M13. \nWhen specifying 8 substitution parameters after the model number the\nlast value represents the number of categories to use in the discrete \napproximation for any continuous distributions.  This must be an integer. ",tempvec.back()); {if(breakonerror) return -1;} } else {ngamcat=atoi((tempvec.back()).c_str());params.pop_back();}}  else {if(nstsize!=7) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M5. \nWas expecting 7 or 8 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} } else ngamcat=10;}}



								else if(mymodel==3  && (nstsize%2)!=0) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M3. \nWas expecting 2K substitution parameters for K categories.\nOrder should be: kappa p_0  p_1  ... p_(k-2)  w_0  w_1  ... w_(k-1)\nFor M0 this would be: kappa w_0\nFor M1a this would be: kappa p_0 w_0 w_1.\nFor M2a this would be: kappa p_0 p_1 w_0 w_1 w_2 etc.\nInstead found an odd number of parameters ("+fh+") after model number",""); {if(breakonerror) return -1;} }


								//else if(mymodel==3  && (nstsize%2)!=0) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M3. \nWas expecting 2K-1 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==4  && nstsize<3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as M4. \nWas expecting at least 2 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								
								else if(mymodel==14 && nstsize!=0) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as ECMrest. \nNot expecting any substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }
								else if(mymodel==15 && nstsize!=0) {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "CODON substitution model set as ECMunrest. \nNot expecting any substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }


								// why is this here?  it says NUCLEOTIDE - is this the same as before? - maybe just used for copying?
								//else if(mymodel==2  && nstsize!=2 && nstsize!=3)   {controlerrorprint2("[MODEL]", name, commands.at(lasttest),  "NUCLEOTIDE substitution model set as K80. \nExpecting 1 or 2 substitution parameters.\nInstead found "+fh+" parameters after model number",""); {if(breakonerror) return -1;} }

								// cout<<"NGAMCAT "<<ngamcat<<endl;
		
								if(mymodel==12 || mymodel==13 || mymodel==2)
								{
									if(params.at(1)>1) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "proportion p0 must be between 0 and 1.  You entered "+tempvec.at(2),""); {if(breakonerror) return -1;} }
									if(params.at(2)>1) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "proportion p1 must be between 0 and 1.  You entered "+tempvec.at(3),""); {if(breakonerror) return -1;} }
									if(mymodel!=12 && params.at(1)+params.at(2)>1) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "sum of proportions p0 and p1 must be between 0 and 1.  \nYou entered "+tempvec.at(2)+" and "+tempvec.at(3),""); {if(breakonerror) return -1;} }

								}
								else if( (mymodel<12 && mymodel>7) || mymodel==6 || mymodel==1)
								{
									if(params.at(1)>1) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "proportion p0 must be between 0 and 1.  You entered "+tempvec.at(2),""); {if(breakonerror) return -1;} }
								}
								else if(mymodel==3 || mymodel==4)
								{
									int thesize=params.size(); 
									if(mymodel==3) thesize=(thesize-1)/2 +1;
									double sum=0;
									for(int gh8=1; gh8<thesize; gh8++) 
									{
										sum+=params.at(gh8);
										if(params.at(gh8)>1) 
										{
											stringstream cd; cd<<gh8-1; string hf=cd.str();
											controlerrorprint2("[MODEL]", name, commands.at(lasttest), "proportion p"+hf+" must be between 0 and 1.  You entered "+tempvec.at(gh8+1),""); 
											{if(breakonerror) return -1;} 
										}			
									}
									if(sum>1) 
									{
										stringstream cd; cd<<sum; string hf=cd.str();
											
										controlerrorprint2("[MODEL]", name, commands.at(lasttest), "sum of proportions must be between 0 and 1.  It was "+hf+".",""); {if(breakonerror) return -1;} 
									}
									
								}
							
							} // end of else if type =3 bracket							

						}
					}// end of lasttest == 7 bracket 	
					else if(lasttest==99)
					{//y1
						// params
						
///////////////////////////////////////////////

						//  disabled bracket now


					}//y1
					else if(lasttest==8)
					{
						// lasttest == 8 sets site specific rates for.   order of commands is --------> rates  pinv alpha ngamcat

						if(type==3)	{controlerrorprint2("[MODEL]", name, commands.at(lasttest), "This command has no effect in a CODON simulation.","");   {if(breakonerror) return -1;} }
						else if(tempvec.size()!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "There should be 3 values following a [rates] command: pinv alpha ngamcat.","");  {if(breakonerror) return -1;} }
						else
						{
							if(AinB('.',tempvec.at(2))) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "The 3rd value representing ngamcat should be an integer",tempvec.at(2));  {if(breakonerror) return -1;} }

							pinv=atof((tempvec.at(0)).c_str());
							alpha=atof((tempvec.at(1)).c_str());
							ngamcat=atoi((tempvec.at(2)).c_str());

							if(pinv>1) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "The 1st value representing pinv should be between 0 and 1.",tempvec.at(0));  {if(breakonerror) return -1;} }
							
						}
					}
					else if(lasttest==9)
					{
						// Crates  ---->    now disabled, was for fixed relative rates for each codon position.

						if(tempvec.size()!=3) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "There should be 3 decimal values following a [Crates] command.","");  {if(breakonerror) return -1;} }
						else
						{
							codonrates[0]=atof((tempvec.at(0)).c_str());
							codonrates[1]=atof((tempvec.at(1)).c_str());
							codonrates[2]=atoi((tempvec.at(2)).c_str());
						}
					}
					else if(lasttest>9)
					{//x1
						// lasttest == 10 is base frequencies being defined
						// lasttest == 11 is insertion base frequenciess being defined (i.e. diff. base freqs being used for creating inserted sequences and for substitution model)
						// lasttest == 12 is root base frequencies being defined  (i.e.  using diff base freqs for creating root sequence than those used by substitution model)

						// base frequency stuff

						int mysize=tempvec.size();
				
						stringstream fd; fd<<mysize; string mysizeS=fd.str();
						bool errorprint=false;
						string mytype, mysizes;

						if(type==3) {if(mysize!=4 && mysize!=12 && mysize!=64) {mytype="CODON"; mysizes="4, 12, or 64"; errorprint=true; }}
						else if(type==2) {if(mysize!=20){mytype="AMINOACID"; mysizes="20"; errorprint=true; }}
						else if(type==1) {if(mysize!=4) {mytype="NUCLEOTIDE";     mysizes="4"; errorprint=true; }}

						if(errorprint) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "There should be "+mysizes+" values for a base frequency command\nwhen simulating "+mytype+" data. You entered "+mysizeS+" values.","");  {if(breakonerror) return -1;} }
						else
						{ //x2
							vector<double> tempbases, tempbases2; 
							for(int kh=0; kh<mysize; kh++) 
							{
								double mytemp=atof( ( tempvec.at(kh) ).c_str() );
								if(mytemp<0) {controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Values for base frequency commands must be positive.","");  {if(breakonerror) return -1;} }
								else tempbases.push_back(mytemp);
							}

							double scaler=0,lastscaler,diff;

							if(type==3 && mysize==4)
							{
								// F1X4 frequencies
								scaler=0;
								for(int gh1q=0; gh1q<mysize; gh1q++) scaler+=tempbases.at(gh1q); 
								for(int gh2q=0; gh2q<mysize; gh2q++) tempbases.at(gh2q)/=scaler;

								diff=scaler-1; if(diff<0) diff=-diff;
								if(diff>0.0001) controlerrorprint2("[MODEL]", name, commands.at(lasttest), "F1X4 base frequencies did not add up to 1 so they have been rescaled.",""); 

								//for(int i1=0; i1<4; i1++) for(int j1=0; j1<4; j1++) for(int k1=0; k1<4; k1++)
								//tempbases2.push_back(tempbases.at(i1)*tempbases.at(j1)*tempbases.at(k1));
								
								//enforcestops(geneticcode,tempbases2);
								//mysize=64;

								tempbases2=tempbases;
								
							}
							else if(type==3 && mysize==12) 
							{
								// F3X4 frequencies
								scaler=0;
								for(int gh1a=0; gh1a<4; gh1a++) scaler+=tempbases.at(gh1a); 
								for(int gh2a=0; gh2a<4; gh2a++) tempbases.at(gh2a)/=scaler;

								diff=scaler-1; if(diff<0) diff=-diff;
								if(diff>0.0001) controlerrorprint2("[MODEL]", name, commands.at(lasttest), "F3X4 base frequencies for the 1st codon position\ndid not add up to 1 so they have been rescaled.",""); 

								scaler=0;
								for(int gh1b=4; gh1b<8; gh1b++) scaler+=tempbases.at(gh1b); 
								for(int gh2b=4; gh2b<8; gh2b++) tempbases.at(gh2b)/=scaler;

								diff=scaler-1; if(diff<0) diff=-diff;
								if(diff>0.0001) controlerrorprint2("[MODEL]", name, commands.at(lasttest), "F3X4 base frequencies for the 2nd codon position\ndid not add up to 1 so they have been rescaled.",""); 

								scaler=0;
								for(int gh1c=8; gh1c<12; gh1c++) scaler+=tempbases.at(gh1c); 
								for(int gh2c=8; gh2c<12; gh2c++) tempbases.at(gh2c)/=scaler;

								diff=scaler-1; if(diff<0) diff=-diff;
								if(diff>0.0001) controlerrorprint2("[MODEL]", name, commands.at(lasttest), "F3X4 base frequencies for the 3rd codon position\ndid not add up to 1 so they have been rescaled.",""); 


								//for(int i1=0; i1<4; i1++) for(int j1=4; j1<8; j1++) for(int k1=8; k1<12; k1++)
								//tempbases2.push_back(tempbases.at(i1)*tempbases.at(j1)*tempbases.at(k1));

								//enforcestops(geneticcode,tempbases2);
								//mysize=64;

								tempbases2=tempbases;

							}
							else if(type==3 && mysize==64)
							{
								//Fcodon frequencies
								scaler=0;
								for(int gh1s=0; gh1s<mysize; gh1s++) scaler+=tempbases.at(gh1s); 
								for(int gh2s=0; gh2s<mysize; gh2s++) tempbases.at(gh2s)/=scaler;

								lastscaler=scaler;
								diff=scaler-1; if(diff<0) diff=-diff;
								if(diff>0.0001) controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Base frequencies did not add up to 1 so they have been rescaled.",""); 
								
								tempbases2=tempbases;
							
							}
							else tempbases2=tempbases;
		
							if(type==1 || type==2)
							{
								scaler=0;
								for(int gh1=0; gh1<mysize; gh1++) scaler+=tempbases2.at(gh1); 
								for(int gh2=0; gh2<mysize; gh2++) tempbases2.at(gh2)/=scaler;

								diff=scaler-1; if(diff<0) diff=-diff;
								if(diff>0.0001) controlerrorprint2("[MODEL]", name, commands.at(lasttest), "Base frequencies did not add up to 1 so they have been rescaled.",""); 
							}
					
							if(lasttest==10) 
							{
								basefreqs=tempbases2;
							}
							else if(lasttest==11) insertfreqs=tempbases2;
							else if(lasttest==12) rootbasefreqs=tempbases2;
						
						}//x2
					}//x1

					tempvec.clear();

				} //end of a2  bracket (long way above!)	

			} //end of a1  bracket ( long long way above!)

		}	// end of for bracket

	} // end of model name else bracket


		totalmodels.push_back(model(totalmodels.size(), type, name, modelnumber, geneticcode, copiedmodel, insertrate, deleterate,  alpha, pinv, ngamcat, codonrates, basefreqs, rootbasefreqs, insertfreqs, params, aamodel, insertmodel, deletemodel ));
		
		if((totalmodels.back()).error==-1) return -1;   // returns error if there was an error resulting from the model constructor in the model class
	
	return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int dealwithsites(vector<string> &block)		
{
	// this function is not used any more as the sites class is disabled now.

	int mymy=block.size();
	string name;
	if(mymy==0) {controlerrorprint2("[SITES]", "", "", "You must specify a name for a sites block. This one is empty.",""); {if(breakonerror) return -1;}  }
	else
	{
		name=block.at(0);
		if(!allAinB(name,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[SITES]", name, "", "First statement in a [SITES] block must be a sites name statement.\nThe name should only contain  ,.;()-_ and alpha-numeric characters.",""); {if(breakonerror) return -1;} }		

		if(checkthere(name,totalsitenames)!=-1)
		{ controlerrorprint2("[SITES]", name, "", "A sites class with this name has been specified twice.",""); {if(breakonerror) return -1;} }
	
	}

	if(mymy==1) {controlerrorprint2("[SITES]", name, "", "This sites block contains no commands.",""); {if(breakonerror) return -1;}  }
	else name=block.at(0);

	vector<double> props;
	vector<string> names;

	if(mymy==2)
	{
		string a=block.at(1),b;
		int commacount=0;
		bool lastcomma=false;
		bool modelnameempty=false;
		
		for(int y=0; y<a.size(); y++) 
		{
			if(lastcomma==true) {if(a[y]==',' || a[y-2]=='[' || a[y]==']') modelnameempty=true;  lastcomma=false;}

			if(a[y]==',' && lastcomma==false) {commacount++; lastcomma=true;}
		}			

		if(type!=1)
		{
			string protcod; if(type==2) protcod="AMINOACID"; else protcod="CODON";
			controlerrorprint2("[SITES]", name, a, "The command C[x,y,z] is used to specify different models or branch\nclasses for each codon position in a NUCLEOTIDE analysis.\nIt cannot be used for a "+protcod+" type analysis.",""); 
			{if(breakonerror) return -1;}  
		}

		if(a[0]!='C' || a[1]!='[' || a[a.size()-1]!=']') 
		{controlerrorprint2("[SITES]", name, a, "Only one statement after the sites class name means you are specifiying\na codon sites class and you must use the following notation:\nC[mymodelname1,mymodelname2,mymodelname3] or\nC[mybranchclass1,mybranchclass2,mybranchclass3]",""); {if(breakonerror) return -1;}  }
		
		
		

		if( commacount!=2)
		{
			stringstream sss; sss<<commacount+1; string cc=sss.str();
			controlerrorprint2("[SITES]", name, a, "You must specify 3 models in a codons site class, you specified "+cc+"\nCorrect format of the statement is:\nC[mymodelname1,mymodelname2,mymodelname3] or\nC[mybranchclass1,mybranchclass2,mybranchclass3]",""); {if(breakonerror) return -1;}  
		}
		
		if(modelnameempty )
		{
			controlerrorprint2("[SITES]", name, a, "This codon sites class contains an empty model name.\nCorrect format of the statement is:\nC[mymodelname1,mymodelname2,mymodelname3] or\nC[mybranchclass1,mybranchclass2,mybranchclass3]",""); {if(breakonerror) return -1;}  
		}

		bool writeon=false;
		for(int y1=0; y1<a.size(); y1++)
		{
			char c=a[y1];

			if(c==',' || c==']') 
			{
				if(allAinB(b,"0.123456789")) { controlerrorprint2("[SITES]", name, "", "This model/branchclass name is a decimal number, is this a mistake?\n\""+a+"\" cannot be used as a model/branchclass name.",""); {if(breakonerror) return -1;} }
				if(!allAinB(b,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[SITES]", name, "", "Model/Branchclass names should only contain ,.;()-_ and alpha-numeric characters.",""); {if(breakonerror) return -1;} }

				names.push_back(b); b="";
			}
			if(writeon && c!=',') b+=c;
			if(c=='[') writeon=true;
			
		}

		// for(int ii=0; ii<names.size(); ii++) cout<<names.at(ii)<<"  "<<endl;

		
	
		totalsitenames.push_back(name); 
		totalsites.push_back(siteclass(name, names,props,totalmodelnames,totalbranchesnames));

	}
	else if(mymy%2!=1){controlerrorprint2("[SITES]", name, "", "Sites blocks must have an even number of statements in pairs of the form:\n'proportion model' or 'proportion branchclass' e.g. 0.2 mymodelname.",""); {if(breakonerror) return -1;}  }
	else if(mymy==3){controlerrorprint2("[SITES]", name, "", "Sites blocks need at least two site classes otherwise\njust evolve a model or a branchclass not a siteclass!",""); {if(breakonerror) return -1;}  }
	else //aaa
	{
		for(int i1=1; i1<mymy; i1++) 
		{
			string a=block.at(i1);
			if(a[0]=='[' && a[a.size()-1]==']')  {controlerrorprint2("[SITES]", name, a, "Sites blocks do not contain other commands.\nHave you made a spelling mistake?",""); {if(breakonerror) return -1;}  }

			if(i1%2==0) 
			{ 
				if(allAinB(a,"0.123456789")) { controlerrorprint2("[SITES]", name, "", "This model/branchclass name is a decimal number, is this a mistake?\n\""+a+"\" cannot be used as a model/branchclass name.",""); {if(breakonerror) return -1;} }
				if(!allAinB(a,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[SITES]", name, "", "Model/Branchclass names should only contain ,.;()-_ and alpha-numeric characters.",""); {if(breakonerror) return -1;} }
			}

			else  {if(!allAinB(a,"0123456789.")) { controlerrorprint2("[SITES]", name, "", "Possible Spelling Mistake...\nProportions should be decimal numbers, control file had:\n"+a+"\nIs this a model name or branch class name?\nIf so, the proportion should come before the name.",""); {if(breakonerror) return -1;} }}
		}
			
		double sum=0; 

		//need to have C[x.y.z] for codon dna  models!

		//and to set up things for putting it all in to a site class!

		//what about the debugging for the whole thing?

		for(int i=1; i<mymy; i=i+2) 
		{
			string a=block.at(i);
			double num=atof(a.c_str());
			sum+=num;
			props.push_back(num);
			names.push_back(block.at(i+1));	
		}

		double diff=sum-1; if(diff<0) diff=-diff;
		if(diff>0.0001) 
		{
			stringstream qwe; qwe<<sum; string blah=qwe.str();
			controlerrorprint2("[SITES]", name, "", "Site class proportions summed to "+blah+" not 1,\nThey have been rescaled so that they sum to 1",""); 
			for(int j=0; j<props.size(); j++) {props.at(j)/=sum;}// cout<<props.at(j)<<" ";} cout<<endl;
			
		}

		for(int k=1; k<props.size(); k++) props.at(k)+=props.at(k-1);
		
		//for(int k1=0; k1<props.size(); k1++) cout<<props.at(k1)<<" "; cout<<endl;

		diff=props.back()-1; if(diff<0) diff=-diff;
		if(diff>0.0001) cout<<"INTERNAL ERROR in rescaling site class proportions."<<endl;

		//for(int ii=0; ii<names.size(); ii++) cout<<names.at(ii)<<"  "<<props.at(ii)<<endl;

		totalsitenames.push_back(name); 
		totalsites.push_back(siteclass(name, names,props,totalmodelnames,totalbranchesnames));

		
		if((totalsites.back()).error==-1) return -1;

	}//aaa

	


	return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int dealwithbranches(vector<string> &block, bool iscodefixed)		
{
	// this function reads in the branch class information and pre-processes it for the branch class constructor

	int mymy=block.size();
	string tree;
	string name;
	if(mymy==0) {controlerrorprint2("[BRANCHES]", "", "", "You must specify a name for a branches block. This one is empty.",""); {if(breakonerror) return -1;}  }
	else
	{
		name=block.at(0);
		if(!allAinB(name,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[BRANCHES]", name, "", "First statement in a [BRANCHES] block must be a branches name statement.\nThe name should only contain  ,.;()-_ and alpha-numeric characters.",""); {if(breakonerror) return -1;} }		

		if(checkthere(name,totalbranchesnames)!=-1)
		{ controlerrorprint2("[BRANCHES]", name, "", "A branches class with this name has been specified twice.",""); {if(breakonerror) return -1;} }
		
	
	}

	if(mymy==1) {controlerrorprint2("[BRANCHES]", name, "", "This branches block contains no tree.",""); {if(breakonerror) return -1;}  }
	else name=block.at(0);

	for(int i=1; i<mymy; i++) 
	{
		string s=block.at(i);
		if(s[0]=='[' && s[s.size()-1]==']')  {controlerrorprint2("[BRANCHES]", name, s, "Branch blocks do not contain other commands.\nHave you made a spelling mistake?",""); {if(breakonerror) return -1;}  }
		else tree+=s;
//		cout<<"branches "<<block.at(i)<<endl;
	}
	
	//cout<<tree<<endl;

	bool treeerror=false;
	char c='q';

	int size=tree.size();

	if(tree[0]!='(' || tree[size-1]!=';' ) //|| tree[size-2]!=')')
	{controlerrorprint2("[BRANCHES]", name, "", "User defined tree does not appear to be a NEWICK format tree (...); ",""); {if(breakonerror) return -1;}  }

	int mybracketleft=0, mybracketright=0, taxaon=0, lengthon=0, rootbranches=0, bracketlevel=0; 

	for(int fg1=0; fg1<size; fg1++)  
	{
		c=tree[fg1]; 
		if(c=='(') {mybracketleft++; bracketlevel++;} 
		if(c==')') {mybracketright++; bracketlevel--;}
		if(bracketlevel==1 && c==',') rootbranches++;
	}
			
	if(mybracketleft!=mybracketright) 	
	{treeerror=true; controlerrorprint2("[BRANCHES]", name, "", "Number of parantheses in guide tree do not match.",""); {if(breakonerror) return -1;}  }

	if(mybracketleft==0)  
	{treeerror=true; controlerrorprint2("[BRANCHES]", name,"", "Branch class guide tree must contain as least one set of parentheses.",""); {if(breakonerror) return -1;}  }

	bool myerror=true;
	//cout<<"W"<<tree<<"W"<<endl;
	int hashcount=0;
	for(int yg=0; yg<size; yg++)
	{
		c=tree[yg];
		if(c=='#') hashcount++;
		
		if(c==',' || c==')')
		{
			if(hashcount==0) {controlerrorprint2("[BRANCHES]", name,"", "Branch class guide tree contains a branch without a model name.\nModel names must have a # character in front of them.",""); {if(breakonerror) return -1;}  }
			if(hashcount> 1) {controlerrorprint2("[BRANCHES]", name,"", "Branch class guide tree contains a branch with more than one model name.\nModel names must be preceded with one single #.",""); {if(breakonerror) return -1;}  }
			hashcount=0;
		}
	}

	// this adds an extra model to the root for the branch class 
	// corresponds to extra branch added for forced insertions in dealwithtrees function
#ifdef checkingindelsubs
	string sd;
	for(int p=tree.size()-1; p>-1; p--) 
	{	
		if(tree[p]=='#') break;
		if(tree[p]!=';' ) sd+=tree[p];
	}
	string temp="("; 
	temp+=tree;
	temp[temp.size()-1]=')';
	temp+="#";
	for(int p2=sd.size()-1; p2>-1; p2--) temp+=sd[p2]; 
	temp+=";";
	tree=temp;
#endif
		
	totalbranchesnames.push_back(name); 

	totalbranches.push_back(branchclass(tree,name,totalmodelnames,iscodefixed));

	if((totalbranches.back()).error==-1) return -1;

	return 0;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
 
vector<unsigned long> treeseeds;   //storage for tree seeds.....
	

class Tree
{
	// this class stores a tree as defined by the user or creates random trees based on user defined settings...

public:
	bool treerandom;
//	bool randomeveryrep;
	bool seedgiven;
	string tree;
	string name;
	int ntaxa, seed, guidetreeroot, last;
	double scaler, scalerd, scalerm, birth, death, sample, mut, max_distance;
	double treelength, treedepth;
	MTRand mymtrand;
	string doctorlengths;
	
	Tree(bool myseedgiven, string mtree, double mscaler, double mscalerd, int mntaxa, double mbirth,double mdeath,double msample, double mmut, int mseed, int mguidetreeroot, string mname, string mydoctor, double mscalerm)
	{
		tree=mtree;
		ntaxa=mntaxa;
		scaler=mscaler;
		scalerd=mscalerd;
		treerandom=false;
		seedgiven=myseedgiven;
		name=mname;
		doctorlengths=mydoctor;
		scalerm=mscalerm;
		max_distance=0;

		if(tree=="RANDOM")
		{
			treerandom=true;
			birth=mbirth;
			death=mdeath;
			sample=msample;
			mut=mmut;
			seed=mseed;
			guidetreeroot=mguidetreeroot;

//			if(seedgiven) {mymtrand.seed(seed);tester.seed(seed);}// cout<<"name "<<name<<" SEEDDDD "<<seed<<endl;} 	
		}


		if(!treerandom)
		{
			if(doctorlengths=="ULTRAMETRIC") tree=randomise_ultrametric(tree);

//			cout<<"Q"<<max_distance<<"Q"<<endl;
			treedepth=getmaxtreedepth(tree);	
//			cout<<"R"<<max_distance<<"R"<<endl<<endl;
			treelength=gettreelength(tree); 

			if(scalerm!=-1) tree=rescaledtree(tree, max_distance, scalerm);
			if(scaler!=-1)  tree=rescaledtree(tree, treelength,   scaler);
			if(scalerd!=-1) tree=rescaledtree(tree, treedepth,    scalerd);////tree=redepthtree(tree, gettreelength(tree), scaler);

			max_distance=0;
			treedepth=getmaxtreedepth(tree);
//			cout<<"S"<<max_distance<<"S"<<endl<<endl;
			treelength=gettreelength(tree); 
//			cout<<"T"<<max_distance<<"T"<<endl<<endl;
			
			if(scaler!=-1)  { double diff=treelength-scaler;    if(diff<0) diff=-diff; if(diff>0.0001) controlerrorprint2("[TREE]", name,"", "The treelength adjustment is not within the limits I expected.\n Please report this and it will be promptly fixed!","");}
			if(scalerd!=-1) { double diff=treedepth-scalerd;    if(diff<0) diff=-diff; if(diff>0.0001) controlerrorprint2("[TREE]", name,"", "The treedepth adjustment is not within the limits I expected.\n Please report this and it will be promptly fixed!","");}
			if(scalerm!=-1) { double diff=max_distance-scalerm; if(diff<0) diff=-diff; if(diff>0.0001) controlerrorprint2("[TREE]", name,"", "The maxdistance adjustment is not within the limits I expected.\n Please report this and it will be promptly fixed!","");}
		}


	}

	////////////
	void newrandom(int rep) //, unsigned long &last)
	{
		// this function should be called at the beginning of simulating each replicate when using random trees.

		/*
		if(rep==1) 
		{
			// [globalseed] overrides [seed] and [randomseed]
			if(globalseed) mymtrand.seed(mtrand1.randInt()); 

			// otherwise if a [TREE] [seed] has been given it is used
			else if(seedgiven) mymtrand.seed(seed); 
			
			// otherwise the random seed is automatically generated.
			else mymtrand.seed(); 
		
		} 	
		else mymtrand.seed(last);
	
		last = mymtrand.randInt();
		
		string temptree=newrandomtree(ntaxa,birth,death,sample,mut,last,guidetreeroot); 

		*/

		if(rep==1) {if(seedgiven) mymtrand.seed(seed);}

		int theseed;  if(seedgiven) theseed=mymtrand.randInt(); else theseed=mtrand1.randInt();

		string temptree=newrandomtree(ntaxa,birth,death,sample,mut,theseed,guidetreeroot); 
		
		// FOR SOME REASON THE CODE FOR RANDOM TREES FROM ZIHENG'S EVOLVER LABELS SPECIES AS 1, 2, 3,... FOR
		// NTAXA <20 BUT AS S1, S2, S3... FOR NTAXA>=20!!  THIS FILTERS THE 'S' OUT.

//		cout<<"Q "<<temptree<<endl;

		tree="";
		for(int yq=0; yq<temptree.size(); yq++) if(temptree[yq]!='S') tree+=temptree[yq];

//		cout<<"W "<<tree<<endl;
	
		if(scaler!=-1) tree=rescaledtree(tree, gettreelength(tree), scaler);

		if(scalerd!=-1) tree=rescaledtree(tree, getmaxtreedepth(tree),scalerd);
			
	//	for(int po=0; po<40; po++) cout<<tree[po]; cout<<endl;

		treelength=gettreelength(tree);
		treedepth=getmaxtreedepth(tree);
	}


private:
	////////////
	
	
	double correctlastdepth(double maxdepth, string originaltree, double depthsofar, vector<string> &originals, vector<string> &changes)
	{
		vector<string> bits;
		string s,t;

		int bracketcount=1, pi;

		if(originaltree[0]=='(')
		{
			//cout<<"Q"<<originaltree<<"W"<<endl;
			for(pi=1; pi<originaltree.size(); pi++)
			{
				char c=originaltree[pi];
				if(c=='(') bracketcount++;
				if(c==')') bracketcount--;

				//cout<<c<<"  "<<bracketcount<<endl;
				if(bracketcount==1 && c==',') {bits.push_back(s); /*cout<<"S"<<s<<"S"<<endl;*/ s=""; continue;}
				if(bracketcount==0) {bits.push_back(s);/* cout<<"SS"<<s<<"SS"<<endl; */break;}

				s+=c;
			}

			//for(int yf=0; yf<bits.size(); yf++) cout<<bits.at(yf)<<endl;
		
//			if(!isitfirstgoo)
//			{
				s="";
				
				pi++;
				if(originaltree[pi]!=':') { controlerrorprint2("[TREE]", name,"", "Something is wrong in the getnextdepthcommand.\n Please report this and it will be promptly fixed!",""); }
		
				
				pi++;
				for(int yg=pi; yg<originaltree.size(); yg++) s+=originaltree[yg];

				/*cout<<"length to add "<<s<<endl;*/

				depthsofar+=atof(s.c_str());

//			}

			for(pi=0; pi<bits.size(); pi++) correctlastdepth(maxdepth, bits.at(pi), depthsofar, originals, changes);
		}	
		else
		{
			//end of the tree!!!

			bool lengthon=false;

			for(pi=0; pi<originaltree.size(); pi++)
			{
				if(lengthon) s+=originaltree[pi];

				if(!lengthon) {	if(originaltree[pi]==':') lengthon=true;}
			}

			//cout<<"END LENGTH "<<s<<endl;

			originals.push_back(s);
			
			double endbranch=atof(s.c_str());

			depthsofar+=endbranch;
			
			double diff=maxdepth-depthsofar;
			stringstream newendbranch;
			newendbranch<<endbranch+diff;
			
			changes.push_back(newendbranch.str());
			
			//cout<<"Depth is "<<depthsofar<<" but we have target of "<<maxdepth<<" which is difference of "<<diff<<endl;
			//cout<<"Old branch was "<<s<<" so now it is "<<newendbranch.str()<<endl;
			//cout<<"******************"<<endl;
		}
		
		return 0;
	}	
	
	string randomise_ultrametric(string originaltree)
	{
		vector<string> onestochange, whattochangeto;
		double diff, depthsofar=0;
		double maxdepth=getmaxtreedepth(originaltree);

		string tree;
		int sizey=originaltree.size()-1;
		if(originaltree[sizey]==';')
		{
			for(int fr=0; fr<sizey; fr++) tree+=originaltree[fr];
			tree+=":0.0";
		}
		else tree=originaltree;
	
		correctlastdepth(maxdepth, tree, depthsofar, onestochange, whattochangeto);
				
		int branch=0;
		string mylength, newtree;
		char c;
			
		tree=originaltree;

		for(int uj=0; uj<tree.size(); uj++)
		{
			c=tree[uj];
			if(c==',' || c==')') 
			{
				branch=0; 
				int pos=-1;
				for(int rf=0; rf<onestochange.size(); rf++)
				{
					//cout<<rf<<" "<<mylength<<"  "<<onestochange.at(rf)<<"  ";
					if(mylength==onestochange.at(rf)) {pos=rf; /*cout<<whattochangeto.at(pos)<<endl;*/ break;}
					//cout<<endl;
				}
				if(pos==-1) newtree+=mylength;	
				else newtree+=whattochangeto.at(pos);
				mylength="";
			}
			if(branch==1) mylength+=c; else newtree+=c;
			if(c==':') branch=1;
		}

	  return newtree;
		
	}
	
	double getmaxtreedepth(string originaltree)
	{
		string ss, newtree;
		double depthsofar=0;
		int sizey=originaltree.size()-1;
		if(originaltree[sizey]==';')
		{
			for(int fr=0; fr<sizey; fr++) newtree+=originaltree[fr];
			newtree+=":0.0";
		}
		else newtree=originaltree;
		
//		if(depths.empty()) controlerrorprint2("[TREE]", name,"", "The treedepth command has misbehaved and the depths vector is empty.\n Please report this and it will be promptly fixed!","");
				
		return getnextdepth(newtree,depthsofar);
	}
	
	double getnextdepth(string originaltree, double depthsofar)
	{
		vector<string> bits;
		string s,t;

		int bracketcount=1, pi;

		if(originaltree[0]=='(')
		{
			//cout<<"Q"<<originaltree<<"W"<<endl;
			for(pi=1; pi<originaltree.size(); pi++)
			{
				char c=originaltree[pi];
				if(c=='(') bracketcount++;
				if(c==')') bracketcount--;

				//cout<<c<<"  "<<bracketcount<<endl;
				if(bracketcount==1 && c==',') {bits.push_back(s); /*cout<<"S"<<s<<"S"<<endl;*/ s=""; continue;}
				if(bracketcount==0) {bits.push_back(s);/* cout<<"SS"<<s<<"SS"<<endl; */break;}

				s+=c;
			}

			//for(int yf=0; yf<bits.size(); yf++) cout<<bits.at(yf)<<endl;
		
//			if(!isitfirstgoo)
//			{
				s="";
				
				pi++;
				if(originaltree[pi]!=':') { controlerrorprint2("[TREE]", name,"", "Something is wrong in the getnextdepthcommand.\n Please report this and it will be promptly fixed!",""); }
		
				
				pi++;
				for(int yg=pi; yg<originaltree.size(); yg++) s+=originaltree[yg];

				/*cout<<"length to add "<<s<<endl;*/

				depthsofar+=atof(s.c_str());

//			}

			vector<double> mydepths, mydistances;
			vector<vector<double> > distances;
			double mymax=0;
			int sizey=bits.size();
 
//			for(pi=0; pi<sizey; pi++) cout<<"R "<<bits.at(pi)<<endl; cout<<"****"<<endl;

// these two lines ended up stuck in a loop on very deep trees!

		//	for(pi=0; pi<sizey; pi++) mydepths.push_back(getnextdepth(bits.at(pi),depthsofar));

		//	for(pi=0; pi<sizey; pi++) mydistances.push_back(getnextdepth(bits.at(pi),0));

			for(pi=0; pi<sizey; pi++) 
			{
				double newdepth=getnextdepth(bits.at(pi),depthsofar);
				mydepths.push_back(newdepth);
				mydistances.push_back(newdepth-depthsofar);
			}

		
			for(pi=0; pi<sizey; pi++) 
			{
	//			cout<<"    "<<mydepths.at(pi)<<endl;
				if(mydepths.at(pi)>mymax) mymax=mydepths.at(pi);
				
				for(int pi2=0; pi2<sizey; pi2++)
				{
					if(pi==pi2) continue;
				
					double x=mydistances.at(pi)+mydistances.at(pi2);
	//				cout<<"Q  "<<mydistances.at(pi)<<"  "<<mydistances.at(pi2)<<"  "<<x<<"  "<<max_distance<<endl;
					if(x>max_distance) max_distance=x;
				}
			}
			
			return mymax;
		}	
		else
		{
			//end of the tree!!!

			bool lengthon=false;

			for(pi=0; pi<originaltree.size(); pi++)
			{
				if(lengthon) s+=originaltree[pi];

				if(!lengthon) {	if(originaltree[pi]==':') lengthon=true;}
			}

			//cout<<"END LENGTH "<<s<<endl;

			depthsofar+=atof(s.c_str());

			return depthsofar;
			//depths.push_back(depthsofar);

		}
		
		return 0;
	}	

	////////////
	double gettreelength(string originaltree)
	{
		// returns tree length of originaltree

		int branch=0;
		string mylength;
		double treelength=0;
		char c;

		for(int uj=0; uj<originaltree.size(); uj++)
		{
			c=originaltree[uj];
			if(c==',' || c==')') {branch=0; double length=atof(mylength.c_str()); treelength+=length; mylength="";}
			if(branch==1) mylength+=c;
			if(c==':') branch=1;
		}

	  return treelength;
	}	 
	/////

	string rescaledtree(string originaltree, double treelength, double newtreelength)
	{
		// this function will rescale a tree according to a desired total tree length.

		int branch=0;
		string mylength, newtree;
		char c;

		for(int uj=0; uj<originaltree.size(); uj++)
		{
			c=originaltree[uj];
			if(c==',' || c==')') 
			{
				branch=0; double length=atof(mylength.c_str()); 
  				length=length*newtreelength/treelength; 
				stringstream fd; fd<<length; string fg=fd.str(); 
				newtree+=fg;
				mylength="";
			}
			if(branch==1) mylength+=c;
			else newtree+=c;
			if(c==':') branch=1;
		}

	  return newtree;
	}
	///////////

};
//////////////////
vector<Tree>	totaltrees;	   // storage for trees!

//////////////////////////

int dealwithtrees(vector<string> &block)		
{
	// this reads in tree information and pre-processes it for the tree class.
	// it also checks whether the tree is inconsistent in any way.

	int seed		  = 1;
	int mysize		  = block.size();	
	int blocksize	  = mysize-1;
	int ntaxa         = 4;
	
	bool seedgiven=false;
	int guidetreeroot = -1;
	double scaler = -1, scalerd=-1, scalerm=-1;

	// change the following to default values used in Yang et al. ?

	// or make it so they must always be entered?

	double birth  = 1;
	double death  = 1;
	double sample = 1;
	double mut    = 1;

	string commandsarray[8]={"[user]","[rooted]","[unrooted]","[seed]","[treedepth]","[treelength]","[branchlengths]","[maxdistance]"};
	
	string tree;
	string doctorlengths="NO";
 
	//rooted or unrooted then: "ntaxa" or "ntaxa birthrate deathrate samplingrate mutationrate"
	
	vector<string> commands(commandsarray,commandsarray+8);

	//	cout<<mysize<<" is block size"<<endl;
	
	//	for(int i=0; i<mysize; i++) cout<<"trees "<<block.at(i)<<endl;
	
	if(mysize==0) {controlerrorprint2("[TREE]", "", "", "No tree name has been given","");  {if(breakonerror) return -1;} }

	string name=block.at(0);

	if(mysize==1) {controlerrorprint2("[TREE]", name, "","The tree block is empty.","");  {if(breakonerror) return -1;} }
	
	if(mysize==2) if(block.at(1)=="[user]") {controlerrorprint2("[TREE]", name, "","There is no tree in the tree block.","");  {if(breakonerror) return -1;} }
	
	string lastbit=block.at(blocksize);
	if(lastbit[0]=='[') {controlerrorprint2("[TREE]", name, lastbit, "The last entry in a [TREE] block cannot be a command.","");  {if(breakonerror) return -1;} }

	bool well[8]={false,false,false,false,false,false,false,false};

	for(int p1=0; p1<mysize; p1++) 
	{
		string ss=block.at(p1);
		for(int p2=0; p2<8; p2++) if(ss==commands.at(p2)) well[p2]=true;
	}

	if(well[0]) for(int p3=1; p3<4; p3++) if(well[p3]) {controlerrorprint2("[TREE]", name, commands.at(p3),"This command cannot be used in a user-defined tree block.",commands.at(p3));  {if(breakonerror) return -1;} }

	bool doultra=false;

	if((well[6])) 
	{
		for(int i1=1; i1<blocksize+1; i1++) 
		{
			string hg=block.at(i1); 
			if(hg=="NON-ULTRAMETRIC" || hg=="ULTRAMETRIC") doultra=true;
		}
	}

	if( !(well[0]) && well[6]) {controlerrorprint2("[TREE]", name, commands.at(6),"This command cannot be used in a random tree block.",commands.at(6));  {if(breakonerror) return -1;} }
	
	
/*

[rooted]
[unrooted] ntaxa birthrate deathrate samplingrate mutationrate
[guidetreeseed]
[treelength]

*/
	string hg=block.at(0), hg1;
//	if(hg[0]!='#' || hg[hg.size()-1]!='#') { controlerrorprint2("[MODEL]", "?", "", "First statement in a [MODEL] block must be a model name statement in the form #modelname#",""); {if(breakonerror) return -1;} }
	if(!allAinB(hg,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[TREE]", name, "", "First statement in a [TREE] block must be a tree name statement.\nThe name should only contain  ,.;()-_ and alpha-numeric characters.",hg); {if(breakonerror) return -1;} }
	else
	{
		if(checkthere(name,totaltreenames)!=-1) { controlerrorprint2("[TREE]", name, "", "A tree with this name has been specified twice.",""); {if(breakonerror) return -1;} }
		else totaltreenames.push_back(name); 
		
	
		vector<string> tempvec;
		
		string hg="";
		int lasttest=-999;
		int mytest=-999;

		for(int i1=1; i1<blocksize+1; i1++) 
		{
			hg=block.at(i1); 
		
			lasttest=mytest;
			if(block.at(i1-1)!="[branchlengths]") mytest=teststring(lasttest, hg, commands, "[TREE]","0123456789.",name);

			if(mytest==-2) {if(breakonerror) return -1;} 
			else 
			{//a1
			
				if(i1==blocksize )  {tempvec.push_back(hg); } //cout<<"H1 "<<hg<<endl;}
				
				if(mytest==lasttest && i1!=blocksize )  {tempvec.push_back(hg);}// cout<<"H1 "<<hg<<endl;}
				else if(i1!=1)
				{//a2
				
					int tempsize=tempvec.size();

					if(tempsize==0) {controlerrorprint2("[TREE]", name, commands.at(lasttest), "No values were found after this command in the control file.",""); {if(breakonerror) return -1;}  }
					
					else if(lasttest>2 && lasttest<8)
					{
						if(tempsize!=1) {/* cout<<tempvec.at(0)<<" "<<tempvec.at(1)<<" "<<tempvec.at(1 )<<" E"<<endl;*/ controlerrorprint2("[TREE]", name, commands.at(lasttest), "There should be only one value after this command in the control file.",""); {if(breakonerror) return -1;}  }
	
						string bit=tempvec.at(0);
						
						if(lasttest==3) 
						{
							if(AinB('.',bit)) {controlerrorprint2("[TREE]", name, commands.at(lasttest), "First value after this command should be an integer.",bit); {if(breakonerror) return -1;}  }

							seedgiven=true;
							seed=(atoi(bit.c_str()));	
							// seed for random tree
							//myinclude.at(8)=true;
						}
						else if(lasttest==5) 
						{
							scaler=atof(bit.c_str());				
							// scaler for treelength				
							// myinclude statement here?
						}
						else if(lasttest==4) 
						{
							scalerd=atof(bit.c_str());				
							// scaler for treedepth				
							// myinclude statement here?
						}
						else if(lasttest==6)
						{
							if(bit=="ULTRAMETRIC" || bit=="NON-ULTRAMETRIC" || bit=="EQUAL") doctorlengths=bit;
							else  {controlerrorprint2("[TREE]", name, commands.at(lasttest), "This command can take the following values\n\"ULTRAMETRIC\"\n\"NON-ULTRAMETRIC\"\n\"EQUAL\"\nYou entered:",bit); {if(breakonerror) return -1;}  }
						}
						else if(lasttest==7) 
						{
							scalerm=atof(bit.c_str());				
							// scaler for treedepth				
							// myinclude statement here?
						}
					}
					else if(lasttest==1 || lasttest==2)
					{
						if(lasttest==1) guidetreeroot=2; //rooted
						else			guidetreeroot=1; //unrooted

						string bit=tempvec.at(0);

						if(AinB('.',bit)) {controlerrorprint2("[TREE]", name, commands.at(lasttest), "First value after this command should be an integer.",bit); {if(breakonerror) return -1;}  }
						
						// myinclude.at(3)=true;
						ntaxa=atoi(bit.c_str());

						if(tempsize!=1 && tempsize!=5)  
						{
							stringstream fg; fg<<tempsize; string fd=fg.str();
							controlerrorprint2("[TREE]", name, commands.at(lasttest), "Found "+fd+" values but was expecting 1 integer, or \n1 integer and 4 decimal values after this command.",""); 
							{if(breakonerror) return -1;}  
						}
						
						
						if(tempsize==5)
						{
							//myinclude.at(4)=true;
							//myinclude.at(5)=true;
							//myinclude.at(6)=true;
							//myinclude.at(7)=true;
							
							birth=atof((tempvec.at(1)).c_str());
							death=atof((tempvec.at(2)).c_str());
							sample=atof((tempvec.at(3)).c_str());
							mut=atof((tempvec.at(4)).c_str());
						}
						guidetreetype=2;

					}
					else 
					{
						// lasttest=0  [user]  user specified tree			
						bool treeerror=false;
						char c='q',c1='q',c2='q',c3='q',c4='q',c5='q',c6='q',c7='q',c8='q',c9='q',c10='q';
						string rawvalue;
						if(tempvec.size()==1) {rawvalue=tempvec.at(0);}// cout<<"RAWVALUE 1 "<<rawvalue<<endl;}
						else 
						{
							for(int jh=0; jh<tempvec.size(); jh++) 
							{
								string bl=tempvec.at(jh);

								//if(bl[0]=='[')	{controlerrorprint2("[TREE]", name, commands.at(lasttest), "No other command can be specified in a tree block with a user-defined tree",bl); {if(breakonerror) return -1;}  }

								rawvalue+=bl;
							}
						}

#ifdef checkingindelsubs

// when checking substitutions for inserted sites, 
// this bit of code adds a very small branch to the 
// root of the tree so that many insertions can
// be "forced" and then these inserted sites can be
// observed to evolve along the tree for testing.

	string temp="("; 
	temp+=rawvalue;
	temp[temp.size()-1]=':';
	temp+="0.0001);";
	rawvalue=temp;
#endif
						int size=rawvalue.size();

						if(rawvalue[0]!='(' || rawvalue[size-1]!=';' || rawvalue[size-2]!=')')
						{controlerrorprint2("[TREE]", name, commands.at(lasttest), "User defined tree does not appear to be a NEWICK format tree (...); ",rawvalue); {if(breakonerror) return -1;}  }

						// if [branchlengths] command has been used then go through tree and 
						if(well[6])
						{
							string rawvalue2, rawvalue3, rawvalue4;
							char c, c1, c2;
							
							// strip lengths from tree if there
							bool skip=false;
							for(int ds2=0; ds2<rawvalue.size(); ds2++)
							{
								c=rawvalue[ds2];
								if(c==':') skip=true;
								if(skip && (c==',' || c==')') ) skip=false;
								if(!skip) rawvalue2+=c;
							}
							
							// add equal or random branch lengths back.
							bool randomise; if(doultra)randomise=true;  else randomise=false;
							//cout<<doctorlengths<<endl;
							rawvalue3+=rawvalue2[0];
							string newlength;

							for(int ds1=1; ds1<rawvalue2.size()-1; ds1++)
							{
								c1=rawvalue2[ds1-1];
								c =rawvalue2[ds1];
								c2=rawvalue2[ds1+1];

								if( (c==')' && c2!=';')  || (c2==')' && c!=')'))// ||  c1!=')') )
								{
									rawvalue3+=c;
									rawvalue3+=':';
									if(randomise)
									{
										double newlengthD=mtrand1();
										stringstream xx; 
										xx<<newlengthD;
										newlength=xx.str();
									}
									else newlength="0.1";
									rawvalue3+=newlength;
								}
								else if(c==',' && c1!=')')// && c1!=')') // || (c!=')' && c1!=')') )
								{
									
									rawvalue3+=':';
									if(randomise)
									{
										double newlengthD=mtrand1();
										stringstream xx; 
										xx<<newlengthD;
										newlength=xx.str();
									}
									else newlength="0.1";
									rawvalue3+=newlength;								
									rawvalue3+=c;
								}
								else rawvalue3+=c;
							}
							rawvalue3+=';';
						
						//	cout<<rawvalue<<endl;
							rawvalue=rawvalue3;
						//	cout<<rawvalue2<<endl;
						//	cout<<rawvalue3<<endl;
						//	cout<<rawvalue4<<endl;
						// 	cout<<rawvalue<<endl;
							size=rawvalue.size();
						}
						
						
						int mybracketleft=0, mybracketright=0, taxaon=0, lengthon=0, rootbranches=0, bracketlevel=0; 

						
						for(int fg1=0; fg1<size; fg1++)  
						{
							c=rawvalue[fg1]; 
							if(c=='(') {mybracketleft++; bracketlevel++;} 
							if(c==')') {mybracketright++; bracketlevel--;}
							if(bracketlevel==1 && c==',') rootbranches++;
						}
			
						if(mybracketleft!=mybracketright) 	
						{cout<<rawvalue<<endl<<mybracketleft<<"  "<<mybracketright<<endl;
						treeerror=true; controlerrorprint2("[TREE]", name, commands.at(lasttest), "Number of parantheses in guide tree do not match.",""); {if(breakonerror) return -1;}  }

						if(mybracketleft==0)  
						{treeerror=true; controlerrorprint2("[TREE]", name, commands.at(lasttest), "Guide tree must contain as least one set of parentheses.",""); {if(breakonerror) return -1;}  }
													
						c1='Q'; c2='Q'; 
						int taxacount=0, coloncount=0, intbranch=0, extbranch=0;
						
						string lasterror, currentlength;

						for(int fg2=0; fg2<size; fg2++) 
						{
							c10=c9; c9=c8; c8=c7; c7=c6; c6=c5; c5=c4; c4=c3; c3=c2; c2=c1;  c1=rawvalue[fg2]; 
							lasterror=""; lasterror+=c10;lasterror+=c9;lasterror+=c8;lasterror+=c7;lasterror+=c6;lasterror+=c5;lasterror+=c4;lasterror+=c3;lasterror+=c2;lasterror+=c1;
						
							if(c1==')' || c1=='(' || c1==',' && (taxaon==2 || lengthon==1) )
							{
								//if(taxaon==1 )	{treeerror=true; controlerrorprint(blocknumber, "There is no branch length in this taxon in the guide tree.",linecount,lasterror); if(breakonerror==1) {breakonerror=2; break;}}
								if(taxaon==2)	{taxaon=3; taxacount++;}
								if(lengthon==1) {lengthon=2;}
							}					
							
							if(lengthon==2 || taxaon==3)
							{
								lengthon=0; taxaon=0;  
								
								int zerowarn=0;
								for(int yg=0; yg<currentlength.size(); yg++)
								{
									c=currentlength[yg];
									if(!AinB(c,"0.")) zerowarn=1;
									if(!AinB(c,"0.123456789:")) {treeerror=true; controlerrorprint2("[TREE]", name, commands.at(lasttest), "A branch length in the guide tree appears non-numerical",""); {if(breakonerror) return -1;}  }
									
								}
							//	if(zerowarn==0 && currentlength!="") {treeerror=true; controlerrorprint2("[TREE]", name, commands.at(lasttest), "This guide tree contains a zero branch length.  This is not necessary.\nPolytomies (of any order) are allowed at root (or any node) of guide tree.",""); {if(breakonerror) return -1;}  }
							}
												
							if(lengthon==1 || taxaon==2) currentlength+=c1;
							
							if(c2==')' && c1!=';' && lengthon==0) 
							{ 
								if(c1!=':') {treeerror=true; controlerrorprint2("[TREE]", name, commands.at(lasttest), "There is no branch length after a right paranthesis in the tree.\nIt happened here \""+lasterror+"\"",""); {if(breakonerror) return -1;}  }
								else {lengthon=1; currentlength="";}
							}
							
							if(c1==')' || c1==',') 
							{
								if(!AinB(c2,"0123456789")) {treeerror=true; controlerrorprint2("[TREE]", name, commands.at(lasttest), "There is a missing branch length in the tree.\nIt happened here \""+lasterror+"\"",""); {if(breakonerror) return -1;}  }
								else {lengthon=1; currentlength="";}
							}
							
							if(c2!=')' && c2!='(' && c2!=',' && taxaon==0) {taxaon=1; }
							if(taxaon==1 && c1==':') {taxaon=2; currentlength=""; }
			
							if(c1==':') {coloncount++; if(c2==')') intbranch++; else extbranch++;}

							if(c2==':' && (c1!='0' && c1!='1' && c1!='2' && c1!='3' && c1!='4' && c1!='5' && c1!='6' && c1!='7' && c1!='8' && c1!='9' ) )
								{treeerror=true; controlerrorprint2("[TREE]", name, commands.at(lasttest), "There is a missing branch length after a : in the guide tree",""); {if(breakonerror) return -1;}  }

						}

						if(!treeerror)
						{	
							tree=rawvalue;
							guidetreetype=1;
							ntaxa=extbranch;
							if(rootbranches==0) ;//cout<<"1 taxon guide tree "<<rawvalue<<" entered. One branch will be simulated."<<endl<<" ";
							else
							{
								if(extbranch==2 || extbranch==3) ;//cout<<"Guide tree of "<<extbranch<<" taxa read in successfully:"<<endl<<" "<<rawvalue<<endl<<" "<<endl<<" ";
								else
								{
									if(intbranch==0) ;//cout<<"Guide tree read in as "<<extbranch<<" taxa star tree:"<<endl<<" "<<rawvalue<<endl<<" "<<endl<<" ";
									else
									{
										int rootvalue=0;
										if(rootbranches==1) rootvalue++;

										
										//cout<<extbranch<<" taxa "; 
										if(coloncount!=2*extbranch-3 +rootvalue) guidetreebinary=false; else guidetreebinary=true;
										if(rootbranches==1) guidetreerooted=true;
										if(rootbranches==2) guidetreerooted=false;

										//cout<<"guide tree with "<<intbranch+extbranch<<" branches read in successfully:"<<endl<<" "<<endl<<" "<<originaltree<<endl<<" "<<endl<<" ";
										//if(coloncount!=2*taxacount-3 ) {cout<<"WARNING: guide tree is not binary."<<endl<<" ";}											
									}
								}
							}
									
						
						}


					} // end of else for user tree
				
					tempvec.clear();

				}// end of i1==1 if/else

			} // end of mytest==-2 else

		}//end of for block

	} // end of initial else

	if(scaler!=-1  && scalerd!=-1) 	{controlerrorprint2("[TREE]", name,"[treedepth]", "You cannot use [treelength] and [treedepth] at the same time!",""); {if(breakonerror) return -1;}  }
	if(scaler!=-1  && scalerm!=-1) 	{controlerrorprint2("[TREE]", name,"[treelength]", "You cannot use [treelength] and [maxdistance] at the same time!",""); {if(breakonerror) return -1;}  }
	if(scalerm!=-1 && scalerd!=-1)  {controlerrorprint2("[TREE]", name,"[treedepth]", "You cannot use [maxdistance] and [treedepth] at the same time!",""); {if(breakonerror) return -1;}  }

	

	if(guidetreetype==2) {tree="RANDOM"; if(scalerm!=-1) {scalerd=0.5*scalerm; scalerm=-1;}} // random trees always ultrametric, so max. pairwise distance i
	totaltrees.push_back(Tree(seedgiven,tree,scaler,scalerd,ntaxa,birth,death,sample,mut,seed,guidetreeroot,name,doctorlengths,scalerm));
	treeseeds.push_back(seed);
	
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class partitionclass
{
	// this class stores the information for simulations that use different partitions
	// A partition may be a branch class, site class or model class.


public:
	string name;
	int ntaxa;
	bool randomtrees;

	vector<int> geneticcodes;
	vector<vector<int> > rootseqints;
	vector<int> blank;
	vector<int> mbsposvec;
	vector<int> rootmodelpos;
	vector<int> rootlengthvec;
//	vector<int> ntaxavec;
	vector<int> treeposvec;

	vector<bool> isitsites;
	vector<bool> isitbranches;
//	vector<bool> isitrandomvec;



	partitionclass(string myname, vector<int> rootlengths, vector<string> &rootseqtxtvec,  vector<string> rootfilenames, vector<int> mmbsposvec, vector<string> mbstypes, vector<string> mbnames, vector<int> mrootmodelpos, vector<int> geneticcodevec, vector<int> mytreeposvec, vector<string> mytreenamevec, int myntaxa, bool isitrandom)
	{
		// constructor

		randomtrees=isitrandom;
		geneticcodes=geneticcodevec;
		name=myname;
		mbsposvec=mmbsposvec;
		rootlengthvec=rootlengths;
		treeposvec=mytreeposvec;

		checktaxaintrees(mytreeposvec);

		ntaxa=myntaxa;
		rootmodelpos=mrootmodelpos;
//		isitrandomvec=myisitrandomvec;
		
		rootseqints.assign(rootseqtxtvec.size(),blank);

		for(int kk=0; kk<rootseqtxtvec.size(); kk++) 
		{
			string mbstype=mbstypes.at(kk);
			
			if(mbstype=="model")    {isitsites.push_back(false); isitbranches.push_back(false);}
			if(mbstype=="branches") {isitsites.push_back(false); isitbranches.push_back(true); }
			if(mbstype=="sites")    
			{
				isitsites.push_back(true);
				siteclass* s=&(totalsites.at(mbsposvec.at(kk)));
				bool now=(*s).therearebranches;
				isitbranches.push_back(now);
			}
			
			makerootseqints(rootseqints.at(kk),rootlengths, rootseqtxtvec.at(kk), kk, mbsposvec.at(kk), rootfilenames, mbstype, mbnames, geneticcodevec);
		}
	}

private:
	
	int checktaxaintrees(vector<int> posvec)
	{
		Tree t=totaltrees.at(posvec.at(0));

		int ntaxa=t.ntaxa;

		for(int p=1; p<posvec.size(); p++)
		{
			Tree t=totaltrees.at(posvec.at(p));
			if(t.ntaxa != ntaxa) 
			{
				stringstream n1, n2; n1<<ntaxa; n2<<t.ntaxa; string s1=n1.str(), s2=n2.str(); string name1=(totaltrees.at(posvec.at(0))).name, name2=(totaltrees.at(posvec.at(p))).name;
				controlerrorprint2("[PARTITIONS]", name, "", "Two trees in this partition group have a different number of taxa.\nTree "+name1+" has "+s1+" taxa.\nTree "+name2+" has "+s2+" taxa.\nThis is not allowed. Please correct.",""); 
				if(breakonerror) return -1; 	
			}

		}
	
		return 0;
	}

	
	int makerootseqints(vector<int> &rootseqint, vector<int> rootlengths, string &rootseqtxt, int kk, int mbspos, vector<string> rootfilenames, string mbstype, vector<string> mbnames, vector<int> geneticcodevec)
	{
		int rootlength=rootlengths.at(kk);
		int myrootlength=rootseqtxt.size();
		int geneticcode=geneticcodevec.at(kk);
		
		string rootfilename=rootfilenames.at(kk);
		string mbname=mbnames.at(kk);

		if(rootseqtxt=="CREATE")
		{
			// CREATION OF ROOT SEQUENCE NOW DONE IN SETUPROOT FUNCTION IN MAIN SKELETON FILE

			/*
			model* m=&(totalmodels.at(mbspos));

			vector<double> rootbasefreqs=(*m).rootbasefreqs;
			vector<double> cumfreqs; cumfreqs.push_back(rootbasefreqs.at(0));
			int rsize=rootbasefreqs.size();
			
			for(int gh1=1; gh1<rsize; gh1++) cumfreqs.push_back(rootbasefreqs.at(gh1)+cumfreqs.at(gh1-1));

			double diff=1-cumfreqs.back(); if(diff<0) diff=-diff;
			if(diff>0.0001) {controlerrorprint2("[PARTITIONS]", name, "", "INTERNAL ERROR: root base freqs do not sum to 1 in partitions class.",""); {if(breakonerror) return -1;} }	

			for(int gf=0; gf<rootlength; gf++) 
			{
				double rand=mtrand1(); 
				for(int ds=0; ds<rsize; ds++) if(rand<cumfreqs.at(ds)) {rootseqint.push_back(ds); break;}
			}

			//for(int as=0; as<rootseqint.size(); as++) cout<<as<<" "<<rootseqint.at(as)<<endl; cout<<endl;	
			*/
		}
		else
		{
				if(type==2)
				{
					string test="ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv";
					if(!allAinB(rootseqtxt,test)) { controlerrorprint2("[PARTITIONS]", name, "", "AMINOACID root sequence in file can only contain following letters:\nARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv",""); {if(breakonerror) return -1;} }
					
					int size=20;

					/*
						//always equal now
					if(rootlength!=myrootlength) 
					{
						stringstream rr1; rr1<<rootlength;   string r1=rr1.str();
						stringstream rr2; rr2<<myrootlength; string r2=rr2.str();
						controlerrorprint2("[PARTITIONS]", name, "", "In partitions command you specified a root length of "+r1+"\nbut in file "+rootfilename+"\nthe root sequence has a length of "+r2+" amino acids.\nThe two numbers should match.  Have you made a mistake?",""); 
						{if(breakonerror) return -1;} 
					}
					*/

					for(int fd0=0; fd0<myrootlength; fd0++)
					{
						char c=rootseqtxt[fd0];
						bool error=true;
						
						for(int fd1=0; fd1<size; fd1++) if(c==test[fd1] || c==test[fd1+size]) {error=false; rootseqint.push_back(fd1); break;}
						
						if(error) { controlerrorprint2("[PARTITIONS]", name, "", "INTERNAL ERROR when making NUCLEOTIDE root sequence.",""); {if(breakonerror) return -1;} }		
					}
				}
				else
				{
					string test="TCAGtcag", bit;
					if(type==1) bit="NUCLEOTIDE"; else bit="CODON";
				
					if(!allAinB(rootseqtxt,test)) { controlerrorprint2("[PARTITIONS]", name, "", bit+" root sequence in file can only contain following letters:\nTCAGtcag",""); {if(breakonerror) return -1;} }		
			
					

					int size=4;
					if(type==1) 
					{
						/*
						//always equal now
						if(rootlength!=myrootlength) 
						{
							stringstream rr1; rr1<<rootlength;   string r1=rr1.str();
							stringstream rr2; rr2<<myrootlength; string r2=rr2.str();
							controlerrorprint2("[PARTITIONS]", name, "", "In partitions command you specified a root length of "+r1+"\nbut in file "+rootfilename+"\nthe root sequence has a length of "+r2+" nucleotide bases.\nThe two numbers should match.  Have you made a mistake?",""); 
							{if(breakonerror) return -1;} 
						}	
						*/
						size=4; 
						for(int fd0=0; fd0<myrootlength; fd0++)
						{
						//	char cc=s[fd0];
						//	stringstream ccc; ccc<<cc; string c=ccc.str();
						
							bool error=true;
							char c=rootseqtxt[fd0];
						
							for(int fd1=0; fd1<size; fd1++) if(c==test[fd1] || c==test[fd1+size]) {error=false; rootseqint.push_back(fd1); break;}
														
							if(error) { controlerrorprint2("[PARTITIONS]", name, "", "INTERNAL ERROR when making "+bit+" root sequence.",""); {if(breakonerror) return -1;} }		
						}
					}
					else
					{ 
						if(myrootlength%3!=0) 
						{
							stringstream rd; rd<<myrootlength; string rdd=rd.str();
							controlerrorprint2("[PARTITIONS]", name, "", "CODON root sequence length is "+rdd+" nucleotides in the file\nbut must be a multiple of 3.",""); 
							{if(breakonerror) return -1;} 
						}		
					
						/*
						if(rootlength!=myrootlength/3) 
						{
							stringstream rr1; rr1<<rootlength;   string r1=rr1.str();
							stringstream rr2; rr2<<myrootlength/3; string r2=rr2.str();
							controlerrorprint2("[PARTITIONS]", name, "", "In partitions command you specified a root length of "+r1+"\nbut in file "+rootfilename+"\nthe root sequence has a length of "+r2+" codons.\nThe two numbers should match.  Have you made a mistake?",""); 
							{if(breakonerror) return -1;} 
						}
						*/

						///////////////////////////////////////////////////////////////////////

						vector<int> notallowed=getstops(geneticcode);

						/*  //OLD WAY
						vector<int> notallowed;

						if(geneticcode!=6 && geneticcode!=12) notallowed.push_back(10);
						if(geneticcode!=6 && geneticcode!=13) notallowed.push_back(11);
						if(geneticcode==0 || geneticcode==6 || geneticcode==9 || geneticcode==10 || geneticcode==13) notallowed.push_back(14);
						if(geneticcode==1) {notallowed.push_back(46); notallowed.push_back(47);}
						*/

						int notsize=notallowed.size();

						stringstream we; we<<geneticcode; string wer=we.str();

						for(int fd0=2; fd0<myrootlength; fd0=fd0+3)
						{
							char c1=rootseqtxt[fd0-2], c2=rootseqtxt[fd0-1], c3=rootseqtxt[fd0];	
							bool error1=true, error2=true, error3=true;
							int i1,i2,i3,tot;

							for(int fd1=0; fd1<size; fd1++) 
							{
								if(c1==test[fd1])      {i1=fd1; error1=false;}
								if(c1==test[fd1+size]) {i1=fd1; error1=false;}
								if(c2==test[fd1])      {i2=fd1; error2=false;}
								if(c2==test[fd1+size]) {i2=fd1; error2=false;}
								if(c3==test[fd1])      {i3=fd1; error3=false;}
								if(c3==test[fd1+size]) {i3=fd1; error3=false;}
							
							}
								
							if(error1 || error2 || error3) { controlerrorprint2("[PARTITIONS]", name, "", "INTERNAL ERROR when making "+bit+" root sequence.",""); {if(breakonerror) return -1;} }		
							
							tot=(i1<<4)+(i2<<2)+i3; //cout<<tot<<endl;

							for(int fd2=0; fd2<notsize; fd2++) if(tot==notallowed.at(fd2)) 
							{
								stringstream m; m<<fd0/3+1; string mm=m.str();
								controlerrorprint2("[PARTITIONS]", name, "", "Root sequence file: "+rootfilename+"\nThis root sequence contains the codon "+c1+c2+c3+" at position "+mm+"\nBut "+mbstype+" class "+mbname+" uses genetic code "+wer+".\nUnder this genetic code "+c1+c2+c3+" is a stop codon and is not allowed.",""); 
								{if(breakonerror) return -1;} 
							}		
														
							rootseqint.push_back(tot); 
							
						}

						myrootlength/=3;

					} // end of type=3 else in type1/type3 else

				} // end of type1/type3 else


		}
		
		return 0;
	}

};

///////////////////////////////////////////////////////////////////////////////////////////////////

vector<partitionclass> totalpartitions;   // storage for partitions

///////////////////////////////////////////////////////////////////////////////////////////////////
int dealwithpartitions(vector<string> &block)	
{
	// this function deals with information in the control file about partitions and pre-processes it for the class constructor

	int mymy=block.size();
	string name;
	if(mymy==0) {controlerrorprint2("[PARTITIONS]", "", "", "You must specify a name for a partitions block. This one is empty.",""); {if(breakonerror) return -1;}  }
	else
	{
		name=block.at(0);
		if(!allAinB(name,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[SITES]", name, "", "First statement in a [PARTITIONS] block must be a name statement.\nThe name should only contain  ,.;()-_ and alpha-numeric characters.",""); {if(breakonerror) return -1;} }		

		if(checkthere(name,totalpartitionnames)!=-1)
		{ controlerrorprint2("[PARTITIONS]", name, "", "A partitions class with this name has been specified twice.",""); {if(breakonerror) return -1;} }
		else totalpartitionnames.push_back(name); 
	
	}

	if(mymy==1) {controlerrorprint2("[PARTITIONS]", name, "", "This partitions block contains no commands.",""); {if(breakonerror) return -1;}  }

	int count=0;
	if(mymy>1)
	{
		string ss1=block.at(1);
		string ss2=block.at(mymy-1);
		int ss2s=ss2.size()-1;
		if(ss1[0]!='[') {controlerrorprint2("[PARTITIONS]", name, "", "First statement in a partitions block must begin with a [. check input.",""); {if(breakonerror) return -1;}  }
		if(ss2[ss2s]!=']') {controlerrorprint2("[PARTITIONS]", name, "", "Last statement in a partitions block must end with a ]. check input.",""); {if(breakonerror) return -1;}  }
	}

	int mbspos;
	string mbstype,tree,treename, rootseqtxt, strippedtree;
	int treepos, rootlength;
	vector<int> ntaxavec, rootlengthvec, mbsposvec;
	vector<string> treenamevec, mbnamevec, mbstypevec, rootseqtxtvec, rootfilenamevec;
	vector<int> treeposvec, rootseqint, rootmodelpos, geneticcodevec;

	bool isitrandom=false, overallrandom=false;
	int geneticcode;

	string mbname;

	for(int k=1; k<mymy; k++)
	{
		string s=block.at(k), t;
		if(s=="") continue;
		count++;

		if(s[s.size()-1]==']') 
		{
			//cout<<"SS "<<s<<endl;
			if( count!=3 && s!="]") 
			{
				//cout<<"Q"<<s<<"Q"<<endl;
				stringstream ds; if(s=="]") count--; ds<<count; string dss=ds.str();
				controlerrorprint2("[PARTITIONS]", name, "", "Statement [..] in a partitions block must contain 3 elements not "+dss+":\n[tree1 model/branchclass/siteclass rootlength] or "+"\n[tree2 model/branchclass/siteclass rootseqfilename]\nConsult manual for more information.",""); {if(breakonerror) return -1;}  
			}
		
			t=s; s=""; for(int gj=0; gj<t.size()-1; gj++) s+=t[gj];
		}

		if(s[0]=='[') 
		{
			count=1;
			if(s=="["){k++; s=block.at(k);}
			else{t=s; s=""; for(int gj=1; gj<t.size(); gj++) s+=t[gj];}
		}

		if(count==1)
		{
			//expect tree
			strippedtree="";

			//cout<<"partitions tree "<<s<<endl;
			if(!allAinB(s,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[SITES]", name, "", "The tree name should only contain  ,.;()-_ and alpha-numeric characters.",""); {if(breakonerror) return -1;} }		

			treename=s;	
			treepos=checkthere(treename,totaltreenames);
			
			//cout<<totaltaxa.size()<<" W "<<totaltreenames.size()<<endl;

			if(treepos==-1) 
			{controlerrorprint2("[PARTITIONS]", name, "", "No tree named \""+treename+"\" has been defined.","");  {if(breakonerror) return -1;} }
			else
			{
				Tree* testtree=&(totaltrees.at(treepos));
		
				tree=(*testtree).tree;
				for(int qw=0; qw<tree.size(); qw++) {char c=tree[qw]; if(c==')' || c=='(' || c==',' || c==';') strippedtree+=c;}
				treeposvec.push_back(treepos);

				treenamevec.push_back(treename);
				ntaxavec.push_back((*testtree).ntaxa);

				isitrandom=(*testtree).treerandom;  if(isitrandom) overallrandom=true;
//				isitrandomvec.push_back(isitrandom);
		
				
		
				// if tree random it cannot be used with branch class .
				// else for this case strip it down to check branchclass tree.
			}
		}
		else if(count==2)
		{
			//expect model/branchclass/siteclass name
			
			//cout<<"partitions mbs name "<<s<<endl;
			if(!allAinB(s,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[SITES]", name, "", "The model/branch class/site class name should only contain  \n,.;()-_ and alpha-numeric characters.",s); {if(breakonerror) return -1;} }		

			mbname=s;
			mbspos=checkthere(mbname,totalmodelnames);
			
			if(mbspos==-1) 
			{
				mbspos=checkthere(mbname,totalbranchesnames); 
				
				if(mbspos==-1) 
				{
					mbspos=checkthere(mbname,totalsitenames); 
				
					if(mbspos==-1) 
					{
						controlerrorprint2("[PARTITIONS]", name, "", "You specified \""+mbname+"\" in this partitions class.\nNo model/branch/sites class with this name has been defined.",""); 
						{if(breakonerror) return -1;} 
					}
					else 
					{
						siteclass* s=&(totalsites.at(mbspos));

						if(isitrandom && (*s).therearebranches) 
						{controlerrorprint2("[PARTITIONS]", name, "", "\""+treename+"\" is a random tree.\n\""+mbname+"\" is a sites class that contains at least one branch class.\nThis combination is not allowed.  Guide tree must be user specified.","");  {if(breakonerror) return -1;} }
					
						else 
						{
							mbstype="sites";			
							geneticcode=(*s).geneticcode;
							string strippedtree2=(*s).mbtree;
							bool anerror=false;

							if(!(*s).therearebranches) {}
							else if(strippedtree.size()!=strippedtree2.size()) anerror=true;
							else for(int we=0; we<strippedtree.size(); we++)
							{
								char c1=strippedtree[we]; char c2=strippedtree2[we];
								if(c1!=c2) {anerror=true; break;}
							}

							if(anerror)	{controlerrorprint2("[PARTITIONS]", name, "", "Guide tree (1) does not match the tree used in \nsites class "+mbname+" (2).\n1) "+strippedtree+"\n2) "+strippedtree2+"\nThe trees must have their taxa/branches written in identical order.",""); {if(breakonerror) return -1;} }

						}
					}
				
				}
				else if(isitrandom) {controlerrorprint2("[PARTITIONS]", name, "","\""+treename+"\" is a random tree.\n\""+mbname+"\" is a branches class.\nThis combination is not allowed.  Guide tree must be user specified.","");  {if(breakonerror) return -1;} }

				else 
				{
					mbstype="branches";		
					branchclass* b=&(totalbranches.at(mbspos));

					geneticcode=(*b).geneticcode;
					string strippedtree2=(*b).baretree;

					bool anerror=false;
					if(strippedtree.size()!=strippedtree2.size()) anerror=true;
					else for(int we=0; we<strippedtree.size(); we++)
					{
						char c1=strippedtree[we]; char c2=strippedtree2[we];
						if(c1!=c2) {anerror=true; break;}
					}
					if(anerror)	{controlerrorprint2("[PARTITIONS]", name, "", "Guide tree (1) does not match the tree used in \nbranches class "+mbname+" (2).\n1) "+strippedtree+"\n2) "+strippedtree2+"\nThe trees must have their taxa/branches written in identical order.",""); {if(breakonerror) return -1;} }
					
				}

			}
			else 
			{
				mbstype="model";
				model* m=&(totalmodels.at(mbspos));
				geneticcode=(*m).geneticcode;
			}
			mbstypevec.push_back(mbstype);
			mbsposvec.push_back(mbspos);
			mbnamevec.push_back(mbname);
			geneticcodevec.push_back(geneticcode);
		}
		else if(count==3)
		{
			// expect either a root length, for use in creating root sequence, or a filename containing a root sequence to be read in.

			if(allAinB(s,"1234567890")) 
			{	
				// i.e. entry is numerical, it is a root length
				rootfilenamevec.push_back("NONE");
				rootseqtxtvec.push_back("CREATE");
				rootlength=atoi(s.c_str());	

				if(rootlength==0) { controlerrorprint2("[PARTITIONS]", name, "", "The root length must be an integer greater than zero.",s); {if(breakonerror) return -1;} }
	
				rootlengthvec.push_back(rootlength);

			}
			else
			{			
				// get root sequence
				ifstream if1;
				if1.open(s.c_str());
				if(!if1.good()){ controlerrorprint2("[PARTITIONS]", name, "", "1) This entry is not a positive integer representing a root length.\n2) No root sequence file with this name exists in the program directory.\nPlease format partition sections of the control file as described in the manual",s); {if(breakonerror) return -1;} }		
				else rootfilenamevec.push_back(s);
	
				string test="TCAGtcag";
				if(type==2) test="ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv";
				char c='\n';

				while(if1.good())
				{
					if(AinB(c,test) && !AinB(c,"\n\r\t ")) rootseqtxt+=c;
					else 
					{
						if(!AinB(c,"\n\r\t "))
						{
							string bit="NUCLEOTIDE";
							if(type==2) bit="AMINOACID";
							if(type==3) bit="CODON";
							string s; s+=c;
							controlerrorprint2("[PARTITIONS]", name, "", bit+" root sequence in file can only contain following letters:\n"+test,s); {if(breakonerror) return -1;} 
						}
					}
					c=if1.get();
				}
				
				rootseqtxtvec.push_back(rootseqtxt);
				rootlengthvec.push_back(rootseqtxt.size());

			} // end of else for if rootseqfilename 
						
		} // end of else bracket for 3rd element.

		/*
				// old way of doing step 3 in two steps....... changed to make root sequence creation the default

		else if(count==3)
		{
			if(allAinB(s,"1234567890")) rootlength=atoi(s.c_str());		
			else { controlerrorprint2("[PARTITIONS]", name, "", "The root length must be an integer greater than zero.",s); {if(breakonerror) return -1;} }
			if(rootlength==0) { controlerrorprint2("[PARTITIONS]", name, "", "The root length must be an integer greater than zero.",s); {if(breakonerror) return -1;} }

			rootlengthvec.push_back(rootlength);
		}
		else if(count==4)
		{
			
			//int ans=checkthere(s,totalmodelnames);
				
			//if(ans!=-1) 
			//{
			//	model* m=&(totalmodels.at(ans));
			//	int geneticcode2=(*m).geneticcode;

			//	if(geneticcode!=geneticcode2)
			//	{ 
			//		string type;
			//		if(mbstype=="model") type="Model";
			//		if(mbstype=="branches") type="Branch class";
			//		if(mbstype=="sites") type="Sites class";
					
			//		stringstream ss1;ss1<<geneticcode; string sss1=ss1.str();
			//		stringstream ss2;ss2<<geneticcode2;string sss2=ss2.str();

			//		controlerrorprint2("[PARTITIONS]", name, "", type+" \""+mbname+"\" has genetic code "+sss1+"\nbut root model \""+s+"\" has genetic code "+sss2,""); {if(breakonerror) return -1;} 
			//	}
			

			//	rootfilenamevec.push_back("NONE");
			//	rootmodelpos.push_back(ans);
			//	rootseqtxtvec.push_back("CREATE");
			//}
			//else
			
			if(s=="CREATE")
			{
				rootfilenamevec.push_back("NONE");
				//rootmodelpos.push_back(ans);
				rootseqtxtvec.push_back("CREATE");
			}
			else
			{			
				// get root sequence
				ifstream if1;
				if1.open(s.c_str());
				if(!if1.good()){ controlerrorprint2("[PARTITIONS]", name, "", "No model with this name has been specified.\nNo root sequence file with this name exists in the program directory.",s); {if(breakonerror) return -1;} }		
				else rootfilenamevec.push_back(s);
	
				string test="TCAGtcag";
				if(type==2) test="ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv";
				char c='\n';

				while(if1.good())
				{
					if(AinB(c,test) && !AinB(c,"\n\r\t ")) rootseqtxt+=c;
					else 
					{
						if(!AinB(c,"\n\r\t "))
						{
							string bit="NUCLEOTIDE";
							if(type==2) bit="AMINOACID";
							if(type==3) bit="CODON";
							string s; s+=c;
							controlerrorprint2("[PARTITIONS]", name, "", bit+" root sequence in file can only contain following letters:\n"+test,s); {if(breakonerror) return -1;} 
						}
					}
					c=if1.get();
				}
				
				rootseqtxtvec.push_back(rootseqtxt);

			} // end of else for if rootseqfilename
						
		} // end of else bracket for 4th element.

		*/




	} //	end of for(int k=1; k<mymy; k++) bracket

	int t1=ntaxavec.at(0), t2;
	for(int ft=1; ft<ntaxavec.size(); ft++) 
	{
		// this bracket simply compares guide trees being used by different partitions in a simulation.
		// They can be different but they need to have the same number of external branches

		t2=t1; t1=ntaxavec.at(ft);

		if(t2!=t1)
		{
			string n1=treenamevec.at(ft);
			string n2=treenamevec.at(ft-1);
			stringstream t1ss,t2ss;
			t1ss<<t1; t2ss<<t2;
			string t1s=t1ss.str();
			string t2s=t2ss.str();
			controlerrorprint2("[PARTITIONS]", name, "", "Guide tree \""+n1+"\" has "+t1s+" taxa.\nGuide tree \""+n2+"\" has "+t2s+" taxa.\nGuide trees in different partitions can differ but the \ntotal number of external branches must be constant.",""); 
			{if(breakonerror) return -1;} 
		}
	}
	
	totalpartitions.push_back(partitionclass(name,rootlengthvec, rootseqtxtvec, rootfilenamevec, mbsposvec, mbstypevec, mbnamevec, rootmodelpos, geneticcodevec,treeposvec,treenamevec,t1,overallrandom));
//	totalpartitions.push_back(partitionclass(treeposvec,treenamevec,rootlengthvec,mbnamevec,mbstypevec,mbsposvec,rootseqtxtvec));


	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class evolve
{
	// real simple class!
	// this is the one that links it all together.
	// Each instance of the class evolve is a standalone simulation of a partition that contains a branch class, or a model or so on....

	public:
	int reps;
	int partitionpos;
	string filenamestub;
	
	evolve(int myreps, int mypartitionpos, string myfilenamestub)
	{
		reps=myreps;						// number of replicates
		partitionpos=mypartitionpos;		// position in totalpartitions of the partition in question
		filenamestub=myfilenamestub;		// filenamestub which is used to create output filenames
	}
};

vector<evolve> totalevolve;		// storage for the evolve class instances.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int dealwithevolve(vector<string> &block)		
{
	//	this function pre-processes the evolve blocks for the class evolve.   pretty self explanatory,.... at least it should be.

	int size=block.size()-1, partitionpos, reps;

	string filenamestub;

	if(size==0) {controlerrorprint2("[EVOLVE]","","","The block is empty.",""); {if(breakonerror) return -1;} }
//	cout<<"size "<<size<<endl;
	if(size%3!=0) {controlerrorprint2("[EVOLVE]","","","The Evolve block must contain a multiple of 3 statements:\nPartitionClassName1  NumberOfReps1  FilenameStub1\nPartitionClassName2  NumberOfReps2  FilenameStub2\nand so on..."/*\nIf you have followed this instruction then the error may be caused\nas there is no white space after the last entry in the control file!\nPlease put a space or a new line after your last entry in this case"*/,""); {if(breakonerror) return -1;} }

	string s1,s2,s3,t1,t2,t3;
	int ans;
	for(int j=2; j<size; j=j+3) 
	{
//		stringstream tt1,tt2,tt3; tt1<<j-2; tt2<<j-1; tt3<<j; t1=tt1.str(); t2=tt2.str(); t3=tt3.str();

		s1=block.at(j-2); s2=block.at(j-1); s3=block.at(j);

		if(!allAinB(s1,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;()-_")) { controlerrorprint2("[EVOLVE]", "", "", "Partition class names should only contain  ,.;()-_ \nand alpha-numeric characters.",s1); {if(breakonerror) return -1;} }		
		ans=checkthere(s1,totalpartitionnames);
		if(ans==-1) {controlerrorprint2("[EVOLVE]","","","No partitions block called "+s1+" has been defined.",""); {if(breakonerror) return -1;} }
		else partitionpos=ans;

		if(!allAinB(s2,"123456780")) { controlerrorprint2("[EVOLVE]", "", "", "Number of repetitions must be an integer value",s2); {if(breakonerror) return -1;} }		
		else reps=atoi(s2.c_str());				
		
		if(!allAinB(s3,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.;/()-_")) { controlerrorprint2("[EVOLVE]", "", "", "Filename stubs should only contain  ,.;()-_ and alpha-numeric characters.",s3); {if(breakonerror) return -1;} }		
		else 
		{
			filenamestub=s3;		
			bool error=false;
			for(int i=0; i<totalevolve.size(); i++)
			{
				if( (totalevolve.at(i)).filenamestub==filenamestub ) error=true;
			}
			if(error) {controlerrorprint2("[EVOLVE]", "", "", "Several lines in the [EVOLVE] block contain the same filename stub.\nThis will cause output files to be overwritten - please correct.",s3); {if(breakonerror) return -1;} }		
	
		}

		// cout<<partitionpos<<" "<<reps<<" "<<filenamestub<<endl;
		totalevolve.push_back(evolve(reps,partitionpos,filenamestub));

	}

	return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

string paupstart, paupmiddle, paupend;

int docontrol()
{		
	// this is the main function that processes the control file and calls the other relevant processing functions

	int isthereanerror=0;
	ifstream if1;
	if1.open(masterfilename.c_str());

	char c1='q', c2='q'; //,c2,c3,c4,c5,c6;
	bool notwhitespace=true;
	string s, newfilename=masterfilename;
	vector<string> sv;


	// get command blocks if they exist
	ifstream paups, paupm, paupf;

	paups.open("paupstart.txt");
	paupm.open("paupmiddle.txt");	
	paupf.open("paupend.txt");

	while(paups.good()) {char c=paups.get(); if(c!='ÿ') paupstart+=c; }//cout<<c<<endl;}   
	while(paupm.good()) {char c=paupm.get(); if(c!='ÿ') paupmiddle+=c; }//cout<<c<<endl;}
	while(paupf.good()) {char c=paupf.get(); if(c!='ÿ') paupend+=c; }//cout<<c<<endl;}

	if(paupstart!="") paupstart+="\n\n";
	if(paupmiddle!="") paupmiddle+="\n"; paupmiddle+="\n";
	if(paupend!="") paupend+="\n\n";

	bool writeon=true;
	string s1,s2,s3;
	if(!if1.good()) 
	{
		bool waste=true;
		for(int qa=0; qa<10; qa++)
		{
			// this prompts for a different control file name if there is no file sharing name of masterfilename that is hardwired in code.

			controlerrorprint2("CONTROL FILE","","","There is no control file. INDELible was looking for file named:\n" +newfilename+"\nPlease enter the correct filename.","");  
			
			cin>>newfilename;

			if1.open(newfilename.c_str());

			if(if1.good()) {waste=false; break;}
		}
		
		if(waste)
		{
			controlerrorprint2("CONTROL FILE","","","There is no control file. INDELible was looking for file named:\n" +newfilename+"\n10 failed attempts - please check your filenames and re-run INDELible.","");  
			isthereanerror=-1;
			return -1;
		}
	}
	

	/*
		// old way of reading in control file 

		while(if1.good())
		{
			// ignore comment lines in command file
			getline(if1,s1); s3="";
			if(s1[0]!='/' && s1[1]!='/') 
			{
				for(int rd=0; rd<s1.size(); rd++) 
				{
					if(s1[rd]=='/' ) break; //&& (s1[rd+1]=='/' || s1[rd+1]=='\n' || s1[rd+1]=='\r')) break;
					else s3+=s1[rd];
				}
				s2+=s3; s2+='\n';
			}
		}

		for(int hb=0; hb<s2.size(); hb++)
		{
			c2=c1; c1=s2[hb];  
			//if(c2=='/' && c1=='*') writeon=false;
			
			if(writeon)
			{	
				if(AinB(c1," \n\r\t")) {if(s!="") sv.push_back(s); s="";}
				else{s+=c1;}
			}

		}
		sv.push_back(s);
	*/


	vector<char> temp1, temp2, temp3;

	temp1.push_back(' ');   // to prevent crashes at *1 below
	while(if1.good()) {char c=if1.get(); temp1.push_back(c); originalcontrol.push_back(c);}          // parse control file character by character into temp1
	temp1.push_back('\n');  // to prevent crashes at *2 below, and also if control file does not end in new line..

	bool WR=true;

	int mymycount=0;

	for(int th1=1; th1<temp1.size()-1; th1++)
	{
		// this for loop will weed out any commented text between /* .... */ in the control file (from temp1 to temp2)

		if(temp1.at(th1) == '/' && temp1.at(th1+1)=='*') {WR=false; mymycount++;}  ///// *1

		if(WR) temp2.push_back(temp1.at(th1));

		if(temp1.at(th1-1)=='*' && temp1.at(th1) == '/') {WR=true; mymycount--;}   ///// *2

		if(mymycount<0)  {controlerrorprint2("CONTROL FILE","","","There is a '*/' statement without a preceding\nand matching '/*' statement! ",""); return -1;}

	}	
	
	temp2.push_back(' ');   WR=true;

	for(int th2=0; th2<temp2.size()-1; th2++) 
	{
		// this for loop should get rid of any remaining lines that begin with double slash //
		// or get rid of the end of any lines that have // ends to them

		char q1=temp2.at(th2);

		if(q1=='/' && temp2.at(th2+1)=='/') WR=false;

		if(WR) temp3.push_back(q1);

		if(q1=='\n' || q1=='\r') WR=true;
	}

	
	// re-sorts information into "word size" string blocks
	for(int th3=0; th3<temp3.size(); th3++) 
	{
		c1=temp3.at(th3);  
		                                     //  again this is to keep backward compatibility with all my old control files
		if(AinB(c1," \n\r\t")) {if(s!="") {if(s=="[basefreq]") sv.push_back("[statefreq]"); else sv.push_back(s);} s="";}
		else{s+=c1;}	
	}
	sv.push_back(s);



	// archaic little function to put "[user]" back in user-defined tree code -----> easier than re-writing the whole of the dealwithtrees function.
	// BUT it doesn't put it back in if it is already there! ---> that means I can use my mountaind of old "validation" control files without having to amend them.

	bool weareonbaby=false;

	for(int sd=0; sd<sv.size()-1; sd++)
	{
		string sq1=sv.at(sd), sq2=sv.at(sd+1);  
		if(sq1=="[TYPE]" || sq1=="[SETTINGS]" || sq1=="[MODEL]" || sq1=="[SITES]" || sq1=="[BRANCHES]" || sq1=="[PARTITIONS]" || sq1=="[EVOLVE]" || sq1== "[BRANCHES*]") weareonbaby=false;
		if(sq1=="[TREE]") weareonbaby=true;
		
		if(weareonbaby && sq2[0]=='(')
		{
			if(sq1 != "[user]" )  {sv.insert(sv.begin()+sd+1,"[user]"); sd++;}
			weareonbaby=false;
		}
 
	}



	// archaic little function to put "M3" back in codon model code -----> easier than re-writing the whole of the type 3 [submodel] command in the dealwithmodel function.
	// BUT it doesn't put it back in if it is already there! ---> that means I can use my mountaind of old "validation" control files without having to amend them.

	if(type==3)
	{
		bool weareonbaby=false;

		for(int sd=0; sd<sv.size()-1; sd++)
		{
			string sq1=sv.at(sd), sq2=sv.at(sd+1);  
			if(sq1=="[submodel]") weareonbaby=true;
			else if(sq1[0]=='[') weareonbaby=false;
			
			if(weareonbaby && sq2[0]=='(')
			{
				if(sq1 != "M3" && sq1 != "M2" && sq1 != "M1" && sq1 != "M0" && sq1 != "ECMrest" && sq1 != "ECMunrest")  {sv.insert(sv.begin()+sd+1,"M3"); sd++;}
				weareonbaby=false;
			}
 
		}

	}	
			
	
		bool settingsblock,modelblock,sitesblock,branchesblock,partitionsblock,evolveblock,doneanymodelsyet;
		settingsblock=modelblock=sitesblock=branchesblock=partitionsblock=evolveblock=doneanymodelsyet=false;
		
		string lastbit;
		vector<string> block;
		int well=0;

		if(sv.at(0)!="[TYPE]") { controlerrorprint2("CONTROL FILE","","","First entry in control file must be a [TYPE] command.",""); return -1;}
		else if(sv.at(1)=="NUCLEOTIDE")		type=1;
		else if(sv.at(1)=="AMINOACID")		type=2;
		else if(sv.at(1)=="CODON")			type=3;
		else {isthereanerror=-1; 	controlerrorprint2("[TYPE]","","","First value for [TYPE] command must be NUCLEOTIDE, AMINOACID, or CODON",""); return -1;}



		
	//this loop will go through and change model names (new system. e.g "HKY") to numbers (old system, e.g. "3") but will ignore old system if used.
	
	string codonmodelnames[16]={"M0","M1","M2","M3","M4xxx","M5xxx","M6xxx","M7xxx","M8xxx","M9xxx","M10xxx","M11xxx","M12xxx","M13xxx","ECMrest","ECMunrest"};

	string nucleotidemodelnames[17]={"JC","F81","K80","HKY","TrNef","TrN","K81","K81uf","TIMef","TIM","TVMef","TVM","SYM","GTR","F84ef","F84","UNREST"};

	string aminoacidmodelnames[17]={"Poisson","JTT","JTT-dcmut","Dayhoff","Dayhoff-dcmut","WAG","mtMAM","mtART","mtREV","rtREV","cpREV","Vt","Blosum","LG","HIVb","HIVw","USER"};

	string option3="M0, M1, M2, M3, ECMrest, ECMunrest.";

	string option1="JC, F81, K80, HKY, TrNef, TrN, K81, K81uf, TIMef,\nTIM, TVMef, TVM, SYM, GTR, F84ef, F84, UNREST.";

	string option2="Poisson, JTT, JTT-dcmut, Dayhoff, Dayhoff-dcmut, WAG, mtMAM,\nmtART, mtREV, rtREV, cpREV, Vt, Blosum, LG, HIVb, HIVw, USER.";


	vector<string> replacenames;   string myname, mybit=sv.at(1), optionstring;

	if(type==3)			{for(int pg=0; pg<16; pg++) replacenames.push_back(codonmodelnames[pg]);       optionstring=option3;}
	else if(type==2)	{for(int pg=0; pg<17; pg++) replacenames.push_back(aminoacidmodelnames[pg]);   optionstring=option2;}
	else if(type==1)	{for(int pg=0; pg<17; pg++) replacenames.push_back(nucleotidemodelnames[pg]);  optionstring=option1;}

	for(int qw=1; qw<sv.size(); qw++)
	{
		string test=sv.at(qw), oldtest=sv.at(qw-1);

		if(oldtest=="[MODEL]") myname=test;

		if(oldtest=="[submodel]")
		{
			if(allAinB(test,"0123456789.")) continue;
			
			bool error=true;
			stringstream sss; string ss;
			for(int pg=0; pg<replacenames.size(); pg++) if(replacenames.at(pg)==test) {sss<<pg; ss=sss.str(); sv.at(qw)=ss; error=false; break;}

			if(type==2 && error) 
			{
				ifstream ig1; ig1.open(test.c_str()); 

				if(!(ig1.good())) {controlerrorprint2("[MODEL]", myname, "[submodel]", test+"\nis not a model name or number and no file of this name exists.\nFor protein models, the entry after a [submodel] command must be\nan integer, model name or a filename.",""); {if(breakonerror) return -1;} }
				else error=false;
				
								
			}

			if(error) 
			{
				string mystring="\nThere is no "+mybit+" substitution model named: "+test+".\n\nYour options are:\n\n"+optionstring+"\n\nor the numerical counterparts. Please consult manual.";
				if(type==2) mystring+="\nThere is also no file named "+mybit+" in the INDELible directory.";
				controlerrorprint2("[MODEL]",myname,"[submodel]",mystring,""); 
				return -1;
			}

		}
	}
/*
		string blah=sv.back(),blah2;
		if(blah.size()!=1)
		{
			//to prevent parsing error if there is no space, horizontal tab or new line after last entry in control file.
			for(int ghg=0; ghg<blah.size()-1; ghg++) blah2+=blah[ghg];

			sv.back()=blah2; sv.push_back("\n");
		}
*/		

		// checks which algorithm people want to use
		if(sv.at(2)=="1") oldmethod=false;		//old model
		else if(sv.at(2)=="2") oldmethod=true; //new model

		else {isthereanerror=-1; 	controlerrorprint2("[TYPE]","","","Second value for [TYPE] command must 1 or 2 for choice of algorithm",""); return -1;}


		int mytype=-1;

		vector<string> rough, settingsV, evolveV;
		vector<vector<string> > modelsV, sitesV, branchesV, treesV, partitionsV;
		vector<bool> branchesCV;


	//	string thisbit=sv.at(2);
		string thisbit;
		string MAIN[8]={"[SETTINGS]","[MODEL]","[SITES]","[BRANCHES]","[TREE]","[PARTITIONS]","[BRANCHES*]","[EVOLVE]"};

		bool doit=true;
		int numberofsettingsblocks=0;
		int numberofmodelsblocks=0;

	/*		
		if(well!=-1)
		{
			if(thisbit=="[SETTINGS]") mytype=0;
			else if (thisbit=="[MODEL]") mytype=1;
			else cout<<"ERROR: After the [TYPE] setting should come a [SETTINGS] or [MODEL] block."<<endl;
		}

		if(mytype==0 || mytype==1) 
	*/
		if(sv.at(3)!="[MODEL]" && sv.at(3)!="[SETTINGS]") {isthereanerror=-1;	controlerrorprint2("[TYPE]","","","There is some text after the [TYPE] command before the first block.\n Is this a mistake? Was expecting a [MODEL] or [SETTINGS] block",sv.at(3)); return -1;}


		if(isthereanerror!=-1) for(int i=3; i<sv.size(); i++) 
		{
			// if there is no error then this function collects together blocks of same type etc.

			thisbit=sv.at(i);

			doit=true;

			for(int j=0; j<8; j++) 
			{
				if(thisbit==MAIN[j]) 
				{
					
					if(mytype==0)		{settingsV=rough; numberofsettingsblocks++;}
	  				else if(mytype==1)	{modelsV.push_back(rough); numberofmodelsblocks++;}
					else if(mytype==2)	sitesV.push_back(rough);
					else if(mytype==3)	{branchesV.push_back(rough); branchesCV.push_back(true);}
					else if(mytype==4)	treesV.push_back(rough);
					else if(mytype==5)	partitionsV.push_back(rough);
					else if(mytype==6)	{branchesV.push_back(rough); branchesCV.push_back(false);}

					mytype=j;
					rough.clear();
					doit=false;
					break;
				}
				

			}

			if(doit) rough.push_back(thisbit);	
		}

		if(well!=-1)
		{
			// final error checking

			if(numberofsettingsblocks>1) {controlerrorprint2("[SETTINGS]","","","ERROR: There was more than one [SETTINGS] block.","");   return -1;}
			if(numberofmodelsblocks==0) {controlerrorprint2("[MODEL]","","","ERROR: No [MODEL] blocks were specified.","");   return -1;}


			if(mytype!=7) {controlerrorprint2("[EVOLVE]","","","ERROR: Control file should end with an [EVOLVE] block","");    return -1;}
			else evolveV=rough;
			rough.clear();
			sv.clear();
		}

	
	int p; //,q;

	// settings blocks first!
	if(isthereanerror!=-1 && numberofsettingsblocks>0) dealwithsettings(settingsV);

	// runs through commands in correct order, and skips if an error occurs anywhere
	for(p=0; p<modelsV.size();     p++) { if(isthereanerror==-1 ) return -1; else isthereanerror=dealwithmodel(modelsV.at(p)); }
	for(p=0; p<branchesV.size();   p++) { if(isthereanerror==-1 ) return -1; else isthereanerror=dealwithbranches(branchesV.at(p), branchesCV.at(p)); }	
	for(p=0; p<sitesV.size();      p++) { if(isthereanerror==-1 ) return -1; else isthereanerror=dealwithsites(sitesV.at(p)); }
	for(p=0; p<treesV.size();      p++) { if(isthereanerror==-1 ) return -1; else isthereanerror=dealwithtrees(treesV.at(p)); }
	for(p=0; p<partitionsV.size(); p++) { if(isthereanerror==-1 ) return -1; else isthereanerror=dealwithpartitions(partitionsV.at(p)); }
	if(isthereanerror==-1 ) return -1; else isthereanerror=dealwithevolve(evolveV);


	
//	cout<<"HERE"<<endl;
	

//	cout<<"HERE"<<endl;
/*
	// for testing model copying and assigning etc

	for(p=0; p<totalmodels.size(); p++)    
	{
		model m=totalmodels.at(p);
		cout<<endl<<"MODEL "<<m.name<<endl<<endl;

		cout<<m.modelnumber<<" modelnumber"<<endl;
		cout<<m.geneticcode<<" geneticcode"<<endl;
		cout<<m.indeltosub<<" indeltosub"<<endl;
		cout<<m.instodel<<" instodel"<<endl;
		cout<<m.inlength<<" inlength"<<endl;
		cout<<m.dellength<<" dellength"<<endl;

		cout<<m.alpha<<" alpha"<<endl;
		cout<<m.pinv<<" pinv"<<endl;
		cout<<m.ngamcat<<" ngamcat"<<endl;

		vector<double> bits; int gg;
		double sum=0;

		sum=0; bits=m.rootbasefreqs; cout<<"rootbasefreqs ";
		for(gg=0;gg<bits.size();gg++) {sum+=bits.at(gg); cout<<bits.at(gg)<<" ";}//<<endl;}
		cout<<sum<<endl;
		sum=0; bits=m.insertfreqs; cout<<"insertfreqs ";
		for(gg=0;gg<bits.size();gg++) {sum+=bits.at(gg); cout<<bits.at(gg)<<" ";}//endl;}
		cout<<sum<<endl;
		sum=0; bits=m.basefreqs; cout<<"basefreqs ";
		for(gg=0;gg<bits.size();gg++) {sum+=bits.at(gg); cout<<bits.at(gg)<<" ";}//endl;}
		cout<<sum<<endl;

	}
*/

	/*
	// for printing of information

	for(q=0; q<settingsV.size(); q++)	cout<<"settings "<<q<<" "<<settingsV.at(q)<<endl;
   
	for(p=0; p<modelsV.size(); p++)		for(q=0; q<(modelsV.at(p)).size(); q++) cout<<"models "<<p<<" "<<q<<" "<<(modelsV.at(p)).at(q)<<endl;
	
	for(p=0; p<sitesV.size(); p++)		for(q=0; q<(sitesV.at(p)).size(); q++) cout<<"sites "<<p<<" "<<q<<" "<<(sitesV.at(p)).at(q)<<endl;
    
	for(p=0; p<branchesV.size(); p++)	for(q=0; q<(branchesV.at(p)).size(); q++) cout<<"branches "<<p<<" "<<q<<" "<<(branchesV.at(p)).at(q)<<endl;
	
	for(p=0; p<partitionsV.size(); p++)	for(q=0; q<(partitionsV.at(p)).size(); q++) cout<<"partitions "<<p<<" "<<q<<" "<<(partitionsV.at(p)).at(q)<<endl;

	for(q=0; q<evolveV.size(); q++)		cout<<"evolve "<<q<<" "<<evolveV.at(q)<<endl;
    */

	printf("\n");
  return isthereanerror;
}				


// include code below if you want to run standalone tests of the code in here without linking to main skeleton program
//double main(int argc, char* argv[])
//{
//	return 0;
//}	
	
