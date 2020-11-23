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




// almost all the functions in this file were borrowed from PAML (by Ziheng Yang).

// the code is for calculation of P(t) from Q

// for calculation of random gamma rate variation

// and for producing random trees through the birth/death process.


#include <math.h>
#include "randoms.cpp"
#include <string>
#include <sstream>
#include <vector>
#include "float.h"
#include <iostream>
#include <fstream>
#include <time.h>

#define square(a) ((a)*(a))
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))

extern double myrand;
extern int idum;
		
using namespace std;

#define FOR(i,n) for(i=0; i<n; i++)

int identity (double x[], int n)
{ int i,j;  FOR (i,n)  { FOR(j,n)  x[i*n+j]=0;  x[i*n+i]=1; }  return (0); }

int zero (double x[], int n)
{ int i; FOR (i,n) x[i]=0; return (0);}

int NPMatUVRoot=0; 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int matby (double a[], double b[], double c[], int n,int m,int k)
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b  */
{
   int i,j,i1;
   double t;
   FOR (i,n)  FOR(j,k) {
      for (i1=0,t=0; i1<m; i1++) t+=a[i*m+i1]*b[i1*k+j];
      c[i*k+j] = t;
   }
   return (0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//int matexp (double Q[], double t, int n, /*int TimeSquare, */double space[])
vector<vector<double> > matexp (vector<vector<double> > Qvec, vector<double> &basefreqs,  double t)
{

	// modified by WF

	int TimeSquare=30,  n=basefreqs.size(), j;
	vector<vector<double> > Pmat;
	vector<double> row;
/* 
Method of calculating P(t) from Q using the repeated matrix squaring that is implemented in the PAML package (Yang, 1997). 

  Yang, Z. (1997) PAML: a program for package for phylogenetic analysis by maximum likelihood. CABIOS 15: 555-556.

  Original notes from ZY follows:
  
-----------------------------------

  This calculates the matrix exponential P(t) = exp(t*Q).
   Input: Q[] has the rate matrix, and t is the time or branch length.
          TimeSquare is the number of times the matrix is squared and should 
          be from 5 to 31.
   Output: Q[] has the transition probability matrix, that is P(Qt).
   space[n*n]: required working space.

      P(t) = (I + Qt/m + (Qt/m)^2/2)^m, with m = 2^TimeSquare.

   T[it=0] is the current matrix, and T[it=1] is the squared result matrix,
   used to avoid copying matrices.
   Use an even TimeSquare to avoid one round of matrix copying.
*/




		int it, i;

		double Q[4096];  double space[8320]; double pi[64];


	   	for(i=0; i<n; i++) 
		{
			pi[i]=basefreqs.at(i);

			for(j=0; j<n; j++) Q[(i*n)+j]=(Qvec.at(i)).at(j); 
		}
	   


	   double *T[2];

	//   if(TimeSquare<2 || TimeSquare>31) error2("TimeSquare not good");

	   T[0]=Q; T[1]=space;
	   for(i=0; i<n*n; i++)  T[0][i] = ldexp( Q[i]*t, -TimeSquare );

	   matby (T[0], T[0], T[1], n, n, n);
	   for(i=0; i<n*n; i++)  T[0][i] += T[1][i]/2;
	   for(i=0; i<n; i++)  T[0][i*n+i] ++;

	   for(i=0,it=0; i<TimeSquare; i++) {
		  it = !it;
		  matby (T[1-it], T[1-it], T[it], n, n, n);
	   }

	   if(it==1) for(i=0;i<n*n;i++) Q[i]=T[1][i];

		for(i=0; i<n; i++) 
		{
			for(j=0; j<n; j++) row.push_back(Q[(i*n)+j]);
			for(j=1; j<n; j++) row.at(j)+=row.at(j-1);
			double diff=row.at(n-1)-1; if(diff<0) diff=-diff;
			if(diff>0.00001) {cout<<"matexp row didn't sum to 1 "<<diff<<endl; for(j=0; j<n;j++) cout<<row.at(j)<<"\t"; cout<<endl;cout<<endl;}
			Pmat.push_back(row);
			row.clear();
		}


   return Pmat;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EigenSort(double d[], double U[], int n)
{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
*/
   int k,j,i;
   double p;

   for (i=0;i<n-1;i++) {
      p=d[k=i];
      for (j=i+1;j<n;j++)
         if (d[j] >= p) p=d[k=j];
      if (k != i) {
         d[k]=d[i];
         d[i]=p;
         for (j=0;j<n;j++) {
            p=U[j*n+i];
            U[j*n+i]=U[j*n+k];
            U[j*n+k]=p;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HouseholderRealSym(double a[], int n, double d[], double e[])
{
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix 
   a[n*n] into a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigends), and e[] the off-diagonal.
*/
   int m,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>=1;i--) {
      m=i-1;
      h=scale=0;
      if (m > 0) {
         for (k=0;k<=m;k++)
            scale += fabs(a[i*n+k]);
         if (scale == 0)
            e[i]=a[i*n+m];
         else {
            for (k=0;k<=m;k++) {
               a[i*n+k] /= scale;
               h += a[i*n+k]*a[i*n+k];
            }
            f=a[i*n+m];
            g=(f >= 0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i*n+m]=f-g;
            f=0;
            for (j=0;j<=m;j++) {
               a[j*n+i]=a[i*n+j]/h;
               g=0;
               for (k=0;k<=j;k++)
                  g += a[j*n+k]*a[i*n+k];
               for (k=j+1;k<=m;k++)
                  g += a[k*n+j]*a[i*n+k];
               e[j]=g/h;
               f += e[j]*a[i*n+j];
            }
            hh=f/(h*2);
            for (j=0;j<=m;j++) {
               f=a[i*n+j];
               e[j]=g=e[j]-hh*f;
               for (k=0;k<=j;k++)
                  a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
            }
         }
      } 
      else
         e[i]=a[i*n+m];
      d[i]=h;
   }
   d[0]=e[0]=0;

   /* Get eigenvectors */
   for (i=0;i<n;i++) {
      m=i-1;
      if (d[i]) {
         for (j=0;j<=m;j++) {
            g=0;
            for (k=0;k<=m;k++)
               g += a[i*n+k]*a[k*n+j];
            for (k=0;k<=m;k++)
               a[k*n+j] -= g*a[k*n+i];
         }
      }
      d[i]=a[i*n+i];
      a[i*n+i]=1;
      for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.  
   d[] is the diagonal (eigenvalues), e[] is the off-diagonal
   z[n*n]: as input should have the identity matrix to get the eigen solution of the 
   tridiagonal matrix, or the output from HouseholderRealSym() to get the 
   eigen solution to the original real symmetric matrix.
   z[n*n]: has the orthogonal matrix as output

   Adapted from routine tqli in Numerical Recipes in C, with reference to
   LAPACK fortran code.
   Ziheng Yang, May 2001
*/
   int m,j,iter,niter=30, status=0, i,k;
   double s,r,p,g,f,dd,c,b, aa,bb;

   for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
   for (j=0;j<n;j++) {
      iter=0;
      do {
         for (m=j;m<n-1;m++) {
            dd=fabs(d[m])+fabs(d[m+1]);
            if (fabs(e[m])+dd == dd) break;  /* ??? */
         }
         if (m != j) {
            if (iter++ == niter) {
               status=-1;
               break;
            }
            g=(d[j+1]-d[j])/(2*e[j]);

            /* r=pythag(g,1); */

            if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
            else                r=sqrt(1+g*g);

            g=d[m]-d[j]+e[j]/(g+SIGN(r,g));
            s=c=1;
            p=0;
            for (i=m-1;i>=j;i--) {
               f=s*e[i];
               b=c*e[i];

               /*  r=pythag(f,g);  */
               aa=fabs(f); bb=fabs(g);
               if(aa>bb)       { bb/=aa;  r=aa*sqrt(1+bb*bb); }
               else if(bb==0)             r=0;
               else            { aa/=bb;  r=bb*sqrt(1+aa*aa); }

               e[i+1]=r;
               if (r == 0) {
                  d[i+1] -= p;
                  e[m]=0;
                  break;
               }
               s=f/r;
               c=g/r;
               g=d[i+1]-p;
               r=(d[i]-g)*s+2*c*b;
               d[i+1]=g+(p=s*r);
               g=c*r-b;
               for (k=0;k<n;k++) {
                  f=z[k*n+i+1];
                  z[k*n+i+1]=s*z[k*n+i]+c*f;
                  z[k*n+i]=c*z[k*n+i]-s*f;
               }
            }
            if (r == 0 && i >= j) continue;
            d[j]-=p; e[j]=g; e[m]=0;
         }
      } while (m != j);
   }
   return(status);
}

#undef SIGN

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int eigenRealSym(double A[], int n, double Root[], double work[])
{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return, 
   A has the right vectors and Root has the eigenvalues. 
   work[n] is the working space.
   The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(), 
   and then using the QL algorithm with implicit shifts.  

   Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
   Ziheng Yang, 23 May 2001
*/
   int status=0;
   HouseholderRealSym(A, n, Root, work);
   status=EigenTridagQLImplicit(Root, work, n, A);
   EigenSort(Root, A, n);

   return(status);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int eigenQREV (double Q[], double pi[], int n, 
               double Root[], double U[], double V[], double spacesqrtpi[])
{
/* 
   This finds the eigen solution of the rate matrix Q for a time-reversible 
   Markov process, using the algorithm for a real symmetric matrix.
   Rate matrix Q = S * diag{pi} = U * diag{Root} * V, 
   where S is symmetrical, all elements of pi are positive, and U*V = I.
   space[n] is for storing sqrt(pi).

   [U 0] [Q_0 0] [U^-1 0]    [Root  0]
   [0 I] [0   0] [0    I]  = [0     0]

   Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
*/

   int i,j, inew, jnew, nnew, status;
   double *pi_sqrt=spacesqrtpi, mysmall=1e-100;

   for(j=0,nnew=0; j<n; j++)
      if(pi[j]>mysmall)
         pi_sqrt[nnew++]=sqrt(pi[j]);

   /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */

   if(nnew==n) {
      for(i=0; i<n; i++)
         for(j=0,U[i*n+i] = Q[i*n+i]; j<i; j++)
            U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);

      status=eigenRealSym(U, n, Root, V);
      for(i=0;i<n;i++) for(j=0;j<n;j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
      for(i=0;i<n;i++) for(j=0;j<n;j++)  U[i*n+j] /= pi_sqrt[i];
   }
   else {
      for(i=0,inew=0; i<n; i++) {
         if(pi[i]>mysmall) {
            for(j=0,jnew=0; j<i; j++) 
               if(pi[j]>mysmall) {
                  U[inew*nnew+jnew] = U[jnew*nnew+inew] 
                                    = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
                  jnew++;
               }
            U[inew*nnew+inew] = Q[i*n+i];
            inew++;
         }
      }

      status=eigenRealSym(U, nnew, Root, V);

      for(i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
         Root[i] = (pi[i]>mysmall ? Root[inew--] : 0);
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
         if(pi[i]>mysmall) {
            for(j=n-1,jnew=nnew-1; j>=0; j--)
               if(pi[j]>mysmall) {
                  V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
                  jnew--;
               }
               else 
                  V[i*n+j] = (i==j);
            inew--;
         }
         else 
            for(j=0; j<n; j++)  V[i*n+j] = (i==j);
      }
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
         if(pi[i]>mysmall) {
            for(j=n-1,jnew=nnew-1;j>=0;j--)
               if(pi[j]>mysmall) {
                  U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[inew];
                  jnew--;
               }
               else 
                  U[i*n+j] = (i==j);
            inew--;
         }
         else 
            for(j=0;j<n;j++)
               U[i*n+j] = (i==j);
      }
   }

/*   This routine works on P(t) as well as Q. */
/*
   if(fabs(Root[0])>1e-10 && noisy) printf("Root[0] = %.5e\n",Root[0]);
   Root[0]=0; 
*/
   return(status);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int PMatUVRoot (double P[], double t, int n, double U[], double V[], double Root[])
{
/* P(t) = U * exp{Root*t} * V
*/
   int i,j,k;
   double expt, uexpt, *pP;
   double smallp = 0;

   NPMatUVRoot++;
   if (t<-0.1) printf ("\nt = %.5f in PMatUVRoot", t);
   if (t<1e-100) { identity (P, n); return(0); }
   for (k=0,zero(P,n*n); k<n; k++)
      for (i=0,pP=P,expt=exp(t*Root[k]); i<n; i++)
         for (j=0,uexpt=U[i*n+k]*expt; j<n; j++)
            *pP++ += uexpt*V[k*n+j];

   for(i=0;i<n*n;i++)  if(P[i]<smallp)  P[i]=0;

#if (DEBUG>=5)
      if (testTransP(P,n)) {
         printf("\nP(%.6f) err in PMatUVRoot.\n", t);
         exit(-1);
      }
#endif

   return (0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//int PMatQRev(double Q[], double pi[], double t, int n, double space[])
vector<vector<double> > PMatQRev(vector<vector<double> > Qvec, vector<double> &basefreqs, double t) 
{
/*
 modified by WF

  Calculates P(t) from Q using the eigenvalues and eigenroots of the rate matrix Q (Yang, 1995).  

  Yang, Z. (1995) On the general reversible Markov-process model of nucleotide substitution: a reply to Saccone et al. Journal of Molecular Evolution 41:254-255.


  Original notes from ZY follows:

  -------------------------------

   This calculates P(t) = exp(Q*t), where Q is the rate matrix for a 
   time-reversible Markov process.

   Q[] or P[] has the rate matrix as input, and P(t) in return.
   space[n*n*2+n*2]  n*n*2+n*2=n*2*(n+1)  
*/
	int i, j, n=basefreqs.size();

	vector<vector<double> > Pmat;
	vector<double> row;


			double Q[4096];  double space[8320]; double pi[64];

	   	for(i=0; i<n; i++) 
		{
			pi[i]=basefreqs.at(i);

			for(j=0; j<n; j++) Q[(i*n)+j]=(Qvec.at(i)).at(j); 
		}

		double *U=space, *V=U+n*n, *Root=V+n*n, *spacesqrtpi=Root+n;

		eigenQREV(Q, pi, n, Root, U, V, spacesqrtpi);
		PMatUVRoot(Q, t, n, U, V, Root);
				
		for(i=0; i<n; i++) 
		{
			for(j=0; j<n; j++) row.push_back(Q[(i*n)+j]);
			for(j=1; j<n; j++) row.at(j)+=row.at(j-1);
			double diff=row.at(n-1)-1; if(diff<0) diff=-diff;
			if(diff>0.00001) {cout<<"Pmat row didn't sum to 1 "<<diff<<endl; for(j=0; j<n;j++) cout<<row.at(j)<<"\t"; cout<<endl;cout<<endl;}
			Pmat.push_back(row);
			row.clear();
		}


	return Pmat;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////






/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
//            ALL THE CODE BELOW IS FOR FOR PRODUCING RANDOM LABELLED HISTORIES                //
//			  (AND BRANCH LENGTHS FROM BIRTH DEATH PROCESS) UNTIL STATED OTHERWISE             //
//            IT HAS BEEN USED WITH KIND PERMISSION FROM ZIHENG YANG'S PAML SUITE              //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

				int LASTROUND;


				#define FPN(file) fputc('\n', file)
				#define F0 stdout
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				#define NS            5000
				#define NBRANCH       (NS*2-2)
				#define MAXNSONS      20
				#define LSPNAME       50
				#define NCODE         64
				#define NCATG         40
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////

				void error2 (string message)
				{ cout<<"\nError thrown by PAML code:\n"<<message<<endl; exit(-1); }

				static unsigned int z_rndu=137;
				static unsigned int w_rndu=123456757;

				void SetSeed (unsigned int seed)
				{
				   if(sizeof(int) != 4) 
					  puts("oh-oh, we are in trouble.  int not 32-bit?");
				   z_rndu=170*(seed%178)+137;
				   w_rndu = seed*127773;
				}


				double rndu (void)
				{
					/* From Ripley (1987, page 46). 32-bit integer assumed */
					   w_rndu = w_rndu*69069+1;
					   return ldexp((double)w_rndu, -32);
				}


				/////////////////////////////////////////////////////////////////////////////////////////////////////
				struct TREEB {
				   int nbranch, nnode, root, branches[NBRANCH][2];
				}  tree;
				struct TREEN {
				   int father, nson, sons[MAXNSONS], ibranch;
				   double branch, age, omega, label, *conP;
				   char *nodeStr, fossil;
				}  *nodes;
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				void NodeToBranchSub (int inode)
				{
				   int i, ison;

				   for(i=0; i<nodes[inode].nson; i++) {
					  tree.branches[tree.nbranch][0] = inode;
					  tree.branches[tree.nbranch][1] = ison = nodes[inode].sons[i];
					  nodes[ison].ibranch = tree.nbranch++;
					  if(nodes[ison].nson>0)  NodeToBranchSub(ison);
				   }
				}

				void NodeToBranch (void)
				{
				   tree.nbranch=0;
				   NodeToBranchSub (tree.root);
				   if(tree.nnode != tree.nbranch+1)
					  error2("nnode != nbranch + 1?");
				}

				void ClearNode (int inode)
				{
					// not sure what ZY's comment below means!!

				/* a source of confusion. Try not to use this routine.
				*/
				   nodes[inode].father=nodes[inode].ibranch=-1;
				   nodes[inode].nson=0;
				   nodes[inode].branch=nodes[inode].age=0;
				   /* nodes[inode].label=0; */
				   /* nodes[inode].branch=0; clear node structure only, not branch lengths */
				   /* FOR (i, com.ns) nodes[inode].sons[i]=-1; */
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				struct CommonInfo {
				   char *z[2*NS-1], spname[NS][LSPNAME+1], daafile[96], cleandata, readpattern;
				   int ns, ls, npatt, np, ntime, ncode, clock, rooted, model, icode;
				   int seqtype, *pose, ncatG, NSsites, fix_omega;
				   double *fpatt, kappa, omega, alpha, pi[64], *conP, daa[20*20];
				   double freqK[NCATG], rK[NCATG];
				   char *siteID;    /* used if ncatG>1 */
				   double *siterates, omega_fix;   /* rates for gamma or omega for site or branch-site models */
				   double *omegaBS, *QfactorBS;     /* omega IDs for branch-site models */
				}  com;
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				int RandomLHistory (int rooted, double space[])
				{
				// Function from evolver.c in PAML
				/* random coalescence tree, with each labeled history having equal probability.
				   interior nodes are numbered ns, ns+1, ..., 2*ns-1-!rooted
				*/
				   int ns=com.ns, i, j, it=0, *nodea=(int*)space;
				   double t;

				   for (i=0; i<2*ns-1-!rooted; i++) ClearNode(i);

				   for (i=0; i<ns; i++) nodea[i]=i;
				   for (i=ns,t=0; i>(1+!rooted); i--) {
					  nodes[it=2*ns-i].nson=2;
					  j=(int)(i*rndu()); 
					  nodes[nodea[j]].father=it; nodes[it].sons[0]=nodea[j];
					  nodea[j]=nodea[i-1];
					  j=(int)((i-1)*rndu()); 
					  nodes[nodea[j]].father=it; nodes[it].sons[1]=nodea[j];
					  nodea[j]=it;
					  if (!rooted && i==3) {
						 nodes[it].nson++; 
						 nodes[nodea[1-j]].father=it; nodes[it].sons[2]=nodea[1-j];
					  }
				   }
				   tree.root=it;  tree.nnode=ns*2-1-!rooted;
				   NodeToBranch();
				   return (0);
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				void BranchLengthBD(int rooted, double birth, double death, double sample, double mut)
				{
					// Function from evolver.c in PAML

				/* Generate random branch lengths (nodes[].branch) using the birth and
				   death process with species sampling, or the Yule (coalescent?) process
				   if sample=0, when only parameter mut is used.
				   Note: older interior nodes have larger node numbers, so root is at
				   node com.ns*2-2 with time t[ns-2], while the youngest node is at 
				   node com.ns with time t[0].  When unrooted=0, the root is removed with
				   branch lengths adjusted.
				   This works with the tree generated from RandomLHistory().
				*/
				   int i,j, it, imin,fixt0=1;
				   double la=birth, mu=death, rho=sample, tmin, r, t[NS-1];
				   double phi, eml, y;

				   if (sample==0)  /* coalescent model.  Check this!!!  */
					  for (i=com.ns,y=0; i>1; i--) 
						  nodes[com.ns*2-i].age=y += -log(rndu())/(i*(i-1.)/2.)*mut/2;
				   else  {         /* BD with sampling */
					  if (fixt0) t[com.ns-2]=1;
					  if (fabs(la-mu)>1e-6) {
						 eml=exp(mu-la);  phi=(rho*la*(eml-1)+(mu-la)*eml)/(eml-1);
						 for (i=0; i<com.ns-1-(fixt0); i++) {
						   r=rndu(); t[i]=log((phi-r*rho*la)/(phi-r*rho*la+r*(la-mu)))/(mu-la);
					   }
					  }
					  else  
						 for (i=0; i<com.ns-1-(fixt0); i++) 
							{ r=rndu();  t[i]=r/(1+la*rho*(1-r)); }
					  /* bubble sort */
					  for (i=0; i<com.ns-1-1; i++) {
						 for (j=i+1,tmin=t[i],imin=i; j<com.ns-1; j++) 
							if (tmin>t[j]) { tmin=t[j]; imin=j; }
						 t[imin]=t[i];  t[i]=tmin;
					  }
					  for (i=com.ns; i>1; i--) nodes[com.ns*2-i].age=t[com.ns-i]*mut;
				   }
				   FOR (i,com.ns) nodes[i].age=0;
				   for (i=0; i<tree.nnode; i++) 
					  if (i!=tree.root) 
						 nodes[i].branch=nodes[nodes[i].father].age-nodes[i].age;
				   if (!rooted) {
					  it=nodes[tree.root].sons[2];
					  nodes[it].branch =
						 2*nodes[2*com.ns-2].age-nodes[tree.root].age-nodes[it].age;
				   }
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				FILE* fout;

				//enum {PrBranch=1, PrNodeNum=2, PrLabel=4, PrAge=8} OutTreeOptions;
				const int PrBranch=1; 
				const int PrNodeNum=2; 
				const int PrLabel=4; 
				const int PrAge=8;
				
				string OutSubTreeN (string &subtree, FILE *fout, int inode, int spnames, int printopt) //, char *labelfmt)
				{
					// Function from evolver.c in PAML - modified by W.F.

					subtree+="(";
					string tempstring;
				   int i,ison;

				   //fputc ('(', fout);
				   for(i=0; i<nodes[inode].nson; i++) {
					  ison=nodes[inode].sons[i];
					  if(nodes[ison].nson==0) {
						 if(spnames) {
							 if(printopt&PrNodeNum) {stringstream fd1; fd1<<ison+1; tempstring=fd1.str(); subtree+=tempstring;}//fprintf(fout, "%d_",ison+1);
							subtree+=com.spname[ison]; //fprintf(tempout,"%s",com.spname[ison]);
						 }
						 else 
						 {stringstream fd2; fd2<<ison+1; tempstring=fd2.str(); subtree+=tempstring;} //fprintf(tempout,"%d",ison+1);
					  }
					  else
						 OutSubTreeN(subtree, fout, ison, spnames, printopt); //, labelfmt);

					  if((printopt & PrNodeNum) && nodes[ison].nson) 
						 {stringstream fd4; fd4<<ison+1; tempstring=fd4.str(); subtree+=tempstring;}//fprintf(tempout," %d ",ison+1);   /* fprintf(fout,"%d",ison+1-sptree.nspecies);  */
					  if((printopt & PrLabel) && nodes[ison].label>0)
						///////////////HERE
					  {stringstream fd5; fd5<<nodes[ison].label; tempstring=fd5.str(); subtree+=tempstring;}  //subtree+=nodes[ison].label; //  {stringstream fd5; fd5<<nodes[ison].label; tempstring=fd5.str(); subtree+=tempstring;}  //fprintf(tempout, labelfmt, nodes[ison].label);
					  if((printopt & PrAge) && nodes[ison].age) 
						{stringstream fd6; fd6<<nodes[ison].age; tempstring=fd6.str(); subtree+=tempstring;} //fprintf(tempout, " @%.3f", nodes[ison].age);


				/*  Add branch labels to be read by Rod Page's TreeView. */


					  if((printopt & PrLabel) && nodes[ison].nodeStr) 
						 subtree+=nodes[ison].nodeStr; //fprintf(tempout," '%s'", nodes[ison].nodeStr);

					  if(printopt & PrBranch)  {stringstream fd7; fd7<<nodes[ison].branch; tempstring=fd7.str(); subtree+=':'; subtree+=tempstring;} //fprintf(tempout,": %.6f", nodes[ison].branch);

					  if(i<nodes[inode].nson-1) subtree+=","; //fprintf(tempout,", ");
				   }
				   subtree+=")"; //fputc (')', fout);
				 //  cout<<endl<<"WWAS"<<endl<<subtree<<endl<<"WWAS"<<endl;

				   return subtree;
				}





				string OutaTreeN (FILE *fout, int spnames, int printopt)
				{
					// Function from evolver.c in PAML - modified by W.F.

				/* 
				   print the current tree.
				*/
					string subtree, tempstring;
				   int i, intlabel=1;
				   //char* labelfmt[2]={"#%.5f", "#%.0f"};

				   if(printopt & PrLabel) {
					  for(i=0;i<tree.nnode;i++) 
						 if(nodes[i].label-(int)nodes[i].label != 0) intlabel=0;
				   }

				   OutSubTreeN(subtree, fout,tree.root,spnames,printopt); //labelfmt[intlabel]);
				   if(printopt&PrNodeNum) {stringstream fd8; fd8<<tree.root+1; tempstring=fd8.str(); subtree+=tempstring;} //fprintf(tempout," %d ", tree.root+1);
				   if((printopt & PrLabel) && nodes[tree.root].label>0) 
					  {stringstream fd5; fd5<<nodes[tree.root].label; tempstring=fd5.str(); subtree+=tempstring;} //subtree+=nodes[tree.root].label;//fprintf(tempout, labelfmt[intlabel], nodes[tree.root].label);
				   if((printopt & PrAge)  && nodes[tree.root].age) 
					  {stringstream fd10; fd10<<nodes[tree.root].age; tempstring=fd10.str(); subtree+=tempstring;}//fprintf(tempout, " @%.3f", nodes[tree.root].age);

				   if((printopt&PrBranch) && nodes[tree.root].branch>1e-8)
					  {stringstream fd11; fd11<<nodes[tree.root].branch; tempstring=fd11.str(); subtree+=":"; subtree+=tempstring;}//fprintf(tempout,": %.6f", nodes[tree.root].branch);

				   subtree+=";"; //fputc(';',tempout);
				//   cout<<"AQA"<<subtree<<"AQA"<<endl;
				   return subtree;
				}


				///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				string newrandomtree(int ntaxa, double birth, double death, double sample, double mut, int randomseed, int option)
				{
					// Function written by W.F. adapting code from the main body of evolver.c in PAML
					
					// produces random trees using the birth-death process

					// Yang and Rannala. (1997) Bayesian phylogenetic inference using DNA sequences: a Markov Chain Monte Carlo Method. Mol. Biol. Evol. 14:717-724.  



					//ofstream tempout("temptree.txt");
					//ns number of species
					// option should be 1 for unrooted and 2 for rooted
					double *space;
					int ntree=1; // number of trees
					int i=randomseed; //randomseed
					com.ns=ntaxa;
					string randomtree;

						 if(com.ns>NS) error2 ("Too many species.  Raise NS.");
						 if((space=(double*)malloc(10000*sizeof(double)))==NULL) error2("oom");

					 int    rooted=!(option%2);
					 int mynewseed= (int)(time(NULL))*(clock());
						 SetSeed(i==-1? mynewseed:i); //(int)time(NULL):i);  //this is what ZY had - makes trees the same if created in quick succession
						int BD=1; //use birth death process for branch lengths

		//				cout<<(i==-1? mynewseed:i)<<"\t";
   
							if(com.ns<3) error2("no need to do this?");
							i=(com.ns*2-1)*sizeof(struct TREEN);
							if((nodes=(struct TREEN*)malloc(i))==NULL) error2("oom");
   

							for(i=0; i<com.ns; i++)          /* default spname */
							   sprintf(com.spname[i],"S%d",i+1);

							for(i=0;i<ntree;i++) {
							   RandomLHistory (rooted, space);
							   BranchLengthBD (1, birth, death, sample, mut);
							   if(com.ns<20&&ntree<10) { randomtree=OutaTreeN(F0,0,BD); }
							   else{randomtree=OutaTreeN(fout,1,BD);}
							   //cout<<" Random Guide Tree set as: "<<endl<<" ";
							  // cout<<randomtree<<endl<<endl;   
							   //FPN(fout);
							}

					free(space);  free(nodes);

					return randomtree;
				}


/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
//            ALL THE CODE BELOW IS FOR GAMMA RATE VARIATION UNTIL STATED OTHERWISE            //
//            IT HAS BEEN USED WITH KIND PERMISSION FROM ZIHENG YANG'S PAML SUITE              //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////
					

				double LnGamma (double alpha);
				double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
				double PointNormal (double prob);
				double PointChi2 (double prob, double v);

					
				double rndgamma1 (double s);
				double rndgamma2 (double s);

				double rndgamma (double s)
				{
					// continuous gamma called with this function

					// 	Yang, Z. (1993) Maximum likelihood estimation of phylogeny from DNA sequences when substitution rates differ over sites. Molecular Biology and Evolution 10:1396-1401.


					double	r=0.0;
					myrand=mtrand1();
					if (s <= 0.0)      
						return 0;
					else if (s < 1.0)  
						r = rndgamma1 (s);
					else if (s > 1.0)  
						r = rndgamma2 (s);
					else           
						r =- log(myrand);
					return (r);
				}


				double rndgamma1 (double s)
				{

					double			r, x=0.0;
					double mysmall=1e-37;
					double w;
					static double	a, p, uf, ss=10.0, d;
					
					if (s!=ss) 
						{
						a  = 1.0-s;
						p  = a/(a+s*exp(-a));
						uf = p*pow(mysmall/a,s);
						d  = a*log(a);
						ss = s;
						}
					for (;;) 
						{
						r = mtrand1();
						if (r > p)        
							x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
						else if (r>uf)  
							x = a*pow(r/p,1/s), w=x;
						else            
							return (0.0);
						r = mtrand1();
						if (1.0-r <= w && r > 0.0)
							if (r*(w+1.0) >= 1.0 || -log(r) <= w)  
								continue;
						break;
						}
					return (x);
				}


				double rndgamma2 (double s)
				{

					double			r ,d, f, g, x;
					static double	b, h, ss=0;
					
					if (s!=ss) 
						{
						b  = s-1.0;
						h  = sqrt(3.0*s-0.75);
						ss = s;
						}
					for (;;) 
						{
						r = mtrand1();
						g = r-r*r;
						f = (r-0.5)*h/sqrt(g);
						x = b+f;
						if (x <= 0.0) 
							continue;
						r = mtrand1();
						d = 64*r*r*g*g*g;
						if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))  
							break;
						}
					return (x);
				}


double CDFNormal (double x)
{
//  Hill ID  (1973)  The normal integral.  Applied Statistics, 22:424-427.
//    Algorithm AS 66.
//    adapted by Z. Yang, March 1994.  Hill's routine is quite bad, and I 
//    haven't consulted 
//      Adams AG  (1969)  Algorithm 39.  Areas under the normal curve.
//      Computer J. 12: 197-198.

    int invers=0;
    double p, limit=10, t=1.28, y=x*x/2;

    if (x<0) {  invers=1;  x*=-1; }
    if (x>limit)  return (invers?0:1);
    if (x<1.28)  
       p = .5 - x * (    .398942280444 - .399903438504 * y
                   /(y + 5.75885480458 - 29.8213557808
		   /(y + 2.62433121679 + 48.6959930692
		   /(y + 5.92885724438))));
    else 
       p = 0.398942280385 * exp(-y) /
           (x - 3.8052e-8 + 1.00000615302 /
           (x + 3.98064794e-4 + 1.98615381364 /
           (x - 0.151679116635 + 5.29330324926 /
           (x + 4.8385912808 - 15.1508972451 /
           (x + 0.742380924027 + 30.789933034 /
           (x + 3.99019417011))))));

    return  invers ? p : 1-p;
}


				double LnGamma (double alpha)
				{
				// returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
				// Stirling's formula is used for the central polynomial part of the procedure.
				// Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
				// Communications of the Association for Computing Machinery, 9:684
				
				   double x=alpha, f=0, z;

				   if (x<7) {
					  f=1;  z=x-1;
					  while (++z<7)  f*=z;
					  x=z;   f=-log(f);
				   }
				   z = 1/(x*x);
				   return  f + (x-0.5)*log(x) - x + .918938533204673 
					  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
						   +.083333333333333)/x;  
				}


double LnGammaFunction (double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;

   if (x<7) {
       f=1;  z=x-1;
       while (++z<7)  f*=z;
       x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*log(x) - x + .918938533204673 
          + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
               +.083333333333333)/x;  
}


				double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
				{
				/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
					   limit of the integration and alpha is the shape parameter.
				   returns (-1) if in error
				   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
				   (1) series expansion     if (alpha>x || x<=1)
				   (2) continued fraction   otherwise
				   RATNEST FORTRAN by
				   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
				   19: 285-287 (AS32)
				*/
				   int i;
				   double p=alpha, g=ln_gamma_alpha;
				   double accurate=1e-8, overflow=1e30;
				   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

				   if (x==0) return (0);
				   if (x<0 || p<=0) return (-1);

				   factor=exp(p*log(x)-x-g);   
				   if (x>1 && x>=p) goto l30;
				   /* (1) series expansion */
				   gin=1;  term=1;  rn=p;
				 l20:
				   rn++;
				   term*=x/rn;   gin+=term;

				   if (term > accurate) goto l20;
				   gin*=factor/p;
				   goto l50;
				 l30:
				   /* (2) continued fraction */
				   a=1-p;   b=a+x+1;  term=0;
				   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
				   gin=pn[2]/pn[3];
				 l32:
				   a++;  b+=2;  term++;   an=a*term;
				   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
				   if (pn[5] == 0) goto l35;
				   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
				   if (dif>accurate) goto l34;
				   if (dif<=accurate*rn) goto l42;
				 l34:
				   gin=rn;
				 l35:
				   for (i=0; i<4; i++) pn[i]=pn[i+2];
				   if (fabs(pn[4]) < overflow) goto l32;
				   for (i=0; i<4; i++) pn[i]/=overflow;
				   goto l32;
				 l42:
				   gin=1-factor*gin;

				 l50:
				   return (gin);
				}


				/* functions concerning the CDF and percentage points of the gamma and
				   Chi2 distribution
				*/
				double PointNormal (double prob)
				{
				/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
				   returns (-9999) if in error
				   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
				   Applied Statistics 22: 96-97 (AS70)

				   Newer methods:
					 Wichura MJ (1988) Algorithm AS 241: the percentage points of the
					   normal distribution.  37: 477-484.
					 Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
					   points of the normal distribution.  26: 118-121.

				*/
				   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
				   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
				   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
				   double y, z=0, p=prob, p1;

				   p1 = (p<0.5 ? p : 1-p);
				   if (p1<1e-20) return (-9999);

				   y = sqrt (log(1/(p1*p1)));   
				   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
				   return (p<0.5 ? -z : z);
				}


				double PointChi2 (double prob, double v)
				{
				/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
				   returns -1 if in error.   0.000002<prob<0.999998
				   RATNEST FORTRAN by
					   Best DJ & Roberts DE (1975) The percentage points of the 
					   Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
				   Converted into C by Ziheng Yang, Oct. 1993.
				*/
				   double e=.5e-6, aa=.6931471805, p=prob, g;
				   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

				   if (p<.000002 || p>.999998 || v<=0) return (-1);

				   g = LnGamma (v/2);
				   xx=v/2;   c=xx-1;
				   if (v >= -1.24*log(p)) goto l1;

				   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
				   if (ch-e<0) return (ch);
				   goto l4;
				l1:
				   if (v>.32) goto l3;
				   ch=0.4;   a=log(1-p);
				l2:
				   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
				   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
				   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
				   if (fabs(q/ch-1)-.01 <= 0) goto l4;
				   else                       goto l2;
  
				l3: 
				   x=PointNormal (p);
				   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
				   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
				l4:
				   q=ch;   p1=.5*ch;
				   if ((t=IncompleteGamma (p1, xx, g))<0) {
					  return (-1);
				   }
				   p2=p-t;
				   t=p2*exp(xx*aa+g+p1-c*log(ch));   
				   b=t/ch;  a=0.5*t-b*c;

				   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
				   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
				   s3=(210+a*(462+a*(707+932*a)))/2520;
				   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
				   s5=(84+264*a+c*(175+606*a))/2520;
				   s6=(120+c*(346+127*c))/5040;
				   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
				   if (fabs(q/ch-1) > e) goto l4;

				   return (ch);
				}


				#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

	//DiscreteGamma(GDfreq,GDrates,(*m).pinv, (*m).alpha,(*m).alpha,(*m).ngamcat,medianORmean);  //makes discrete gamma frequencies etc..   

	//Discrete Gamma + Proportion invariant
			
	//if(mtrand1()>(*m).pinv)  rate=GDrates.at( (int)( mtrand1() * ((*m).ngamcat) ) )/(1-((*m).pinv));  else rate=0;

		
				int DiscreteGamma (vector<double> &cumfreqK, vector<double> &rK, double pinv,
					double alfa, double beta, int K, int median)
				{

				/* 
				
					Yang, Z. (1994) Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. Journal of Molecular Evolution 39:306-314.
	
					discretization of gamma distribution with equal proportions in each category
				*/
					vector<double> freqK;
				   cumfreqK.assign(K,0);
				   freqK.assign(K,0);
				   rK.assign(K,0);
				   int i;
				   double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

				   if (median==1) {
					  for (i=0; i<K; i++) rK.at(i)=PointGamma((i*2.0+1)*gap05, alfa, beta);
					  for (i=0,t=0; i<K; i++) t+=rK.at(i);
					  for (i=0; i<K; i++)     rK.at(i)*=factor/t;
				   }
				   else {
					  lnga1=LnGamma(alfa+1);
					  for (i=0; i<K-1; i++)
					 freqK.at(i)=PointGamma((i+1.0)/K, alfa, beta);
					  for (i=0; i<K-1; i++)
					 freqK.at(i)=IncompleteGamma(freqK.at(i)*beta, alfa+1, lnga1);
					  
					  rK.at(0) = freqK[0]*factor/(1-pinv);
					  rK.at(K-1) = (1-freqK.at(K-2))*factor/(1-pinv);
					  for (i=1; i<K-1; i++)  rK.at(i) = (freqK.at(i)-freqK.at(i-1))*factor/(1-pinv);
				   }
				   for (i=0; i<K; i++) freqK.at(i)=(1.0 - pinv)/K;

				   if(pinv>0) {rK.push_back(0); freqK.push_back(pinv); cumfreqK.push_back(pinv); }
				   

				   cumfreqK.at(0)=freqK.at(0);

				   for(i=1; i<freqK.size(); i++) {cumfreqK.at(i)=freqK.at(i)+cumfreqK.at(i-1); }

				   double diff=cumfreqK.back()-1; if(diff<0) diff=-diff;
				   if(diff>0.00001) cout<<endl<<endl<<"diff is "<<diff<<"   ERROR IN DISCRETEGAMMA CUMFREQS"<<endl<<endl;
//			cout<<"HERE"<<endl;
				   return (0);
				}

					
				double rndbeta(double p, double q)
				{
					double gamma1, gamma2;
					gamma1=rndgamma(p);
					gamma2=rndgamma(q);
					return gamma1/(gamma1+gamma2);
				}

/////////////////////////////////////////////////////////////////////////////////////////////////////

	// code below is to do with codon models discrete rates... used by scipting function to set up M3

/////////////////////////////////////////////////////////////////////////////////////////////////////
				
const int NSnneutral=1;
const int NSpselection=2;
const int NSdiscrete=3;
const int NSfreqs=4;
const int NSgamma=5;
const int NS2gamma=6;
const int NSbeta=7;
const int NSbetaw=8;
const int NSbetagamma=9;
const int NSbeta1gamma=10;
const int NSbeta1normal=11;
const int NS02normal=12;
const int NS3normal=13;


double PDFBeta(double x, double p, double q)
{
/* Returns pdf of beta(p,q)
*/
   double y, mysmall=1e-20;

   if(x<mysmall || x>1-mysmall) 
      error2("bad x in PDFbeta");

   y = (p-1)*log(x) + (q-1)*log(1-x);
   y-= LnGamma(p)+LnGamma(q)-LnGamma(p+q);

   return(exp(y));
}

double CDFBeta(double x, double pin, double qin, double lnbeta)
{
/* Returns distribution function of the standard form of the beta distribution, 
   that is, the incomplete beta ratio I_x(p,q).

   lnbeta is log of the complete beta function; provide it if known,
   and otherwise use 0.

   This is called from InverseCDFBeta() in a root-finding loop.

    This routine is a translation into C of a Fortran subroutine
    by W. Fullerton of Los Alamos Scientific Laboratory.
    Bosten and Battiste (1974).
    Remark on Algorithm 179, CACM 17, p153, (1974).
*/
   double ans, c, finsum, p, ps, p1, q, term, xb, xi, y, mysmall=1e-15;
   int n, i, ib;
   static double eps = 0, alneps = 0, sml = 0, alnsml = 0;

   if(x<mysmall)        return 0;
   else if(x>1-mysmall) return 1;
   if(pin<=0 || qin<=0)  { 
      printf("p=%.4f q=%.4f: parameter outside range in CDFBeta",pin,qin); 
      return (-1); 
   }

   if (eps == 0) {/* initialize machine constants ONCE */
      eps = pow((double)FLT_RADIX, -(double)DBL_MANT_DIG);
      alneps = log(eps);
      sml = DBL_MIN;
      alnsml = log(sml);
   }
   y = x;  p = pin;  q = qin;

    /* swap tails if x is greater than the mean */
   if (p / (p + q) < x) {
      y = 1 - y;
      p = qin;
      q = pin;
   }

   if(lnbeta==0) lnbeta=LnGamma(p)+LnGamma(q)-LnGamma(p+q);

   if ((p + q) * y / (p + 1) < eps) {  /* tail approximation */
      ans = 0;
      xb = p * log(max2(y, sml)) - log(p) - lnbeta;
      if (xb > alnsml && y != 0)
         ans = exp(xb);
      if (y != x || p != pin)
      ans = 1 - ans;
   }
   else {
      /* evaluate the infinite sum first.  term will equal */
      /* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */
      ps = q - floor(q);
      if (ps == 0)
         ps = 1;

      xb=LnGamma(ps)+LnGamma(p)-LnGamma(ps+p);
      xb = p * log(y) - xb - log(p);

      ans = 0;
      if (xb >= alnsml) {
         ans = exp(xb);
         term = ans * p;
         if (ps != 1) {
            n = (int)max2(alneps/log(y), 4.0);
         for(i=1 ; i<= n ; i++) {
            xi = i;
            term = term * (xi - ps) * y / xi;
            ans = ans + term / (p + xi);
         }
      }
   }

   /* evaluate the finite sum. */
   if (q > 1) {
      xb = p * log(y) + q * log(1 - y) - lnbeta - log(q);
      ib = (int) (xb/alnsml);  if(ib<0) ib=0;
      term = exp(xb - ib * alnsml);
      c = 1 / (1 - y);
      p1 = q * c / (p + q - 1);

      finsum = 0;
      n = (int) q;
      if (q == (double)n)
         n = n - 1;
         for(i=1 ; i<=n ; i++) {
            if (p1 <= 1 && term / eps <= finsum)
               break;
            xi = i;
            term = (q - xi + 1) * c * term / (p + q - xi);
            if (term > 1) {
               ib = ib - 1;
               term = term * sml;
            }
            if (ib == 0)
               finsum = finsum + term;
         }
         ans = ans + finsum;
      }
      if (y != x || p != pin)
         ans = 1 - ans;
      if(ans>1) ans=1;
      if(ans<0) ans=0;
   }
   return ans;
}

double InverseCDFBeta(double prob, double p, double q, double lnbeta)
{
/* This calculates the inverseCDF of the beta distribution

   Cran, G. W., K. J. Martin and G. E. Thomas (1977).
   Remark AS R19 and Algorithm AS 109, Applied Statistics, 26(1), 111-114.
   Remark AS R83 (v.39, 309-310) and correction (v.40(1) p.236).

   My own implementation of the algorithm did not bracket the variable well.  
   This version is Adpated from the pbeta and qbeta routines from 
   "R : A Computer Language for Statistical Data Analysis".  It fails for 
   extreme values of p and q as well, although it seems better than my 
   previous version.
   Ziheng Yang, May 2001
*/
   double fpu=3e-308, acu_min=1e-300, lower=fpu, upper=1-2.22e-16;
   /* acu_min>= fpu: Minimal value for accuracy 'acu' which will depend on (a,p); */
   int swap_tail, i_pb, i_inn, niterations=2000;
   double a, adj, g, h, pp, prev=0, qq, r, s, t, tx=0, w, y, yprev;
   double acu, xinbta;

   if(prob<0 || prob>1 || p<0 || q<0) error2("out of range in InverseCDFBeta");

   /* define accuracy and initialize */
   xinbta = prob;

   /* test for admissibility of parameters */
   if(p<0 || q<0 || prob<0 || prob>1)  error2("beta par err");
   if (prob == 0 || prob == 1)
      return prob;

   if(lnbeta==0) lnbeta=LnGamma(p)+LnGamma(q)-LnGamma(p+q);

   /* change tail if necessary;  afterwards   0 < a <= 1/2    */
   if (prob <= 0.5) {
      a = prob;   pp = p; qq = q; swap_tail = 0;
   }
   else {
      a = 1. - prob; pp = q; qq = p; swap_tail = 1;
   }

   /* calculate the initial approximation */
   r = sqrt(-log(a * a));
   y = r - (2.30753+0.27061*r)/(1.+ (0.99229+0.04481*r) * r);

   if (pp > 1. && qq > 1.) {
      r = (y * y - 3.) / 6.;
      s = 1. / (pp*2. - 1.);
      t = 1. / (qq*2. - 1.);
      h = 2. / (s + t);
      w = y * sqrt(h + r) / h - (t - s) * (r + 5./6. - 2./(3.*h));
      xinbta = pp / (pp + qq * exp(w + w));
   }
   else {
      r = qq*2.;
      t = 1. / (9. * qq);
      t = r * pow(1. - t + y * sqrt(t), 3.);
      if (t <= 0.)
         xinbta = 1. - exp((log((1. - a) * qq) + lnbeta) / qq);
      else {
         t = (4.*pp + r - 2.) / t;
         if (t <= 1.)
            xinbta = exp((log(a * pp) + lnbeta) / pp);
         else
            xinbta = 1. - 2./(t+1.);
      }
   }

   /* solve for x by a modified newton-raphson method, using CDFBeta */
   r = 1. - pp;
   t = 1. - qq;
   yprev = 0.;
   adj = 1.;


   
/* Changes made by Ziheng to fix a bug in qbeta()
   qbeta(0.25, 0.143891, 0.05) = 3e-308   wrong (correct value is 0.457227)
*/
   if(xinbta<=lower || xinbta>=upper)  xinbta=(a+.5)/2;

   /* Desired accuracy should depend on (a,p)
    * This is from Remark .. on AS 109, adapted.
    * However, it's not clear if this is "optimal" for IEEE double prec.
    * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));
    * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
    * ---- i.e.,  "new acu" = sqrt(old acu)
    */
   acu = pow(10., -13. - 2.5/(pp * pp) - 0.5/(a * a));
   acu = max2(acu, acu_min);

   for (i_pb=0; i_pb<niterations; i_pb++) {
      y = CDFBeta(xinbta, pp, qq, lnbeta);
      y = (y - a) *
         exp(lnbeta + r * log(xinbta) + t * log(1. - xinbta));
      if (y * yprev <= 0)
         prev = max2(fabs(adj),fpu);
      for (i_inn=0,g=1; i_inn<niterations; i_inn++) {
         adj = g * y;
         if (fabs(adj) < prev) {
            tx = xinbta - adj; /* trial new x */
            if (tx >= 0. && tx <= 1.) {
               if (prev <= acu || fabs(y) <= acu)   goto L_converged;
               if (tx != 0. && tx != 1.)  break;
            }
         }
         g /= 3.;
      }
      if (fabs(tx-xinbta)<fpu) 
         goto L_converged;
      xinbta = tx;
      yprev = y;
   }
//   if(!PAML_RELEASE) 
//    printf("\nInverseCDFBeta(%.2f, %.5f, %.5f) = %.6e\t%d rounds\n", 
//     prob,p,q, (swap_tail ? 1. - xinbta : xinbta), niterations);

   L_converged:
   return (swap_tail ? 1. - xinbta : xinbta);
}






double LineSearch (double(*fun)(double x),double *f,double *x0,double xb[2],double step, double e)
{
/* linear search using quadratic interpolation 

   From Wolfe M. A.  1978.  Numerical methods for unconstrained
   optimization: An introduction.  Van Nostrand Reinhold Company, New York.
   pp. 62-73.
   step is used to find the bracket (a1,a2,a3)

   This is the same routine as LineSearch2(), but I have not got time 
   to test and improve it properly.  Ziheng note, September, 2002
*/
   int ii=0, maxround=100, i;
   double factor=2, step1, percentUse=0;
   double a0,a1,a2,a3,a4=-1,a5,a6, f0,f1,f2,f3,f4=-1,f5,f6;

/* find a bracket (a1,a2,a3) with function values (f1,f2,f3)
   so that a1<a2<a3 and f2<f1 and f2<f3
*/

   if(step<=0) return(*x0);
   a0=a1=a2=a3=f0=f1=f2=f3=-1;
   if(*x0<xb[0]||*x0>xb[1]) 
      error2("err LineSearch: x0 out of range");
   f2=f0=fun(a2=a0=*x0);
   step1=min2(step,(a0-xb[0])/4);
   step1=max2(step1,e);
   for(i=0,a1=a0,f1=f0; ; i++) {
      a1-=(step1*=factor); 
      if(a1>xb[0]) {
         f1=fun(a1);
         if(f1>f2)  break;
         else {
            a3=a2; f3=f2; a2=a1; f2=f1;
         }
      }
      else {
         a1=xb[0];  f1=fun(a1);
         if(f1<=f2) { a2=a1; f2=f1; }
         break;
      }

      /* if(noisy>2) printf("\ta = %.6f\tf = %.6f %5d\n", a2, f2, NFunCall);
      */

   }

   if(i==0) { /* *x0 is the best point during the previous search */
      step1=min2(step,(xb[1]-a0)/4);
      for(i=0,a3=a2,f3=f2; ; i++) {
         a3+=(step1*=factor); 
         if(a3<xb[1]) {
            f3=fun(a3);
            if(f3>f2)  break;
            else 
               { a1=a2; f1=f2; a2=a3; f2=f3; }
         }
         else {
            a3=xb[1];  f3=fun(a3);
            if(f3<f2) { a2=a3; f2=f3; }
            break;
         }

//	 if(noisy>2) printf("\ta = %.6f\tf = %.6f %5d\n", a3, f3, NFunCall);

      }
   }

   /* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
   for (ii=0; ii<maxround; ii++) {
      /* a4 is the minimum from the parabola over (a1,a2,a3)  */

      if (a1>a2+1e-99 || a3<a2-1e-99 || f2>f1+1e-99 || f2>f3+1e-99) /* for linux */
         { printf("\npoints out of order (ii=%d)!",ii+1); break; }

      a4 = (a2-a3)*f1+(a3-a1)*f2+(a1-a2)*f3;
      if (fabs(a4)>1e-100)
         a4=((a2*a2-a3*a3)*f1+(a3*a3-a1*a1)*f2+(a1*a1-a2*a2)*f3)/(2*a4);
      if (a4>a3 || a4<a1)  a4=(a1+a2)/2;  /* out of range */
      else                 percentUse++;
      f4=fun(a4);

      /*
      if (noisy>2) printf("\ta = %.6f\tf = %.6f %5d\n", a4, f4, NFunCall);
      */

      if (fabs(f2-f4)*(1+fabs(f2))<=e && fabs(a2-a4)*(1+fabs(a2))<=e)  break;

      if (a1<=a4 && a4<=a2) {    /* fig 2.2.10 */
         if (fabs(a2-a4)>.2*fabs(a1-a2)) {
            if (f1>=f4 && f4<=f2) { a3=a2; a2=a4;  f3=f2; f2=f4; }
            else { a1=a4; f1=f4; }
         }
         else {
            if (f4>f2) {
               a5=(a2+a3)/2; f5=fun(a5);
               if (f5>f2) { a1=a4; a3=a5;  f1=f4; f3=f5; }
               else       { a1=a2; a2=a5;  f1=f2; f2=f5; }
            }
            else {
               a5=(a1+a4)/2; f5=fun(a5);
               if (f5>=f4 && f4<=f2)
                  { a3=a2; a2=a4; a1=a5;  f3=f2; f2=f4; f1=f5; }
               else {
                  a6=(a1+a5)/2; f6=fun(a6);
                  if (f6>f5)
                       { a1=a6; a2=a5; a3=a4;  f1=f6; f2=f5; f3=f4; }
                  else { a2=a6; a3=a5;  f2=f6; f3=f5; }
               }
            }
         }
      }
      else {                     /* fig 2.2.9 */
         if (fabs(a2-a4)>.2*fabs(a2-a3)) {
            if (f2>=f4 && f4<=f3) { a1=a2; a2=a4;  f1=f2; f2=f4; }
            else                  { a3=a4; f3=f4; }
         }
         else {
            if (f4>f2) {
               a5=(a1+a2)/2; f5=fun(a5);
               if (f5>f2) { a1=a5; a3=a4;  f1=f5; f3=f4; }
               else       { a3=a2; a2=a5;  f3=f2; f2=f5; }
            }
            else {
               a5=(a3+a4)/2; f5=fun(a5);
               if (f2>=f4 && f4<=f5)
                  { a1=a2; a2=a4; a3=a5;  f1=f2; f2=f4; f3=f5; }
               else {
                  a6=(a3+a5)/2; f6=fun(a6);
                  if (f6>f5)
                      { a1=a4; a2=a5; a3=a6;  f1=f4; f2=f5; f3=f6; }
                  else { a1=a5; a2=a6;  f1=f5; f2=f6; }
               }
            }
         }
      }
   }   /*  for (ii) */
   if (f2<=f4)  { *f=f2; a4=a2; }
   else           *f=f4;

   return (*x0=(a4+a2)/2);
}





static double prob_InverseCDF, *par_InverseCDF;
static double (*cdf_InverseCDF)(double x,double par[]);
double diff_InverseCDF(double x);

double diff_InverseCDF(double x)
{
// This is the difference between the given p and the CDF(x), the 
// objective function to be minimized.

   double px=(*cdf_InverseCDF)(x,par_InverseCDF);
   return(square(prob_InverseCDF-px));
}


double InverseCDF(double(*cdf)(double x,double par[]),
       double p,double x,double par[],double xb[2])
{
// Use x for initial value if in range

//   int noisy0=noisy;
   double sdiff,step=min2(0.05,(xb[1]-xb[0])/100), e=1e-15;

//   noisy=0;
   prob_InverseCDF=p;  par_InverseCDF=par; cdf_InverseCDF=cdf;
   if(x<=xb[0]||x>=xb[1]) x=.5;
   LineSearch(diff_InverseCDF, &sdiff, &x, xb, step, e);
//   noisy=noisy0;

   return(x);
}

#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGammaFunction(alpha))


int matout (FILE *fout, double x[], int n, int m)
{
   int i,j;
   for (i=0,FPN(fout); i<n; i++,FPN(fout)) 
      FOR(j,m) fprintf(fout," %11.6f", x[i*m+j]);
   return (0);
}

double CDFdN_dS(double x,double p[])
{
/* This calculates the CDF of the continuous dN/dS distribution over sites, 
   to be used as argument to the routine InverseCDF().  When the distribution
   has spikes, the spikes are ignored in this routine, and the scaling
   is done outside this routine, for example, in DiscreteNSsites().
   All parameters (par) for the w distribution are passed to this routine, 
   although some (p0 for the spike at 0) are not used in this routine.  
   Parameters are arranged in the following order:

      NSgamma (2):       alpha, beta
      NS2gamma (4):      p0, alpha1, beta1, alpha2 (=beta2)
      NSbeta (2):        p_beta, q_beta
      NSbetaw (4):       p0, p_beta, q_beta, w (if !com.fix_omega, not used here)
      NSbetagamma (5):   p0, p_beta, q_beta, alpha, beta
      NSbeta1gamma (5):  p0, p_beta, q_beta, alpha, beta (1+gamma)
      NSbeta1normal (5): p0, p_beta, q_beta, mu, s (normal>1)
      NS02normal (5):    p0, p1, mu2, s1, s2 (s are sigma's)
      NS3normal (6):     p0, p1, mu2, s0, s1, s2 (s are sigma's)

   Parameters p0 & p1 are transformed if (!LASTROUND)

*/
   double cdf=-1;
   double z, f[3],mu[3]={0,1,2},sig[3]; /* 3normal: mu0=0 fixed. mu2 estimated */

   switch(com.NSsites) {
   case(NSgamma):  cdf=CDFGamma(x,p[0],p[1]);   break;
   case(NS2gamma): 
      cdf=p[0] *CDFGamma(x,p[1],p[2])+(1-p[0])*CDFGamma(x,p[3],p[3]);  break;
   case(NSbeta):   cdf=CDFBeta(x,p[0],p[1],0);  break;
   case(NSbetaw):  cdf=CDFBeta(x,p[1],p[2],0);  break;
   case(NSbetagamma):
      cdf=p[0]*CDFBeta(x,p[1],p[2],0)+(1-p[0])*CDFGamma(x,p[3],p[4]);  break;

   case(NSbeta1gamma):
      if(x<=1) cdf=p[0]*CDFBeta(x,p[1],p[2],0);
      else     cdf=p[0]+(1-p[0])*CDFGamma(x-1,p[3],p[4]);
      break;
   case(NSbeta1normal):
      if(x<=1) cdf=p[0]*CDFBeta(x,p[1],p[2],0);
      else {
         cdf=CDFNormal((p[3]-1)/p[4]);
         if(cdf<1e-9) {
            matout(F0,p,1,5);;
            printf("PHI(%.6f)=%.6f\n",(p[3]-1)/p[4],cdf);  getchar();
         }
         cdf=p[0]+(1-p[0])*(1- CDFNormal((p[3]-x)/p[4])/cdf);
      }
      break;
   case(NS02normal):
      mu[2]=p[2]; sig[1]=p[3]; sig[2]=p[4];
      f[1]=p[1];  f[2]=1-f[1];
      cdf = 1 - f[1]* CDFNormal(-(x-mu[1])/sig[1])/CDFNormal(mu[1]/sig[1])
              - f[2]* CDFNormal(-(x-mu[2])/sig[2])/CDFNormal(mu[2]/sig[2]);
      break;
   case(NS3normal):
      mu[2]=p[2]; sig[0]=p[3]; sig[1]=p[4]; sig[2]=p[5];

      if(LASTROUND) { f[0]=p[0]; f[1]=p[1]; }
      else          { z=(f[0]=exp(p[0]))+(f[1]=exp(p[1]))+1; f[0]/=z; f[1]/=z;}
      f[2]=1-f[0]-f[1];
      cdf = 1 - f[0]* 2*CDFNormal(-x/sig[0])
              - f[1]* CDFNormal(-(x-mu[1])/sig[1])/CDFNormal(mu[1]/sig[1])
              - f[2]* CDFNormal(-(x-mu[2])/sig[2])/CDFNormal(mu[2]/sig[2]);
      break;
   }
   return(cdf);
}
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


/*

int DiscreteNSsites(vector<double> &params, int &ncatG, int &NSsites, vector<double> &output)//double a, double b, double c,double d) //double par[])
{
/* This discretizes the continuous distribution for dN/dS ratios among sites
   and calculates freqK[] and rK[], using the median method.
   par[] contains all paras in the w distribution.  par[0] is the 
   proportion of beta if (NSsites==betaw), or the proportion of w=0 if 
   (NSsites=NS02normal).
   This routine uses NSsites, ncatG, freqK, rK.
   betaw has ncatG-1 site classes in the beta distribution, and 02normal 
   has ncatG-1 site classes in the mixed normal distribution.
   See the function CDFdN_dS() for definitions of parameters.
*/
/*
//	double par[4]={a,b,c,d};
int thesize=params.size();
double *par; par=new double[thesize];
double *rK; rK=new double[thesize];
double *freqK; freqK=new double[thesize];

for(int gv=1; gv<params.size(); gv++) par[gv]=params.at(gv);
for( gv=1; gv<params.size(); gv++) cout<<" WEWE "<<par[gv]<<endl;

cout<<ncatG<<" EWEWE "<<NSsites<<endl;
	int status=0, j,off, K=ncatG-(NSsites==NSbetaw || NSsites==NS02normal);

double xb[2]={1e-7,99};  // bounds for omega.  
   int K1=6, K2=4, UseK1K2=0;
   double p01=0, p,w0, lnbeta;

   if(NSsites==NSbeta || NSsites==NSbetaw) xb[1]=1;

#ifdef NSSITES_K1_K2_CLASSES
 //  cout<<"A"<<endl;
   if((NSsites==NSgamma||NSsites==NS2gamma||NSsites>=NSbetagamma)){
      K2=max2(K2,K/3);  K1=K-K2;  UseK1K2=1;
      p01=CDFdN_dS(1.,par);

      // printf("\nK:%3d%3d\t p01=%9.5f\n",K1,K2,p01); 
      FOR(j,K) {
         if(j<K1) { p=p01*(j*2.+1)/(2.*K1); w0=p; }
         else     { p=p01+(1-p01)*((j-K1)*2.+1)/(2.*K2); w0=1.01+(j-K1)/K2; }
         rK[j]=InverseCDF(CDFdN_dS,p,w0,par,xb);
         freqK[j]=(j<K1 ? p01/K1 : (1-p01)/K2); thread
      }
   }
#endif

   if(!UseK1K2) { // this is currently used 
   cout<<"B"<<endl;
      if(NSsites==NSbeta || NSsites==NSbetaw) {
         off=(NSsites==NSbetaw);  // par[0] is proportion for beta for M8 
         lnbeta=LnGamma(par[off])+LnGamma(par[off+1])-LnGamma(par[off]+par[off+1]);
         for(j=0; j<K; j++) {
            p=(j*2.+1)/(2.*K);
            rK[j]=InverseCDFBeta(p, par[off], par[off+1], lnbeta);
         }
      }
      else {
   cout<<"C"<<endl;
         FOR(j,K) {
            p=(j*2.+1)/(2.*K);
            w0=.01+j/K; if(rK[j]) w0=(w0+rK[j])/2;
            rK[j]=InverseCDF(CDFdN_dS,p,w0,par,xb);
			cout<<"WDWD "<<InverseCDF(CDFdN_dS,p,w0,par,xb)<<endl;
         }
      }
      FOR(j,K) freqK[j]=1./K;
   }

   if(NSsites==NSbetaw) {
    //if(!fix_omega) rK[ncatG-1]=par[3];
	   rK[ncatG-1]=par[3];
    //else               rK[ncatG-1]=omega_fix;
      freqK[K]=1-par[0]; FOR(j,K) freqK[j]*=par[0];
   }
   if(NSsites==NS02normal) {
      for(j=K-1;j>=0;j--) // shift to right by 1 to make room for spike at 0
         { rK[j+1]=rK[j]; freqK[j+1]=freqK[j];  }
      rK[0]=0;  freqK[0]=par[0];
      for(j=1;j<K+1;j++) freqK[j]*=(1-par[0]);
   }

   if(NSsites>=NSgamma){
      if(!status && NSsites==NSbeta) 
         for(j=1;j<ncatG;j++) if(rK[j]+1e-7<rK[j-1]) status=1;

      if(status) {
         printf("\nwarning: DiscreteNSsites\nparameters: ");
         FOR(j,(NSsites==7?2:4)) printf(" %12.6f", par[j]);  FPN(F0);
         FOR(j,ncatG)            printf("%13.5f", freqK[j]);  FPN(F0);
         FOR(j,ncatG)            printf("%13.5e", rK[j]);  FPN(F0);
      }
   }

   output.clear();
   cout<<"WOWOWOWOW "<<ncatG<<endl;

   for(int jk=0; jk<ncatG; jk++) 
   {output.push_back(rK[jk]); cout<<rK[jk]<<" WEGFQEGEG"<<endl;}
   
   return(0);
}*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// note from WF - not really used any more as options for models M4-M13 have been disabled.


int DiscreteNSsites(double par[], int ngamcat, int model, vector<double> &output)
{
/* 
Yang, Z., Nielsen, R., Goldman, N. and Pedersen, A-M. K. (2000) Codon-Substitution Models for Heterogeneous Selection Pressure at Amino Acid Sites. Genetics 155: 431-439 (May 2000)

  
  
   This discretizes the continuous distribution for dN/dS ratios among sites
   and calculates freqK[] and rK[], using the median method.
   par[] contains all paras in the w distribution.  par[0] is the 
   proportion of beta if (com.NSsites==betaw), or the proportion of w=0 if 
   (com.NSsites=NS02normal).
   This routine uses com.NSsites, com.ncatG, com.freqK, com.rK.
   betaw has com.ncatG-1 site classes in the beta distribution, and 02normal 
   has com.ncatG-1 site classes in the mixed normal distribution.
   See the function CDFdN_dS() for definitions of parameters.
*/
//	cout<<"!!!!!!!!!! !"<<endl;
	com.ncatG=ngamcat;
	com.NSsites=model;

   int status=0, j,off, K=com.ncatG-(com.NSsites==NSbetaw || com.NSsites==NS02normal);
   //cout<<"KKKKKKKKKKKKKKKKKKK   "<<K<<"  "<<model<<endl;
   double xb[2]={1e-7,99};  /* bounds for omega.  */
   int K1=6, K2=4, UseK1K2=0;
   double p01=0, p,w0, lnbeta;

   if(com.NSsites==NSbeta || com.NSsites==NSbetaw) xb[1]=1;

#ifdef NSSITES_K1_K2_CLASSES
   if((com.NSsites==NSgamma||com.NSsites==NS2gamma||com.NSsites>=NSbetagamma)){
      K2=max2(K2,K/3);  K1=K-K2;  UseK1K2=1;
      p01=CDFdN_dS(1.,par);

      /* printf("\nK:%3d%3d\t p01=%9.5f\n",K1,K2,p01); */
      FOR(j,K) {
         if(j<K1) { p=p01*(j*2.+1)/(2.*K1); w0=p; }
         else     { p=p01+(1-p01)*((j-K1)*2.+1)/(2.*K2); w0=1.01+(j-K1)/K2; }
         com.rK[j]=InverseCDF(CDFdN_dS,p,w0,par,xb);
         com.freqK[j]=(j<K1 ? p01/K1 : (1-p01)/K2); thread
      }
   }
#endif

   if(!UseK1K2) { /* this is currently used */
      if(com.NSsites==NSbeta || com.NSsites==NSbetaw) {
         off=(com.NSsites==NSbetaw);  /* par[0] is proportion for beta for M8 */
         lnbeta=LnGamma(par[off])+LnGamma(par[off+1])-LnGamma(par[off]+par[off+1]);
         for(j=0; j<K; j++) {
            p=(j*2.+1)/(2.*K);
            com.rK[j]=InverseCDFBeta(p, par[off], par[off+1], lnbeta);
         }
      }
      else {
         FOR(j,K) {
            p=(j*2.+1)/(2.*K);
            w0=.01+j/K; if(com.rK[j]) w0=(w0+com.rK[j])/2;
            com.rK[j]=InverseCDF(CDFdN_dS,p,w0,par,xb);
         }
      }
      FOR(j,K) com.freqK[j]=1./K;
   }

   if(com.NSsites==NSbetaw) {
      if(!com.fix_omega) com.rK[com.ncatG-1]=par[3];
      else               com.rK[com.ncatG-1]=com.omega_fix;
      com.freqK[K]=1-par[0]; FOR(j,K) com.freqK[j]*=par[0];
   }
   if(com.NSsites==NS02normal) {
      for(j=K-1;j>=0;j--) /* shift to right by 1 to make room for spike at 0*/
         { com.rK[j+1]=com.rK[j]; com.freqK[j+1]=com.freqK[j];  }
      com.rK[0]=0;  com.freqK[0]=par[0];
      for(j=1;j<K+1;j++) com.freqK[j]*=(1-par[0]);
   }

   if(com.NSsites>=NSgamma){
      if(!status && com.NSsites==NSbeta) 
         for(j=1;j<com.ncatG;j++) if(com.rK[j]+1e-7<com.rK[j-1]) status=1;

      if(status) {
         printf("\nwarning: DiscreteNSsites\nparameters: ");
         FOR(j,(com.NSsites==7?2:4)) printf(" %12.6f", par[j]);  FPN(F0);
         FOR(j,com.ncatG)            printf("%13.5f", com.freqK[j]);  FPN(F0);
         FOR(j,com.ncatG)            printf("%13.5e", com.rK[j]);  FPN(F0);
      }
   }

   for(int i=0; i<ngamcat; i++) output.push_back(com.rK[i]);
   return(0);
}

///////////////////
// old testing of routines.......  from Jan 08 I think ....  delete?

/*
double main(int argc, char* argv[])
{	
vector<double> params;
vector<double> output;
params.push_back(0.2);
params.push_back(0.2);
params.push_back(0.2);
params.push_back(0.2);
params.push_back(0.2);

//int ng=10;
//int blah=5;
//DiscreteNSsites(params,ng , blah, output);

	double parameters[4]=
	//{3,3,3,3}; //
	{0.383,0.967,1.452,0.283};
	
	//int DiscreteNSsites(double par[], int ncatG, int NSsites, double rK[], double freqK[])
DiscreteNSsites(parameters, 10, 6, output);

   for(int i=0; i<10; i++) cout<<"qwdewf "<<output.at(i)<<endl;


	//	parameters);
cout<<"EEEEEEEEEEEEEEEEE"<<endl;
//	cout<<com.rK[0]<<"\t"<<com.rK[1]<<"\t"<<com.rK[2]<<"\t"<<com.rK[3]<<"\t"<<com.rK[4]<<"\t"<<com.rK[5]<<"\t"<<com.rK[6]<<"\t"<<com.rK[7]<<"\t"<<com.rK[8]<<"\t"<<com.rK[9]<<endl;


	return 0;
}

/*
double main(int argc, char* argv[])
{					
	vector<double> freq,rK, freq1,rK1;
	double alfa, beta;
	int K=10, median=0, ig;
	
	alfa=0.967; beta=1.452;
	DiscreteGamma (freq, rK,alfa, beta, K, median);

	for( ig=0; ig<K; ig++) cout<<rK.at(ig)<<endl;

	cout<<"*****************"<<endl ;
	alfa=beta=0.283;
	DiscreteGamma (freq1, rK1,alfa, beta, K, median);

	for( ig=0; ig<K; ig++) cout<<rK1.at(ig)<<endl;

	cout<<"*****************"<<endl;

	for(ig=0; ig<K; ig++) cout<<0.383*rK.at(ig)+0.617*rK1.at(ig)<<endl;

	cout<<"*****************"<<endl ;
	alfa=0.967*0.383+0.617*0.283; beta=1.452*0.383+0.283*0.617;
	DiscreteGamma (freq1, rK1,alfa, beta, K, median);

	for( ig=0; ig<K; ig++) cout<<rK1.at(ig)<<endl;


	
	


	com.NSsites=5;
	double parameters[4]={3,3,3,3}; //0.383,0.967,1.452,0.283};

	com.ncatG=10;
	DiscreteNSsites(parameters);
cout<<"EEEEEEEEEEEEEEEEE"<<endl;
	cout<<com.rK[0]<<"\t"<<com.rK[1]<<"\t"<<com.rK[2]<<"\t"<<com.rK[3]<<"\t"<<com.rK[4]<<"\t"<<com.rK[5]<<"\t"<<com.rK[6]<<"\t"<<com.rK[7]<<"\t"<<com.rK[8]<<"\t"<<com.rK[9]<<endl;
cout<<"EEEEEEEEEEEEEEEEE"<<endl;
	DiscreteGamma (freq1, rK1,3, 3, 10, 0);

	for( ig=0; ig<10; ig++) cout<<rK1.at(ig)<<"\t";

	int a4,b4,c4,d4;
	double a,b,c,d,a1,b1,c1,d1;
	a1=0.383;
	b1=0.967;
	c1=1.452;
	d1=0.283;
/*
	ofstream rout("bits.txt");
	
	a=a1; b=b1; c=c1; d=d1; //a=0.303; //a=0.38;
a=a1;a+=0.0005;
	for(a4=0; a4<101; a4++)
	{b=b1;b+=0.0005;

		a-=0.00001;
	for(b4=0; b4<101; b4++)
	{
		c=c1;c+=0.0005;

		b-=0.00001;
	
	for(c4=0; c4<101; c4++)
	{
d=d1;d+=0.0005;

		c-=0.00001;

	for(d4=0; d4<101; d4++)
	{
		d-=0.00001;
		double diff;
	
		DiscreteNSsites(a,b,c,d); //parameters);
		//for(ig=0; ig<10; ig++) cout<<com.rK[ig]<<endl;
		double T[10]={0.0003,0.0135,0.0598,0.1424,0.2621,0.4267,0.6569,1.0037,1.6282,3.5598};
		double x0=com.rK[0]; diff=T[0]-x0; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x1=com.rK[1]; diff=T[1]-x1; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x2=com.rK[2]; diff=T[2]-x2; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x3=com.rK[3]; diff=T[3]-x3; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x4=com.rK[4]; diff=T[4]-x4; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x5=com.rK[5]; diff=T[5]-x5; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x6=com.rK[6]; diff=T[6]-x6; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x7=com.rK[7]; diff=T[7]-x7; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x8=com.rK[8]; diff=T[8]-x8; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		double x9=com.rK[9]; diff=T[9]-x9; if(diff<0) diff=-diff; if(diff>0.00005) continue;
		
		rout<<a4<<"\t"<<b4<<"\t"<<c4<<"\t"<<d4<<"\t"<<com.rK[0]<<"\t"<<com.rK[1]<<"\t"<<com.rK[2]<<"\t"<<com.rK[3]<<"\t"<<com.rK[4]<<"\t"<<com.rK[5]<<"\t"<<com.rK[6]<<"\t"<<com.rK[7]<<"\t"<<com.rK[8]<<"\t"<<com.rK[9]<<endl;
		cout<<a4<<"\t"<<b4<<"\t"<<c4<<"\t"<<d4<<"\t"<<com.rK[0]<<"\t"<<com.rK[1]<<"\t"<<com.rK[2]<<"\t"<<com.rK[3]<<"\t"<<com.rK[4]<<"\t"<<com.rK[5]<<"\t"<<com.rK[6]<<"\t"<<com.rK[7]<<"\t"<<com.rK[8]<<"\t"<<com.rK[9]<<endl;

	}}}}
*/
/*	return 0;
}
*/
