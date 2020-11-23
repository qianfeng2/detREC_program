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



// uniform random numbers are produced using the code in this file.
#include "MersenneTwister.h"

MTRand mtrand1;
MTRand mtrand2;
	
// needed for time monitoring

#include <stdio.h>
#include <stdlib.h>
/*
// time monitoring stuff works as below:

clock_t start, finish;
   double  duration;

 // Measure the duration of an event. 
   printf( "Time to do %ld empty loops is ", i );
   start = clock();
   while( i-- )
      ;
   finish = clock();
   duration = (double)(finish - start) / CLOCKS_PER_SEC;
   printf( "%2.1f seconds\n", duration );
*/

#include <time.h>        
#include <math.h>
//////////////////////////////////////////////////////////////

double myrand;
int idum;


// exponential random number generation
double expmyrand(double myrate)
{
	double dum;
	do
	dum=mtrand1();
	while (dum ==0.0);
	return -(log(dum)/myrate);
}


// geometric random number generation
int oldgeomyrand(double p) 
{
	if(p==1) return 1;
	else
	{
		double dum;
		do
		dum=mtrand1();
		while (dum ==0.0);
		return int(log(dum)/log(1-p));
	}
}

// alternative parameterisation - for testing
int geomyrand(double q)  
{
	if(q==0) return 1;
	else
	{
		double dum;
		do
		dum=mtrand1();
		while (dum ==0.0);
		return int(log(dum)/log(q));
	}
}

// negative binomial distribution
int randnegbin(int r, double q)
{
	int u=0;
	while(r--)
		u += geomyrand(q);
	return u+1;
}

// alternative parameterisation - for testing
int oldrandnegbin(int r, double q)
{
	double p=1-q;
	int u=0;
	while(r--)
		u += oldgeomyrand(p);
	return u+1;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Zipfian random number generation

#define H(x,q1,q2,v) exp(q1*log(v+x))*q2

#define H1(x,q1,q2,v) -v+exp(q2*log((1-q)*x))

unsigned int imax = ~0;
	
	
double Zq1,  Zq2,  ZHx0,  Zs,  ZHimax;

int Zipf(double q, double v)
{
	// rejection-inversion method of Hörmann and Derflinger (1996) 
	
	// Hörmann, W. and G. Derflinger, G. (1996) Rejection-inversion to generate variates from monotone discrete distributions. ACM Transactions on Modelling and Computer Simulation, 6(3), 169-184 

	// produces random zipf number from distribution with following unnormalised probability function:

    // p_k = (v + k)^(-q)  where k=0,1,... and q>1, v>0

	// at least twice as fast as method of Devroye Luc (1986) used in DAWG - see below
	
	double U, X, K;

	do
	{
		U = mtrand1();

		U = ZHimax + ( U * ( ZHx0 - ZHimax ) );

		X = H1(U,Zq1,Zq2,v);

		K = X + 0.5;

		if( K - X <= Zs) return int(K+1);

		else if ( U >= H(K+0.5,Zq1,Zq2,v) - exp(-1 * log(v + K) * q) ) return int(K+1);
		
	} while(true); 
}
#define newZipf(a) Zipf(a,1)
  


// algorithm from DAWG (Cartwright, 2005).
// Draw from Zipf distribution, with parameter a > 1.0
// Devroye Luc (1986) Non-uniform random variate generation.
// Springer-Verlag: Berlin. p551

int oldZipf(double a)
{
	double b = pow(2.0, a-1.0);
	double x,t;
	do {
	 x = floor(pow(mtrand1(), -1.0/(a-1.0)));
	 t = pow(1.0+1.0/x, a-1.0);
	} while( mtrand1()*x*(t-1.0)*b > t*(b-1.0));
	return (int)x;
}

