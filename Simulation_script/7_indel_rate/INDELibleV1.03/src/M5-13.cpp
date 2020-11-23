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

#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include "float.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>

using namespace std;




/*  
  //script to calculate discrete M3 values for M4 to M13

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

				int LASTROUND;

void error2 (char * message)
				{ printf("\nError: %s.\n", message); exit(-1); }

#define FOR(i,n) for(i=0; i<n; i++)
#define square(a) ((a)*(a))
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
				#define NS            5000
				#define NBRANCH       (NS*2-2)
				#define MAXNSONS      20
				#define LSPNAME       50
				#define NCODE         64
				#define NCATG         40
					#define FPN(file) fputc('\n', file)
				#define F0 stdout

	
				double LnGamma (double alpha);
				double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);

				double PointNormal (double prob);
				double PointChi2 (double prob, double v);

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
double PDFBeta(double x, double p, double q)
{
/* Returns pdf of beta(p,q)
*/
   double y, small=1e-20;

   if(x<small || x>1-small) 
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
   double ans, c, finsum, p, ps, p1, q, term, xb, xi, y, small=1e-15;
   int n, i, ib;
   static double eps = 0, alneps = 0, sml = 0, alnsml = 0;

   if(x<small)        return 0;
   else if(x>1-small) return 1;
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


int DiscreteNSsites(double par[], int ngamcat, int model, vector<double> &output, vector<double> &freqs)
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
 
   for(int pp=0; pp<ngamcat; pp++) freqs.push_back(com.freqK[pp]);

   return(0);
}





/*  
  //script to calculate discrete M3 values for M4 to M13

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



bool itis(string test, string check)
{
	bool result=true;

	for(int p=0; p<test.size(); p++)
	{
		char c1=test[p];
		bool minitest=false;
		
		for(int y=0; y<check.size(); y++)
		{
			char c2=check[y];
			if(c2==c1) {minitest=true; break;}
		}

		if(minitest==false) {result=false; break;}
	}

	return result;

}



int main(int argc, char* argv[])
{



/*
if(modelnumber==4)
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

				}

*/


	

/*  
  //script to calculate discrete M3 values for M4 to M13

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
	ofstream of1; of1.open("M5-13_output.txt");


while(true)
{

	double mypar[11]={0,0,0,0,0,0,0,0,0,0,0};
	int	nparams[11]={2,4,2,4,5,5,5,5,6,3,4};
	vector<string> modelnames;
	modelnames.push_back("gamma");
	modelnames.push_back("2gamma");
	modelnames.push_back("beta");
	modelnames.push_back("beta&w");
	modelnames.push_back("beta&gamma");
	modelnames.push_back("beta&1+gamma");
	modelnames.push_back("beta&1>normal");
	modelnames.push_back("0&2normal");
	modelnames.push_back("3normal");
	modelnames.push_back("M8a:beta&w=1");
	modelnames.push_back("M8a:beta&w>=1");


	vector<vector<string> > paramlist; vector<string> blah;
		blah.push_back("alpha"); blah.push_back("beta");  paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("alpha1"); blah.push_back("beta1"); blah.push_back("alpha2 (=beta2)"); paramlist.push_back(blah); blah.clear();
		blah.push_back("p_beta"); blah.push_back("q_beta"); paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("p_beta"); blah.push_back("q_beta"); blah.push_back("w");     paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("p_beta"); blah.push_back("q_beta"); blah.push_back("alpha"); blah.push_back("beta");     paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("p_beta"); blah.push_back("q_beta"); blah.push_back("alpha"); blah.push_back("beta"); paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("p_beta"); blah.push_back("q_beta"); blah.push_back("mu"); blah.push_back("s"); paramlist.push_back(blah); blah.clear();
 		blah.push_back("p0"); blah.push_back("p1"); blah.push_back("mu2"); blah.push_back("s1"); blah.push_back("s2"); paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("p1"); blah.push_back("mu2"); blah.push_back("s0"); blah.push_back("s1"); blah.push_back("s2"); paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("p_beta"); blah.push_back("q_beta"); blah.push_back("w=1 fixed"); paramlist.push_back(blah); blah.clear();
		blah.push_back("p0"); blah.push_back("p_beta"); blah.push_back("q_beta"); blah.push_back("w>=1 estimated"); paramlist.push_back(blah); blah.clear();



	int modelnumber;
	string test="!";

	while(true)
	{
		cout<<endl<<" Please enter an integer between 5 and 13 (e.g. 5 for M5, 6 for M6 etc): "<<flush;
		cin>>test;
		if(itis(test,"0123456789"))
		{
			modelnumber=atoi(test.c_str());

			// model number 5 to 13
			if(modelnumber>4 && modelnumber<14) break;
			else cout<<endl<<"    That integer is outside the accepted range."<<endl;
		}
		else cout<<endl<<"    That is not an integer! "<<endl;
	}

	int mysize=nparams[modelnumber-5];

	cout<<endl<<" Model chosen: M"<<modelnumber<<"  ("<<modelnames.at(modelnumber-5)<<")"<<endl<<endl;
	cout<<" This model has "<<mysize<<" parameters."<<endl<<endl;
	cout<<" The parameters are: "<<endl<<endl;
	//blah=paramlist.at(modelnumber-5); cout<<blah.at(0); for(int yg=1; yg<blah.size(); yg++) cout<<", "<<blah.at(yg);
	blah=paramlist.at(modelnumber-5); for(int yg2=0; yg2<blah.size(); yg2++) cout<<"    Parameter "<<yg2+1<<":\t"<<blah.at(yg2)<<endl;
	cout<<endl;

	vector<double> params;

	cout<<" Please enter them in this order pressing enter after each parameter"<<endl<<endl;

	params.assign(mysize,0);

	for(int pi=0; pi<mysize; pi++)
	{ 
		// get input
		double temp;
		string readin;
		cout<<"    Parameter "<<pi+1<<": "<<flush;
		cin>>readin;




		// check string does not contain illegal characters
		string check="0123456789.";

		if( modelnumber==11 && pi==3) check+="-";
	
		if(  (modelnumber==12 || modelnumber==13) && pi==2) check+="-";

		if(!itis(readin,check))
		{
			if(  (modelnumber==12 || modelnumber==13) && pi==2) cout<<endl<<"      ERROR: this parameter must be a real number. Allowed characters are:"<<endl<<"      -.0123456789"<<endl;
				
			else if( modelnumber==11 && pi==3) cout<<endl<<"      ERROR: this parameter must be a real number. Allowed characters are:"<<endl<<"      -.0123456789"<<endl;
	
			else  cout<<endl<<"      ERROR: this parameter must be a decimal number. Allowed characters are:"<<endl<<"      .0123456789"<<endl;

			cout<<"     Try again."<<endl<<endl;  pi--; continue;
		}
		else temp=atof(readin.c_str());
			
		
		// check value is not incorrect
		if(modelnumber==5 || modelnumber==7) { if(temp==0) {cout<<endl<<"      ERROR: this parameter must be greater than zero. Try again."<<endl<<endl; pi--; continue;}}
		
		else
		{
			if(  (modelnumber==12 || modelnumber==13) && pi==2) {}  //mu_2
				
			else if( modelnumber==11 && pi==3) {}  // mu

			else if( modelnumber==8 && pi==3) {}  // omega 
	
			else
			{

				if(pi==0) {if(temp>1) {cout<<endl<<"      ERROR: proportion p0 must be between 0 and 1.  You entered "<<temp<<endl<<"      Try again."<<endl<<endl; pi--; continue;}}

				else if(pi==1) 
				{
					if(modelnumber==12 || modelnumber==13)
					{
						if(temp>1) {cout<<endl<<"      ERROR: proportion p1 must be between 0 and 1.  You entered "<<temp<<endl<<"      Try again."<<endl<<endl; pi--; continue;}
						
						if(temp+params.at(0)>1) {cout<<endl<<"      ERROR: sum of proportions p0 and p1 must be between 0 and 1.\n    You entered "<<temp+params.at(0)<<endl<<"      Try again."<<endl<<endl; pi--; continue;}
					}
					else if(temp==0) {cout<<endl<<"      ERROR: this parameter must be greater than zero. Try again."<<endl<<endl; pi--; continue;}

				}
				else if(temp==0) {cout<<endl<<"      ERROR: this parameter must be greater than zero. Try again."<<endl<<endl; pi--; continue;}
			}			

		}

		
		params.at(pi)=temp;

	}
		
	int ngamcat;  test="";

	cout<<endl<<endl;
	cout<<"  I also need to know how many categories you want to use for"<<endl;
	cout<<"   the discrete approximation to the continuous distribution."<<endl<<endl;

	while(true)
	{
		cout<<endl<<" Please enter an integer between 4 and 40: "<<flush;
		cin>>test;
		if(itis(test,"0123456789"))
		{
			ngamcat=atoi(test.c_str());

			if(ngamcat>3 && ngamcat<41) break;
			else cout<<endl<<"    That integer is outside the accepted range."<<endl;
		}
		else cout<<endl<<"    That is not an integer! "<<endl;
	}
	
	
	
	cout<<endl<<endl<<"******************************************************"<<endl<<endl;
	


	vector<double> output2, cumfreqs, freqs2;
	for(int iu=0; iu<params.size(); iu++) {mypar[iu]=params.at(iu); }  
	
	DiscreteNSsites(mypar, ngamcat, modelnumber, output2, freqs2);

	vector<string> output,freqs; for(int tg=0; tg<ngamcat; tg++)
	{
		stringstream sd1; sd1<<output2.at(tg); output.push_back(sd1.str());
		stringstream sd2; sd2<<freqs2.at(tg); freqs.push_back(sd2.str());
	}

	/*
					if(modelnumber==12 || modelnumber==8)
					{
						if(modelnumber==12)
						{
							vector<double> temp=output; output.clear(); output.push_back(0); for(int y=0; y<temp.size(); y++) output.push_back(temp.at(y));
					 					
							cumfreqs.push_back(params.at(0));
							for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back((1-params.at(0))/double(ngamcat));
						}
						else
						{
							vector<double> temp=output; output.clear(); output.push_back(params.at(3)); for(int y=0; y<temp.size(); y++) output.push_back(temp.at(y));
						 					
							cumfreqs.push_back(1-params.at(0));
							for(int hfd=0; hfd<ngamcat; hfd++) cumfreqs.push_back(params.at(0)/double(ngamcat));
						}
					}
					else cumfreqs.assign(ngamcat,1/double(ngamcat));
*/
	cout<<"Your [submodel] command for INDELible has been output in the file:"<<endl<<"M5-13_output.txt"<<endl<<endl;

	of1<<"/*"<<endl;
	of1<<"  M3 approximation to M"<<modelnumber<<"  ("<<modelnames.at(modelnumber-5)<<")  with "<<ngamcat<<" categories used in discrete approximation"<<endl<<endl;
	of1<<" You chose the following values:"<<endl<<endl;
	blah=paramlist.at(modelnumber-5); for(int yg1=0; yg1<blah.size(); yg1++) of1<<"    Parameter "<<yg1+1<<":\t"<<blah.at(yg1)<<"\t"<<params.at(yg1)<<endl;
	of1<<"*/"<<endl<<endl;

	of1<<"  [submodel]\tkappa  // use your own value for kappa here!"<<endl;

	of1<<"\t\t"; for(int s=0; s<freqs.size(); s++) 
	{
		string s1=output.at(s), s2=freqs.at(s); 
		if(s==freqs.size()-1) s2="";
		int diff=s1.size()-s2.size();  

		of1<<s2<<" "; 

		if(diff>0) for(int fd=0; fd<diff; fd++) of1<<" ";
	}
	of1<<"  // proportions"<<endl;

	of1<<"\t\t"; for(int t=0; t<output.size(); t++) 
	{
		string s1=output.at(t), s2=freqs.at(t); 
		
		int diff=s2.size()-s1.size();  

		of1<<s1<<" "; 

		if(diff>0) for(int fd=0; fd<diff; fd++) of1<<" ";
	}
	of1<<"  // omega values"<<endl;
	

	
	

//	cout<<"              "kappa"<<endl;
	of1<<endl<<endl<<"/*********************************************/"<<endl<<endl;



	
	}


  //script to calculate discrete M3 values for M4 to M13
/*
	 5  "gamma"         2:    alpha, beta
	 6  "2gamma"        4:    p0, alpha1,beta1, alpha2=beta2
	 7  "beta"          2:    p_beta, q_beta
	 8  "beta&w"        4:    p0, p_beta, q_beta, w estimated
	 9  "beta&gamma"    5:    p0, p_beta, q_beta, alpha, beta
	10  "beta&1+gamma"  5:    p0, p_beta, q_beta, alpha, beta (1+gamma used)
	11  "beta&1>normal" 5:    p0, p_beta, q_beta, mu, s    (normal truncated w>1)
	12  "0&2normal"     5:    p0, p1, mu2, s1, s2
	13  "3normal"       6:    p0, p1, mu2, s0, s1, s2
	14  "M8a:beta&w=1"  3:    p0, p_beta, q_beta, w=1 fixed
	15  "M8a:beta&w>=1" 4:    p0, p_beta, q_beta, w>=1 estimated
*/



				

	return 0;
}