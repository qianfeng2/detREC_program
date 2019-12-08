#include "tools.h"

void nrerror(error_text)
char error_text[];
{
	void exit();

	fprintf(stderr,"\n\nRun-time error\n\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"\n..now exiting to system...\n");
	exit(1);
}



float *vector(nl,nh)
int nl,nh;
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}



float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	float **m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

char **cmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
        int i;
	char **m;

        m=(char **)malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
        if (!m) nrerror("allocation failure 1 in cmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(char *)malloc((unsigned) (nch-ncl+1)*sizeof(char));
                if (!m[i]) nrerror("allocation failure 2 in cmatrix()");
                m[i] -= ncl;
        }
        return m; 
}


void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
	free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
	free((char*) (v+nl));
}



void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) if ((m[i]+ncl) != NULL) free((char*) (m[i]+ncl));
	if ((m+nrl) != NULL) free((char*) (m+nrl));
}

        
void free_cmatrix(m,nrl,nrh,ncl,nch)
char **m;
int nrl,nrh,ncl,nch;
{
        int i;
 
        for(i=nrh;i>=nrl;i--) if ((m[i]+ncl) != NULL) free((char*) (m[i]+ncl));
        if ((m+nrl) != NULL) free((char*) (m+nrl));
}



int mini(i,j) 
int i,j;
{
        if (i<j) return i;
        else return j; 
}

int maxi(i,j) 
int i,j;
{
	if (i>j) return i;
	else return j;
}


float minf(f1, f2)
     float f1, f2;
{
  if (f1<f2) return (float) f1;
  else return (float) f2;
}

float maxf(f1, f2)
     float f1, f2;
{
  if (f1>f2) return (float) f1;
  else return (float) f2;
}

double lnfac(i) 
int i;
{

  int j;
  double cp=0.0;

  for (j=2;j<=i;j++) cp += log((double) j);
  return cp;
}


float minc(l1,l2,ls) 
float l1,l2,ls;
{
        
        double d;
        if ((d= (l2-l1)) >  ls/2) d = ls-l2+l1;
        return (float) d; 
}


void pswap(pt,s1,s2) 
int *pt,s1,s2;
{
                 
        int tmp;
                 
        tmp = pt[s2];
        pt[s2]=pt[s1];
        pt[s1]=tmp; 
}


double lognC2(n,a)
int n, a;
{
	int i;
	double x=0;

	if (a>n) nrerror("Error in lognC2 (1)");
	if ((n<2)||(a==0)||(a==n)) return 0.0;

	for (i=a+1;i<=n;i++) x += (double) log(i);
	for (i=2;i<=n-a;i++) x -= (double) log(i);

	if (x>0) return x;
	else {
		nrerror("Error in logCn2 (2)");
		return 0;
	}
}

double lognC4(n,a,b,c,d)
int n,a,b,c,d;
{
	int i;
	double x=0;

	if (a+b+c+d != n) nrerror("Error in logCn4 (1)");
	if ((n<2)||(a==n)||(b==n)||(c==n)||(d==n)) return 0.0;

	for (i=a+1;i<=n;i++) x += (double) log(i);
	for (i=2;i<=b;i++) x -= (double) log(i);
	for (i=2;i<=c;i++) x -= (double) log(i);	
	for (i=2;i<=d;i++) x -= (double) log(i);

	if (x>0) return x;
	else {
		nrerror("Error in lognC4 (2)");
		return 0;
	}
}


void sort(array, ne)
float *array;
int ne;
{
	int pass, i;
	float tmp;

	for (pass=1;pass<=ne;pass++)
		for (i=1;i<ne;i++) 
			if (array[i+1]<array[i]) {
				tmp = array[i];
				array[i]=array[i+1];
				array[i+1]=tmp;
			}
}


/*Routine to set seed from clock*/

long setseed(void) {
	time_t lt;
	lt=time(NULL);
	return lt;
}


/*C routine random number generator from
Numerical Recipes in C: Press et al*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(void) {
	int j;
	long k;
	extern long *idum;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7; j>=0; j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
			
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j]=*idum;
	if (iy < 1) iy += IMM1;
/*	if (((float) temp = (float) AM*iy) > RNMX) return (float)  RNMX;*/
	if ((temp = AM*iy) > RNMX) return RNMX;
	else return (float) temp;
}


