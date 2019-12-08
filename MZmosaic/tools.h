#if !defined TOOLS_H
#define TOOLS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

float *vector();
float **matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
char **cmatrix();
void free_vector();
void free_dvector();
void free_ivector();
void free_matrix();
void free_dmatrix();
void free_imatrix();
void free_cmatrix();
void nrerror();
int mini();
int maxi();
float minc();
float minf();
float maxf();
double lnfac();
void pswap();
double lognC2();
double lognC4();
void sort();
long setseed();
float ran2();

#endif


