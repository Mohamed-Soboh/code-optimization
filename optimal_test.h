#pragma once


typedef double (*FUN_PTR)(double[]); // pointer to function 
typedef void (*GRAD_FUN_PTR)(double grad[], double x[]); // pointer 
typedef int VECTOR_CONVERGENCE_TEST(double arr[], int n, double epsilon);


FUN_PTR objective_function;

FUN_PTR global_h[NMAX];

double grad_vector[NMAX];
double xnm1[NMAX];
double gtemp[NMAX];
double lambda[NMAX];
double grad_vector[NMAX];
double diff[NMAX];
double delta_h[NMAX][NMAX];
double FLH[NMAX][NMAX];
double Reduced_F[NMAX][NMAX];

int vector_n;
int global_n;
int global_m;
double epsilon = 0.001;


int naux;
double Aaux[NMAX][NMAX];
int basis[NMAX];
int basis_size;
double y[NMAX];
double E[NMAX][NMAX];
double Et[NMAX][NMAX];
double temp[NMAX][NMAX];
double B[NMAX][NMAX];
double Eaux[NMAX];
double W[NMAX][NMAX];
double precision_epsilon = 0.000001;