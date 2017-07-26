#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifdef TASK_10_07
#define TASK_10_07

int usage(FILE* fl);
int evc_memsize_10_07(int n);
int sim_memsize_10_07(int n);
void print_matr(int n, double *A);
void print_vec(int n, double *B);


double norma_10_07(int n, double *a);
double u_10_07(int i, int j, double *x);
void multip_matr_10_07(int n, double *A, double *B, double *C);
void precision_10_07(int n, double *A, double precision);
double norma_matr_10_07(int n, double *A);
void copy_matr_10_07(int n, double *A, double *B);
void transpose_matr_10_07(int n, double *A);
void exhaustion_10_07(int n, double *A, double epsilon);
double shift_10_07(int n, double *A, int *m);
int check_q_10_07(int n, double *Q, double eps);
int chek_count_step_10_07(int y, int max_iterations);
void acceleration(int n, double *Q, double *R, double epsilon);


int sim_10_07(int n, double* A, double* tmp, double precision);
int evc_10_07(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision);

#endif TASK_10_07
