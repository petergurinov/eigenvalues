#include "task_10_07.h"

extern char d, dd, e;

//tmp memory
//n*n - T
//n - X
//n*n - U

static double norma_10_07(int n, double *a)
{
	int i;
	double norma=0;

	if (d)
		printf("norma: sqrt( ");
	for(i=0; i < n; i++) {
		if (d)
			printf("%1.2lg^2 + ", a[i]);
		norma += a[i]*a[i];
	}
	norma = sqrt(norma);
	if (d)
		printf(") = %1.2lg\n", norma);
	return norma;
}

static double u_10_07(int i, int j, double *x)//elements of matr U[i,j]
{
	if (i == j)
		return 1-2*x[i]*x[j];
	return 0-2*x[i]*x[j];
}

void multip_matr_10_07(int n, double *A, double *B, double *C)//multiplay matr A*B and write answer to C
{
	int i, j, k;
	for(i=0; i < n; i++) {
		for(j=0; j < n; j++) {
			C[i*n+j] = 0;
			for(k=0; k < n; k++) {
				C[i*n+j] += (A[i*n+k] * B[k*n+j]);
			}
		}
	}
	return;
}

void precision_10_07(int n, double *A, double precision)
{
	int i;
	for(i=0; i < n*n; i++)
		if (fabs(A[i]) < precision)
			A[i] = 0;
}

int sim_10_07(int n, double* A, double* tmp, double precision)
{
	int i, j, k;
	double *T = tmp;
	double *X = T + n*n;
	double *U = X + n;
	double norma;
	if (d)
		printf("Start algorithm to almost triangular form!\n");
	for(i=0; i < n-1; i++) { //by column
		for(k=0; k < i+1; k++)
			X[k] = 0;
		for(k=i+1; k < n; k++)
			X[k] = A[k*n+i];
		norma = norma_10_07(n, X);//norma of a_i
		X[i+1] -= norma;
		if (dd) {
			printf("^x_i ");
			print_vec(n, X);
		}
		norma = norma_10_07(n, X);
		if (norma == 0)
			continue;
		for(k=0; k < n; k++)
			X[k] /= norma;
		if (dd) {
			printf("~x_i ");
			print_vec(n, X);
		}
		if (dd) {
			printf("U ");
		}
		for(k=0; k < n; k++) {
			for(j=0; j < n; j++) {
				U[k*n+j] = u_10_07(k, j, X);
			}
		}
		if (dd) {
			print_matr(n, U);
			printf("T = U*A ");
		}
		multip_matr_10_07(n, U, A, T);
		precision_10_07(n, T, precision);
/*		transpose_matr_10_07(n, U);*/
		if (dd) {
			print_matr(n, T);
			printf("A = T*U ");
		}
		multip_matr_10_07(n, T, U, A);
		precision_10_07(n, A, precision);
		if (dd) {
			print_matr(n, A);
			printf("======\n");
		}
	}
	if (d) {
		printf("A ");
		print_matr(n, A);
		printf("========\n");
	}

	return 0;
}
