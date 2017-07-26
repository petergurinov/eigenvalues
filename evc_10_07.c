#include "task_10_07.h"

extern char d, dd, e;

//tmp memory
//n*n - Q
//n - X
//n*n - R
//n*n - U

int cmp(const void *x, const void *y)
{
	double xx = *(double*)x, yy = *(double*)y;
	if (xx < yy) return -1;
	if (xx > yy) return  1;
	return 0;
}

static double norma_10_07(int n, double *a)
{
	int i;
	double norma=0;

	if (dd) {
		printf("norma: sqrt( ");
	}
	for(i=0; i < n; i++) {
		if (dd) {
			printf("%1.2lg^2 + ", a[i]);
		}
		norma += a[i]*a[i];
	}
	norma = sqrt(norma);
	if (dd) {
		printf(") = %1.2lg\n", norma);
	}
	return norma;
}

double norma_matr_10_07(int n, double *A)
{
	int i, j;
	double norma, max=0;
	for(i=0; i < n; i++) {
		norma = 0;
		for(j=0; j < n; j++) {
			norma += fabs(A[i*n+j]);
		}
		if (max < norma) {
			max = norma;
		}
	}
	return max;
}

static double u_10_07(int i, int j, double *x)//elements of matr U[i,j]
{
	if (i == j)
		return 1-2*x[i]*x[j];
	return 0-2*x[i]*x[j];
}


void copy_matr_10_07(int n, double *A, double *B) //copy matr from B to A
{
	int i;
	for(i=0; i < n*n; i++)
		A[i] = B[i];
}

void transpose_matr_10_07(int n, double *A) //get matr A and make one Transpose
{
	int i, j;
	double t;
	for(i=0; i < n; i++) {
		for(j=0; j < i; j++) {
			t = A[i*n+j];
			A[i*n+j] = A[j*n+i];
			A[j*n+i] = t;
		}
	}
}

void exhaustion_10_07(int n, double *A, double epsilon)
{
	int i;
	double norma;
	norma = norma_matr_10_07(n, A);
	for(i=0; i < n - 1; i++) {
		if (fabs(A[(i + 1)*n+i]) < epsilon*norma)
			A[(i + 1)*n+i] = 0;
	}
}

double shift_10_07(int n, double *A, int *l,double epsilon)
{
	int i;
	double s_k;

	if ((fabs(A[(*l)*n+(*l)-1]) > epsilon ) && ((*l) >= 0)) {
		s_k = A[(*l)*n+(*l)];
		(*l)--;
	} else {
		s_k = 0;
	}

	
	for(i=0; i < n; i++)
		A[i*n+i] -= s_k;

	return s_k;
}

void acceleration(int n, double *Q, double *R, double epsilon)
{
	int i, k, check;
	for(i=1; i < n; i++) {
		if (fabs(fabs(Q[i*n+i]) - 1) < epsilon) {
			check = 1;
			for(k=0; k < i; k++) {
				if ( fabs(Q[i*n+k]) > epsilon || fabs(Q[k*n+i]) > epsilon )
					check *= 0;
			}
			if (check) {
				for(k=0; k < n; k++) {
					if (k != i) {
						R[i*n+k] = 0;
						R[k*n+i] = 0;
					}
				}
			}
		}
	}
	return;
}

int check_q_10_07(int n, double *Q, double eps) //1 - true; 0 - false
{
	int i, j;

	for(i=0; i < n; i++) {
		for(j=0; j < n; j++) {
			if (fabs(fabs(Q[i*n+j]) - (i == j?1:0)) > eps)
				return 0;
		}
	}
	return 1;
}

int chek_count_step_10_07(int y, int max_iterations)
{
	if ((max_iterations == 0) && (y > 4000)) return 0;
	if (max_iterations == 0) return 1;
	if (max_iterations > y) return 1;
	return 0;
}

int evc_10_07(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision)
{
	int i, j, k, y=1, l;
	double *Q = tmp;
	double *X = Q + n*n;
	double *R = X + n;
	double *U = R + n*n;
	double norma;
	double s_k = 0;
	l = n-1;
	if (d) {
		printf("Start QR algorithm for finding the eigenvalues!\n");
	}
	do { //find Q and R
		for(i=0; i < n; i++) { //Q matr
			for(j=0; j < n; j++) {
				if (i == j) { Q[i*n+j] = 1; } else { Q[i*n+j] = 0; }
			}
		}
		for(i=0; i < n; i++) { //R matr
			for(j=0; j < n; j++) {
				R[i*n+j] = A[i*n+j];
			}
		}
		exhaustion_10_07(n, A, epsilon);
		s_k = shift_10_07(n, A, &l, epsilon);
		if (dd) {
			printf("l = %d; s_k = %1.3lg; A after shift -", l, s_k);
			print_matr(n, A);
			printf("==\n");
		}
		for(k=0; k < n; k++) { //by column
			for(i=0; i < k; i++)
				X[i] = 0;
			for(i=k; i < n; i++)
				X[i] = A[i*n+k];
			norma = norma_10_07(n, X);//norma of a_k
			X[k] -= norma;
			if (fabs(X[k]) < precision) {
				X[k] = 0;
			}
			if (dd) {
				printf("^x_k ");
				print_vec(n, X);
			}
			norma = norma_10_07(n, X);
			if (fabs(norma) < precision) {
				if (dd) {
					printf("norma x == %1.3lg\n", norma);
				}
				continue;
			}
			for(i=0; i < n; i++) {
				X[i] /= norma;
				if (fabs(X[i]) < precision) {
					X[i] = 0;
				}
			}
			if (dd) {
				printf("~x_k ");
				print_vec(n, X);
				printf("U = ");
			}
			for(i=0; i < n; i++) {
				for(j=0; j < n; j++) {
					U[i*n+j] = u_10_07(i, j, X);
				}
			}
			precision_10_07(n, U, precision);
			if (dd) {
				print_matr(n, U);
				printf("Q = ");
				print_matr(n, Q);
				printf("Q <- R <- Q*U = ");
			}
			multip_matr_10_07(n, Q, U, R);
			precision_10_07(n, R, precision);
			copy_matr_10_07(n, Q, R);
			if (dd) {
				print_matr(n, Q);
				printf("A = ");
				print_matr(n, A);
				printf("A <- R <- U*A = ");
			}
			multip_matr_10_07(n, U, A, R);
			precision_10_07(n, R, precision);
			copy_matr_10_07(n, A, R);
			if (dd) {
				print_matr(n, A);
				printf("====\n");
			}
		}
/*		transpose_matr_10_07(n, Q);*/
		if (d) {
			printf("step: %d: \n", y);
			printf("Q = ");
			print_matr(n, Q);
			printf("R = ");
			print_matr(n, R);
		}
		acceleration(n, Q, R, epsilon);
		if (dd) {
			printf("after accel R = ");
			print_matr(n, R);
		}
		multip_matr_10_07(n, Q, R, A);
		precision_10_07(n, A, precision);
		if (dd) {
			printf("A = Q*R ");
			print_matr(n, A);
			printf("====\n");
		}
		multip_matr_10_07(n, R, Q, A);
		precision_10_07(n, A, precision);
		for(i=0; i < n; i++) { //for shift_10_07
			A[i*n+i] += s_k;
		}
		precision_10_07(n, A, precision);
		if (d) {
			printf("A = R*Q ");
			print_matr(n, A);
			if (dd) {
				printf("l = %d; s_k = %1.3lg; A after shift +", l, s_k);
				print_matr(n, A);
			}
			printf("===========\n");
		}
		if (l < 0) break;
		y++;
	} while (!check_q_10_07(n, Q, epsilon) && chek_count_step_10_07(y-1, max_iterations) );

	for(i=0; i < n; i++)
		E[i] = A[i*n+i];

	qsort(E, n, sizeof(double), cmp);

	if (d) {
		printf("E vector:\n");
		print_vec(n, E);
		printf("===========\n");
	}


	return 0;
}
