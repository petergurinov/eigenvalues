#include "task_10_07.h"

#define DEFAULT_INPUT_FILENAME "10_07_in.txt"
#define DEFAULT_OUTPUT_FILENAME "10_07_out.txt"
#define DEFAULT_D_OPTION 0	/*0-OFF; 1-ON*/
#define DEFAULT_DD_OPTION 0
#define DEFAULT_E_OPTION 0
#define DEFAULT_P_OPTION 0
#define DEFAULT_T_OPTION 0
#define DEFAULT_PREC_OPTION 1e-14
#define DEFAULT_EPSILON_OPTION 1e-10
#define DEFAULT_MAXITER_OPTION 0

char d = DEFAULT_D_OPTION, dd = DEFAULT_DD_OPTION, e = DEFAULT_E_OPTION;

int usage(FILE* fl)
{
	fprintf(fl, "Usage: lss [input_filename] [output_filename] [options]\n");
	fprintf(fl, "Where options include:\n");
	fprintf(fl, "	-d		print debug messages [default OFF]\n");
	fprintf(fl, "	-dd		print deep debug messages [default OFF]\n");
	fprintf(fl, "	-e		print errors [default OFF]\n");
	fprintf(fl, "	-p		print matrix [default OFF]\n");
	fprintf(fl, "	-t		print execution time [default OFF]\n");
	fprintf(fl, "	-h, -?		print this help and exit\n");
	fprintf(fl, "		-prec=<num>	precision [default - 1e-14]\n");
	fprintf(fl, "		-eps=<num>	'epsilon' [default - 1e-10]\n");
	fprintf(fl, "		-max_iter=<num>	limit number of iterations [default - 0, i.e. not limit]\n");
	fprintf(fl, "Default input_filename value is %s\n", DEFAULT_INPUT_FILENAME);
	fprintf(fl, "Default output_filename value is %s\n\n",  DEFAULT_OUTPUT_FILENAME);
	return 0;
}

int evc_memsize_10_07(int n)
{
	return (sizeof(double))*(n*n + n + n*n + n*n);
}

int sim_memsize_10_07(int n)
{
	return (sizeof(double))*(n*n + n + n*n);
}

void print_matr(int n, double *A)
{
	int i, j;
/*	printf("N: %d\n", n);*/
	printf("Matr:\n");
	printf("{ ");
	for (i = 0; i < n; i++) {
		if (i != 0) { printf("  { "); } else { printf("{ "); }
		for (j = 0; j < n; j++) {
			printf("%1.2lg",A[i*n+j]);
			if (j < n-1) { printf(", "); }
		}
		if (i < n-1) { printf("}, \n"); } else { printf("}"); }
	}
	printf(" }\n");
	return;
}

void print_vec(int n, double *B)
{
	int i;
/*	printf("N: %d\n", n);*/
	printf("Vec:\n");
	for (i = 0; i < n; i++)
		printf("%1.2lg ",B[i]);
	printf("\n");
	return;
}

int main(int argc, char **argv)
{
	FILE *fin;
	FILE *fout;
	char *input_filename = DEFAULT_INPUT_FILENAME;
	char *output_filename = DEFAULT_OUTPUT_FILENAME;
	char p = DEFAULT_P_OPTION, t = DEFAULT_T_OPTION;
	double prec = DEFAULT_PREC_OPTION, eps = DEFAULT_EPSILON_OPTION;
	int max_iter = DEFAULT_MAXITER_OPTION;
	
	int i, j, answ_evc, answ_sim, max;
	time_t tm_b, tm_e;
	
	int n;
	double *A = NULL, *E = NULL, *tmp = NULL;

/*	parsing command line*/

	if(argc == 1) { /*work with default options*/
		d = DEFAULT_D_OPTION;
		dd = DEFAULT_DD_OPTION;
		e = DEFAULT_E_OPTION;
		p = DEFAULT_P_OPTION;
		t = DEFAULT_T_OPTION;
		prec = DEFAULT_PREC_OPTION;
		eps = DEFAULT_PREC_OPTION;
		max_iter = DEFAULT_MAXITER_OPTION;
		input_filename = DEFAULT_INPUT_FILENAME;
		output_filename = DEFAULT_OUTPUT_FILENAME;
	}
	if(argc >= 2) { /*parse options*/
		j = 0;/*count for filename parameters*/
		for (i = 1; i < argc; i++) {
			if ((argv[i][0] != '-') && (i < 4)) { /*work with filename*/
				if (j == 0) { /*input filename*/
					input_filename = argv[i];
					j++;
					continue;
				} else if (j == 1) { /*output filename*/
					output_filename = argv[i];
					j++;
					continue;
				} else { /*wrong parameter*/
					fprintf(stderr, "Error: Wrong parameters string - too many non-\'minus\' prefixed parameters!\n\n");
					usage(stderr);
					return 2;
				}
			} else { /*work with options*/
				if ((argv[i][1] == 'm') && (argv[i][2] == 'a') && (argv[i][3] == 'x') && (argv[i][4] == '_') && (argv[i][5] == 'i') && (argv[i][6] == 't') && (argv[i][7] == 'e') && (argv[i][8] == 'r') && (argv[i][9] == '=')) {
					max_iter = atoi(&argv[i][10]);
					continue;
				} else if ((argv[i][1] == 'e') && (argv[i][2] == 'p') && (argv[i][3] == 's') && (argv[i][4] == '=')) {
					eps = atof(&argv[i][5]);
					continue;
				} else if ((argv[i][1] == 'd') && (argv[i][2] == 'd')) {
					dd = 1;
					d = 1;
					continue;
				} else if((argv[i][1] == 'p') && (argv[i][2] == 'r') && (argv[i][3] == 'e') && (argv[i][4] == 'c') && (argv[i][5] == '=')) {
					prec = atof(&argv[i][6]);
					continue;
				}
				if (argv[i][2] != 0) {
					fprintf(stderr, "Error: Unknown option \'%s\'!\n", argv[i]);
					usage(stderr);
					return 2;
				} else if ((argv[i][1] == 'h') || (argv[i][1] == '?')) {
					usage(stdout);
					return 6;
				} else if (argv[i][1] == 'd') {
					d = 1;
					continue;
				} else if (argv[i][1] == 't') {
					t = 1;
					continue;
				} else if (argv[i][1] == 'e') {
					e = 1;
					continue;
				} else if (argv[i][1] == 'p') {
					p = 1;
					continue;
				} else {
					fprintf(stderr, "Error: Unknown option \'%s\'!\n", argv[i]);
					usage(stderr);
					return 2;
				}
			}
		}
	}
/*	if (d) {*/
/*		printf("Parsing command line success!\n");*/
/*		printf("Parametr:	input_filename: %s", input_file_name);*/
/*		printf("		output_filename: %s", output_file_name);*/
/*	}*/


/*	work with files*/

	fin = fopen(input_filename, "r");
	if (fin == NULL) {
		if (e)
			fprintf(stderr, "Error: Couldn\'t open file %s for reading!\n", input_filename);
		return 3;
	}
	if (d)
		printf("Open input file %s success!\n", input_filename);

	fout = fopen(output_filename, "w");
	if (fout == NULL) {
		if (e)
			fprintf(stderr, "Error: Couldn\'t open file %s for writing!\n", output_filename);
		fclose(fin);
		return 4;
	}
	if (d)
		printf("Open output file %s success!\n", output_filename);

/*	reading from input file and allocate memory*/

	fscanf(fin, "%d", &n);
	if (n < 1) {
		if (e)
			fprintf(stderr, "Error: Wrong in file %s, parameter \'n\' = %d, it should be > 1!\n", input_filename, n);
		fclose(fin);
		fclose(fout);
		return 5;
	}
/*	allocate memory*/
	A = (double *)malloc(sizeof(double)*(n*n));
	E = (double *)malloc(sizeof(double)*n);
	if (evc_memsize_10_07(n) > sim_memsize_10_07(n)) {
		tmp = (double*)malloc(evc_memsize_10_07(n));
	} else {
		tmp = (double*)malloc(sim_memsize_10_07(n));
	}
	if (d)
		printf("Allocate memory success!\n");
/*	reading from file*/
	for (i = 0; i < n*n; i++){
	 	fscanf(fin, "%lf", &A[i]);
	}
	if (d)
		printf("Reading information from file %s success!\n", input_filename);
	if ((p==1) || (d==1))
		print_matr(n, A);

/*	execution*/

	tm_b = time(NULL);
	answ_sim = sim_10_07(n, A, tmp, prec);
	if (answ_sim == 0)
		answ_evc = evc_10_07(n, max_iter, eps, A, E, tmp, prec);
	tm_e = time(NULL);
	if (d)
		printf("Execution complite!\n");
	if (t)
		printf("Execution time: %1.9lf\n", difftime(tm_b, tm_e));
	if (p)
		print_matr(n, A);

/*	writing to file*/

	if (answ_evc == 0) {
		fprintf(fout,"%d\n", n);
		for (i = 0; i < n; i++)
		{
			fprintf(fout, "%1.9lf\n", E[i]);
		}
	} else {
		fprintf(fout, "0");
	}
	if (d)
		printf("Writing answer to file %s success!\n", output_filename);

/*	closing all files*/

	fclose(fin);
	fclose(fout);
	if (d)
		printf("Closing file success!\n");

/*	free memory*/

	if (A != NULL) free(A);
	if (E != NULL) free(E);
	if (tmp != NULL) free(tmp);
	if (d)
		printf("Free memory success!\n");

	return answ_evc;
}
