#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#define EPSILON 0.000001

int compareDouble(double a, double b) {
	double diff = a - b;
	if(diff >= 0)
		return (diff < EPSILON);
	else 
		return (-diff < EPSILON);
}

void matInitNull(double *mat, int n) {
	int i, j;
	for(i = 0; i < n * n; i++) 
		mat[i] = 0;
}

void matInitId(double *mat, int n) {
	int i, j;
	for(i = 0; i < n; i++) 
		for (j = 0; j < n; j++)
			if(i == j)
				mat[i * n + j] = 1;
			else
				mat[i * n + j] = 0;
}

void matInitA(double *mat, int n) {
	int i, j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			mat[i * n + j] = (double) i / (j + 1);
		}
	}
}

void matInitB(double *mat, int n) {
	int i, j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			mat[i * n + j] = (double)(i + 1)/(double)(n * (j + 1));
		}
	}
}

void matPrint(double *mat, int n) {
	int i,j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			printf("%.2f ", mat[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void matMult(double *A, double *B, double *mat, int n) {
	int i,j,k;
	for(i = 0; i < n; i++)	
		for(j = 0; j < n; j++)	
			for(k = 0; k < n; k++)
				mat[i * n + j] += A[i * n + k] * B[k * n + j];
	
}

void matMultCannon(double *A, double *B, double *mat, int n, int numProcs, int myId) {
	assert ((n % (int) sqrt((double) numProcs)) == 0);

	// MPI_Sendrecv_replace();
	printf("TODO");
}

int matEquals(double *a, double *b, int n){
	int i,j,k;
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			if( !compareDouble(a[i * n + j], b[i * n + j]) ) {
				printf("a:%f -  b: %f\n", a[i * n + j], b[i * n + j]);
				return 0;
			}
	return 1;
}

int main(int argc, char** argv) {
	int numProcs, myId;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	int n = 1000;
	double num_flops = 2 * pow(n,3) - pow(n,2);
	double mil = 1000000;

	double* A = malloc(n*n*sizeof(double));
	double* B = malloc(n*n*sizeof(double));
	double* C = malloc(n*n*sizeof(double));

	matInitId(A,n);
	matInitB(B,n);

	double startTime = MPI_Wtime();
	matMult(A,B,C,n);
	double diff_time = MPI_Wtime() - startTime;
	printf("Equal?: %d\n",	matEquals(B,C,n));
	printf("Time:     %.2f MFLOPS\n", num_flops / (diff_time * mil));

	startTime = MPI_Wtime();
	matMatMult(n,A,B,C);
	diff_time = MPI_Wtime() - startTime;
	printf("Bib-Time: %.2f MFLOPS\n", num_flops / (diff_time * mil));
	

	MPI_Finalize();		
	return 0;
} 
