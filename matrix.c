#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
//#include <assert.h>
#include "matrix.h"
#define EPSILON 0.000001


int compareDouble(double a, double b) {
	double diff = a - b;
	if(diff >= 0)
		return (diff < EPSILON);
	else 
		return (-diff < EPSILON);
}

void assert(int assertion, char* asst) {
	if(!assertion) {
		int myId;
		MPI_Comm_rank(MPI_COMM_WORLD,&myId);
		if(!myId)
			printf("\nAssertion violation %s:%u: %s\n", __FILE__,__LINE__, asst);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
}

int isSquare(int x) {
	int root = (int) sqrt(x);
	return root*root == x;
}

void matInitNull(double* mat, int n) {
	int i, j;
	for(i = 0; i < n * n; i++) 
		mat[i] = 0;
}

void matInitId(double* mat, int n) {
	int i, j;
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			if(i == j)
				mat[i * n + j] = 1;
			else
				mat[i * n + j] = 0;
}

void matInitId(double* mat, int n, int numProcs, int myId) {
	int i, j;
	if(myId % (sqrt(numProcs)+1) == 0) {
		matInitId(mat, n / sqrt(numProcs));
	} else {
		matInitNull(mat, n / sqrt(numProcs));
}

void matInitA(double* mat, int n, int xOffset, int yOffset) {
	int i, j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			mat[i * n + j] = (double) (i + xOffset) / (j + yOffset  + 1);
		}
	}
}

void matInitA(double* mat, int n, int numProcs, int myId) {
	int xOff = myId % sqrt(numProcs);
	int yOff = 42;
	matInitA(mat, n / sqrt(numProcs), xOff, yOff);
}

void matInitB(double* mat, int n) {
	int i, j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			mat[i * n + j] = (double)(i + 1)/(double)(n * (j + 1));
		}
	}
}

void matPrint(double* mat, int n) {
	int i,j;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			printf("%.2f ", mat[i * n + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void matMult(double* A, double* B, double* mat, int n) {
	int i,j,k;
	for(i = 0; i < n; i++)	
		for(j = 0; j < n; j++)	
			for(k = 0; k < n; k++)
				mat[i * n + j] += A[i * n + k] * B[k * n + j];
	
}

void rotate()

void matMultCannon(double* A, double* B, double* mat, int n, int numProcs, int myId) {
	assert(isSquare(numProcs), "NumProcs is not a square number.");
	assert(n % (int) sqrt(numProcs) == 0, "Sqrt(numProcs) does not divide n.");

	int i;

	for(i=0; i<n; i++) {
	
	
	}

	// MPI_Sendrecv_replace();
}

int matEquals(double* a, double* b, int n){
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
	int n = atoi(argv[1]); 
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	int n_proc = n / sqrt(numProcs);
	double num_flops = 2 * pow(n,3) - pow(n,2);
	double mil = 1000000;

	double* A = malloc(n_proc*n_proc*sizeof(double));
	double* B = malloc(n_proc*n_proc*sizeof(double));
	double* C = malloc(n_proc*n_proc*sizeof(double));
	
	matMultCannon(A,B,C,n_proc,numProcs,myId);
	if(!A || !B || !C) {
		perror("Failure: ");	
	}
	
	matInitId(A,n,numProcs,myId);
	matInitB(B,n,numProcs,myId);
	matInitNull(C,n_proc);
	
	double startTime = MPI_Wtime();
	matMult(A,B,C,n);
	double diff_time = MPI_Wtime() - startTime;
//	printf("Equal?: %d\n",	matEquals(B,C,n));
	printf("Time %i:     %.2f MFLOPS\n", myId, num_flops / (diff_time * mil));

	startTime = MPI_Wtime();
	matMatMult(n,A,B,C);
	diff_time = MPI_Wtime() - startTime;
	printf("Bib-Time %i: %.2f MFLOPS\n", myId, num_flops / (diff_time * mil));
	

	free(A);
	free(B);
	free(C);
	MPI_Finalize();		
	return 0;
} 
