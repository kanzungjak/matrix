#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "matrix.h"
#define EPSILON 0.000001

int myId;
int numProcs;
int blocksPerDim;

void assert(int assertion, char* errMsg) {
	if(!assertion) {
		if(!myId)
			printf("\nAssertion violation %s:%u: %s\n", __FILE__,__LINE__, errMsg);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
}

int nProc(int n) { // Width (and height) for each process
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (n / blocksPerDim);
}

int xPos() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (myId % blocksPerDim);
}

int yPos() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (myId / blocksPerDim);
}

int compareDouble(double a, double b) {
	double diff = a - b;
	if(diff >= 0)
		return (diff < EPSILON);
	else 
		return (-diff < EPSILON);
}

int isSquare(int x) {
	int root = (int) sqrt(x);
	return root * root == x;
}

// Null Matrix (local)
void matLocInitNull(double* mat, int nProc) {
	int i;
	for(i = 0; i < nProc * nProc; i++) 
		mat[i] = 0;
}

// Null matrix (global)
void matInitNull(double* mat, int n) {
	matLocInitNull(mat, nProc(n));
}

// Identity matrix (local)
void matLocInitId(double* mat, int nProc) {
	int i, j;
	for(i = 0; i < nProc; i++)
		for(j = 0; j < nProc; j++)
			if(i == j)
				mat[i * nProc + j] = 1;
			else
				mat[i * nProc + j] = 0;
}

// Identity Matrix (global)
void matInitId(double* mat, int n) {
	if(xPos() == yPos()) {
		matLocInitId(mat, nProc(n));
	}
	else {
		matLocInitNull(mat, nProc(n));
	}
}

// A[i,j] = i / (i + j) (local)
void matLocInitA(double* mat, int nProc, int xOff, int yOff) {
	int i, j;
	xOff *= nProc;
	yOff *= nProc;
	for(i = 0; i < nProc; i++) {
		for(j = 0; j < nProc; j++) {
			mat[i * nProc + j] = (double) (i + xOff) / (j + yOff  + 1);
		}
	}
}

// A[i,j] = i / (i + j) (global)
void matInitA(double* mat, int n) {
	matLocInitA(mat, nProc(n), xPos(), yPos());
}


// A[i,j] = (i + 1) / (n * (j + 1)) (local)
void matLocInitB(double* mat, int nProc, int xOff, int yOff) {
	int i, j;
	xOff *= nProc;
	yOff *= nProc;
	for(i = 0; i < nProc; i++) {
		for(j = 0; j < nProc; j++) {
			mat[i * nProc + j] = (double) (i+1 + xOff) / (nProc * (j + yOff + 1));
		}
	}
}

// A[i,j] = (i + 1) / (n * (j + 1)) (global)
void matInitB(double* mat, int n) {
	matLocInitB(mat, nProc(n), xPos(), yPos());
}

// Print entries of a matrix, not a good idea for big matrices
void matPrint(double* mat, int n) {
	int i,j;
	double* complete = malloc(n * n * sizeof(double));
	if(numProcs > 1) {
		/*MPI_Gather(mat,    nProc(n) * nProc(n), MPI_DOUBLE, 
			   complete, nProc(n) * nProc(n), MPI_DOUBLE,
			   0, MPI_COMM_WORLD);*/
		MPI_Comm line, row;
		MPI_Comm_split(MPI_COMM_WORLD, yPos(), myId, &line);
		for(i = 0; i < nProc(n); i++) {
			MPI_Gather(&mat[i * nProc(n)], nProc(n), MPI_DOUBLE,
		    		   &complete[i * n], nProc(n), MPI_DOUBLE,
				   0, line);
		}
		MPI_Comm_free(&line);

		MPI_Comm_split(MPI_COMM_WORLD, xPos(), myId, &row);
		MPI_Gather(complete, n * nProc(n), MPI_DOUBLE,
	     	     	   complete, n * nProc(n), MPI_DOUBLE,
			   0, row);  // Hmm, besser, aber immer noch nicht ganz richtig glaubich ...
		MPI_Comm_free(&row);
	} else {
		complete = mat;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if(!myId) {
		for(i = 0; i < n; i++) {
			for(j = 0; j < n; j++) {
				if (complete[i * n + j] < 10)
					printf("%.2f ", complete[i * n + j]);
				else printf("E.00 ");
			}
			printf("\n");
		}
		printf("\n");
	}
}

// Naive matrix multiplication
void matMult(double* A, double* B, double* mat, int n) {
	int i,j,k;
	// Wollen wir das wirklich so machen?
	// Muss man doch irgendwie zwischen den Prozessen kommunizieren ?! -> sehr oft sendrecv vermutl.
	// Momentan machen wir eher das Folgende (mat) für myProcs = 4:
	// mat = |A0*B0|A1*B1| =!= A*B = |A0*B0+A1*B2|A0*B1+A1*B3|
	//       |A2*B2|A3*B3|           |A2*B0+A3*B2|A2*B1+A3*B3|
	for(i = 0; i < nProc(n); i++)
		 for(j = 0; j < nProc(n); j++)
			for(k = 0; k < nProc(n); k++)
				mat[i * nProc(n) + j] += A[i * nProc(n) + k] * B[k * nProc(n) + j];
}

void matAdd(double* A, double* B, double* mat, int n) {
	int i;
	for(i = 0; i < nProc(n) * nProc(n); i++)	
		mat[i] = A[i] + B[i];
}

void matMultAdd(double* A, double* B, double* mat, int n ) {
	double* C = malloc(nProc(n) * nProc(n) * sizeof(double));
	matMult(A, B, C, n);
	matAdd(mat, C, mat, n);
	free(C);
}

int isLeftMost() {
	return xPos() == 0;
}

int isRightMost() {
	return (myId + 1) % blocksPerDim == 0;
}

void rotateLeft(double* mat, int n) {
	int left_nb = isLeftMost() ? (myId + blocksPerDim - 1) : (myId - 1);
	int right_nb = isRightMost() ? (myId - blocksPerDim + 1) : (myId + 1);
	MPI_Status status;
	MPI_Sendrecv_replace(mat, nProc(n)*nProc(n), MPI_DOUBLE,
			     left_nb, 37, right_nb, MPI_ANY_TAG,
			     MPI_COMM_WORLD, &status);
}

void rotateUp(double* mat, int n) {
	int upper_nb = (myId - blocksPerDim + numProcs) % numProcs;
	int lower_nb = (myId + blocksPerDim) % numProcs;
	MPI_Status status;
	MPI_Sendrecv_replace(mat, nProc(n)*nProc(n), MPI_DOUBLE,
			     upper_nb, 38, lower_nb, MPI_ANY_TAG,
			     MPI_COMM_WORLD, &status);
}

void matMultCannon(double* A, double* B, double* mat, int n) {
	assert(isSquare(numProcs), "numProcs is not a square number.");
	assert(n % blocksPerDim == 0, "Sqrt(numProcs) does not divide n.");

	int i;
	matInitNull(mat, n);

	for(i = 0; i < yPos(); i++){ // y-dim
		rotateLeft(A, n);
	}
	for(i = 0; i < xPos(); i++){ // x-dim
		rotateUp(B, n);	
	}

	for (i = 0; i < blocksPerDim; i++) {
		matMultAdd(A, B, mat, n);
		rotateLeft(A, n);
		rotateUp(B, n);
	}
}


//0: First entry, which is unequal
//1: A = B
int matEquals(double* A, double* B, int n){
	int i,j,k;
	for(i = 0; i < nProc(n); i++)
		for(j = 0; j < nProc(n); j++)
			if( !compareDouble(A[i * nProc(n) + j], B[i * nProc(n) + j]) ) {
				printf("a: %f -  b: %f index: (%d,%d)\n", A[i * nProc(n) + j], B[i * nProc(n) + j], i, j);
				return 0;
			}
	return 1;
}

int main(int argc, char** argv) {
	int n = atoi(argv[1]); // Width (and height) of a matrix
//	printf ("%d\n", argc);
	if (argc < 2) {
		printf("Usage: matrix 'length'");
		exit(EXIT_FAILURE);
	}	

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	
	blocksPerDim = (int) sqrt(numProcs);

	double numFlops = 2 * pow(n,3) - pow(n,2); // Number of flops for naive matrix-multiplication
	

	double* A = malloc(nProc(n) * nProc(n) * sizeof(double));
	double* B = malloc(nProc(n) * nProc(n) * sizeof(double));
	double* C = malloc(nProc(n) * nProc(n) * sizeof(double));
	double* D = malloc(nProc(n) * nProc(n) * sizeof(double));
	double* E = malloc(nProc(n) * nProc(n) * sizeof(double));
	double* F = malloc(nProc(n) * nProc(n) * sizeof(double));
	if(!A || !B || !C || !D || !E) {
		perror("Failure: ");	
	}
	
	matInitA(A, n);
	matInitB(B, n);
	matInitNull(C, n);
	matInitId(D, n);
	matInitNull(E, n);
	matInitNull(F, n);

	//double startTime = MPI_Wtime();
	matMult(D, D, F, n); 
	matMultCannon(D, D, C, n);
	//double diff_time = MPI_Wtime() - startTime;
//	printf("Equal?: %d\n",	matEquals(B,C,n));
	//printf("Time %i:     %.2f MFLOPS\n", myId, numFlops / (diff_time * 1000000.0));

	//startTime = MPI_Wtime();
//	matMatMult(nProc(n), D, D, E);
	//diff_time = MPI_Wtime() - startTime;
	//printf("Bib-Time %i: %.2f MFLOPS\n", myId, numFlops / (diff_time * 1000000.0));
	printf("%d Eins\n", myId);
	//MPI_Barrier(MPI_COMM_WORLD);

	if(!myId) printf("Null\n");
	matPrint(E, n);
	printf("%d Zwei\n", myId);
	//MPI_Barrier(MPI_COMM_WORLD);
	if(!myId) printf("B\n");
	matPrint(B, n);
	if(!myId) printf("A\n");
	matPrint(A, n);
	if(!myId) printf("Id\n");
	matPrint(D, n);

	if(!myId) printf("Cannon\n");
	matPrint(C, n);
//	if(!myId) printf("Library multiplication\n");
//	matPrint(E, n);
	if(!myId) printf("Naive multiplication\n");
	matPrint(F, n);

	//MPI_Barrier(MPI_COMM_WORLD);

	printf("\n%d Equals? %d \n",myId, matEquals(C, F, n));

	free(A);
	free(B);
	free(C);
	free(D);
	free(E);
	free(F);
	MPI_Finalize();		
	return 0;
}
