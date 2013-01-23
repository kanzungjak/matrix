#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "matrix.h"
#define EPSILON 0.000001

int myId, numProcs, blocksPerDim;

// Assert with error message
void assert(int assertion, char* errMsg) {
	if(!assertion) {
		if(!myId)
			printf("\nAssertion violation %s:%u: %s\n", __FILE__,__LINE__, errMsg);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
}

// Width (and height) for each PE
int nProc(int n) {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (n / blocksPerDim);
}

// Position of local matrix block in total matrix (horizontally from left to right)
int xPos() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (myId % blocksPerDim);
}

// Position of local matrix block in total matrix (vertically from top to bottom)
int yPos() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (myId / blocksPerDim);
}

// Is the local matrix block on the left side?
int isLeftMost() {
	return xPos() == 0;
}

// Is the local matrix block on the right side?
int isRightMost() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (myId + 1) % blocksPerDim == 0;
}

// To which PE does the block to the left belong?
int leftNb() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (isLeftMost() ? (myId + blocksPerDim - 1) : (myId - 1));
}

// To which PE does the block to the right belong?
int rightNb() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	return (isRightMost() ? (myId - blocksPerDim + 1) : (myId + 1));
}

// To which PE does the block above belong?
int upperNb() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	assert(numProcs > 0, "numProcs is not calculated yet.");
	return ((myId - blocksPerDim + numProcs) % numProcs);
}

// To which PE does the block below belong?
int lowerNb() {
	assert(blocksPerDim > 0, "blocksPerDim is not calculated yet.");
	assert(numProcs > 0, "numProcs is not calculated yet.");
	return ((myId + blocksPerDim) % numProcs);
}

// Since floating point values are rounded, we do not test on exact equality
int compareDouble(double a, double b) {
	double diff = a - b;
	if(diff >= 0)
		return (diff < EPSILON);
	else 
		return (-diff < EPSILON);
}

// Is the number a square number?
int isSquare(int x) {
	int root = (int) sqrt(x);
	return root * root == x;
}

// Null Matrix (local)
void matLocInitNull(double* mat, int nLoc) {
	int i;
	for(i = 0; i < nLoc * nLoc; i++) 
		mat[i] = 0;
}

// Null matrix (global)
void matInitNull(double* mat, int n) {
	matLocInitNull(mat, nProc(n));
}

// Identity matrix (local)
void matLocInitId(double* mat, int nLoc) {
	int i, j;
	for(i = 0; i < nLoc; i++)
		for(j = 0; j < nLoc; j++)
			if(i == j)
				mat[i * nLoc + j] = 1;
			else
				mat[i * nLoc + j] = 0;
}

// Identity Matrix (global)
void matInitId(double* mat, int n) {
	if(xPos() == yPos())
		matLocInitId(mat, nProc(n));
	else 
		matLocInitNull(mat, nProc(n));
}

// A[i,j] = i / (j + 1) (local)
void matLocInitA(double* mat, int nLoc, int xOff, int yOff) {
	int i, j;
	xOff *= nLoc;
	yOff *= nLoc;
	for(i = 0; i < nLoc; i++)
		for(j = 0; j < nLoc; j++)
			mat[i * nLoc + j] = (double) (i + yOff) / (j + xOff  + 1);	
}

// A[i,j] = i / (j + 1) (global)
void matInitA(double* mat, int n) {
	matLocInitA(mat, nProc(n), xPos(), yPos());
}


// B[i,j] = (i + 1) / (n * (j + 1)) (local)
void matLocInitB(double* mat, int nLoc, int xOff, int yOff) {
	int i, j;
	xOff *= nLoc;
	yOff *= nLoc;
	for(i = 0; i < nLoc; i++) 
		for(j = 0; j < nLoc; j++) 
			mat[i * nLoc + j] = (double) (i+1 + yOff) / (nLoc * blocksPerDim * (j + xOff + 1));	
}

// B[i,j] = (i + 1) / (n * (j + 1)) (global)
void matInitB(double* mat, int n) {
	matLocInitB(mat, nProc(n), xPos(), yPos());
}

// Print entries of whole matrix, not a good idea for big matrices
void matPrint(double* mat, int n) {
	int i,j, nLoc;
	double* complete = malloc(n * n * sizeof(double));
	matLocInitNull(complete, n);
	if(numProcs > 1) {
		nLoc = nProc(n);
		MPI_Comm line, row;
		MPI_Comm_split(MPI_COMM_WORLD, yPos(), myId, &line);
		for(i = 0; i < nLoc; i++) {
			MPI_Gather(&mat[i * nLoc], nLoc, MPI_DOUBLE,
		   		   &complete[i * n], nLoc, MPI_DOUBLE,
				   0, line);
		}
		MPI_Comm_free(&line);

		MPI_Comm_split(MPI_COMM_WORLD, xPos(), myId, &row);

		if( !yPos() )
			MPI_Gather(MPI_IN_PLACE, n * nLoc, MPI_DOUBLE,
				   complete, n * nLoc, MPI_DOUBLE,
				   0, row);
		else
			MPI_Gather(complete, n * nLoc, MPI_DOUBLE,
				   complete, n * nLoc, MPI_DOUBLE,
				   0, row);
		MPI_Comm_free(&row);
	} else {
		complete = mat;
	}
	MPI_Barrier(MPI_COMM_WORLD);
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

	if(!complete)
		free(complete);
}

// Naive matrix multiplication: mat = A * B (local)
void matLocMult(double* A, double* B, double* mat, int nLoc) {
	int i, j, k;
	for(i = 0; i < nLoc; i++)
		for(j = 0; j < nLoc; j++)
			for(k = 0; k < nLoc; k++)
				mat[i * nLoc + j] += A[i * nLoc + k] * B[k * nLoc + j];
}

// Naive matrix sum: mat = A + B (local)
void matLocAdd(double* A, double* B, double* mat, int nLoc) {
	int i;
	for(i = 0; i < nLoc * nLoc; i++)
		mat[i] = A[i] + B[i];
}

// Naive matrix sum: mat = A + B (global)
void matAdd(double* A, double* B, double* mat, int n) {
	matLocAdd(A, B, mat, nProc(n));
}

// Naive matrix multiplication and sum: mat = mat + A * B (local)
void matLocMultAdd(double* A, double* B, double* mat, int nLoc) {
	double* C = malloc(nLoc * nLoc * sizeof(double));
	matLocInitNull(C, nLoc);
	matLocMult(A, B, C, nLoc);
	matLocAdd(mat, C, mat, nLoc);
	free(C);
}

// Block rotation to the left
void locRotateLeft(double* mat, int nLoc) {
	int left, right;
	MPI_Status status;
	
	left = leftNb();
	right = rightNb();
	
	MPI_Sendrecv_replace(mat, nLoc * nLoc, MPI_DOUBLE,
			     left, 37, right, MPI_ANY_TAG,
			     MPI_COMM_WORLD, &status);
}

// Block rotation upwards
void locRotateUp(double* mat, int nLoc) {
	int upper, lower;
	MPI_Status status;
	
	upper = upperNb();
	lower = lowerNb();
	
	MPI_Sendrecv_replace(mat, nLoc * nLoc, MPI_DOUBLE,
			     upper, 38, lower, MPI_ANY_TAG,
			     MPI_COMM_WORLD, &status);
}

// Good Cannon Matrix Multiplication
void matMultCannon(double* A, double* B, double* mat, int n) {
	assert(isSquare(numProcs), "numProcs is not a square number.");
	assert(n % blocksPerDim == 0, "Sqrt(numProcs) does not divide n.");

	int i, nLoc;
	nLoc = nProc(n);
	matInitNull(mat, n);

	for(i = 0; i < yPos(); i++){ // y-dim
		locRotateLeft(A, nLoc);
	}
	for(i = 0; i < xPos(); i++){ // x-dim
		locRotateUp(B, nLoc);	
	}
	for (i = 0; i < blocksPerDim; i++) {
		matMatMultAdd(nLoc, A, B, mat);
		locRotateLeft(A, nLoc);
		locRotateUp(B, nLoc);
	}
}

// Naive (concerning local multiplicatiuon and sum) Cannon Matrix Multiplication
void matMultNaiveCannon(double* A, double* B, double* mat, int n) {
	assert(isSquare(numProcs), "numProcs is not a square number.");
	assert(n % blocksPerDim == 0, "Sqrt(numProcs) does not divide n.");

	int i, nLoc;
	nLoc = nProc(n);
	matInitNull(mat, n);

	for(i = 0; i < yPos(); i++){ // y-dim
		locRotateLeft(A, nLoc);
	}
	for(i = 0; i < xPos(); i++){ // x-dim
		locRotateUp(B, nLoc);	
	}
	for (i = 0; i < blocksPerDim; i++) {
		matLocMultAdd(A, B, mat, nLoc);
		locRotateLeft(A, nLoc);
		locRotateUp(B, nLoc);
	}
}

// Are all values of our matrix block equal? (local)
// 1 if yes, 0 if not (we might print out the first unequal value we find)
int matLocEquals(double* A, double* B, int nLoc){
	int i, j, k;
	for(i = 0; i < nLoc; i++)
		for(j = 0; j < nLoc; j++)
			if( !compareDouble(A[i * nLoc + j], B[i * nLoc + j]) ) {
				//printf("a: %f -  b: %f index: (%d,%d)\n", A[i * nLoc + j], B[i * nLoc + j], i, j);
				return 0;
			}
	return 1;
}

// Are all matrix values equal? (global)
// 1 if yes, 0 if not (each process might print out the first unequal value it finds)
int matEquals(double* A, double* B, int n){
	int e, nLoc, eq;
	nLoc = nProc(n);
	e = matLocEquals(A, B, nLoc);
	MPI_Reduce(&e, &eq, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	return eq;
}

int main(int argc, char** argv) {
	int n, nLoc, eq; //length of matrix, local length for a PE, return value of matEquals (for root)
	double numFlops, startTime, diffTime; //estimated number of FLOPS, timing stuff
	double *A, *B, *C; //matrices

	n = atoi(argv[1]); // Width (and height) of a matrix
	if (argc < 2) {
		printf("Usage: matrix 'length'");
		exit(EXIT_FAILURE);
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);

	blocksPerDim = (int) sqrt(numProcs); //number of PEs for one dimension
	nLoc = nProc(n);

	numFlops = 2 * pow(n,3) - pow(n,2); // Number of flops for naive matrix multiplication


	A = malloc(nLoc * nLoc * sizeof(double));
	matInitNull(A, n);
	B = malloc(nLoc * nLoc * sizeof(double));
	matInitNull(B, n);
	C = malloc(nLoc * nLoc * sizeof(double));
	matInitNull(C, n);
	if(!A || !B || !C) {
		perror("Failure: ");	
	}

	/*Good Cannon Matrix Multiplication*/
	matInitA(A, n);
	matInitB(B, n);
	matInitNull(C, n);
	if (!myId)
		printf("Good Cannon Matrix Multiplication:\n");
	startTime = MPI_Wtime();
	matMultCannon(A, B, C, n);
	diffTime = MPI_Wtime() - startTime;

	/*if (!myId)
		printf("Result:\n");
	matPrint(C, n);*/

	if (!myId) {
		printf("Operations per Time: %.2f MFLOPS\n", numFlops / (diffTime * 1000000.0));
		printf("Time: %.2f seconds\n", diffTime);
	}

	matInitA(A, n);;
	eq = matEquals(C, A, n);
	if (!myId && eq)
		printf("Correct Result!\n");
	else if (!myId)
		printf("Incorrect Result!\n");

	if (!myId)
		printf("---------------------------------------------------------------\n");

	/*Naive Cannon Matrix Multiplication*/
	matInitA(A, n);
	matInitB(B, n);
	matInitNull(C, n);
	if (!myId)
		printf("Naive Cannon Matrix Multiplication:\n");
	startTime = MPI_Wtime();	
	matMultNaiveCannon(A, B, C, n);
	diffTime = MPI_Wtime() - startTime;

	/*if (!myId)
		printf("Result:\n");
	matPrint(C, n);*/

	if (!myId) {
		printf("Operations per Time: %.2f MFLOPS\n", numFlops / (diffTime * 1000000.0));
		printf("Time: %.2f seconds\n", diffTime);
	}

	matInitA(A, n);
	eq = matEquals(C, A, n);
	if (!myId && eq)
		printf("Correct Result!\n");
	else if (!myId)
		printf("Incorrect Result!\n");

	free(A);
	free(B);
	free(C);
	MPI_Finalize();
	return 0;
}
