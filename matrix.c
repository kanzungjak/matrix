#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "matrix.h"
#define EPSILON 0.000001

int myId;
int numProcs;

int compareDouble(double a, double b) {
	double diff = a - b;
	if(diff >= 0)
		return (diff < EPSILON);
	else 
		return (-diff < EPSILON);
}

void assert(int assertion, char* errMsg) {
	if(!assertion) {
		int myId;
		MPI_Comm_rank(MPI_COMM_WORLD,&myId);
		if(!myId)
			printf("\nAssertion violation %s:%u: %s\n", __FILE__,__LINE__, errMsg);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
}

int isSquare(int x) {
	int root = (int) sqrt(x);
	return root*root == x;
}

//Null Matrix
void matInitNull(double* mat, int n) {
	int i, j;
	for(i = 0; i < n * n; i++) 
		mat[i] = 0;
}

//Identity matrix (local)
void matLocInitId(double* mat, int n) {
	int i, j;
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			if(i == j)
				mat[i * n + j] = 1;
			else
				mat[i * n + j] = 0;
}

//Identity Matrix (global)
void matInitId(double* mat, int n) {
	if((myId % ((int) (sqrt(numProcs)+1))) == 0)  //on the diagonal?
		matLocInitId(mat, n / sqrt(numProcs));
	else 
		matInitNull(mat, n / sqrt(numProcs));
}

//a_ij = i/(i+j) (local)
void matLocInitA(double* mat, int n, int xOffset, int yOffset) {
	int i, j;
	xOffset *= n;
	yOffset *= n;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			mat[i * n + j] = (double) (i + xOffset) / (j + yOffset  + 1);
		}
	}
}

//a_ij = i/(i+j) (global)
void matInitA(double* mat, int n) {
	int xOff = myId % ((int) sqrt(numProcs));
	int yOff = myId / sqrt(numProcs);
	matLocInitA(mat, n / sqrt(numProcs), xOff, yOff);
}


//a_ij = (i+1)/n*(j+1) (local)
void matLocInitB(double* mat, int n, int xOffset, int yOffset) {
	int i, j;
	xOffset *= n;
	yOffset *= n;
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			mat[i * n + j] = (double) (i+1 + xOffset) / (n*(j + yOffset  + 1));
		}
	}
}

//a_ij = (i+1)/n*(j+1) (global)
void matInitB(double* mat, int n) {
	int xOff = myId % ((int) sqrt(numProcs));
	int yOff = myId / sqrt(numProcs);
	matLocInitB(mat, n / sqrt(numProcs), xOff, yOff);
}


//print entries of a matrix, not a good idea for big matrices
void matPrint(double* mat, int n, int numProcs, int myId) {
	int i,j,k,l,m;
	int blocks_per_dim = (int) sqrt(numProcs);
	if(numProcs > 1) {
		double* complete = malloc(n*n*numProcs*sizeof(double));
		//MPI_Gather(mat,      n*n, MPI_DOUBLE, 
		//	   complete, n*n, MPI_DOUBLE,
		//	   0, MPI_COMM_WORLD);

		MPI_Status status;
		for(k = 0; k < blocks_per_dim; k++) {
			for(j = 0; j < n; j++) {
				for(i = 0; i < blocks_per_dim; i++) {
					l = i + k * blocks_per_dim;
					m = i * n + (k + 1) * j * n * blocks_per_dim;
					MPI_Sendrecv(&mat[j], n,  MPI_DOUBLE, 0, i,
				     	     	     &complete[m], n, MPI_DOUBLE, l, i,
				     	     	     MPI_COMM_WORLD, &status);
				}
			}
		}

		n *= blocks_per_dim;
		if(!myId) {
			for(i = 0; i < n; i++) {
				for(j = 0; j < n; j++) {
					printf("%.2f ", complete[i * n + j]);
				}
				printf("\n");
			}
			printf("\n");
		}
	} else {
		if(!myId) {
			for(i = 0; i < n; i++) {
				for(j = 0; j < n; j++) {
					printf("%.2f ", mat[i * n + j]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}
}

//naive matrix multiplication
void matMult(double* A, double* B, double* mat, int n) {
	int i,j,k;
	for(i = 0; i < n; i++)	
		for(j = 0; j < n; j++)	
			for(k = 0; k < n; k++)
				mat[i * n + j] += A[i * n + k] * B[k * n + j];
	
}

void rotate_left(double* mat, int n) {
	int i,j;
	double* tmp;
	int left_nb;
	tmp = malloc(n*sizeof(double)); //one row
	//save first row temporally
	for(i = 0; i < n; i++ )
		tmp[i] = mat[i*n];
	for(i = 0; i < n; i++) { //all columns
		for(j = 1; j < n-1; j++) { //all rows except first and last one (which will be sent (received) to (from) others PEs)
			mat[i*n + j] = mat[i*n+j + 1];
		}
	}

	//MPI_Sendrecv_replace(void* buf, int count, MPI_Datatype datatype,
	//		     int dest, int sendtag, int source, int recvtag, MPI_Comm comm,
        // 	     	     MPI_Status *status);
	MPI_Status status;
	for(i = 0; i < n; i++) {
		MPI_Sendrecv_replace(&mat[i*n], 1, MPI_DOUBLE,
			   	     left_nb, MPI_ANY_TAG, myId, MPI_ANY_TAG, 
				     MPI_COMM_WORLD, &status);
	}

//	free(tmp);
}

void rotate_up() {

}

void matMultCannon(double* A, double* B, double* mat, int n, int numProcs, int myId) {
	assert(isSquare(numProcs), "numProcs is not a square number.");
	assert(n % (int) sqrt(numProcs) == 0, "Sqrt(numProcs) does not divide n.");

	int i,j,k;
	i = 42;
	j = i;
	int index_mod;
	for(k = 0; k < n*n; k++) { //bis n? wohl eher n*n
//		c[i,j] = c[i,j] + a[i,(i+j+k)%n] * b[(i+j+k)%n,j];
		index_mod = i+j+k % n;
		mat[i*n + j] += A[i*n + index_mod]*B[index_mod*n + j]; 	
	}

	// MPI_Sendrecv_replace();
}


//0: first entry, which is unequal
//1: A = B
int matEquals(double* a, double* b, int n){
	int i,j,k;
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			if( !compareDouble(a[i * n + j], b[i * n + j]) ) {
				printf("a:%f -  b: %f index: (%d,%d)\n", a[i*n + j], b[i*n + j], i, j);
				return 0;
			}
	return 1;
}

int main(int argc, char** argv) {
	int n = atoi(argv[1]); //width (and height) of a matrix

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);

	int n_proc = n / sqrt(numProcs); //size of a block, for 1 processor
	double num_flops = 2 * pow(n,3) - pow(n,2); //number of flops for naive matrix-multiplication

	double* A = malloc(n_proc*n_proc*sizeof(double));
	double* B = malloc(n_proc*n_proc*sizeof(double));
	double* C = malloc(n_proc*n_proc*sizeof(double));
	double* D = malloc(n_proc*n_proc*sizeof(double));
	if(!A || !B || !C || !D) {
		perror("Failure: ");	
	}
	
	matInitId(D, n);
	matInitA(A, n);
	matInitB(B, n);
	matInitNull(C, n_proc);
	
	//matMultCannon(A,B,C,n_proc,numProcs,myId);
	
	//double startTime = MPI_Wtime();
	//matMult(A,B,C,n);
	//double diff_time = MPI_Wtime() - startTime;
//	printf("Equal?: %d\n",	matEquals(B,C,n));
	//printf("Time %i:     %.2f MFLOPS\n", myId, num_flops / (diff_time * 1000000.0));

	//startTime = MPI_Wtime();
	//matMatMult(n,A,B,C);
	//diff_time = MPI_Wtime() - startTime;
	//printf("Bib-Time %i: %.2f MFLOPS\n", myId, num_flops / (diff_time * 1000000.0));

	if(!myId) printf("null\n");
	matPrint(C, n_proc, numProcs, myId);
	if(!myId) printf("B\n");
	matPrint(B, n_proc, numProcs, myId);
	if(!myId) printf("A\n");
	matPrint(A, n_proc, numProcs, myId);
	if(!myId) printf("Id\n");
	matPrint(D, n_proc, numProcs, myId);

	free(A);
	free(B);
	free(C);
	free(D);
	MPI_Finalize();		
	return 0;
} 
