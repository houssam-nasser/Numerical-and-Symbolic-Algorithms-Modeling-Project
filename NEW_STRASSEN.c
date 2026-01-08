#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define THRESHOLD 1024  // Threshold size for switching to naive multiplication
#define DEBUG 0      // Set to 0 to disable debug prints

// Allocate memory for a matrix stored as a 1D array
double *allocate_matrix(int size) {
    return (double *)malloc(size * size * sizeof(double));
}

// Free allocated memory for a 1D array matrix
void free_matrix(double *matrix) {
    free(matrix);
}

// Add two matrices (1D arrays)
void add_matrix(int size, int row_length, const double *A, const double *B, double *C) {
    assert(size > 0 && row_length >= size && "Invalid matrix dimensions");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int index = i * row_length + j;
            assert(index < size * row_length && "Index out of bounds in add_matrix");
            C[index] = A[index] + B[index];
        }
    }
}

// Subtract two matrices (1D arrays)
void subtract_matrix(int size, int row_length, const double *A, const double *B, double *C) {
    assert(size > 0 && row_length >= size && "Invalid matrix dimensions");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int index = i * row_length + j;
            assert(index < size * row_length && "Index out of bounds in subtract_matrix");
            C[index] = A[index] - B[index];
        }
    }
}

// Naive multiplication for two matrices (1D arrays)
void naive_mult(int size, int row_length_A, const double *A, int row_length_B, const double *B, int row_length_C, double *C) {
    assert(size > 0 && row_length_A >= size && row_length_B >= size && row_length_C >= size && "Invalid matrix dimensions");
    int unroll_factor = (size >= 64) ? 4 : 2;  // Dynamically adjust unroll factor based on matrix size
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double sum = 0.0;
            for (int k = 0; k < size; k += unroll_factor) {
                for (int u = 0; u < unroll_factor && (k + u) < size; u++) {
                    assert(i * row_length_A + k + u < size * row_length_A && "Index out of bounds in naive_mult A");
                    assert((k + u) * row_length_B + j < size * row_length_B && "Index out of bounds in naive_mult B");
                    sum += A[i * row_length_A + k + u] * B[(k + u) * row_length_B + j];
                }
            }
            assert(i * row_length_C + j < size * row_length_C && "Index out of bounds in naive_mult C");
            C[i * row_length_C + j] += sum;
        }
    }
}

// Function to find the next power of 2
int next_power_of_two(int n) {
    int power = 1;
    while (power < n) {
        power *= 2;
    }
    return power;
}

// Function to pad a matrix to the next power of 2, stored in 1D array
void pad_matrix(const double *original, double *padded, int old_size, int new_size) {
    for (int i = 0; i < new_size; i++) {
        for (int j = 0; j < new_size; j++) {
            if (i < old_size && j < old_size)
                padded[i * new_size + j] = original[i * old_size + j];  // Copy original value
            else
                padded[i * new_size + j] = 0.0;  // Pad with zeros
        }
    }
}

// Temporary matrices reused across recursive calls
static double *temp1 = NULL;
static double *temp2 = NULL;
static double *q1 = NULL;
static double *q2 = NULL;
static double *q3 = NULL;
static double *q4 = NULL;
static double *q5 = NULL;
static double *q6 = NULL;
static double *q7 = NULL;

// Strassen Algorithm
void strassen_mult(double *M, double *N, double *R, int size, int stride) {
    assert(size > 0 && stride >= size && "Invalid matrix dimensions");
    if (DEBUG) {
        printf("Strassen_mult: size = %d\n", size);
        fflush(stdout);
    }

    // Base case: Use naive multiplication for small matrices
    if (size <= THRESHOLD) {
        if (DEBUG) {
            printf("Using naive multiplication for size: %d\n", size);
        }
        naive_mult(size, stride, M, stride, N, stride, R);
        return;
    }

    int newSize = size / 2;

    // Pointers to submatrices (blocks)
    double *a = M;
    double *b = M + newSize;
    double *c = M + newSize * stride;
    double *d = M + newSize * stride + newSize;

    double *x = N;
    double *y = N + newSize;
    double *z = N + newSize * stride;
    double *t = N + newSize * stride + newSize;

    double *r11 = R;
    double *r12 = R + newSize;
    double *r21 = R + newSize * stride;
    double *r22 = R + newSize * stride + newSize;

    // Ensure temporary matrices are allocated only once
    if (!temp1) {
        temp1 = (double *)calloc(stride * stride, sizeof(double));
        temp2 = (double *)calloc(stride * stride, sizeof(double));
        q1 = (double *)calloc(stride * stride, sizeof(double));
        q2 = (double *)calloc(stride * stride, sizeof(double));
        q3 = (double *)calloc(stride * stride, sizeof(double));
        q4 = (double *)calloc(stride * stride, sizeof(double));
        q5 = (double *)calloc(stride * stride, sizeof(double));
        q6 = (double *)calloc(stride * stride, sizeof(double));
        q7 = (double *)calloc(stride * stride, sizeof(double));
        if (!temp1 || !temp2 || !q1 || !q2 || !q3 || !q4 || !q5 || !q6 || !q7) {
            fprintf(stderr, "Memory allocation failed.\n");
            exit(EXIT_FAILURE);
        }
    }

    // q1 = a * (x + z)
    add_matrix(newSize, stride, x, z, temp2);
    strassen_mult(a, temp2, q1, newSize, stride);

    // q2 = d * (y + t)
    add_matrix(newSize, stride, y, t, temp2);
    strassen_mult(d, temp2, q2, newSize, stride);

    // q3 = (d - a) * (z - y)
    subtract_matrix(newSize, stride, d, a, temp1);
    subtract_matrix(newSize, stride, z, y, temp2);
    strassen_mult(temp1, temp2, q3, newSize, stride);

    // q4 = (b - d) * (z + t)
    subtract_matrix(newSize, stride, b, d, temp1);
    add_matrix(newSize, stride, z, t, temp2);
    strassen_mult(temp1, temp2, q4, newSize, stride);

    // q5 = (b - a) * z
    subtract_matrix(newSize, stride, b, a, temp1);
    strassen_mult(temp1, z, q5, newSize, stride);

    // q6 = (c - a) * (x + y)
    subtract_matrix(newSize, stride, c, a, temp1);
    add_matrix(newSize, stride, x, y, temp2);
    strassen_mult(temp1, temp2, q6, newSize, stride);

    // q7 = (c - d) * y
    subtract_matrix(newSize, stride, c, d, temp1);
    strassen_mult(temp1, y, q7, newSize, stride);

    // r11 = q1 + q5
    add_matrix(newSize, stride, q1, q5, r11);

    // r12 = q1 + q3 + q6 - q7
    add_matrix(newSize, stride, q1, q3, temp1);
    add_matrix(newSize, stride, temp1, q6, temp2);
    subtract_matrix(newSize, stride, temp2, q7, r12);

    // r21 = q2 + q3 + q4 - q5
    add_matrix(newSize, stride, q2, q3, temp1);
    add_matrix(newSize, stride, temp1, q4, temp2);
    subtract_matrix(newSize, stride, temp2, q5, r21);

    // r22 = q2 + q7
    add_matrix(newSize, stride, q2, q7, r22);

    if (DEBUG) {
        printf("Strassen_mult: completed size = %d\n", size);
        fflush(stdout);
    }
}

// Free reused temporary matrices
void free_temp_matrices() {
    free(temp1);
    free(temp2);
    free(q1);
    free(q2);
    free(q3);
    free(q4);
    free(q5);
    free(q6);
    free(q7);
}

int main() {
    int size;
    printf("Enter the size of the matrices (NxN): ");
    scanf("%d", &size);

    // Calculate the next power of 2 greater than or equal to size
    int padded_size = next_power_of_two(size);
    printf("Adjusted matrix size after padding: %d x %d\n", padded_size, padded_size);

    // Allocate memory for matrices in 1D format
    double *A = (double *)malloc(size * size * sizeof(double));        // Original matrix A
    double *B = (double *)malloc(size * size * sizeof(double));        // Original matrix B
    double *C = (double *)malloc(size * size * sizeof(double));        // Resultant matrix C
    double *A_padded = (double *)calloc(padded_size * padded_size, sizeof(double)); // Padded A
    double *B_padded = (double *)calloc(padded_size * padded_size, sizeof(double)); // Padded B
    double *C_padded = (double *)calloc(padded_size * padded_size, sizeof(double)); // Padded C

    // Check memory allocation
    if (!A || !B || !C || !A_padded || !B_padded || !C_padded) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize matrices A and B with random values
    srand(time(NULL));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A[i * size + j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;  // Random values between -1 and 1
            B[i * size + j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
    }

    printf("Matrix A:\n");
    //printSqMatrix(A, size, size);

    printf("Matrix B:\n");
    //printSqMatrix(B, size, size);

    // Pad matrices to the next power of 2
    pad_matrix(A, A_padded, size, padded_size);
    pad_matrix(B, B_padded, size, padded_size);

    // Perform matrix multiplication using Strassen's Algorithm
    clock_t start = clock();
    strassen_mult(A_padded, B_padded, C_padded, padded_size, padded_size);
    clock_t end = clock();

    // Display the resultant matrix
    printf("Resultant Matrix C:\n");/*
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%8.4f ", C[i * size + j]);
        }
        printf("\n");
    }*/

    printf("Time taken for Strassen's Algorithm: %.6f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    // Free allocated memory
    free(A);
    free(B);
    free(C);
    free(A_padded);
    free(B_padded);
    free(C_padded);

    // Free temporary matrices
    free_temp_matrices();

    return 0;
}
