#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Function to allocate memory for a matrix
double **allocate_matrix(int size)
{
    double **matrix = (double **)malloc(size * sizeof(double *));
    for (int i = 0; i < size; i++)
    {
        matrix[i] = (double *)malloc(size * sizeof(double));
    }
    return matrix;
}

// Function to free memory of a matrix
void free_matrix(double **matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

// Function for naive matrix multiplication
void naive_matrix_mult(double **A, double **B, double **C, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            C[i][j] = 0.0; // Initialize result matrix element
            for (int k = 0; k < size; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main()
{
    int size;

    printf("Enter matrix size (must be a power of 2): ");
    scanf("%d", &size);

    // Allocate matrices
    double **M = allocate_matrix(size);
    double **N = allocate_matrix(size);
    double **R = allocate_matrix(size);

    // Measure matrix initialization time
    clock_t init_start = clock();
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            M[i][j] = (0.3) * i + (0.7) * j;
            N[i][j] = (0.5) * i - (0.9) * j;
        }
    }
    clock_t init_end = clock();

    // Measure matrix multiplication time
    clock_t mult_start = clock();
    naive_matrix_mult(M, N, R, size);
    clock_t mult_end = clock();

    // Calculate and print times
    double init_time = (double)(init_end - init_start) / CLOCKS_PER_SEC;
    double mult_time = (double)(mult_end - mult_start) / CLOCKS_PER_SEC;

    printf("Matrix initialization time: %.6f seconds\n", init_time);
    printf("Matrix multiplication time: %.6f seconds\n", mult_time);

    // Free allocated memory
    free_matrix(M, size);
    free_matrix(N, size);
    free_matrix(R, size);

    return 0;
}
