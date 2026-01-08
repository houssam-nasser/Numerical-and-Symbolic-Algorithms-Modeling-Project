/*
Group 03

Hani Abdallah - 21400302
Houssam Eddine Jamil Nasser - 21400407
Tan Viet Nguyen - 21400381

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// LU Decomposition function - using gaussian elimination
void LU_Decomposition(double **A, int n, double **L, double **U)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double factor = U[j][i] / U[i][i];
            L[j][i] = factor;
            for (int k = i; k < n; k++)
            {
                U[j][k] -= factor * U[i][k];
            }
        }
    }
}

// Print a 2D matrix
void printSqMatrix(const char *name, double **matrix, int n)
{
    printf("%s:\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.3f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function that request user to input Matrix elements
void RequestInput(const char *name, double **matrix, int n)
{
    printf("Input Matrix %s Elements:\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%s[%d][%d]=", name, i, j);
            scanf("%lf", &matrix[i][j]);
        }
    }
}

// free allocated memory of matrices
void freeMatrix(double **matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

// allocate memory for matrices
double **allocateMatrix(int n)
{
    double **matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        matrix[i] = (double *)malloc(n * sizeof(double));
    }
    return matrix;
}

// Build an Identity matrix
void makeIdentity(double **matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                matrix[i][j] = 1.0;
            }
            else
            {
                matrix[i][j] = 0.0;
            }
        }
    }
}

// Copy Matrix A to Matrix B
void copyMatrix(double **original_matrix, double **copied_matrix, int n)
{

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copied_matrix[i][j] = original_matrix[i][j];
        }
    }
}

int main()
{
    clock_t start_time, end_time;
    int n;

    // Step 1: Get dimensions for Matrix A
    printf("\nChoose Matrix Dimension for the square matrix: ");
    scanf("%d", &n);

    // Error condition
    if (n <= 0)
    {
        printf("Matrix Size can only be larger than zero, the program will exit...");
        exit(1);
    }

    // Step 2: Allocate matrices memory
    double **A = allocateMatrix(n);
    double **L = allocateMatrix(n);
    double **U = allocateMatrix(n);

    // Step 3: Asking the user to input Matrix A and Matrix B elements
    RequestInput("A", A, n);

    // Step 4: Building Matrix L as Identity Matrix
    makeIdentity(L, n);

    // Step 5: Copying Matrix A to Matrix U to preserve Matrix A before decomposition
    copyMatrix(A, U, n);

    start_time = clock();
    // Step 6: Perform LU decomposition
    LU_Decomposition(A, n, L, U);
    end_time = clock();

    printf("Time taken for LU Decomposition Algorithm: %.6f seconds\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);


    // Step 7: Print results
    printSqMatrix("\nMatrix (A)", A, n);
    printSqMatrix("\nLower Triangular Matrix (L)", L, n);
    printSqMatrix("\nUpper Triangular Matrix (U)", U, n);

    // Step 8: Free memory
    freeMatrix(A, n);
    freeMatrix(L, n);
    freeMatrix(U, n);

    return 0;
}
