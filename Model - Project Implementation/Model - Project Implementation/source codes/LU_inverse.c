/*
Group 03

Hani Abdallah - 21400302
Houssam Eddine Jamil Nasser - 21400407
Tan Viet Nguyen - 21400381

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// Function to allocate memory for a square matrix
double **allocateMatrix(int n)
{
    double **matrix = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        matrix[i] = malloc(n * sizeof(double));
    }
    return matrix;
}

// Free allocated matrix memory
void freeMatrix(double **matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

// Function to print a matrix
void printMatrix(const char *name, double **matrix, int n)
{
    printf("%s:\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%10.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// LU Decomposition
void LU_Decomposition(double **A, int n, double **L, double **U)
{
    for (int i = 0; i < n; i++)
    {
        // Diagonal of L
        L[i][i] = 1.0;

        // Upper triangular U
        for (int j = i; j < n; j++)
        {
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++)
            {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }

        // Lower triangular L
        for (int j = i + 1; j < n; j++)
        {
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++)
            {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
}

// Forward substitution
void forwardSubstitution(double **L, double *b, double *y, int n)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = b[i];
        for (int j = 0; j < i; j++)
        {
            y[i] -= L[i][j] * y[j];
        }
    }
}

// Backward substitution
void backwardSubstitution(double **U, double *y, double *x, int n)
{
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++)
        {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

// Invert matrix using LU decomposition
void invertMatrix(double **A, double **A_inverse, int n)
{
    double **L = allocateMatrix(n);
    double **U = allocateMatrix(n);
    LU_Decomposition(A, n, L, U);

    double *b = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));
    double *x = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        // Set up identity matrix column
        for (int j = 0; j < n; j++)
        {
            b[j] = (i == j) ? 1.0 : 0.0;
        }

        // Solve Ly = b
        forwardSubstitution(L, b, y, n);

        // Solve Ux = y
        backwardSubstitution(U, y, x, n);

        // Store solution in the inverse matrix
        for (int j = 0; j < n; j++)
        {
            A_inverse[j][i] = x[j];
        }
    }

    freeMatrix(L, n);
    freeMatrix(U, n);
    free(b);
    free(y);
    free(x);
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

int main()
{

    clock_t start_time, end_time;
    double cpu_time;

    int n;
    // Step 1: Get dimensions for Matrix A
    printf("\nChoose Matrix Dimension for the square matrix: ");
    scanf("%d", &n);

    // Error condition 1
    if (n <= 0)
    {
        printf("Matrix Size can only be larger than zero, the program will exit...");
        exit(1);
    }

    // Step 2: Allocate matrices memory
    double **A = allocateMatrix(n);
    double **A_inverse = allocateMatrix(n);

    // Note: You can either enter the matrix input manually, or choose to fill the matrices automatically. But make sure to comment "RequestInput" function calls, and to uncomment the alternative random input option.
    // Step 3: Asking the user to input Matrix A and Matrix B elements
    RequestInput("A", A, n);

    // This is an alternative random input option, but make sure to comment "RequestInput" function calls.
    /*for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = (i == j) ? 1.0 : ((double)rand() / RAND_MAX);
        }
    }*/

    start_time = clock();
    // Step 4: Perform LU inversion
    invertMatrix(A, A_inverse, n);
    end_time = clock();

    // Step 5: Print results
    printMatrix("Matrix A", A, n);
    printMatrix("Matrix A Inverse", A_inverse, n);

    printf("Time taken for LU Inversion Algorithm: %.6f seconds\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

    // Step 6: Free memory
    freeMatrix(A, n);
    freeMatrix(A_inverse, n);

    return 0;
}
