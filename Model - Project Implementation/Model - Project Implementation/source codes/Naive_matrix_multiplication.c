/*
Group 03

Hani Abdallah - 21400302
Houssam Eddine Jamil Nasser - 21400407
Tan Viet Nguyen - 21400381

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printSqMatrix(const char *name, double **matrix, int rows, int cols)
{
    printf("%s Matrix:\n", name);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%.3f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function that request user to input Matrix elements
void RequestInput(const char *name, double **matrix, int rows, int columns)
{
    printf("Input Matrix %s Elements:\n", name);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
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

double **allocateMatrix(int rows, int columns)
{
    double **matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = (double *)malloc(columns * sizeof(double));
    }
    return matrix;
}

void Naive_matrix_multiplication(double **A, double **B, double **R, int rows_A, int cols_A, int cols_B)
{
    for (int i = 0; i < rows_A; i++)
    {
        for (int j = 0; j < cols_B; j++)
        {
            R[i][j] = 0.0;
            for (int k = 0; k < cols_A; k++)
            {
                R[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main()
{
    clock_t start_time, end_time;
    int rows_A;
    int cols_A;
    int rows_B;
    int cols_B;

    // Step 1: Get dimensions for Matrix A and Matrix B
    printf("\nChoose Matrix A Rows: ");
    scanf("%d", &rows_A);
    printf("\nChoose Matrix A columns: ");
    scanf("%d", &cols_A);
    printf("\nChoose Matrix B Rows: ");
    scanf("%d", &rows_B);
    printf("\nChoose Matrix B columns: ");
    scanf("%d", &cols_B);

    int rows_C = rows_A;
    int cols_C = cols_B;

    // Validate dimensions
    // Error condition 1
    if (rows_A <= 0 || cols_A <= 0 || rows_B <= 0 || cols_B <= 0)
    {
        printf("Matrix dimensions can only be larger than zero, the program will exit...");
        exit(1);
    }

    // Error condition 2
    if (cols_A != rows_B)
    {
        printf("The columns of Matrix A should be equal to rows of Matrix B, the program will exit...");
        exit(1);
    }

    // Step 2: Allocate matrices memory
    double **A = allocateMatrix(rows_A, cols_A);
    double **B = allocateMatrix(rows_B, cols_B);
    double **C = allocateMatrix(rows_C, cols_C);

    // Note: You can either enter the matrix input manually, or choose to fill the matrices automatically. But make sure to comment "RequestInput" function calls, and to uncomment the alternative random input option.
    // Step 3: Asking the user to input Matrix A and Matrix B elements
    RequestInput("A", A, rows_A, cols_A);
    RequestInput("B", B, rows_B, cols_B);

    // Step 3: Initialize Matrices A and B - **Alternative random input**
    /*for (int i = 0; i < rows_A; i++)
    {
        for (int j = 0; j < cols_A; j++)
        {
            A[i][j] = (-0.33) * i + (0.051) * j;
        }
    }

    for (int i = 0; i < rows_B; i++)
    {
        for (int j = 0; j < cols_B; j++)
        {
            B[i][j] = (0.036) * i - (1.9047) * j;
        }
    }*/
    start_time = clock();
    // Step 4: Perform Strassen's Multiplication
    Naive_matrix_multiplication(A, B, C, rows_A, cols_A, cols_B);
    end_time = clock();

    // Step 5: Print results
    printSqMatrix("\nA", A, rows_A, cols_A);
    printSqMatrix("\nB", B, rows_B, cols_B);
    printSqMatrix("\nC", C, rows_C, cols_C);

    printf("\nTime taken for Naive multiplication Algorithm: %.6f seconds\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

    // Step 6: Free memory
    freeMatrix(A, rows_A);
    freeMatrix(B, rows_B);
    freeMatrix(C, rows_C);

    return 0;
}
