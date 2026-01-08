#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// allocate memory for matrices
double **allocate_matrix(int size)
{
    double **matrix = (double **)malloc(size * sizeof(double *));
    for (int i = 0; i < size; i++)
    {
        matrix[i] = (double *)malloc(size * sizeof(double));
    }
    return matrix;
}

// free allocated memory of matrices
void free_matrix(double **matrix, int size)
{
    for (int i = 0; i < size; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

// add two matrices
void add_matrix(double **A, double **B, double **C, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// subtract two matrices
void subtract_matrix(double **A, double **B, double **C, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// Function to print 2D matrix
void printSqMatrix(const char *name, double **matrix, int n)
{
    printf("%s Matrix:\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function to find the next power of 2
int next_power_of_two(int n)
{
    int power = 1;
    while (power < n)
    {
        power *= 2;
    }
    return power;
}

 //Function to pad the matrix
double **pad_matrix(double **original, int old_rows, int old_cols, int new_size)
{
    double **padded_matrix = allocate_matrix(new_size);

    for (int i = 0; i < new_size; i++)
    {
        for (int j = 0; j < new_size; j++)
        {
            if (i < old_rows && j < old_cols)
                padded_matrix[i][j] = original[i][j]; // Copy original value
            else
                padded_matrix[i][j] = 0.0; // Pad with zeros
        }
    }
    return padded_matrix;
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

// Strassen's algorithm:
void strassen_mult(double **M, double **N, double **R, int size)
{

    if (size == 1)
    {
        R[0][0] = M[0][0] * N[0][0];
        return;
    }

    int newSize = size / 2;

    // submatrices allocation
    double **a = allocate_matrix(newSize);
    double **b = allocate_matrix(newSize);
    double **c = allocate_matrix(newSize);
    double **d = allocate_matrix(newSize);
    double **x = allocate_matrix(newSize);
    double **y = allocate_matrix(newSize);
    double **z = allocate_matrix(newSize);
    double **t = allocate_matrix(newSize);
    double **q1 = allocate_matrix(newSize);
    double **q2 = allocate_matrix(newSize);
    double **q3 = allocate_matrix(newSize);
    double **q4 = allocate_matrix(newSize);
    double **q5 = allocate_matrix(newSize);
    double **q6 = allocate_matrix(newSize);
    double **q7 = allocate_matrix(newSize);
    double **r11 = allocate_matrix(newSize);
    double **r12 = allocate_matrix(newSize);
    double **r21 = allocate_matrix(newSize);
    double **r22 = allocate_matrix(newSize);

    double **temp1 = allocate_matrix(newSize);
    double **temp2 = allocate_matrix(newSize);

    // M and N submatrices (Blocks)
    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            a[i][j] = M[i][j];                     // M11
            b[i][j] = M[i][j + newSize];           // M12
            c[i][j] = M[i + newSize][j];           // M21
            d[i][j] = M[i + newSize][j + newSize]; // M22
            x[i][j] = N[i][j];                     // N11
            y[i][j] = N[i][j + newSize];           // N12
            z[i][j] = N[i + newSize][j];           // N21
            t[i][j] = N[i + newSize][j + newSize]; // N22
        }
    }

    // q1 = a * (x + z)
    add_matrix(x, z, temp2, newSize);
    strassen_mult(a, temp2, q1, newSize); // recursive call

    // q2 = d * (y + t)
    add_matrix(y, t, temp2, newSize);
    strassen_mult(d, temp2, q2, newSize); // recursive call

    // q3 = (d - a) * (z - y)
    subtract_matrix(d, a, temp1, newSize);
    subtract_matrix(z, y, temp2, newSize);
    strassen_mult(temp1, temp2, q3, newSize); // recursive call

    // q4 = (b - d) * (z + t)
    subtract_matrix(b, d, temp1, newSize);
    add_matrix(z, t, temp2, newSize);
    strassen_mult(temp1, temp2, q4, newSize); // recursive call

    // q5 = (b - a) * z
    subtract_matrix(b, a, temp1, newSize);
    strassen_mult(temp1, z, q5, newSize); // recursive call

    // q6 = (c - a) * (x + y)
    subtract_matrix(c, a, temp1, newSize);
    add_matrix(x, y, temp2, newSize);
    strassen_mult(temp1, temp2, q6, newSize); // recursive call

    // q7 = (c - d) * y
    subtract_matrix(c, d, temp1, newSize);
    strassen_mult(temp1, y, q7, newSize); // recursive call

    // r11 = q1 + q5
    add_matrix(q1, q5, r11, newSize);

    // r12 = q2 + q3 + q4 - q5
    add_matrix(q2, q3, temp1, newSize);
    add_matrix(temp1, q4, temp2, newSize);
    subtract_matrix(temp2, q5, r12, newSize);

    // r21 = q1 + q3 + q6 - q7
    add_matrix(q1, q3, temp1, newSize);
    add_matrix(temp1, q6, temp2, newSize);
    subtract_matrix(temp2, q7, r21, newSize);

    // r22 = q2 + q7
    add_matrix(q2, q7, r22, newSize);

    // result matrix R
    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            R[i][j] = r11[i][j];
            R[i][j + newSize] = r12[i][j];
            R[i + newSize][j] = r21[i][j];
            R[i + newSize][j + newSize] = r22[i][j];
        }
    }

    // free allocated memory
    free_matrix(a, newSize);
    free_matrix(b, newSize);
    free_matrix(c, newSize);
    free_matrix(d, newSize);
    free_matrix(x, newSize);
    free_matrix(y, newSize);
    free_matrix(z, newSize);
    free_matrix(t, newSize);
    free_matrix(r11, newSize);
    free_matrix(r12, newSize);
    free_matrix(r21, newSize);
    free_matrix(r22, newSize);
    free_matrix(q1, newSize);
    free_matrix(q2, newSize);
    free_matrix(q3, newSize);
    free_matrix(q4, newSize);
    free_matrix(q5, newSize);
    free_matrix(q6, newSize);
    free_matrix(q7, newSize);
    free_matrix(temp1, newSize);
    free_matrix(temp2, newSize);
}

int main()
{
    int rows_A, rows_B, cols_A, cols_B;

    // Step 1: Get dimensions for Matrix A and Matrix B
    printf("Enter dimensions for Matrix A (rows_A cols_A): ");
    scanf("%d %d", &rows_A, &cols_A);

    printf("Enter dimensions for Matrix B (cols_A cols_B): ");
    scanf("%d %d", &rows_B, &cols_B);

    // Validate dimensions
    if (rows_A <= 0 || cols_A <= 0 || cols_B <= 0)
    {
        printf("Matrix dimensions must be greater than zero. Exiting...\n");
        return 1;
    }

        if (rows_B != cols_A)
    {
        printf("Choose correct rows and columns. Exiting...\n");
        return 1;
    }

    // Step 2: Find the next power of 2 for padding
    int max_dim = rows_A > cols_A ? rows_A : cols_A; // Start with largest dimension of A
    if (cols_B > max_dim)
        max_dim = cols_B; // Compare with B's columns
    int padded_size = next_power_of_two(max_dim);

    // Step 3: Allocate original matrices with their true sizes
    double **A = allocate_matrix(rows_A);
    for (int i = 0; i < rows_A; i++)
        A[i] = (double *)malloc(cols_A * sizeof(double));

    double **B = allocate_matrix(rows_B);
    for (int i = 0; i < rows_B; i++)
        B[i] = (double *)malloc(cols_B * sizeof(double));


    // Step 4: Initialize Matrices A and B
    for (int i = 0; i < rows_A; i++)
    {
        for (int j = 0; j < cols_A; j++)
        {
            A[i][j] = ((double)rand() / RAND_MAX) * i + ((double)rand() / RAND_MAX) * j;
        }
    }

    for (int i = 0; i < rows_B; i++)
    {
        for (int j = 0; j < cols_B; j++)
        {
            
            B[i][j] = ((double)rand() / RAND_MAX) * i - ((double)rand() / RAND_MAX) * j;
        }
    }

    // Step 5: Pad matrices to the square size (padded_size x padded_size)
    double **A_padded = pad_matrix(A, rows_A, cols_A, padded_size);
    double **B_padded = pad_matrix(B, rows_B, cols_B, padded_size);

    double **R_padded = allocate_matrix(padded_size);

    // Step 6: Perform Strassen's Multiplication
    strassen_mult(A_padded, B_padded, R_padded, padded_size);

    // Step 7: Print results
    printf("\nMatrix A (Padded):\n");
    printSqMatrix("A", A_padded, padded_size);

    printf("\nMatrix B (Padded):\n");
    printSqMatrix("B", B_padded, padded_size);

    printf("\nResult Matrix (R):\n");
    printSqMatrix("R", R_padded, padded_size);

    // Step 8: Free memory
    for (int i = 0; i < rows_A; i++)
        free(A[i]);
    free(A);

    for (int i = 0; i < cols_A; i++)
        free(B[i]);
    free(B);

    free_matrix(A_padded, padded_size);
    free_matrix(B_padded, padded_size);
    free_matrix(R_padded, padded_size);

    return 0;
}
