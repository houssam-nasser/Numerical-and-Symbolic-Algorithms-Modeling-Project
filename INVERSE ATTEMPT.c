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

// Function to print 2D matrix
void printSqMatrix(const char *name, double **matrix, int n)
{
    printf("%s Matrix:\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.4f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

// subtraction of two matrices
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

// addition of two matrices
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
 double **pad_matrix(double **original, int old_size, int new_size)
{
    double **padded_matrix = allocate_matrix(new_size);
    for (int i = 0; i < new_size; i++)
    {
        for (int j = 0; j < new_size; j++)
        {
            if (i < old_size && j < old_size)
                padded_matrix[i][j] = original[i][j]; // Copy original value
            else
                padded_matrix[i][j] = (i == j) ? 1.0 : 0.01; // Pad with non zeros
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
    double **a = allocate_matrix(newSize); // M11
    double **b = allocate_matrix(newSize); // M12
    double **c = allocate_matrix(newSize); // M21
    double **d = allocate_matrix(newSize); // M22
    double **x = allocate_matrix(newSize); // N11
    double **y = allocate_matrix(newSize); // N12
    double **z = allocate_matrix(newSize); // N21
    double **t = allocate_matrix(newSize); // N22
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
            a[i][j] = M[i][j];
            b[i][j] = M[i][j + newSize];
            c[i][j] = M[i + newSize][j];
            d[i][j] = M[i + newSize][j + newSize];
            x[i][j] = N[i][j];
            y[i][j] = N[i][j + newSize];
            z[i][j] = N[i + newSize][j];
            t[i][j] = N[i + newSize][j + newSize];
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

// Strassen's Matrix Inversion
void strassen_inversion(double **A, double **A_inv, int size)
{
    if (size == 1)
    {
        if (A[0][0] == 0)
        {
            fprintf(stderr, "Matrix is singular and cannot be inverted.\n");
            exit(1);
        }
        A_inv[0][0] = 1.0 / A[0][0];
        return;
    }

    int newSize = size / 2;

    double **a = allocate_matrix(newSize);
    double **b = allocate_matrix(newSize);
    double **c = allocate_matrix(newSize);
    double **d = allocate_matrix(newSize);
    double **e = allocate_matrix(newSize);
    double **z = allocate_matrix(newSize);
    double **t = allocate_matrix(newSize);
    double **x = allocate_matrix(newSize);
    double **y = allocate_matrix(newSize);

    double **temp1 = allocate_matrix(newSize);
    double **temp2 = allocate_matrix(newSize);

    // Split A into submatrices
    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            a[i][j] = A[i][j];                     // A11
            b[i][j] = A[i][j + newSize];           // A12
            c[i][j] = A[i + newSize][j];           // A21
            d[i][j] = A[i + newSize][j + newSize]; // A22
        }
    }

    // e = a^-1
    strassen_inversion(a, e, newSize); // recursive call

    // z = d - c * e * b
    strassen_mult(e, b, temp1, newSize);
    strassen_mult(c, temp1, temp2, newSize);
    subtract_matrix(d, temp2, z, newSize);

    // t = z^-1
    strassen_inversion(z, t, newSize); // recursive call

    // y = -e * b * t
    strassen_mult(b, t, temp1, newSize);
    strassen_mult(e, temp1, y, newSize);
    for (int i = 0; i < newSize; i++)
    { // y = -(e * b * t)
        for (int j = 0; j < newSize; j++)
        {
            y[i][j] = -y[i][j];
        }
    }

    // z = -t * c * e
    strassen_mult(c, e, temp1, newSize);
    strassen_mult(t, temp1, z, newSize);
    for (int i = 0; i < newSize; i++)
    { // z = -(t * c * e)
        for (int j = 0; j < newSize; j++)
        {
            z[i][j] = -z[i][j];
        }
    }

    // x = e + e * b * t * c * e
    strassen_mult(b, t, temp1, newSize);     // temp1 = b * t
    strassen_mult(temp1, c, temp2, newSize); // temp2 = (b * t) * c
    strassen_mult(e, temp2, temp1, newSize); // temp1 = e * ((b * t) * c)
    strassen_mult(temp1, e, temp2, newSize); // temp2 = (e * ((b * t) * c)) * e
    add_matrix(e, temp2, x, newSize);        // x = e + (e * (b * t * c))

    // Combine x, y, z, t, into A_inv
    for (int i = 0; i < newSize; i++)
    {
        for (int j = 0; j < newSize; j++)
        {
            A_inv[i][j] = x[i][j];                     // A_inv11
            A_inv[i][j + newSize] = y[i][j];           // A_inv12
            A_inv[i + newSize][j] = z[i][j];           // A_inv21
            A_inv[i + newSize][j + newSize] = t[i][j]; // A_inv22
        }
    }

    // Free allocated memory
    free_matrix(a, newSize);
    free_matrix(b, newSize);
    free_matrix(c, newSize);
    free_matrix(d, newSize);
    free_matrix(e, newSize);
    free_matrix(z, newSize);
    free_matrix(t, newSize);
    free_matrix(x, newSize);
    free_matrix(y, newSize);
    free_matrix(temp1, newSize);
    free_matrix(temp2, newSize);
}

int main()
{
    int n;

    printf("\nChoose Matrix Dimension: ");
    scanf("%d", &n);

    // Error condition (matrix size must be > 0 and a power of 2)

    if (n <= 0)
    {
        printf("Matrix Size can only be larger than zero. The program will exit...");
        exit(1);
    }

    int padded_size = next_power_of_two(n);

    // Allocate matrices
    double **A = allocate_matrix(n);

    // Request to input Matrices entries from user
    //RequestInput("A", A, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = (i == j) ? 1.0 : ((double)rand() / RAND_MAX);

        }
    }

     double **A_padded = pad_matrix(A, n, padded_size);
     
    double **A_inv = allocate_matrix(padded_size);

    // Perform Strassen's Inversion ALgo
    strassen_inversion(A, A_inv, padded_size);

    printSqMatrix("\nA", A, padded_size);
    printSqMatrix("\nInverse of A", A_inv, padded_size);

    free_matrix(A, padded_size);
    free_matrix(A_inv, padded_size);

    return 0;
}