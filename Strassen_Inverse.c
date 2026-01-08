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

// Invert a 1x1 matrix
void invert_1x1(double **A, double **A_inv)
{
    if (A[0][0] == 0)
    {
        fprintf(stderr, "Matrix is singular and cannot be inverted.\n");
        exit(EXIT_FAILURE);
    }
    A_inv[0][0] = 1.0 / A[0][0];
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

bool is_power_of_two(int n)
{
    return (n > 0) && ((n & (n - 1)) == 0);
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
        invert_1x1(A, A_inv);
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

int main() // QUICK TEST . GO DOWN FOR THE MANUAL TEST
{
    int size = 2; // test size (must be a power of 2)

    // Allocate matrices
    double **A = allocate_matrix(size);
    double **A_inv = allocate_matrix(size);
    double **C = allocate_matrix(size);

    /// test initialization
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i][j] = 1.34 + i + j;
        }
    }

    // Perform Strassen's Inversion ALgo
    strassen_inversion(A, A_inv, size);



    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                C[i][j] += A[i][k] * A_inv[k][j];
            }
        }
    }

    // result
    printf("Inverse of C:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%f ", C[i][j]);
        }
        printf("\n");
    }

    free_matrix(A, size);
    free_matrix(A_inv, size);
    free_matrix(C, size);

    return 0;
}

/*
int main() // MANUAL TEST
{
    int size;

    printf("Enter size of matrix (must be a power of 2): ");
    scanf("%d", &size);

    // Validate the size
    if (size <= 0 || !is_power_of_two(size))
    {
        printf("Error: Size must be a positive power of 2.\n");
        return EXIT_FAILURE;
    }

    double **A = allocate_matrix(size);
    double **A_inv = allocate_matrix(size);

    // Input matrix M
    printf("Enter elements of A:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("A[%d][%d] = ", i, j);
            scanf("%lf", &A[i][j]);
        }
    }

    // Perform Strassen's Inversion
    strassen_inversion(A, A_inv, size);

    // Result
    printf("Inverse of A:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%f ", A_inv[i][j]);
        }
        printf("\n");
    }

    // Free allocated memory
    free_matrix(A, size);
    free_matrix(A_inv, size);

    return EXIT_SUCCESS;
}
*/
