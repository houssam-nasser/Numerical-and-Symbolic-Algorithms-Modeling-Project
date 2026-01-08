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
            a[i][j] = M[i][j];                       // M11    
            b[i][j] = M[i][j + newSize];             // M12                
            c[i][j] = M[i + newSize][j];             // M21                
            d[i][j] = M[i + newSize][j + newSize];   // M22                        
            x[i][j] = N[i][j];                       // N11    
            y[i][j] = N[i][j + newSize];             // N12                
            z[i][j] = N[i + newSize][j];             // N21                
            t[i][j] = N[i + newSize][j + newSize];   // N22                        
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

int main()   // QUICK TEST . GO DOWN FOR THE MANUAL TEST
{
    int size = 10; //test size (must be a power of 2)

    //allocate M and N
    double **M = allocate_matrix(size);
    double **N = allocate_matrix(size);
    double **R = allocate_matrix(size);

    //test initialization
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            M[i][j] = ((double)rand() / RAND_MAX) * i + ((double)rand() / RAND_MAX) * j;
            N[i][j] = ((double)rand() / RAND_MAX) * i - ((double)rand() / RAND_MAX) * j;
        }
    }

    // Perform Strassen's Algo
    strassen_mult(M, N, R, size);

    //result
    printf("R:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%f ", R[i][j]);
        }
        printf("\n");
    }

    free_matrix(M, size);
    free_matrix(N, size);
    free_matrix(R, size);

    return 0;
}

/* int main()      // MANUAL TEST 
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

    //allocate M and N
    double **M = allocate_matrix(size);
    double **N = allocate_matrix(size);
    double **R = allocate_matrix(size);

    //matrix M
    printf("Enter elements of matrix M:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("M[%d][%d]=", i,j);
            scanf("%lf", &M[i][j]);
        }
    }

    //matrix N
    printf("Enter elements of matrix N:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("N[%d][%d]=", i,j);
            scanf("%lf", &N[i][j]);
        }
    }

    // Perform Strassen's Algo
    strassen_mult(M, N, R, size);

    //result
    printf("Matrix R:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%f ", R[i][j]); 
        }
        printf("\n");
    }

    free_matrix(M, size);
    free_matrix(N, size);
    free_matrix(R, size);

    return EXIT_SUCCESS;
} */
