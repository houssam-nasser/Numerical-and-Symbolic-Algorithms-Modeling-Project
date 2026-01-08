#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS 1
#define FAILURE 0
enum ERROR_CODES
{
    OK,
    INCORRECT_MATRIX,
    CALCULATION_ERROR
};

typedef struct matrix_struct
{
    double **matrix;
    int rows;
    int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result)
{
    if (result == NULL)
    {
        return INCORRECT_MATRIX;
    }

    int matrix_status = OK;

    if (rows < 1 || columns < 1)
    {
        matrix_status = INCORRECT_MATRIX;
        result->matrix = NULL;
    }
    else
    {
        result->rows = rows;
        result->columns = columns;
        result->matrix = (double **)calloc(rows + rows * columns, sizeof(double));
        if (result->matrix == NULL)
        {
            matrix_status = INCORRECT_MATRIX;
        }
        else
        {
            // allocating one block of memory for everything at once
            double *start = (double *)(result->matrix + rows);
            // indexing our matrix
            int iteration_limit = rows <= columns ? columns : rows;
            for (int i = 0; i < iteration_limit; i++)
            {
                result->matrix[i] = start + i * columns;
            }
        }
    }

    return matrix_status;
}

int s21_validate_matrix(int matrix_amount, matrix_t *A, ...)
{
    if (A == NULL || A->matrix == NULL || A->rows < 1 || A->columns < 1)
    {
        return INCORRECT_MATRIX;
    }

    int matrices_status = OK;

    va_list matrix_list;
    va_start(matrix_list, A);
    for (int i = 0; i < matrix_amount - 1 && matrices_status == OK; i++)
    {
        matrix_t *current_matrix = va_arg(matrix_list, matrix_t *);
        if (current_matrix == NULL || current_matrix->matrix == NULL ||
            current_matrix->rows < 1 || current_matrix->columns < 1)
        {
            matrices_status = INCORRECT_MATRIX;
        }
    }
    va_end(matrix_list);

    return matrices_status;
}

double s21_mult_matrix_res(int i, int j, matrix_t *A, matrix_t *B)
{
    if (A == NULL || B == NULL)
    {
        return 0;
    }

    double res = 0;

    // Burroughs reference
    for (int k = 0; k < B->rows; k++)
    {
        res += A->matrix[i][k] * B->matrix[k][j];
    }

    return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result)
{
    if (A == NULL || A->matrix == NULL || result == NULL)
    {
        return INCORRECT_MATRIX;
    }
    else if (!isfinite(number))
    {
        return CALCULATION_ERROR;
    }

    int error_code = OK;
    error_code = s21_create_matrix(A->rows, A->columns, result);
    if (!error_code)
    {
        for (int i = 0; i < A->rows && error_code == OK; i++)
        {
            for (int j = 0; j < A->columns && error_code == OK; j++)
            {
                result->matrix[i][j] = A->matrix[i][j] * number;
                if (!isfinite(result->matrix[i][j]))
                {
                    error_code = CALCULATION_ERROR;
                }
            }
        }
    }

    return error_code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result)
{
    if (s21_validate_matrix(2, A, B) || result == NULL)
    {
        return INCORRECT_MATRIX;
    }
    else if (A->columns != B->rows)
    {
        return CALCULATION_ERROR;
    }

    int error_code = OK;
    error_code = s21_create_matrix(A->rows, B->columns, result);
    if (!error_code)
    {
        for (int i = 0; i < result->rows && error_code == OK; i++)
        {
            for (int j = 0; j < result->columns && error_code == OK; j++)
            {
                result->matrix[i][j] = s21_mult_matrix_res(i, j, A, B);
                if (!isfinite(result->matrix[i][j]))
                {
                    error_code = CALCULATION_ERROR;
                }
            }
        }
    }

    return error_code;
}

int s21_is_matrix_same_size(int matrix_amount, matrix_t *A, ...)
{
    if (A == NULL)
    {
        return FAILURE;
    }

    int dimensions_comp_status = SUCCESS;
    int rows_to_compare = A->rows;
    int columns_to_compare = A->columns;

    va_list matrix_list;
    va_start(matrix_list, A);
    for (int i = 0; i < matrix_amount - 1 && dimensions_comp_status == SUCCESS;
         i++)
    {
        matrix_t *current_matrix = va_arg(matrix_list, matrix_t *);
        if (current_matrix == NULL || current_matrix->rows != rows_to_compare ||
            current_matrix->columns != columns_to_compare)
        {
            dimensions_comp_status = FAILURE;
        }
    }
    va_end(matrix_list);

    return dimensions_comp_status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result)
{
    if (s21_validate_matrix(2, A, B) || result == NULL)
    {
        return INCORRECT_MATRIX;
    }
    else if (!s21_is_matrix_same_size(2, A, B))
    {
        return CALCULATION_ERROR;
    }

    int error_code = OK;
    error_code = s21_create_matrix(A->rows, A->columns, result);
    if (!error_code)
    {
        for (int i = 0; i < A->rows && error_code == OK; i++)
        {
            for (int j = 0; j < A->columns && error_code == OK; j++)
            {
                result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
                if (!isfinite(result->matrix[i][j]))
                {
                    error_code = CALCULATION_ERROR;
                }
            }
        }
    }

    return error_code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result)
{
    if (s21_validate_matrix(2, A, B) || result == NULL)
    {
        return INCORRECT_MATRIX;
    }
    else if (!s21_is_matrix_same_size(2, A, B))
    {
        return CALCULATION_ERROR;
    }

    int error_code = OK;
    error_code = s21_create_matrix(A->rows, A->columns, result);
    if (!error_code)
    {
        for (int i = 0; i < A->rows && error_code == OK; i++)
        {
            for (int j = 0; j < A->columns && error_code == OK; j++)
            {
                result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
                if (!isfinite(result->matrix[i][j]))
                {
                    error_code = CALCULATION_ERROR;
                }
            }
        }
    }

    return error_code;
}

void get_identity_matrix(int n, matrix_t *result)
{
    s21_create_matrix(n, n, result);
    for (int i = 0; i < n; i++)
    {
        result->matrix[i][i] = 1;
    }
}

void s21_print_matrix(matrix_t *A)
{
    if (A != NULL && A->matrix != NULL)
    {
        for (int i = 0; i < A->rows; i++)
        {
            for (int j = 0; j < A->columns; j++)
            {
                printf("%f\t", A->matrix[i][j]);
            }
            printf("\n");
        }
    }
}

void merge_matrices(matrix_t *first_matrix, matrix_t *second_matrix,
                    matrix_t *third_matrix, matrix_t *fourth_matrix, matrix_t *filled_matrix)
{
    // int columns = first_matrix->columns+second_matrix->columns;
    // int rows = first_matrix->rows+third_matrix->rows;

    for (int i = 0; i < first_matrix->rows; i++)
    {
        for (int j = 0; j < first_matrix->columns; j++)
        {
            filled_matrix->matrix[i][j] = first_matrix->matrix[i][j];
        }
    }

    for (int i = 0; i < second_matrix->rows; i++)
    {
        for (int j = 0; j < second_matrix->columns; j++)
        {
            filled_matrix->matrix[i][first_matrix->columns + j] = second_matrix->matrix[i][j];
        }
    }

    for (int i = 0; i < third_matrix->rows; i++)
    {
        for (int j = 0; j < third_matrix->columns; j++)
        {
            filled_matrix->matrix[third_matrix->rows + i][j] = third_matrix->matrix[i][j];
        }
    }

    for (int i = 0; i < fourth_matrix->rows; i++)
    {
        for (int j = 0; j < fourth_matrix->columns; j++)
        {
            filled_matrix->matrix[fourth_matrix->rows + i][fourth_matrix->columns + j] = fourth_matrix->matrix[i][j];
        }
    }
}

void split_matrix_into_quadrants(matrix_t *A,
                                 matrix_t *top_left, matrix_t *top_right,
                                 matrix_t *bottom_left, matrix_t *bottom_right)
{
    int size = A->rows;
    int split_point = size / 2;

    for (int i = 0; i < split_point; i++)
    {
        for (int j = 0; j < split_point; j++)
        {
            top_left->matrix[i][j] = A->matrix[i][j];
        }
    }

    for (int i = split_point; i < size; i++)
    {
        for (int j = 0; j < split_point; j++)
        {
            bottom_left->matrix[i - split_point][j] = A->matrix[i][j];
        }
    }

    for (int i = 0; i < split_point; i++)
    {
        for (int j = split_point; j < size; j++)
        {
            top_right->matrix[i][j - split_point] = A->matrix[i][j];
        }
    }

    for (int i = split_point; i < size; i++)
    {
        for (int j = split_point; j < size; j++)
        {
            bottom_right->matrix[i - split_point][j - split_point] = A->matrix[i][j];
        }
    }
}

void inverse(matrix_t *A, matrix_t *inv_A)
{
    if (A->columns == 1 && A->rows == 1)
    {
        if (A->matrix[0][0] == 0)
        {
            printf("Cannot devide by 0!!\n");
        }
        else
        {
            inv_A->matrix[0][0] = 1 / A->matrix[0][0];
        }
    }
    else if (A->rows == 2 && A->columns == 2)
    {
        matrix_t A_11 = {};
        s21_create_matrix(1, 1, &A_11);
        A_11.matrix[0][0] = A->matrix[0][0];

        matrix_t A_12 = {};
        s21_create_matrix(1, 1, &A_12);
        A_12.matrix[0][0] = A->matrix[0][1];

        matrix_t A_21 = {};
        s21_create_matrix(1, 1, &A_21);
        A_21.matrix[0][0] = A->matrix[1][0];

        matrix_t A_22 = {};
        s21_create_matrix(1, 1, &A_22);
        A_22.matrix[0][0] = A->matrix[1][1];

        matrix_t A_11_inv = {};
        s21_create_matrix(1, 1, &A_11_inv);
        inverse(&A_11, &A_11_inv);

        matrix_t A_21_A_11_inv = {};
        s21_create_matrix(1, 1, &A_21_A_11_inv);

        s21_mult_matrix(&A_21, &A_11_inv, &A_21_A_11_inv);

        matrix_t A_21_A_11_inv_A_12 = {};
        s21_create_matrix(1, 1, &A_21_A_11_inv_A_12);

        s21_mult_matrix(&A_21_A_11_inv, &A_12, &A_21_A_11_inv_A_12);

        matrix_t S_22 = {};
        s21_create_matrix(1, 1, &S_22);

        s21_sub_matrix(&A_22, &A_21_A_11_inv_A_12, &S_22);

        matrix_t S_22_inv = {};
        s21_create_matrix(1, 1, &S_22_inv);
        inverse(&S_22, &S_22_inv);

        matrix_t A_12_S_22_inv = {};
        s21_create_matrix(1, 1, &A_12_S_22_inv);
        s21_mult_matrix(&A_12, &S_22_inv, &A_12_S_22_inv);

        matrix_t A_12_S_22_inv_A_21 = {};
        s21_create_matrix(1, 1, &A_12_S_22_inv_A_21);

        s21_mult_matrix(&A_12_S_22_inv, &A_21, &A_12_S_22_inv_A_21);

        matrix_t A_12_S_22_inv_A_21_A_11_inv = {};
        s21_create_matrix(1, 1, &A_12_S_22_inv_A_21_A_11_inv);

        s21_mult_matrix(&A_12_S_22_inv_A_21, &A_11_inv, &A_12_S_22_inv_A_21_A_11_inv);

        A_12_S_22_inv_A_21_A_11_inv.matrix[0][0] += 1;

        matrix_t B_11 = {};
        s21_create_matrix(1, 1, &B_11);
        s21_mult_matrix(&A_11_inv, &A_12_S_22_inv_A_21_A_11_inv, &B_11);

        matrix_t A_11_inv_A_12 = {};
        s21_create_matrix(1, 1, &A_11_inv_A_12);
        s21_mult_matrix(&A_11_inv, &A_12, &A_11_inv_A_12);
        matrix_t B_12 = {};
        s21_create_matrix(1, 1, &B_12);
        s21_mult_matrix(&A_11_inv_A_12, &S_22_inv, &B_12);

        B_12.matrix[0][0] *= -1;

        matrix_t S_22_inv_A_21 = {};
        s21_create_matrix(1, 1, &S_22_inv_A_21);
        s21_mult_matrix(&S_22_inv, &A_21, &S_22_inv_A_21);

        matrix_t B_21 = {};
        s21_create_matrix(1, 1, &B_21);
        s21_mult_matrix(&S_22_inv_A_21, &A_11_inv, &B_21);
        B_21.matrix[0][0] *= -1;

        matrix_t B_22 = {};
        s21_create_matrix(1, 1, &B_22);

        B_22.matrix[0][0] = S_22_inv.matrix[0][0];

        merge_matrices(&B_11, &B_12, &B_21, &B_22, inv_A);
    }

    else
    {
        int split_point = A->rows / 2;
        matrix_t A_11 = {};
        s21_create_matrix(split_point, split_point, &A_11);

        matrix_t A_12 = {};
        s21_create_matrix(split_point, split_point, &A_12);

        matrix_t A_21 = {};
        s21_create_matrix(split_point, split_point, &A_21);

        matrix_t A_22 = {};
        s21_create_matrix(split_point, split_point, &A_22);

        split_matrix_into_quadrants(A, &A_11, &A_12, &A_21, &A_22);

        matrix_t id = {};
        s21_create_matrix(A_11.rows, A_11.rows, &id);
        get_identity_matrix(A_11.rows, &id);

        matrix_t A_11_inv = {};
        s21_create_matrix(split_point, split_point, &A_11_inv);

        inverse(&A_11, &A_11_inv);

        matrix_t A_21_A_11_inv = {};
        s21_create_matrix(split_point, split_point, &A_21_A_11_inv);

        s21_mult_matrix(&A_21, &A_11_inv, &A_21_A_11_inv);

        matrix_t A_21_A_11_inv_A_12 = {};
        s21_create_matrix(split_point, split_point, &A_21_A_11_inv_A_12);

        s21_mult_matrix(&A_21_A_11_inv, &A_12, &A_21_A_11_inv_A_12);

        matrix_t S_22 = {};
        s21_create_matrix(split_point, split_point, &S_22);

        s21_sub_matrix(&A_22, &A_21_A_11_inv_A_12, &S_22);

        matrix_t S_22_inv = {};
        s21_create_matrix(split_point, split_point, &S_22_inv);
        inverse(&S_22, &S_22_inv);

        matrix_t A_12_S_22_inv = {};
        s21_create_matrix(split_point, split_point, &A_12_S_22_inv);
        s21_mult_matrix(&A_12, &S_22_inv, &A_12_S_22_inv);

        matrix_t A_12_S_22_inv_A_21 = {};
        s21_create_matrix(split_point, split_point, &A_12_S_22_inv_A_21);

        s21_mult_matrix(&A_12_S_22_inv, &A_21, &A_12_S_22_inv_A_21);

        matrix_t A_12_S_22_inv_A_21_A_11_inv = {};
        s21_create_matrix(split_point, split_point, &A_12_S_22_inv_A_21_A_11_inv);

        s21_mult_matrix(&A_12_S_22_inv_A_21, &A_11_inv, &A_12_S_22_inv_A_21_A_11_inv);

        matrix_t A_12_S_22_inv_A_21_A_11_inv_id = {};
        s21_create_matrix(split_point, split_point, &A_12_S_22_inv_A_21_A_11_inv_id);

        s21_sum_matrix(&A_12_S_22_inv_A_21_A_11_inv, &id, &A_12_S_22_inv_A_21_A_11_inv_id);

        matrix_t B_11 = {};
        s21_create_matrix(split_point, split_point, &B_11);
        s21_mult_matrix(&A_11_inv, &A_12_S_22_inv_A_21_A_11_inv_id, &B_11);

        matrix_t A_11_inv_A_12 = {};
        s21_create_matrix(split_point, split_point, &A_11_inv_A_12);
        s21_mult_matrix(&A_11_inv, &A_12, &A_11_inv_A_12);
        matrix_t B_12d = {};
        s21_create_matrix(split_point, split_point, &B_12d);
        s21_mult_matrix(&A_11_inv_A_12, &S_22_inv, &B_12d);

        matrix_t B_12 = {};
        s21_create_matrix(split_point, split_point, &B_12);
        s21_mult_number(&B_12d, -1, &B_12);

        matrix_t S_22_inv_A_21 = {};
        s21_create_matrix(split_point, split_point, &S_22_inv_A_21);
        s21_mult_matrix(&S_22_inv, &A_21, &S_22_inv_A_21);

        matrix_t B_21d = {};
        s21_create_matrix(split_point, split_point, &B_21d);
        s21_mult_matrix(&S_22_inv_A_21, &A_11_inv, &B_21d);

        matrix_t B_21 = {};
        s21_create_matrix(split_point, split_point, &B_21);
        s21_mult_number(&B_21d, -1, &B_21);
        // B_21.matrix[0][0] *= -1;

        matrix_t B_22 = {};
        s21_create_matrix(split_point, split_point, &B_22);

        for (int i = 0; i < split_point; i++)
        {
            for (int j = 0; j < split_point; j++)
            {
                B_22.matrix[i][j] = S_22_inv.matrix[i][j];
            }
        }
        merge_matrices(&B_11, &B_12, &B_21, &B_22, inv_A);
    }
}

// You can use the below to make dynamic input!!

// Function that request user to input Matrix elements
void RequestInput(const char *name, matrix_t *A, int n)
{
    printf("Input Matrix %s Elements:\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%s[%d][%d]=", name, i, j);
            scanf("%lf", &A->matrix[i][j]);
        }
    }
}

// free allocated memory of matrices
void freeMatrix(matrix_t *A, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(A->matrix[i]);
    }
    free(A);
}

int main()
{

    int n;

    // Asking the user to input the matrix dimensions
    printf("\nChoose Matrix Dimension: ");
    scanf("%d", &n);

    // Error condition
    if (!is_power_of_two(n) || n <= 0)
    {
        printf("Matrix Size can only be larger than zero, and a power of 2. The program will exit...");
        exit(1);
    }

    matrix_t A = {};

    // Asking the user to input Matrix A elements
    // RequestInput("A", A, n);

    /*
    A.matrix[0][0] = 3;
    A.matrix[1][0] = 5;
    A.matrix[2][0] = 6;
    A.matrix[3][0] = 5;
    A.matrix[4][0] = 4;
    A.matrix[5][0] = 3;
    A.matrix[6][0] = 2;
    A.matrix[7][0] = 1;

    A.matrix[0][1] = 3;
    A.matrix[1][1] = 4;
    A.matrix[2][1] = 5;
    A.matrix[3][1] = 6;
    A.matrix[4][1] = 7;
    A.matrix[5][1] = 8;
    A.matrix[6][1] = 9;
    A.matrix[7][1] = 43;

    A.matrix[0][2] = 2;
    A.matrix[1][2] = 3;
    A.matrix[2][2] = 4;
    A.matrix[3][2] = 5;
    A.matrix[4][2] = 6;
    A.matrix[5][2] = 7;
    A.matrix[6][2] = 4;
    A.matrix[7][2] = 9;

    A.matrix[0][3] = 2;
    A.matrix[1][3] = 31;
    A.matrix[2][3] = 2;
    A.matrix[3][3] = 3;
    A.matrix[4][3] = 46;
    A.matrix[5][3] = 7;
    A.matrix[6][3] = 4;
    A.matrix[7][3] = 49;

    A.matrix[0][4] = 3;
    A.matrix[1][4] = 2;
    A.matrix[2][4] = 3;
    A.matrix[3][4] = 5;
    A.matrix[4][4] = 6;
    A.matrix[5][4] = 8;
    A.matrix[6][4] = 49;
    A.matrix[7][4] = 7;

    A.matrix[0][5] = 2;
    A.matrix[1][5] = 3;
    A.matrix[2][5] = 1;
    A.matrix[3][5] = 3;
    A.matrix[4][5] = 2;
    A.matrix[5][5] = 3;
    A.matrix[6][5] = 4;
    A.matrix[7][5] = 5;

    A.matrix[0][6] = 2;
    A.matrix[1][6] = 3;
    A.matrix[2][6] = 1;
    A.matrix[3][6] = 3;
    A.matrix[4][6] = 2;
    A.matrix[5][6] = 4;
    A.matrix[6][6] = 5;
    A.matrix[7][6] = 6;

    A.matrix[0][7] = 2;
    A.matrix[1][7] = 3;
    A.matrix[2][7] = 1;
    A.matrix[3][7] = 3;
    A.matrix[4][7] = 1;
    A.matrix[5][7] = 2;
    A.matrix[6][7] = 3;
    A.matrix[7][7] = 6;
    */

    s21_create_matrix(2, 2, &A);

    printf("Matrix A:\n");
    s21_print_matrix(&A);

    matrix_t inv_A = {};
    s21_create_matrix(2, 2, &inv_A);

    inverse(&A, &inv_A);

    printf("\n Matrix inverse (A^-1): ");
    s21_print_matrix(&inv_A);

    printf("\n");

    matrix_t A_inv_A = {};

    s21_create_matrix(2, 2, &A_inv_A);

    s21_mult_matrix(&A, &inv_A, &A_inv_A);

    // Added here
    // freeMatrix(A, n);
    // freeMatrix(inv_A, n);

    return 0;
}
