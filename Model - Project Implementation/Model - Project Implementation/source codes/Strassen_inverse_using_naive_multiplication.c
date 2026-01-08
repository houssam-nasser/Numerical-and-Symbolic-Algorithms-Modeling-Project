/*
Group 03

Hani Abdallah - 21400302
Houssam Eddine Jamil Nasser - 21400407
Tan Viet Nguyen - 21400381

*/
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

int create_matrix(int rows, int columns, matrix_t *result)
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

int validate_matrix(int matrix_amount, matrix_t *A, ...)
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

double mult_matrix_res(int i, int j, matrix_t *A, matrix_t *B)
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

int mult_number(matrix_t *A, double number, matrix_t *result)
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
  error_code = create_matrix(A->rows, A->columns, result);
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

int mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result)
{
  if (validate_matrix(2, A, B) || result == NULL)
  {
    return INCORRECT_MATRIX;
  }
  else if (A->columns != B->rows)
  {
    return CALCULATION_ERROR;
  }

  int error_code = OK;
  error_code = create_matrix(A->rows, B->columns, result);
  if (!error_code)
  {
    for (int i = 0; i < result->rows && error_code == OK; i++)
    {
      for (int j = 0; j < result->columns && error_code == OK; j++)
      {
        result->matrix[i][j] = mult_matrix_res(i, j, A, B);
        if (!isfinite(result->matrix[i][j]))
        {
          error_code = CALCULATION_ERROR;
        }
      }
    }
  }

  return error_code;
}

int is_matrix_same_size(int matrix_amount, matrix_t *A, ...)
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

int sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result)
{
  if (validate_matrix(2, A, B) || result == NULL)
  {
    return INCORRECT_MATRIX;
  }
  else if (!is_matrix_same_size(2, A, B))
  {
    return CALCULATION_ERROR;
  }

  int error_code = OK;
  error_code = create_matrix(A->rows, A->columns, result);
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

int sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result)
{
  if (validate_matrix(2, A, B) || result == NULL)
  {
    return INCORRECT_MATRIX;
  }
  else if (!is_matrix_same_size(2, A, B))
  {
    return CALCULATION_ERROR;
  }

  int error_code = OK;
  error_code = create_matrix(A->rows, A->columns, result);
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
  create_matrix(n, n, result);
  for (int i = 0; i < n; i++)
  {
    result->matrix[i][i] = 1;
  }
}

void print_matrix(matrix_t *A)
{
  if (A != NULL && A->matrix != NULL)
  {
    for (int i = 0; i < A->rows; i++)
    {
      for (int j = 0; j < A->columns; j++)
      {
        printf("%.4f\t", A->matrix[i][j]);
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

float randomFloat()
{
  double randnumber = rand();
  return (float)(randnumber) / (float)(randnumber + rand());
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
    create_matrix(1, 1, &A_11);
    A_11.matrix[0][0] = A->matrix[0][0];

    matrix_t A_12 = {};
    create_matrix(1, 1, &A_12);
    A_12.matrix[0][0] = A->matrix[0][1];

    matrix_t A_21 = {};
    create_matrix(1, 1, &A_21);
    A_21.matrix[0][0] = A->matrix[1][0];

    matrix_t A_22 = {};
    create_matrix(1, 1, &A_22);
    A_22.matrix[0][0] = A->matrix[1][1];

    matrix_t A_11_inv = {};
    create_matrix(1, 1, &A_11_inv);
    inverse(&A_11, &A_11_inv);

    matrix_t A_21_A_11_inv = {};
    create_matrix(1, 1, &A_21_A_11_inv);

    mult_matrix(&A_21, &A_11_inv, &A_21_A_11_inv);

    matrix_t A_21_A_11_inv_A_12 = {};
    create_matrix(1, 1, &A_21_A_11_inv_A_12);

    mult_matrix(&A_21_A_11_inv, &A_12, &A_21_A_11_inv_A_12);

    matrix_t S_22 = {};
    create_matrix(1, 1, &S_22);

    sub_matrix(&A_22, &A_21_A_11_inv_A_12, &S_22);

    matrix_t S_22_inv = {};
    create_matrix(1, 1, &S_22_inv);
    inverse(&S_22, &S_22_inv);

    matrix_t A_12_S_22_inv = {};
    create_matrix(1, 1, &A_12_S_22_inv);
    mult_matrix(&A_12, &S_22_inv, &A_12_S_22_inv);

    matrix_t A_12_S_22_inv_A_21 = {};
    create_matrix(1, 1, &A_12_S_22_inv_A_21);

    mult_matrix(&A_12_S_22_inv, &A_21, &A_12_S_22_inv_A_21);

    matrix_t A_12_S_22_inv_A_21_A_11_inv = {};
    create_matrix(1, 1, &A_12_S_22_inv_A_21_A_11_inv);

    mult_matrix(&A_12_S_22_inv_A_21, &A_11_inv, &A_12_S_22_inv_A_21_A_11_inv);

    A_12_S_22_inv_A_21_A_11_inv.matrix[0][0] += 1;

    matrix_t B_11 = {};
    create_matrix(1, 1, &B_11);
    mult_matrix(&A_11_inv, &A_12_S_22_inv_A_21_A_11_inv, &B_11);

    matrix_t A_11_inv_A_12 = {};
    create_matrix(1, 1, &A_11_inv_A_12);
    mult_matrix(&A_11_inv, &A_12, &A_11_inv_A_12);
    matrix_t B_12 = {};
    create_matrix(1, 1, &B_12);
    mult_matrix(&A_11_inv_A_12, &S_22_inv, &B_12);

    B_12.matrix[0][0] *= -1;

    matrix_t S_22_inv_A_21 = {};
    create_matrix(1, 1, &S_22_inv_A_21);
    mult_matrix(&S_22_inv, &A_21, &S_22_inv_A_21);

    matrix_t B_21 = {};
    create_matrix(1, 1, &B_21);
    mult_matrix(&S_22_inv_A_21, &A_11_inv, &B_21);
    B_21.matrix[0][0] *= -1;

    matrix_t B_22 = {};
    create_matrix(1, 1, &B_22);

    B_22.matrix[0][0] = S_22_inv.matrix[0][0];

    merge_matrices(&B_11, &B_12, &B_21, &B_22, inv_A);
  }

  else
  {
    int split_point = A->rows / 2;
    matrix_t A_11 = {};
    create_matrix(split_point, split_point, &A_11);

    matrix_t A_12 = {};
    create_matrix(split_point, split_point, &A_12);

    matrix_t A_21 = {};
    create_matrix(split_point, split_point, &A_21);

    matrix_t A_22 = {};
    create_matrix(split_point, split_point, &A_22);

    split_matrix_into_quadrants(A, &A_11, &A_12, &A_21, &A_22);

    matrix_t id = {};
    create_matrix(A_11.rows, A_11.rows, &id);
    get_identity_matrix(A_11.rows, &id);

    matrix_t A_11_inv = {};
    create_matrix(split_point, split_point, &A_11_inv);

    inverse(&A_11, &A_11_inv);

    matrix_t A_21_A_11_inv = {};
    create_matrix(split_point, split_point, &A_21_A_11_inv);

    mult_matrix(&A_21, &A_11_inv, &A_21_A_11_inv);

    matrix_t A_21_A_11_inv_A_12 = {};
    create_matrix(split_point, split_point, &A_21_A_11_inv_A_12);

    mult_matrix(&A_21_A_11_inv, &A_12, &A_21_A_11_inv_A_12);

    matrix_t S_22 = {};
    create_matrix(split_point, split_point, &S_22);

    sub_matrix(&A_22, &A_21_A_11_inv_A_12, &S_22);

    matrix_t S_22_inv = {};
    create_matrix(split_point, split_point, &S_22_inv);
    inverse(&S_22, &S_22_inv);

    matrix_t A_12_S_22_inv = {};
    create_matrix(split_point, split_point, &A_12_S_22_inv);
    mult_matrix(&A_12, &S_22_inv, &A_12_S_22_inv);

    matrix_t A_12_S_22_inv_A_21 = {};
    create_matrix(split_point, split_point, &A_12_S_22_inv_A_21);

    mult_matrix(&A_12_S_22_inv, &A_21, &A_12_S_22_inv_A_21);

    matrix_t A_12_S_22_inv_A_21_A_11_inv = {};
    create_matrix(split_point, split_point, &A_12_S_22_inv_A_21_A_11_inv);

    mult_matrix(&A_12_S_22_inv_A_21, &A_11_inv, &A_12_S_22_inv_A_21_A_11_inv);

    matrix_t A_12_S_22_inv_A_21_A_11_inv_id = {};
    create_matrix(split_point, split_point, &A_12_S_22_inv_A_21_A_11_inv_id);

    sum_matrix(&A_12_S_22_inv_A_21_A_11_inv, &id, &A_12_S_22_inv_A_21_A_11_inv_id);

    matrix_t B_11 = {};
    create_matrix(split_point, split_point, &B_11);
    mult_matrix(&A_11_inv, &A_12_S_22_inv_A_21_A_11_inv_id, &B_11);

    matrix_t A_11_inv_A_12 = {};
    create_matrix(split_point, split_point, &A_11_inv_A_12);
    mult_matrix(&A_11_inv, &A_12, &A_11_inv_A_12);
    matrix_t B_12d = {};
    create_matrix(split_point, split_point, &B_12d);
    mult_matrix(&A_11_inv_A_12, &S_22_inv, &B_12d);

    matrix_t B_12 = {};
    create_matrix(split_point, split_point, &B_12);
    mult_number(&B_12d, -1, &B_12);

    matrix_t S_22_inv_A_21 = {};
    create_matrix(split_point, split_point, &S_22_inv_A_21);
    mult_matrix(&S_22_inv, &A_21, &S_22_inv_A_21);

    matrix_t B_21d = {};
    create_matrix(split_point, split_point, &B_21d);
    mult_matrix(&S_22_inv_A_21, &A_11_inv, &B_21d);

    matrix_t B_21 = {};
    create_matrix(split_point, split_point, &B_21);
    mult_number(&B_21d, -1, &B_21);
    // B_21.matrix[0][0] *= -1;

    matrix_t B_22 = {};
    create_matrix(split_point, split_point, &B_22);

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
int main()
{

  clock_t start_time, end_time;
  double cpu_time;

  matrix_t A = {};
  int n;
  printf("\nChoose Matrix Dimension for the square matrix: ");
  scanf("%d", &n);
  printf("\n");
  int next_size = pow(2, ceil(log(n) / log(2)));
  int adjusted_size = next_size - n;

  if (adjusted_size == 0)
  {

    create_matrix(n, n, &A);

    // Note: You can either enter the matrix input manually, or choose to fill the matrices automatically. But make sure to comment "RequestInput" function calls, and to uncomment the alternative random input option.
    RequestInput("A", &A, n);
    // This is an alternative random input option, but make sure to comment "RequestInput" function calls.
    /*for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        A.matrix[i][j] = randomFloat();
      }
    }*/

    printf("Generate elements...\n");
    printf("\n");
    printf("Matrix A:\n");
    print_matrix(&A);
    printf("\n");
    matrix_t inv_A = {};
    create_matrix(n, n, &inv_A);

    start_time = clock();
    // Perform inversion
    inverse(&A, &inv_A);
    end_time = clock();

    printf("Matrix A inverse: \n");
    print_matrix(&inv_A);

    printf("\n");
    printf("Time taken for Strassen's Algorithm: %.6f seconds\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

    matrix_t A_inv_A = {};

    create_matrix(n, n, &A_inv_A);

    mult_matrix(&A, &inv_A, &A_inv_A);
  }
  else
  {

    create_matrix(next_size, next_size, &A);

    // Note: You can either enter the matrix input manually, or choose to fill the matrices automatically. But make sure to comment "RequestInput" function calls, and to uncomment the alternative random input option.
    RequestInput("A", &A, n);

    // This is an alternative random input option, but make sure to comment "RequestInput" function calls.
    /*for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
      {
        A.matrix[i][j] = randomFloat();
      }
    }*/

    printf("Generate elements and Padding Elements...\n");

    // padding
    for (int i = n; i < next_size; i++)
    {
      for (int j = n; j < next_size; j++)
      {
        if (i == j)
        {
          A.matrix[i][j] = 1;
        }
        else
        {
          A.matrix[i][j] = 0;
        }
      }
    }

    printf("Matrix A:\n");
    print_matrix(&A);
    printf("\n");
    matrix_t inv_A = {};
    create_matrix(next_size, next_size, &inv_A);

    start_time = clock();
    inverse(&A, &inv_A);
    end_time = clock();

    matrix_t orig_A = {};
    create_matrix(n, n, &orig_A);

    matrix_t orig_inv_A = {};
    create_matrix(n, n, &orig_inv_A);

    // store original matrix
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        orig_A.matrix[i][j] = A.matrix[i][j];
        orig_inv_A.matrix[i][j] = inv_A.matrix[i][j];
      }
    }

    printf("Matrix A inverse:\n");
    print_matrix(&orig_inv_A);

    printf("\n");
    printf("Time taken for Strassen's Inversion using Naive multiplication Algorithm: %.6f seconds\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

    matrix_t A_inv_A = {};

    create_matrix(n, n, &A_inv_A);
    mult_matrix(&orig_A, &orig_inv_A, &A_inv_A);
  }
}