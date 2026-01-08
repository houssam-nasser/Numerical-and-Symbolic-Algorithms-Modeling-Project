# Matrix Manipulation Programs
This repository contains C programs for various matrix manipulation operations including LU decomposition, matrix multiplication, and matrix inversion using different algorithms.

## Compilation
To compile all the programs, run the following command:
make all

## Clean up
To clean all the compiled programs, run the following command:
make clean

## Changing matrices input

// Example 

// Note: You can either enter the matrix input manually, or choose to fill the matrices automatically. But make sure to comment "RequestInput" function calls and to uncomment the alternative random input option.
// Step 3: Asking the user to input Matrix A and Matrix B elements
RequestInput("A", A, size);

// This is an alternative random input option, but make sure to comment "RequestInput" function calls.
/*for (int i = 0; i < size; i++)
 {
 for (int j = 0; j < size; j++)
 {
  A[i][j] = (i == j) ? 1.0 : ((double)rand() / RAND_MAX);
  }
}*/
