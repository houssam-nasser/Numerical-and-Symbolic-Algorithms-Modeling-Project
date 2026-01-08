all: LU_decomposition.c LU_inverse.c Naive_matrix_multiplication.c Strassen_inverse_using_naive_multiplication.c Strassen_inverse_using_strassen_multiplication.c Strassen_multiplication.c
	gcc -O3 -o LU_decomposition LU_decomposition.c -lm
	gcc -O3 -o Naive_matrix_multiplication Naive_matrix_multiplication.c -lm
	gcc -O3 -o Strassen_multiplication Strassen_multiplication.c -lm
	gcc -O3 -o LU_inverse LU_inverse.c -lm
	gcc -O3 -o Strassen_inverse_using_strassen_multiplication Strassen_inverse_using_strassen_multiplication.c -lm
	gcc -O3 -o Strassen_inverse_using_naive_multiplication Strassen_inverse_using_naive_multiplication.c -lm

clean:
	rm -f LU_decomposition Naive_matrix_multiplication Strassen_multiplication LU_inverse Strassen_inverse_using_strassen_multiplication Strassen_inverse_using_naive_multiplication