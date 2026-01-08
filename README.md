# Numerical and Symbolic Algorithms Modeling

## Project Overview

The goal of the project is to design, implement, and evaluate core numerical linear algebra algorithms, with a particular focus on correctness, performance, and algorithmic complexity.

The project includes multiple implementations of matrix operations, comparisons between classical and advanced algorithms, and an experimental performance evaluation.

---

## Project Structure

The project is organized around three main algorithmic components:

1. **Matrix Multiplication**
2. **LU Decomposition (Gaussian Elimination)**
3. **Matrix Inversion**

Each component includes both implementation and benchmarking code to analyze computational performance.

---

## Implemented Algorithms

### 1. Matrix Multiplication

Two different matrix multiplication approaches are implemented and compared.

#### 1.1 Naive Matrix Multiplication

* Classical triple-loop algorithm
* Time complexity: **O(n³)** for two n × n matrices
* Used as a baseline for correctness and performance comparison

#### 1.2 Strassen Matrix Multiplication

* Divide-and-conquer algorithm based on matrix block decomposition
* Asymptotic complexity: **O(n^log₂7)**
* Recursive implementation with a threshold-based base case
* Demonstrates reduced arithmetic complexity for large matrices

---

### 2. LU Decomposition (Gaussian Elimination)

* Decomposes a matrix **A = L × U**, where:

  * **L** is a lower triangular matrix
  * **U** is an upper triangular matrix
* Based on Gaussian elimination
* Time complexity: **O(n³)**
* Serves as a foundation for matrix inversion and system solving

---

### 3. Matrix Inversion

Multiple approaches to matrix inversion are implemented:

#### 3.1 LU-based Matrix Inversion

* Computes the inverse using LU decomposition
* Relies on forward and backward substitution

#### 3.2 Strassen-based Inversion (with Strassen multiplication)

* Uses block matrix inversion formulas
* Employs Strassen multiplication internally

#### 3.3 Strassen-based Inversion (with Naive multiplication)

* Same block inversion strategy
* Uses classical multiplication for comparison

---

## Performance Evaluation

A detailed benchmarking study is conducted to compare:

* Naive vs Strassen matrix multiplication
* Different inversion strategies
* LU decomposition performance

Metrics analyzed include:

* Execution time
* Scalability with matrix size
* Algorithmic crossover points

Experimental results are presented and discussed in the accompanying report.

---

## Build and Execution

The project is written in **C** and compiled using a standard C compiler.

Typical compilation example:

```bash
gcc -O2 -o model main.c
```

Execution:

```bash
./model
```

(Exact compilation and execution commands may vary depending on the selected module.)

---

## Key Learning Outcomes

* Practical implementation of numerical linear algebra algorithms
* Comparison between asymptotically optimal and classical methods
* Performance measurement and analysis
* Memory management and recursive algorithm design in C

---

## Notes

* This repository is part of a **group academic project**.
* The code is intended for educational and experimental purposes.
* Results and conclusions are discussed in the final project report.

---

## License

This project is provided for academic use only. Redistribution or reuse should properly credit all contributors.
