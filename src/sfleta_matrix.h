#ifndef SRC_sfleta_MATRIX_H_
#define SRC_sfleta_MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#define SUCCESS 1
#define FAILURE 0

typedef enum {
    CORRECT_MATRIX = 0,
    INCORRECT_MATRIX = 1,
    IDENTITY_MATRIX = 2,
    ZERO_MATRIX = 3
} matrix_type_t;


typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
    matrix_type_t matrix_type;
} matrix_t;

matrix_t sfleta_create_matrix(int rows, int columns);
void sfleta_remove_matrix(matrix_t *A);
int sfleta_eq_matrix(matrix_t *A, matrix_t *B);
matrix_t sfleta_sum_matrix(matrix_t *A, matrix_t *B);
matrix_t sfleta_sub_matrix(matrix_t *A, matrix_t *B);
matrix_t sfleta_mult_number(matrix_t *A, double number);
matrix_t sfleta_mult_matrix(matrix_t *A, matrix_t *B);
matrix_t sfleta_transpose(matrix_t *A);
matrix_t sfleta_calc_complements(matrix_t *A);
double sfleta_determinant(matrix_t *A);
matrix_t sfleta_inverse_matrix(matrix_t *A);

////// internal_functions //////
matrix_t getSubMatrix(matrix_t *src, int rowForDelete, int columnForDelete);
matrix_type_t GetMatrixType(matrix_t *A);
bool isIdentityMatrix(matrix_t *A);
bool isZeroMatrix(matrix_t *A);
int isEquale(double src1, double src2);

#endif  // SRC_sfleta_MATRIX_H_
