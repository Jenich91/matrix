#include <check.h>
#include <stdio.h>
#include <stdlib.h>

#include "sfleta_matrix.h"
void generate_matrix(matrix_t *A, double start, double step) {
    double value = start;
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
            A->matrix[i][j] = value;
            value = value + step;
        }
    }
    A->matrix_type = CORRECT_MATRIX;
    if (start == 0. && step == 0.) {
    A->matrix_type = ZERO_MATRIX;
    }
    if (start == 1. && step == 0. && A->columns == 1 && A->rows == 1) {
        A->matrix_type = IDENTITY_MATRIX;
    }
}

void fill_matrix_from_arr(matrix_t *dst, int columns, double src[][columns]) {
    for (int i = 0; i < dst->rows; i++) {
        for (int j = 0; j < dst->columns; j++) {
            dst->matrix[i][j] = src[i][j];
        }
    }

    dst->matrix_type = GetMatrixType(dst);
}

START_TEST(test_create_and_remove_matrix) {
    int rows = 5;
    int columns = 5;
    matrix_t A = sfleta_create_matrix(rows, columns);
    ck_assert_int_eq(A.rows, rows);
    ck_assert_int_eq(A.columns, columns);
    ck_assert_int_eq(A.matrix_type, ZERO_MATRIX);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
        ck_assert_double_eq(A.matrix[i][j], 0);
        }
    }
    sfleta_remove_matrix(&A);

    rows = 1;
    columns = 1;
    A = sfleta_create_matrix(rows, columns);
    A.matrix[0][0] = 5;
    A.matrix_type = GetMatrixType(&A);
    ck_assert_int_eq(A.rows, rows);
    ck_assert_int_eq(A.columns, columns);
    ck_assert_int_eq(A.matrix_type, CORRECT_MATRIX);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
        ck_assert_double_eq(A.matrix[i][j], 5);
        }
    }
    sfleta_remove_matrix(&A);

    rows = 1;
    columns = 1;
    A = sfleta_create_matrix(rows, columns);
    A.matrix[0][0] = 1;
    A.matrix_type = GetMatrixType(&A);
    ck_assert_int_eq(A.rows, rows);
    ck_assert_int_eq(A.columns, columns);
    ck_assert_int_eq(A.matrix_type, IDENTITY_MATRIX);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
        ck_assert_double_eq(A.matrix[i][j], 1);
        }
    }
    sfleta_remove_matrix(&A);

    rows = 0;
    columns = 0;
    A = sfleta_create_matrix(rows, columns);
    ck_assert_int_eq(A.rows, rows);
    ck_assert_int_eq(A.columns, columns);
    ck_assert_int_eq(A.matrix_type, INCORRECT_MATRIX);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
        ck_assert_double_eq(A.matrix[i][j], 0);
        }
    }
    sfleta_remove_matrix(&A);

    rows = -0;
    columns = -1;
    A = sfleta_create_matrix(rows, columns);
    ck_assert_int_eq(A.matrix_type, INCORRECT_MATRIX);
    sfleta_remove_matrix(&A);
}
END_TEST

START_TEST(test_eq_matrix) {
    double arr1[3][3] = {
                        {1, 0, -3},
                        {0, 1, 2},
                        {0, 1, 5}
                        };

    double arr2[3][3] = {
                        {1, 0, -3},
                        {0, 1, 2},
                        {0, 1, 5}
                        };


    matrix_t testMatrix1 = sfleta_create_matrix(3, 3);
    fill_matrix_from_arr(&testMatrix1, 3, arr1);
    matrix_t testMatrix2 = sfleta_create_matrix(3, 3);
    fill_matrix_from_arr(&testMatrix2, 3, arr2);
    ck_assert_int_eq(sfleta_eq_matrix(&testMatrix1, &testMatrix2), SUCCESS);
    sfleta_remove_matrix(&testMatrix1);
    sfleta_remove_matrix(&testMatrix2);

    double arr3[2][3] = {
                        {1, 0, -3},
                        {0, 1, 2}
                        };

    double arr4[2][4] = {
                        {1, 0, -3, 5},
                        {0, 1, 2, 4}
                        };


    matrix_t testMatrix3 = sfleta_create_matrix(2, 3);
    fill_matrix_from_arr(&testMatrix3, 3, arr3);
    matrix_t testMatrix4 = sfleta_create_matrix(2, 4);
    fill_matrix_from_arr(&testMatrix4, 4, arr4);
    ck_assert_int_eq(sfleta_eq_matrix(&testMatrix3, &testMatrix4), FAILURE);
    sfleta_remove_matrix(&testMatrix3);
    sfleta_remove_matrix(&testMatrix4);

        double arr5[2][3] = {
                        {1, 0, -3},
                        {0, 1, 2}
                        };

    double arr6[2][3] = {
                        {1, 0, -3},
                        {0, 666, 2}
                        };

    matrix_t testMatrix5 = sfleta_create_matrix(2, 3);
    fill_matrix_from_arr(&testMatrix5, 3, arr5);
    matrix_t testMatrix6 = sfleta_create_matrix(2, 3);
    fill_matrix_from_arr(&testMatrix6, 3, arr6);
    ck_assert_int_eq(sfleta_eq_matrix(&testMatrix5, &testMatrix6), FAILURE);
    sfleta_remove_matrix(&testMatrix5);
    sfleta_remove_matrix(&testMatrix6);

        double arr7[2][3] = {
                        {1, 0, -3},
                        {0, 1, 2}
                        };

    matrix_t testMatrix7 = sfleta_create_matrix(2, 3);
    fill_matrix_from_arr(&testMatrix7, 3, arr7);
    matrix_t testMatrix8 = sfleta_create_matrix(2, -3);
    ck_assert_int_eq(sfleta_eq_matrix(&testMatrix7, &testMatrix8), FAILURE);
    sfleta_remove_matrix(&testMatrix7);
    sfleta_remove_matrix(&testMatrix8);
}
END_TEST

START_TEST(test_sum_matrix) {
    int rows = 2;
    int columns = 2;
    double a = 999999.999999;
    double b = -1.999998;
    double c = 999998.000001;
    int type = CORRECT_MATRIX;
    matrix_t A = sfleta_create_matrix(rows, columns);
    matrix_t B = sfleta_create_matrix(rows, columns);
    matrix_t C = sfleta_create_matrix(rows, columns);

    generate_matrix(&A, a, 0);
    generate_matrix(&B, b, 0);
    generate_matrix(&C, c, 0);

    matrix_t C_test = sfleta_sum_matrix(&A, &B);
    ck_assert_int_eq(sfleta_eq_matrix(&C, &C_test), SUCCESS);
    ck_assert_int_eq(C_test.matrix_type, type);

    sfleta_remove_matrix(&A);
    sfleta_remove_matrix(&B);
    sfleta_remove_matrix(&C);
    sfleta_remove_matrix(&C_test);
}
END_TEST

START_TEST(test_sub_matrix) {
    int rows = 100;
    int columns = 300;
    double a = 1;
    double b = 0.999999;
    double c = 0.000001;
    int type = CORRECT_MATRIX;
    matrix_t A = sfleta_create_matrix(rows, columns);
    matrix_t B = sfleta_create_matrix(rows, columns);
    matrix_t C = sfleta_create_matrix(rows, columns);

    generate_matrix(&A, a, 0);
    generate_matrix(&B, b, 0);
    generate_matrix(&C, c, 0);

    matrix_t C_test = sfleta_sub_matrix(&A, &B);
    ck_assert_int_eq(sfleta_eq_matrix(&C, &C_test), SUCCESS);
    ck_assert_int_eq(C_test.matrix_type, type);

    sfleta_remove_matrix(&A);
    sfleta_remove_matrix(&B);
    sfleta_remove_matrix(&C);
    sfleta_remove_matrix(&C_test);
}
END_TEST

START_TEST(test_mult_number) {
    int rows = 2000;
    int columns = 20;
    double a = 1;
    double b = 500;
    int type = CORRECT_MATRIX;
    matrix_t A = sfleta_create_matrix(rows, columns);
    matrix_t B = sfleta_create_matrix(rows, columns);

    generate_matrix(&A, a, 0);
    generate_matrix(&B, a * b, 0);

    matrix_t B_test = sfleta_mult_number(&A, b);
    ck_assert_int_eq(sfleta_eq_matrix(&B, &B_test), SUCCESS);
    ck_assert_int_eq(B_test.matrix_type, type);

    sfleta_remove_matrix(&A);
    sfleta_remove_matrix(&B);
    sfleta_remove_matrix(&B_test);
}
END_TEST

START_TEST(test_mult_matrix) {
    int rows = 100;
    int columns = 200;
    matrix_t A = sfleta_create_matrix(rows, columns);
    matrix_t B = sfleta_create_matrix(columns, rows);
    matrix_t C = sfleta_create_matrix(rows, rows);

    generate_matrix(&A, 1, 0);
    generate_matrix(&B, 1, 0);
    generate_matrix(&C, A.columns, 0);

    matrix_t C_test = sfleta_mult_matrix(&A, &B);
    ck_assert_int_eq(C.rows, A.rows);
    ck_assert_int_eq(C.columns, B.columns);
    ck_assert_int_eq(sfleta_eq_matrix(&C, &C_test), SUCCESS);

    sfleta_remove_matrix(&A);
    sfleta_remove_matrix(&B);
    sfleta_remove_matrix(&C);
    sfleta_remove_matrix(&C_test);
}
END_TEST

START_TEST(test_transpose) {
    int rows = 1000;
    int columns = 2000;
    double a = 1111111;
    double b = 000001;
    matrix_t A = sfleta_create_matrix(rows, columns);
    generate_matrix(&A, a, b);

    matrix_t B = sfleta_transpose(&A);
    matrix_t C = sfleta_transpose(&B);

    ck_assert_int_eq(sfleta_eq_matrix(&A, &C), SUCCESS);
    ck_assert_int_eq(A.matrix_type, C.matrix_type);
    ck_assert_int_eq(A.rows, B.columns);
    ck_assert_int_eq(A.columns, B.rows);

    sfleta_remove_matrix(&A);
    sfleta_remove_matrix(&B);
    sfleta_remove_matrix(&C);
}
END_TEST

START_TEST(test_determinant) {
    int rows = 2;
    int value = 100.;
    double result = 20000.;
    matrix_t A = sfleta_create_matrix(rows, rows);
    double eps = 1e-7;
    for (int i = 0; i < A.rows; i++) {
    A.matrix[i][i] = value;
    }
    ck_assert_double_le(sfleta_determinant(&A) - result, eps);
    sfleta_remove_matrix(&A);
}
END_TEST

START_TEST(test_calc_complements_and_inverse_matrix) {
    int rows = 4;
    double start = 5;
    double step = 0.1;
    matrix_t A = sfleta_create_matrix(rows, rows);
    generate_matrix(&A, start, step);
    A.matrix[0][0] = 0;
    A.matrix[rows - 1][rows - 1] = A.matrix[rows - 1][rows - 1] * 10;
    matrix_t I = sfleta_inverse_matrix(&A);
    matrix_t M = sfleta_mult_matrix(&A, &I);
    ck_assert_int_eq(M.matrix_type, IDENTITY_MATRIX);
    sfleta_remove_matrix(&A);
    sfleta_remove_matrix(&I);
    sfleta_remove_matrix(&M);

    matrix_t B = sfleta_create_matrix(1, 1);
    B.matrix[0][0] = 0;
    I = sfleta_inverse_matrix(&B);
    ck_assert_int_eq(I.matrix_type, INCORRECT_MATRIX);
    sfleta_remove_matrix(&B);
    sfleta_remove_matrix(&I);

    B = sfleta_create_matrix(1, 1);
    B.matrix[0][0] = 5;
    I = sfleta_inverse_matrix(&B);
    M = sfleta_mult_matrix(&B, &I);
    ck_assert_int_eq(M.matrix_type, IDENTITY_MATRIX);
    sfleta_remove_matrix(&B);
    sfleta_remove_matrix(&I);
    sfleta_remove_matrix(&M);
}
END_TEST

Suite *sfleta_suite(void) {
    Suite *suite = suite_create("Matrix");
    TCase *tcase_core = tcase_create("Equelity");

    tcase_add_test(tcase_core, test_create_and_remove_matrix);
    tcase_add_test(tcase_core, test_eq_matrix);
    tcase_add_test(tcase_core, test_sum_matrix);
    tcase_add_test(tcase_core, test_sub_matrix);
    tcase_add_test(tcase_core, test_mult_number);
    tcase_add_test(tcase_core, test_mult_matrix);
    tcase_add_test(tcase_core, test_transpose);
    tcase_add_test(tcase_core, test_determinant);
    tcase_add_test(tcase_core, test_calc_complements_and_inverse_matrix);
    suite_add_tcase(suite, tcase_core);
    return suite;
}

int main(void) {
    Suite *suite = sfleta_suite();
    SRunner *suite_runner = srunner_create(suite);
    srunner_run_all(suite_runner, CK_NORMAL);
    int failed_count = srunner_ntests_failed(suite_runner);
    srunner_free(suite_runner);
    return (failed_count != 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}
