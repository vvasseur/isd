/******************************************************************************
 *
 *                 M4RI: Linear Algebra over GF(2)
 *
 *    Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
 *    Copyright (C) 2007,2008 Martin Albrecht <malb@informatik.uni-bremen.de>
 *    Copyright (C) 2020 Valentin Vasseur
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  version 2 or higher.
 *
 *    This code is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    General Public License for more details.
 *
 *  The full text of the GPL is available at:
 *
 *                  http://www.gnu.org/licenses/
 ******************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H
#include <stddef.h>
#include <stdint.h>

#define MAX_K 7

typedef uint64_t word_t;
typedef word_t **matrix_t;

#define WORD_SIZE (8 * sizeof(word_t))

matrix_t matrix_alloc(size_t rows, size_t cols);
void matrix_reset(matrix_t M, size_t rows, size_t cols);
void matrix_free(matrix_t M, size_t rows);
void matrix_swap_rows(matrix_t M, size_t i, size_t j);
void matrix_swap_cols(matrix_t M, size_t i, size_t j, size_t rows);
void matrix_alloc_gray_code(int ***rev, int ***diff);
void matrix_free_gray_code(int **rev, int **diff);
void matrix_build_gray_code(int **rev, int **diff);
size_t matrix_gauss_submatrix(matrix_t M, size_t r, size_t c, size_t rows,
                              size_t cols, size_t k);
size_t matrix_echelonize_partial(matrix_t M, size_t rows, size_t cols, size_t k,
                                 size_t rstop, uint64_t *xor_rows, int **rev,
                                 int **diff);
size_t matrix_opt_k(size_t a, size_t b);
#endif /* MATRIX_H */
