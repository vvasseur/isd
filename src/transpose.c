/*
   Copyright (c) 2021 Valentin Vasseur

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to
   deal in the Software without restriction, including without limitation the
   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
   sell copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
   IN THE SOFTWARE
*/
#include "transpose.h"

#ifdef __aarch64__
// https://github.com/kunpengcompute/AvxToNeon
#include "../../AvxToNeon/avx2neon.h"
#else
#include <immintrin.h>
#endif

/* Transpose a bit-matrix using the vpmovmskb instruction.
 *
 * See Bitshuffle - https://github.com/kiyo-masui/bitshuffle (MIT)
 * Copyright (c) 2014 Kiyoshi Masui (kiyo@physics.ubc.ca) */
void matrix_transpose(matrix_t At, const matrix_t A, const size_t nrows,
                      const size_t ncols) {
  __m256i vec_x;

  int8_t in[32] __attribute__((aligned(32)));
  for (size_t i = 0; i < (nrows + 31) / 32; i += 1) {
    for (size_t j = 0; j < (ncols + 7) / 8; j += 1) {
      for (size_t k = 0; k < 32; ++k) {
        in[k] = (((int8_t **)A)[k + i * 32])[j];
      }

      vec_x = _mm256_load_si256((__m256i*)in);
      for (size_t k = 8; k-- > 0;) {
        int32_t hi;
        hi = _mm256_movemask_epi8(vec_x);
        vec_x = _mm256_slli_epi64(vec_x, 1);
        (((int32_t **)At)[k + j * 8])[i] = hi;
      }
    }
  }
}

/* Reverse the rows of a matrix.
 *
 * Modify the matrix in place. */
void matrix_reverse_rows(matrix_t A, const size_t nrows) {
  for (size_t i = 0; i < nrows / 2; ++i) {
    word_t *swp = A[i];
    A[i] = A[nrows - 1 - i];
    A[nrows - 1 - i] = swp;
  }
}

/* It is equivalent to transposing then reversing the columns or, equivalently,
 * reversing the rows then transposing.
 *
 * The inverse of this function is 'matrix_transpose_rev_rows'. */
void matrix_transpose_rev_cols(matrix_t At, matrix_t A, const size_t nrows,
                               const size_t ncols) {
  matrix_reverse_rows(A, nrows);
  matrix_transpose(At, A, nrows, ncols);
}

/* It is equivalent to transposing then reversing the rows or, equivalently,
 * reversing the columns then transposing.
 *
 * The inverse of this function is 'matrix_transpose_rev_cols'. */
void matrix_transpose_rev_rows(matrix_t At, matrix_t A, const size_t nrows,
                               const size_t ncols) {
  matrix_transpose(At, A, nrows, ncols);
  matrix_reverse_rows(At, ncols);
}
